# STATCOM mitigation on a D-Suite feeder (snapshot)
# - Baseline PF -> weakest bus
# - Add STATCOM (reactive-only generator)
# - Sweep STATCOM rating and solve OPF
# - add inverter-loss branches using rosetta_distribution_opf (RPMD)

import Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate()

using PowerModelsDistribution
using JuMP
using Ipopt
using DataFrames
using CSV
using Plots

const PMD = PowerModelsDistribution
PMD.silence!()

# -----------------------------
# 0) User edits
# -----------------------------
# ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"
ROOT = joinpath(@__DIR__, "..", "..")
FEEDER = "spd_s"

master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", FEEDER, "master_scaled.dss")

OUTDIR = joinpath(ROOT, "results", "statcom_mitigation", FEEDER)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

LOAD_ALPHA = 1.0
VMIN_PU = 0.94
VMAX_PU = 1.10

# sweep sizes (kVAr per phase as a human-friendly knob)
QCAP_KVAR_LIST = [0.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0]

USE_RPMD_LOSSES = false

# -----------------------------
# Helpers
# -----------------------------
count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0
normalize_bus(s::AbstractString) = lowercase(split(String(s), ".")[1])

function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng
    for (_, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        for k in ("pd_nom","qd_nom","pd","qd")
            if haskey(ld, k)
                v = ld[k]
                ld[k] = v isa Number ? alpha*v : alpha .* v
            end
        end
    end
    return eng
end

function set_bus_voltage_limits!(math; vmin=0.94, vmax=1.10)
    for (_, bus_any) in math["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vmin"] = vmin .* ones(length(bus["vmin"]))
        bus["vmax"] = vmax .* ones(length(bus["vmax"]))
    end
end

# convert kVAr -> per-unit-ish "model" units when settings exist
function kvar_to_model_q(math::Dict{String,Any}, kvar::Real)
    settings = get(math, "settings", Dict{String,Any}())
    sbase = get(settings, "sbase", missing)
    psf   = get(settings, "power_scale_factor", missing)

    if (sbase isa Number) && (psf isa Number) && isfinite(sbase) && isfinite(psf) && sbase > 0 && psf > 0
        sbase_va = Float64(sbase) * Float64(psf)
        q_va = Float64(kvar) * 1000.0
        return q_va / sbase_va
    end
    return Float64(kvar)
end

function run_pf(math::Dict{String,Any})
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
    return PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)
end

function run_opf(math::Dict{String,Any})
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes", "max_iter"=>100000)
    return PMD.solve_mc_opf(math, PMD.IVRUPowerModel, ipopt)
end

function extract_bus_vm(soln::Dict{String,Any}, math::Dict{String,Any})
    sol_bus = soln["solution"]["bus"]
    rows = NamedTuple[]
    for (bus_id_any, bus_any) in math["bus"]
        bus_id = string(bus_id_any)
        haskey(sol_bus, bus_id) || continue

        sb = sol_bus[bus_id]
        vm =
            haskey(sb, "vm") ? sb["vm"] :
            (haskey(sb,"vr") && haskey(sb,"vi")) ? sqrt.(sb["vr"].^2 .+ sb["vi"].^2) :
            nothing
        vm === nothing && continue

        eng_name =
            haskey(bus_any, "name") ? normalize_bus(string(bus_any["name"])) :
            bus_id

        vmin = minimum(vm[1:min(3,length(vm))])
        vmax = maximum(vm[1:min(3,length(vm))])

        push!(rows, (bus_id=bus_id, bus=eng_name, vmin_pu=vmin, vmax_pu=vmax))
    end
    df = DataFrame(rows)
    sort!(df, :vmin_pu)
    return df
end

function count_violations(bus_df::DataFrame; vmin=0.94, vmax=1.10)
    n_under = count(bus_df.vmin_pu .< vmin)
    n_over  = count(bus_df.vmax_pu .> vmax)
    return n_under, n_over
end

# PMD/InfrastructureModels expects numeric-like string ids for components
function next_numeric_id(existing::Dict{String,Any}; start::Int=900001)
    ids = Int[]
    for k in keys(existing)
        x = tryparse(Int, k)
        x === nothing || push!(ids, x)
    end
    mx = isempty(ids) ? (start - 1) : maximum(ids)
    return string(max(mx + 1, start))
end

# Add STATCOM as a reactive-only gen at a specific *math* bus_id
function add_statcom!(math::Dict{String,Any}, bus_id::String; qcap_model::Float64)
    haskey(math, "gen") || (math["gen"] = Dict{String,Any}())

    # pick a template generator that already has the right schema
    # (usually the source gen exists)
    template_id = first(keys(math["gen"]))
    template = deepcopy(math["gen"][template_id])

    bus = math["bus"][bus_id]::Dict{String,Any}

    # determine how many phases to connect (use 3 if available)
    nconn =
        haskey(bus, "terminals") ? min(3, length(bus["terminals"])) :
        haskey(bus, "vmin") ? min(3, length(bus["vmin"])) :
        3

    gen_id = next_numeric_id(math["gen"])

    # overwrite the important bits
    template["name"] = "statcom"
    template["gen_bus"] = bus["index"]
    template["connections"] = collect(1:nconn)

    template["gen_status"] = 1

    # reactive-only device
    template["pmin"] = fill(0.0, nconn)
    template["pmax"] = fill(0.0, nconn)
    template["qmin"] = fill(-qcap_model, nconn)
    template["qmax"] = fill(qcap_model, nconn)

    # optional: kill any weird leftovers if template had them
    if haskey(template, "pg"); template["pg"] = fill(0.0, nconn); end
    if haskey(template, "qg"); template["qg"] = fill(0.0, nconn); end

    # ensure index exists if PMD expects it
    template["index"] = tryparse(Int, gen_id) === nothing ? template["index"] : parse(Int, gen_id)

    math["gen"][gen_id] = template

    return gen_id
end

# -----------------------------
# Main
# -----------------------------
println("Parsing: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

println("eng counts: buses=", count_dict(eng0,"bus"),
        " lines=", count_dict(eng0,"line"),
        " loads=", count_dict(eng0,"load"),
        " transformers=", count_dict(eng0,"transformer"))

scale_loads!(eng0, LOAD_ALPHA)

math0 = PMD.transform_data_model(eng0; multinetwork=false, kron_reduce=true, phase_project=true)
set_bus_voltage_limits!(math0; vmin=VMIN_PU, vmax=VMAX_PU)

pf0 = run_pf(deepcopy(math0))
println("Baseline PF status: ", pf0["termination_status"])

bus_pf = extract_bus_vm(pf0, math0)
n_under0, n_over0 = count_violations(bus_pf; vmin=VMIN_PU, vmax=VMAX_PU)

weak = bus_pf[1, :]
weak_bus_id = string(weak.bus_id)
println("Weakest bus: ", weak.bus, " | vmin_pu=", weak.vmin_pu)

CSV.write(joinpath(TBLDIR, "baseline_bus_vm_sorted.csv"), bus_pf)

rows = NamedTuple[]
for qkvar in QCAP_KVAR_LIST
    math = deepcopy(math0)

    q_model = kvar_to_model_q(math, qkvar)
    gen_id = add_statcom!(math, weak_bus_id; qcap_model=q_model)

    opf = run_opf(math)
    status = string(opf["termination_status"])

    bus_opf = extract_bus_vm(opf, math)
    n_under, n_over = count_violations(bus_opf; vmin=VMIN_PU, vmax=VMAX_PU)

    push!(rows, (
        feeder = FEEDER,
        load_alpha = LOAD_ALPHA,
        statcom_bus = weak.bus,
        statcom_bus_id = weak_bus_id,
        statcom_gen_id = gen_id,
        qcap_kvar_per_phase = Float64(qkvar),
        termination_status = status,
        vmin_pu = minimum(bus_opf.vmin_pu),
        vmax_pu = maximum(bus_opf.vmax_pu),
        n_under = n_under,
        n_over = n_over
    ))
end

df = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "statcom_sweep_summary.csv"), df)

p1 = plot(df.qcap_kvar_per_phase, df.vmin_pu, marker=:circle,
    xlabel="STATCOM rating (kVAr per phase)",
    ylabel="Minimum bus voltage (pu)",
    title="STATCOM sweep | vmin (OPF) | $FEEDER | alpha=$LOAD_ALPHA",
    legend=false
)
hline!([VMIN_PU], linestyle=:dash, color=:red)
savefig(p1, joinpath(FIGDIR, "statcom_sweep_vmin.png"))

p2 = plot(df.qcap_kvar_per_phase, df.vmax_pu, marker=:circle,
    xlabel="STATCOM rating (kVAr per phase)",
    ylabel="Maximum bus voltage (pu)",
    title="STATCOM sweep | vmax (OPF) | $FEEDER | alpha=$LOAD_ALPHA",
    legend=false
)
hline!([VMAX_PU], linestyle=:dash, color=:red)
savefig(p2, joinpath(FIGDIR, "statcom_sweep_vmax.png"))

p3 = plot(df.qcap_kvar_per_phase, df.n_under, marker=:circle,
    xlabel="STATCOM rating (kVAr per phase)",
    ylabel="Under-voltage buses (count)",
    title="STATCOM sweep | under count | $FEEDER",
    legend=false
)
savefig(p3, joinpath(FIGDIR, "statcom_sweep_n_under.png"))

p4 = plot(df.qcap_kvar_per_phase, df.n_over, marker=:circle,
    xlabel="STATCOM rating (kVAr per phase)",
    ylabel="Over-voltage buses (count)",
    title="STATCOM sweep | over count | $FEEDER",
    legend=false
)
savefig(p4, joinpath(FIGDIR, "statcom_sweep_n_over.png"))

println("Baseline violations: under=$n_under0 over=$n_over0")
println("Saved to: ", OUTDIR)
println("Done.")
