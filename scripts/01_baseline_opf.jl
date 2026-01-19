import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PowerModelsDistribution
using JuMP
using Ipopt
using DataFrames
using CSV
using Plots
using Statistics

const PMD = PowerModelsDistribution
PMD.silence!()

# ============================================================
# Baseline snapshot OPF
# with LOAD SCALING
# ============================================================

ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"
NET  = "spd_s"   # USE spd_s so something actually happens

master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_scaled.dss")

OUTDIR = joinpath(ROOT, "results", "baseline_opf", NET)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

VBASE_LN = 230.0
VMIN_PU  = 0.94
VMAX_PU  = 1.10
VMIN_V   = VMIN_PU * VBASE_LN
VMAX_V   = VMAX_PU * VBASE_LN

LOAD_ALPHA = 2.0   # deliberately stress the network

# -----------------------------
# Helpers
# -----------------------------
count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0

function pick_source_bus_eng(eng)::String
    haskey(eng["bus"], "sourcebusz") && return "sourcebusz"
    haskey(eng["bus"], "sourcebusZ") && return "sourcebusZ"
    haskey(eng["bus"], "sourcebus")  && return "sourcebus"
    return first(keys(eng["bus"]))
end

function bus_eng_name(bus_id::String, bus_data::Dict{String,Any})::String
    if haskey(bus_data, "name")
        return lowercase(string(bus_data["name"]))
    end
    if haskey(bus_data, "source_id")
        sid = lowercase(string(bus_data["source_id"]))
        occursin("bus.", sid) && return split(sid, "bus.")[end]
    end
    return lowercase(bus_id)
end

function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng
    for (_, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        for k in ("pd", "qd")
            if haskey(ld, k)
                v = ld[k]
                ld[k] = v isa Number ? alpha * v : alpha .* v
            end
        end
    end
    return eng
end

function solved_bus_vm_volts_keyed_by_eng(result, math; vbase_ln=230.0)
    sol_bus = result["solution"]["bus"]
    out = Dict{String,Dict{String,Any}}()

    for (bus_id_any, bus_data_any) in math["bus"]
        bus_id = string(bus_id_any)
        haskey(sol_bus, bus_id) || continue
        sb = sol_bus[bus_id]

        vm_pu =
            haskey(sb, "vm") ? sb["vm"] :
            haskey(sb, "vr") ? sqrt.(sb["vr"].^2 .+ sb["vi"].^2) :
            nothing
        vm_pu === nothing && continue

        vmV = vm_pu .* vbase_ln
        eng = bus_eng_name(bus_id, bus_data_any)

        out[eng] = Dict(
            "vma" => vmV[1],
            "vmb" => vmV[2],
            "vmc" => vmV[3]
        )
    end
    return out
end

function voltage_stats_pu(buses; vbase_ln=230.0)
    vmins = [min(b["vma"], b["vmb"], b["vmc"]) / vbase_ln for b in values(buses)]
    sort!(vmins)
    return (
        min = minimum(vmins),
        q05 = vmins[clamp(Int(ceil(0.05 * length(vmins))),1,length(vmins))],
        median = vmins[clamp(Int(ceil(0.5  * length(vmins))),1,length(vmins))],
        q95 = vmins[clamp(Int(ceil(0.95 * length(vmins))),1,length(vmins))]
    )
end

# -----------------------------
# 1) Parse model
# -----------------------------
println("Parsing: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

println("eng counts: buses=", count_dict(eng0,"bus"),
        " lines=", count_dict(eng0,"line"),
        " loads=", count_dict(eng0,"load"),
        " transformers=", count_dict(eng0,"transformer"))

scale_loads!(eng0, LOAD_ALPHA)

SOURCE_BUS = pick_source_bus_eng(eng0)
println("SOURCE_BUS = ", SOURCE_BUS)

# -----------------------------
# 2) Build OPF model
# -----------------------------
ipopt = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 0,
    "sb" => "yes"
)

math = PMD.transform_data_model(
    eng0;
    multinetwork=false,
    kron_reduce=true,
    phase_project=true
)

# Give OPF voltage bounds (important)
for (_, bus) in math["bus"]
    bus["vmin"] = VMIN_PU .* ones(length(bus["vmin"]))
    bus["vmax"] = VMAX_PU .* ones(length(bus["vmax"]))
end

println("Running OPF (IVRUPowerModel)...")
result = PMD.solve_mc_opf(math, PMD.IVRUPowerModel, ipopt)

println("OPF status: ", result["termination_status"])
println("Objective value: ", get(result, "objective", missing))

# -----------------------------
# 3) Extract results
# -----------------------------
buses = solved_bus_vm_volts_keyed_by_eng(result, math; vbase_ln=VBASE_LN)
stats = voltage_stats_pu(buses; vbase_ln=VBASE_LN)

println("Voltage stats (pu): ",
    "min=", round(stats.min, digits=4),
    " q05=", round(stats.q05, digits=4),
    " median=", round(stats.median, digits=4),
    " q95=", round(stats.q95, digits=4)
)

# -----------------------------
# 4) Plots
# -----------------------------
va = [b["vma"] for b in values(buses)]
vb = [b["vmb"] for b in values(buses)]
vc = [b["vmc"] for b in values(buses)]

bins = (VMIN_V-5):0.5:(VMAX_V+5)
p = histogram(va; bins, label="A")
histogram!(p, vb; bins, label="B")
histogram!(p, vc; bins, label="C")
xlabel!(p, "Voltage (V)")
ylabel!(p, "Count")
title!(p, "Baseline OPF voltage histogram | $NET | alpha=$LOAD_ALPHA")

savefig(p, joinpath(FIGDIR, "baseline_opf_voltage_hist.png"))

# -----------------------------
# 5) Tables
# -----------------------------
rows = NamedTuple[]
for (bus, b) in buses
    push!(rows, (
        bus = bus,
        vA = b["vma"],
        vB = b["vmb"],
        vC = b["vmc"],
        vmin = min(b["vma"], b["vmb"], b["vmc"]),
        vmax = max(b["vma"], b["vmb"], b["vmc"])
    ))
end

df = DataFrame(rows)
sort!(df, :vmin)
CSV.write(joinpath(TBLDIR, "baseline_opf_bus_voltages.csv"), df)

meta = DataFrame([
    (net=NET, alpha=LOAD_ALPHA,
     vmin_pu=stats.min, vq05_pu=stats.q05,
     vmed_pu=stats.median, vq95_pu=stats.q95)
])
CSV.write(joinpath(TBLDIR, "baseline_opf_summary.csv"), meta)

println("Saved results to: ", OUTDIR)
println("Done.")
