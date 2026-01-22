# Compare baseline vs PV using time-series scaling, but only on selected timesteps.
#
# Steps:
# - Read the baseline PF metrics CSV from 01
# - Pick "PV-stress" timesteps (high PV proxy, low load alpha)
# - For each selected timestep:
#     - Run baseline PF
#     - Add PV (as negative load) in a PMD-safe way
#     - Run PF again
#     - Save comparison metrics + voltage plots

import Pkg
# Pkg.activate(joinpath(@__DIR__, "..", ".."))
# Pkg.instantiate()

using PowerModelsDistribution
using DataFrames
using Dates
using Statistics
using CSV
using Plots
using Ipopt
using Random
using InfrastructureModels
using JuMP
using LinearAlgebra: diag, diagm

const PMD = PowerModelsDistribution
const IM = InfrastructureModels
PMD.silence!()

# using rosetta_distribution_opf
# const RDOPF = rosetta_distribution_opf
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")


## --------------------------------------------------
# 0) Settings
# --------------------------------------------------

# ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"
ROOT = joinpath(@__DIR__, "../..")

include(joinpath(ROOT, "src/read_functions.jl"))

NET_3w = "spd_s"
NET_4w = "spd_s_4w"
master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4w, "master_scaled.dss")

YEAR = 2023
STRIDE = 4

BASELINE_OUTDIR = joinpath(ROOT, "results", "time_series", "baseline_pf", NET_3w, "year=$(YEAR)_stride=$(STRIDE)")
BASELINE_CSV    = joinpath(BASELINE_OUTDIR, "tables", "timeseries_baseline_pf_metrics.csv")

# PV settings
PV_BUS = "AUTO"             # "AUTO" picks a far-ish non-source bus
PV_KW_PER_PHASE = 3.0
PV_PHASES = [1]             # typical single-phase PV for homes

# How many timesteps we re-simulate
N_PV_STRESS = 30
RANDOM_SEED = 7

# Reporting limits (same style as 01)
VMIN_PU = 0.9
VMAX_PU = 1.10
VBASE_LN = 230.0

# Base scaling for loads before timestep-specific alpha
LOAD_ALPHA_BASE = 2.5

OUTDIR = joinpath(
    ROOT,
    "results",
    "time_series",
    "pv_pf",
    NET,
    "year=$(YEAR)_stride=$(STRIDE)",
    "pvkw=$(PV_KW_PER_PHASE)_ph=$(join(PV_PHASES, '_'))_K=$(N_PV_STRESS)"
)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR)
mkpath(TBLDIR)

## --------------------------------------------------
# 1) Same helpers as 01: scale loads, solve PF, metrics
# --------------------------------------------------

function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng
    for (_, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        for k in ("pd_nom", "qd_nom", "pd", "qd")
            if haskey(ld, k)
                v = ld[k]
                ld[k] = v isa Number ? alpha * v : alpha .* v
            end
        end
    end
    return eng
end

function solve_pf(eng::Dict{String,Any})
    # EXACTLY like 01
    math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")
    pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)
    return pf, math
end

function pf_metrics(pf::Dict{String,Any}; vmin_pu=0.9, vmax_pu=1.10)
    sol_bus = pf["solution"]["bus"]
    vmins = Float64[]
    vmaxs = Float64[]

    for sb in values(sol_bus)
        vm =
            haskey(sb, "vm") ? sb["vm"][1:3] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"][1:3].^2 .+ sb["vi"][1:3].^2) :
            nothing
        vm === nothing && continue

        push!(vmins, minimum(vm))
        push!(vmaxs, maximum(vm))
    end

    if isempty(vmins)
        return (status=string(pf["termination_status"]), vmin=NaN, vmax=NaN, n_under=0, n_over=0)
    end

    vmin_sys = minimum(vmins)
    vmax_sys = maximum(vmaxs)
    n_under = count(x -> x < vmin_pu, vmins)
    n_over  = count(x -> x > vmax_pu, vmaxs)

    return (status=string(pf["termination_status"]), vmin=vmin_sys, vmax=vmax_sys, n_under=n_under, n_over=n_over)
end

function print_baseline_load_summary(eng::Dict{String,Any})
    total_pd = 0.0
    total_qd = 0.0
    nloads = 0

    if haskey(eng, "load")
        for (_, ld_any) in eng["load"]
            ld = ld_any::Dict{String,Any}
            if haskey(ld, "pd")
                total_pd += sum(ld["pd"])
            elseif haskey(ld, "pd_nom")
                total_pd += sum(ld["pd_nom"])
            end
            if haskey(ld, "qd")
                total_qd += sum(ld["qd"])
            elseif haskey(ld, "qd_nom")
                total_qd += sum(ld["qd_nom"])
            end
            nloads += 1
        end
    end

    println("\nBaseline load summary from engineering model")
    println("  number of loads = ", nloads)
    println("  total P (kW)    = ", total_pd)
    println("  total Q (kvar)  = ", total_qd)
end

# --------------------------------------------------
# 2) Plot helpers for voltage profiles
# --------------------------------------------------

function make_lines_df_from_eng(eng::Dict{String,Any})
    rows = NamedTuple[]
    if haskey(eng, "line")
        for (_, ln_any) in eng["line"]
            ln = ln_any::Dict{String,Any}
            f_bus = ln["f_bus"]
            t_bus = ln["t_bus"]
            f_phases = ln["f_connections"]
            length_km = (haskey(ln, "length") ? ln["length"] : 0.0) / 1000.0
            push!(rows, (Bus1=f_bus, Bus2=t_bus, phases=f_phases, length_km=length_km))
        end
    elseif haskey(eng, "branch")
        for (_, br_any) in eng["branch"]
            br = br_any::Dict{String,Any}
            f_bus = br["f_bus"]
            t_bus = br["t_bus"]
            f_phases = br["f_connections"]
            length_km = (haskey(br, "length") ? br["length"] : 0.0) / 1000.0
            push!(rows, (Bus1=f_bus, Bus2=t_bus, phases=f_phases, length_km=length_km))
        end
    else
        error("No line or branch data found")
    end
    return DataFrame(rows)
end

function pick_source_bus(eng::Dict{String,Any})
    if haskey(eng, "bus")
        names = Set(keys(eng["bus"]))
        if "sourcebusz" in names
            return "sourcebusz"
        elseif "sourcebus" in names
            return "sourcebus"
        end
    end
    return "sourcebus"
end

function compute_bus_distances(lines_df::DataFrame; source_bus::String)
    adj = Dict{String, Vector{Tuple{String, Float64}}}()
    for r in eachrow(lines_df)
        b1 = r.Bus1
        b2 = r.Bus2
        len = r.length_km
        push!(get!(adj, b1, Tuple{String,Float64}[]), (b2, len))
        push!(get!(adj, b2, Tuple{String,Float64}[]), (b1, len))
    end

    dist = Dict{String,Float64}(source_bus => 0.0)
    queue = [source_bus]
    while !isempty(queue)
        u = popfirst!(queue)
        for (v, w) in get(adj, u, Tuple{String,Float64}[])
            if !haskey(dist, v)
                dist[v] = dist[u] + w
                push!(queue, v)
            end
        end
    end
    return dist
end

function solved_bus_vm_volts_dual(pf::Dict{String,Any}, math::Dict{String,Any}; vbase_ln::Float64)
    sol_bus = pf["solution"]["bus"]

    buses_by_name = Dict{String, Dict{String,Any}}()
    buses_by_id   = Dict{String, Dict{String,Any}}()

    for (bus_id, bus_data_any) in math["bus"]
        bus_data = bus_data_any::Dict{String,Any}
        name = string(bus_data["name"])

        haskey(sol_bus, bus_id) || continue
        sb = sol_bus[bus_id]

        vm_pu =
            haskey(sb, "vm") ? sb["vm"] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"].^2 .+ sb["vi"].^2) :
            nothing
        vm_pu === nothing && continue

        vmV = vm_pu .* vbase_ln

        entry = Dict(
            "vma" => [vmV[1]],
            "vmb" => length(vmV) >= 2 ? [vmV[2]] : [vmV[1]],
            "vmc" => length(vmV) >= 3 ? [vmV[3]] : [vmV[1]]
        )

        buses_by_name[name] = deepcopy(entry)
        buses_by_id[string(bus_id)] = deepcopy(entry)
    end

    return buses_by_name, buses_by_id
end

function plot_voltage_along_feeder_snap(buses_dict::Dict{String,Dict{String,Any}}, lines_df::DataFrame;
        t::Int=1, Vthreshold::Float64=1000.0, vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0)

    p = plot(legend=false)
    xlabel!("Distance from reference bus (km)")
    ylabel!("Voltage magnitude P-N (V)")
    title!("Voltage along feeder")

    colors = Dict(1=>:blue, 2=>:red, 3=>:black)

    for r in eachrow(lines_df)
        b1 = r.Bus1
        b2 = r.Bus2
        phases = r.phases

        if !(haskey(buses_dict, b1) && haskey(buses_dict, b2))
            continue
        end
        if !(haskey(buses_dict[b1], "distance") && haskey(buses_dict[b2], "distance"))
            continue
        end

        for ph in phases
            vm_f = ph == 1 ? buses_dict[b1]["vma"][t] :
                   ph == 2 ? buses_dict[b1]["vmb"][t] :
                             buses_dict[b1]["vmc"][t]
            vm_t = ph == 1 ? buses_dict[b2]["vma"][t] :
                   ph == 2 ? buses_dict[b2]["vmb"][t] :
                             buses_dict[b2]["vmc"][t]

            if vm_f < Vthreshold && vm_t < Vthreshold
                plot!(
                    [buses_dict[b1]["distance"], buses_dict[b2]["distance"]],
                    [vm_f, vm_t],
                    color=colors[ph],
                    marker=:circle,
                    markersize=1
                )
            end
        end
    end

    maxdist = maximum(bus["distance"] for bus in values(buses_dict) if haskey(bus, "distance"))
    plot!([0, maxdist], [vminV, vminV], linestyle=:dash, color=:red)
    plot!([0, maxdist], [vmaxV, vmaxV], linestyle=:dash, color=:red)

    return p
end

function plot_voltage_histogram_snap(buses_dict::Dict{String,Dict{String,Any}};
        t::Int=1, Vthreshold::Float64=1000.0, vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0)

    phase_a = Float64[]
    phase_b = Float64[]
    phase_c = Float64[]

    for (_, bus_data) in buses_dict
        if haskey(bus_data, "vma") && bus_data["vma"][t] < Vthreshold
            push!(phase_a, bus_data["vma"][t])
        end
        if haskey(bus_data, "vmb") && bus_data["vmb"][t] < Vthreshold
            push!(phase_b, bus_data["vmb"][t])
        end
        if haskey(bus_data, "vmc") && bus_data["vmc"][t] < Vthreshold
            push!(phase_c, bus_data["vmc"][t])
        end
    end

    bins = (vminV-5):0.5:(vmaxV+5)

    p = histogram(phase_a; bins=bins, label="phase a")
    histogram!(phase_b; bins=bins, label="phase b")
    histogram!(phase_c; bins=bins, label="phase c")

    xlabel!("Voltage magnitude (V)")
    ylabel!("Count")
    title!("Voltage histogram")

    return p
end

function plot_voltage_combined_snap(buses_dict::Dict{String,Dict{String,Any}}, lines_df::DataFrame;
        t::Int=1, Vthreshold::Float64=1000.0, vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0)

    p1 = plot_voltage_along_feeder_snap(buses_dict, lines_df; t=t, Vthreshold=Vthreshold, vminV=vminV, vmaxV=vmaxV)
    p2 = plot_voltage_histogram_snap(buses_dict; t=t, Vthreshold=Vthreshold, vminV=vminV, vmaxV=vmaxV)
    return plot(p1, p2, layout=(1,2))
end

# --------------------------------------------------
# 3) PV-stress selection (proxy PV + low alpha)
# --------------------------------------------------

# quick PV proxy: only used for ranking timesteps
function pv_shape_simple(time_vec::Vector{DateTime})
    pv = zeros(Float64, length(time_vec))
    for (i, ts) in enumerate(time_vec)
        h = hour(ts) + minute(ts)/60
        x = (h - 6.0)/12.0
        pv[i] = (0.0 <= x <= 1.0) ? sin(pi*x)^2 : 0.0
    end
    m = maximum(pv)
    return m > 0 ? pv ./ m : pv
end

# --------------------------------------------------
# 4) PV injection: PMD-safe negative load
# --------------------------------------------------

# Find an existing load at pv_bus. If found, subtract PV without breaking vector shapes.
# If not found, create a new load with REQUIRED keys, including `model`.
function add_pv_negative_load!(eng::Dict{String,Any}, pv_bus::String, pv_kw_per_phase::Float64, pv_phases::Vector{Int})
    haskey(eng, "load") || (eng["load"] = Dict{String,Any}())

    # Try find an existing load on that bus
    for (lid, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        haskey(ld, "bus") || continue
        string(ld["bus"]) == pv_bus || continue

        # Identify which P field exists
        pd_key = haskey(ld, "pd") ? "pd" : (haskey(ld, "pd_nom") ? "pd_nom" : nothing)
        pd_key === nothing && continue

        pd_val = ld[pd_key]

        # If scalar, subtract total PV
        if pd_val isa Number
            ld[pd_key] = Float64(pd_val) - pv_kw_per_phase * length(pv_phases)
            return ("modified_existing_load", string(lid))
        end

        # If vector, subtract per-phase only if we can
        if pd_val isa AbstractVector
            pd_new = Float64.(pd_val)

            # safest mapping: only subtract on indices that exist
            for ph in pv_phases
                if 1 <= ph <= length(pd_new)
                    pd_new[ph] -= pv_kw_per_phase
                end
            end

            ld[pd_key] = pd_new
            return ("modified_existing_load", string(lid))
        end
    end

    # Otherwise create a new negative load (MUST include model)
    new_id = "pv_neg_load_$(pv_bus)"

    eng["load"][new_id] = Dict{String,Any}(
        "bus" => pv_bus,
        "connections" => pv_phases,
        "phases" => length(pv_phases),
        "pd" => fill(-pv_kw_per_phase, length(pv_phases)),
        "qd" => fill(0.0, length(pv_phases)),
        "model" => "constant_power",
        "status" => 1
    )

    return ("created_new_load", new_id)
end

# pick PV bus: farthest non-source bus (using dist)
function pick_pv_bus_auto(eng::Dict{String,Any}, dist::Dict{String,Float64})
    bad = Set(["sourcebus", "sourcebusz", "SourceBus", "SourceBusZ"])
    best_bus = ""
    best_d = -Inf
    for (b, d) in dist
        if lowercase(b) in Set(lowercase.(collect(bad)))
            continue
        end
        if d > best_d
            best_d = d
            best_bus = b
        end
    end
    if best_bus != ""
        return best_bus
    end
    # fallback: first non-source bus key
    for b in keys(eng["bus"])
        if lowercase(string(b)) in Set(lowercase.(collect(bad)))
            continue
        end
        return string(b)
    end
    return first(keys(eng["bus"])) |> string
end

## --------------------------------------------------
# 5) Load baseline CSV and select timesteps
# --------------------------------------------------

isfile(BASELINE_CSV) || error("Baseline CSV not found: $(BASELINE_CSV)")

df0 = CSV.read(BASELINE_CSV, DataFrame)
df0 = df0[df0.pf_status .== "LOCALLY_SOLVED", :]

df0.pv_pu = pv_shape_simple(df0.time)

# Stress score: high PV proxy / low load alpha
df0.score = df0.pv_pu ./ (df0.alpha_t .+ 1e-6)

df_sel = sort(df0, :score, rev=true)[1:min(N_PV_STRESS, nrow(df0)), :]
df_sel = unique(df_sel, :timestep)
df_sel = sort(df_sel, :time)

println("\nSelected PV stress timesteps: ", nrow(df_sel))
println("Baseline CSV: ", BASELINE_CSV)
println("Score range: min=", minimum(df_sel.score), " max=", maximum(df_sel.score))

# --------------------------------------------------
# 6) Parse feeder once and prep plot objects
# --------------------------------------------------

println("\nParsing feeder: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!, reduce_lines!])
math0 = PMD.transform_data_model(eng0; kron_reduce=false, phase_project=false)
PMD.add_start_vrvi!(math0)
pf0 = PMD.solve_mc_opf(math0, PMD.IVRENPowerModel, ipopt)
m0 = pf_metrics(pf0; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)


scale_loads!(eng0, LOAD_ALPHA_BASE)
print_baseline_load_summary(eng0)

lines_df = make_lines_df_from_eng(eng0)
source_bus = pick_source_bus(eng0)
dist = compute_bus_distances(lines_df; source_bus=source_bus)
println("\nDistance reference bus: ", source_bus)

pv_bus = (PV_BUS == "AUTO") ? pick_pv_bus_auto(eng0, dist) : PV_BUS
println("\nChosen PV bus: ", pv_bus)

## --------------------------------------------------
# 7) Loop selected timesteps: baseline vs PV
# --------------------------------------------------

rows = NamedTuple[]

for (rank, r) in enumerate(eachrow(df_sel[1:1,:]))
    k = Int(r.timestep)
    a = Float64(r.alpha_t)

    # ---- baseline at timestep k ----
    eng_base = deepcopy(eng0)
    scale_loads!(eng_base, a)
    math_base = PMD.transform_data_model(eng_base; kron_reduce=false, phase_project=false)
    PMD.add_start_vrvi!(math_base)
    pf_base = PMD.solve_mc_opf(math_base, PMD.IVRENPowerModel, ipopt)
    # pf_base, math_base = solve_pf(eng_base)
    m_base = pf_metrics(pf_base; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ---- PV case at timestep k ----
    eng_pv = deepcopy(eng0)
    scale_loads!(eng_pv, a)
    pv_action, pv_load_id = add_pv_negative_load!(eng_pv, pv_bus, PV_KW_PER_PHASE, PV_PHASES)
    math_pv = PMD.transform_data_model(eng_pv; kron_reduce=false, phase_project=false)
    PMD.add_start_vrvi!(math_pv)
    pf_pv = PMD.solve_mc_opf(math_pv, PMD.IVRENPowerModel, ipopt)
    # pf_pv, math_pv = solve_pf(eng_pv)
    m_pv = pf_metrics(pf_pv; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ---- STATCOM case at timestep k ----
    eng_stc = deepcopy(eng0)
    scale_loads!(eng_stc, a)
    pv_action, pv_load_id = add_pv_negative_load!(eng_stc, pv_bus, PV_KW_PER_PHASE, PV_PHASES)
    math_stc = PMD.transform_data_model(eng_stc; kron_reduce=false, phase_project=false)
    ### add inverter lossy branches
    stc_ids = [string(i) for i=1:10] # what is stc id? [i for (i, gen) in data_math["gen"] if !occursin("source", gen["name"])]
    for gen_id in stc_ids
        add_inverter_losses!(math_stc, gen_id; three_wire=false, c_rating_a=30*ones(3))
    end
    for (i, bus) in math_stc["bus"]
        bus["vm_vuf_max"] = 0.02
    end
    
    PMD.add_start_vrvi!(math_stc)
    model_stc = PMD.instantiate_mc_model(math_stc, PMD.IVRENPowerModel, build_mc_opf_mx_Rhea)
    # data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 30] / Ibase   # here you can change STATCOM ratings
    pf_stc = PMD.optimize_model!(model_stc, optimizer=ipopt)

    m_stc = pf_metrics(pf_stc; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)


    push!(rows, (
        net = NET,
        year = YEAR,
        rank = rank,
        timestep = k,
        time = r.time,
        alpha_t = a,
        pv_pu_proxy = Float64(r.pv_pu),
        score = Float64(r.score),
        pv_bus = pv_bus,
        pv_kw_per_phase = PV_KW_PER_PHASE,
        pv_phases = join(PV_PHASES, ","),
        pv_action = pv_action,
        pv_load_id = pv_load_id,
        base_status = m_base.status,
        base_vmin_pu = m_base.vmin,
        base_vmax_pu = m_base.vmax,
        pv_status = m_pv.status,
        pv_vmin_pu = m_pv.vmin,
        pv_vmax_pu = m_pv.vmax,
        delta_vmin_pu = m_pv.vmin - m_base.vmin,
        delta_vmax_pu = m_pv.vmax - m_base.vmax
    ))

    # ---- plots (same mapping logic as 01) ----
    buses_base_by_name, buses_base_by_id = solved_bus_vm_volts_dual(pf_base, math_base; vbase_ln=VBASE_LN)
    buses_pv_by_name,   buses_pv_by_id   = solved_bus_vm_volts_dual(pf_pv,   math_pv;   vbase_ln=VBASE_LN)

    hits_id_base   = count(b -> haskey(buses_base_by_id, b), keys(dist))
    hits_name_base = count(b -> haskey(buses_base_by_name, b), keys(dist))
    buses_base = hits_id_base >= hits_name_base ? buses_base_by_id : buses_base_by_name

    hits_id_pv   = count(b -> haskey(buses_pv_by_id, b), keys(dist))
    hits_name_pv = count(b -> haskey(buses_pv_by_name, b), keys(dist))
    buses_pv = hits_id_pv >= hits_name_pv ? buses_pv_by_id : buses_pv_by_name

    for (bus, d) in dist
        haskey(buses_base, bus) && (buses_base[bus]["distance"] = d)
        haskey(buses_pv,   bus) && (buses_pv[bus]["distance"] = d)
    end

    vminV = VMIN_PU * VBASE_LN
    vmaxV = VMAX_PU * VBASE_LN

    p_base = plot_voltage_combined_snap(buses_base, lines_df; t=1, Vthreshold=1000.0, vminV=vminV, vmaxV=vmaxV)
    savefig(p_base, joinpath(FIGDIR, "rank=$(rank)_timestep=$(k)_baseline.png"))

    p_pv = plot_voltage_combined_snap(buses_pv, lines_df; t=1, Vthreshold=1000.0, vminV=vminV, vmaxV=vmaxV)
    savefig(p_pv, joinpath(FIGDIR, "rank=$(rank)_timestep=$(k)_pv.png"))

    println("Saved rank=", rank, " t=", k,
        " alpha_t=", a,
        " base_vmin=", m_base.vmin,
        " pv_vmin=", m_pv.vmin,
        " delta_vmin=", (m_pv.vmin - m_base.vmin)
    )
end

df_out = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps.csv"), df_out)

# Summary plots
p1 = scatter(df_out.time, df_out.delta_vmin_pu;
    xlabel="Time",
    ylabel="PV minus baseline vmin (pu)",
    title="PV impact on vmin | selected timesteps | $(NET)",
    legend=false
)
savefig(p1, joinpath(FIGDIR, "delta_vmin_selected.png"))

p2 = scatter(df_out.time, df_out.delta_vmax_pu;
    xlabel="Time",
    ylabel="PV minus baseline vmax (pu)",
    title="PV impact on vmax | selected timesteps | $(NET)",
    legend=false
)
savefig(p2, joinpath(FIGDIR, "delta_vmax_selected.png"))

println("\nSaved results to: ", OUTDIR)
println("Done.")
