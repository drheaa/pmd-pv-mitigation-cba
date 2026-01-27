# Compare baseline vs PV vs STATCOM using time-series scaling,
# but only on selected "PV-stress" timesteps.

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
const IM  = InfrastructureModels
PMD.silence!()

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

# --------------------------------------------------
# 0) Settings
# --------------------------------------------------

ROOT = joinpath(@__DIR__, "../..")

# This file is assumed to define:
# - reduce_lines!
# - add_inverter_losses!
# - build_mc_opf_mx_Rhea
include(joinpath(ROOT, "src/read_functions.jl"))

NET_3W = "spd_s"
NET_4W = "spd_s_4w"
master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_scaled.dss")

YEAR   = 2023
STRIDE = 4

BASELINE_OUTDIR = joinpath(ROOT, "results", "time_series", "baseline_pf", NET_3W, "year=$(YEAR)_stride=$(STRIDE)")
BASELINE_CSV    = joinpath(BASELINE_OUTDIR, "tables", "timeseries_baseline_pf_metrics.csv")

# PV settings
PV_BUS          = "AUTO"     # "AUTO" picks farthest non-source bus
PV_KW_PER_PHASE = 10.0       # requested (visibility)
PV_PHASES       = [1]        # typical single-phase PV

# How many timesteps to re-simulate
N_PV_STRESS  = 30
RANDOM_SEED  = 7
MAX_LOOP     = 100           # safety cap (keep)
SKIP_NIGHT_PV = true         # skip pv_kw_eff ~ 0 to avoid boring solves

# Limits
VMIN_PU   = 0.90
VMAX_PU   = 1.10
VBASE_LN  = 230.0

# Base scaling before timestep-specific alpha
LOAD_ALPHA_BASE = 2.5

OUTDIR = joinpath(
    ROOT, "results", "time_series", "pv_pf",
    NET_4W,
    "year=$(YEAR)_stride=$(STRIDE)",
    "pvkw=$(PV_KW_PER_PHASE)_ph=$(join(PV_PHASES, '_'))_K=$(N_PV_STRESS)"
)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

# --------------------------------------------------
# 1) Helpers: scaling, PF metrics
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

function pf_metrics(pf::Dict{String,Any}; vmin_pu=0.9, vmax_pu=1.10)
    sol_bus = pf["solution"]["bus"]
    vmins = Float64[]
    vmaxs = Float64[]

    for sb in values(sol_bus)
        vm =
            haskey(sb, "vm") ? sb["vm"][1:min(3, length(sb["vm"]))] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"][1:min(3, length(sb["vr"]))].^2 .+ sb["vi"][1:min(3, length(sb["vi"]))].^2) :
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
# 2) Plot helpers: feeder distance voltage snapshots
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

function plot_voltage_along_feeder_snap(
    buses_dict::Dict{String,Dict{String,Any}},
    lines_df::DataFrame;
    t::Int=1, Vthreshold::Float64=1000.0,
    vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0
)
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

function plot_voltage_histogram_snap(
    buses_dict::Dict{String,Dict{String,Any}};
    t::Int=1, Vthreshold::Float64=1000.0,
    vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0
)
    phase_a = Float64[]; phase_b = Float64[]; phase_c = Float64[]
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

function plot_voltage_combined_snap(
    buses_dict::Dict{String,Dict{String,Any}},
    lines_df::DataFrame;
    t::Int=1, Vthreshold::Float64=1000.0,
    vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0
)
    p1 = plot_voltage_along_feeder_snap(buses_dict, lines_df; t=t, Vthreshold=Vthreshold, vminV=vminV, vmaxV=vmaxV)
    p2 = plot_voltage_histogram_snap(buses_dict; t=t, Vthreshold=Vthreshold, vminV=vminV, vmaxV=vmaxV)
    return plot(p1, p2, layout=(1,2))
end

# --------------------------------------------------
# 3) PV-stress selection
# --------------------------------------------------

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

function add_pv_negative_load!(eng::Dict{String,Any}, pv_bus::String, pv_kw_per_phase::Float64, pv_phases::Vector{Int})
    haskey(eng, "load") || (eng["load"] = Dict{String,Any}())

    for (lid, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        haskey(ld, "bus") || continue
        string(ld["bus"]) == pv_bus || continue

        pd_key = haskey(ld, "pd") ? "pd" : (haskey(ld, "pd_nom") ? "pd_nom" : nothing)
        pd_key === nothing && continue

        pd_val = ld[pd_key]

        if pd_val isa Number
            ld[pd_key] = Float64(pd_val) - pv_kw_per_phase * length(pv_phases)
            return ("modified_existing_load", string(lid))
        end

        if pd_val isa AbstractVector
            pd_new = Float64.(pd_val)
            for ph in pv_phases
                if 1 <= ph <= length(pd_new)
                    pd_new[ph] -= pv_kw_per_phase
                end
            end
            ld[pd_key] = pd_new
            return ("modified_existing_load", string(lid))
        end
    end

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

function pick_pv_bus_auto(eng::Dict{String,Any}, dist::Dict{String,Float64})
    bad = Set(["sourcebus", "sourcebusz", "SourceBus", "SourceBusZ"])
    best_bus = ""; best_d = -Inf
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
    for b in keys(eng["bus"])
        if lowercase(string(b)) in Set(lowercase.(collect(bad)))
            continue
        end
        return string(b)
    end
    return first(keys(eng["bus"])) |> string
end

# --------------------------------------------------
# 5) Feeder-head current extraction (source branch + pu->A)
# --------------------------------------------------

function _phasor_from(obj::Dict{String,Any}, rkey::String, ikey::String, idx::Int)
    (haskey(obj, rkey) && haskey(obj, ikey)) || return nothing
    r = obj[rkey]; im = obj[ikey]
    (r isa AbstractVector && im isa AbstractVector) || return nothing
    (length(r) >= idx && length(im) >= idx) || return nothing
    return complex(Float64(r[idx]), Float64(im[idx]))
end

function _seq_components(Aa::Complex, Ab::Complex, Ac::Complex)
    a = cis(2pi / 3)
    A0 = (Aa + Ab + Ac) / 3
    A1 = (Aa + a * Ab + a^2 * Ac) / 3
    A2 = (Aa + a^2 * Ab + a * Ac) / 3
    return (A0=A0, A1=A1, A2=A2)
end

function current_base_A(math::Dict{String,Any}; vbase_ln::Float64)
    # PMD typically stores sbase_default in MVA
    sbase_mva = try
        math["settings"]["sbase_default"]
    catch
        1.0
    end
    S = Float64(sbase_mva) * 1e6   # VA (3φ base)
    return S / (3.0 * vbase_ln)    # A per phase (S = 3*Vln*I)
end

function pick_source_branch_id(pf::Dict{String,Any}, math::Dict{String,Any}, source_bus_name::String)
    sol = pf["solution"]
    haskey(sol, "branch") || return nothing
    isempty(sol["branch"]) && return nothing

    # source bus id in math
    src_id = nothing
    for (bid, b_any) in math["bus"]
        b = b_any::Dict{String,Any}
        if string(b["name"]) == source_bus_name
            src_id = bid
            break
        end
    end
    src_id === nothing && return first(keys(sol["branch"]))

    if haskey(math, "branch")
        for (br_id, br_any) in math["branch"]
            br = br_any::Dict{String,Any}
            if br["f_bus"] == src_id || br["t_bus"] == src_id
                if haskey(sol["branch"], br_id)
                    return br_id
                elseif haskey(sol["branch"], string(br_id))
                    return string(br_id)
                end
            end
        end
    end

    return first(keys(sol["branch"]))
end

function feeder_head_current_metrics(
    pf::Dict{String,Any},
    math::Dict{String,Any},
    source_bus_name::String;
    vbase_ln::Float64=230.0
)
    sol = pf["solution"]
    if !haskey(sol, "branch") || isempty(sol["branch"])
        return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)
    end

    bid = pick_source_branch_id(pf, math, source_bus_name)
    bid === nothing && return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)

    br = sol["branch"][bid]
    keypairs = [
        ("cr_fr", "ci_fr"),
        ("cfr",   "cfi"),
        ("cr",    "ci"),
        ("ctr",   "cti"),
    ]

    IbaseA = current_base_A(math; vbase_ln=vbase_ln)

    for (rk, ik) in keypairs
        Ia = _phasor_from(br, rk, ik, 1)
        Ib = _phasor_from(br, rk, ik, 2)
        Ic = _phasor_from(br, rk, ik, 3)
        In = _phasor_from(br, rk, ik, 4)   # neutral if present

        if Ia !== nothing && Ib !== nothing && Ic !== nothing
            ma = abs(Ia); mb = abs(Ib); mc = abs(Ic)
            mn = (In === nothing) ? NaN : abs(In)

            # Heuristic: likely per-unit if small
            if maximum([ma, mb, mc]) < 10.0
                Ia *= IbaseA; Ib *= IbaseA; Ic *= IbaseA
                ma *= IbaseA; mb *= IbaseA; mc *= IbaseA
                if isfinite(mn); mn *= IbaseA; end
                if In !== nothing; In *= IbaseA; end
            end

            seq = _seq_components(Ia, Ib, Ic)
            I0 = abs(seq.A0)
            I1 = abs(seq.A1)
            I2 = abs(seq.A2)
            r  = I1 > 1e-12 ? I2 / I1 : NaN

            return (Ia=ma, Ib=mb, Ic=mc, In=mn, Imean=(ma+mb+mc)/3, Imax=max(ma,mb,mc), I0=I0, I2=I2, I2_over_I1=r)
        end
    end

    return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)
end

function duration_curve(values::AbstractVector{<:Real}; high_is_worse::Bool=true)
    v = Float64.(values)
    v = filter(isfinite, v)
    isempty(v) && return (p=Float64[], x=Float64[])
    x = sort(v; rev=high_is_worse)
    p = collect(range(0, 100; length=length(x)))
    return (p=p, x=x)
end

# --------------------------------------------------
# 6) Load baseline CSV and select PV-stress timesteps
# --------------------------------------------------

isfile(BASELINE_CSV) || error("Baseline CSV not found: $(BASELINE_CSV)")

df0 = CSV.read(BASELINE_CSV, DataFrame)
df0 = df0[df0.pf_status .== "LOCALLY_SOLVED", :]

df0.pv_pu = pv_shape_simple(df0.time)
df0.score = df0.pv_pu ./ (df0.alpha_t .+ 1e-6)

df_sel = sort(df0, :score, rev=true)
df_sel = unique(df_sel, :timestep)
df_sel = sort(df_sel, :score, rev=true)
df_sel = df_sel[1:min(N_PV_STRESS, nrow(df_sel)), :]
df_sel = sort(df_sel, :time)

println("\nSelected PV-stress timesteps: ", nrow(df_sel))
println("Baseline CSV: ", BASELINE_CSV)
println("Score range: min=", minimum(df_sel.score), " max=", maximum(df_sel.score))

# --------------------------------------------------
# 7) Parse feeder once; prep distances & choose PV bus
# --------------------------------------------------

println("\nParsing feeder: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!, reduce_lines!])

scale_loads!(eng0, LOAD_ALPHA_BASE)
print_baseline_load_summary(eng0)

lines_df   = make_lines_df_from_eng(eng0)
source_bus = pick_source_bus(eng0)
dist       = compute_bus_distances(lines_df; source_bus=source_bus)
println("\nDistance reference bus: ", source_bus)

pv_bus = (PV_BUS == "AUTO") ? pick_pv_bus_auto(eng0, dist) : PV_BUS
println("Chosen PV bus: ", pv_bus)

# --------------------------------------------------
# 8) One-step stress test (before full loop)
# --------------------------------------------------

r_test = first(eachrow(sort(df_sel, :score, rev=true)))
k_test = Int(r_test.timestep)
a_test = Float64(r_test.alpha_t)
pv_test = PV_KW_PER_PHASE * Float64(r_test.pv_pu)

println("\nSTRESS TEST")
println("  k=", k_test, " time=", r_test.time, " score=", r_test.score, " alpha=", a_test,
        " pv_pu=", r_test.pv_pu, " pv_kw_eff=", pv_test)

# --------------------------------------------------
# 9) Main loop over df_sel (ranked rows)
# --------------------------------------------------

t_vec = DateTime[]

# Per-phase currents (A)
base_Ia = Float64[]; base_Ib = Float64[]; base_Ic = Float64[]; base_In = Float64[]
pv_Ia   = Float64[]; pv_Ib   = Float64[]; pv_Ic   = Float64[]; pv_In   = Float64[]
stc_Ia  = Float64[]; stc_Ib  = Float64[]; stc_Ic  = Float64[]; stc_In  = Float64[]

base_Imean = Float64[]; base_Imax = Float64[]
pv_Imean   = Float64[]; pv_Imax   = Float64[]
stc_Imean  = Float64[]; stc_Imax  = Float64[]

base_I0 = Float64[]; base_I2 = Float64[]; base_I2r = Float64[]
pv_I0   = Float64[]; pv_I2   = Float64[]; pv_I2r   = Float64[]
stc_I0  = Float64[]; stc_I2  = Float64[]; stc_I2r  = Float64[]

# For Fig-7 style bands
base_Imin = Float64[]; pv_Imin = Float64[]; stc_Imin = Float64[]
base_Ipmax = Float64[]; pv_Ipmax = Float64[]; stc_Ipmax = Float64[]

rows = NamedTuple[]
df_loop = df_sel[1:min(MAX_LOOP, nrow(df_sel)), :]

for (rank, r) in enumerate(eachrow(df_loop))
    k = Int(r.timestep)
    a = Float64(r.alpha_t)
    pv_pu = Float64(r.pv_pu)
    pv_kw_eff = PV_KW_PER_PHASE * pv_pu

    if SKIP_NIGHT_PV && pv_kw_eff < 1e-6
        continue
    end

    # ---- baseline ----
    eng_base = deepcopy(eng0)
    scale_loads!(eng_base, a)
    math_base = PMD.transform_data_model(eng_base; kron_reduce=false, phase_project=false)
    PMD.add_start_vrvi!(math_base)
    pf_base = PMD.solve_mc_opf(math_base, PMD.IVRENPowerModel, ipopt)
    m_base = pf_metrics(pf_base; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ---- PV ----
    eng_pv = deepcopy(eng0)
    scale_loads!(eng_pv, a)
    pv_action, pv_load_id = add_pv_negative_load!(eng_pv, pv_bus, pv_kw_eff, PV_PHASES)
    math_pv = PMD.transform_data_model(eng_pv; kron_reduce=false, phase_project=false)
    PMD.add_start_vrvi!(math_pv)
    pf_pv = PMD.solve_mc_opf(math_pv, PMD.IVRENPowerModel, ipopt)
    m_pv = pf_metrics(pf_pv; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ---- STATCOM ----
    eng_stc = deepcopy(eng0)
    scale_loads!(eng_stc, a)
    pv_action2, pv_load_id2 = add_pv_negative_load!(eng_stc, pv_bus, pv_kw_eff, PV_PHASES)
    math_stc = PMD.transform_data_model(eng_stc; kron_reduce=false, phase_project=false)

    # Add inverter losses / constraints (your existing approach)
    stc_ids = [string(i) for i=1:10]
    for gen_id in stc_ids
        add_inverter_losses!(math_stc, gen_id; three_wire=false, c_rating_a=30*ones(3))
    end
    for (_, bus_any) in math_stc["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vm_vuf_max"] = 0.02
    end

    PMD.add_start_vrvi!(math_stc)
    model_stc = PMD.instantiate_mc_model(math_stc, PMD.IVRENPowerModel, build_mc_opf_mx_Rhea)
    # data_math_conv["gen"]["1"]["c_rating"] = [30 ; 30 ; 30 ; 30] / Ibase   # here you can change STATCOM ratings
    pf_stc = PMD.optimize_model!(model_stc, optimizer=ipopt)
    m_stc = pf_metrics(pf_stc; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ---- verify the branch keys you’re extracting ----
    if rank == 1
        bid = pick_source_branch_id(pf_base, math_base, source_bus)
        println("\nDEBUG feeder-head branch id = ", bid)
        if bid !== nothing
            println("DEBUG branch solution keys: ", collect(keys(pf_base["solution"]["branch"][bid])))
        end
        println("DEBUG Ibase(A) = ", current_base_A(math_base; vbase_ln=VBASE_LN))
    end

    # ---- feeder-head currents (A) ----
    Ibase = feeder_head_current_metrics(pf_base, math_base, source_bus; vbase_ln=VBASE_LN)
    Ipv   = feeder_head_current_metrics(pf_pv,   math_pv,   source_bus; vbase_ln=VBASE_LN)
    Istc  = feeder_head_current_metrics(pf_stc,  math_stc,  source_bus; vbase_ln=VBASE_LN)

    # Store time axis
    push!(t_vec, r.time)

    # Per-phase currents
    push!(base_Ia, Ibase.Ia); push!(base_Ib, Ibase.Ib); push!(base_Ic, Ibase.Ic); push!(base_In, Ibase.In)
    push!(pv_Ia,   Ipv.Ia);   push!(pv_Ib,   Ipv.Ib);   push!(pv_Ic,   Ipv.Ic);   push!(pv_In,   Ipv.In)
    push!(stc_Ia,  Istc.Ia);  push!(stc_Ib,  Istc.Ib);  push!(stc_Ic,  Istc.Ic);  push!(stc_In,  Istc.In)

    # Mean/Max
    push!(base_Imean, Ibase.Imean); push!(base_Imax, Ibase.Imax)
    push!(pv_Imean,   Ipv.Imean);   push!(pv_Imax,   Ipv.Imax)
    push!(stc_Imean,  Istc.Imean);  push!(stc_Imax,  Istc.Imax)

    # Min/Max across phases (for range-bands)
    base_min = minimum([Ibase.Ia, Ibase.Ib, Ibase.Ic])
    base_max = maximum([Ibase.Ia, Ibase.Ib, Ibase.Ic])
    pv_min   = minimum([Ipv.Ia,   Ipv.Ib,   Ipv.Ic])
    pv_max   = maximum([Ipv.Ia,   Ipv.Ib,   Ipv.Ic])
    stc_min  = minimum([Istc.Ia,  Istc.Ib,  Istc.Ic])
    stc_max  = maximum([Istc.Ia,  Istc.Ib,  Istc.Ic])

    push!(base_Imin, base_min); push!(base_Ipmax, base_max)
    push!(pv_Imin,   pv_min);   push!(pv_Ipmax,   pv_max)
    push!(stc_Imin,  stc_min);  push!(stc_Ipmax,  stc_max)

    # Sequence indicators
    push!(base_I0, Ibase.I0); push!(base_I2, Ibase.I2); push!(base_I2r, Ibase.I2_over_I1)
    push!(pv_I0,   Ipv.I0);   push!(pv_I2,   Ipv.I2);   push!(pv_I2r,   Ipv.I2_over_I1)
    push!(stc_I0,  Istc.I0);  push!(stc_I2,  Istc.I2);  push!(stc_I2r,  Istc.I2_over_I1)

    push!(rows, (
        net = NET_4W,
        year = YEAR,
        rank = rank,
        timestep = k,
        time = r.time,
        alpha_t = a,
        pv_pu_proxy = pv_pu,
        score = Float64(r.score),
        pv_bus = pv_bus,
        pv_kw_per_phase = PV_KW_PER_PHASE,
        pv_kw_eff = pv_kw_eff,
        pv_phases = join(PV_PHASES, ","),
        pv_action = pv_action,
        pv_load_id = pv_load_id,
        base_status = m_base.status,
        base_vmin_pu = m_base.vmin,
        base_vmax_pu = m_base.vmax,
        pv_status = m_pv.status,
        pv_vmin_pu = m_pv.vmin,
        pv_vmax_pu = m_pv.vmax,
        stc_status = m_stc.status,
        stc_vmin_pu = m_stc.vmin,
        stc_vmax_pu = m_stc.vmax,
        delta_pv_vmin_pu = m_pv.vmin - m_base.vmin,
        delta_pv_vmax_pu = m_pv.vmax - m_base.vmax,
        delta_stc_vmin_pu = m_stc.vmin - m_base.vmin,
        delta_stc_vmax_pu = m_stc.vmax - m_base.vmax
    ))

    # ---- snapshot voltage panels ----
    vminV = VMIN_PU * VBASE_LN
    vmaxV = VMAX_PU * VBASE_LN

    buses_base_by_name, buses_base_by_id = solved_bus_vm_volts_dual(pf_base, math_base; vbase_ln=VBASE_LN)
    buses_pv_by_name,   buses_pv_by_id   = solved_bus_vm_volts_dual(pf_pv,   math_pv;   vbase_ln=VBASE_LN)
    buses_stc_by_name,  buses_stc_by_id  = solved_bus_vm_volts_dual(pf_stc,  math_stc;  vbase_ln=VBASE_LN)

    hits_id_base   = count(b -> haskey(buses_base_by_id, b), keys(dist))
    hits_name_base = count(b -> haskey(buses_base_by_name, b), keys(dist))
    buses_base = hits_id_base >= hits_name_base ? buses_base_by_id : buses_base_by_name

    hits_id_pv   = count(b -> haskey(buses_pv_by_id, b), keys(dist))
    hits_name_pv = count(b -> haskey(buses_pv_by_name, b), keys(dist))
    buses_pv = hits_id_pv >= hits_name_pv ? buses_pv_by_id : buses_pv_by_name

    hits_id_stc   = count(b -> haskey(buses_stc_by_id, b), keys(dist))
    hits_name_stc = count(b -> haskey(buses_stc_by_name, b), keys(dist))
    buses_stc = hits_id_stc >= hits_name_stc ? buses_stc_by_id : buses_stc_by_name

    for (bus, d) in dist
        haskey(buses_base, bus) && (buses_base[bus]["distance"] = d)
        haskey(buses_pv,   bus) && (buses_pv[bus]["distance"]   = d)
        haskey(buses_stc,  bus) && (buses_stc[bus]["distance"]  = d)
    end

    p_base = plot_voltage_combined_snap(buses_base, lines_df; t=1, Vthreshold=1000.0, vminV=vminV, vmaxV=vmaxV)
    savefig(p_base, joinpath(FIGDIR, "rank=$(rank)_timestep=$(k)_baseline.png"))

    p_pv = plot_voltage_combined_snap(buses_pv, lines_df; t=1, Vthreshold=1000.0, vminV=vminV, vmaxV=vmaxV)
    savefig(p_pv, joinpath(FIGDIR, "rank=$(rank)_timestep=$(k)_pv.png"))

    p_stc = plot_voltage_combined_snap(buses_stc, lines_df; t=1, Vthreshold=1000.0, vminV=vminV, vmaxV=vmaxV)
    savefig(p_stc, joinpath(FIGDIR, "rank=$(rank)_timestep=$(k)_statcom.png"))

    p_panel = plot(p_base, p_pv, p_stc, layout=(1,3))
    savefig(p_panel, joinpath(FIGDIR, "rank=$(rank)_timestep=$(k)_panel_base_pv_statcom.png"))

    println("Saved rank=", rank,
        " k=", k,
        " alpha=", a,
        " pv_kw_eff=", round(pv_kw_eff, digits=3),
        " base_vmin=", round(m_base.vmin, digits=6),
        " pv_vmin=", round(m_pv.vmin, digits=6),
        " stc_vmin=", round(m_stc.vmin, digits=6)
    )
end

# --------------------------------------------------
# 10) Save tables
# --------------------------------------------------

df_out = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps.csv"), df_out)


curr_df = DataFrame(
    time=t_vec,
    base_Ia=base_Ia, base_Ib=base_Ib, base_Ic=base_Ic, base_In=base_In,
    base_Imean=base_Imean, base_Imax=base_Imax, base_Imin=base_Imin, base_Ipmax=base_Ipmax,
    base_I0=base_I0, base_I2=base_I2, base_I2r=base_I2r,

    pv_Ia=pv_Ia, pv_Ib=pv_Ib, pv_Ic=pv_Ic, pv_In=pv_In,
    pv_Imean=pv_Imean, pv_Imax=pv_Imax, pv_Imin=pv_Imin, pv_Ipmax=pv_Ipmax,
    pv_I0=pv_I0, pv_I2=pv_I2, pv_I2r=pv_I2r,

    stc_Ia=stc_Ia, stc_Ib=stc_Ib, stc_Ic=stc_Ic, stc_In=stc_In,
    stc_Imean=stc_Imean, stc_Imax=stc_Imax, stc_Imin=stc_Imin, stc_Ipmax=stc_Ipmax,
    stc_I0=stc_I0, stc_I2=stc_I2, stc_I2r=stc_I2r
)
CSV.write(joinpath(TBLDIR, "substation_current_comparison_selected.csv"), curr_df)

# --------------------------------------------------
# 11) Summary “impact” plots (PV and STATCOM vs baseline)
# --------------------------------------------------

p1 = scatter(df_out.time, df_out.delta_pv_vmin_pu;
    xlabel="Time",
    ylabel="PV minus baseline vmin (pu)",
    title="PV impact on vmin | selected timesteps | $(NET_4W)",
    legend=false
)
savefig(p1, joinpath(FIGDIR, "delta_vmin_selected_pv.png"))

p2 = scatter(df_out.time, df_out.delta_pv_vmax_pu;
    xlabel="Time",
    ylabel="PV minus baseline vmax (pu)",
    title="PV impact on vmax | selected timesteps | $(NET_4W)",
    legend=false
)
savefig(p2, joinpath(FIGDIR, "delta_vmax_selected_pv.png"))

p3 = scatter(df_out.time, df_out.delta_stc_vmin_pu;
    xlabel="Time",
    ylabel="STATCOM minus baseline vmin (pu)",
    title="STATCOM impact on vmin | selected timesteps | $(NET_4W)",
    legend=false
)
savefig(p3, joinpath(FIGDIR, "delta_vmin_selected_statcom.png"))

p4 = scatter(df_out.time, df_out.delta_stc_vmax_pu;
    xlabel="Time",
    ylabel="STATCOM minus baseline vmax (pu)",
    title="STATCOM impact on vmax | selected timesteps | $(NET_4W)",
    legend=false
)
savefig(p4, joinpath(FIGDIR, "delta_vmax_selected_statcom.png"))

# --------------------------------------------------
# 12) Paper-style plots (Fig 6–9 vibe)
# --------------------------------------------------

# Fig 6 vibe: phase currents (baseline)
pI_base = plot(t_vec, base_Ia; label="phase a", xlabel="Time", ylabel="Substation current (A)",
    title="Substation phase currents | baseline | PV-stress timesteps")
plot!(t_vec, base_Ib; label="phase b")
plot!(t_vec, base_Ic; label="phase c")
plot!(t_vec, base_Imean; label="mean", linewidth=2)
savefig(pI_base, joinpath(FIGDIR, "fig6_like_currents_baseline.png"))

# PV
pI_pv = plot(t_vec, pv_Ia; label="phase a", xlabel="Time", ylabel="Substation current (A)",
    title="Substation phase currents | PV | PV-stress timesteps")
plot!(t_vec, pv_Ib; label="phase b")
plot!(t_vec, pv_Ic; label="phase c")
plot!(t_vec, pv_Imean; label="mean", linewidth=2)
savefig(pI_pv, joinpath(FIGDIR, "fig6_like_currents_pv.png"))

# STATCOM
pI_stc = plot(t_vec, stc_Ia; label="phase a", xlabel="Time", ylabel="Substation current (A)",
    title="Substation phase currents | STATCOM | PV-stress timesteps")
plot!(t_vec, stc_Ib; label="phase b")
plot!(t_vec, stc_Ic; label="phase c")
plot!(t_vec, stc_Imean; label="mean", linewidth=2)
savefig(pI_stc, joinpath(FIGDIR, "fig6_like_currents_statcom.png"))

# Fig 7 vibe: mean + range band (min/max across phases) baseline vs STATCOM
# ribbon expects distance from y to bounds (upper/lower deltas)
base_ribbon = (base_Imean .- base_Imin, base_Ipmax .- base_Imean)
stc_ribbon  = (stc_Imean  .- stc_Imin,  stc_Ipmax  .- stc_Imean)

p_band = plot(t_vec, base_Imean; label="Mean before", xlabel="Time", ylabel="Substation current (A)",
    title="Feeder-head phase current range | baseline vs STATCOM", ribbon=base_ribbon, fillalpha=0.15)
plot!(t_vec, stc_Imean; label="Mean after", ribbon=stc_ribbon, fillalpha=0.15)
savefig(p_band, joinpath(FIGDIR, "fig7_like_mean_with_range_band_base_vs_statcom.png"))

# Fig 8: duration curve of max phase current
dc_base = duration_curve(base_Imax; high_is_worse=true)
dc_pv   = duration_curve(pv_Imax;   high_is_worse=true)
dc_stc  = duration_curve(stc_Imax;  high_is_worse=true)

p_dc = plot(dc_base.p, dc_base.x; label="baseline", xlabel="Load percentile (%)", ylabel="Max phase current (A)",
    title="Duration curve: max phase current | baseline vs PV vs STATCOM")
plot!(dc_pv.p, dc_pv.x; label="PV")
plot!(dc_stc.p, dc_stc.x; label="STATCOM")
savefig(p_dc, joinpath(FIGDIR, "fig8_like_duration_i_max_phase_base_pv_statcom.png"))

# Fig 9 vibe: phase + neutral (only if neutral is present)
hasN = any(isfinite.(base_In)) || any(isfinite.(pv_In)) || any(isfinite.(stc_In))

if hasN
    pN = plot(t_vec, base_Ia; label="A (base)", xlabel="Time", ylabel="Current (A)",
        title="Phase + neutral currents | baseline vs STATCOM (neutral if available)")
    plot!(t_vec, base_Ib; label="B (base)")
    plot!(t_vec, base_Ic; label="C (base)")
    plot!(t_vec, base_In; label="N (base)", linewidth=2)

    plot!(t_vec, stc_Ia; label="A (stc)", linestyle=:dash)
    plot!(t_vec, stc_Ib; label="B (stc)", linestyle=:dash)
    plot!(t_vec, stc_Ic; label="C (stc)", linestyle=:dash)
    plot!(t_vec, stc_In; label="N (stc)", linestyle=:dash, linewidth=2)

    savefig(pN, joinpath(FIGDIR, "fig9_like_phase_and_neutral_base_vs_statcom.png"))
end

# Fig 11 vibe (simple): jittered scatter for I0 and I2/I1 per case + median lines
function jitter(x::Real, j::Real=0.08)
    return x + (rand() - 0.5) * 2j
end

Random.seed!(RANDOM_SEED)

x_base = [jitter(1.0) for _ in 1:length(base_I2r)]
x_pv   = [jitter(2.0) for _ in 1:length(pv_I2r)]
x_stc  = [jitter(3.0) for _ in 1:length(stc_I2r)]

p_seq = scatter(x_base, base_I2r; label="baseline", xlabel="", ylabel="I2/I1", title="Negative sequence ratio (I2/I1) | distribution across timesteps")
scatter!(x_pv,   pv_I2r;  label="PV")
scatter!(x_stc,  stc_I2r; label="STATCOM")
hline!([median(filter(isfinite, base_I2r))]; label="median base", linestyle=:dash)
hline!([median(filter(isfinite, stc_I2r))];  label="median stc", linestyle=:dash)
xticks!([1,2,3], ["baseline","PV","STATCOM"])
savefig(p_seq, joinpath(FIGDIR, "fig11_like_i2_over_i1_scatter.png"))

p_i0 = scatter(x_base, base_I0; label="baseline", xlabel="", ylabel="I0 (A)", title="Zero-sequence magnitude (I0) | distribution across timesteps")
scatter!(x_pv,   pv_I0;  label="PV")
scatter!(x_stc,  stc_I0; label="STATCOM")
xticks!([1,2,3], ["baseline","PV","STATCOM"])
savefig(p_i0, joinpath(FIGDIR, "fig11_like_i0_scatter.png"))

println("\nSaved results to: ", OUTDIR)
println("Tables:")
println("  - ", joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps.csv"))
println("  - ", joinpath(TBLDIR, "substation_current_comparison_selected.csv"))
println("Done.")
