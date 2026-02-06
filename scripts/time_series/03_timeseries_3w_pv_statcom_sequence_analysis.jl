# Time-series comparison on selected PV-stress timesteps and a
# contiguous 48-hour window for paper-style current plots:
#
# (1) Baseline (NO PV, NO STATCOM)
# (2) PV only (PV as negative load)
# (3) PV + STATCOM (OPF w inverter gens)
#
# Non-negotiables enforced:
# - Baseline and PV never parse master_scaled.dss
# - STATCOM scenario only parses master_scaled.dss
# - STATCOM integration path remains OPF + inverter gens workflow
# ================================================================

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
using JuMP
using LinearAlgebra
using Random

const PMD = PowerModelsDistribution
PMD.silence!()

# --------------------------------------------------
# 0) Settings
# --------------------------------------------------
ROOT = joinpath(@__DIR__, "../..")

# Required helper functions
include(joinpath(ROOT, "src/read_functions.jl"))

NET_4W = "spd_s_4w"

# MASTER_BASELINE = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w.dss")
# MASTER_STATCOM = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_scaled.dss")

MASTER_BASELINE = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w_new.dss")
# MASTER_BASELINE = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w.dss")

MASTER_BASELINE_3w = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", "spd_s", "master_scaled_new.dss")

YEAR = 2023
STRIDE = 4

BASELINE_OUTDIR = joinpath(
    ROOT, "results", "time_series", "baseline_pf", NET_4W,
    "year=$(YEAR)_stride=$(STRIDE)_K=17520"
)
BASELINE_CSV = joinpath(BASELINE_OUTDIR, "tables", "timeseries_baseline_pf_metrics.csv")

# PV settings
PV_BUS = "AUTO"
PV_KW_PER_PHASE = 10.0          # increase if PV impact remains too small
PV_PHASES = [1]                 # single-phase PV on phase A by default

# PV-stress
N_PV_STRESS = 30
SKIP_NIGHT_PV = true
RANDOM_SEED = 42

# 48-hour window for paper-style current plots (Fig6/Fig7/Fig9)
WINDOW_HOURS = 48
WINDOW_STEPS = Int(WINDOW_HOURS * 60 ÷ 30) + 1  # 30-min baseline sampling
WINDOW_CENTER = "BEST_STRESS"                   # "BEST_STRESS" or "FIRST_STRESS"

# Limits and conversion
VMIN_PU = 0.90
VMAX_PU = 1.10
VBASE_LN = 230.0

# Baseline engineering scaling applied once, then alpha(t)
LOAD_ALPHA_BASE = 1.5  # scale baseline loads by this factor once at parsing and should be between 1.5 to 2.0

# Solver
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

OUTDIR = joinpath(
    ROOT, "results", "time_series", "pv_pf", NET_4W,
    "year=$(YEAR)_stride=$(STRIDE)",
    "pvkw=$(PV_KW_PER_PHASE)_ph=$(join(PV_PHASES, '_'))_K=$(N_PV_STRESS)"
)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

## --------------------------------------------------
# 1) Helpers: scaling, metrics, PV proxy, PV injection
# --------------------------------------------------
function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng  # early exit if no loads

    for (_, ld_any) in eng["load"]  # loop through every load object
        ld = ld_any::Dict{String,Any}  # treat each load as a dictionary
        for k in ("pd_nom", "qd_nom", "pd", "qd")  # for each power-related key, scale if it exists
            if haskey(ld, k)  # only scale if the load actually has that field
                v = ld[k]  # get the current value
                ld[k] = v isa Number ? alpha * v : alpha .* v  # scale numbers vs vectors correctly
            end
        end
    end
    return eng  # return the mutated dictionary
end

function add_load_q_from_pf!(eng::Dict{String,Any}; pf::Float64=0.95)
    # add reactive power to loads based on specified power factor
    haskey(eng, "load") || return eng
    pf = clamp(pf, 0.5, 0.999999)  # avoid invalid power factor values
    q_factor = tan(acos(pf))       # Q = |P| * tan(arccos(pf))

    println("DEBUG add_load_q_from_pf!: pf=", pf, " q_factor=", tan(acos(pf)))

    for (_, ld_any) in eng["load"]  # loop through every load object
        ld = ld_any::Dict{String,Any}

        # Choose the active P key (pd preferred, otherwise pd_nom)
        pkey = haskey(ld, "pd") ? "pd" : (haskey(ld, "pd_nom") ? "pd_nom" : nothing)
        pkey === nothing && continue

        # Write matching Q key
        qkey = (pkey == "pd") ? "qd" : "qd_nom"

        P = ld[pkey]
        if P isa Number
            ld[qkey] = abs(Float64(P)) * q_factor
        elseif P isa AbstractVector
            Pvec = Float64.(P)
            ld[qkey] = abs.(Pvec) .* q_factor
        end
    end
    return eng
end

function sum_loads_eng(eng::Dict{String,Any})
    # sum up loads in kW and kVAR from eng dictionary
    P = 0.0
    Q = 0.0
    for ld_any in values(get(eng, "load", Dict{String,Any}()))
        ld = ld_any::Dict{String,Any}

        # prefer *_nom in eng
        if haskey(ld, "pd_nom")
            v = ld["pd_nom"]
            P += v isa Number ? Float64(v) : sum(Float64.(v))
        elseif haskey(ld, "pd")
            v = ld["pd"]
            P += v isa Number ? Float64(v) : sum(Float64.(v))
        end

        if haskey(ld, "qd_nom")
            v = ld["qd_nom"]
            Q += v isa Number ? Float64(v) : sum(Float64.(v))
        elseif haskey(ld, "qd")
            v = ld["qd"]
            Q += v isa Number ? Float64(v) : sum(Float64.(v))
        end
    end
    return (P_kW = P, Q_kvar = Q)
end

function pf_metrics(pf::Dict{String,Any}; vmin_pu=0.90, vmax_pu=1.10)
    sol_bus = pf["solution"]["bus"]  # grab the solved bus results
    vmins = Float64[]
    vmaxs = Float64[]

    for sb in values(sol_bus)
        vm = haskey(sb, "vm") ? sb["vm"][1:min(3, length(sb["vm"]))] :
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

function pv_shape_simple(time_vec::Vector{DateTime})
    # simple normalized PV shape over the day
    pv = zeros(Float64, length(time_vec))
    for (i, ts) in enumerate(time_vec)
        h = hour(ts) + minute(ts)/60
        x = (h - 6.0)/12.0
        pv[i] = (0.0 <= x <= 1.0) ? sin(pi*x)^2 : 0.0
    end
    m = maximum(pv)
    return m > 0 ? pv ./ m : pv
end

function add_pv_negative_load!(eng::Dict{String,Any}, pv_bus::String, pv_kw_per_phase::Float64, pv_phases::Vector{Int})
    haskey(eng, "load") || (eng["load"] = Dict{String,Any}())

    # modifying an existing load at pv_bus
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

    # Fallback: create a new negative load
    load_dict_template = eng["load"][first(keys(eng["load"]))]::Dict{String,Any}
    new_id = "pv_neg_load_$(pv_bus)"
    eng["load"][new_id] = load_dict_template
    eng["load"][new_id]["bus"] = pv_bus
    eng["load"][new_id]["connections"] = pv_phases
    eng["load"][new_id]["phases"] = length(pv_phases)
    eng["load"][new_id]["pd"] = fill(-pv_kw_per_phase, length(pv_phases))
    eng["load"][new_id]["qd"] = fill(0.0, length(pv_phases))
    @show eng["load"][new_id]
    return ("created_new_load", new_id)
end

function kw_to_w!(eng::Dict{String,Any})
    # convert all kW/kVar load ratings to W/Var
    haskey(eng, "load") || return eng

    for ld_any in values(eng["load"])
        ld = ld_any::Dict{String,Any}
        for k in ("pd_nom", "qd_nom", "pd", "qd")
            haskey(ld, k) || continue
            v = ld[k]
            ld[k] = v isa Number ? 1000.0 * Float64(v) : 1000.0 .* Float64.(v)
        end
    end
    return eng
end

# --------------------------------------------------
# 2) Helpers: PV bus AUTO selection (distance-based)
# --------------------------------------------------
function make_lines_df_from_eng(eng::Dict{String,Any})
    # extract lines/branches into a DataFrame
    rows = NamedTuple[]
    if haskey(eng, "line")
        for (_, ln_any) in eng["line"]
            ln = ln_any::Dict{String,Any}
            push!(rows, (Bus1=string(ln["f_bus"]), Bus2=string(ln["t_bus"]),
                         length_km=(get(ln, "length", 0.0) / 1000.0)))
        end
    elseif haskey(eng, "branch")
        for (_, br_any) in eng["branch"]
            br = br_any::Dict{String,Any}
            push!(rows, (Bus1=string(br["f_bus"]), Bus2=string(br["t_bus"]),
                         length_km=(get(br, "length", 0.0) / 1000.0)))
        end
    else
        error("No line or branch data found")
    end
    return DataFrame(rows)
end

function pick_source_bus_name(eng::Dict{String,Any})
    if haskey(eng, "bus")
        names = Set(string.(keys(eng["bus"])))
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
        push!(get!(adj, r.Bus1, Tuple{String,Float64}[]), (r.Bus2, r.length_km))
        push!(get!(adj, r.Bus2, Tuple{String,Float64}[]), (r.Bus1, r.length_km))
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

function pick_pv_bus_auto(eng::Dict{String,Any}, dist::Dict{String,Float64})
    bad = Set(["sourcebus", "sourcebusz", "SourceBus", "SourceBusZ"])
    best_bus = ""
    best_d = -Inf

    for (b, d) in dist
        if lowercase(b) in lowercase.(collect(bad))
            continue
        end
        if d > best_d
            best_d = d
            best_bus = b
        end
    end

    best_bus != "" && return best_bus

    for b in keys(eng["bus"])
        if lowercase(string(b)) in lowercase.(collect(bad))
            continue
        end
        return string(b)
    end

    return first(keys(eng["bus"])) |> string
end

# --------------------------------------------------
# 3) Helpers: robust feeder-head current extraction and metrics
# --------------------------------------------------
function current_base_A(math::Dict{String,Any}; vbase_ln::Float64)
    sbase = get(get(math, "settings", Dict{String,Any}()), "sbase_default", 1.0)
    S = Float64(sbase)
    return S / (3.0 * vbase_ln)
end

function _phasor_from(br::Dict{String,Any}, rkey::String, ikey::String, idx::Int)
    (haskey(br, rkey) && haskey(br, ikey)) || return nothing
    r = br[rkey]; im = br[ikey]
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

function _branch_current_phasors(br::Dict{String,Any})
    keypairs = [
        ("cr_fr", "ci_fr"),
        ("cfr", "cfi"),
        ("cr", "ci"),
        ("ctr", "cti"),
    ]

    for (rk, ik) in keypairs
        Ia = _phasor_from(br, rk, ik, 1)
        Ib = _phasor_from(br, rk, ik, 2)
        Ic = _phasor_from(br, rk, ik, 3)
        In = _phasor_from(br, rk, ik, 4)
        if Ia !== nothing && Ib !== nothing && Ic !== nothing
            return (Ia=Ia, Ib=Ib, Ic=Ic, In=In, rk=rk, ik=ik)
        end
    end
    return nothing
end

function _find_source_bus_id(math::Dict{String,Any}, source_bus_name::String)
    haskey(math, "bus") || return nothing
    for (bid, b_any) in math["bus"]
        b = b_any::Dict{String,Any}
        if haskey(b, "name") && string(b["name"]) == source_bus_name
            return bid
        end
    end
    return nothing
end

function _pick_head_branch_id(pf::Dict{String,Any}, math::Dict{String,Any}, source_bus_name::String)
    sol = pf["solution"]
    haskey(sol, "branch") || return nothing
    isempty(sol["branch"]) && return nothing

    # Priority 1: explicit source impedance line
    if haskey(math, "branch")
        for (br_id, br_any) in math["branch"]
            brm = br_any::Dict{String,Any}
            sid = get(brm, "source_id", "")
            if occursin("line.sourceZ", string(sid)) || occursin("sourceZ", lowercase(string(sid)))
                if haskey(sol["branch"], br_id)
                    return br_id
                elseif haskey(sol["branch"], string(br_id))
                    return string(br_id)
                end
            end
        end
    end

    # Priority 2: branch connected to source bus with largest 3-phase current
    src_id = _find_source_bus_id(math, source_bus_name)
    best = nothing
    best_mag = -Inf

    if src_id !== nothing && haskey(math, "branch")
        for (br_id, br_any) in math["branch"]
            brm = br_any::Dict{String,Any}
            connected = (brm["f_bus"] == src_id || brm["t_bus"] == src_id)
            connected || continue

            sol_id = haskey(sol["branch"], br_id) ? br_id :
                     (haskey(sol["branch"], string(br_id)) ? string(br_id) : nothing)
            sol_id === nothing && continue

            brs = sol["branch"][sol_id]
            ph = _branch_current_phasors(brs)
            ph === nothing && continue

            mag = maximum([abs(ph.Ia), abs(ph.Ib), abs(ph.Ic)])
            if mag > best_mag
                best_mag = mag
                best = sol_id
            end
        end
        best !== nothing && return best
    end

    # Priority 3: any branch with largest 3-phase current
    for (br_id, brs_any) in sol["branch"]
        brs = brs_any::Dict{String,Any}
        ph = _branch_current_phasors(brs)
        ph === nothing && continue
        mag = maximum([abs(ph.Ia), abs(ph.Ib), abs(ph.Ic)])
        if mag > best_mag
            best_mag = mag
            best = br_id
        end
    end

    return best === nothing ? first(keys(sol["branch"])) : best
end

function feeder_head_current_metrics(pf::Dict{String,Any}, math::Dict{String,Any}, source_bus_name::String; vbase_ln::Float64=230.0)
    sol = pf["solution"]
    if !haskey(sol, "branch") || isempty(sol["branch"])
        return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)
    end

    bid = _pick_head_branch_id(pf, math, source_bus_name)
    bid === nothing && return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)

    br = sol["branch"][bid]
    ph = _branch_current_phasors(br)
    ph === nothing && return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)

    Ia, Ib, Ic, In = ph.Ia, ph.Ib, ph.Ic, ph.In
    ma, mb, mc = abs(Ia), abs(Ib), abs(Ic)
    mn = (In === nothing) ? NaN : abs(In)

    if maximum([ma, mb, mc]) < 10.0
        IbaseA = current_base_A(math; vbase_ln=vbase_ln)
        Ia *= IbaseA; Ib *= IbaseA; Ic *= IbaseA
        ma *= IbaseA; mb *= IbaseA; mc *= IbaseA
        if isfinite(mn)
            mn *= IbaseA
        end
        if In !== nothing
            In *= IbaseA
        end
    end

    seq = _seq_components(Ia, Ib, Ic)
    I0 = abs(seq.A0)
    I1 = abs(seq.A1)
    I2 = abs(seq.A2)
    r = I1 > 1e-12 ? I2 / I1 : NaN

    return (Ia=ma, Ib=mb, Ic=mc, In=mn, Imean=(ma+mb+mc)/3, Imax=max(ma,mb,mc), I0=I0, I2=I2, I2_over_I1=r)
end

function sum_loads_math_pu(math::Dict{String,Any})
    Ppu = 0.0; Qpu = 0.0
    for ld_any in values(get(math, "load", Dict{String,Any}()))
        ld = ld_any::Dict{String,Any}
        if haskey(ld, "pd")
            v = ld["pd"]
            Ppu += v isa Number ? Float64(v) : sum(Float64.(v))
        end
        if haskey(ld, "qd")
            v = ld["qd"]
            Qpu += v isa Number ? Float64(v) : sum(Float64.(v))
        end
    end
    return (P_pu=Ppu, Q_pu=Qpu)
end

function sbase_math(math::Dict{String,Any})
    if haskey(math, "settings") && haskey(math["settings"], "sbase_default")
        return Float64(math["settings"]["sbase_default"])
    end
    if haskey(math, "baseMVA")
        return Float64(math["baseMVA"])
    end
    return NaN
end

# --------------------------------------------------
# 4) Helpers: voltage sequence extraction for all buses
# --------------------------------------------------
function bus_voltage_phasors(pf::Dict{String,Any})
    sol_bus = pf["solution"]["bus"]
    out = Dict{String, Tuple{ComplexF64,ComplexF64,ComplexF64}}()

    for (bid, sb_any) in sol_bus
        sb = sb_any::Dict{String,Any}

        if haskey(sb, "vr") && haskey(sb, "vi")
            vr = sb["vr"]; vi = sb["vi"]
            if vr isa AbstractVector && vi isa AbstractVector && length(vr) >= 3 && length(vi) >= 3
                Va = complex(Float64(vr[1]), Float64(vi[1]))
                Vb = complex(Float64(vr[2]), Float64(vi[2]))
                Vc = complex(Float64(vr[3]), Float64(vi[3]))
                out[string(bid)] = (Va, Vb, Vc)
            end
        elseif haskey(sb, "vm") && haskey(sb, "va")
            vm = sb["vm"]; va = sb["va"]
            if vm isa AbstractVector && va isa AbstractVector && length(vm) >= 3 && length(va) >= 3
                Va = Float64(vm[1]) * cis(Float64(va[1]))
                Vb = Float64(vm[2]) * cis(Float64(va[2]))
                Vc = Float64(vm[3]) * cis(Float64(va[3]))
                out[string(bid)] = (Va, Vb, Vc)
            end
        end
    end

    return out
end

function voltage_sequence_rows(pf::Dict{String,Any}, time::DateTime, scenario::String)
    V = bus_voltage_phasors(pf)
    rows = NamedTuple[]
    for (bus_id, (Va,Vb,Vc)) in V
        seq = _seq_components(Va,Vb,Vc)
        V0 = abs(seq.A0)
        V1 = abs(seq.A1)
        V2 = abs(seq.A2)
        r = V1 > 1e-12 ? (V2 / V1) : NaN
        push!(rows, (time=time, scenario=scenario, bus_id=bus_id, V0=V0, V1=V1, V2=V2, V2_over_V1=r))
    end
    return rows
end

function bus_vm_abc(pf::Dict{String,Any})
    sol_bus = pf["solution"]["bus"]

    ks = collect(keys(sol_bus))
    keys_sorted = sort(ks; lt = (a,b) -> begin
        sa = string(a); sb = string(b)
        ia = tryparse(Int, sa)
        ib = tryparse(Int, sb)
        if ia !== nothing && ib !== nothing
            return ia < ib
        else
            return sa < sb
        end
    end)

    vm_a = Float64[]; vm_b = Float64[]; vm_c = Float64[]
    for bid in keys_sorted
        bus = sol_bus[bid]::Dict{String,Any}
        if haskey(bus, "vm")
            vm = bus["vm"]
            push!(vm_a, Float64(vm[1]))
            push!(vm_b, Float64(vm[2]))
            push!(vm_c, Float64(vm[3]))
        else
            vr = bus["vr"]; vi = bus["vi"]
            push!(vm_a, hypot(Float64(vr[1]), Float64(vi[1])))
            push!(vm_b, hypot(Float64(vr[2]), Float64(vi[2])))
            push!(vm_c, hypot(Float64(vr[3]), Float64(vi[3])))
        end
    end
    return (a=vm_a, b=vm_b, c=vm_c, keys=keys_sorted)
end

function branch_current_sequence_rows(pf::Dict{String,Any}, time::DateTime, scenario::String)
    rows = NamedTuple[]
    for (br_id, br_any) in pf["solution"]["branch"]
        br = br_any::Dict{String,Any}
        ph = _branch_current_phasors(br)
        ph === nothing && continue
        Ia, Ib, Ic = ph.Ia, ph.Ib, ph.Ic
        seq = _seq_components(Ia, Ib, Ic)
        I0 = abs(seq.A0)
        I1 = abs(seq.A1)
        I2 = abs(seq.A2)
        r = I1 > 1e-12 ? I2 / I1 : NaN
        push!(rows, (time=time, scenario=scenario, branch_id=string(br_id), I0=I0, I1=I1, I2=I2, I2_over_I1=r))
    end
    return rows
end

# --------------------------------------------------
# 5) Helpers: plotting utilities (paper-style)
# --------------------------------------------------
function duration_curve(values::AbstractVector{<:Real}; high_is_worse::Bool=true)
    v = Float64.(values)
    v = filter(isfinite, v)
    isempty(v) && return (p=Float64[], x=Float64[])
    x = sort(v; rev=high_is_worse)
    p = collect(range(0, 100; length=length(x)))
    return (p=p, x=x)
end

function plot_inputs(df_sel::DataFrame, outpath::String)
    p1 = plot(df_sel.time, df_sel.alpha_t; xlabel="Time", ylabel="alpha(t)", label="alpha(t)", title="Input load scaling alpha(t)")
    p2 = plot(df_sel.time, df_sel.pv_pu; xlabel="Time", ylabel="pv_pu_proxy", label="pv_pu_proxy", title="Input PV proxy pv_pu_proxy(t)")
    p3 = plot(df_sel.time, df_sel.pv_kw_eff; xlabel="Time", ylabel="PV kW per phase", label="pv_kw_eff", title="Effective PV injection per phase (kW)")
    plt = plot(p1, p2, p3; layout=(3,1), size=(900,900))
    savefig(plt, outpath)
end

function fig6_like_currents(df::DataFrame, prefix::String, outpath::String)
    p = plot(df.time, df[!, Symbol(prefix*"_Ia")]; label="phase a", xlabel="Time", ylabel="Current (A)", title="Substation phase currents | $(prefix)")
    plot!(df.time, df[!, Symbol(prefix*"_Ib")]; label="phase b")
    plot!(df.time, df[!, Symbol(prefix*"_Ic")]; label="phase c")
    if hasproperty(df, Symbol(prefix*"_In"))
        plot!(df.time, df[!, Symbol(prefix*"_In")]; label="neutral")
    end
    plot!(df.time, df[!, Symbol(prefix*"_Imean")]; label="mean", linewidth=2)
    savefig(p, outpath)
end

function fig7_like_range_band(df::DataFrame, base_prefix::String, after_prefix::String, outpath::String)
    base_min = min.(df[!, Symbol(base_prefix*"_Ia")], df[!, Symbol(base_prefix*"_Ib")], df[!, Symbol(base_prefix*"_Ic")])
    base_max = max.(df[!, Symbol(base_prefix*"_Ia")], df[!, Symbol(base_prefix*"_Ib")], df[!, Symbol(base_prefix*"_Ic")])
    after_min = min.(df[!, Symbol(after_prefix*"_Ia")], df[!, Symbol(after_prefix*"_Ib")], df[!, Symbol(after_prefix*"_Ic")])
    after_max = max.(df[!, Symbol(after_prefix*"_Ia")], df[!, Symbol(after_prefix*"_Ib")], df[!, Symbol(after_prefix*"_Ic")])

    p = plot(df.time, df[!, Symbol(base_prefix*"_Imean")]; label="Mean before", xlabel="Time", ylabel="Substation current (A)", title="Substation current range band | before vs after")
    plot!(df.time, df[!, Symbol(after_prefix*"_Imean")]; label="Mean after", linewidth=2)

    plot!(df.time, base_min; label="Before (min)", linestyle=:dot)
    plot!(df.time, base_max; label="Before (max)", linestyle=:dot)

    plot!(df.time, after_min; label="After (min)", linestyle=:dash)
    plot!(df.time, after_max; label="After (max)", linestyle=:dash)

    savefig(p, outpath)
end

function fig8_like_duration(df::DataFrame, outpath::String)
    dc_base = duration_curve(df.base_Imax; high_is_worse=true)
    dc_pv   = duration_curve(df.pv_Imax; high_is_worse=true)
    dc_stc  = duration_curve(df.stc_Imax; high_is_worse=true)

    p = plot(dc_base.p, dc_base.x; label="Baseline", xlabel="Load percentile (%)", ylabel="Max phase current (A)", title="Load duration curve of maximum phase current (substation)")
    plot!(dc_pv.p, dc_pv.x; label="PV")
    plot!(dc_stc.p, dc_stc.x; label="PV + STATCOM")

    savefig(p, outpath)
end

function fig11_like_voltage_seq_distributions(df_seq::DataFrame, outdir::String)
    mkpath(outdir)

    function scat_med(ycol::Symbol, title_str::String, outpath::String)
        scenarios = ["baseline", "pv", "statcom"]
        xs = Float64[]
        ys = Float64[]
        for (i, sc) in enumerate(scenarios)
            sub = df_seq[df_seq.scenario .== sc, :]
            y = sub[!, ycol]
            y = y[isfinite.(y)]
            append!(ys, y)
            append!(xs, fill(i, length(y)) .+ (rand(length(y)) .- 0.5) .* 0.15)
        end

        p = scatter(xs, ys; markersize=2, xlabel="", ylabel=String(ycol), title=title_str, legend=false)
        xticks!(p, (1:3, scenarios))

        for (i, sc) in enumerate(scenarios)
            sub = df_seq[df_seq.scenario .== sc, :]
            y = sub[!, ycol]
            y = y[isfinite.(y)]
            isempty(y) && continue
            med = median(y)
            plot!(p, [i-0.25, i+0.25], [med, med]; linewidth=3, label=false)
        end

        savefig(p, outpath)
    end

    scat_med(:V0, "Zero-sequence voltage magnitude | distribution across buses and timesteps",
             joinpath(outdir, "fig11_like_V0_scatter_median.png"))

    scat_med(:V2_over_V1, "Negative sequence ratio (V2/V1) | distribution across buses and timesteps",
             joinpath(outdir, "fig11_like_V2_over_V1_scatter_median.png"))
end

## --------------------------------------------------
# 6) Load baseline CSV and select PV-stress timesteps
# --------------------------------------------------
isfile(BASELINE_CSV) || error("Baseline CSV not found: $(BASELINE_CSV)")

df0 = CSV.read(BASELINE_CSV, DataFrame)
df0 = df0[df0.pf_status .== "LOCALLY_SOLVED", :]

# PV proxy computed here (answers the pv_pu_proxy question)
df0.pv_pu = pv_shape_simple(df0.time)
df0.score = df0.pv_pu ./ (df0.alpha_t .+ 1e-6)

# add effective PV injection for ALL timesteps (df0)
df0.pv_kw_eff = PV_KW_PER_PHASE .* df0.pv_pu

df_sel = sort(df0, :score, rev=true)
df_sel = unique(df_sel, :timestep)
df_sel = df_sel[1:min(N_PV_STRESS, nrow(df_sel)), :]
df_sel = sort(df_sel, :time)

if SKIP_NIGHT_PV
    df_sel = df_sel[df_sel.pv_pu .> 1e-8, :]
end

nrow(df_sel) == 0 && error("No PV-stress timesteps selected after filters")

df_sel.pv_kw_eff = PV_KW_PER_PHASE .* df_sel.pv_pu
println("\nSelected PV-stress timesteps: ", nrow(df_sel))

# Choose center timestep for 48h window
center_row = WINDOW_CENTER == "FIRST_STRESS" ? df_sel[1, :] : df_sel[argmax(df_sel.score), :]
center_k = Int(center_row.timestep)

# Build 48h window from baseline CSV time axis
k_list = collect(df0.timestep)
center_idx = findfirst(==(center_k), k_list)
center_idx === nothing && error("Center timestep not found in baseline CSV")

half = Int(floor(WINDOW_STEPS ÷ 2))
i1 = max(1, center_idx - half)
i2 = min(nrow(df0), i1 + WINDOW_STEPS - 1)

df_win = df0[i1:i2, :]
df_win.pv_pu = pv_shape_simple(df_win.time)
df_win.pv_kw_eff = PV_KW_PER_PHASE .* df_win.pv_pu

println("48h window rows: ", nrow(df_win), " | center timestep = ", center_k)

# --------------------------------------------------
# 7) Rank 1 PV-stress timestep detailed sequence analysis
# --------------------------------------------------
rank=1
r = df_sel[1, :]
k = Int(r.timestep)
a = Float64(r.alpha_t)
pv_pu = Float64(r.pv_pu)
pv_kw_eff = PV_KW_PER_PHASE * pv_pu
load_multiplier = 50

eng_base = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
scale_loads!(eng_base, a)
math_base = PMD.transform_data_model(eng_base; kron_reduce=true, phase_project=true)

for (i, bus) in math_base["bus"]
    bus["vmin"][1:3] .= 0.9
    bus["vmax"][1:3] .= 1.1
end
for (i, load) in math_base["load"]
    load["pd"] *= load_multiplier
end

pf_base = PMD.solve_mc_pf(math_base, PMD.IVRUPowerModel, ipopt)
m_base = pf_metrics(pf_base; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

# PV only
eng_pv = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
scale_loads!(eng_pv, a)
pv_load_id = first(keys(eng_pv["load"])) |> string
eng_pv["load"][pv_load_id]["pd_nom"] .-= pv_kw_eff
pv_bus = string(eng_pv["load"][pv_load_id]["bus"])

math_pv = PMD.transform_data_model(eng_pv; kron_reduce=true, phase_project=true)
for (_, bus) in math_pv["bus"]
    bus["vmin"][1:3] .= 0.9
    bus["vmax"][1:3] .= 1.1
end
for (i, load) in math_pv["load"]
    load["pd"] *= load_multiplier
end

pf_pv = PMD.solve_mc_pf(math_pv, PMD.IVRUPowerModel, ipopt)
m_pv = pf_metrics(pf_pv; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

# PV + STATCOM
eng_stc = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
scale_loads!(eng_stc, a)
pv_load_id = first(keys(eng_stc["load"])) |> string
eng_stc["load"][pv_load_id]["pd_nom"] .-= pv_kw_eff
pv_bus = string(eng_stc["load"][pv_load_id]["bus"])

math_stc = PMD.transform_data_model(eng_stc; kron_reduce=true, phase_project=true)
for (i, bus) in math_stc["bus"]
    bus["vmin"][1:3] .= 0.9
    bus["vmax"][1:3] .= 1.1
end
for (i, load) in math_stc["load"]
    load["pd"] *= load_multiplier
end

# Inverter gen ids expected in master_scaled.dss
stc_ids = [string(i) for i in 1:10]
for gen_id in stc_ids
    add_inverter_losses!(math_stc, gen_id; three_wire=true, c_rating_a=30*ones(3))
end

# Voltage unbalance constraint field used by the OPF builder for (_, bus_any) in math_stc["bus"]
for (_, bus_any) in math_stc["bus"]
    bus = bus_any::Dict{String,Any}
    bus["vm_vuf_max"] = 0.02
end

PMD.add_start_vrvi!(math_stc)
model = PMD.instantiate_mc_model(math_stc, PMD.IVRUPowerModel, build_mc_opf_mx_3w_Rhea);
pf_stc = PMD.optimize_model!(model, optimizer=ipopt)
m_stc = pf_metrics(pf_stc; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

# Extract and print feeder head current magnitudes
# Manually specify reference branch id for feeder head if known
# ref_branch_id = 176
# ref_branch_current_magnitude_base = abs.(pf_base["solution"]["branch"]["$ref_branch_id"]["cr_fr"] .+ im .* pf_base["solution"]["branch"]["$ref_branch_id"]["ci_fr"])
# ref_branch_current_magnitude_pv = abs.(pf_pv["solution"]["branch"]["$ref_branch_id"]["cr_fr"] .+ im .* pf_pv["solution"]["branch"]["$ref_branch_id"]["ci_fr"])
# ref_branch_current_magnitude_stc = abs.(pf_stc["solution"]["branch"]["$ref_branch_id"]["cr_fr"] .+ im .* pf_stc["solution"]["branch"]["$ref_branch_id"]["ci_fr"])

# pick source bus name (based on common naming conventions)
source_bus = pick_source_bus_name(math_base)

# pick feeder-head / substation branch id automatically
ref_branch_id = _pick_head_branch_id(pf_base, math_base, source_bus)
println("Using reference branch: ", ref_branch_id)

# branch current magnitude extraction function
function branch_current_magnitude(pf::Dict{String,Any}, bid)
    br = pf["solution"]["branch"][string(bid)]
    cr = br["cr_fr"]
    ci = br["ci_fr"]
    return abs.(cr .+ im .* ci)
end

# extract reference branch current magnitudes
ref_branch_current_magnitude_base = branch_current_magnitude(pf_base, ref_branch_id)
ref_branch_current_magnitude_pv = branch_current_magnitude(pf_pv, ref_branch_id)
ref_branch_current_magnitude_stc = branch_current_magnitude(pf_stc, ref_branch_id)

println("\nFeeder head current magnitudes (A):")
println(" Base max: ", maximum(ref_branch_current_magnitude_base))
println(" PV max: ", maximum(ref_branch_current_magnitude_pv))
println(" STC max: ", maximum(ref_branch_current_magnitude_stc))

# --------------------------------------------------
# Detailed sequence voltage analysis plots
# --------------------------------------------------
# Sequence voltage extraction for all buses
v_neg_base = Float64[]
v_zero_base = Float64[]
v_neg_pv = Float64[]
v_zero_pv = Float64[]
v_neg_stc = Float64[]
v_zero_stc = Float64[]

# Helper: convert vm/va to vr/vi
function vm_va_to_rectangular(vm::AbstractVector, va::AbstractVector)
    vr = vm .* cos.(va)
    vi = vm .* sin.(va)
    return (vr=vr, vi=vi)
end

for (i, bus) in pf_base["solution"]["bus"]
    # base results
    if haskey(bus, "vm")
        vr, vi = vm_va_to_rectangular(bus["vm"], bus["va"])
    else
        vr, vi = bus["vr"], bus["vi"]
    end
    _, _, v_seq_m = get_sequence_components(vr[1:3] + im .* vi[1:3])
    push!(v_neg_base, v_seq_m[3])
    push!(v_zero_base, v_seq_m[1])

    # pv results
    bus_res_pv = pf_pv["solution"]["bus"][i]
    if haskey(bus_res_pv, "vm")
        vr, vi = vm_va_to_rectangular(bus_res_pv["vm"], bus_res_pv["va"])
    else
        vr, vi = bus_res_pv["vr"], bus_res_pv["vi"]
    end
    _, _, v_seq_m = get_sequence_components(vr[1:3] + im .* vi[1:3])
    push!(v_neg_pv, v_seq_m[3])
    push!(v_zero_pv, v_seq_m[1])

    # statcom results
    bus_res_stc = pf_stc["solution"]["bus"][i]
    if haskey(bus_res_stc, "vm")
        vr, vi = vm_va_to_rectangular(bus_res_stc["vm"], bus_res_stc["va"])
    else
        vr, vi = bus_res_stc["vr"], bus_res_stc["vi"]
    end
    _, _, v_seq_m = get_sequence_components(vr[1:3] + im .* vi[1:3])
    push!(v_neg_stc, v_seq_m[3])
    push!(v_zero_stc, v_seq_m[1])
end

# Zero-sequence and negative-sequence voltage scatter plots
v_zero_plot = plot(v_zero_base; xlabel="bus index", label="base", title="Zero-sequence voltage magnitude", seriestype=:scatter)
plot!(v_zero_pv; label="pv", seriestype=:scatter)
plot!(v_zero_stc; label="stc", seriestype=:scatter)

v_neg_plot = plot(v_neg_base; xlabel="bus index", label="base", title="Negative-sequence voltage magnitude", seriestype=:scatter)
plot!(v_neg_pv; label="pv", seriestype=:scatter)
plot!(v_neg_stc; label="stc", seriestype=:scatter)

# Save zero-sequence scatter
savefig(v_zero_plot, joinpath(FIGDIR, "seq_V0_rank$(rank)_k$(k).png"))

# Save negative-sequence scatter
savefig(v_neg_plot, joinpath(FIGDIR, "seq_V2_rank$(rank)_k$(k).png"))

# Bus voltage magnitude per phase plots
vm_base = bus_vm_abc(pf_base)
vm_pv = bus_vm_abc(pf_pv)
vm_stc = bus_vm_abc(pf_stc)

# currently doing Phase A only as PV_PHASES = [1]
pA = plot(vm_base.a; label="base A", xlabel="Bus index", ylabel="Voltage (pu)", title="Bus voltage magnitude phase A | k=$(k)")
plot!(vm_pv.a; label="pv A")
plot!(vm_stc.a; label="stc A")
savefig(pA, joinpath(FIGDIR, "vm_phaseA_rank$(rank)_k$(k).png"))

# --------------------------------------------------
# 8) Main loop: detailed analysis on ALL selected PV-stress timesteps
# --------------------------------------------------
rows = NamedTuple[]
seq_rows = NamedTuple[]
br_seq_rows = NamedTuple[]

for (rank, r) in enumerate(eachrow(df_sel))   # iterate over selected PV-stress timesteps
# for (rank, r) in enumerate(eachrow(df0))   # iterate over ALL baseline timesteps (for sanity check and to fill seq_rows/br_seq_rows for all timesteps)
# for (rank, r) in enumerate(eachrow(df_win))  # iterate over 48h window timesteps around center PV-stress timestep
    println("\n==============================")
    println("Running rank = ", rank, " / ", nrow(df_sel))
    println("timestep = ", r.timestep, " | time = ", r.time)
    println("==============================")

    k = Int(r.timestep)
    a = Float64(r.alpha_t)
    pv_kw_eff = Float64(r.pv_kw_eff)
    load_multiplier = 50

    # ==================================================
    # BASELINE
    # ==================================================
    eng_base = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
    scale_loads!(eng_base, a)
    math_base = PMD.transform_data_model(eng_base; kron_reduce=true, phase_project=true)

    for (_, bus) in math_base["bus"]
        bus["vmin"][1:3] .= VMIN_PU
        bus["vmax"][1:3] .= VMAX_PU
    end
    for (_, load) in math_base["load"]
        load["pd"] *= load_multiplier
    end

    pf_base = PMD.solve_mc_pf(math_base, PMD.IVRUPowerModel, ipopt)
    m_base = pf_metrics(pf_base; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ==================================================
    # PV ONLY
    # ==================================================
    eng_pv = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
    scale_loads!(eng_pv, a)
    pv_load_id = first(keys(eng_pv["load"])) |> string
    eng_pv["load"][pv_load_id]["pd_nom"] .-= pv_kw_eff
    pv_bus = string(eng_pv["load"][pv_load_id]["bus"])

    math_pv = PMD.transform_data_model(eng_pv; kron_reduce=true, phase_project=true)
    for (_, bus) in math_pv["bus"]
        bus["vmin"][1:3] .= VMIN_PU
        bus["vmax"][1:3] .= VMAX_PU
    end
    for (_, load) in math_pv["load"]
        load["pd"] *= load_multiplier
    end

    pf_pv = PMD.solve_mc_pf(math_pv, PMD.IVRUPowerModel, ipopt)
    m_pv = pf_metrics(pf_pv; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ==================================================
    # PV + STATCOM
    # ==================================================
    eng_stc = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
    scale_loads!(eng_stc, a)
    pv_load_id = first(keys(eng_stc["load"])) |> string
    eng_stc["load"][pv_load_id]["pd_nom"] .-= pv_kw_eff

    math_stc = PMD.transform_data_model(eng_stc; kron_reduce=true, phase_project=true)

    for (_, bus_any) in math_stc["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vmin"][1:3] .= VMIN_PU
        bus["vmax"][1:3] .= VMAX_PU
        bus["vm_vuf_max"] = 0.02
    end
    for (_, load) in math_stc["load"]
        load["pd"] *= load_multiplier
    end

    stc_ids = [string(i) for i in 1:10]
    for gen_id in stc_ids
        add_inverter_losses!(math_stc, gen_id; three_wire=true, c_rating_a=30 .* ones(3))
    end

    PMD.add_start_vrvi!(math_stc)
    model = PMD.instantiate_mc_model(math_stc, PMD.IVRUPowerModel, build_mc_opf_mx_3w_Rhea)
    pf_stc = PMD.optimize_model!(model, optimizer=ipopt)
    m_stc = pf_metrics(pf_stc; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ==================================================
    # SUBSTATION / HEAD CURRENTS
    # ==================================================
    source_bus = pick_source_bus_name(eng_base)
    Ibase = feeder_head_current_metrics(pf_base, math_base, source_bus; vbase_ln=VBASE_LN)
    Ipv   = feeder_head_current_metrics(pf_pv,   math_pv,   source_bus; vbase_ln=VBASE_LN)
    Istc  = feeder_head_current_metrics(pf_stc,  math_stc,  source_bus; vbase_ln=VBASE_LN)

    # ==================================================
    # VOLTAGE + BRANCH SEQUENCES (Rahmat’s request)
    # ==================================================
    append!(seq_rows, voltage_sequence_rows(pf_base, r.time, "baseline"))
    append!(seq_rows, voltage_sequence_rows(pf_pv, r.time, "pv"))
    append!(seq_rows, voltage_sequence_rows(pf_stc, r.time, "statcom"))

    append!(br_seq_rows, branch_current_sequence_rows(pf_base, r.time, "baseline"))
    append!(br_seq_rows, branch_current_sequence_rows(pf_pv, r.time, "pv"))
    append!(br_seq_rows, branch_current_sequence_rows(pf_stc, r.time, "statcom"))

    # ==================================================
    # MAIN OUTPUT ROW
    # ==================================================
    push!(rows, (
        rank = rank,
        timestep = k,
        time = r.time,
        alpha_t = a,
        pv_bus = pv_bus,
        pv_kw_eff = pv_kw_eff,

        base_vmin = m_base.vmin,
        base_vmax = m_base.vmax,

        pv_vmin = m_pv.vmin,
        pv_vmax = m_pv.vmax,

        stc_vmin = m_stc.vmin,
        stc_vmax = m_stc.vmax,

        d_pv_vmin = m_pv.vmin - m_base.vmin,
        d_pv_vmax = m_pv.vmax - m_base.vmax,

        d_stc_vmin = m_stc.vmin - m_base.vmin,
        d_stc_vmax = m_stc.vmax - m_base.vmax,

        base_Imax = Ibase.Imax,
        pv_Imax = Ipv.Imax,
        stc_Imax = Istc.Imax,

        base_I2r = Ibase.I2_over_I1,
        pv_I2r = Ipv.I2_over_I1,
        stc_I2r = Istc.I2_over_I1
    ))

    println("Completed rank ", rank,
            " | base_vmin=", round(m_base.vmin, digits=4),
            " pv_vmin=", round(m_pv.vmin, digits=4),
            " stc_vmin=", round(m_stc.vmin, digits=4))
end

df_out = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps.csv"), df_out)
println("\nSaved table: ", joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps.csv"))

df_seq = DataFrame(seq_rows)
CSV.write(joinpath(TBLDIR, "voltage_sequence_components_selected_timesteps.csv"), df_seq)
println("Saved table: ", joinpath(TBLDIR, "voltage_sequence_components_selected_timesteps.csv"))

df_br_seq = DataFrame(br_seq_rows)
CSV.write(joinpath(TBLDIR, "branch_current_sequence_components_selected_timesteps.csv"), df_br_seq)
println("Saved table: ", joinpath(TBLDIR, "branch_current_sequence_components_selected_timesteps.csv"))

#### sanity check
println("\n--- sanity checks ---")
@show nrow(df_out)
@show nrow(df_seq)
@show combine(groupby(df_seq, :scenario), nrow => :N)
@show sum(isfinite.(df_seq.V0))
sum(isfinite.(df_seq.V2))
sum(isfinite.(df_seq.V2_over_V1))
println("----------------------\n")

# --------------------------------------------------
# analysis plots
# --------------------------------------------------
function duration_curve(values::AbstractVector{<:Real}; high_is_worse::Bool=true)
    v = Float64.(values)
    v = v[isfinite.(v)]
    isempty(v) && return (p=Float64[], x=Float64[])
    x = sort(v; rev=high_is_worse)
    p = collect(range(0, 100; length=length(x)))
    return (p=p, x=x)
end

function duration_curve_from_df(df::DataFrame, col::Symbol, scenario::String; high_is_worse::Bool=true)
    sub = df[df.scenario .== scenario, col]
    return duration_curve(sub; high_is_worse=high_is_worse)
end

df_out = sort(df_out, :time)

# df_out must contain: time, base_vmin, base_vmax, pv_vmin, pv_vmax, stc_vmin, stc_vmax
p = plot(; xlabel="Time", ylabel="Voltage magnitude (pu)", title="Voltage magnitude range (min–max) across time")

# Baseline band
plot!(p, df_out.time, df_out.base_vmin; label="baseline vmin")
plot!(p, df_out.time, df_out.base_vmax; label="baseline vmax", fillrange=df_out.base_vmin, fillalpha=0.15)

# PV band
plot!(p, df_out.time, df_out.pv_vmin; label="pv vmin")
plot!(p, df_out.time, df_out.pv_vmax; label="pv vmax", fillrange=df_out.pv_vmin, fillalpha=0.15)

# STATCOM band
plot!(p, df_out.time, df_out.stc_vmin; label="statcom vmin")
plot!(p, df_out.time, df_out.stc_vmax; label="statcom vmax", fillrange=df_out.stc_vmin, fillalpha=0.15)

savefig(p, joinpath(FIGDIR, "voltage_minmax_band_timeseries.png"))

# --- duration curves across ALL selected timesteps (uses df_seq) ---
dcV0_base = duration_curve_from_df(df_seq, :V0, "baseline")
dcV0_pv   = duration_curve_from_df(df_seq, :V0, "pv")
dcV0_stc  = duration_curve_from_df(df_seq, :V0, "statcom")

pV0 = plot(dcV0_base.p, dcV0_base.x; label="baseline", xlabel="Percentile (%)", ylabel="|V0| (pu)",
           title="Zero-sequence voltage duration curve (all selected timesteps)")
plot!(pV0, dcV0_pv.p, dcV0_pv.x; label="pv")
plot!(pV0, dcV0_stc.p, dcV0_stc.x; label="statcom")
savefig(pV0, joinpath(FIGDIR, "dc_V0_selected_timesteps.png"))

dcV2_base = duration_curve_from_df(df_seq, :V2, "baseline")
dcV2_pv   = duration_curve_from_df(df_seq, :V2, "pv")
dcV2_stc  = duration_curve_from_df(df_seq, :V2, "statcom")

pV2 = plot(dcV2_base.p, dcV2_base.x; label="baseline", xlabel="Percentile (%)", ylabel="|V2| (pu)",
           title="Negative-sequence voltage duration curve (all selected timesteps)")
plot!(pV2, dcV2_pv.p, dcV2_pv.x; label="pv")
plot!(pV2, dcV2_stc.p, dcV2_stc.x; label="statcom")
savefig(pV2, joinpath(FIGDIR, "dc_V2_selected_timesteps.png"))

# Min voltage vs time
p_vmin = plot(df_out.time, df_out.base_vmin; label="baseline", xlabel="Time", ylabel="Vmin (pu)",
              title="Minimum bus voltage over time", legend=:bottomleft)
plot!(p_vmin, df_out.time, df_out.pv_vmin; label="pv")
plot!(p_vmin, df_out.time, df_out.stc_vmin; label="statcom")
savefig(p_vmin, joinpath(FIGDIR, "timeseries_vmin.png"))

# Max voltage vs time
p_vmax = plot(df_out.time, df_out.base_vmax; label="baseline", xlabel="Time", ylabel="Vmax (pu)",
              title="Maximum bus voltage over time", legend=:bottomleft)
plot!(p_vmax, df_out.time, df_out.pv_vmax; label="pv")
plot!(p_vmax, df_out.time, df_out.stc_vmax; label="statcom")
savefig(p_vmax, joinpath(FIGDIR, "timeseries_vmax.png"))

# mean lines + shaded min–max bands + zoom inse
# min/max across phases for the band
phase_min(Ia, Ib, Ic) = min(Ia, Ib, Ic)
phase_max(Ia, Ib, Ic) = max(Ia, Ib, Ic)

# 48-hour window: compute feeder-head currents (baseline vs STATCOM) for paper-style Fig.7
rows_win = NamedTuple[]

for r in eachrow(df_sel)
    k = Int(r.timestep)
    a = Float64(r.alpha_t)
    pv_kw_eff = Float64(r.pv_kw_eff)
    load_multiplier = 50

    # =========================
    # BASELINE
    # =========================
    eng_base = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
    scale_loads!(eng_base, a)
    math_base = PMD.transform_data_model(eng_base; kron_reduce=true, phase_project=true)

    for (_, bus) in math_base["bus"]
        bus["vmin"][1:3] .= VMIN_PU
        bus["vmax"][1:3] .= VMAX_PU
    end
    for (_, load) in math_base["load"]
        load["pd"] *= load_multiplier
    end

    pf_base = PMD.solve_mc_pf(math_base, PMD.IVRUPowerModel, ipopt)

    # =========================
    # PV + STATCOM (your "after")
    # =========================
    eng_stc = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])
    scale_loads!(eng_stc, a)
    pv_load_id = first(keys(eng_stc["load"])) |> string
    eng_stc["load"][pv_load_id]["pd_nom"] .-= pv_kw_eff

    math_stc = PMD.transform_data_model(eng_stc; kron_reduce=true, phase_project=true)

    for (_, bus_any) in math_stc["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vmin"][1:3] .= VMIN_PU
        bus["vmax"][1:3] .= VMAX_PU
        bus["vm_vuf_max"] = 0.02
    end
    for (_, load) in math_stc["load"]
        load["pd"] *= load_multiplier
    end

    for gen_id in [string(i) for i in 1:10]
        add_inverter_losses!(math_stc, gen_id; three_wire=true, c_rating_a=30 .* ones(3))
    end

    PMD.add_start_vrvi!(math_stc)
    model = PMD.instantiate_mc_model(math_stc, PMD.IVRUPowerModel, build_mc_opf_mx_3w_Rhea)
    pf_stc = PMD.optimize_model!(model, optimizer=ipopt)

    # =========================
    # feeder-head current metrics
    # =========================
    source_bus = pick_source_bus_name(math_base)  # <-- use math_base here (more consistent)
    Ibase = feeder_head_current_metrics(pf_base, math_base, source_bus; vbase_ln=VBASE_LN)
    Istc  = feeder_head_current_metrics(pf_stc,  math_stc,  source_bus; vbase_ln=VBASE_LN)

    push!(rows_win, (
        time = r.time,
        mean_before = Ibase.Imean,
        min_before = phase_min(Ibase.Ia, Ibase.Ib, Ibase.Ic),
        max_before = phase_max(Ibase.Ia, Ibase.Ib, Ibase.Ic),
        mean_after = Istc.Imean,
        min_after = phase_min(Istc.Ia, Istc.Ib, Istc.Ic),
        max_after = phase_max(Istc.Ia, Istc.Ib, Istc.Ic),
    ))
end

df_win_curr = DataFrame(rows_win)
df_win_curr = sort(df_win_curr, :time)
CSV.write(joinpath(TBLDIR, "currents_48h_before_after.csv"), df_win_curr)
println("Saved: ", joinpath(TBLDIR, "currents_48h_before_after.csv"))

function plot_substation_current_envelope_48h(df::DataFrame; outpath::String,
    ylab::String="Feeder F2 Substation Current (A)",
    title_str::String="Substation range of phase currents before and after STATCOM",
    inset_hours::Tuple{Int,Int}=(5,7))

    # Main plot
    p = plot(df.time, df.mean_after; label="Mean after", linewidth=2,
             xlabel="Time", ylabel=ylab, title=title_str)
    plot!(df.time, df.mean_before; label="Mean before", linewidth=2)

    # shaded bands:
    plot!(df.time, df.max_after; label="After STATCOM (range)",
          fillrange=df.min_after, fillalpha=0.20, linewidth=0)
    plot!(df.time, df.max_before; label="Before STATCOM (range)",
          fillrange=df.min_before, fillalpha=0.20, linewidth=0)

    # ---------- inset zoom (simple time-of-day window) ----------
    h1, h2 = inset_hours
    mask = [(hour(t) >= h1 && hour(t) <= h2) for t in df.time]
    dfz = df[mask, :]

    pin = plot(dfz.time, dfz.mean_after; label=false, linewidth=2, title="", xlabel="", ylabel="")
    plot!(pin, dfz.time, dfz.mean_before; label=false, linewidth=2)
    plot!(pin, dfz.time, dfz.max_after; fillrange=dfz.min_after, fillalpha=0.20, linewidth=0, label=false)
    plot!(pin, dfz.time, dfz.max_before; fillrange=dfz.min_before, fillalpha=0.20, linewidth=0, label=false)

    # Inset placement 
    plot!(p,
        inset_subplots = [
            (1, bbox(0.62, 0.62, 0.35, 0.30, :left, :bottom))
        ]
    )

    # Draw into subplot 2 (the inset)
    plot!(p, dfz.time, dfz.mean_after;
        subplot=2, label=false, linewidth=2)

    plot!(p, dfz.time, dfz.mean_before;
        subplot=2, label=false, linewidth=2)

    plot!(p, dfz.time, dfz.max_after;
        subplot=2, fillrange=dfz.min_after,
        fillalpha=0.20, linewidth=0, label=false)

    plot!(p, dfz.time, dfz.max_before;
        subplot=2, fillrange=dfz.min_before,
        fillalpha=0.20, linewidth=0, label=false)


    savefig(p, outpath)
end

outpath = "/home/auc009/test_plot.png"

plot_substation_current_envelope_48h(df_win_curr;
    outpath=joinpath(FIGDIR, "substation_current_envelope_before_after_statcom_48h.png"))

# savefig(p, joinpath(FIGDIR, "substation_current_envelope_before_after_statcom_48h.png"))