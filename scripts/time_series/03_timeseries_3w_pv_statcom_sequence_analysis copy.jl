# ================================================================
# Time-series comparison on selected PV-stress timesteps and a
# contiguous 48-hour window for paper-style current plots:
#
# (1) Baseline (NO PV, NO STATCOM)
# (2) PV only (PV as negative load)
# (3) PV + STATCOM (OPF w inverter gens)
#
# Added for Rahmat:
# - Plot current on objective/reference branch 176 (time series + duration curve)
# - Network-wide max branch phase current over time + duration curve
# - Network-wide current sequence improvement: p95(I2/I1) across branches over time
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
# MASTER_STATCOM  = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_scaled.dss")

MASTER_BASELINE   = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w_new.dss")
# MASTER_BASELINE = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w.dss")

MASTER_BASELINE_3w = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", "spd_s", "master_scaled_new.dss")

YEAR   = 2023
STRIDE = 4

BASELINE_OUTDIR = joinpath(
    ROOT, "results", "time_series", "baseline_pf", NET_4W,
    "year=$(YEAR)_stride=$(STRIDE)_K=17520"
)
BASELINE_CSV = joinpath(BASELINE_OUTDIR, "tables", "timeseries_baseline_pf_metrics.csv")

# PV settings
PV_BUS = "AUTO"
PV_KW_PER_PHASE = 10.0
PV_PHASES = [1]

# PV-stress selection
N_PV_STRESS = 30
SKIP_NIGHT_PV = true
RANDOM_SEED = 42

# 48-hour window
WINDOW_HOURS = 48
WINDOW_STEPS = Int(WINDOW_HOURS * 60 ÷ 30) + 1
WINDOW_CENTER = "BEST_STRESS"  # "BEST_STRESS" or "FIRST_STRESS"

# Voltage limits + conversion
VMIN_PU = 0.90
VMAX_PU = 1.10
VBASE_LN = 230.0

# Baseline engineering scaling applied once, then alpha(t)
LOAD_ALPHA_BASE = 1.5

# Objective/reference branch
REF_BRANCH_ID = "176"   # <- Rahmat’s objective branch

# Solver
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

OUTDIR = joinpath(
    ROOT, "results", "time_series", "pv_pf", "3w",
    "year=$(YEAR)_stride=$(STRIDE)",
    "pvkw=$(PV_KW_PER_PHASE)_ph=$(join(PV_PHASES, '_'))_K=$(N_PV_STRESS)"
)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

# --------------------------------------------------
# 1) Helpers: scaling, metrics, PV proxy
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

function pf_metrics(pf::Dict{String,Any}; vmin_pu=0.90, vmax_pu=1.10)
    sol_bus = pf["solution"]["bus"]
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
    isempty(vmins) && return (status=string(pf["termination_status"]), vmin=NaN, vmax=NaN, n_under=0, n_over=0)
    vmin_sys = minimum(vmins)
    vmax_sys = maximum(vmaxs)
    n_under = count(x -> x < vmin_pu, vmins)
    n_over  = count(x -> x > vmax_pu, vmaxs)
    return (status=string(pf["termination_status"]), vmin=vmin_sys, vmax=vmax_sys, n_under=n_under, n_over=n_over)
end

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
# 2) Helpers: current extraction + sequences (NETWORK METRICS)
# --------------------------------------------------
function _phasor_from(br::Dict{String,Any}, rkey::String, ikey::String, idx::Int)
    (haskey(br, rkey) && haskey(br, ikey)) || return nothing
    r = br[rkey]; im = br[ikey]
    (r isa AbstractVector && im isa AbstractVector) || return nothing
    (length(r) >= idx && length(im) >= idx) || return nothing
    return complex(Float64(r[idx]), Float64(im[idx]))
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

function _seq_components(Aa::Complex, Ab::Complex, Ac::Complex)
    a = cis(2pi/3)
    A0 = (Aa + Ab + Ac)/3
    A1 = (Aa + a*Ab + a^2*Ac)/3
    A2 = (Aa + a^2*Ab + a*Ac)/3
    return (A0=A0, A1=A1, A2=A2)
end

"""
Return max phase current magnitude for a specific branch id (string).
If branch id not present / currents not found -> NaN.
"""
function branch_Imax(pf::Dict{String,Any}, branch_id::String)
    sol = pf["solution"]
    haskey(sol, "branch") || return NaN
    brdict = sol["branch"]
    if !haskey(brdict, branch_id)
        return NaN
    end
    br = brdict[branch_id]::Dict{String,Any}
    ph = _branch_current_phasors(br)
    ph === nothing && return NaN
    return maximum([abs(ph.Ia), abs(ph.Ib), abs(ph.Ic)])
end

"""
Network-wide maximum of branch max-phase-current (per timestep).
"""
function network_max_branch_Imax(pf::Dict{String,Any})
    sol = pf["solution"]
    haskey(sol, "branch") || return NaN
    isempty(sol["branch"]) && return NaN
    vals = Float64[]
    for (_, br_any) in sol["branch"]
        br = br_any::Dict{String,Any}
        ph = _branch_current_phasors(br)
        ph === nothing && continue
        push!(vals, maximum([abs(ph.Ia), abs(ph.Ib), abs(ph.Ic)]))
    end
    isempty(vals) ? NaN : maximum(vals)
end

"""
Network-wide p95 of current unbalance ratio I2/I1 across branches (per timestep).
"""
function network_p95_I2_over_I1(pf::Dict{String,Any})
    sol = pf["solution"]
    haskey(sol, "branch") || return NaN
    isempty(sol["branch"]) && return NaN
    ratios = Float64[]
    for (_, br_any) in sol["branch"]
        br = br_any::Dict{String,Any}
        ph = _branch_current_phasors(br)
        ph === nothing && continue
        Ia, Ib, Ic = ph.Ia, ph.Ib, ph.Ic
        seq = _seq_components(Ia, Ib, Ic)
        I1 = abs(seq.A1)
        I2 = abs(seq.A2)
        r = I1 > 1e-12 ? (I2 / I1) : NaN
        isfinite(r) && push!(ratios, r)
    end
    isempty(ratios) ? NaN : quantile(ratios, 0.95)
end

# --------------------------------------------------
# 3) Duration curve helper
# --------------------------------------------------
function duration_curve(values::AbstractVector{<:Real}; high_is_worse::Bool=true)
    v = Float64.(values)
    v = v[isfinite.(v)]
    isempty(v) && return (p=Float64[], x=Float64[])
    x = sort(v; rev=high_is_worse)
    p = collect(range(0, 100; length=length(x)))
    return (p=p, x=x)
end

# --------------------------------------------------
# 4) Load baseline CSV and select PV-stress timesteps
# --------------------------------------------------
isfile(BASELINE_CSV) || error("Baseline CSV not found: $(BASELINE_CSV)")

df0 = CSV.read(BASELINE_CSV, DataFrame)
df0 = df0[df0.pf_status .== "LOCALLY_SOLVED", :]

df0.pv_pu = pv_shape_simple(df0.time)
df0.score = df0.pv_pu ./ (df0.alpha_t .+ 1e-6)
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

# 48h window around the chosen center
center_row = WINDOW_CENTER == "FIRST_STRESS" ? df_sel[1, :] : df_sel[argmax(df_sel.score), :]
center_k = Int(center_row.timestep)
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
# 5) Main loop: solve + extract Rahmat metrics
# --------------------------------------------------
rows = NamedTuple[]
for (rank, r) in enumerate(eachrow(df_sel))
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

    for (_, bus_any) in math_base["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vmin"][1:3] .= VMIN_PU
        bus["vmax"][1:3] .= VMAX_PU
    end
    for (_, load_any) in math_base["load"]
        load = load_any::Dict{String,Any}
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
    for (_, bus_any) in math_pv["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vmin"][1:3] .= VMIN_PU
        bus["vmax"][1:3] .= VMAX_PU
    end
    for (_, load_any) in math_pv["load"]
        load = load_any::Dict{String,Any}
        load["pd"] *= load_multiplier
    end

    pf_pv = PMD.solve_mc_pf(math_pv, PMD.IVRUPowerModel, ipopt)
    m_pv = pf_metrics(pf_pv; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ==================================================
    # PV + STATCOM (OPF w inverter gens)
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
    for (_, load_any) in math_stc["load"]
        load = load_any::Dict{String,Any}
        load["pd"] *= load_multiplier
    end

    for gen_id in [string(i) for i in 1:10]
        add_inverter_losses!(math_stc, gen_id; three_wire=true, c_rating_a=30 .* ones(3))
    end

    PMD.add_start_vrvi!(math_stc)
    model = PMD.instantiate_mc_model(math_stc, PMD.IVRUPowerModel, build_mc_opf_mx_3w_Rhea)
    pf_stc = PMD.optimize_model!(model, optimizer=ipopt)
    m_stc = pf_metrics(pf_stc; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # ==================================================
    # Rahmat metrics:
    # (A) Branch 176 current (objective)
    # (B) Network-wide max branch phase current
    # (C) Network-wide p95(I2/I1) across branches
    # ==================================================
    ref176_base = branch_Imax(pf_base, REF_BRANCH_ID)
    ref176_pv   = branch_Imax(pf_pv,   REF_BRANCH_ID)
    ref176_stc  = branch_Imax(pf_stc,  REF_BRANCH_ID)

    netI_base = network_max_branch_Imax(pf_base)
    netI_pv   = network_max_branch_Imax(pf_pv)
    netI_stc  = network_max_branch_Imax(pf_stc)

    netI2r_base = network_p95_I2_over_I1(pf_base)
    netI2r_pv   = network_p95_I2_over_I1(pf_pv)
    netI2r_stc  = network_p95_I2_over_I1(pf_stc)

    # ==================================================
    # Output row
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
        pv_vmin   = m_pv.vmin,
        pv_vmax   = m_pv.vmax,
        stc_vmin  = m_stc.vmin,
        stc_vmax  = m_stc.vmax,

        ref176_Imax_base = ref176_base,
        ref176_Imax_pv   = ref176_pv,
        ref176_Imax_stc  = ref176_stc,

        net_Imax_base = netI_base,
        net_Imax_pv   = netI_pv,
        net_Imax_stc  = netI_stc,

        net_I2r_p95_base = netI2r_base,
        net_I2r_p95_pv   = netI2r_pv,
        net_I2r_p95_stc  = netI2r_stc
    ))

    println("Done rank ", rank,
            " | ref176 base/pv/stc = ",
            round(ref176_base, digits=4), " / ",
            round(ref176_pv, digits=4), " / ",
            round(ref176_stc, digits=4),
            " | netI base/pv/stc = ",
            round(netI_base, digits=4), " / ",
            round(netI_pv, digits=4), " / ",
            round(netI_stc, digits=4),
            " | p95(I2/I1) base/pv/stc = ",
            round(netI2r_base, digits=4), " / ",
            round(netI2r_pv, digits=4), " / ",
            round(netI2r_stc, digits=4)
    )
end

df_out = DataFrame(rows)
df_out = sort(df_out, :time)

CSV.write(joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps_with_currents.csv"), df_out)
println("\nSaved table: ", joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps_with_currents.csv"))

# --------------------------------------------------
# 6) Plots Rahmat asked for
# --------------------------------------------------

# --- Branch 176 time series ---
p176 = plot(df_out.time, df_out.ref176_Imax_base; label="baseline",
            xlabel="Time", ylabel="Branch $(REF_BRANCH_ID) Imax (A)",
            title="Objective branch current (branch $(REF_BRANCH_ID)) over time", legend=:bottomleft)
plot!(p176, df_out.time, df_out.ref176_Imax_pv;  label="pv")
plot!(p176, df_out.time, df_out.ref176_Imax_stc; label="statcom")
savefig(p176, joinpath(FIGDIR, "timeseries_ref176_Imax.png"))

# --- Branch 176 duration curve ---
dc176_base = duration_curve(df_out.ref176_Imax_base; high_is_worse=true)
dc176_pv   = duration_curve(df_out.ref176_Imax_pv;   high_is_worse=true)
dc176_stc  = duration_curve(df_out.ref176_Imax_stc;  high_is_worse=true)

p176dc = plot(dc176_base.p, dc176_base.x; label="baseline",
              xlabel="Percentile (%)", ylabel="Branch $(REF_BRANCH_ID) Imax (A)",
              title="Duration curve: objective branch current ($(REF_BRANCH_ID))")
plot!(p176dc, dc176_pv.p,  dc176_pv.x;  label="pv")
plot!(p176dc, dc176_stc.p, dc176_stc.x; label="statcom")
savefig(p176dc, joinpath(FIGDIR, "dc_ref176_Imax.png"))

# --- Network-wide max branch current time series ---
pNet = plot(df_out.time, df_out.net_Imax_base; label="baseline",
            xlabel="Time", ylabel="Network max branch phase current (A)",
            title="Network-wide maximum branch phase current over time", legend=:bottomleft)
plot!(pNet, df_out.time, df_out.net_Imax_pv;  label="pv")
plot!(pNet, df_out.time, df_out.net_Imax_stc; label="statcom")
savefig(pNet, joinpath(FIGDIR, "timeseries_network_Imax.png"))

# --- Network-wide max branch current duration curve ---
dcNet_base = duration_curve(df_out.net_Imax_base; high_is_worse=true)
dcNet_pv   = duration_curve(df_out.net_Imax_pv;   high_is_worse=true)
dcNet_stc  = duration_curve(df_out.net_Imax_stc;  high_is_worse=true)

pNetdc = plot(dcNet_base.p, dcNet_base.x; label="baseline",
              xlabel="Percentile (%)", ylabel="Network max branch phase current (A)",
              title="Duration curve: network max branch phase current")
plot!(pNetdc, dcNet_pv.p,  dcNet_pv.x;  label="pv")
plot!(pNetdc, dcNet_stc.p, dcNet_stc.x; label="statcom")
savefig(pNetdc, joinpath(FIGDIR, "dc_network_Imax.png"))

# --- Network-wide I2/I1 p95 time series ---
pI2 = plot(df_out.time, df_out.net_I2r_p95_base; label="baseline",
           xlabel="Time", ylabel="I2/I1 (p95 across branches)",
           title="Unbalance current ratio across network (95th percentile)", legend=:bottomleft)
plot!(pI2, df_out.time, df_out.net_I2r_p95_pv;  label="pv")
plot!(pI2, df_out.time, df_out.net_I2r_p95_stc; label="statcom")
savefig(pI2, joinpath(FIGDIR, "timeseries_network_I2r_p95.png"))

# (Optional) duration curve of p95(I2/I1) across timesteps
dcI2_base = duration_curve(df_out.net_I2r_p95_base; high_is_worse=true)
dcI2_pv   = duration_curve(df_out.net_I2r_p95_pv;   high_is_worse=true)
dcI2_stc  = duration_curve(df_out.net_I2r_p95_stc;  high_is_worse=true)

pI2dc = plot(dcI2_base.p, dcI2_base.x; label="baseline",
             xlabel="Percentile (%)", ylabel="p95(I2/I1)",
             title="Duration curve: p95 unbalance current ratio across network")
plot!(pI2dc, dcI2_pv.p,  dcI2_pv.x;  label="pv")
plot!(pI2dc, dcI2_stc.p, dcI2_stc.x; label="statcom")
savefig(pI2dc, joinpath(FIGDIR, "dc_network_I2r_p95.png"))

# Keep voltage plots if you want, but Rahmat said not interesting:
pVmin = plot(df_out.time, df_out.base_vmin; label="baseline",
             xlabel="Time", ylabel="Vmin (pu)", title="Minimum bus voltage over time", legend=:bottomleft)
plot!(pVmin, df_out.time, df_out.pv_vmin;  label="pv")
plot!(pVmin, df_out.time, df_out.stc_vmin; label="statcom")
savefig(pVmin, joinpath(FIGDIR, "timeseries_vmin.png"))

pVmax = plot(df_out.time, df_out.base_vmax; label="baseline",
             xlabel="Time", ylabel="Vmax (pu)", title="Maximum bus voltage over time", legend=:bottomleft)
plot!(pVmax, df_out.time, df_out.pv_vmax;  label="pv")
plot!(pVmax, df_out.time, df_out.stc_vmax; label="statcom")
savefig(pVmax, joinpath(FIGDIR, "timeseries_vmax.png"))

println("\nSaved figures to: ", FIGDIR)
println("Done.")
