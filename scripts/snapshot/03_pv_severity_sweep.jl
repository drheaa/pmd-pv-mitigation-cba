# scripts/03_pv_severity_sweep.jl
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PowerModelsDistribution
using JuMP
using Ipopt
using DataFrames
using CSV
using Plots

const PMD = PowerModelsDistribution
PMD.silence!()

# ============================================================
# 03) PV severity sweep (snapshot PF)
#
# Part A: Sweep PV size at a fixed location (far-end bus)
#         Record vmin/vmax + violation counts.
#
# Part B: Compare clustered vs distributed placements (same PV size)
#         - Clustered: PV1 + PV2 at far-end bus
#         - Distributed: PV1 at far-end, PV2 at mid bus
# ============================================================

# -----------------------------
# 0) Network + paths
# -----------------------------
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"
NET  = "spd_s"  # change to spd_r/spd_u/spm_r/... when you loop later

master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_scaled.dss")

OUTDIR = joinpath(ROOT, "results", "pv_severity_sweep", NET)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

# -----------------------------
# 1) Study settings
# -----------------------------
LOAD_ALPHA = 1.0

VBASE_LN = 230.0
VMIN_PU  = 0.94
VMAX_PU  = 1.10

# Part A: PV sweep sizes (kW, single-phase PV modeled as negative load)
PV_SWEEP_KW = collect(0.0:1.0:15.0)

# Part B: pick one PV size for placement comparison
PV2_KW = 5.0

# -----------------------------
# Helpers
# -----------------------------
count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0

normalize_bus(s::AbstractString) = lowercase(split(String(s), ".")[1])

function pick_source_bus_eng(eng)::String
    haskey(eng["bus"], "sourcebusz") && return "sourcebusz"
    haskey(eng["bus"], "sourcebusZ") && return "sourcebusZ"
    haskey(eng["bus"], "sourcebus")  && return "sourcebus"
    return first(keys(eng["bus"]))
end

function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng
    for (_, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        for k in ("pd","qd","pd_nom","qd_nom","kW","kvar","kw","kvar")
            if haskey(ld, k)
                v = ld[k]
                ld[k] = v isa Number ? alpha * v : alpha .* v
            end
        end
    end
    return eng
end

# ---------- distances (lines + transformer) ----------
function make_edges_from_lines(eng)
    rows = NamedTuple[]
    haskey(eng, "line") || return DataFrame(rows)
    for (id, ln) in eng["line"]
        b1 = normalize_bus(string(get(ln, "f_bus", "")))
        b2 = normalize_bus(string(get(ln, "t_bus", "")))
        (isempty(b1) || isempty(b2)) && continue
        len = Float64(get(ln, "length", 0.0)) / 1000.0
        if !isfinite(len) || len < 0
            len = 0.0
        end
        push!(rows, (Bus1=b1, Bus2=b2, length_km=len, kind="line"))
    end
    return DataFrame(rows)
end

function make_edges_from_transformers(eng)
    rows = NamedTuple[]
    haskey(eng, "transformer") || return DataFrame(rows)
    for (id, tx) in eng["transformer"]
        buses =
            haskey(tx, "bus")   ? tx["bus"] :
            haskey(tx, "buses") ? tx["buses"] :
            nothing
        buses === nothing && continue
        length(buses) < 2 && continue
        b1 = normalize_bus(string(buses[1]))
        b2 = normalize_bus(string(buses[2]))
        push!(rows, (Bus1=b1, Bus2=b2, length_km=0.0, kind="transformer"))
    end
    return DataFrame(rows)
end

function compute_bus_distances(edges_df::DataFrame; source_bus::String)
    adj = Dict{String, Vector{Tuple{String, Float64}}}()
    for r in eachrow(edges_df)
        push!(get!(adj, r.Bus1, Tuple{String,Float64}[]), (r.Bus2, r.length_km))
        push!(get!(adj, r.Bus2, Tuple{String,Float64}[]), (r.Bus1, r.length_km))
    end

    dist = Dict{String,Float64}(source_bus => 0.0)
    q = [source_bus]
    while !isempty(q)
        u = popfirst!(q)
        for (v,w) in get(adj, u, [])
            if !haskey(dist, v)
                dist[v] = dist[u] + w
                push!(q, v)
            end
        end
    end
    return dist
end

# ---------- PV placement on customer loads ----------
"""
Collect candidate PV attachment points:
- only "load_*" (customer loads)
- single-phase only (phases == 1)
- must have an active power key we can modify (pd/pd_nom/kW/kw/p/p_kw)
Returns NamedTuples with (load_id, bus, phase).
"""
function collect_customer_1ph_load_points(eng::Dict{String,Any})
    pts = NamedTuple[]
    haskey(eng, "load") || return pts

    for (lid_any, ld_any) in eng["load"]
        lid = String(lid_any)
        startswith(lid, "load_") || continue

        ld = ld_any::Dict{String,Any}
        haskey(ld, "bus") || continue

        # single-phase only
        nph = get(ld, "phases", 1)
        (nph isa Number && Int(nph) == 1) || continue

        bus_raw  = String(ld["bus"])
        bus_base = normalize_bus(bus_raw)

        # phase from "bus.2" if present, else 1
        parts = split(bus_raw, ".")
        ph = (length(parts) >= 2 && tryparse(Int, parts[2]) !== nothing) ? clamp(parse(Int, parts[2]), 1, 3) : 1

        # active power key availability
        has_p =
            haskey(ld, "pd") || haskey(ld, "pd_nom") ||
            haskey(ld, "kW") || haskey(ld, "kw") ||
            haskey(ld, "p")  || haskey(ld, "p_kw")

        has_p || continue

        push!(pts, (load_id=lid, bus=bus_base, phase=ph))
    end

    return pts
end

"""
Subtract pv_kw from the load's active power (single-phase PV as negative load).
Tries keys in a safe order.
"""
function apply_pv_to_load!(eng::Dict{String,Any}, load_id::String; pv_kw::Float64)
    ld = eng["load"][load_id]::Dict{String,Any}

    function sub_kw!(key::String)
        v = ld[key]
        if v isa Number
            ld[key] = Float64(v) - pv_kw
            return true
        elseif v isa AbstractVector
            vv = Float64.(v)
            vv[1] -= pv_kw
            ld[key] = vv
            return true
        end
        return false
    end

    for k in ("pd", "pd_nom", "kW", "kw", "p", "p_kw")
        if haskey(ld, k) && sub_kw!(k)
            return nothing
        end
    end

    error("Load $load_id has no usable active power key. Keys=$(collect(keys(ld)))")
end

# ---------- PF solve + metrics ----------
function run_pf_metrics(eng::Dict{String,Any})
    math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)
    ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

    pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)
    status = pf["termination_status"]

    sol_bus = pf["solution"]["bus"]

    vmins = Float64[]
    vmaxs = Float64[]

    for sb in values(sol_bus)
        vm =
            haskey(sb, "vm") ? sb["vm"] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"].^2 .+ sb["vi"].^2) :
            nothing
        vm === nothing && continue

        push!(vmins, minimum(vm))
        push!(vmaxs, maximum(vm))
    end

    isempty(vmins) && return (status=status, vmin=NaN, vmax=NaN, n_under=0, n_over=0)

    vmin_sys = minimum(vmins)
    vmax_sys = maximum(vmaxs)

    n_under = count(x -> x < VMIN_PU, vmins)
    n_over  = count(x -> x > VMAX_PU, vmaxs)

    return (status=status, vmin=vmin_sys, vmax=vmax_sys, n_under=n_under, n_over=n_over)
end

# ============================================================
# Main
# ============================================================
println("Parsing: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

println("eng counts: buses=", count_dict(eng0,"bus"),
        " lines=", count_dict(eng0,"line"),
        " loads=", count_dict(eng0,"load"),
        " transformers=", count_dict(eng0,"transformer"))

# match your baseline habit
eng0["voltage_source"]["source"]["rs"] .= 0
eng0["voltage_source"]["source"]["xs"] .= 0

scale_loads!(eng0, LOAD_ALPHA)

SOURCE_BUS = normalize_bus(pick_source_bus_eng(eng0))

edges = vcat(make_edges_from_lines(eng0), make_edges_from_transformers(eng0))
dist  = compute_bus_distances(edges; source_bus=SOURCE_BUS)

pts = collect_customer_1ph_load_points(eng0)
isempty(pts) && error("No customer 1ϕ loads found (load_*) with active power keys.")

# rank points by distance
pts_ranked = [(; p..., d=get(dist, p.bus, Inf)) for p in pts]
sort!(pts_ranked; by = x -> x.d)

# pick mid and far based on ranked distances
mid_i = clamp(Int(round(length(pts_ranked) / 2)), 1, length(pts_ranked))
far   = pts_ranked[end]
mid   = pts_ranked[mid_i]

println("Picked PV attachment points (customer load buses, 1ϕ):")
println("  mid: load=$(mid.load_id) bus=$(mid.bus) phase=$(mid.phase) d=$(mid.d) km")
println("  far: load=$(far.load_id) bus=$(far.bus) phase=$(far.phase) d=$(far.d) km")

# ------------------------------------------------------------
# Baseline PF
# ------------------------------------------------------------
base = run_pf_metrics(deepcopy(eng0))
println("Baseline PF status: ", base.status,
        " | vmin=", base.vmin, " vmax=", base.vmax,
        " n_under=", base.n_under, " n_over=", base.n_over)

# ============================================================
# Part A: Sweep PV size at FAR
# ============================================================
rowsA = NamedTuple[]
for pv_kw in PV_SWEEP_KW
    eng = deepcopy(eng0)
    if pv_kw > 0
        apply_pv_to_load!(eng, far.load_id; pv_kw=pv_kw)
    end
    r = run_pf_metrics(eng)

    push!(rowsA, (
        net = NET,
        location = "far",
        pv_bus = far.bus,
        pv_load_id = far.load_id,
        pv_phase = far.phase,
        pv_kw = pv_kw,
        load_alpha = LOAD_ALPHA,
        status = string(r.status),
        vmin_pu = r.vmin,
        vmax_pu = r.vmax,
        n_under = r.n_under,
        n_over  = r.n_over
    ))
end

dfA = DataFrame(rowsA)
CSV.write(joinpath(TBLDIR, "pv_size_sweep_far_pf.csv"), dfA)

pA1 = plot(dfA.pv_kw, dfA.vmin_pu, marker=:circle,
           xlabel="PV size (kW, 1ϕ)", ylabel="Minimum bus voltage (pu)",
           title="PV size sweep at far-end (PF) | $NET | alpha=$LOAD_ALPHA",
           label="vmin")
hline!([VMIN_PU], linestyle=:dash, color=:red, label="vmin limit")
savefig(pA1, joinpath(FIGDIR, "pv_sweep_far_vmin.png"))

pA2 = plot(dfA.pv_kw, dfA.n_over, marker=:circle,
           xlabel="PV size (kW, 1ϕ)", ylabel="Over-voltage buses (count)",
           title="PV size sweep at far-end (PF) | $NET | alpha=$LOAD_ALPHA",
           label="n_over")
savefig(pA2, joinpath(FIGDIR, "pv_sweep_far_over_count.png"))

# ============================================================
# Part B: Clustered vs Distributed (same PV size)
# ============================================================
function run_two_pv_case(case_name::String; pv_kw::Float64)
    eng = deepcopy(eng0)

    if case_name == "clustered"
        # PV1 + PV2 both on far load (subtract twice)
        apply_pv_to_load!(eng, far.load_id; pv_kw=pv_kw)
        apply_pv_to_load!(eng, far.load_id; pv_kw=pv_kw)
        pv1_bus = far.bus
        pv2_bus = far.bus
    elseif case_name == "distributed"
        # PV1 on far, PV2 on mid
        apply_pv_to_load!(eng, far.load_id; pv_kw=pv_kw)
        apply_pv_to_load!(eng, mid.load_id; pv_kw=pv_kw)
        pv1_bus = far.bus
        pv2_bus = mid.bus
    else
        error("Unknown case: $case_name")
    end

    r = run_pf_metrics(eng)

    return (
        net = NET,
        case = case_name,
        pv_kw_each = pv_kw,
        pv1_bus = pv1_bus,
        pv2_bus = pv2_bus,
        pv1_load_id = (case_name=="distributed" ? far.load_id : far.load_id),
        pv2_load_id = (case_name=="distributed" ? mid.load_id : far.load_id),
        load_alpha = LOAD_ALPHA,
        status = string(r.status),
        vmin_pu = r.vmin,
        vmax_pu = r.vmax,
        n_under = r.n_under,
        n_over  = r.n_over
    )
end

row_clustered    = run_two_pv_case("clustered"; pv_kw=PV2_KW)
row_distributed  = run_two_pv_case("distributed"; pv_kw=PV2_KW)

dfB = DataFrame([row_clustered, row_distributed])
CSV.write(joinpath(TBLDIR, "pv_clustered_vs_distributed_pf.csv"), dfB)

pB = bar(dfB.case, dfB.vmax_pu,
         xlabel="Case", ylabel="Maximum bus voltage (pu)",
         title="Clustered vs Distributed (PF) | $NET | PV=$(PV2_KW) kW each | alpha=$LOAD_ALPHA",
         label="vmax")
hline!([VMAX_PU], linestyle=:dash, color=:red, label="vmax limit")
savefig(pB, joinpath(FIGDIR, "pv_clustered_vs_distributed_vmax.png"))

println("Saved to: ", OUTDIR)
println("Done.")
