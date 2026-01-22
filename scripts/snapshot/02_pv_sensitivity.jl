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

# ============================================================
# 02) PV sensitivity (OPF snapshot)
# - Fixed 1ϕ PV size
# - Move PV to near/mid/far load buses
# - PV modeled as negative load on an existing customer load
# ============================================================

ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"
NET  = "spd_s"
master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_scaled.dss")

LOAD_ALPHA = 1.0
PV_KW      = 50.0     # kW (single-phase)

VBASE_LN = 230.0
VMIN_PU  = 0.94
VMAX_PU  = 1.10

RUN_TAG = "opf_pv_sensitivity_alpha$(LOAD_ALPHA)_pv$(PV_KW)kW_1ph"
OUTDIR = joinpath(ROOT, "results", "pv_sensitivity_opf", NET, RUN_TAG)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0

# -----------------------------
# name normalization
# -----------------------------
normalize_bus(s::AbstractString) = lowercase(split(String(s), ".")[1])

function pick_source_bus_eng(eng)::String
    haskey(eng["bus"], "sourcebusz") && return "sourcebusz"
    haskey(eng["bus"], "sourcebusZ") && return "sourcebusZ"
    haskey(eng["bus"], "sourcebus")  && return "sourcebus"
    return first(keys(eng["bus"]))
end

# -----------------------------
# distances
# -----------------------------
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

# -----------------------------
# load scaling
# -----------------------------
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

# -----------------------------
# choose PV placement targets 
# -----------------------------
function collect_customer_load_points(eng::Dict{String,Any})
    pts = NamedTuple[]
    haskey(eng, "load") || return pts

    for (lid, ld_any) in eng["load"]
        id = String(lid)
        startswith(id, "load_") || continue 
        ld = ld_any::Dict{String,Any}
        haskey(ld, "bus") || continue

        bus_raw = String(ld["bus"])
        bus_base = normalize_bus(bus_raw)

        # phase: from "12345.2" if present; else 1
        parts = split(bus_raw, ".")
        ph = (length(parts) >= 2 && tryparse(Int, parts[2]) !== nothing) ? clamp(parse(Int, parts[2]),1,3) : 1

        # single-phase only (for this experiment)
        nph = get(ld, "phases", 1)
        (nph isa Number && Int(nph) == 1) || continue

        push!(pts, (load_id=id, bus=bus_base, phase=ph))
    end
    return pts
end

function apply_pv_to_load!(eng::Dict{String,Any}, load_id::String; pv_kw::Float64)
    ld = eng["load"][load_id]::Dict{String,Any}

    # subtract pv_kw from scalar or first element of vector
    function sub_kw!(ld::Dict{String,Any}, key::String)
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

    # D-Suite style
    if haskey(ld, "pd_nom") && sub_kw!(ld, "pd_nom")
        return nothing
    end

    # PMD style (some feeders / conversions use this)
    if haskey(ld, "pd") && sub_kw!(ld, "pd")
        return nothing
    end

    # OpenDSS style (rare here but harmless to try)
    if haskey(ld, "kW") && sub_kw!(ld, "kW")
        return nothing
    end
    if haskey(ld, "kw") && sub_kw!(ld, "kw")
        return nothing
    end

    error("Load $load_id has no active power key (pd_nom/pd/kW/kw). Keys = $(collect(keys(ld)))")
end

# -----------------------------
# OPF runner
# -----------------------------
function run_opf(eng::Dict{String,Any})
    math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)

    for (_, bus) in math["bus"]
        bus["vmin"] = VMIN_PU .* ones(length(bus["vmin"]))
        bus["vmax"] = VMAX_PU .* ones(length(bus["vmax"]))
    end

    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
    res = PMD.solve_mc_opf(math, PMD.IVRUPowerModel, ipopt)

    sol = res["solution"]["bus"]
    vmins = Float64[]
    vmaxs = Float64[]
    for sb in values(sol)
        vm = haskey(sb,"vm") ? sb["vm"] : sqrt.(sb["vr"].^2 .+ sb["vi"].^2)
        push!(vmins, minimum(vm))
        push!(vmaxs, maximum(vm))
    end

    return (
        status = string(res["termination_status"]),                 # <-- FIXED
        objective = Float64(get(res, "objective", NaN)),
        vmin = minimum(vmins),
        vmax = maximum(vmaxs),
        n_under = count(<(VMIN_PU), vmins),
        n_over  = count(>(VMAX_PU), vmaxs)
    )
end

# ============================================================
# main
# ============================================================
println("Parsing: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

println("eng counts: buses=", count_dict(eng0,"bus"),
        " lines=", count_dict(eng0,"line"),
        " loads=", count_dict(eng0,"load"),
        " transformers=", count_dict(eng0,"transformer"))

scale_loads!(eng0, LOAD_ALPHA)

SOURCE_BUS = normalize_bus(pick_source_bus_eng(eng0))

edges = vcat(make_edges_from_lines(eng0), make_edges_from_transformers(eng0))
dist  = compute_bus_distances(edges; source_bus=SOURCE_BUS)

pts = collect_customer_load_points(eng0)
isempty(pts) && error("No customer 1ϕ loads found (load_*).")

# rank by distance
pts2 = [(; p..., d=get(dist, p.bus, Inf)) for p in pts]
sort!(pts2; by = x -> x.d)

near = pts2[1]
mid  = pts2[clamp(Int(round(length(pts2)/2)), 1, length(pts2))]
far  = pts2[end]

println("PV buses (distance-based, customer load buses only):")
println("  near_source: $(near.bus) (d=$(near.d) km) load_id=$(near.load_id)")
println("  mid:         $(mid.bus) (d=$(mid.d) km)  load_id=$(mid.load_id)")
println("  far:         $(far.bus) (d=$(far.d) km)  load_id=$(far.load_id)")

cases = [("near_source", near), ("mid", mid), ("far", far)]

# baseline
base = run_opf(deepcopy(eng0))
println("Baseline OPF status: $(base.status) | objective: $(base.objective)")

rows = NamedTuple[]
push!(rows, (net=NET, case="baseline", pv_bus="(none)", pv_kw=0.0, pv_phase=1,
             load_alpha=LOAD_ALPHA, status=base.status, objective=base.objective,
             vmin_pu=base.vmin, vmax_pu=base.vmax, n_under=base.n_under, n_over=base.n_over))

labels = ["baseline"]
vmins  = [base.vmin]
vmaxs  = [base.vmax]

for (label, p) in cases
    eng = deepcopy(eng0)
    apply_pv_to_load!(eng, p.load_id; pv_kw=PV_KW)
    res = run_opf(eng)

    push!(rows, (net=NET, case=label, pv_bus=p.bus, pv_kw=PV_KW, pv_phase=p.phase,
                 load_alpha=LOAD_ALPHA, status=res.status, objective=res.objective,
                 vmin_pu=res.vmin, vmax_pu=res.vmax, n_under=res.n_under, n_over=res.n_over))

    push!(labels, label)
    push!(vmins, res.vmin)
    push!(vmaxs, res.vmax)
end

df = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "pv_sensitivity_opf_summary.csv"), df)

# plots
p1 = bar(labels, vmins; ylabel="Minimum bus voltage (pu)",
         title="PV sensitivity (OPF snapshot) | $NET | PV=$(PV_KW) kW (1ϕ)")
hline!([VMIN_PU, VMAX_PU], linestyle=:dash, color=:red, label="limits")
savefig(p1, joinpath(FIGDIR, "pv_sensitivity_vmin.png"))

p2 = bar(labels, vmaxs; ylabel="Maximum bus voltage (pu)",
         title="PV sensitivity (OPF snapshot) | $NET | PV=$(PV_KW) kW (1ϕ)")
hline!([VMIN_PU, VMAX_PU], linestyle=:dash, color=:red, label="limits")
savefig(p2, joinpath(FIGDIR, "pv_sensitivity_vmax.png"))

println("Saved results to: ", OUTDIR)
println("Done.")
