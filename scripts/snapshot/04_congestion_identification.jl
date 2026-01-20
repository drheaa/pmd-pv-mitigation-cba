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
# 04) Congestion identification (snapshot PF)
#
# A) Voltage congestion: vmin/vmax + under/over counts + weakest buses
# B) Branch stress: top |S| branches, and loading % if rating exists
#
# Note: Many D-Suite feeders don't include thermal ratings, so we rank by |S|.
# ============================================================

# -----------------------------
# 0) Network + paths
# -----------------------------
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"
NET  = "spd_s"   # change if needed

master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_scaled.dss")

OUTDIR = joinpath(ROOT, "results", "congestion_identification", NET)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

# -----------------------------
# 1) Parameters (edit)
# -----------------------------
LOAD_ALPHA = 1.0

VMIN_PU = 0.94
VMAX_PU = 1.10

TOP_K_BRANCHES = 15

# -----------------------------
# Helpers
# -----------------------------
count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0

normalize_bus(s::AbstractString) = lowercase(split(String(s), ".")[1])

function pick_source_bus_eng(eng)::String
    haskey(eng, "bus") || return ""
    haskey(eng["bus"], "sourcebusz") && return "sourcebusz"
    haskey(eng["bus"], "sourcebusZ") && return "sourcebusZ"
    haskey(eng["bus"], "sourcebus")  && return "sourcebus"
    return first(keys(eng["bus"]))
end

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

function make_edges_from_lines(eng)
    rows = NamedTuple[]
    haskey(eng, "line") || return DataFrame(rows)
    for (id, ln) in eng["line"]
        b1 = normalize_bus(string(get(ln, "f_bus", "")))
        b2 = normalize_bus(string(get(ln, "t_bus", "")))
        (isempty(b1) || isempty(b2)) && continue
        len = Float64(get(ln, "length", 0.0)) / 1000.0
        isfinite(len) || (len = 0.0)
        push!(rows, (Bus1=b1, Bus2=b2, length_km=len, kind="line", edge_id=string(id)))
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
        push!(rows, (Bus1=b1, Bus2=b2, length_km=0.0, kind="transformer", edge_id=string(id)))
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

function run_pf(eng::Dict{String,Any})
    math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
    pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)
    return (pf=pf, math=math)
end

function extract_bus_vm_pu(pf, math)
    sol_bus = pf["solution"]["bus"]
    out = NamedTuple[]

    for (bus_id_any, bus_data_any) in math["bus"]
        bus_id = string(bus_id_any)
        haskey(sol_bus, bus_id) || continue
        sb = sol_bus[bus_id]

        vm =
            haskey(sb, "vm") ? sb["vm"] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"].^2 .+ sb["vi"].^2) :
            nothing
        vm === nothing && continue

        eng_name =
            haskey(bus_data_any, "name") ? normalize_bus(string(bus_data_any["name"])) :
            bus_id

        vmin = minimum(vm[1:min(3,length(vm))])
        vmax = maximum(vm[1:min(3,length(vm))])
        vA = length(vm) >= 1 ? vm[1] : NaN
        vB = length(vm) >= 2 ? vm[2] : NaN
        vC = length(vm) >= 3 ? vm[3] : NaN

        push!(out, (bus=eng_name, vA_pu=vA, vB_pu=vB, vC_pu=vC, vmin_pu=vmin, vmax_pu=vmax))
    end

    return DataFrame(out)
end

# robust: magnitude of apparent power (supports scalar or vector)
function s_mag(p, q)
    if p isa Number && q isa Number
        return sqrt(Float64(p)^2 + Float64(q)^2)
    elseif p isa AbstractVector && q isa AbstractVector
        n = min(length(p), length(q))
        if n == 0
            return NaN
        end
        pv = Float64.(p[1:n])
        qv = Float64.(q[1:n])
        return maximum(sqrt.(pv.^2 .+ qv.^2))
    else
        return NaN
    end
end

# robust: rating might be scalar OR vector OR missing
function rating_to_scalar(rate)
    rate === missing && return missing
    if rate isa Number
        r = Float64(rate)
        return (isfinite(r) && r > 0) ? r : missing
    elseif rate isa AbstractVector
        rr = Float64.(rate)
        rr = rr[isfinite.(rr) .& (rr .> 0)]
        isempty(rr) && return missing
        return maximum(rr)  # choose max rating if per-phase
    else
        return missing
    end
end

function extract_branch_stress(pf, math)
    sol = get(pf["solution"], "branch", Dict{String,Any}())
    haskey(math, "branch") || return DataFrame()

    rows = NamedTuple[]
    for (br_id_any, br_data_any) in math["branch"]
        br_id = string(br_id_any)
        haskey(sol, br_id) || continue

        sb = sol[br_id]
        br = br_data_any::Dict{String,Any}

        f_bus = haskey(br, "f_bus") ? normalize_bus(string(br["f_bus"])) : ""
        t_bus = haskey(br, "t_bus") ? normalize_bus(string(br["t_bus"])) : ""

        pfk = haskey(sb, "pf") ? sb["pf"] :
              haskey(sb, "p_fr") ? sb["p_fr"] :
              haskey(sb, "p") ? sb["p"] : missing

        qfk = haskey(sb, "qf") ? sb["qf"] :
              haskey(sb, "q_fr") ? sb["q_fr"] :
              haskey(sb, "q") ? sb["q"] : missing

        s_from = (pfk === missing || qfk === missing) ? NaN : s_mag(pfk, qfk)

        rate_raw =
            haskey(br, "rate_a") ? br["rate_a"] :
            haskey(br, "rate")   ? br["rate"]   :
            missing

        rate = rating_to_scalar(rate_raw)

        loading_pct =
            (rate === missing || !isfinite(rate) || rate <= 0 || !isfinite(s_from)) ?
            missing :
            100.0 * (s_from / rate)

        push!(rows, (
            branch_id = br_id,
            f_bus = f_bus,
            t_bus = t_bus,
            s_from = s_from,
            rate = rate,
            loading_pct = loading_pct
        ))
    end

    df = DataFrame(rows)
    sort!(df, :s_from, rev=true)
    return df
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
dist = compute_bus_distances(edges; source_bus=SOURCE_BUS)

res = run_pf(eng0)
pf = res.pf
math = res.math

println("PF status: ", pf["termination_status"])

# -----------------------------
# A) voltage congestion
# -----------------------------
bus_df = extract_bus_vm_pu(pf, math)
bus_df.distance_km = [haskey(dist, b) ? dist[b] : missing for b in bus_df.bus]
bus_df.under_vmin = bus_df.vmin_pu .< VMIN_PU
bus_df.over_vmax  = bus_df.vmax_pu .> VMAX_PU

sort!(bus_df, :vmin_pu)
CSV.write(joinpath(TBLDIR, "bus_voltages_sorted_by_vmin.csv"), bus_df)

n_under = count(bus_df.under_vmin)
n_over  = count(bus_df.over_vmax)

println("Voltage violations: under=", n_under, " over=", n_over)
println("Weakest bus: ", bus_df.bus[1], " vmin_pu=", bus_df.vmin_pu[1])

k = min(15, nrow(bus_df))
p_vweak = bar(bus_df.bus[1:k], bus_df.vmin_pu[1:k];
    xlabel="Bus",
    ylabel="Min phase voltage (pu)",
    title="Weakest buses (PF) | $NET | alpha=$LOAD_ALPHA",
    xrotation=45,
    legend=false
)
hline!([VMIN_PU], linestyle=:dash, color=:red)
savefig(p_vweak, joinpath(FIGDIR, "weakest_buses_vmin.png"))

# -----------------------------
# B) branch stress congestion
# -----------------------------
br_df = extract_branch_stress(pf, math)
CSV.write(joinpath(TBLDIR, "branch_stress_ranked.csv"), br_df)

if nrow(br_df) > 0
    k2 = min(TOP_K_BRANCHES, nrow(br_df))
    top = br_df[1:k2, :]
    labels = ["$(top.f_bus[i])->$(top.t_bus[i])" for i in 1:k2]

    p_br = bar(labels, top.s_from;
        xlabel="Branch",
        ylabel="|S_from| (model units)",
        title="Top stressed branches (PF) | $NET | alpha=$LOAD_ALPHA",
        xrotation=45,
        legend=false
    )
    savefig(p_br, joinpath(FIGDIR, "top_branch_stress.png"))

    if any(.!ismissing.(top.loading_pct))
        vals = [ismissing(x) ? 0.0 : Float64(x) for x in top.loading_pct]
        p_pct = bar(labels, vals;
            xlabel="Branch",
            ylabel="Loading (%)",
            title="Branch loading % (if ratings exist) | $NET",
            xrotation=45,
            legend=false
        )
        savefig(p_pct, joinpath(FIGDIR, "top_branch_loading_pct.png"))
    end
else
    println("No branch solution data found in pf[\"solution\"][\"branch\"].")
end

# -----------------------------
# summary
# -----------------------------
summary = DataFrame([
    (net=NET,
     load_alpha=LOAD_ALPHA,
     pf_status=string(pf["termination_status"]),
     n_buses=nrow(bus_df),
     vmin_pu=minimum(bus_df.vmin_pu),
     vmax_pu=maximum(bus_df.vmax_pu),
     n_under=n_under,
     n_over=n_over,
     top_branch_s = (nrow(br_df) > 0 ? br_df.s_from[1] : missing),
     top_branch_id = (nrow(br_df) > 0 ? br_df.branch_id[1] : missing))
])
CSV.write(joinpath(TBLDIR, "congestion_summary.csv"), summary)

println("Saved to: ", OUTDIR)
println("Done.")
