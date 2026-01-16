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
# Baseline snapshot PF (no PV, no STATCOM) + topology distances
# with LOAD SCALING
# ============================================================

# -----------------------------
# 0) Choose network
# -----------------------------
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"
NET  = "spd_r"  # spd_s, spd_u, spm_r, spm_s, spm_u

master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_scaled.dss")

OUTDIR = joinpath(ROOT, "results", "baseline", NET)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

VBASE_LN = 230.0
VMIN_PU  = 0.94
VMAX_PU  = 1.10
VMIN_V   = VMIN_PU * VBASE_LN
VMAX_V   = VMAX_PU * VBASE_LN

# -----------------------------
# 0.5) Load scaling controls
# -----------------------------
# Option A: fixed alpha (simple, start here)
LOAD_ALPHA = 1.0

# Option B: sweep alphas and pick one automatically
DO_ALPHA_SWEEP = false
ALPHAS = 0.5:0.1:3.0

# "Target" for what the baseline to look like.
# This is a loose heuristic, not a law of physics.
TARGET_Q05_PU = 0.98   # 5th percentile of min-phase voltage, in pu
TARGET_MIN_PU = 0.96   # absolute minimum voltage in pu
PICK_POLICY = :closest_q05  # :closest_q05 or :closest_min

# -----------------------------
# Helpers
# -----------------------------
count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0

"""
Pick a sensible starting bus for BFS in the engineering namespace.
For this dataset, transformer uses sourcebusZ, and PMD often normalizes to lowercase.
"""
function pick_source_bus_eng(eng)::String
    if haskey(eng, "bus")
        haskey(eng["bus"], "sourcebusz") && return "sourcebusz"
        haskey(eng["bus"], "sourcebusZ") && return "sourcebusZ"
        haskey(eng["bus"], "sourcebus")  && return "sourcebus"
    end
    return first(keys(eng["bus"]))
end

"""
Extract an "engineering-like" name for a bus from the math model bus data.
We need eng-like names because eng["line"] uses those bus labels (e.g. 9089776, sourcebusz).
"""
function bus_eng_name(bus_id::String, bus_data::Dict{String,Any})::String
    if haskey(bus_data, "name")
        return lowercase(string(bus_data["name"]))
    end

    if haskey(bus_data, "source_id")
        sid = lowercase(string(bus_data["source_id"]))
        if occursin("bus.", sid)
            parts = split(sid, "bus.")
            if length(parts) >= 2
                return strip(parts[end])
            end
        end
    end

    return lowercase(bus_id)
end

"""
Create edges from eng["line"].
Each edge has (Bus1, Bus2, length_km, kind, id).
"""
function make_edges_from_lines(eng)::DataFrame
    rows = NamedTuple[]
    if !haskey(eng, "line")
        return DataFrame(rows)
    end

    for (id, ln) in eng["line"]
        b1 = lowercase(string(get(ln, "f_bus", "")))
        b2 = lowercase(string(get(ln, "t_bus", "")))
        len_km = (get(ln, "length", 0.0)) / 1000.0

        (isempty(b1) || isempty(b2)) && continue

        push!(rows, (
            edge_id = string(id),
            Bus1 = b1,
            Bus2 = b2,
            length_km = len_km,
            kind = "line"
        ))
    end

    return DataFrame(rows)
end

"""
Create edges from eng["transformer"].

In this dataset, transformer is the "bridge" from sourcebusz to LV bus.
We assign 0 km because physical distance is not given here.
"""
function make_edges_from_transformers(eng)::DataFrame
    rows = NamedTuple[]
    if !haskey(eng, "transformer")
        return DataFrame(rows)
    end

    for (id, tx) in eng["transformer"]
        buses =
            haskey(tx, "bus")   ? tx["bus"] :
            haskey(tx, "buses") ? tx["buses"] :
            nothing

        buses === nothing && continue
        length(buses) < 2 && continue

        b1 = lowercase(string(buses[1]))
        b2 = lowercase(string(buses[2]))

        push!(rows, (
            edge_id = string(id),
            Bus1 = b1,
            Bus2 = b2,
            length_km = 0.0,
            kind = "transformer"
        ))
    end

    return DataFrame(rows)
end

"""
Build adjacency and compute shortest-path distances from source_bus.
This is BFS style, summing line lengths.
"""
function compute_bus_distances(edges_df::DataFrame; source_bus::String)
    adj = Dict{String, Vector{Tuple{String, Float64}}}()

    for r in eachrow(edges_df)
        b1 = r.Bus1
        b2 = r.Bus2
        w  = r.length_km
        push!(get!(adj, b1, Tuple{String,Float64}[]), (b2, w))
        push!(get!(adj, b2, Tuple{String,Float64}[]), (b1, w))
    end

    dist = Dict{String,Float64}(source_bus => 0.0)
    queue = [source_bus]

    while !isempty(queue)
        u = popfirst!(queue)
        for (v, w) in get(adj, u, [])
            if !haskey(dist, v)
                dist[v] = dist[u] + w
                push!(queue, v)
            end
        end
    end

    return dist
end

"""
Extract solved bus voltages in VOLTS keyed by engineering bus name (lowercased).
"""
function solved_bus_vm_volts_keyed_by_eng(pf, math; vbase_ln=230.0)
    sol_bus = pf["solution"]["bus"]
    out = Dict{String, Dict{String,Any}}()

    for (bus_id_any, bus_data_any) in math["bus"]
        bus_id = string(bus_id_any)
        haskey(sol_bus, bus_id) || continue

        sb = sol_bus[bus_id]

        vm_pu =
            haskey(sb, "vm") ? sb["vm"] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"].^2 .+ sb["vi"].^2) :
            nothing
        vm_pu === nothing && continue

        vmV = vm_pu .* vbase_ln

        bus_data = bus_data_any::Dict{String,Any}
        eng_name = bus_eng_name(bus_id, bus_data)

        out[eng_name] = Dict(
            "vma" => vmV[1],
            "vmb" => vmV[2],
            "vmc" => vmV[3]
        )
    end

    return out
end

"""
Plot voltage along feeder, using bus "distance" and bus phase voltages.
"""
function plot_voltage_along_feeder(buses, edges_df; vmin=0.94*230, vmax=1.10*230)
    p = plot(
        legend=false,
        xlabel="Distance from source (km)",
        ylabel="Voltage (V)",
        title="Voltage along feeder (snapshot)"
    )

    colors = Dict(1=>:blue, 2=>:red, 3=>:black)

    line_edges = edges_df[edges_df.kind .== "line", :]

    drawn = 0
    for r in eachrow(line_edges)
        b1, b2 = r.Bus1, r.Bus2
        if !(haskey(buses, b1) && haskey(buses, b2))
            continue
        end
        if !(haskey(buses[b1], "distance") && haskey(buses[b2], "distance"))
            continue
        end

        for ph in (1,2,3)
            v1 = ph==1 ? buses[b1]["vma"] : ph==2 ? buses[b1]["vmb"] : buses[b1]["vmc"]
            v2 = ph==1 ? buses[b2]["vma"] : ph==2 ? buses[b2]["vmb"] : buses[b2]["vmc"]

            plot!(
                p,
                [buses[b1]["distance"], buses[b2]["distance"]],
                [v1, v2],
                color=colors[ph],
                linewidth=1,
                marker=:circle,
                markersize=1
            )
            drawn += 1
        end
    end

    maxdist = maximum(get(b, "distance", 0.0) for b in values(buses))
    plot!(p, [0,maxdist], [vmin,vmin], linestyle=:dash, color=:red)
    plot!(p, [0,maxdist], [vmax,vmax], linestyle=:dash, color=:red)

    return p
end

function plot_voltage_histogram(buses; vmin=0.94*230, vmax=1.10*230)
    va = [b["vma"] for b in values(buses)]
    vb = [b["vmb"] for b in values(buses)]
    vc = [b["vmc"] for b in values(buses)]
    bins = (vmin-2):0.5:(vmax+2)

    p = histogram(va; bins, label="A")
    histogram!(p, vb; bins, label="B")
    histogram!(p, vc; bins, label="C")
    xlabel!(p, "Voltage (V)")
    ylabel!(p, "Count")
    title!(p, "Voltage histogram (snapshot)")
    return p
end

# -----------------------------
# Load scaling helpers
# -----------------------------
"""
Scale all loads in the engineering model by alpha.
Multiplies both active (pd) and reactive (qd).
Handles scalars or vectors.
"""
function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng

    for (id, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}

        for key in ("pd", "qd")
            if haskey(ld, key)
                v = ld[key]
                if v isa Number
                    ld[key] = alpha * v
                elseif v isa AbstractVector
                    ld[key] = alpha .* v
                end
            end
        end
    end

    return eng
end

"""
Compute simple voltage summary stats in pu from 'buses' Dict.
We use min of (A,B,C) per bus.
"""
function voltage_stats_pu(buses::Dict{String,Dict{String,Any}}; vbase_ln=230.0)
    vmins_pu = Float64[]
    for b in values(buses)
        vminV = min(b["vma"], b["vmb"], b["vmc"])
        push!(vmins_pu, vminV / vbase_ln)
    end

    # Guard for empty
    if isempty(vmins_pu)
        return (min=NaN, q05=NaN, median=NaN, q95=NaN)
    end

    sort!(vmins_pu)
    q(p) = vmins_pu[clamp(Int(ceil(p * length(vmins_pu))), 1, length(vmins_pu))]

    return (
        min = minimum(vmins_pu),
        q05 = q(0.05),
        median = q(0.50),
        q95 = q(0.95)
    )
end

# -----------------------------
# 1) Parse engineering model
# -----------------------------
println("Parsing: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

println("eng counts: buses=", count_dict(eng0,"bus"),
        " lines=", count_dict(eng0,"line"),
        " loads=", count_dict(eng0,"load"),
        " transformers=", count_dict(eng0,"transformer"))

eng0["voltage_source"]["source"]["rs"] .= 0
eng0["voltage_source"]["source"]["xs"] .= 0

SOURCE_BUS = pick_source_bus_eng(eng0)
println("Chosen SOURCE_BUS for distances: ", SOURCE_BUS)

# -----------------------------
# 2) Build edge list and distances (from eng0, not scaled)
# -----------------------------
edges_lines = make_edges_from_lines(eng0)
edges_tx    = make_edges_from_transformers(eng0)
edges = vcat(edges_lines, edges_tx)

println("edge counts: lines=", nrow(edges_lines), " transformers=", nrow(edges_tx), " total=", nrow(edges))

dist = compute_bus_distances(edges; source_bus=lowercase(SOURCE_BUS))
println("Reachable buses from ", lowercase(SOURCE_BUS), ": ", length(dist))

# -----------------------------
# 3) Choose alpha (fixed or sweep)
# -----------------------------
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

function run_pf_with_alpha(eng_template, alpha)
    eng = deepcopy(eng_template)
    scale_loads!(eng, alpha)

    math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)
    pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)

    buses = solved_bus_vm_volts_keyed_by_eng(pf, math; vbase_ln=VBASE_LN)
    stats = voltage_stats_pu(buses; vbase_ln=VBASE_LN)

    return (pf=pf, math=math, buses=buses, stats=stats)
end

chosen_alpha = LOAD_ALPHA
sweep_df = DataFrame()

if DO_ALPHA_SWEEP
    rows = NamedTuple[]
    best_alpha = nothing
    best_score = Inf

    for a in ALPHAS
        res = run_pf_with_alpha(eng0, a)
        st = res.stats

        # Skip failed cases
        if isnan(st.min) || isnan(st.q05)
            continue
        end

        score =
            PICK_POLICY == :closest_min ? abs(st.min - TARGET_MIN_PU) :
            abs(st.q05 - TARGET_Q05_PU)

        push!(rows, (alpha=a, vmin_pu=st.min, vq05_pu=st.q05, vmed_pu=st.median, vq95_pu=st.q95, score=score))

        if score < best_score
            best_score = score
            best_alpha = a
        end
    end

    sweep_df = DataFrame(rows)
    sort!(sweep_df, :score)

    if nrow(sweep_df) == 0
        error("Alpha sweep produced no valid results. Check solve status / data.")
    end

    chosen_alpha = best_alpha
    println("Chosen alpha from sweep = ", chosen_alpha)

    CSV.write(joinpath(TBLDIR, "baseline_alpha_sweep.csv"), sweep_df)
end

println("Using LOAD_ALPHA = ", chosen_alpha)

# -----------------------------
# 4) Run PF with chosen alpha
# -----------------------------
println("Solving unbalanced PF (IVRUPowerModel) with load scaling...")
res = run_pf_with_alpha(eng0, chosen_alpha)
pf   = res.pf
math = res.math
buses = res.buses

println("PF status: ", pf["termination_status"])
println("solution bus count: ", haskey(pf["solution"], "bus") ? length(pf["solution"]["bus"]) : 0)

st = res.stats
println("Voltage stats (min-phase per bus, pu): min=", round(st.min, digits=4),
        " q05=", round(st.q05, digits=4),
        " median=", round(st.median, digits=4),
        " q95=", round(st.q95, digits=4))

# -----------------------------
# 5) Attach distances
# -----------------------------
matched = 0
for (bus, d) in dist
    if haskey(buses, bus)
        buses[bus]["distance"] = d
        matched += 1
    end
end
println("distance matched for ", matched, " buses")

# -----------------------------
# 6) Plots
# -----------------------------
p1 = plot_voltage_along_feeder(buses, edges; vmin=VMIN_V, vmax=VMAX_V)
p2 = plot_voltage_histogram(buses; vmin=VMIN_V, vmax=VMAX_V)
p  = plot(p1, p2, layout=(1,2), size=(1100,450))

savefig(p1, joinpath(FIGDIR, "baseline_voltage_along_alpha_$(chosen_alpha).png"))
savefig(p2, joinpath(FIGDIR, "baseline_voltage_hist_alpha_$(chosen_alpha).png"))
savefig(p,  joinpath(FIGDIR, "baseline_voltage_combined_alpha_$(chosen_alpha).png"))

println("Saved figures to: ", FIGDIR)

# -----------------------------
# 7) CSV tables
# -----------------------------
rows = NamedTuple[]
for (bus, b) in buses
    push!(rows, (
        bus = bus,
        distance_km = get(b, "distance", missing),
        vA = b["vma"],
        vB = b["vmb"],
        vC = b["vmc"],
        vmin = min(b["vma"], b["vmb"], b["vmc"]),
        vmax = max(b["vma"], b["vmb"], b["vmc"])
    ))
end

df = DataFrame(rows)
sort!(df, :vmin)
CSV.write(joinpath(TBLDIR, "baseline_bus_voltages_sorted_by_vmin_alpha_$(chosen_alpha).csv"), df)

# Save the chosen alpha and headline stats
meta = DataFrame([
    (net=NET, alpha=chosen_alpha, vmin_pu=st.min, vq05_pu=st.q05, vmed_pu=st.median, vq95_pu=st.q95)
])
CSV.write(joinpath(TBLDIR, "baseline_alpha_chosen.csv"), meta)

println("Saved tables to: ", TBLDIR)
println("Done.")
