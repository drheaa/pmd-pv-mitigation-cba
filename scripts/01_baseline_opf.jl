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
# Baseline snapshot PF (no PV, no STATCOM) + topology distances
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
VMIN_V = 0.94 * VBASE_LN
VMAX_V = 1.10 * VBASE_LN

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

        # Skip empty
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

In this dataset, the transformer is the only "bridge" from sourcebusz to the LV bus
(e.g. buses=(sourcebusz 9089776)). If we ignore transformers, BFS from source dies.
We assign a tiny length (0.0 km) because the file doesn't give physical distance.
"""
function make_edges_from_transformers(eng)::DataFrame
    rows = NamedTuple[]
    if !haskey(eng, "transformer")
        return DataFrame(rows)
    end

    for (id, tx) in eng["transformer"]
        # PMD transformer representation typically stores buses in tx["bus"] or tx["buses"]
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
            length_km = 0.0,        # transformer doesn't have a line length in DSS here
            kind = "transformer"
        ))
    end

    return DataFrame(rows)
end

"""
Build adjacency and compute shortest-path distances from source_bus.
Because this is radial-ish, BFS works fine, but we use a simple queue and sum lengths.
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

function plot_voltage_along_feeder(buses, edges_df; vmin=0.94*230, vmax=1.10*230)
    p = plot(
        legend=false,
        xlabel="Distance from source (km)",
        ylabel="Voltage (V)",
        title="Voltage along feeder (snapshot)"
    )
    colors = Dict(1=>:blue, 2=>:red, 3=>:black)

    # Only draw physical line segments (not transformers)
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

    if drawn == 0
        println("WARNING: no segments drawn in voltage-vs-distance plot.")
        println("This usually means distances still didn't attach, or bus keys mismatch.")
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
# 1) Parse engineering model
# -----------------------------
println("Parsing: ", master_dss)
eng = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

println("eng counts: buses=", count_dict(eng,"bus"),
        " lines=", count_dict(eng,"line"),
        " loads=", count_dict(eng,"load"),
        " transformers=", count_dict(eng,"transformer"))

eng["settings"]["sbase_default"] = 1
eng["voltage_source"]["source"]["rs"] .= 0
eng["voltage_source"]["source"]["xs"] .= 0

SOURCE_BUS = pick_source_bus_eng(eng)
println("Chosen SOURCE_BUS for distances: ", SOURCE_BUS)

# -----------------------------
# 2) Transform to math model
# -----------------------------
math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)

println("math counts: buses=", count_dict(math,"bus"),
        " loads=", count_dict(math,"load"))

# -----------------------------
# 3) Solve PF snapshot
# -----------------------------
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 1, "sb" => "yes")

println("Solving unbalanced PF (IVRUPowerModel)...")
pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)

println("PF status: ", pf["termination_status"])
println("solution bus count: ", haskey(pf["solution"], "bus") ? length(pf["solution"]["bus"]) : 0)

# -----------------------------
# 4) Build edge list (lines + transformers) and compute distances
# -----------------------------
edges_lines = make_edges_from_lines(eng)
edges_tx    = make_edges_from_transformers(eng)

edges = vcat(edges_lines, edges_tx)

println("edge counts: lines=", nrow(edges_lines), " transformers=", nrow(edges_tx), " total=", nrow(edges))

dist = compute_bus_distances(edges; source_bus=lowercase(SOURCE_BUS))
println("Reachable buses from ", lowercase(SOURCE_BUS), ": ", length(dist))

# -----------------------------
# 5) Extract voltages (keyed by eng names) and attach distances
# -----------------------------
buses = solved_bus_vm_volts_keyed_by_eng(pf, math; vbase_ln=VBASE_LN)
println("extracted bus count (keyed by eng names): ", length(buses))

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

savefig(p1, joinpath(FIGDIR, "baseline_voltage_along.png"))
savefig(p2, joinpath(FIGDIR, "baseline_voltage_hist.png"))
savefig(p,  joinpath(FIGDIR, "baseline_voltage_combined.png"))

println("Saved figures to: ", FIGDIR)

# -----------------------------
# 7) CSV table
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
CSV.write(joinpath(TBLDIR, "baseline_bus_voltages_sorted_by_vmin.csv"), df)

println("Saved table to: ", TBLDIR)
println("Done.")
    