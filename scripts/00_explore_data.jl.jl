import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PowerModelsDistribution
using DataFrames
using CSV
using Plots
using Parquet

const PMD = PowerModelsDistribution
PMD.silence!()

# -----------------------------
# 0) Paths
# -----------------------------
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"

NETWORK = "spd_r"  # spd_r spd_s spd_u spm_r spm_s spm_u

master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NETWORK, "master_scaled.dss")
xy_csv     = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NETWORK, "opendss_xy_$(NETWORK)_scaled.csv")

YEAR = 2023
parquet_file = joinpath(ROOT, "data/raw/dsuite_load_profiles/profiles", "year=$(YEAR)", "49661634_0.parquet")

# outputs
figdir = joinpath(ROOT, "results", "data_analysis")
mkpath(figdir)

# -----------------------------
# 1) Parse network (engineering model)
# -----------------------------
println("\n--- Parsing network ---")
println("Master DSS: ", master_dss)

eng = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!, PMD.reduce_lines!])

count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0

println("\n--- Component counts (engineering model) ---")
println("buses         : ", count_dict(eng, "bus"))
println("lines         : ", count_dict(eng, "line"))
println("branches      : ", count_dict(eng, "branch"))
println("loads         : ", count_dict(eng, "load"))
println("generators    : ", count_dict(eng, "generator"))

println("transformers  : ", count_dict(eng, "transformer"))
println("switches      : ", count_dict(eng, "switch"))
println("shunts        : ", count_dict(eng, "shunt"))
println("capacitors    : ", count_dict(eng, "capacitor"))

# -----------------------------
# 2) Device presence (best-effort scan by name)
# -----------------------------
function find_named(eng, component_key, patterns::Vector{String})
    if !haskey(eng, component_key)
        return String[]
    end
    hits = String[]
    for (id, comp) in eng[component_key]
        nm = lowercase(get(comp, "name", string(id)))
        if any(p -> occursin(p, nm), patterns)
            push!(hits, "$(id):$(get(comp, "name", ""))")
        end
    end
    return hits
end

pv_hits      = find_named(eng, "generator", ["pv", "solar"])
statcom_hits = find_named(eng, "generator", ["statcom", "dstatcom"])
sop_hits     = find_named(eng, "branch",    ["sop", "soft", "openpoint"])

println("\n--- Device presence (best-effort name scan) ---")
println("PV-like generators found      : ", isempty(pv_hits) ? "none" : pv_hits)
println("STATCOM-like generators found : ", isempty(statcom_hits) ? "none" : statcom_hits)
println("SOP-like branches found       : ", isempty(sop_hits) ? "none" : sop_hits)

# -----------------------------
# 3) Load → bus mapping (first 10)
# -----------------------------
println("--- Load → bus mapping (first 10 loads) ---")
if haskey(eng, "load")
    i = 0
    for (id, ld) in eng["load"]
        i += 1
        bus = get(ld, "bus", missing)
        conns = get(ld, "connections", missing)
        pd = get(ld, "pd", missing)
        qd = get(ld, "qd", missing)
        println("load ", id, "  bus=", bus, "  phases=", conns, "  pd=", pd, "  qd=", qd)
        i >= 10 && break
    end
else
    println("No loads found.")
end

# -----------------------------
# 4) Topology plot using XY + lines
#    Your XY file has NO header and is: bus,x,y
# -----------------------------
println("\n--- Plotting topology (xy + lines) ---")
println("XY CSV: ", xy_csv)

if isfile(xy_csv) && haskey(eng, "bus") && haskey(eng, "line")

    # Read CSV with NO header
    xy = CSV.read(xy_csv, DataFrame; header=false)
    if ncol(xy) < 3
        error("XY CSV has < 3 columns. Found columns: $(names(xy))")
    end

    # Rename first 3 columns to bus/x/y
    rename!(xy, Dict(names(xy)[1] => :bus, names(xy)[2] => :x, names(xy)[3] => :y))

    # bus should be string keys
    xy.bus = string.(xy.bus)
    xy.x   = Float64.(xy.x)
    xy.y   = Float64.(xy.y)

    println("XY rows: ", nrow(xy), " | columns: ", names(xy)[1:3])

    # Map bus_key -> (x,y)
    bus_xy = Dict{String,Tuple{Float64,Float64}}()
    for r in eachrow(xy)
        bus_xy[r.bus] = (r.x, r.y)
    end

    # Build mapping: bus_id <-> bus_name
    id_to_name = Dict{String,String}()
    name_to_id = Dict{String,String}()

    for (bid, b) in eng["bus"]
        bid_s = string(bid)
        nm_s  = string(get(b, "name", bid))
        id_to_name[bid_s] = nm_s
        name_to_id[nm_s]  = bid_s
    end

    # Debug samples
    println("Sample XY bus key        : ", first(keys(bus_xy)))
    bid0, b0 = first(collect(eng["bus"]))
    println("Sample eng bus id/name   : ", string(bid0), " / ", string(get(b0, "name", bid0)))
    lid0, ln0 = first(collect(eng["line"]))
    println("Sample line f_bus/t_bus  : ", string(get(ln0, "f_bus", "")), " / ", string(get(ln0, "t_bus", "")))

    # Resolve bus key to match bus_xy keys
    # tries: raw -> id->name -> name->id
    function resolve_xy_key(raw::String)
        if haskey(bus_xy, raw)
            return raw
        end
        if haskey(id_to_name, raw) && haskey(bus_xy, id_to_name[raw])
            return id_to_name[raw]
        end
        if haskey(name_to_id, raw) && haskey(bus_xy, name_to_id[raw])
            return name_to_id[raw]
        end
        return nothing
    end

    # Collect bus points that exist in xy
    xs = Float64[]
    ys = Float64[]
    for (bid, b) in eng["bus"]
        raw = string(get(b, "name", bid))  # prefer name
        k = resolve_xy_key(raw)
        if k !== nothing
            x,y = bus_xy[k]
            push!(xs, x); push!(ys, y)
        else
            # try by id too
            k2 = resolve_xy_key(string(bid))
            if k2 !== nothing
                x,y = bus_xy[k2]
                push!(xs, x); push!(ys, y)
            end
        end
    end

    p = scatter(xs, ys; markersize=2, legend=false,
                title="Topology: $(NETWORK)", xlabel="x", ylabel="y")

    # Draw lines
    drawn = 0
    skipped = 0

    for (lid, ln) in eng["line"]
        f_raw = string(get(ln, "f_bus", ""))
        t_raw = string(get(ln, "t_bus", ""))

        f_key = resolve_xy_key(f_raw)
        t_key = resolve_xy_key(t_raw)

        # If f_bus/t_bus are numeric IDs but stored as names in eng, try mapping
        if f_key === nothing
            f_key = resolve_xy_key(get(id_to_name, f_raw, f_raw))
        end
        if t_key === nothing
            t_key = resolve_xy_key(get(id_to_name, t_raw, t_raw))
        end

        if f_key !== nothing && t_key !== nothing
            x1,y1 = bus_xy[f_key]
            x2,y2 = bus_xy[t_key]
            plot!(p, [x1,x2], [y1,y2]; linewidth=1, label=false)
            drawn += 1
        else
            skipped += 1
        end
    end

    display(p)
    savepath = joinpath(figdir, "topology_$(NETWORK).png")
    savefig(p, savepath)

    println("Lines drawn: ", drawn, " | skipped: ", skipped)
    println("Saved: ", savepath)

else
    println("Skipping topology plot.")
    println("xy exists? ", isfile(xy_csv), " | has bus? ", haskey(eng, "bus"), " | has line? ", haskey(eng, "line"))
end

# -----------------------------
# 5) Read load profile parquet with DuckDB (robust)
# -----------------------------
println("\n--- Reading load profile parquet (DuckDB) ---")
println("Parquet: ", parquet_file)

if isfile(parquet_file)
    using DuckDB
    using DBInterface

    con = DBInterface.connect(DuckDB.DB, ":memory:")

    # Row count (fast)
    nrows_df = DBInterface.execute(con,
        "SELECT COUNT(*) AS n FROM read_parquet('$(parquet_file)')"
    ) |> DataFrame
    nrows = nrows_df.n[1]
    println("Row count: ", nrows)

    # Schema (column names + types)
    schema_df = DBInterface.execute(con,
        "DESCRIBE SELECT * FROM read_parquet('$(parquet_file)')"
    ) |> DataFrame

    println("Column count: ", nrow(schema_df))
    println("\nFirst 30 columns (name :: type):")
    for r in eachrow(schema_df[1:min(30, nrow(schema_df)), :])
        println("  ", r.column_name, " :: ", r.column_type)
    end

    # Pick first numeric column to plot (DuckDB numeric types)
    numeric_mask = occursin.(r"(INTEGER|BIGINT|DOUBLE|FLOAT|REAL|DECIMAL)", schema_df.column_type)
    numcols = schema_df.column_name[numeric_mask]

    if isempty(numcols)
        println("\nNo numeric columns detected from schema types.")
    else
        col = numcols[1]
        println("\nPlotting sample numeric column: ", col)

        # Pull first ~1 week if the table is 30-min resolution (336 points), else first 500
        n = min(nrows, 500)
        sample_df = DBInterface.execute(con,
            "SELECT \"$(col)\" AS y FROM read_parquet('$(parquet_file)') LIMIT $(n)"
        ) |> DataFrame

        # convert to Float64 where possible
        y = Float64.(sample_df.y)

        pprof = plot(1:length(y), y; legend=false, xlabel="row index", ylabel=col,
                     title="Parquet sample: $(col) (first $(length(y)) rows)")
        display(pprof)

        savepath = joinpath(figdir, "profile_sample_year$(YEAR)_$(col).png")
        savefig(pprof, savepath)
        println("Saved: ", savepath)
    end

    DBInterface.close!(con)
else
    println("Parquet file not found at: ", parquet_file)
end


println("\nDone.")
