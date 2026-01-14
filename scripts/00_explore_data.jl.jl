import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using PowerModelsDistribution
using DataFrames
using CSV
using Plots
using DuckDB
using DBInterface
using Statistics

const PMD = PowerModelsDistribution
PMD.silence!()

# -----------------------------
# 0) Paths + configuration
# -----------------------------
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"

NETWORKS = ["spd_r", "spd_s", "spd_u", "spm_r", "spm_s", "spm_u"]

PROFILES_ROOT = joinpath(ROOT, "data/raw/dsuite_load_profiles/profiles")
NET_ROOT      = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1")

OUTDIR = joinpath(ROOT, "results", "data_analysis")
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR)
mkpath(TBLDIR)

MAKE_TOPOLOGY_PLOTS = true
MAKE_PROFILE_SAMPLE_PLOT = true

# -----------------------------
# Helpers
# -----------------------------
count_dict(d, key) = haskey(d, key) ? length(d[key]) : 0

function find_named(eng, component_key, patterns::Vector{String})
    if !haskey(eng, component_key)
        return String[]
    end
    hits = String[]
    for (id, comp) in eng[component_key]
        nm = lowercase(string(get(comp, "name", id)))
        if any(p -> occursin(p, nm), patterns)
            push!(hits, string(get(comp, "name", id)))
        end
    end
    sort!(unique(hits))
    return hits
end

# resolve XY bus key differences (id vs name)
function build_bus_maps(eng)
    id_to_name = Dict{String,String}()
    name_to_id = Dict{String,String}()
    for (bid, b) in eng["bus"]
        bid_s = string(bid)
        nm_s  = string(get(b, "name", bid))
        id_to_name[bid_s] = nm_s
        name_to_id[nm_s]  = bid_s
    end
    return id_to_name, name_to_id
end

function resolve_xy_key(raw::String, bus_xy, id_to_name, name_to_id)
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

# -----------------------------
# 1) Loop over networks: parse + counts + topology plot + CSV
# -----------------------------
network_rows = NamedTuple[]

for net in NETWORKS
    println("\n==============================")
    println("NETWORK: ", net)
    println("==============================")

    master_dss = joinpath(NET_ROOT, net, "master_scaled.dss")
    xy_csv     = joinpath(NET_ROOT, net, "opendss_xy_$(net)_scaled.csv")

    if !isfile(master_dss)
        println("Missing master DSS: ", master_dss)
        continue
    end

    println("Parsing: ", master_dss)
    eng = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!, PMD.reduce_lines!])

    nb = count_dict(eng, "bus")
    nl = count_dict(eng, "line")
    nbr = count_dict(eng, "branch")
    nld = count_dict(eng, "load")
    ngen = count_dict(eng, "generator")
    ntx = count_dict(eng, "transformer")

    pv_hits      = find_named(eng, "generator", ["pv", "solar"])
    statcom_hits = find_named(eng, "generator", ["statcom", "dstatcom"])
    sop_hits     = find_named(eng, "branch",    ["sop", "soft", "openpoint"])

    push!(network_rows, (
        network=net,
        master_dss=master_dss,
        xy_csv=xy_csv,
        buses=nb,
        lines=nl,
        branches=nbr,
        loads=nld,
        generators=ngen,
        transformers=ntx,
        pv_like_generators=join(pv_hits, ";"),
        statcom_like_generators=join(statcom_hits, ";"),
        sop_like_branches=join(sop_hits, ";")
    ))

    # ---- Topology plot ----
    if MAKE_TOPOLOGY_PLOTS && isfile(xy_csv) && haskey(eng, "bus") && haskey(eng, "line")
        println("Plotting topology...")

        xy = CSV.read(xy_csv, DataFrame; header=false)
        rename!(xy, Dict(names(xy)[1] => :bus, names(xy)[2] => :x, names(xy)[3] => :y))
        xy.bus = string.(xy.bus)
        xy.x   = Float64.(xy.x)
        xy.y   = Float64.(xy.y)

        bus_xy = Dict{String,Tuple{Float64,Float64}}()
        for r in eachrow(xy)
            bus_xy[r.bus] = (r.x, r.y)
        end

        id_to_name, name_to_id = build_bus_maps(eng)

        xs = Float64[]
        ys = Float64[]
        for (bid, b) in eng["bus"]
            raw = string(get(b, "name", bid))
            k = resolve_xy_key(raw, bus_xy, id_to_name, name_to_id)
            if k === nothing
                k = resolve_xy_key(string(bid), bus_xy, id_to_name, name_to_id)
            end
            if k !== nothing
                x,y = bus_xy[k]
                push!(xs, x); push!(ys, y)
            end
        end

        p = scatter(xs, ys; markersize=2, legend=false,
                    title="Topology: $(net)", xlabel="x", ylabel="y")

        drawn = 0
        skipped = 0
        for (lid, ln) in eng["line"]
            f_raw = string(get(ln, "f_bus", ""))
            t_raw = string(get(ln, "t_bus", ""))

            f_key = resolve_xy_key(f_raw, bus_xy, id_to_name, name_to_id)
            t_key = resolve_xy_key(t_raw, bus_xy, id_to_name, name_to_id)

            if f_key === nothing
                f_key = resolve_xy_key(get(id_to_name, f_raw, f_raw), bus_xy, id_to_name, name_to_id)
            end
            if t_key === nothing
                t_key = resolve_xy_key(get(id_to_name, t_raw, t_raw), bus_xy, id_to_name, name_to_id)
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

        savepath = joinpath(FIGDIR, "topology_$(net).png")
        savefig(p, savepath)
        println("Saved: ", savepath, " (lines drawn=", drawn, ", skipped=", skipped, ")")
    else
        println("Skipping topology plot for ", net)
    end
end

network_df = DataFrame(network_rows)
CSV.write(joinpath(TBLDIR, "network_summary.csv"), network_df)
println("\nSaved network summary CSV: ", joinpath(TBLDIR, "network_summary.csv"))

# -----------------------------
# 2) Loop over parquet years: schema + stats + CSV
# -----------------------------
println("\n==============================")
println("PROFILES: scanning ", PROFILES_ROOT)
println("==============================")

year_dirs = filter(d -> isdir(joinpath(PROFILES_ROOT, d)) && occursin("year=", d), readdir(PROFILES_ROOT))
sort!(year_dirs)

profile_rows = NamedTuple[]

for yd in year_dirs
    year = replace(yd, "year=" => "")
    yr_path = joinpath(PROFILES_ROOT, yd)

    parquet_files = filter(f -> endswith(lowercase(f), ".parquet"), readdir(yr_path))
    if isempty(parquet_files)
        println("No parquet in ", yr_path)
        continue
    end

    # dataset seems to have a single parquet per year folder
    for pq in parquet_files
        parquet_file = joinpath(yr_path, pq)
        println("\n--- YEAR ", year, " | ", pq, " ---")

        con = DBInterface.connect(DuckDB.DB, ":memory:")

        nrows_df = DBInterface.execute(con,
            "SELECT COUNT(*) AS n FROM read_parquet('$(parquet_file)')"
        ) |> DataFrame
        nrows = nrows_df.n[1]

        schema_df = DBInterface.execute(con,
            "DESCRIBE SELECT * FROM read_parquet('$(parquet_file)')"
        ) |> DataFrame
        ncols = nrow(schema_df)

        # numeric columns
        numeric_mask = occursin.(r"(INTEGER|BIGINT|DOUBLE|FLOAT|REAL|DECIMAL)", schema_df.column_type)
        numcols = schema_df.column_name[numeric_mask]

        push!(profile_rows, (
            year=parse(Int, year),
            parquet=parquet_file,
            rows=nrows,
            cols=ncols,
            numeric_cols=length(numcols)
        ))

        # ---- compute per-column stats in chunks (avoid loading 3300 cols at once) ----
        # We’ll output a CSV with mean/min/max/std for each numeric column.
        stats_out = NamedTuple[]
        chunk_size = 50  # adjust if you want faster vs memory

        println("Computing stats for ", length(numcols), " numeric columns (chunk_size=", chunk_size, ")")

        for i in 1:chunk_size:length(numcols)
            cols_chunk = numcols[i:min(i+chunk_size-1, length(numcols))]

            select_list = join(["\"$(c)\" AS \"$(c)\"" for c in cols_chunk], ", ")
            q = "SELECT $(select_list) FROM read_parquet('$(parquet_file)')"

            df = DBInterface.execute(con, q) |> DataFrame

            for c in cols_chunk
                v = df[!, Symbol(c)]
                # skip missings if any
                vv = skipmissing(v)
                # DuckDB sometimes returns Union{Missing,Float64} etc
                vals = collect(Float64, vv)

                push!(stats_out, (
                    year=parse(Int, year),
                    column=c,
                    mean=mean(vals),
                    std=std(vals),
                    min=minimum(vals),
                    max=maximum(vals)
                ))
            end
        end

        stats_df = DataFrame(stats_out)
        stats_csv = joinpath(TBLDIR, "profile_stats_year$(year).csv")
        CSV.write(stats_csv, stats_df)
        println("Saved: ", stats_csv)

        # ---- optional sample plot (first numeric column only) ----
        if MAKE_PROFILE_SAMPLE_PLOT && !isempty(numcols)
            col = numcols[1]
            n = min(nrows, 500)
            sample_df = DBInterface.execute(con,
                "SELECT \"$(col)\" AS y FROM read_parquet('$(parquet_file)') LIMIT $(n)"
            ) |> DataFrame
            y = Float64.(sample_df.y)

            pprof = plot(1:length(y), y; legend=false, xlabel="row index", ylabel=col,
                         title="Parquet sample: $(col) (first $(length(y)) rows)")
            savepath = joinpath(FIGDIR, "profile_sample_year$(year)_$(col).png")
            savefig(pprof, savepath)
            println("Saved: ", savepath)
        end

        DBInterface.close!(con)
    end
end

profile_df = DataFrame(profile_rows)
CSV.write(joinpath(TBLDIR, "parquet_summary.csv"), profile_df)
println("\nSaved parquet summary CSV: ", joinpath(TBLDIR, "parquet_summary.csv"))

println("\nDone. Outputs in:")
println("  ", OUTDIR)
