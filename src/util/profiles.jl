"""
find_parquet_for_year(profiles_root, year)

Returns one parquet path under profiles_root/year=YYYY.
This is a stable convention in the D-Suite profile bundle.
"""
function find_parquet_for_year(profiles_root::String, year::Int)
    year_dir = joinpath(profiles_root, "year=$(year)")
    isdir(year_dir) || error("Year folder not found: $year_dir")

    files = filter(f -> endswith(lowercase(f), ".parquet"), readdir(year_dir))
    isempty(files) && error("No parquet files found in: $year_dir")

    return joinpath(year_dir, files[1])
end

"""
numeric_columns(schema_df)

DuckDB DESCRIBE provides column_type strings.
This selects numeric-like columns for profile selection.
"""
function numeric_columns(schema_df::DataFrame)
    mask = occursin.(r"(INTEGER|BIGINT|DOUBLE|FLOAT|REAL|DECIMAL)", schema_df.column_type)
    return Vector{String}(schema_df.column_name[mask])
end

"""
pick_customer_column_auto(con, parquet_file, cols)

Chooses a profile column that has variance and nonzero values.
A small sample is used for speed.
"""
function pick_customer_column_auto(con, parquet_file::String, cols::Vector{String})
    best_col = nothing
    best_score = -Inf

    for c in cols
        if c == "hh" || c == "year"
            continue
        end

        q = "SELECT \"$(c)\" AS y FROM read_parquet('$(parquet_file)') LIMIT 2000"
        df = DBInterface.execute(con, q) |> DataFrame

        vals = Float64[]
        for v in df.y
            v === missing && continue
            try
                push!(vals, Float64(v))
            catch
            end
        end
        isempty(vals) && continue

        nonzero_frac = count(x -> abs(x) > 1e-9, vals) / length(vals)
        score = Statistics.var(vals) * nonzero_frac

        if score > best_score
            best_score = score
            best_col = c
        end
    end

    best_col === nothing && error("No suitable customer column found")
    return best_col
end

"""
build_time_axis(nrows; start_date, step_minutes)

Builds a DateTime axis for half-hourly, hourly, and similar data.
"""
function build_time_axis(nrows::Int; start_date::Date, step_minutes::Int)
    t0 = DateTime(start_date)
    dt = Minute(step_minutes)
    return [t0 + (i - 1) * dt for i in 1:nrows]
end

"""
load_profile_column(parquet_file, customer_col)

Reads a single column from parquet using DuckDB, returns Vector{Float64}.
"""
function load_profile_column(parquet_file::String, customer_col::String)
    con = DBInterface.connect(DuckDB.DB, ":memory:")
    prof_df = DBInterface.execute(con, "SELECT \"$(customer_col)\" AS y FROM read_parquet('$(parquet_file)')") |> DataFrame
    DBInterface.close!(con)
    return Float64.(prof_df.y)
end

"""
build_alpha(y; alpha_min, alpha_max)

Returns alpha(t) = y / mean(y), clipped.
This keeps the alpha definition consistent across all scripts.
"""
function build_alpha(y::Vector{Float64}; alpha_min::Float64, alpha_max::Float64)
    μ = mean(y)
    μ > 0 || error("Mean of profile column is non-positive, cannot build alpha")
    a = y ./ μ
    return clamp.(a, alpha_min, alpha_max)
end
