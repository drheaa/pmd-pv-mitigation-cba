# This script is only about understanding the customer load profile data.
# Nothing here touches the network model or PowerModelsDistribution yet.
#
# The goal is to answer these basic questions:
# - Does the load profile change over the day?
# - Does it change over the year?
# - Are the values reasonable and not obviously broken?
#
# If this script makes sense and the plots look sane,
# then the data is safe to use for time-series power flow later.

# --------------------------------------------------
# 0) Environment setup
# --------------------------------------------------

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

# Packages needed for this script
using DuckDB          # to read parquet files easily
using DBInterface     # interface for database queries
using DataFrames      # tables
using Statistics      # mean, std, etc
using Dates           # time handling
using Plots           # plotting
using CSV             # saving tables

# --------------------------------------------------
# 1) User settings
# --------------------------------------------------

# Root of the project, same as snapshot scripts
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"

# Folder that contains the D-Suite customer load profiles
PROFILES_ROOT = joinpath(ROOT, "data/raw/dsuite_load_profiles/profiles")

# Choose one year of data to inspect
YEAR = 2023

# Choose which customer column to inspect
# If set to "AUTO", the script automatically picks a reasonable one
CUSTOMER_COL = "AUTO"

# Starting date for building a time axis
# This does not need to be exact for sanity checking
START_DATE = Date(YEAR, 1, 1)

# Expected time resolution of the data in minutes
# D-Suite profiles are usually 30-minute resolution
STEP_MINUTES = 30

# How many days to plot for the raw time-series example
SAMPLE_DAYS = 14
SAMPLE_START_DAY = 1   # 1 means start from January 1

# Output folders
OUTDIR = joinpath(ROOT, "results", "time_series", "profile_sanity", "year=$(YEAR)")
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR)
mkpath(TBLDIR)

# --------------------------------------------------
# 2) Helper functions
# --------------------------------------------------

# Find the parquet file corresponding to the chosen year
function find_parquet_for_year(profiles_root::String, year::Int)
    year_dir = joinpath(profiles_root, "year=$(year)")
    isdir(year_dir) || error("Year folder not found")

    files = filter(f -> endswith(lowercase(f), ".parquet"), readdir(year_dir))
    isempty(files) && error("No parquet files found")

    # Usually there is only one parquet file per year
    return joinpath(year_dir, files[1])
end

# Identify numeric columns in the parquet file
# Customer loads are stored as numeric columns
function numeric_columns(schema_df::DataFrame)
    mask = occursin.(r"(INTEGER|BIGINT|DOUBLE|FLOAT|REAL|DECIMAL)", schema_df.column_type)
    return Vector{String}(schema_df.column_name[mask])
end

# Automatically pick one customer column that looks like a real load
# This explicitly skips metadata columns such as "hh" and "year"
function pick_customer_column_auto(con, parquet_file::String, cols::Vector{String})
    best_col = nothing
    best_score = -Inf

    for c in cols

        # Skip known non-load metadata columns
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
                # skip values that cannot be converted cleanly
            end
        end

        isempty(vals) && continue

        # Heuristic:
        # - variance ensures the series is not constant
        # - nonzero fraction avoids mostly-zero columns
        nonzero_frac = count(x -> abs(x) > 1e-9, vals) / length(vals)
        score = var(vals) * nonzero_frac

        if score > best_score
            best_score = score
            best_col = c
        end
    end

    best_col === nothing && error("No suitable customer load column found")
    return best_col
end


# Build a simple time axis assuming fixed time steps
function build_time_axis(nrows::Int; start_date::Date, step_minutes::Int)
    t0 = DateTime(start_date)
    dt = Minute(step_minutes)
    return [t0 + (i-1)*dt for i in 1:nrows]
end

# Simple quantile helper
function quantile_approx(x::Vector{Float64}, q::Float64)
    isempty(x) && return NaN
    xs = sort(x)
    idx = clamp(Int(ceil(q * length(xs))), 1, length(xs))
    return xs[idx]
end

# --------------------------------------------------
# 3) Load parquet file and inspect schema
# --------------------------------------------------

parquet_file = find_parquet_for_year(PROFILES_ROOT, YEAR)
println("Using parquet file: ", parquet_file)

# Create an in-memory DuckDB database
con = DBInterface.connect(DuckDB.DB, ":memory:")

# Count number of rows
nrows_df = DBInterface.execute(con,
    "SELECT COUNT(*) AS n FROM read_parquet('$(parquet_file)')"
) |> DataFrame
nrows = Int(nrows_df.n[1])
println("Number of rows: ", nrows)

# Inspect schema to find numeric columns
schema_df = DBInterface.execute(con,
    "DESCRIBE SELECT * FROM read_parquet('$(parquet_file)')"
) |> DataFrame
numcols = numeric_columns(schema_df)
println("Numeric columns found: ", length(numcols))

isempty(numcols) && error("No numeric columns found")

# Choose customer column
customer_col =
    CUSTOMER_COL != "AUTO" ? CUSTOMER_COL :
    pick_customer_column_auto(con, parquet_file, numcols)

println("Chosen customer column: ", customer_col)

# --------------------------------------------------
# 4) Read the customer load time-series
# --------------------------------------------------

df = DBInterface.execute(con,
    "SELECT \"$(customer_col)\" AS y FROM read_parquet('$(parquet_file)')"
) |> DataFrame

# Convert to Float64 vector and track missing values
y = Vector{Union{Missing,Float64}}(undef, nrows)
missing_count = 0

for i in 1:nrows
    v = df.y[i]
    if v === missing
        y[i] = missing
        missing_count += 1
    else
        y[i] = Float64(v)
    end
end

# Remove missing values for statistics
y_nomiss = [Float64(v) for v in y if v !== missing]

# --------------------------------------------------
# 5) Basic sanity statistics
# --------------------------------------------------

stats = DataFrame([
    (
        year = YEAR,
        customer = customer_col,
        rows = nrows,
        missing = missing_count,
        min = minimum(y_nomiss),
        p05 = quantile_approx(y_nomiss, 0.05),
        median = quantile_approx(y_nomiss, 0.50),
        p95 = quantile_approx(y_nomiss, 0.95),
        max = maximum(y_nomiss),
        mean = mean(y_nomiss),
        std = std(y_nomiss)
    )
])

CSV.write(joinpath(TBLDIR, "profile_stats.csv"), stats)
println(stats)

# --------------------------------------------------
# 6) Build time axis
# --------------------------------------------------

t = build_time_axis(nrows; start_date=START_DATE, step_minutes=STEP_MINUTES)

# --------------------------------------------------
# 7) Plot raw load for a short window
# --------------------------------------------------

steps_per_day = Int(round(24*60 / STEP_MINUTES))
i0 = (SAMPLE_START_DAY - 1) * steps_per_day + 1
i1 = min(i0 + SAMPLE_DAYS*steps_per_day - 1, nrows)

p_raw = plot(
    t[i0:i1],
    [v === missing ? NaN : Float64(v) for v in y[i0:i1]],
    xlabel="Time",
    ylabel="Load",
    title="Raw customer load sample"
)
savefig(p_raw, joinpath(FIGDIR, "raw_timeseries_sample.png"))

# --------------------------------------------------
# 8) Diurnal profile
# --------------------------------------------------

slot_sum = zeros(Float64, steps_per_day)
slot_cnt = zeros(Int, steps_per_day)

for i in 1:nrows
    v = y[i]
    v === missing && continue
    slot = ((i-1) % steps_per_day) + 1
    slot_sum[slot] += Float64(v)
    slot_cnt[slot] += 1
end

slot_mean = [slot_cnt[s] > 0 ? slot_sum[s]/slot_cnt[s] : NaN for s in 1:steps_per_day]
hours = collect(0:(24/steps_per_day):(24 - 24/steps_per_day))

p_diurnal = plot(
    hours,
    slot_mean,
    xlabel="Hour of day",
    ylabel="Average load",
    title="Diurnal load profile",
    marker=:circle,
    markersize=2
)
savefig(p_diurnal, joinpath(FIGDIR, "diurnal_profile.png"))

# --------------------------------------------------
# 9) Seasonal trend
# --------------------------------------------------

ndays = Int(ceil(nrows / steps_per_day))
daily_mean = fill(NaN, ndays)
dates = [START_DATE + Day(d-1) for d in 1:ndays]

for d in 1:ndays
    a = (d-1)*steps_per_day + 1
    b = min(d*steps_per_day, nrows)
    vals = [Float64(y[i]) for i in a:b if y[i] !== missing]
    daily_mean[d] = isempty(vals) ? NaN : mean(vals)
end

p_season = plot(
    dates,
    daily_mean,
    xlabel="Date",
    ylabel="Daily mean load",
    title="Seasonal load variation"
)
savefig(p_season, joinpath(FIGDIR, "seasonal_profile.png"))

# --------------------------------------------------
# 10) Distribution check
# --------------------------------------------------

p_hist = histogram(
    y_nomiss,
    bins=50,
    xlabel="Load",
    ylabel="Count",
    title="Load distribution"
)
savefig(p_hist, joinpath(FIGDIR, "load_histogram.png"))

DBInterface.close!(con)

println("Finished profile sanity check")
println("Results saved to:")
println(OUTDIR)
