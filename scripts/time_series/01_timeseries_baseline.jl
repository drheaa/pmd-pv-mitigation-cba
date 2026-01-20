# Build a first time-series baseline for one feeder:
# - Read one real customer profile column from parquet
# - Convert it into a time scaling factor alpha(t)
# - Scale all feeder loads by alpha(t)
# - Run PF at each timestep (or every Nth timestep)
# - Save voltage metrics over time
#
# This is a stepping stone:
# Later, replace "global alpha(t)" with full customer to load allocation.

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

using PowerModelsDistribution
using JuMP
using Ipopt
using DuckDB
using DBInterface
using DataFrames
using Dates
using Statistics
using CSV
using Plots

const PMD = PowerModelsDistribution
PMD.silence!()

# --------------------------------------------------
# 0) Settings that control the run
# --------------------------------------------------

# Project root, same style as snapshot scripts
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"

# Choose one feeder that already works in snapshot scripts
NET = "spd_s"
master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_scaled.dss")

# Load profile folder and one year
PROFILES_ROOT = joinpath(ROOT, "data/raw/dsuite_load_profiles/profiles")
YEAR = 2023

# Pick one customer column from the parquet file
# "AUTO" means choose automatically, but skip metadata columns like hh and year
CUSTOMER_COL = "AUTO"

# Time axis assumptions
# D-Suite profiles are typically 30 minute resolution
START_DATE = Date(YEAR, 1, 1)
STEP_MINUTES = 30

# Voltage thresholds used only for reporting
VMIN_PU = 0.94
VMAX_PU = 1.10

# Base scaling already used in snapshot work
# This is separated from time-series scaling
LOAD_ALPHA_BASE = 1.0

# STRIDE controls how many timesteps are skipped
# STRIDE=1 runs all 17520 steps
# STRIDE=2 runs every hour
# STRIDE=4 runs every 2 hours
STRIDE = 4

# Alpha clipping keeps the first run stable
# This avoids extreme stress while testing the pipeline
ALPHA_MIN = 0.3
ALPHA_MAX = 2.5

# Output folders
OUTDIR = joinpath(ROOT, "results", "time_series", "baseline_pf", NET, "year=$(YEAR)_stride=$(STRIDE)")
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR)
mkpath(TBLDIR)

# --------------------------------------------------
# 1) Helper functions for reading parquet files
# --------------------------------------------------

# Find the parquet file for the chosen year
function find_parquet_for_year(profiles_root::String, year::Int)
    year_dir = joinpath(profiles_root, "year=$(year)")
    isdir(year_dir) || error("Year folder not found")
    files = filter(f -> endswith(lowercase(f), ".parquet"), readdir(year_dir))
    isempty(files) && error("No parquet files found")
    return joinpath(year_dir, files[1])
end

# DuckDB schema output includes column_name and column_type
# This function keeps only numeric columns
function numeric_columns(schema_df::DataFrame)
    mask = occursin.(r"(INTEGER|BIGINT|DOUBLE|FLOAT|REAL|DECIMAL)", schema_df.column_type)
    return Vector{String}(schema_df.column_name[mask])
end

# Automatically pick a customer load column
# This skips metadata columns and picks something with variance and nonzero values
function pick_customer_column_auto(con, parquet_file::String, cols::Vector{String})
    best_col = nothing
    best_score = -Inf

    for c in cols
        # Skip known metadata columns
        if c == "hh" || c == "year"
            continue
        end

        # Read only a small sample for speed
        q = "SELECT \"$(c)\" AS y FROM read_parquet('$(parquet_file)') LIMIT 2000"
        df = DBInterface.execute(con, q) |> DataFrame

        vals = Float64[]
        for v in df.y
            v === missing && continue
            try
                push!(vals, Float64(v))
            catch
                # Skip values that cannot convert cleanly
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

# Build a DateTime axis using a start date and constant step size
function build_time_axis(nrows::Int; start_date::Date, step_minutes::Int)
    t0 = DateTime(start_date)
    dt = Minute(step_minutes)
    return [t0 + (i-1)*dt for i in 1:nrows]
end

# --------------------------------------------------
# 2) Helper functions for scaling loads and running PF
# --------------------------------------------------

# Multiply all loads in the engineering model by alpha
# This matches the snapshot scripts style
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

# Run PF and extract simple system level metrics
# This returns vmin, vmax, and violation counts
function run_pf_metrics(eng::Dict{String,Any}; vmin_pu=0.94, vmax_pu=1.10)
    # Convert engineering model to math model for the solver
    math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)

    # Use IPOPT with silent output
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

    # Solve power flow using an IVR unbalanced model
    pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)

    # Extract bus voltage magnitudes from the solution
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

    if isempty(vmins)
        return (status=string(pf["termination_status"]), vmin=NaN, vmax=NaN, n_under=0, n_over=0)
    end

    vmin_sys = minimum(vmins)
    vmax_sys = maximum(vmaxs)

    n_under = count(x -> x < vmin_pu, vmins)
    n_over  = count(x -> x > vmax_pu, vmaxs)

    return (status=string(pf["termination_status"]), vmin=vmin_sys, vmax=vmax_sys, n_under=n_under, n_over=n_over)
end

# --------------------------------------------------
# 3) Read one customer profile and turn it into alpha(t)
# --------------------------------------------------

parquet_file = find_parquet_for_year(PROFILES_ROOT, YEAR)
println("Parquet file: ", parquet_file)

con = DBInterface.connect(DuckDB.DB, ":memory:")

# Count rows, expected 17520 for half-hour data across a year
nrows_df = DBInterface.execute(con, "SELECT COUNT(*) AS n FROM read_parquet('$(parquet_file)')") |> DataFrame
nrows = Int(nrows_df.n[1])
println("Profile rows: ", nrows)

# Read schema and pick numeric columns
schema_df = DBInterface.execute(con, "DESCRIBE SELECT * FROM read_parquet('$(parquet_file)')") |> DataFrame
numcols = numeric_columns(schema_df)

# Pick customer column
customer_col =
    CUSTOMER_COL != "AUTO" ? CUSTOMER_COL :
    pick_customer_column_auto(con, parquet_file, numcols)

println("Chosen customer column: ", customer_col)

# Read the chosen column into memory
prof_df = DBInterface.execute(con, "SELECT \"$(customer_col)\" AS y FROM read_parquet('$(parquet_file)')") |> DataFrame
DBInterface.close!(con)

y = Float64.(prof_df.y)

# Convert y(t) into alpha(t) by dividing by the mean
# This makes average alpha close to 1
y_mean = mean(y)
alpha = y ./ y_mean

# Clip alpha to keep the first run stable
alpha = clamp.(alpha, ALPHA_MIN, ALPHA_MAX)

# Build time axis for the saved output table and plots
t = build_time_axis(nrows; start_date=START_DATE, step_minutes=STEP_MINUTES)

# --------------------------------------------------
# 4) Parse the feeder once and prepare a baseline engineering model
# --------------------------------------------------

println("Parsing feeder: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

# Apply the base scaling if needed
scale_loads!(eng0, LOAD_ALPHA_BASE)

# --------------------------------------------------
# 5) Time loop: apply alpha(t), run PF, store metrics
# --------------------------------------------------

rows = NamedTuple[]

for k in 1:STRIDE:nrows
    # Make a fresh copy of the feeder data for this timestep
    eng = deepcopy(eng0)

    # Apply time-varying scaling for this timestep
    scale_loads!(eng, alpha[k])

    # Run PF and collect metrics
    r = run_pf_metrics(eng; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # Store one row for this timestep
    push!(rows, (
        net = NET,
        year = YEAR,
        timestep = k,
        time = t[k],
        customer_col = customer_col,
        alpha_t = alpha[k],
        pf_status = r.status,
        vmin_pu = r.vmin,
        vmax_pu = r.vmax,
        n_under = r.n_under,
        n_over = r.n_over
    ))

    # Print progress sometimes, not every step
    if k % (STRIDE * 200) == 1
        println("Progress ", k, " / ", nrows, " | vmin=", r.vmin, " vmax=", r.vmax)
    end
end

df = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "timeseries_baseline_pf_metrics.csv"), df)

# --------------------------------------------------
# 6) Plots that show how voltages change across the year
# --------------------------------------------------

p1 = plot(df.time, df.vmin_pu;
    xlabel="Time",
    ylabel="System vmin (pu)",
    title="Time-series baseline PF | vmin | $(NET) | year=$(YEAR) | stride=$(STRIDE)",
    legend=false
)
hline!([VMIN_PU], linestyle=:dash, color=:red)
savefig(p1, joinpath(FIGDIR, "ts_vmin.png"))

p2 = plot(df.time, df.vmax_pu;
    xlabel="Time",
    ylabel="System vmax (pu)",
    title="Time-series baseline PF | vmax | $(NET) | year=$(YEAR) | stride=$(STRIDE)",
    legend=false
)
hline!([VMAX_PU], linestyle=:dash, color=:red)
savefig(p2, joinpath(FIGDIR, "ts_vmax.png"))

p3 = plot(df.time, df.n_under;
    xlabel="Time",
    ylabel="Under-voltage bus count",
    title="Time-series baseline PF | under count | $(NET)",
    legend=false
)
savefig(p3, joinpath(FIGDIR, "ts_n_under.png"))

p4 = plot(df.time, df.n_over;
    xlabel="Time",
    ylabel="Over-voltage bus count",
    title="Time-series baseline PF | over count | $(NET)",
    legend=false
)
savefig(p4, joinpath(FIGDIR, "ts_n_over.png"))

println("Saved results to: ", OUTDIR)

println("\nSummary checks")
println("alpha_t min = ", minimum(df.alpha_t), " | max = ", maximum(df.alpha_t))
println("alpha hits ALPHA_MIN count = ", count(==(ALPHA_MIN), df.alpha_t))
println("alpha hits ALPHA_MAX count = ", count(==(ALPHA_MAX), df.alpha_t))

println("vmin over run: min = ", minimum(df.vmin_pu), " | median = ", median(df.vmin_pu), " | max = ", maximum(df.vmin_pu))
println("vmax over run: min = ", minimum(df.vmax_pu), " | median = ", median(df.vmax_pu), " | max = ", maximum(df.vmax_pu))

println("timesteps with any under-voltage buses = ", count(>(0), df.n_under))
println("timesteps with any over-voltage buses  = ", count(>(0), df.n_over))


println("Done.")
