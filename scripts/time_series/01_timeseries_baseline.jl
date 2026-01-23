# ================================================================
# 01_timeseries_baseline.jl
# ================================================================
# Time-series baseline PF for one feeder using a global alpha(t).
# Utilities are loaded from src/util/PMDStudyUtils.jl.
# ================================================================

import Pkg
# Pkg.activate(joinpath(@__DIR__, "..", ".."))
# Pkg.instantiate()

using PowerModelsDistribution
using DataFrames
using Dates
using Statistics
using CSV
using Plots
using DuckDB
using DBInterface

const PMD = PowerModelsDistribution
PMD.silence!()

# Local module include for scripts
include(joinpath(@__DIR__, "..", "..", "src", "util", "PMDStudyUtils.jl"))
using .PMDStudyUtils

# --------------------------------------------------
# 0) Run settings
# --------------------------------------------------

ROOT = joinpath(@__DIR__, "../..")

NET = "spd_s_4w"
master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_baseline_4w.dss")

PROFILES_ROOT = joinpath(ROOT, "data/raw/dsuite_load_profiles/profiles")
YEAR = 2023
CUSTOMER_COL = "AUTO"

START_DATE   = Date(YEAR, 1, 1)
STEP_MINUTES = 30

VMIN_PU = 0.90
VMAX_PU = 1.10

LOAD_ALPHA_BASE = 2.5
STRIDE = 4

ALPHA_MIN = 0.3
ALPHA_MAX = 2.5

MAX_STEPS = 17520
N_WORST = 5

OUTDIR = joinpath(ROOT, "results", "time_series", "baseline_pf", NET, "year=$(YEAR)_stride=$(STRIDE)_K=$(MAX_STEPS)")
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR)
mkpath(TBLDIR)

# --------------------------------------------------
# 1) Load profile, build alpha(t), build time axis
# --------------------------------------------------

parquet_file = PMDStudyUtils.find_parquet_for_year(PROFILES_ROOT, YEAR)
println("Parquet file: ", parquet_file)

con = DBInterface.connect(DuckDB.DB, ":memory:")
schema_df = DBInterface.execute(con, "DESCRIBE SELECT * FROM read_parquet('$(parquet_file)')") |> DataFrame
numcols = PMDStudyUtils.numeric_columns(schema_df)

customer_col =
    CUSTOMER_COL != "AUTO" ? CUSTOMER_COL :
    PMDStudyUtils.pick_customer_column_auto(con, parquet_file, numcols)
DBInterface.close!(con)

println("Chosen customer column: ", customer_col)

y = load_profile_column(parquet_file, customer_col)
alpha = build_alpha(y; alpha_min=ALPHA_MIN, alpha_max=ALPHA_MAX)
time_axis = build_time_axis(length(y); start_date=START_DATE, step_minutes=STEP_MINUTES)

# --------------------------------------------------
# 2) Parse feeder once, apply baseline scaling once, cache base
# --------------------------------------------------

println("\nParsing feeder: ", master_dss)
eng = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

cache0 = cache_load_bases(eng)
apply_load_scaling!(eng, cache0, LOAD_ALPHA_BASE)

cache_base = cache_load_bases(eng)

print_load_summary(eng; title="Baseline load summary after LOAD_ALPHA_BASE")
sanity_check_no_mitigation_loads(eng)

# --------------------------------------------------
# 3) Choose PF representation
# --------------------------------------------------
# For explicit 4-wire studies, switch the solver call to solve_pf_4w_explicit.
# If explicit 4-wire errors, keep projected 3-wire as a temporary reference.
#
# Solver function handle is set here so later scripts stay consistent.
solve_pf = solve_pf_3w_projected
# solve_pf = solve_pf_4w_explicit

# --------------------------------------------------
# 4) Time loop
# --------------------------------------------------

rows = NamedTuple[]
last_k = min(length(alpha), MAX_STEPS)

pf_debug = nothing
math_debug = nothing
k_debug = nothing

for k in 1:STRIDE:last_k
    apply_load_scaling!(eng, cache_base, alpha[k])

    pf, math = solve_pf(eng)

    if pf_debug === nothing
        pf_debug = pf
        math_debug = math
        k_debug = k
        save_debug_snapshot(OUTDIR, k, pf, math)
    end

    mv = pf_voltage_metrics(pf; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)
    vuf = pf_system_vuf_max(pf)
    vseq = pf_system_vseq_max(pf)
    p_kw = pf_transformer_or_slack_p_kw(pf)
    Ih = feeder_head_currents(pf)

    push!(rows, (
        net = NET,
        year = YEAR,
        timestep = k,
        time = time_axis[k],
        customer_col = customer_col,
        alpha_t = alpha[k],
        pf_status = mv.status,
        vmin_pu = mv.vmin,
        vmax_pu = mv.vmax,
        vuf_max = vuf,
        v0_max = vseq.V0_max,
        v2_max = vseq.V2_max,
        v2_over_v1_max = vseq.V2_over_V1_max,
        tx_p_kw = p_kw,
        i_a = Ih.Ia,
        i_b = Ih.Ib,
        i_c = Ih.Ic,
        i_max_phase = Ih.Imax,
        i_mean_phase = Ih.Imean,
        i0_mag = Ih.I0,
        i2_mag = Ih.I2,
        i2_over_i1 = Ih.I2_over_I1,
        n_under = mv.n_under,
        n_over = mv.n_over
    ))

    if k % (STRIDE * 200) == 1
        println("Progress ", k, " / ", length(alpha),
            " | vmin=", mv.vmin,
            " vmax=", mv.vmax,
            " vuf=", vuf,
            " p_kw=", p_kw
        )
    end
end

df = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "timeseries_baseline_pf_metrics.csv"), df)
println("\nSaved table: ", joinpath(TBLDIR, "timeseries_baseline_pf_metrics.csv"))

# --------------------------------------------------
# 5) Summary checks
# --------------------------------------------------

println("\nSummary checks")
println("alpha_t min = ", minimum(df.alpha_t), " | max = ", maximum(df.alpha_t))
println("alpha hits ALPHA_MIN count = ", count(==(ALPHA_MIN), df.alpha_t))
println("alpha hits ALPHA_MAX count = ", count(==(ALPHA_MAX), df.alpha_t))

println("vmin over run: min = ", minimum(df.vmin_pu),
    " | median = ", median(df.vmin_pu),
    " | max = ", maximum(df.vmin_pu)
)
println("vmax over run: min = ", minimum(df.vmax_pu),
    " | median = ", median(df.vmax_pu),
    " | max = ", maximum(df.vmax_pu)
)
println("vuf_max over run: min = ", minimum(df.vuf_max),
    " | median = ", median(df.vuf_max),
    " | max = ", maximum(df.vuf_max)
)

println("timesteps with any under-voltage buses = ", count(>(0), df.n_under))
println("timesteps with any over-voltage buses  = ", count(>(0), df.n_over))

if nrow(df) >= 5
    println("corr(alpha_t, vmin_pu) = ", cor(df.alpha_t, df.vmin_pu))
    if all(isfinite.(df.tx_p_kw))
        println("corr(alpha_t, tx_p_kw) = ", cor(df.alpha_t, df.tx_p_kw))
    end
end

# --------------------------------------------------
# 6) Research-ready plots
# --------------------------------------------------

plot_voltage_envelope(
    df;
    vmin_limit=VMIN_PU,
    vmax_limit=VMAX_PU,
    title_str="Voltage envelope | $(NET) | year=$(YEAR) | stride=$(STRIDE)",
    outpath=joinpath(FIGDIR, "ts_voltage_envelope.png")
)

plot_series(df.time, df.vuf_max;
    xlabel="Time",
    ylabel="VUF (max across buses)",
    title_str="Voltage unbalance (VUF) | $(NET) | year=$(YEAR) | stride=$(STRIDE)",
    outpath=joinpath(FIGDIR, "ts_vuf_over_time.png")
)

plot_series(df.time, df.tx_p_kw;
    xlabel="Time",
    ylabel="Active power proxy (kW)",
    title_str="Substation active power proxy | $(NET) | year=$(YEAR) | stride=$(STRIDE)",
    outpath=joinpath(FIGDIR, "ts_tx_p_kw_over_time.png")
)

plot_scatter(df.alpha_t, df.vmin_pu;
    xlabel="alpha(t)",
    ylabel="system vmin (pu)",
    title_str="Sanity: alpha(t) vs vmin",
    outpath=joinpath(FIGDIR, "sanity_alpha_vs_vmin.png")
)

plot_scatter(df.alpha_t, df.vuf_max;
    xlabel="alpha(t)",
    ylabel="system VUF max",
    title_str="Sanity: alpha(t) vs VUF max",
    outpath=joinpath(FIGDIR, "sanity_alpha_vs_vuf.png")
)

# Duration curves
dc_vmin = duration_curve(df.vmin_pu; high_is_worse=false)
dc_vuf  = duration_curve(df.vuf_max; high_is_worse=true)

plot_duration_curve(dc_vmin.pct, dc_vmin.x;
    xlabel="Percentile (%)",
    ylabel="vmin (pu)",
    title_str="Duration curve: system vmin",
    outpath=joinpath(FIGDIR, "duration_vmin.png")
)

plot_duration_curve(dc_vuf.pct, dc_vuf.x;
    xlabel="Percentile (%)",
    ylabel="VUF max",
    title_str="Duration curve: system VUF max",
    outpath=joinpath(FIGDIR, "duration_vuf.png")
)

write_duration_curves_csv(joinpath(TBLDIR, "duration_curves.csv"), df)

# Worst timestep diagnostics
worst = select_worst_timesteps(df, N_WORST)
CSV.write(joinpath(TBLDIR, "worst_timesteps.csv"), worst)

for r in eachrow(worst)
    k = Int(r.timestep)
    PMDStudyUtils.apply_load_scaling!(eng, cache_base, alpha[k])
    pf_k, _ = solve_pf(eng)

    vt = PMDStudyUtils.bus_voltage_table(pf_k)
    CSV.write(joinpath(TBLDIR, "worst_timestep_bus_voltages_k=$(k).csv"), vt)

    PMDStudyUtils.plot_worst_voltage_profile(vt;
        title_str="Worst timestep voltage profile | k=$(k) | vmin=$(round(r.vmin_pu, digits=4)) | alpha=$(round(r.alpha_t, digits=3))",
        outpath=joinpath(FIGDIR, "worst_voltage_profile_k=$(k).png")
    )
end

println("\nSaved results to: ", OUTDIR)
println("Debug snapshot timestep k = ", k_debug)
println("Done.")
