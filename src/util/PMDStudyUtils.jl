module PMDStudyUtils

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
using Serialization

const PMD = PowerModelsDistribution

include("profiles.jl")
include("load_scaling.jl")
include("sequences.jl")
include("pf_solve.jl")
include("metrics.jl")
include("paper_outputs.jl")
include("plots.jl")
include("debug.jl")

export
    # profiles
    find_parquet_for_year,
    numeric_columns,
    pick_customer_column_auto,
    build_time_axis,
    load_profile_column,
    build_alpha,

    # load scaling
    cache_load_bases,
    apply_load_scaling!,
    print_load_summary,
    sanity_check_no_mitigation_loads,

    # sequences
    seq_components,
    seq_magnitudes,

    # pf solve
    solve_pf_3w_projected,
    solve_pf_4w_explicit,

    # metrics
    pf_voltage_metrics,
    pf_system_vuf_max,
    pf_system_vseq_max,
    pf_transformer_or_slack_p_kw,
    feeder_head_currents,
    bus_voltage_table,

    # paper outputs
    duration_curve,
    write_duration_curves_csv,
    select_worst_timesteps,

    # plots
    plot_voltage_envelope,
    plot_series,
    plot_scatter,
    plot_duration_curve,
    plot_worst_voltage_profile,

    # debug
    save_debug_snapshot

end
