# Build a first time-series baseline for one feeder:
# 1) Read one real customer profile column from parquet
# 2) Convert it into a time scaling factor alpha(t)
# 3) Scale all feeder loads by alpha(t)
# 4) Run PF at each timestep (or every Nth timestep)
# 5) Save voltage metrics over time
#
# EXTRA ACTIONS IN THIS VERSION
# 1) Print baseline load magnitude from the parsed engineering model
# 2) Create voltage profile plots for the worst timesteps (lowest vmin_pu)
# 3) Create a simple topology plot if bus coordinate CSV exists
#
# NOTE
# This is still a "global alpha(t)" baseline.
# Later work can map different customer profiles to different loads and phases.

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

# Project root
ROOT = "/mnt/c/Users/auc009/OneDrive - CSIRO/Documents/power-models-distribution/pmd_pv_experiments"

# Choose one feeder
NET = "spd_u"
master_dss = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "master_scaled.dss")

# Optional bus coordinates file for topology plotting
# Many D-Suite feeders include a file like opendss_xy_spd_r_scaled.csv
buscoords_csv = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET, "opendss_xy_$(NET)_scaled.csv")

# Load profile folder and one year
PROFILES_ROOT = joinpath(ROOT, "data/raw/dsuite_load_profiles/profiles")
YEAR = 2023

# Pick one customer column from the parquet file
# "AUTO" means choose automatically, but skip metadata columns like hh and year
CUSTOMER_COL = "AUTO"

# Time axis assumptions
START_DATE = Date(YEAR, 1, 1)
STEP_MINUTES = 30

# Voltage thresholds used only for reporting (per unit)
VMIN_PU = 0.94
VMAX_PU = 1.10

# Nominal LN voltage for converting pu to volts in plots
VBASE_LN = 230.0

# Base scaling already used in snapshot work
# This is separated from time-series scaling
LOAD_ALPHA_BASE = 2.5

# STRIDE controls how many timesteps are skipped
# STRIDE=1 runs all 17520 steps
# STRIDE=2 runs every hour
# STRIDE=4 runs every 2 hours
STRIDE = 4

# Alpha clipping keeps the first run stable
ALPHA_MIN = 0.3
ALPHA_MAX = 2.5

# Number of worst timesteps to plot
N_WORST_PLOTS = 5

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
    isdir(year_dir) || error("Year folder not found: $year_dir")
    files = filter(f -> endswith(lowercase(f), ".parquet"), readdir(year_dir))
    isempty(files) && error("No parquet files found in: $year_dir")
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

# Solve PF and return both result and math model
# Returning math allows voltage plots using the same solved state
function solve_pf(eng::Dict{String,Any})
    math = PMD.transform_data_model(eng; multinetwork=false, kron_reduce=true, phase_project=true)
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")
    pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)
    return pf, math
end

# Extract system level metrics from PF result
function pf_metrics(pf::Dict{String,Any}; vmin_pu=0.94, vmax_pu=1.10)
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
# 3) Diagnostic: baseline load magnitude in PMD
# --------------------------------------------------

function print_baseline_load_summary(eng::Dict{String,Any})
    total_pd = 0.0
    total_qd = 0.0
    nloads = 0

    if haskey(eng, "load")
        for (_, ld_any) in eng["load"]
            ld = ld_any::Dict{String,Any}

            if haskey(ld, "pd")
                total_pd += sum(ld["pd"])
            elseif haskey(ld, "pd_nom")
                total_pd += sum(ld["pd_nom"])
            end

            if haskey(ld, "qd")
                total_qd += sum(ld["qd"])
            elseif haskey(ld, "qd_nom")
                total_qd += sum(ld["qd_nom"])
            end

            nloads += 1
        end
    end

    println("\nBaseline load summary from engineering model")
    println("  number of loads = ", nloads)
    println("  total P (kW)    = ", total_pd)
    println("  total Q (kvar)  = ", total_qd)
end

# --------------------------------------------------
# 4) Plot helpers for topology and voltage profiles
# --------------------------------------------------

# Choose a reference source bus name
# Some feeders include sourcebus and sourcebusz as separate buses
function pick_source_bus(eng::Dict{String,Any})
    if haskey(eng, "bus")
        names = Set(keys(eng["bus"]))
        if "sourcebusz" in names
            return "sourcebusz"
        elseif "sourcebus" in names
            return "sourcebus"
        end
    end
    return "sourcebus"
end

# Build a DataFrame of line connections from the engineering model
function make_lines_df_from_eng(eng::Dict{String,Any})
    rows = NamedTuple[]

    if haskey(eng, "line")
        for (_, ln_any) in eng["line"]
            ln = ln_any::Dict{String,Any}
            f_bus = ln["f_bus"]
            t_bus = ln["t_bus"]
            f_phases = ln["f_connections"]
            t_phases = ln["t_connections"]
            @assert length(f_phases) == length(t_phases)

            length_km = (haskey(ln, "length") ? ln["length"] : 0.0) / 1000.0

            push!(rows, (Bus1=f_bus, Bus2=t_bus, phases=f_phases, length_km=length_km))
        end
    elseif haskey(eng, "branch")
        for (_, br_any) in eng["branch"]
            br = br_any::Dict{String,Any}
            f_bus = br["f_bus"]
            t_bus = br["t_bus"]
            f_phases = br["f_connections"]
            t_phases = br["t_connections"]
            @assert length(f_phases) == length(t_phases)

            length_km = (haskey(br, "length") ? br["length"] : 0.0) / 1000.0

            push!(rows, (Bus1=f_bus, Bus2=t_bus, phases=f_phases, length_km=length_km))
        end
    else
        error("No line or branch data found in engineering model")
    end

    return DataFrame(rows)
end

# Compute distances (km) from a reference bus using a simple graph walk
function compute_bus_distances(lines_df::DataFrame; source_bus::String)
    adj = Dict{String, Vector{Tuple{String, Float64}}}()

    for r in eachrow(lines_df)
        b1 = r.Bus1
        b2 = r.Bus2
        len = r.length_km
        push!(get!(adj, b1, Tuple{String,Float64}[]), (b2, len))
        push!(get!(adj, b2, Tuple{String,Float64}[]), (b1, len))
    end

    dist = Dict{String,Float64}(source_bus => 0.0)
    queue = [source_bus]

    while !isempty(queue)
        u = popfirst!(queue)
        for (v, w) in get(adj, u, Tuple{String,Float64}[])
            if !haskey(dist, v)
                dist[v] = dist[u] + w
                push!(queue, v)
            end
        end
    end

    return dist
end

# Extract solved bus voltages in volts.
# Returns two dictionaries:
# 1) buses_by_name: keyed by math bus "name"
# 2) buses_by_id: keyed by math bus id (the key in math["bus"])
#
# Reason: distance keys come from engineering line data, and sometimes match bus ids,
# sometimes match bus names. This keeps both so the rest of the script can choose.
function solved_bus_vm_volts_dual(pf::Dict{String,Any}, math::Dict{String,Any}; vbase_ln::Float64)
    sol_bus = pf["solution"]["bus"]

    buses_by_name = Dict{String, Dict{String,Any}}()
    buses_by_id   = Dict{String, Dict{String,Any}}()

    for (bus_id, bus_data_any) in math["bus"]
        bus_data = bus_data_any::Dict{String,Any}
        name = string(bus_data["name"])

        if !haskey(sol_bus, bus_id)
            continue
        end

        sb = sol_bus[bus_id]

        vm_pu =
            haskey(sb, "vm") ? sb["vm"] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"].^2 .+ sb["vi"].^2) :
            nothing

        vm_pu === nothing && continue

        vmV = vm_pu .* vbase_ln

        entry = Dict(
            "vma" => [vmV[1]],
            "vmb" => length(vmV) >= 2 ? [vmV[2]] : [vmV[1]],
            "vmc" => length(vmV) >= 3 ? [vmV[3]] : [vmV[1]]
        )

        buses_by_name[name] = deepcopy(entry)
        buses_by_id[string(bus_id)] = deepcopy(entry)
    end

    return buses_by_name, buses_by_id
end

# Plot voltage magnitude along feeder (distance vs volts)
function plot_voltage_along_feeder_snap(buses_dict::Dict{String,Dict{String,Any}}, lines_df::DataFrame;
        t::Int=1, Vthreshold::Float64=1000.0, vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0)

    p = plot(legend=false)
    xlabel!("Distance from reference bus (km)")
    ylabel!("Voltage magnitude P-N (V)")
    title!("Voltage along feeder")

    colors = Dict(1=>:blue, 2=>:red, 3=>:black)

    for r in eachrow(lines_df)
        b1 = r.Bus1
        b2 = r.Bus2
        phases = r.phases

        if !(haskey(buses_dict, b1) && haskey(buses_dict, b2))
            continue
        end
        if !(haskey(buses_dict[b1], "distance") && haskey(buses_dict[b2], "distance"))
            continue
        end

        for ph in phases
            vm_f = ph == 1 ? buses_dict[b1]["vma"][t] :
                   ph == 2 ? buses_dict[b1]["vmb"][t] :
                             buses_dict[b1]["vmc"][t]

            vm_t = ph == 1 ? buses_dict[b2]["vma"][t] :
                   ph == 2 ? buses_dict[b2]["vmb"][t] :
                             buses_dict[b2]["vmc"][t]

            if vm_f < Vthreshold && vm_t < Vthreshold
                plot!(
                    [buses_dict[b1]["distance"], buses_dict[b2]["distance"]],
                    [vm_f, vm_t],
                    color=colors[ph],
                    marker=:circle,
                    markersize=1
                )
            end
        end
    end

    maxdist = maximum(bus["distance"] for bus in values(buses_dict) if haskey(bus, "distance"))
    plot!([0, maxdist], [vminV, vminV], linestyle=:dash, color=:red)
    plot!([0, maxdist], [vmaxV, vmaxV], linestyle=:dash, color=:red)

    return p
end

# Plot voltage histogram at a snapshot
function plot_voltage_histogram_snap(buses_dict::Dict{String,Dict{String,Any}};
        t::Int=1, Vthreshold::Float64=1000.0, vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0)

    phase_a = Float64[]
    phase_b = Float64[]
    phase_c = Float64[]

    for (_, bus_data) in buses_dict
        if haskey(bus_data, "vma") && bus_data["vma"][t] < Vthreshold
            push!(phase_a, bus_data["vma"][t])
        end
        if haskey(bus_data, "vmb") && bus_data["vmb"][t] < Vthreshold
            push!(phase_b, bus_data["vmb"][t])
        end
        if haskey(bus_data, "vmc") && bus_data["vmc"][t] < Vthreshold
            push!(phase_c, bus_data["vmc"][t])
        end
    end

    bins = (vminV-5):0.5:(vmaxV+5)

    p = histogram(phase_a; bins=bins, label="phase a")
    histogram!(phase_b; bins=bins, label="phase b")
    histogram!(phase_c; bins=bins, label="phase c")

    xlabel!("Voltage magnitude (V)")
    ylabel!("Count")
    title!("Voltage histogram")

    return p
end

# Combine voltage profile and histogram in one figure
function plot_voltage_combined_snap(buses_dict::Dict{String,Dict{String,Any}}, lines_df::DataFrame;
        t::Int=1, Vthreshold::Float64=1000.0, vminV::Float64=0.94*230.0, vmaxV::Float64=1.10*230.0)

    p1 = plot_voltage_along_feeder_snap(buses_dict, lines_df; t=t, Vthreshold=Vthreshold, vminV=vminV, vmaxV=vmaxV)
    p2 = plot_voltage_histogram_snap(buses_dict; t=t, Vthreshold=Vthreshold, vminV=vminV, vmaxV=vmaxV)
    return plot(p1, p2, layout=(1,2))
end

# Load bus coordinate CSV (if available) and plot topology
function plot_topology_if_coords(lines_df::DataFrame, coords_csv::String; outpath::String)
    isfile(coords_csv) || return nothing

    coords = CSV.read(coords_csv, DataFrame)

    # Try common column names used by OpenDSS exports
    # If these differ, rename here based on actual CSV header
    col_bus = findfirst(n -> lowercase(String(n)) in ["bus", "busname", "name"], names(coords))
    col_x   = findfirst(n -> lowercase(String(n)) in ["x", "xm", "easting"], names(coords))
    col_y   = findfirst(n -> lowercase(String(n)) in ["y", "ym", "northing"], names(coords))

    (col_bus === nothing || col_x === nothing || col_y === nothing) && return nothing

    bus_col = String(names(coords)[col_bus])
    x_col   = String(names(coords)[col_x])
    y_col   = String(names(coords)[col_y])

    xy = Dict{String,Tuple{Float64,Float64}}()
    for r in eachrow(coords)
        b = string(r[bus_col])
        x = Float64(r[x_col])
        y = Float64(r[y_col])
        xy[b] = (x, y)
    end

    p = plot(legend=false, xlabel="x", ylabel="y", title="Network topology")

    for r in eachrow(lines_df)
        b1 = r.Bus1
        b2 = r.Bus2
        if haskey(xy, b1) && haskey(xy, b2)
            x1, y1 = xy[b1]
            x2, y2 = xy[b2]
            plot!([x1, x2], [y1, y2])
        end
    end

    savefig(p, outpath)
    return p
end

# --------------------------------------------------
# 5) Read one customer profile and turn it into alpha(t)
# --------------------------------------------------

parquet_file = find_parquet_for_year(PROFILES_ROOT, YEAR)
println("Parquet file: ", parquet_file)

con = DBInterface.connect(DuckDB.DB, ":memory:")

nrows_df = DBInterface.execute(con, "SELECT COUNT(*) AS n FROM read_parquet('$(parquet_file)')") |> DataFrame
nrows = Int(nrows_df.n[1])
println("Profile rows: ", nrows)

schema_df = DBInterface.execute(con, "DESCRIBE SELECT * FROM read_parquet('$(parquet_file)')") |> DataFrame
numcols = numeric_columns(schema_df)

customer_col =
    CUSTOMER_COL != "AUTO" ? CUSTOMER_COL :
    pick_customer_column_auto(con, parquet_file, numcols)

println("Chosen customer column: ", customer_col)

prof_df = DBInterface.execute(con, "SELECT \"$(customer_col)\" AS y FROM read_parquet('$(parquet_file)')") |> DataFrame
DBInterface.close!(con)

y = Float64.(prof_df.y)

# alpha(t) is a normalized version of the customer profile
# mean alpha is about 1 before clipping
y_mean = mean(y)
alpha = y ./ y_mean

# clip alpha for stability in early runs
alpha = clamp.(alpha, ALPHA_MIN, ALPHA_MAX)

t = build_time_axis(nrows; start_date=START_DATE, step_minutes=STEP_MINUTES)

# --------------------------------------------------
# 6) Parse feeder once and run baseline diagnostic
# --------------------------------------------------

println("Parsing feeder: ", master_dss)
eng0 = PMD.parse_file(master_dss, transformations=[PMD.transform_loops!])

# Apply base scaling if needed
scale_loads!(eng0, LOAD_ALPHA_BASE)

# Print baseline load magnitude so it is obvious what the solver is seeing
print_baseline_load_summary(eng0)

# Prepare line data and distances for voltage profile plots
lines_df = make_lines_df_from_eng(eng0)
source_bus = pick_source_bus(eng0)
dist = compute_bus_distances(lines_df; source_bus=source_bus)
println("\nDistance reference bus: ", source_bus)

# Optional topology plot if coordinate file exists
topo_out = joinpath(FIGDIR, "topology.png")
plot_topology_if_coords(lines_df, buscoords_csv; outpath=topo_out)

# --------------------------------------------------
# 7) Time loop: apply alpha(t), run PF, store metrics
# --------------------------------------------------

rows = NamedTuple[]

for k in 1:STRIDE:nrows
    # fresh copy per timestep, so scaling does not accumulate
    eng = deepcopy(eng0)

    # apply time-varying scaling for this timestep
    scale_loads!(eng, alpha[k])

    # solve PF and extract metrics
    pf, math = solve_pf(eng)
    r = pf_metrics(pf; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

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

    # occasional progress print
    if k % (STRIDE * 200) == 1
        println("Progress ", k, " / ", nrows, " | vmin=", r.vmin, " vmax=", r.vmax)
    end
end

df = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "timeseries_baseline_pf_metrics.csv"), df)

# --------------------------------------------------
# 8) Basic time-series plots (system level)
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

# --------------------------------------------------
# 9) Worst timestep voltage profile plots
# --------------------------------------------------

# Pick worst timesteps by lowest vmin_pu
df_sorted = sort(df, :vmin_pu)
nplot = min(N_WORST_PLOTS, nrow(df_sorted))

println("\nCreating voltage profile plots for worst timesteps:")
println("  count = ", nplot)

for i in 1:nplot
    row = df_sorted[i, :]
    k = Int(row.timestep)
    tt = row.time
    a = Float64(row.alpha_t)

    # rebuild the PF case for this timestep
    eng = deepcopy(eng0)
    scale_loads!(eng, a)

    pf, math = solve_pf(eng)

    # extract bus voltages in volts
    buses_by_name, buses_by_id = solved_bus_vm_volts_dual(pf, math; vbase_ln=VBASE_LN)

    # Attach distances using whichever keying matches the dist dictionary best
    hits_id = count(k -> haskey(buses_by_id, k), keys(dist))
    hits_name = count(k -> haskey(buses_by_name, k), keys(dist))

    buses_dict = hits_id >= hits_name ? buses_by_id : buses_by_name

    for (bus, d) in dist
        if haskey(buses_dict, bus)
            buses_dict[bus]["distance"] = d
        end
    end

println("Distance mapping hits: by_id=", hits_id, " by_name=", hits_name)


    vminV = VMIN_PU * VBASE_LN
    vmaxV = VMAX_PU * VBASE_LN

    p_comb = plot_voltage_combined_snap(
        buses_dict,
        lines_df;
        t=1,
        Vthreshold=1000.0,
        vminV=vminV,
        vmaxV=vmaxV
    )

    # filename includes rank and timestep
    fname = "worst_voltage_combined_rank=$(i)_timestep=$(k).png"
    savefig(p_comb, joinpath(FIGDIR, fname))

    println("  saved rank=", i, " timestep=", k, " time=", tt, " alpha_t=", a, " vmin_pu=", row.vmin_pu)
end

# --------------------------------------------------
# 10) Small summary printout
# --------------------------------------------------

println("\nSummary checks")
println("alpha_t min = ", minimum(df.alpha_t), " | max = ", maximum(df.alpha_t))
println("alpha hits ALPHA_MIN count = ", count(==(ALPHA_MIN), df.alpha_t))
println("alpha hits ALPHA_MAX count = ", count(==(ALPHA_MAX), df.alpha_t))

println("vmin over run: min = ", minimum(df.vmin_pu), " | median = ", median(df.vmin_pu), " | max = ", maximum(df.vmin_pu))
println("vmax over run: min = ", minimum(df.vmax_pu), " | median = ", median(df.vmax_pu), " | max = ", maximum(df.vmax_pu))

println("timesteps with any under-voltage buses = ", count(>(0), df.n_under))
println("timesteps with any over-voltage buses  = ", count(>(0), df.n_over))

println("\nSaved results to: ", OUTDIR)
println("Done.")
