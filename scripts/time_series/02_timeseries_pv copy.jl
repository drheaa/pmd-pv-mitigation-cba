# Time-series comparison on selected PV-stress timesteps and a
# contiguous 48-hour window for paper-style current plots:
#
#   (1) Baseline (NO PV, NO STATCOM)         -> master_baseline_4w.dss
#   (2) PV only (PV as negative load)        -> master_baseline_4w.dss
#   (3) PV + STATCOM (OPF w inverter gens)   -> master_scaled.dss
#
# Non-negotiables enforced:
# - Baseline and PV never parse master_scaled.dss
# - STATCOM scenario only parses master_scaled.dss
# - STATCOM integration path remains OPF + inverter gens workflow
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
using Ipopt
using JuMP
using LinearAlgebra
using Random

const PMD = PowerModelsDistribution
PMD.silence!()

# --------------------------------------------------
# 0) Settings
# --------------------------------------------------

ROOT = joinpath(@__DIR__, "../..")

# Required helper functions (non-negotiable STATCOM workflow support)
include(joinpath(ROOT, "src/read_functions.jl"))

NET_4W = "spd_s_4w"

MASTER_BASELINE = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w.dss")
MASTER_STATCOM  = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_scaled.dss")

YEAR   = 2023
STRIDE = 4

BASELINE_OUTDIR = joinpath(ROOT, "results", "time_series", "baseline_pf", NET_4W, "year=$(YEAR)_stride=$(STRIDE)_K=17520")
BASELINE_CSV    = joinpath(BASELINE_OUTDIR, "tables", "timeseries_baseline_pf_metrics.csv")

# PV settings
PV_BUS          = "AUTO"
PV_KW_PER_PHASE = 10.0          # increase if PV impact remains too small
PV_PHASES       = [1]           # single-phase PV on phase A by default

# PV-stress
N_PV_STRESS    = 30
SKIP_NIGHT_PV  = true
RANDOM_SEED    = 42

# 48-hour window for paper-style current plots (Fig6/Fig7/Fig9)
WINDOW_HOURS   = 48
WINDOW_STEPS   = Int(WINDOW_HOURS * 60 ÷ 30) + 1   # 30-min baseline sampling
WINDOW_CENTER  = "BEST_STRESS"  # "BEST_STRESS" or "FIRST_STRESS"

# Limits and conversion
VMIN_PU   = 0.90
VMAX_PU   = 1.10
VBASE_LN  = 230.0

# Baseline engineering scaling applied once, then alpha(t)
LOAD_ALPHA_BASE = 1.5    # scale baseline loads by this factor once at parsing and should be between 1.5 to 2.0

# Solver
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

OUTDIR = joinpath(
    ROOT, "results", "time_series", "pv_pf",
    NET_4W,
    "year=$(YEAR)_stride=$(STRIDE)",
    "pvkw=$(PV_KW_PER_PHASE)_ph=$(join(PV_PHASES, '_'))_K=$(N_PV_STRESS)"
)
FIGDIR = joinpath(OUTDIR, "figures")
TBLDIR = joinpath(OUTDIR, "tables")
mkpath(FIGDIR); mkpath(TBLDIR)

## --------------------------------------------------
# 1) Helpers: scaling, metrics, PV proxy, PV injection
# --------------------------------------------------

function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng     # early exit if no loads
    for (_, ld_any) in eng["load"]        # loop through every load object
        ld = ld_any::Dict{String,Any}     # treat each load as a dictionary
        for k in ("pd_nom", "qd_nom", "pd", "qd")   # for each power-related key, scale if it exists
            if haskey(ld, k)     # only scale if the load actually has that field
                v = ld[k]        # get the current value
                ld[k] = v isa Number ? alpha * v : alpha .* v   # scale numbers vs vectors correctly
            end
        end
    end
    return eng     # return the mutated dictionary
end

function add_load_q_from_pf!(eng::Dict{String,Any}; pf::Float64=0.95)    # add reactive power to loads based on specified power factor    
    haskey(eng, "load") || return eng
    pf = clamp(pf, 0.5, 0.999999)    # avoid invalid power factor values
    q_factor = tan(acos(pf))  # Q = |P| * tan(arccos(pf))
    println("DEBUG add_load_q_from_pf!: pf=", pf, " q_factor=", tan(acos(pf)))

    for (_, ld_any) in eng["load"]      # loop through every load object
        ld = ld_any::Dict{String,Any}

        # Choose the active P key (pd preferred, otherwise pd_nom)
        pkey = haskey(ld, "pd") ? "pd" : (haskey(ld, "pd_nom") ? "pd_nom" : nothing)
        pkey === nothing && continue

        # Write matching Q key
        qkey = (pkey == "pd") ? "qd" : "qd_nom"

        P = ld[pkey]
        if P isa Number    # scalar active power
            ld[qkey] = abs(Float64(P)) * q_factor   # set reactive power accordingly
        elseif P isa AbstractVector
            Pvec = Float64.(P)
            ld[qkey] = abs.(Pvec) .* q_factor    # set reactive power vector accordingly
        end
    end

    return eng
end

function sum_loads_eng(eng::Dict{String,Any}) 
    # sum up loads in kW and kVAR from eng dictionary
    P = 0.0
    Q = 0.0

    for ld_any in values(get(eng, "load", Dict{String,Any}()))  # loop through each load
        ld = ld_any::Dict{String,Any}

        # prefer *_nom in eng
        if haskey(ld, "pd_nom")
            v = ld["pd_nom"]
            P += v isa Number ? Float64(v) : sum(Float64.(v))   # sum up kW
        elseif haskey(ld, "pd")
            v = ld["pd"]
            P += v isa Number ? Float64(v) : sum(Float64.(v))   # sum up kW
        end

        if haskey(ld, "qd_nom")
            v = ld["qd_nom"]
            Q += v isa Number ? Float64(v) : sum(Float64.(v))   # sum up kVAR
        elseif haskey(ld, "qd")
            v = ld["qd"]
            Q += v isa Number ? Float64(v) : sum(Float64.(v))   # sum up kVAR
        end
    end

    return (P_kW = P, Q_kvar = Q)
end

function pf_metrics(pf::Dict{String,Any}; vmin_pu=0.90, vmax_pu=1.10)
    sol_bus = pf["solution"]["bus"]     # grab the solved bus results
    vmins = Float64[]
    vmaxs = Float64[]      # prepare to collect min/max voltages
    
    for sb in values(sol_bus)      # loop through each bus solution
        vm =
            haskey(sb, "vm") ? sb["vm"][1:min(3, length(sb["vm"]))] :   # try vm field first if it exists and use that
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"][1:min(3, length(sb["vr"]))].^2 .+ sb["vi"][1:min(3, length(sb["vi"]))].^2) :   # otherwise compute from vr and vi
            nothing    # if neither exists, skip this bus
        vm === nothing && continue  # skip bus if no voltage info available
        push!(vmins, minimum(vm))  
        push!(vmaxs, maximum(vm))   # store this bus’s min/max across phases
    end

    if isempty(vmins)    # no voltage data found at all 
        return (status=string(pf["termination_status"]), vmin=NaN, vmax=NaN, n_under=0, n_over=0)  # early exit
    end

    vmin_sys = minimum(vmins)
    vmax_sys = maximum(vmaxs)    # system-wide min/max voltages
    n_under = count(x -> x < vmin_pu, vmins)
    n_over  = count(x -> x > vmax_pu, vmaxs)    # count how many buses violate limits

    return (status=string(pf["termination_status"]), vmin=vmin_sys, vmax=vmax_sys, n_under=n_under, n_over=n_over)   # return metrics as a named tuple
end

function pv_shape_simple(time_vec::Vector{DateTime})   # simple normalized PV shape over the day
    pv = zeros(Float64, length(time_vec))              # prepare output vector
    for (i, ts) in enumerate(time_vec)                 # loop through each timestamp
        h = hour(ts) + minute(ts)/60                   # compute decimal hour
        x = (h - 6.0)/12.0                             # map 6am-6pm to 0-1
        pv[i] = (0.0 <= x <= 1.0) ? sin(pi*x)^2 : 0.0  # simple sin^2 shape between 6am and 6pm
    end
    m = maximum(pv)                                    # normalize to peak of 1
    return m > 0 ? pv ./ m : pv                        # return normalized PV shape
end

function add_pv_negative_load!(eng::Dict{String,Any}, pv_bus::String, pv_kw_per_phase::Float64, pv_phases::Vector{Int})
    haskey(eng, "load") || (eng["load"] = Dict{String,Any}())
    
    # modifying an existing load at pv_bus
    for (lid, ld_any) in eng["load"]   
        ld = ld_any::Dict{String,Any}
        haskey(ld, "bus") || continue
        string(ld["bus"]) == pv_bus || continue    # only consider loads at the PV bus and skip loads without a bus, or not on the PV bus

        pd_key = haskey(ld, "pd") ? "pd" : (haskey(ld, "pd_nom") ? "pd_nom" : nothing)   # find pd or pd_nom key
        pd_key === nothing && continue    # skip loads without pd or pd_nom

        pd_val = ld[pd_key]   # get the current pd value

        if pd_val isa Number
            ld[pd_key] = Float64(pd_val) - pv_kw_per_phase * length(pv_phases)    # subtract total PV from scalar load
            return ("modified_existing_load", string(lid))
        end

        if pd_val isa AbstractVector
            pd_new = Float64.(pd_val)   # make a copy to modify
            for ph in pv_phases    # subtract PV from specified phases
                if 1 <= ph <= length(pd_new)
                    pd_new[ph] -= pv_kw_per_phase   # subtract PV from this phase
                end
            end
            ld[pd_key] = pd_new   # update the load with modified vector
            return ("modified_existing_load", string(lid))
        end
    end

    # Fallback: create a new negative load
    load_dict_template = eng["load"][first(keys(eng["load"]))]::Dict{String,Any}   # get a template load dictionary
    new_id = "pv_neg_load_$(pv_bus)"    # unique load ID
    eng["load"][new_id] = load_dict_template
    eng["load"][new_id]["bus"] = pv_bus
    eng["load"][new_id]["connections"] = pv_phases
    eng["load"][new_id]["phases"] = length(pv_phases)
    eng["load"][new_id]["pd"] = fill(-pv_kw_per_phase, length(pv_phases))
    eng["load"][new_id]["qd"] = fill(0.0, length(pv_phases))
    
    @show eng["load"][new_id]
    return ("created_new_load", new_id)
end

function kw_to_w!(eng::Dict{String,Any})    # convert all kW/kVar load ratings to W/Var
    haskey(eng, "load") || return eng      # early exit if no loads
    for ld_any in values(eng["load"])      # loop through every load object
        ld = ld_any::Dict{String,Any}      # treat each load as a dictionary

        for k in ("pd_nom", "qd_nom", "pd", "qd")       # for each power-related key 
            haskey(ld, k) || continue
            v = ld[k]                                   # get the current value
            ld[k] = v isa Number ? 1000.0 * Float64(v) : 1000.0 .* Float64.(v)     # convert numbers vs vectors correctly
        end
    end
    return eng
end

# --------------------------------------------------
# 2) Helpers: PV bus AUTO selection (distance-based)
# --------------------------------------------------

function make_lines_df_from_eng(eng::Dict{String,Any})               #extract lines/branches into a DataFrame
    rows = NamedTuple[]                                              
    if haskey(eng, "line")
        for (_, ln_any) in eng["line"]                               # loop through each line
            ln = ln_any::Dict{String,Any}
            push!(rows, (Bus1=string(ln["f_bus"]), Bus2=string(ln["t_bus"]),    # extract from and to buses
                         length_km=(get(ln, "length", 0.0) / 1000.0)))          # extract length in km (default 0.0 if missing)
        end
    elseif haskey(eng, "branch")                          # fallback to branches if no lines
        for (_, br_any) in eng["branch"]                  # loop through each branch
            br = br_any::Dict{String,Any}
            push!(rows, (Bus1=string(br["f_bus"]), Bus2=string(br["t_bus"]),          # extract from and to buses
                         length_km=(get(br, "length", 0.0) / 1000.0)))                # extract length in km (default 0.0 if missing)
        end
    else
        error("No line or branch data found")
    end
    return DataFrame(rows)
end

function pick_source_bus_name(eng::Dict{String,Any})       # pick source bus name from common conventions
    if haskey(eng, "bus")
        names = Set(string.(keys(eng["bus"])))             # collect all bus names
        if "sourcebusz" in names                
            return "sourcebusz"
        elseif "sourcebus" in names
            return "sourcebus"
        end
    end
    return "sourcebus"
end

function compute_bus_distances(lines_df::DataFrame; source_bus::String)         # compute shortest-path distances from source bus though Breadth-First Search
    adj = Dict{String, Vector{Tuple{String, Float64}}}()                        # adjacency list of the network which is undirected and shows connections between buses
    for r in eachrow(lines_df)
        push!(get!(adj, r.Bus1, Tuple{String,Float64}[]), (r.Bus2, r.length_km))         
        push!(get!(adj, r.Bus2, Tuple{String,Float64}[]), (r.Bus1, r.length_km))      # undirected graph
    end

    dist = Dict{String,Float64}(source_bus => 0.0)                 # distances from source bus
    queue = [source_bus]                       # BFS queue
    while !isempty(queue)                      # BFS loop
        u = popfirst!(queue)                   # dequeue
        for (v, w) in get(adj, u, Tuple{String,Float64}[])         # explore neighbors
            if !haskey(dist, v)
                dist[v] = dist[u] + w                              # update distance to neighbor
                push!(queue, v)                                    # enqueue neighbor
            end
        end
    end
    return dist
end

function pick_pv_bus_auto(eng::Dict{String,Any}, dist::Dict{String,Float64})       # pick PV bus automatically based on distance from source bus, avoiding common source bus names
    bad = Set(["sourcebus", "sourcebusz", "SourceBus", "SourceBusZ"])              # bad bus names to avoid
    best_bus = ""; best_d = -Inf                                                   # initialize best bus and distance

    for (b, d) in dist
        if lowercase(b) in lowercase.(collect(bad))         
            continue
        end
        if d > best_d
            best_d = d
            best_bus = b
        end
    end

    best_bus != "" && return best_bus                    # return best bus if found

    for b in keys(eng["bus"])                            # fallback: first non-bad bus
        if lowercase(string(b)) in lowercase.(collect(bad))           
            continue
        end
        return string(b)
    end

    return first(keys(eng["bus"])) |> string              # last resort: first bus in the dictionary
end

# --------------------------------------------------
# 3) Helpers: robust feeder-head current extraction and metrics
# --------------------------------------------------

function current_base_A(math::Dict{String,Any}; vbase_ln::Float64)              # compute current base in amps from math dictionary and line-to-neutral voltage base
    sbase = get(get(math, "settings", Dict{String,Any}()), "sbase_default", 1.0)    # default to 1.0 if missing
    S = Float64(sbase)                                                              
    # sbase_default in PMD math is in VA (commonly), but can be W scaling-based.
    # The existing heuristic converts only when magnitudes look per-unit.
    return S / (3.0 * vbase_ln)                                            # Ibase = Sbase / (sqrt(3) * Vbase_LL) = Sbase / (3 * Vbase_LN)
end

function _phasor_from(br::Dict{String,Any}, rkey::String, ikey::String, idx::Int)    #extract a complex phasor from real/imag keys at given phase index
    (haskey(br, rkey) && haskey(br, ikey)) || return nothing
    r = br[rkey]; im = br[ikey]                               # get real and imaginary parts
    (r isa AbstractVector && im isa AbstractVector) || return nothing     # ensure both are vectors
    (length(r) >= idx && length(im) >= idx) || return nothing             # ensure index is valid
    return complex(Float64(r[idx]), Float64(im[idx]))          # return complex phasor
end

function _seq_components(Aa::Complex, Ab::Complex, Ac::Complex)         # compute sequence components from phase phasors
    a = cis(2pi / 3)                                  # operator a = e^(j120°)
    A0 = (Aa + Ab + Ac) / 3                    # zero-sequence
    A1 = (Aa + a * Ab + a^2 * Ac) / 3          # positive-sequence
    A2 = (Aa + a^2 * Ab + a * Ac) / 3          # negative-sequence
    return (A0=A0, A1=A1, A2=A2)
end

function _branch_current_phasors(br::Dict{String,Any})         # extract branch current phasors robustly from various possible key namings
    keypairs = [
        ("cr_fr", "ci_fr"),
        ("cfr",   "cfi"),
        ("cr",    "ci"),
        ("ctr",   "cti"),
    ]                                           # possible key name pairs for real and imaginary parts of current
    for (rk, ik) in keypairs
        Ia = _phasor_from(br, rk, ik, 1)        # extract phase A current
        Ib = _phasor_from(br, rk, ik, 2)        # extract phase B current
        Ic = _phasor_from(br, rk, ik, 3)        # extract phase C current
        In = _phasor_from(br, rk, ik, 4)        # extract neutral current (if present)
        if Ia !== nothing && Ib !== nothing && Ic !== nothing
            return (Ia=Ia, Ib=Ib, Ic=Ic, In=In, rk=rk, ik=ik)        # return phasors if successfully extracted
        end
    end
    return nothing
end

function _find_source_bus_id(math::Dict{String,Any}, source_bus_name::String)      # find source bus ID from its name
    haskey(math, "bus") || return nothing
    for (bid, b_any) in math["bus"]                # loop through all buses
        b = b_any::Dict{String,Any}
        if haskey(b, "name") && string(b["name"]) == source_bus_name       # match by name
            return bid
        end
    end
    return nothing
end

function _pick_head_branch_id(pf::Dict{String,Any}, math::Dict{String,Any}, source_bus_name::String)   # pick the feeder head branch ID based on priorities
    sol = pf["solution"]
    haskey(sol, "branch") || return nothing
    isempty(sol["branch"]) && return nothing

    # Priority 1: the explicit source impedance line if present
    if haskey(math, "branch")
        for (br_id, br_any) in math["branch"]           # loop through all branches
            brm = br_any::Dict{String,Any}
            sid = get(brm, "source_id", "")             # get source_id field if it exists
            if occursin("line.sourceZ", string(sid)) || occursin("sourceZ", lowercase(string(sid)))         # check for sourceZ indication
                if haskey(sol["branch"], br_id)         # check if this branch is in the solution
                    return br_id
                elseif haskey(sol["branch"], string(br_id))   # check string version
                    return string(br_id)
                end
            end
        end
    end

    # Priority 2: branch connected to the source bus with largest 3-phase current magnitude
    src_id = _find_source_bus_id(math, source_bus_name)       # find source bus ID
    best = nothing   
    best_mag = -Inf

    if src_id !== nothing && haskey(math, "branch")
        for (br_id, br_any) in math["branch"]
            brm = br_any::Dict{String,Any}
            connected = (brm["f_bus"] == src_id || brm["t_bus"] == src_id)     # check if branch is connected to source bus
            connected || continue

            sol_id = haskey(sol["branch"], br_id) ? br_id : (haskey(sol["branch"], string(br_id)) ? string(br_id) : nothing)   # find corresponding solution branch ID
            sol_id === nothing && continue

            brs = sol["branch"][sol_id]         # get branch solution
            ph = _branch_current_phasors(brs)   # extract current phasors
            ph === nothing && continue          # skip if unable to extract currents

            mag = maximum([abs(ph.Ia), abs(ph.Ib), abs(ph.Ic)])    # compute largest phase current magnitude
            if mag > best_mag                  # update best if larger
                best_mag = mag  
                best = sol_id
            end
        end
        best !== nothing && return best
    end

    # Priority 3: any branch with largest 3-phase current magnitude
    for (br_id, brs_any) in sol["branch"]         # loop through all branches in solution
        brs = brs_any::Dict{String,Any}
        ph = _branch_current_phasors(brs)         # extract current phasors
        ph === nothing && continue                # skip if unable to extract currents
        mag = maximum([abs(ph.Ia), abs(ph.Ib), abs(ph.Ic)])    # compute largest phase current magnitude
        if mag > best_mag               # update best if larger
            best_mag = mag
            best = br_id
        end
    end

    return best === nothing ? first(keys(sol["branch"])) : best
end

function feeder_head_current_metrics(pf::Dict{String,Any}, math::Dict{String,Any}, source_bus_name::String; vbase_ln::Float64=230.0)    # extract feeder head current metrics robustly
    sol = pf["solution"]
    if !haskey(sol, "branch") || isempty(sol["branch"])
        return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)   # early exit if no branch solution
    end

    bid = _pick_head_branch_id(pf, math, source_bus_name)
    bid === nothing && return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)  # early exit if unable to pick branch ID

    br = sol["branch"][bid]
    ph = _branch_current_phasors(br)   # extract current phasors
    ph === nothing && return (Ia=NaN, Ib=NaN, Ic=NaN, In=NaN, Imean=NaN, Imax=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)   # early exit if unable to extract currents

    Ia, Ib, Ic, In = ph.Ia, ph.Ib, ph.Ic, ph.In   # get phasors
    ma, mb, mc = abs(Ia), abs(Ib), abs(Ic)   # magnitudes
    mn = (In === nothing) ? NaN : abs(In)    # neutral magnitude if present

    # Per-unit to amps conversion heuristic
    if maximum([ma, mb, mc]) < 10.0
        IbaseA = current_base_A(math; vbase_ln=vbase_ln)   # compute current base in amps
        Ia *= IbaseA; Ib *= IbaseA; Ic *= IbaseA           # scale phasors to amps
        ma *= IbaseA; mb *= IbaseA; mc *= IbaseA           # scale magnitudes to amps
        if isfinite(mn)     # scale neutral if present
            mn *= IbaseA
        end
        if In !== nothing   # scale neutral phasor if present
            In *= IbaseA
        end
    end

    seq = _seq_components(Ia, Ib, Ic)     # compute sequence components
    I0 = abs(seq.A0)                      # zero-sequence magnitude
    I1 = abs(seq.A1)                      # positive-sequence magnitude
    I2 = abs(seq.A2)                      # negative-sequence magnitude
    r  = I1 > 1e-12 ? I2 / I1 : NaN       # negative-to-positive sequence ratio

    return (Ia=ma, Ib=mb, Ic=mc, In=mn, Imean=(ma+mb+mc)/3, Imax=max(ma,mb,mc), I0=I0, I2=I2, I2_over_I1=r)    # return metrics as named tuple
end

function sum_loads_math_pu(math::Dict{String,Any})   # sum up loads in per-unit from math dictionary
    Ppu = 0.0; Qpu = 0.0    # initialize sums
    for ld_any in values(get(math, "load", Dict{String,Any}()))    # loop through each load
        ld = ld_any::Dict{String,Any}

        if haskey(ld, "pd")
            v = ld["pd"]
            Ppu += v isa Number ? Float64(v) : sum(Float64.(v))     # sum up per-unit active power
        end
        if haskey(ld, "qd")
            v = ld["qd"]
            Qpu += v isa Number ? Float64(v) : sum(Float64.(v))     # sum up per-unit reactive power
        end
    end
    return (P_pu=Ppu, Q_pu=Qpu)     # return sums as named tuple
end

function sbase_math(math::Dict{String,Any})   # extract sbase from math dictionary
    # common places sbase appears
    if haskey(math, "settings") && haskey(math["settings"], "sbase_default")
        return Float64(math["settings"]["sbase_default"])    # return sbase_default if present
    end
    if haskey(math, "baseMVA")
        return Float64(math["baseMVA"])     # return baseMVA if present
    end
    return NaN
end

# --------------------------------------------------
# 4) Helpers: voltage sequence extraction for all buses
# --------------------------------------------------

function bus_voltage_phasors(pf::Dict{String,Any})
    sol_bus = pf["solution"]["bus"]
    out = Dict{String, Tuple{ComplexF64,ComplexF64,ComplexF64}}()

    for (bid, sb_any) in sol_bus
        sb = sb_any::Dict{String,Any}
        if haskey(sb, "vr") && haskey(sb, "vi")
            vr = sb["vr"]; vi = sb["vi"]
            if vr isa AbstractVector && vi isa AbstractVector && length(vr) >= 3 && length(vi) >= 3
                Va = complex(Float64(vr[1]), Float64(vi[1]))
                Vb = complex(Float64(vr[2]), Float64(vi[2]))
                Vc = complex(Float64(vr[3]), Float64(vi[3]))
                out[string(bid)] = (Va, Vb, Vc)
            end
        elseif haskey(sb, "vm") && haskey(sb, "va")
            vm = sb["vm"]; va = sb["va"]
            if vm isa AbstractVector && va isa AbstractVector && length(vm) >= 3 && length(va) >= 3
                Va = Float64(vm[1]) * cis(Float64(va[1]))
                Vb = Float64(vm[2]) * cis(Float64(va[2]))
                Vc = Float64(vm[3]) * cis(Float64(va[3]))
                out[string(bid)] = (Va, Vb, Vc)
            end
        end
    end

    return out
end

function voltage_sequence_rows(pf::Dict{String,Any}, time::DateTime, scenario::String)
    V = bus_voltage_phasors(pf)
    rows = NamedTuple[]
    for (bus_id, (Va,Vb,Vc)) in V
        seq = _seq_components(Va,Vb,Vc)
        V0 = abs(seq.A0)
        V1 = abs(seq.A1)
        V2 = abs(seq.A2)
        r  = V1 > 1e-12 ? (V2 / V1) : NaN
        push!(rows, (time=time, scenario=scenario, bus_id=bus_id, V0=V0, V2=V2, V2_over_V1=r))
    end
    return rows
end

# --------------------------------------------------
# 5) Helpers: plotting utilities (paper-style)
# --------------------------------------------------

function duration_curve(values::AbstractVector{<:Real}; high_is_worse::Bool=true)
    v = Float64.(values)
    v = filter(isfinite, v)
    isempty(v) && return (p=Float64[], x=Float64[])
    x = sort(v; rev=high_is_worse)
    p = collect(range(0, 100; length=length(x)))
    return (p=p, x=x)
end

function plot_inputs(df_sel::DataFrame, outpath::String)
    p1 = plot(df_sel.time, df_sel.alpha_t; xlabel="Time", ylabel="alpha(t)", label="alpha(t)",
              title="Input load scaling alpha(t)")
    p2 = plot(df_sel.time, df_sel.pv_pu; xlabel="Time", ylabel="pv_pu_proxy", label="pv_pu_proxy",
              title="Input PV proxy pv_pu_proxy(t)")
    p3 = plot(df_sel.time, df_sel.pv_kw_eff; xlabel="Time", ylabel="PV kW per phase", label="pv_kw_eff",
              title="Effective PV injection per phase (kW)")
    plt = plot(p1, p2, p3; layout=(3,1), size=(900,900))
    savefig(plt, outpath)
end

function fig6_like_currents(df::DataFrame, prefix::String, outpath::String)
    p = plot(df.time, df[!, Symbol(prefix*"_Ia")]; label="phase a", xlabel="Time", ylabel="Current (A)",
             title="Substation phase currents | $(prefix)")
    plot!(df.time, df[!, Symbol(prefix*"_Ib")]; label="phase b")
    plot!(df.time, df[!, Symbol(prefix*"_Ic")]; label="phase c")
    if hasproperty(df, Symbol(prefix*"_In"))
        plot!(df.time, df[!, Symbol(prefix*"_In")]; label="neutral")
    end
    plot!(df.time, df[!, Symbol(prefix*"_Imean")]; label="mean", linewidth=2)
    savefig(p, outpath)
end

function fig7_like_range_band(df::DataFrame, base_prefix::String, after_prefix::String, outpath::String)
    base_min = min.(df[!, Symbol(base_prefix*"_Ia")], df[!, Symbol(base_prefix*"_Ib")], df[!, Symbol(base_prefix*"_Ic")])
    base_max = max.(df[!, Symbol(base_prefix*"_Ia")], df[!, Symbol(base_prefix*"_Ib")], df[!, Symbol(base_prefix*"_Ic")])
    after_min = min.(df[!, Symbol(after_prefix*"_Ia")], df[!, Symbol(after_prefix*"_Ib")], df[!, Symbol(after_prefix*"_Ic")])
    after_max = max.(df[!, Symbol(after_prefix*"_Ia")], df[!, Symbol(after_prefix*"_Ib")], df[!, Symbol(after_prefix*"_Ic")])

    p = plot(df.time, df[!, Symbol(base_prefix*"_Imean")]; label="Mean before", xlabel="Time", ylabel="Substation current (A)",
             title="Substation current range band | before vs after")
    plot!(df.time, df[!, Symbol(after_prefix*"_Imean")]; label="Mean after", linewidth=2)

    plot!(df.time, base_min; label="Before (min)", linestyle=:dot)
    plot!(df.time, base_max; label="Before (max)", linestyle=:dot)
    plot!(df.time, after_min; label="After (min)", linestyle=:dash)
    plot!(df.time, after_max; label="After (max)", linestyle=:dash)

    savefig(p, outpath)
end

function fig8_like_duration(df::DataFrame, outpath::String)
    dc_base = duration_curve(df.base_Imax; high_is_worse=true)
    dc_pv   = duration_curve(df.pv_Imax;   high_is_worse=true)
    dc_stc  = duration_curve(df.stc_Imax;  high_is_worse=true)

    p = plot(dc_base.p, dc_base.x; label="Baseline", xlabel="Load percentile (%)", ylabel="Max phase current (A)",
             title="Load duration curve of maximum phase current (substation)")
    plot!(dc_pv.p, dc_pv.x; label="PV")
    plot!(dc_stc.p, dc_stc.x; label="PV + STATCOM")
    savefig(p, outpath)
end

function fig11_like_voltage_seq_distributions(df_seq::DataFrame, outdir::String)
    mkpath(outdir)

    # jittered scatter + median line per scenario, for V0 and V2/V1
    function scat_med(ycol::Symbol, title_str::String, outpath::String)
        scenarios = ["baseline", "pv", "statcom"]
        xs = Float64[]
        ys = Float64[]
        for (i, sc) in enumerate(scenarios)
            sub = df_seq[df_seq.scenario .== sc, :]
            y = sub[!, ycol]
            y = y[isfinite.(y)]
            append!(ys, y)
            append!(xs, fill(i, length(y)) .+ (rand(length(y)) .- 0.5) .* 0.15)
        end

        p = scatter(xs, ys; markersize=2, xlabel="", ylabel=String(ycol), title=title_str, legend=false)
        xticks!(p, (1:3, scenarios))

        for (i, sc) in enumerate(scenarios)
            sub = df_seq[df_seq.scenario .== sc, :]
            y = sub[!, ycol]
            y = y[isfinite.(y)]
            isempty(y) && continue
            med = median(y)
            plot!(p, [i-0.25, i+0.25], [med, med]; linewidth=3, label=false)
        end

        savefig(p, outpath)
    end

    scat_med(:V0, "Zero-sequence voltage magnitude | distribution across buses and timesteps",
             joinpath(outdir, "fig11_like_V0_scatter_median.png"))

    scat_med(:V2_over_V1, "Negative sequence ratio (V2/V1) | distribution across buses and timesteps",
             joinpath(outdir, "fig11_like_V2_over_V1_scatter_median.png"))
end

## --------------------------------------------------
# 6) Load baseline CSV and select PV-stress timesteps
# --------------------------------------------------
    
isfile(BASELINE_CSV) || error("Baseline CSV not found: $(BASELINE_CSV)")   # check file existence of baseline CSV which contains alpha_t and pf_status

df0 = CSV.read(BASELINE_CSV, DataFrame)           # load baseline CSV
df0 = df0[df0.pf_status .== "LOCALLY_SOLVED", :]  # keep only locally solved cases

# PV proxy computed here (answers the pv_pu_proxy question)
df0.pv_pu = pv_shape_simple(df0.time)               # compute simple PV shape proxy to get pv_pu_proxy
df0.score = df0.pv_pu ./ (df0.alpha_t .+ 1e-6)      # score = pv_pu_proxy / alpha_t which favors high PV and low load

df_sel = sort(df0, :score, rev=true)                # sort by score descending
df_sel = unique(df_sel, :timestep)                  # keep only unique timesteps
df_sel = df_sel[1:min(N_PV_STRESS, nrow(df_sel)), :]  # keep only top N_PV_STRESS timesteps
df_sel = sort(df_sel, :time)                        # sort back by time

if SKIP_NIGHT_PV                                    # filter out nighttime PV cases if specified
    df_sel = df_sel[df_sel.pv_pu .> 1e-8, :]
end

nrow(df_sel) == 0 && error("No PV-stress timesteps selected after filters")   # error if no timesteps left

df_sel.pv_kw_eff = PV_KW_PER_PHASE .* df_sel.pv_pu         # compute effective PV kW per phase

println("\nSelected PV-stress timesteps: ", nrow(df_sel))    

# Choose center timestep for 48h window
center_row = WINDOW_CENTER == "FIRST_STRESS" ? df_sel[1, :] : df_sel[argmax(df_sel.score), :]     # pick first stress or highest-score stress
center_k   = Int(center_row.timestep)            # get center timestep for 48h window

# Build 48h window from baseline CSV time axis
k_list = collect(df0.timestep)                 # list of all timesteps in baseline CSV
center_idx = findfirst(==(center_k), k_list)            # find index of center timestep
center_idx === nothing && error("Center timestep not found in baseline CSV")          # error if not found

half = Int(floor(WINDOW_STEPS ÷ 2))             # half window size
i1 = max(1, center_idx - half)                  # start index of window
i2 = min(nrow(df0), i1 + WINDOW_STEPS - 1)      # end index of window
df_win = df0[i1:i2, :]                          # extract 48h window
df_win.pv_pu = pv_shape_simple(df_win.time)     # compute simple PV shape proxy for window
df_win.pv_kw_eff = PV_KW_PER_PHASE .* df_win.pv_pu    # compute effective PV kW per phase for window

println("48h window rows: ", nrow(df_win), " | center timestep = ", center_k)    

# Save input plots for stress selection
plot_inputs(df_sel, joinpath(FIGDIR, "inputs_alpha_pvproxy_pveff_selected.png"))

# --------------------------------------------------
# 7) Parse feeders ONCE (baseline-clean and statcom-enabled)
# --------------------------------------------------

println("\nParsing BASELINE feeder: ", MASTER_BASELINE)
eng_base0 = PMD.parse_file(MASTER_BASELINE, transformations=[PMD.transform_loops!, reduce_lines!])   # parse baseline feeder file and reduce lines

@show haskey(eng_base0, "load")
@show length(get(eng_base0, "load", Dict()))





########################################################################
########################################################################
## ######################################################################
MASTER_BASELINE = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w_new.dss")
# MASTER_BASELINE = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", NET_4W, "master_baseline_4w.dss")
MASTER_BASELINE_3w = joinpath(ROOT, "data/raw/dsuite_networks_scaled_v1.1", "spd_s", "master_scaled_new.dss")

rank=1
r = df_sel[1, :]
k = Int(r.timestep)
a = Float64(r.alpha_t)      # load scaling factor
pv_pu = Float64(r.pv_pu)    # effective per-phase PV kW injection
pv_kw_eff = PV_KW_PER_PHASE * pv_pu    # effective PV kW per phase

load_multiplier = 50     

eng_base = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])   # parse baseline feeder file and reduce lines
scale_loads!(eng_base, a)
math_base = PMD.transform_data_model(eng_base; kron_reduce=true, phase_project=true)
for (i, bus) in math_base["bus"]
    bus["vmin"][1:3] .= 0.9
    bus["vmax"][1:3] .= 1.1
end
for (i, load) in math_base["load"]
    load["pd"] *= load_multiplier
end
pf_base = PMD.solve_mc_pf(math_base, PMD.IVRUPowerModel, ipopt)
m_base = pf_metrics(pf_base; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)


# PV only
eng_pv = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])   # parse baseline feeder file and reduce lines
scale_loads!(eng_pv, a)
pv_bus = first(keys(eng_pv["load"])) |> string    # pick first bus as PV bus
eng_pv["load"][pv_bus]["pd_nom"] .-= pv_kw_eff
math_pv = PMD.transform_data_model(eng_pv; kron_reduce=true, phase_project=true)
for (_, bus) in math_pv["bus"]
    bus["vmin"][1:3] .= 0.9
    bus["vmax"][1:3] .= 1.1
end
for (i, load) in math_pv["load"]
    load["pd"] *= load_multiplier
end
pf_pv = PMD.solve_mc_pf(math_pv, PMD.IVRUPowerModel, ipopt)
m_pv = pf_metrics(pf_pv; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)



include(joinpath(ROOT, "src/read_functions.jl"))

# PV + STATCOM
eng_stc = PMD.parse_file(MASTER_BASELINE_3w, transformations=[PMD.transform_loops!, reduce_lines!])   # parse baseline feeder file and reduce lines
scale_loads!(eng_stc, a)
pv_bus = first(keys(eng_stc["load"])) |> string    # pick first bus as PV bus
eng_stc["load"][pv_bus]["pd_nom"] .-= pv_kw_eff
math_stc = PMD.transform_data_model(eng_stc; kron_reduce=true, phase_project=true)
for (i, bus) in math_stc["bus"]
    bus["vmin"][1:3] .= 0.9
    bus["vmax"][1:3] .= 1.1
end
for (i, load) in math_stc["load"]
    load["pd"] *= load_multiplier
end
# Inverter gen ids expected in master_scaled.dss
stc_ids = [string(i) for i in 1:10]
for gen_id in stc_ids
    add_inverter_losses!(math_stc, gen_id; three_wire=true, c_rating_a=30*ones(3))
end
# Voltage unbalance constraint field used by the OPF builder
for (_, bus_any) in math_stc["bus"]
    bus = bus_any::Dict{String,Any}
    bus["vm_vuf_max"] = 0.02
end
PMD.add_start_vrvi!(math_stc)
model = PMD.instantiate_mc_model(math_stc, PMD.IVRUPowerModel, build_mc_opf_mx_3w_Rhea);
pf_stc = PMD.optimize_model!(model, optimizer=ipopt)
m_stc = pf_metrics(pf_stc; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)


ref_branch_id = 176
ref_branch_current_magnitude_base = abs.(pf_base["solution"]["branch"]["$ref_branch_id"]["cr_fr"] .+ im .* pf_base["solution"]["branch"]["$ref_branch_id"]["ci_fr"])
ref_branch_current_magnitude_pv = abs.(pf_pv["solution"]["branch"]["$ref_branch_id"]["cr_fr"] .+ im .* pf_pv["solution"]["branch"]["$ref_branch_id"]["ci_fr"])
ref_branch_current_magnitude_stc = abs.(pf_stc["solution"]["branch"]["$ref_branch_id"]["cr_fr"] .+ im .* pf_stc["solution"]["branch"]["$ref_branch_id"]["ci_fr"])

v_neg_base = []
v_zero_base = []
v_neg_pv = []
v_zero_pv = []
v_neg_stc = []
v_zero_stc = []

function vm_va_to_rectangular(vm::AbstractVector, va::AbstractVector)
    vr = vm .* cos.(va)
    vi = vm .* sin.(va)
    return (vr=vr, vi=vi)
end

for (i, bus) in pf_base["solution"]["bus"]
    # base results
    if haskey(bus, "vm")
        vr, vi = vm_va_to_rectangular(bus["vm"], bus["va"])
    else
        vr, vi = bus["vr"], bus["vi"]
    end
    _, _, v_seq_m = get_sequence_components(vr[1:3] + im .* vi[1:3])
    push!(v_neg_base, v_seq_m[3])
    push!(v_zero_base, v_seq_m[1])

    # pv results
    bus_res_pv = pf_pv["solution"]["bus"][i]
    if haskey(bus_res_pv, "vm")
        vr, vi = vm_va_to_rectangular(bus_res_pv["vm"], bus_res_pv["va"])
    else
        vr, vi = bus_res_pv["vr"], bus_res_pv["vi"]
    end
    _, _, v_seq_m = get_sequence_components(vr[1:3] + im .* vi[1:3])
    push!(v_neg_pv, v_seq_m[3])
    push!(v_zero_pv, v_seq_m[1])

    # statcom results
    bus_res_stc = pf_stc["solution"]["bus"][i]
    if haskey(bus_res_stc, "vm")
        vr, vi = vm_va_to_rectangular(bus_res_stc["vm"], bus_res_stc["va"])
    else
        vr, vi = bus_res_stc["vr"], bus_res_stc["vi"]
    end
    _, _, v_seq_m = get_sequence_components(vr[1:3] + im .* vi[1:3])
    push!(v_neg_stc, v_seq_m[3])
    push!(v_zero_stc, v_seq_m[1])
end

v_zero_plot = plot(v_zero_base; xlabel="bus index", label="base", title="Zero-sequence voltage magnitude", seriestype=:scatter)
plot!(v_zero_pv; label="pv", seriestype=:scatter)
plot!(v_zero_stc; label="stc", seriestype=:scatter)

v_neg_plot = plot(v_neg_base; xlabel="bus index", label="base", title="Negative-sequence voltage magnitude", seriestype=:scatter)
plot!(v_neg_pv; label="pv", seriestype=:scatter)
plot!(v_neg_stc; label="stc", seriestype=:scatter)






## ######################################################################
########################################################################
########################################################################

function sum_eng_PQ(eng)
    P = 0.0; Q = 0.0
    for ld in values(get(eng, "load", Dict()))
        p = get(ld, "pd_nom", get(ld, "pd", 0.0))
        q = get(ld, "qd_nom", get(ld, "qd", 0.0))
        P += p isa Number ? Float64(p) : sum(Float64.(p))
        Q += q isa Number ? Float64(q) : sum(Float64.(q))
    end
    return (P=P, Q=Q)
end

@show sum_eng_PQ(eng_base0)

@show sum_loads_eng(eng_base0)  # show total loads before scaling

@show haskey(eng_base0, "voltage_source")
@show keys(get(eng_base0, "voltage_source", Dict()))

# eng_base0["voltage_source"]["source"]["rs"] .= 0
# eng_base0["voltage_source"]["source"]["xs"] .= 0

eng_base_template = deepcopy(eng_base0)    # template to use for each timestep
scale_loads!(eng_base_template, LOAD_ALPHA_BASE)     # scale loads by base alpha
add_load_q_from_pf!(eng_base_template; pf=0.95)      # add some Q load based on power factor

@show sum_loads_eng(eng_base_template) # should now show P ~ 74*1.5 and Q > 0

# kw_to_w!(eng_base_template)   # convert kW loads to W for better numerical conditioning

# math_pf0 = PMD.transform_data_model(eng_base_template; kron_reduce=true, phase_project=true)
# # PMD.add_start_vrvi!(math_pf0)
# pf0 = PMD.solve_mc_pf(math_pf0, PMD.IVRUPowerModel, ipopt)
# @show pf0["termination_status"]
# @show pf_metrics(pf0)

math_base0 = PMD.transform_data_model(eng_base_template; kron_reduce=false, phase_project=false)    # transform to math model without reduction
length(math_base0["load"]) == 0 && error("No loads found in baseline feeder")                       # check loads exist

@show haskey(math_base0, "load")
@show length(get(math_base0, "load", Dict()))

ld = first(values(math_base0["load"]))
@show keys(ld)
@show get(ld, "pd", nothing)
@show get(ld, "qd", nothing)

function sum_math_PQ(math)
    P = 0.0; Q = 0.0
    for ld in values(get(math, "load", Dict()))
        p = get(ld, "pd", 0.0)
        q = get(ld, "qd", 0.0)
        P += p isa Number ? Float64(p) : sum(Float64.(p))
        Q += q isa Number ? Float64(q) : sum(Float64.(q))
    end
    return (P=P, Q=Q)
end

@show sum_math_PQ(math_base0)
@show math_base0["settings"]["sbase_default"]

for (_, bus_any) in math_base0["bus"]     # set voltage limits on all buses
    bus = bus_any::Dict{String,Any}
    bus["vmin"] .= 0.94
    bus["vmax"] .= 1.10
end

S = sbase_math(math_base0)                     # extract sbase from math model
pq = sum_loads_math_pu(math_base0)             # sum loads in per-unit from math model
@show S pq                                     # show sbase and per-unit loads 
@show (P_kW_est = pq.P_pu * S, Q_kvar_est = pq.Q_pu * S)  # estimated loads in kW and kVAR from math model

@show haskey(math_base0, "voltage_source")
@show keys(get(math_base0, "voltage_source", Dict()))

pf_builtin = PMD.solve_mc_pf(math_base0, PMD.IVRUPowerModel, ipopt)
@show pf_builtin["termination_status"]
@show pf_metrics(pf_builtin)

# unique(length.(get.(values(math_base0["bus"]), "vmin", [])))
# unique(length.(get.(values(math_base0["bus"]), "vmax", [])))

lines_df   = make_lines_df_from_eng(eng_base0)
source_bus = pick_source_bus_name(eng_base0)
dist       = compute_bus_distances(lines_df; source_bus=source_bus)

pv_bus = (PV_BUS == "AUTO") ? pick_pv_bus_auto(eng_base0, dist) : PV_BUS   # pick PV bus automatically or use specified one
println("Distance reference bus: ", source_bus)   
println("Chosen PV bus: ", pv_bus)

println("\nParsing STATCOM feeder: ", MASTER_STATCOM)
eng_stc0 = PMD.parse_file(MASTER_STATCOM, transformations=[PMD.transform_loops!, reduce_lines!])     # parse statcom feeder file and reduce lines
scale_loads!(eng_stc0, LOAD_ALPHA_BASE)     # scale loads by base alpha to match baseline and avoid unbalanced comparisons
eng_stc0["load"]    

# --------------------------------------------------
# 8) Main loop on PV-stress timesteps (tables + voltage sequences)
# --------------------------------------------------

rows = NamedTuple[]
seq_rows = NamedTuple[]

# for (rank, r) in enumerate(eachrow(df_sel))
    rank=1
    r = df_sel[1, :]
    k = Int(r.timestep)
    a = Float64(r.alpha_t)      # load scaling factor
    pv_pu = Float64(r.pv_pu)    # effective per-phase PV kW injection
    pv_kw_eff = PV_KW_PER_PHASE * pv_pu    # effective PV kW per phase

    # Baseline
    eng_base = deepcopy(eng_base_template)
    scale_loads!(eng_base, a)                     # scale loads by alpha_t

    total_kw = sum(sum(ld["pd_nom"]) for ld in values(eng_base["load"]))
    total_kw

    ld = first(values(eng_base["load"]))
    @show keys(ld)
    @show get(ld, "pd_nom", nothing)
    @show get(ld, "pd", nothing)
    @show get(ld, "qd_nom", nothing)
    @show get(ld, "qd", nothing)

    # kw_to_w!(eng_base)  # convert kW loads to W for better numerical conditioning

    math_base = PMD.transform_data_model(eng_base; kron_reduce=false, phase_project=false)
    
    for (_, bus) in math_base["bus"]
        bus["vmin"] .= 0.9
        bus["vmax"] .= 1.1
    end

    total_pu = total_kw * 1000 / math_base["settings"]["sbase_default"]
    total_pu
    
    pq = sum_loads_math_pu(math_base)
    pq.P_pu 

    @show pq
    @show (P_seen_kW = pq.P_pu * math_base["settings"]["sbase_default"] / 1000)

    for (i, load) in math_base["load"]
        load["pd"] *= 10
    end

    PMD.add_start_vrvi!(math_base)      # telling the solver where “normal” is before asking it to find the truth
    # pf_base = PMD.solve_mc_opf(math_base, PMD.IVRENPowerModel, ipopt)
    model = PMD.instantiate_mc_model(math_base, PMD.IVRENPowerModel, build_mc_pf_mx_Rhea)
    pf_base = PMD.optimize_model!(model, optimizer=ipopt)
    # pf_base, math_base = solve_pf_case(eng_base, ipopt)
    pf_base = PMD.solve_mc_pf(math_base, PMD.IVRENPowerModel, ipopt)
    pf_base = PMD.compute_mc_pf(math_base, explicit_neutral=true, kron_reduce=false, phase_project=false)   # ensure all quantities computed
    pf
    m_base = pf_metrics(pf_base; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    pf_base["solution"]["bus"] |> length
    pf_base["objective"]
    pf_base["solution"]["bus"]["1"]  # or any key that exists

    math_base["settings"]["sbase_default"]

    @show pf_base["termination_status"]
    @show sum_loads_eng(eng_base)
    @show feeder_head_current_metrics(pf_base, math_base, source_bus; vbase_ln=VBASE_LN)


    # PV only
    eng_pv = deepcopy(eng_base_template)
    scale_loads!(eng_pv, a)
    # kw_to_w!(eng_pv) # convert kW loads to W for better numerical conditioning
    pv_action, pv_load_id = add_pv_negative_load!(eng_pv, pv_bus, pv_kw_eff, PV_PHASES)
    math_pv = PMD.transform_data_model(eng_pv; kron_reduce=false, phase_project=false)
    for (_, bus) in math_pv["bus"]
        bus["vmin"] .= 0.94
        bus["vmax"] .= 1.1
    end
    PMD.add_start_vrvi!(math_pv)
    # pf_pv = PMD.solve_mc_opf(math_pv, PMD.IVRENPowerModel, ipopt)
    model = PMD.instantiate_mc_model(math_pv, PMD.IVRENPowerModel, build_mc_pf_mx_Rhea)
    pf_pv = PMD.optimize_model!(model, optimizer=ipopt)
    # pf_pv, math_pv = solve_pf_case(eng_pv, ipopt)
    m_pv = pf_metrics(pf_pv; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # PV + STATCOM
    eng_stc = deepcopy(eng_stc0)
    scale_loads!(eng_stc, a)
    # pv_action2, pv_load_id2 = add_pv_negative_load!(eng_stc, pv_bus, pv_kw_eff, PV_PHASES)
    math_stc = PMD.transform_data_model(eng_stc; kron_reduce=false, phase_project=false)
    for (_, bus) in math_stc["bus"]
        bus["vmin"] .= 0.9
        bus["vmax"] .= 1.1
    end
    # Inverter gen ids expected in master_scaled.dss
    stc_ids = [string(i) for i in 1:10]
    for gen_id in stc_ids
        add_inverter_losses!(math_stc, gen_id; three_wire=true, c_rating_a=30*ones(3))
    end
    # Voltage unbalance constraint field used by the OPF builder
    for (_, bus_any) in math_stc["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vm_vuf_max"] = 0.02
    end
    PMD.add_start_vrvi!(math_stc)
    model = PMD.instantiate_mc_model(math_stc, PMD.IVRENPowerModel, build_mc_opf_mx_Rhea)
    pf_stc = PMD.optimize_model!(model, optimizer=ipopt)
    
    # pf_stc, math_stc = solve_statcom_case(eng_stc, ipopt)
    m_stc = pf_metrics(pf_stc; vmin_pu=VMIN_PU, vmax_pu=VMAX_PU)

    # Debug generator counts: baseline and PV should not contain 1..10
    if rank == 1
        ngen_base = haskey(math_base, "gen") ? length(math_base["gen"]) : 0
        ngen_pv   = haskey(math_pv,   "gen") ? length(math_pv["gen"])   : 0
        ngen_stc  = haskey(math_stc,  "gen") ? length(math_stc["gen"])  : 0
        println("\nDEBUG gen counts: base=$(ngen_base) pv=$(ngen_pv) stc=$(ngen_stc)")
        if haskey(math_stc, "gen")
            println("DEBUG statcom gen ids sample: ", collect(keys(math_stc["gen"]))[1:min(12, ngen_stc)])
        end
    end

    # Feeder head currents (Fix A: robust head branch selection)
    Ibase = feeder_head_current_metrics(pf_base, math_base, source_bus; vbase_ln=VBASE_LN)
    Ipv   = feeder_head_current_metrics(pf_pv,   math_pv,   source_bus; vbase_ln=VBASE_LN)
    Istc  = feeder_head_current_metrics(pf_stc,  math_stc,  source_bus; vbase_ln=VBASE_LN)

    # Voltage sequences for all buses (Rahmat requirement)
    append!(seq_rows, voltage_sequence_rows(pf_base, r.time, "baseline"))
    append!(seq_rows, voltage_sequence_rows(pf_pv,   r.time, "pv"))
    append!(seq_rows, voltage_sequence_rows(pf_stc,  r.time, "statcom"))

    push!(rows, (
        net = NET_4W,
        year = YEAR,
        rank = rank,
        timestep = k,
        time = r.time,
        alpha_t = a,
        pv_pu_proxy = pv_pu,
        score = Float64(r.score),
        pv_bus = pv_bus,
        pv_kw_per_phase = PV_KW_PER_PHASE,
        pv_kw_eff = pv_kw_eff,
        pv_phases = join(PV_PHASES, ","),
        pv_action = pv_action,
        pv_load_id = pv_load_id,

        base_status = m_base.status,
        base_vmin_pu = m_base.vmin,
        base_vmax_pu = m_base.vmax,

        pv_status = m_pv.status,
        pv_vmin_pu = m_pv.vmin,
        pv_vmax_pu = m_pv.vmax,

        stc_status = m_stc.status,
        stc_vmin_pu = m_stc.vmin,
        stc_vmax_pu = m_stc.vmax,

        delta_pv_vmin_pu = m_pv.vmin - m_base.vmin,
        delta_pv_vmax_pu = m_pv.vmax - m_base.vmax,
        delta_stc_vmin_pu = m_stc.vmin - m_base.vmin,
        delta_stc_vmax_pu = m_stc.vmax - m_base.vmax,

        base_Ia = Ibase.Ia, base_Ib = Ibase.Ib, base_Ic = Ibase.Ic, base_In = Ibase.In,
        base_Imean = Ibase.Imean, base_Imax = Ibase.Imax,
        base_I0 = Ibase.I0, base_I2 = Ibase.I2, base_I2r = Ibase.I2_over_I1,

        pv_Ia = Ipv.Ia, pv_Ib = Ipv.Ib, pv_Ic = Ipv.Ic, pv_In = Ipv.In,
        pv_Imean = Ipv.Imean, pv_Imax = Ipv.Imax,
        pv_I0 = Ipv.I0, pv_I2 = Ipv.I2, pv_I2r = Ipv.I2_over_I1,

        stc_Ia = Istc.Ia, stc_Ib = Istc.Ib, stc_Ic = Istc.Ic, stc_In = Istc.In,
        stc_Imean = Istc.Imean, stc_Imax = Istc.Imax,
        stc_I0 = Istc.I0, stc_I2 = Istc.I2, stc_I2r = Istc.I2_over_I1
    ))

    println("Saved rank=", rank,
        " k=", k,
        " alpha=", round(a, digits=3),
        " pv_kw_eff=", round(pv_kw_eff, digits=3),
        " base_vmin=", round(m_base.vmin, digits=6),
        " pv_vmin=", round(m_pv.vmin, digits=6),
        " stc_vmin=", round(m_stc.vmin, digits=6),
        " | base_Imax=", round(Ibase.Imax, digits=3),
        " pv_Imax=", round(Ipv.Imax, digits=3),
        " stc_Imax=", round(Istc.Imax, digits=3)
    )
# end

## ##########
using Plots
vm_base_a = [abs.(bus["vr"][1] + im*bus["vi"][1]) for (i, bus) in pf_base["solution"]["bus"]]
vm_base_b = [abs.(bus["vr"][2] + im*bus["vi"][2]) for (i, bus) in pf_base["solution"]["bus"]]
vm_base_c = [abs.(bus["vr"][3] + im*bus["vi"][3]) for (i, bus) in pf_base["solution"]["bus"]]

vm_pv_a = [abs.(bus["vr"][1] + im*bus["vi"][1]) for (i, bus) in pf_pv["solution"]["bus"]]
vm_pv_b = [abs.(bus["vr"][2] + im*bus["vi"][2]) for (i, bus) in pf_pv["solution"]["bus"]]
vm_pv_c = [abs.(bus["vr"][3] + im*bus["vi"][3]) for (i, bus) in pf_pv["solution"]["bus"]]

vm_stc_a = [abs.(bus["vr"][1] + im*bus["vi"][1]) for (i, bus) in pf_stc["solution"]["bus"]]
vm_stc_b = [abs.(bus["vr"][2] + im*bus["vi"][2]) for (i, bus) in pf_stc["solution"]["bus"]]
vm_stc_c = [abs.(bus["vr"][3] + im*bus["vi"][3]) for (i, bus) in pf_stc["solution"]["bus"]]

plot(vm_base_a, label="base", xlabel="Bus index", color=:blue, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")
plot!(vm_base_b, label="base", xlabel="Bus index", color=:blue, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")
plot!(vm_base_c, label="base", xlabel="Bus index", color=:blue, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")

plot!(vm_pv_a, label="pv", xlabel="Bus index", color=:green, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")
plot!(vm_pv_b, label="pv", xlabel="Bus index", color=:green, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")
plot!(vm_pv_c, label="pv", xlabel="Bus index", color=:green, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")

plot!(vm_stc_a, label="stc", xlabel="Bus index", color=:red, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")
plot!(vm_stc_b, label="stc", xlabel="Bus index", color=:red, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")
plot!(vm_stc_c, label="stc", xlabel="Bus index", color=:red, ylabel="Voltage magnitude (V)", title="Voltage magnitudes at timestep $(k)")



### 
# - increase the load pd, and also set qd to a fixed power factor (e.g., 0.9)
# - calculate all voltage sequences at all buses for the three cases (x_seq_m[1:3] -> zero sequence, positive sequence, negative sequence)
# - calculate all branch flow sequences at all branches for the three cases (x_seq_f[1:3] -> zero sequence, positive sequence, negative sequence)
# - compare them

c_base = [get_sequence_components(branch["cr_fr"][1:3]+im*branch["ci_fr"][1:3]) for (i, branch) in pf_stc["solution"]["branch"]]
c_pv = [get_sequence_components(branch["cr_fr"][1:3]+im*branch["ci_fr"][1:3]) for (i, branch) in pf_stc["solution"]["branch"]]
c_stc = [get_sequence_components(branch["cr_fr"][1:3]+im*branch["ci_fr"][1:3]) for (i, branch) in pf_stc["solution"]["branch"]]


##################



df_out = DataFrame(rows)
CSV.write(joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps.csv"), df_out)
println("\nSaved table: ", joinpath(TBLDIR, "timeseries_pv_pf_selected_timesteps.csv"))

curr_df = select(df_out,
    :time,
    :base_Ia, :base_Ib, :base_Ic, :base_In, :base_Imean, :base_Imax,
    :pv_Ia,   :pv_Ib,   :pv_Ic,   :pv_In,   :pv_Imean,   :pv_Imax,
    :stc_Ia,  :stc_Ib,  :stc_Ic,  :stc_In,  :stc_Imean,  :stc_Imax
)
CSV.write(joinpath(TBLDIR, "substation_current_comparison_selected.csv"), curr_df)
println("Saved table: ", joinpath(TBLDIR, "substation_current_comparison_selected.csv"))

df_seq = DataFrame(seq_rows)
CSV.write(joinpath(TBLDIR, "voltage_sequences_all_buses_all_timesteps.csv"), df_seq)
println("Saved table: ", joinpath(TBLDIR, "voltage_sequences_all_buses_all_timesteps.csv"))

# --------------------------------------------------
# 10) Paper-style plots on PV-stress set
# --------------------------------------------------

# Delta voltage scatter plots (selected timesteps)
p1 = scatter(df_out.time, df_out.delta_pv_vmax_pu; xlabel="Time", ylabel="PV - baseline vmax (pu)",
    title="PV impact on vmax | selected timesteps", legend=false)
savefig(p1, joinpath(FIGDIR, "delta_vmax_selected_pv.png"))

p2 = scatter(df_out.time, df_out.delta_stc_vmax_pu; xlabel="Time", ylabel="STATCOM - baseline vmax (pu)",
    title="STATCOM impact on vmax | selected timesteps", legend=false)
savefig(p2, joinpath(FIGDIR, "delta_vmax_selected_statcom.png"))

p3 = scatter(df_out.time, df_out.delta_pv_vmin_pu; xlabel="Time", ylabel="PV - baseline vmin (pu)",
    title="PV impact on vmin | selected timesteps", legend=false)
savefig(p3, joinpath(FIGDIR, "delta_vmin_selected_pv.png"))

p4 = scatter(df_out.time, df_out.delta_stc_vmin_pu; xlabel="Time", ylabel="STATCOM - baseline vmin (pu)",
    title="STATCOM impact on vmin | selected timesteps", legend=false)
savefig(p4, joinpath(FIGDIR, "delta_vmin_selected_statcom.png"))

# Fig8-like duration curve on selected timesteps (max phase current)
fig8_like_duration(df_out, joinpath(FIGDIR, "fig8_like_duration_i_max_selected.png"))

# Fig11-like distributions for voltage sequences
fig11_like_voltage_seq_distributions(df_seq, FIGDIR)

# --------------------------------------------------
# 11) Evaluate scenarios over the 48-hour window (Fig6/Fig7/Fig9)
# --------------------------------------------------

win_rows = NamedTuple[]

for r in eachrow(df_win)
    k = Int(r.timestep)
    a = Float64(r.alpha_t)
    pv_pu = Float64(r.pv_pu)
    pv_kw_eff = PV_KW_PER_PHASE * pv_pu

    # Baseline
    eng_base = deepcopy(eng_base0)
    scale_loads!(eng_base, a)
    pf_base, math_base = solve_pf_case(eng_base, ipopt)
    Ibase = feeder_head_current_metrics(pf_base, math_base, source_bus; vbase_ln=VBASE_LN)

    # PV only
    eng_pv = deepcopy(eng_base0)
    scale_loads!(eng_pv, a)
    add_pv_negative_load!(eng_pv, pv_bus, pv_kw_eff, PV_PHASES)
    pf_pv, math_pv = solve_pf_case(eng_pv, ipopt)
    Ipv = feeder_head_current_metrics(pf_pv, math_pv, source_bus; vbase_ln=VBASE_LN)

    # PV + STATCOM
    eng_stc = deepcopy(eng_stc0)
    scale_loads!(eng_stc, a)
    add_pv_negative_load!(eng_stc, pv_bus, pv_kw_eff, PV_PHASES)
    pf_stc, math_stc = solve_statcom_case(eng_stc, ipopt)
    Istc = feeder_head_current_metrics(pf_stc, math_stc, source_bus; vbase_ln=VBASE_LN)

    push!(win_rows, (
        time = r.time,
        timestep = k,
        alpha_t = a,
        pv_pu = pv_pu,
        pv_kw_eff = pv_kw_eff,

        base_Ia = Ibase.Ia, base_Ib = Ibase.Ib, base_Ic = Ibase.Ic, base_In = Ibase.In,
        base_Imean = Ibase.Imean, base_Imax = Ibase.Imax,

        pv_Ia = Ipv.Ia, pv_Ib = Ipv.Ib, pv_Ic = Ipv.Ic, pv_In = Ipv.In,
        pv_Imean = Ipv.Imean, pv_Imax = Ipv.Imax,

        stc_Ia = Istc.Ia, stc_Ib = Istc.Ib, stc_Ic = Istc.Ic, stc_In = Istc.In,
        stc_Imean = Istc.Imean, stc_Imax = Istc.Imax
    ))
end

df_win_out = DataFrame(win_rows)
CSV.write(joinpath(TBLDIR, "window_48h_substation_currents.csv"), df_win_out)
println("Saved table: ", joinpath(TBLDIR, "window_48h_substation_currents.csv"))

# Fig6-like for each scenario
fig6_like_currents(df_win_out, "base", joinpath(FIGDIR, "fig6_like_substation_currents_baseline_48h.png"))
fig6_like_currents(df_win_out, "pv",   joinpath(FIGDIR, "fig6_like_substation_currents_pv_48h.png"))
fig6_like_currents(df_win_out, "stc",  joinpath(FIGDIR, "fig6_like_substation_currents_statcom_48h.png"))

# Fig7-like: before vs after range band (baseline vs STATCOM)
fig7_like_range_band(df_win_out, "base", "stc", joinpath(FIGDIR, "fig7_like_rangeband_baseline_vs_statcom_48h.png"))

# Fig9-like: phase + neutral overlay baseline vs STATCOM
p9a = plot(df_win_out.time, df_win_out.base_Ia; label="Phase a (before)", xlabel="Time", ylabel="Current (A)",
           title="Phase and neutral currents | before vs after")
plot!(df_win_out.time, df_win_out.stc_Ia; label="Phase a (after)")
plot!(df_win_out.time, df_win_out.base_Ib; label="Phase b (before)")
plot!(df_win_out.time, df_win_out.stc_Ib; label="Phase b (after)")
plot!(df_win_out.time, df_win_out.base_Ic; label="Phase c (before)")
plot!(df_win_out.time, df_win_out.stc_Ic; label="Phase c (after)")
if any(isfinite.(df_win_out.base_In)) || any(isfinite.(df_win_out.stc_In))
    plot!(df_win_out.time, df_win_out.base_In; label="Neutral (before)")
    plot!(df_win_out.time, df_win_out.stc_In; label="Neutral (after)")
end
savefig(p9a, joinpath(FIGDIR, "fig9_like_phase_neutral_baseline_vs_statcom_48h.png"))

# Fig8-like duration curve on 48h window
fig8_like_duration(df_win_out, joinpath(FIGDIR, "fig8_like_duration_i_max_window_48h.png"))

println("\nSaved results to: ", OUTDIR)
println("Done.")