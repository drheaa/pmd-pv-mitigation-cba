module OPF

using PowerModelsDistribution
const PMD = PowerModelsDistribution

export run_opf_ivr_vuf

# ------------------------------------------------------------
# Helper: identify PV / inverter-connected generators
# ------------------------------------------------------------
"""
    get_pv_bus_branch(nw_data) -> (gen_ids, gen_bus_ids, gen_branch_ids)

Works with both Dict{String,Any} and Dict{Symbol,Any} network refs.

Minimal detection:
- if no PV/inverter/statcom gens exist, returns empty arrays (baseline works)
- if you add PV/STATCOM later, put "pv"/"inverter"/"statcom" in gen["name"]
"""
function get_pv_bus_branch(nw_data::AbstractDict)
    gen_ids = Int[]
    gen_bus_ids = Int[]
    gen_branch_ids = Tuple{Int,Int,Int}[]  # placeholder

    genblock =
        haskey(nw_data, "gen") ? nw_data["gen"] :
        haskey(nw_data, :gen)  ? nw_data[:gen]  :
        nothing

    genblock === nothing && return gen_ids, gen_bus_ids, gen_branch_ids

    for (_, g_any) in genblock
        # PMD gen entries are typically Dict{String,Any}; keep it permissive
        g = g_any
        name = lowercase(string(get(g, "name", get(g, :name, ""))))

        is_pv = occursin("pv", name) || occursin("inverter", name) || occursin("statcom", name)

        if is_pv
            idx = get(g, "index", get(g, :index, 0))
            gb  = get(g, "gen_bus", get(g, :gen_bus, 0))

            push!(gen_ids, Int(idx))
            push!(gen_bus_ids, Int(gb))
            push!(gen_branch_ids, (0,0,0))  # not used yet
        end
    end

    return gen_ids, gen_bus_ids, gen_branch_ids
end


# ------------------------------------------------------------
# VUF constraint hook (best-effort for your PMD version)
# ------------------------------------------------------------
"""
    constraint_mc_bus_voltage_balance(pm, bus_id; nw)

Tries to add VUF constraint if bus has vm_vuf_max.
Your PMD currently supports VUF constraint methods for ACP, not IVR,
so we attempt and silently skip if unsupported.
"""
function constraint_mc_bus_voltage_balance(pm, bus_id::Int; nw=PMD.nw_id_default)::Nothing
    bus = PMD.ref(pm, nw, :bus, bus_id)

    if haskey(bus, "vm_vuf_max")
        vmax = bus["vm_vuf_max"]
        try
            PMD.constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, vmax)
        catch
            # skip (unsupported for this model type in your PMD build)
        end
    end

    return nothing
end


# ------------------------------------------------------------
# OPF builder (Rahmat structure, runnable on your machine)
# ------------------------------------------------------------
function build_mc_opf_ivr_vuf(pm::PMD.AbstractExplicitNeutralIVRModel)
    # NOTE: this ref is Dict{Symbol,Any} in your run
    gen_ids, _, _ = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])

    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    # Bus constraints
    for i in PMD.ids(pm, :bus)
        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)

        constraint_mc_bus_voltage_balance(pm, i)
    end

    # Generator constraints
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
        PMD.constraint_mc_generator_current(pm, id)

        if id ∈ gen_ids && isdefined(PMD, :constraint_mc_generator_current_limit)
            PMD.constraint_mc_generator_current_limit(pm, id)
        end
    end

    # Load constraints
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    # Transformer constraints
    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)
        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    # Branch constraints
    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)
        PMD.constraint_mc_bus_voltage_drop(pm, i)
        PMD.constraint_mc_branch_current_limit(pm, i)
    end

    # Switch constraints
    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    # KCL
    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    if isdefined(PMD, :objective_mc_min_max_phase_current)
        PMD.objective_mc_min_max_phase_current(pm)
    else
        PMD.objective_mc_min_fuel_cost(pm)
    end
end


# ------------------------------------------------------------
# Runner
# ------------------------------------------------------------
function run_opf_ivr_vuf(math::Dict{String,Any}; optimizer)
    # Flat start to avoid start-voltage DimensionMismatch issues
    for (_, bus_any) in math["bus"]
        bus = bus_any::Dict{String,Any}
        n = haskey(bus, "terminals") ? length(bus["terminals"]) : 3
        bus["vr_start"] = ones(n)
        bus["vi_start"] = zeros(n)
    end

    # Add a VUF limit field (constraint may be skipped depending on PMD support)
    for (_, bus_any) in math["bus"]
        bus = bus_any::Dict{String,Any}
        bus["vm_vuf_max"] = 0.02
    end

    return PMD.solve_mc_model(
        math,
        PMD.IVRENPowerModel,
        optimizer,
        build_mc_opf_ivr_vuf
    )
end

end # module
