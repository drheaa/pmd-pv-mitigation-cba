"""
    constraint_mc_bus_voltage_balance(pm::AbstractUnbalancedACRModel, bus_id::Int; nw=nw_id_default)::Nothing

Template function for bus voltage balance constraints.
"""
function constraint_mc_bus_voltage_balance_Rhea(pm::PMD.AbstractUnbalancedACRModel, bus_id::Int; nw=PMD.nw_id_default)::Nothing
    # @assert(length(PMD.ref(pm, nw, :conductor_ids))==3)

    bus = PMD.ref(pm, nw, :bus, bus_id)
    # constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, 0)

    if haskey(bus, "vm_vuf_max")
        constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    end

    # if haskey(bus, "vm_seq_neg_max")
    #     constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, bus["vm_seq_neg_max"])
    # end

    # if haskey(bus, "vm_seq_pos_max")
    #     constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"])
    # end

    # if haskey(bus, "vm_seq_zero_max")
    #     constraint_mc_bus_voltage_magnitude_zero_sequence(pm, nw, bus_id, bus["vm_seq_zero_max"])
    # end

    # if haskey(bus, "vm_ll_min")|| haskey(bus, "vm_ll_max")
    #     vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : fill(0, 3)
    #     vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : fill(Inf, 3)
    #     PMD.constraint_mc_bus_voltage_magnitude_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    # end
    nothing
end

"""
	function build_mc_opf_mx(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals including reconfigurable inverters
"""
function build_mc_opf_mx_Rhea(pm::PMD.AbstractExplicitNeutralIVRModel)
    inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))


    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    # variable_mc_generator_current(pm)
    # variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)

        constraint_mc_bus_voltage_balance_Rhea(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        if id ∈ gen_ids  # Generators connected with inverter
            constraint_mc_generator_power(pm, id)
            # PMD.constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current_limit(pm, id)
            
        else  # Other generators
            PMD.constraint_mc_generator_power(pm, id)
            PMD.constraint_mc_generator_current(pm, id)
        end
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)
        PMD.constraint_mc_bus_voltage_drop(pm, i)
        PMD.constraint_mc_branch_current_limit(pm, i)
        # PMD.constraint_mc_thermal_limit_from(pm, i)
        # PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    # PMD.objective_mc_min_fuel_cost(pm)
    objective_mc_min_IUF(pm)
    # objective_mc_min_max_phase_current(pm)
    # objective_mc_min_ref_branch_loss(pm)
end



########## pf case ########## 
function build_mc_pf_mx_Rhea(pm::PMD.AbstractExplicitNeutralIVRModel)
    # inverter_branches = [branch["index"] for (i, branch) in pm.data["branch"] if occursin("inverter_branch", branch["name"])]
    # gen_ids, gen_bus_ids, gen_branch_ids = get_pv_bus_branch(pm.ref[:it][:pmd][:nw][0])  # maybe output branch_ids and arcs seperately?
    # branch_gens = Dict(branch_id[1] => gen_ids[i] for (i, branch_id) in enumerate(gen_branch_ids))


    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    # variable_mc_generator_current(pm)
    # variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)

        constraint_mc_bus_voltage_balance_Rhea(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    # for id in PMD.ids(pm, :gen)
    #     if id ∈ gen_ids  # Generators connected with inverter
    #         constraint_mc_generator_power(pm, id)
    #         # PMD.constraint_mc_generator_current(pm, id)
    #         constraint_mc_generator_current(pm, id)
    #         constraint_mc_generator_current_limit(pm, id)
            
    #     else  # Other generators
    #         PMD.constraint_mc_generator_power(pm, id)
    #         PMD.constraint_mc_generator_current(pm, id)
    #     end
    # end

    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
        PMD.constraint_mc_generator_current(pm, id)
    end


    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)
        PMD.constraint_mc_bus_voltage_drop(pm, i)
        PMD.constraint_mc_branch_current_limit(pm, i)
        # PMD.constraint_mc_thermal_limit_from(pm, i)
        # PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    # PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_min_IUF(pm)
    # objective_mc_min_max_phase_current(pm)
    # objective_mc_min_ref_branch_loss(pm)
    JuMP.@objective(pm.model, Min, 0.0)

end



