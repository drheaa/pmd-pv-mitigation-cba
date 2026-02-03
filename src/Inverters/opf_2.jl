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




"""
	function build_mc_opf_mx(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals including reconfigurable inverters
"""
function build_mc_opf_mx_3w_Rhea(pm::PMD.AbstractUnbalancedIVRModel)
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

        # if i in PMD.ids(pm, :ref_buses)
        #     PMD.constraint_mc_voltage_reference(pm, i)
        # end

        # PMD.constraint_mc_voltage_absolute(pm, i)
        # PMD.constraint_mc_voltage_pairwise(pm, i)

        constraint_mc_bus_voltage_balance_Rhea(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        if id ∈ gen_ids  # Generators connected with inverter
            constraint_mc_generator_power_stc(pm, id)
            # PMD.constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current(pm, id)
            constraint_mc_generator_current_limit(pm, id)
            
        else  # Other generators
            constraint_mc_generator_power(pm, id)
            # PMD.constraint_mc_generator_current(pm, id)
        end
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        # PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        # constraint_mc_current_from(pm, i)
        # constraint_mc_current_to(pm, i)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)
        PMD.constraint_mc_bus_voltage_drop(pm, i)
        constraint_mc_branch_current_limit(pm, i)
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
        # PMD.constraint_mc_current_balance(pm, i)
        constraint_mc_current_balance(pm, i)
    end

    # Objective
    # PMD.objective_mc_min_fuel_cost(pm)
    # objective_mc_min_IUF(pm)
    objective_mc_min_max_phase_current(pm)
    # objective_mc_min_ref_branch_loss(pm)
end


function constraint_mc_current_balance(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    bus = PMD.ref(pm, nw, :bus, i)
    bus_arcs = PMD.ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_arcs_sw = PMD.ref(pm, nw, :bus_arcs_conns_switch, i)
    bus_arcs_trans = PMD.ref(pm, nw, :bus_arcs_conns_transformer, i)
    bus_gens = PMD.ref(pm, nw, :bus_conns_gen, i)
    bus_storage = PMD.ref(pm, nw, :bus_conns_storage, i)
    bus_loads = PMD.ref(pm, nw, :bus_conns_load, i)
    bus_shunts = PMD.ref(pm, nw, :bus_conns_shunt, i)

    constraint_mc_current_balance(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_arcs_sw, bus_arcs_trans, bus_gens, bus_storage, bus_loads, bus_shunts)
    nothing
end

function constraint_mc_current_balance(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_sw::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_arcs_trans::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_storage::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    cr    = get(PMD.var(pm, nw),    :cr, Dict()); PMD._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(PMD.var(pm, nw),    :ci, Dict()); PMD._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(PMD.var(pm, nw),   :crd_bus, Dict()); PMD._check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(PMD.var(pm, nw),   :cid_bus, Dict()); PMD._check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(PMD.var(pm, nw),   :crg_bus, Dict()); PMD._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(PMD.var(pm, nw),   :cig_bus, Dict()); PMD._check_var_keys(cig, bus_gens, "imaginary current", "generator")
    crs   = get(PMD.var(pm, nw),   :crs, Dict()); PMD._check_var_keys(crs, bus_storage, "real currentr", "storage")
    cis   = get(PMD.var(pm, nw),   :cis, Dict()); PMD._check_var_keys(cis, bus_storage, "imaginary current", "storage")
    crsw  = get(PMD.var(pm, nw),  :crsw, Dict()); PMD._check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    cisw  = get(PMD.var(pm, nw),  :cisw, Dict()); PMD._check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    crt   = get(PMD.var(pm, nw),   :crt, Dict()); PMD._check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    cit   = get(PMD.var(pm, nw),   :cit, Dict()); PMD._check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")
    Gt, Bt = PMD._build_bus_shunt_matrices(pm, nw, terminals, bus_shunts)

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]

    for (idx, t) in ungrounded_terminals
        JuMP.@constraint(pm.model,
                                      sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vr[u] -Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
        JuMP.@constraint(pm.model,
                                      sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                    + sum(cisw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    + sum(cit[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                    - sum(cis[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vi[u] +Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
    end
end
