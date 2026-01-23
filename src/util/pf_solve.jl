"""
solve_pf_3w_projected(eng)

Projected representation for stability.
Useful for fast baselines and early development.
"""
function solve_pf_3w_projected(eng::Dict{String,Any})
    math = PMD.transform_data_model(
        eng;
        multinetwork=false,
        kron_reduce=true,
        phase_project=true
    )
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")
    pf = PMD.solve_mc_pf(math, PMD.IVRUPowerModel, ipopt)
    return pf, math
end

"""
solve_pf_4w_explicit(eng)

Explicit neutral representation.
Use for neutral current, zero sequence current, and 4-wire device behavior.
"""
function solve_pf_4w_explicit(eng::Dict{String,Any})
    math = PMD.transform_data_model(
        eng;
        multinetwork=false,
        kron_reduce=false,
        phase_project=false
    )
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")
    pf = PMD.solve_mc_pf(math, PMD.IVRENPowerModel, ipopt)
    return pf, math
end
