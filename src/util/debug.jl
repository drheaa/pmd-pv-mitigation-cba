"""
save_debug_snapshot(outdir, k, pf, math)

Saves:
  - debug_pf_k=....jls
  - debug_math_k=....jls
  - debug_solution_keys.txt

This is a reproducibility hook that helps trace missing fields later.
"""
function save_debug_snapshot(outdir::String, k::Int, pf::Dict{String,Any}, math::Dict{String,Any})
    serialize(joinpath(outdir, "debug_pf_k=$(k).jls"), pf)
    serialize(joinpath(outdir, "debug_math_k=$(k).jls"), math)

    open(joinpath(outdir, "debug_solution_keys.txt"), "w") do io
        println(io, "solution keys: ", keys(pf["solution"]))
        if haskey(pf["solution"], "branch") && !isempty(pf["solution"]["branch"])
            bid = first(keys(pf["solution"]["branch"]))
            println(io, "branch id: ", bid)
            println(io, "branch fields: ", keys(pf["solution"]["branch"][bid]))
        end
    end
end
