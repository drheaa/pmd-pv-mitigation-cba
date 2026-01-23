"""
cache_load_bases(eng)

Stores base load values so scaling never accumulates across timesteps.
Returns Dict(load_id => Dict(field => base_value)).
"""
function cache_load_bases(eng::Dict{String,Any})
    cache = Dict{Any, Dict{String,Any}}()
    haskey(eng, "load") || return cache

    for (lid, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        d = Dict{String,Any}()

        for k in ("pd_nom", "qd_nom", "pd", "qd")
            if haskey(ld, k)
                d[k] = deepcopy(ld[k])
            end
        end
        cache[lid] = d
    end
    return cache
end

"""
apply_load_scaling!(eng, cache, scale)

Writes load values as base .* scale for all cached fields.
This is fast and avoids cumulative multiplication bugs.
"""
function apply_load_scaling!(eng::Dict{String,Any}, cache::Dict, scale::Real)
    haskey(eng, "load") || return

    for (lid, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        base = get(cache, lid, nothing)
        base === nothing && continue

        for (k, v0) in base
            if v0 isa Number
                ld[k] = scale * v0
            else
                ld[k] = scale .* v0
            end
        end
    end
end

function print_load_summary(eng::Dict{String,Any}; title::String="Load summary")
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

    println("\n", title)
    println("  number of loads = ", nloads)
    println("  total P (kW)    = ", total_pd)
    println("  total Q (kvar)  = ", total_qd)
end

function sanity_check_no_mitigation_loads(eng::Dict{String,Any})
    if !haskey(eng, "load")
        println("\nMitigation sanity check: no loads table found")
        return
    end

    keys_load = collect(keys(eng["load"]))
    n_statcom = count(k -> occursin("statcom", lowercase(string(k))), keys_load)
    n_sop     = count(k -> occursin("sop", lowercase(string(k))), keys_load)
    n_str     = count(k -> occursin("str", lowercase(string(k))), keys_load)

    println("\nMitigation sanity check (expected 0 for baseline master)")
    println("  loads containing 'statcom' = ", n_statcom)
    println("  loads containing 'sop'     = ", n_sop)
    println("  loads containing 'str'     = ", n_str)
end
