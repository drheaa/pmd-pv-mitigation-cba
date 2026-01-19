module Network

using DataFrames

export scale_loads!,
       pick_source_bus_eng,
       make_edges_from_lines,
       make_edges_from_transformers,
       compute_bus_distances

function scale_loads!(eng::Dict{String,Any}, alpha::Real)
    haskey(eng, "load") || return eng
    for (_, ld_any) in eng["load"]
        ld = ld_any::Dict{String,Any}
        for key in ("pd", "qd")
            if haskey(ld, key)
                v = ld[key]
                if v isa Number
                    ld[key] = alpha * v
                elseif v isa AbstractVector
                    ld[key] = alpha .* v
                end
            end
        end
    end
    return eng
end

function pick_source_bus_eng(eng)::String
    if haskey(eng, "bus")
        haskey(eng["bus"], "sourcebusz") && return "sourcebusz"
        haskey(eng["bus"], "sourcebusZ") && return "sourcebusZ"
        haskey(eng["bus"], "sourcebus")  && return "sourcebus"
    end
    return first(keys(eng["bus"]))
end

function make_edges_from_lines(eng)::DataFrame
    rows = NamedTuple[]
    haskey(eng, "line") || return DataFrame(rows)

    for (id, ln) in eng["line"]
        b1 = lowercase(string(get(ln, "f_bus", "")))
        b2 = lowercase(string(get(ln, "t_bus", "")))
        len_km = get(ln, "length", 0.0) / 1000.0
        (isempty(b1) || isempty(b2)) && continue
        push!(rows, (edge_id=string(id), Bus1=b1, Bus2=b2, length_km=len_km, kind="line"))
    end
    return DataFrame(rows)
end

function make_edges_from_transformers(eng)::DataFrame
    rows = NamedTuple[]
    haskey(eng, "transformer") || return DataFrame(rows)

    for (id, tx) in eng["transformer"]
        buses = haskey(tx,"bus") ? tx["bus"] : haskey(tx,"buses") ? tx["buses"] : nothing
        buses === nothing || length(buses) < 2 && continue
        push!(rows, (
            edge_id=string(id),
            Bus1=lowercase(string(buses[1])),
            Bus2=lowercase(string(buses[2])),
            length_km=0.0,
            kind="transformer"
        ))
    end
    return DataFrame(rows)
end

function compute_bus_distances(edges_df::DataFrame; source_bus::String)
    adj = Dict{String,Vector{Tuple{String,Float64}}}()
    for r in eachrow(edges_df)
        push!(get!(adj, r.Bus1, Tuple{String,Float64}[]), (r.Bus2, r.length_km))
        push!(get!(adj, r.Bus2, Tuple{String,Float64}[]), (r.Bus1, r.length_km))
    end

    dist = Dict(source_bus => 0.0)
    queue = [source_bus]
    while !isempty(queue)
        u = popfirst!(queue)
        for (v,w) in get(adj,u,[])
            haskey(dist,v) && continue
            dist[v] = dist[u] + w
            push!(queue,v)
        end
    end
    return dist
end

end
