module Metrics

export bus_eng_name,
       solved_bus_vm_volts_keyed_by_eng,
       voltage_stats_pu

function bus_eng_name(bus_id::String, bus_data::Dict{String,Any})::String
    if haskey(bus_data, "name")
        return lowercase(string(bus_data["name"]))
    end
    if haskey(bus_data, "source_id")
        sid = lowercase(string(bus_data["source_id"]))
        occursin("bus.", sid) && return split(sid,"bus.")[end]
    end
    return lowercase(bus_id)
end

function solved_bus_vm_volts_keyed_by_eng(pf, math; vbase_ln=230.0)
    sol_bus = pf["solution"]["bus"]
    out = Dict{String,Dict{String,Any}}()

    for (bus_id_any, bus_data_any) in math["bus"]
        bus_id = string(bus_id_any)
        haskey(sol_bus,bus_id) || continue
        sb = sol_bus[bus_id]
        vm_pu = haskey(sb,"vm") ? sb["vm"] : sqrt.(sb["vr"].^2 .+ sb["vi"].^2)
        eng = bus_eng_name(bus_id, bus_data_any)
        out[eng] = Dict("vma"=>vm_pu[1]*vbase_ln,
                        "vmb"=>vm_pu[2]*vbase_ln,
                        "vmc"=>vm_pu[3]*vbase_ln)
    end
    return out
end

function voltage_stats_pu(buses; vbase_ln=230.0)
    vmins = [min(b["vma"],b["vmb"],b["vmc"]) / vbase_ln for b in values(buses)]
    sort!(vmins)
    q(p) = vmins[clamp(Int(ceil(p*length(vmins))),1,length(vmins))]
    return (min=minimum(vmins), q05=q(0.05), median=q(0.5), q95=q(0.95))
end

end
