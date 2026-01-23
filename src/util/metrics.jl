function pf_voltage_metrics(pf::Dict{String,Any}; vmin_pu=0.90, vmax_pu=1.10)
    sol_bus = pf["solution"]["bus"]
    vmins = Float64[]
    vmaxs = Float64[]

    for sb in values(sol_bus)
        vm =
            haskey(sb, "vm") ? sb["vm"][1:min(3, length(sb["vm"]))] :
            (haskey(sb, "vr") && haskey(sb, "vi")) ? sqrt.(sb["vr"][1:min(3, length(sb["vr"]))].^2 .+ sb["vi"][1:min(3, length(sb["vi"]))].^2) :
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

function pf_system_vuf_max(pf::Dict{String,Any})
    sol_bus = pf["solution"]["bus"]
    vufs = Float64[]

    for sb in values(sol_bus)
        haskey(sb, "vr") && haskey(sb, "vi") || continue
        vr = sb["vr"]; vi = sb["vi"]
        (length(vr) >= 3 && length(vi) >= 3) || continue

        Va = complex(vr[1], vi[1])
        Vb = complex(vr[2], vi[2])
        Vc = complex(vr[3], vi[3])

        mags = seq_magnitudes(Va, Vb, Vc)
        isfinite(mags.A2_over_A1) && push!(vufs, mags.A2_over_A1)
    end

    isempty(vufs) && return NaN
    return maximum(vufs)
end

function pf_system_vseq_max(pf::Dict{String,Any})
    sol_bus = pf["solution"]["bus"]
    v0s = Float64[]
    v2s = Float64[]
    rats = Float64[]

    for sb in values(sol_bus)
        haskey(sb, "vr") && haskey(sb, "vi") || continue
        vr = sb["vr"]; vi = sb["vi"]
        (length(vr) >= 3 && length(vi) >= 3) || continue

        Va = complex(vr[1], vi[1])
        Vb = complex(vr[2], vi[2])
        Vc = complex(vr[3], vi[3])

        mags = seq_magnitudes(Va, Vb, Vc)
        isfinite(mags.A0) && push!(v0s, mags.A0)
        isfinite(mags.A2) && push!(v2s, mags.A2)
        isfinite(mags.A2_over_A1) && push!(rats, mags.A2_over_A1)
    end

    return (
        V0_max = isempty(v0s) ? NaN : maximum(v0s),
        V2_max = isempty(v2s) ? NaN : maximum(v2s),
        V2_over_V1_max = isempty(rats) ? NaN : maximum(rats)
    )
end

function pf_transformer_or_slack_p_kw(pf::Dict{String,Any})
    sol = pf["solution"]

    if haskey(sol, "transformer") && !isempty(sol["transformer"])
        vals = Float64[]
        for tr in values(sol["transformer"])
            key = haskey(tr, "pf") ? "pf" : (haskey(tr, "pt") ? "pt" : nothing)
            key === nothing && continue
            v = tr[key]
            if v isa AbstractVector
                push!(vals, sum(Float64.(v[1:min(3, length(v))])))
            elseif v isa Number
                push!(vals, Float64(v))
            end
        end
        !isempty(vals) && return maximum(vals)
    end

    if haskey(sol, "voltage_source") && !isempty(sol["voltage_source"])
        vals = Float64[]
        for vs in values(sol["voltage_source"])
            haskey(vs, "pg") || continue
            v = vs["pg"]
            if v isa AbstractVector
                push!(vals, sum(Float64.(v[1:min(3, length(v))])))
            elseif v isa Number
                push!(vals, Float64(v))
            end
        end
        !isempty(vals) && return maximum(vals)
    end

    if haskey(sol, "gen") && !isempty(sol["gen"])
        vals = Float64[]
        for g in values(sol["gen"])
            haskey(g, "pg") || continue
            v = g["pg"]
            if v isa AbstractVector
                push!(vals, sum(Float64.(v[1:min(3, length(v))])))
            elseif v isa Number
                push!(vals, Float64(v))
            end
        end
        !isempty(vals) && return maximum(vals)
    end

    return NaN
end

function _phasor_from_solution(obj::Dict{String,Any}, rkey::String, ikey::String, idx::Int)
    haskey(obj, rkey) && haskey(obj, ikey) || return nothing
    r = obj[rkey]; i = obj[ikey]
    (length(r) >= idx && length(i) >= idx) || return nothing
    return complex(Float64(r[idx]), Float64(i[idx]))
end

"""
feeder_head_currents(pf)

Extracts phase current magnitudes (A) from the PF solution using branch current fields.
Supported key pairs (most common in PMD):
  - cr_fr, ci_fr  (from side)
  - cr_to, ci_to  (to side)

Selection rule:
  - iterates over all branches with available current phasors
  - chooses the branch with the largest mean phase current as a feeder-head proxy

Returns:
  (Ia, Ib, Ic, Imax, Imean, I0, I2, I2_over_I1)

Notes:
  - This is a proxy for "substation current" until a specific feeder-head branch
    is selected using network topology and the source bus.
"""
function feeder_head_currents(pf::Dict{String,Any})
    sol = pf["solution"]

    if !haskey(sol, "branch") || isempty(sol["branch"])
        return (Ia=NaN, Ib=NaN, Ic=NaN, Imax=NaN, Imean=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)
    end

    # Common PMD branch current key pairs
    pairs = [
        ("cr_fr", "ci_fr"),
        ("cr_to", "ci_to"),
        # optional fallbacks for older variants
        ("cfr", "cfi"),
        ("cr",  "ci"),
        ("ctr", "cti"),
    ]

    best = nothing
    best_score = -Inf

    for br_any in values(sol["branch"])
        br = br_any::Dict{String,Any}

        for (rk, ik) in pairs
            haskey(br, rk) && haskey(br, ik) || continue

            cr = br[rk]
            ci = br[ik]
            (cr isa AbstractVector && ci isa AbstractVector) || continue
            (length(cr) >= 3 && length(ci) >= 3) || continue

            Ia = complex(Float64(cr[1]), Float64(ci[1]))
            Ib = complex(Float64(cr[2]), Float64(ci[2]))
            Ic = complex(Float64(cr[3]), Float64(ci[3]))

            ima = abs(Ia); imb = abs(Ib); imc = abs(Ic)
            imean = (ima + imb + imc) / 3

            if imean > best_score
                mags = seq_magnitudes(Ia, Ib, Ic)
                best_score = imean
                best = (
                    Ia = ima,
                    Ib = imb,
                    Ic = imc,
                    Imax = max(ima, imb, imc),
                    Imean = imean,
                    I0 = mags.A0,
                    I2 = mags.A2,
                    I2_over_I1 = mags.A2_over_A1
                )
            end
        end
    end

    best === nothing && return (Ia=NaN, Ib=NaN, Ic=NaN, Imax=NaN, Imean=NaN, I0=NaN, I2=NaN, I2_over_I1=NaN)
    return best
end

function bus_voltage_table(pf::Dict{String,Any})
    sol_bus = pf["solution"]["bus"]
    rows = NamedTuple[]

    for (bid, sb_any) in sol_bus
        sb = sb_any::Dict{String,Any}

        if haskey(sb, "vm")
            vm = sb["vm"]
            vm3 = vm[1:min(3, length(vm))]
            vmin = minimum(vm3)
            vmax = maximum(vm3)

            va = vm3[1]
            vb = vm3[min(2, length(vm3))]
            vc = vm3[min(3, length(vm3))]

            push!(rows, (bus=string(bid), va=va, vb=vb, vc=vc, vmin=vmin, vmax=vmax))

        elseif haskey(sb, "vr") && haskey(sb, "vi")
            vr = sb["vr"]; vi = sb["vi"]
            if length(vr) >= 3 && length(vi) >= 3
                vm3 = sqrt.(vr[1:3].^2 .+ vi[1:3].^2)
                vmin = minimum(vm3)
                vmax = maximum(vm3)
                push!(rows, (bus=string(bid), va=vm3[1], vb=vm3[2], vc=vm3[3], vmin=vmin, vmax=vmax))
            end
        end
    end

    t = DataFrame(rows)
    sort!(t, :vmin)
    return t
end
