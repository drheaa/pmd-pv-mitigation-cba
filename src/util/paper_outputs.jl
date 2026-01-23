function duration_curve(values::AbstractVector{<:Real}; high_is_worse::Bool=true)
    v = Float64.(values)
    v = filter(isfinite, v)
    isempty(v) && return (pct=Float64[], x=Float64[])
    x = sort(v; rev=high_is_worse)
    pct = collect(range(0, 100; length=length(x)))
    return (pct=pct, x=x)
end

"""
write_duration_curves_csv(outpath, df)

Writes a standard duration curve table.
This is useful for paper figures and consistency across cases.
"""
function write_duration_curves_csv(outpath::String, df::DataFrame)
    dc_vmin = duration_curve(df.vmin_pu; high_is_worse=false)
    dc_vuf  = duration_curve(df.vuf_max; high_is_worse=true)

    n = maximum(length.([dc_vmin.x, dc_vuf.x]))

    function pad(v, n)
        w = Vector{Float64}(undef, n)
        fill!(w, NaN)
        w[1:length(v)] .= v
        return w
    end

    out = DataFrame(
        pct = collect(range(0, 100; length=n)),
        vmin_pu_sorted = pad(sort(filter(isfinite, Float64.(df.vmin_pu)); rev=false), n),
        vuf_sorted = pad(sort(filter(isfinite, Float64.(df.vuf_max)); rev=true), n),
    )

    if :i_max_phase in names(df) && any(isfinite.(df.i_max_phase))
        out.i_max_phase_sorted = pad(sort(filter(isfinite, Float64.(df.i_max_phase)); rev=true), n)
    end

    CSV.write(outpath, out)
    return out
end

"""
select_worst_timesteps(df, n_worst)

Selects n_worst rows by minimum vmin_pu.
"""
function select_worst_timesteps(df::DataFrame, n_worst::Int)
    n = min(n_worst, nrow(df))
    return sort(df, :vmin_pu)[1:n, :]
end
