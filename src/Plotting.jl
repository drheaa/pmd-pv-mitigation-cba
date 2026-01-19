module Plotting

using Plots

export plot_voltage_along_feeder,
       plot_voltage_histogram

"""
Plot voltage along feeder, using bus distances and phase voltages.
Expects:
buses[bus] = Dict("vma","vmb","vmc","distance")
edges_df with Bus1, Bus2, length_km, kind
"""
function plot_voltage_along_feeder(buses, edges_df; vmin=0.94*230, vmax=1.10*230)

    p = plot(
        legend=false,
        xlabel="Distance from source (km)",
        ylabel="Voltage (V)",
        title="Voltage along feeder (snapshot)"
    )

    colors = Dict(1 => :blue, 2 => :red, 3 => :black)

    line_edges = edges_df[edges_df.kind .== "line", :]

    for r in eachrow(line_edges)
        b1, b2 = r.Bus1, r.Bus2

        if !(haskey(buses, b1) && haskey(buses, b2))
            continue
        end
        if !(haskey(buses[b1], "distance") && haskey(buses[b2], "distance"))
            continue
        end

        for ph in (1, 2, 3)
            v1 = ph == 1 ? buses[b1]["vma"] :
                 ph == 2 ? buses[b1]["vmb"] :
                           buses[b1]["vmc"]

            v2 = ph == 1 ? buses[b2]["vma"] :
                 ph == 2 ? buses[b2]["vmb"] :
                           buses[b2]["vmc"]

            plot!(
                p,
                [buses[b1]["distance"], buses[b2]["distance"]],
                [v1, v2],
                color = colors[ph],
                linewidth = 1,
                marker = :circle,
                markersize = 1
            )
        end
    end

    maxdist = maximum(get(b, "distance", 0.0) for b in values(buses))
    plot!(p, [0, maxdist], [vmin, vmin], linestyle=:dash, color=:red)
    plot!(p, [0, maxdist], [vmax, vmax], linestyle=:dash, color=:red)

    return p
end


"""
Histogram of phase voltages across all buses.
"""
function plot_voltage_histogram(buses; vmin=0.94*230, vmax=1.10*230)

    va = [b["vma"] for b in values(buses)]
    vb = [b["vmb"] for b in values(buses)]
    vc = [b["vmc"] for b in values(buses)]

    bins = (vmin - 2):0.5:(vmax + 2)

    p = histogram(va; bins=bins, label="A")
    histogram!(p, vb; bins=bins, label="B")
    histogram!(p, vc; bins=bins, label="C")

    xlabel!(p, "Voltage (V)")
    ylabel!(p, "Count")
    title!(p, "Voltage histogram (snapshot)")

    return p
end

end
