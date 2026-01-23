function plot_voltage_envelope(df::DataFrame; vmin_limit::Float64, vmax_limit::Float64, title_str::String, outpath::String)
    p = plot(df.time, df.vmin_pu;
        xlabel="Time",
        ylabel="Voltage (pu)",
        title=title_str,
        label="vmin"
    )
    plot!(df.time, df.vmax_pu; label="vmax", fillrange=df.vmin_pu, fillalpha=0.15)
    hline!([vmin_limit], linestyle=:dash, color=:red, label="lower limit")
    hline!([vmax_limit], linestyle=:dash, color=:red, label="upper limit")
    savefig(p, outpath)
end

function plot_series(x, y; xlabel::String, ylabel::String, title_str::String, outpath::String)
    p = plot(x, y; xlabel=xlabel, ylabel=ylabel, title=title_str, legend=false)
    savefig(p, outpath)
end

function plot_scatter(x, y; xlabel::String, ylabel::String, title_str::String, outpath::String)
    p = scatter(x, y; xlabel=xlabel, ylabel=ylabel, title=title_str, legend=false)
    savefig(p, outpath)
end

function plot_duration_curve(pct, x; xlabel::String, ylabel::String, title_str::String, outpath::String)
    p = plot(pct, x; xlabel=xlabel, ylabel=ylabel, title=title_str, legend=false)
    savefig(p, outpath)
end

function plot_worst_voltage_profile(volts::DataFrame; title_str::String, outpath::String)
    x = 1:nrow(volts)
    p = plot(x, volts.va; xlabel="Bus rank (sorted by vmin)", ylabel="Voltage (pu)", title=title_str, label="Va")
    plot!(x, volts.vb; label="Vb")
    plot!(x, volts.vc; label="Vc")
    savefig(p, outpath)
end
