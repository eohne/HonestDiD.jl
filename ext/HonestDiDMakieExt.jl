module HonestDiDMakieExt

# Makie.jl rendering for the HonestDiD plot functions (loaded with any Makie
# backend, e.g. CairoMakie / GLMakie). Data prep lives in HonestDiD core.

using HonestDiD
using Makie

const HD = HonestDiD

function _sensitivity(robustResults, originalResults, xlbl, xcol; rescaleFactor, maxX, add_xAxis)
    d = HD._sensitivity_plot_data(robustResults, originalResults, xcol; rescaleFactor=rescaleFactor, maxX=maxX)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; xlabel=xlbl, ylabel="")
    add_xAxis && Makie.hlines!(ax, [0.0]; color=:black, linestyle=:dash)
    for m in unique(d.method)
        idx = d.method .== m
        color = m == "Original" ? :red : Makie.RGBf(0.004, 0.635, 0.851)  # #01a2d9
        mid = (d.lb[idx] .+ d.ub[idx]) ./ 2
        Makie.rangebars!(ax, d.x[idx], d.lb[idx], d.ub[idx]; color=color)
        Makie.scatter!(ax, d.x[idx], mid; color=color, markersize=8, label=m)
    end
    Makie.axislegend(ax; position=:rb)
    return fig
end

HD.createSensitivityPlot(robustResults, originalResults; rescaleFactor=1, maxM=Inf, add_xAxis=true) =
    _sensitivity(robustResults, originalResults, "M", :M; rescaleFactor=rescaleFactor, maxX=maxM, add_xAxis=add_xAxis)

HD.createSensitivityPlot_relativeMagnitudes(robustResults, originalResults; rescaleFactor=1, maxMbar=Inf, add_xAxis=true) =
    _sensitivity(robustResults, originalResults, "Mbar", :Mbar; rescaleFactor=rescaleFactor, maxX=maxMbar, add_xAxis=add_xAxis)

function HD.createEventStudyPlot(betahat; stdErrors=nothing, sigma=nothing,
        numPrePeriods, numPostPeriods, alpha=0.05, timeVec, referencePeriod, useRelativeEventTime=false)
    d = HD._eventstudy_plot_data(betahat; stdErrors=stdErrors, sigma=sigma, numPrePeriods=numPrePeriods,
        numPostPeriods=numPostPeriods, alpha=alpha, timeVec=timeVec, referencePeriod=referencePeriod,
        useRelativeEventTime=useRelativeEventTime)
    fin = .!isnan.(d.se)
    blue = Makie.RGBf(0.004, 0.635, 0.851)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; xlabel="Event time", ylabel="")
    Makie.rangebars!(ax, d.t[fin], d.beta[fin] .- d.z .* d.se[fin], d.beta[fin] .+ d.z .* d.se[fin]; color=blue)
    Makie.scatter!(ax, d.t, d.beta; color=:red, markersize=10)
    return fig
end

end # module
