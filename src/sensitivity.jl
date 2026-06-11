# Sensitivity-analysis orchestrators, ported from R/sensitivityresults.R.
# Results are returned as SensitivityResults tables (R returns tibbles).

# Extract [lb, ub] of the accepted region of a grid CS (R: min/max of accepted grid).
function _ci_from_grid(df)
    accepted = df.grid[df.accept .== 1]
    lb = isempty(accepted) ? Inf : minimum(accepted)
    ub = isempty(accepted) ? -Inf : maximum(accepted)
    return lb, ub
end

"""
    constructOriginalCS(betahat, sigma, numPrePeriods, numPostPeriods; l_vec=basisVector(1, numPostPeriods), alpha=0.05)

The conventional (non-robust) `1 - alpha` confidence interval for
`θ = l_vec' * τ_post`, i.e. the usual normal interval that assumes parallel
trends holds exactly. Use this as the baseline to compare against the robust
results from [`createSensitivityResults`](@ref) /
[`createSensitivityResults_relativeMagnitudes`](@ref), and as the `originalResults`
argument to the plotting functions.

Returns a single-row [`SensitivityResults`](@ref) table with columns `lb`, `ub`,
`method` (`"Original"`) and `Delta`.

# Example
```julia
constructOriginalCS(betahat, sigma, 3, 5)                 # first post-period effect
constructOriginalCS(betahat, sigma, 3, 5; l_vec = fill(0.2, 5))   # average effect
```
"""
function constructOriginalCS(betahat, sigma, numPrePeriods, numPostPeriods;
                             l_vec=basisVector(1, numPostPeriods), alpha=0.05)
    stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
    warnIfNotSymmPSD(sigma)
    bh = _asvec(betahat); lv = _asvec(l_vec)
    post = (numPrePeriods+1):(numPrePeriods+numPostPeriods)
    stdError = sqrt(_scalar(lv' * sigma[post, post] * lv))
    pe = dot(lv, bh[post])
    z = quantile(Normal(), 1 - alpha / 2)
    return SensitivityResults((lb=[pe - z * stdError], ub=[pe + z * stdError],
                              method=["Original"], Delta=Union{String,Missing}[missing]))
end

"""
    createSensitivityResults(betahat, sigma, numPrePeriods, numPostPeriods; kwargs...)

Sensitivity analysis under the smoothness restriction Δ^SD(M): report a robust
confidence interval for each value of the smoothness bound `M`, tracing how the
conclusion changes as one allows the differential trend to be increasingly
non-linear. The smallest `M` at which the interval first includes 0 is the
"breakdown" value for a significant effect.

# Arguments
 * `betahat`, `sigma`, `numPrePeriods`, `numPostPeriods`: the event study (ordered
  pre-periods then post-periods, reference omitted).

# Keyword arguments
 * `Mvec = nothing`: vector of `M` values. If `nothing`, a default grid from 0 to a
  data-driven upper bound ([`DeltaSD_upperBound_Mpre`](@ref)) is used.
 * `method = nothing`: confidence-set construction. `"FLCI"` (fixed-length, the
  default when no shape/sign restriction is given), `"Conditional"`, `"C-F"`
  (conditional-FLCI), or `"C-LF"` (conditional least-favorable).
 * `l_vec = basisVector(1, numPostPeriods)`: the target `θ = l_vec' * τ_post`.
 * `monotonicityDirection`: `"increasing"` / `"decreasing"` to additionally impose
  a shape restriction (uses Δ^SDM; default `method` becomes `"C-F"`).
 * `biasDirection`: `"positive"` / `"negative"` to additionally impose a sign
  restriction (uses Δ^SDB; default `method` becomes `"C-F"`).
 * `alpha = 0.05`: confidence level is `1 - alpha`.
 * `seed = 0`: RNG seed for the least-favorable simulation (only for `"C-LF"`).

# Returns
A [`SensitivityResults`](@ref) table with columns `lb`, `ub`, `method`, `Delta`,
`M` - one row per value of `M`.

# Example
```julia
res  = createSensitivityResults(betahat, sigma, 3, 5; Mvec = range(0, 0.06, length = 7))
orig = constructOriginalCS(betahat, sigma, 3, 5)
# using CairoMakie / Plots:
createSensitivityPlot(res, orig)
```
"""
function createSensitivityResults(betahat, sigma, numPrePeriods, numPostPeriods;
        method=nothing, Mvec=nothing, l_vec=basisVector(1, numPostPeriods),
        monotonicityDirection=nothing, biasDirection=nothing, alpha=0.05, seed=0)
    stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
    warnIfNotSymmPSD(sigma)
    bh = _asvec(betahat); lv = _asvec(l_vec)

    if Mvec === nothing
        if numPrePeriods == 1
            Mvec = collect(range(0.0, sqrt(sigma[1, 1]); length=10))
        else
            Mub = DeltaSD_upperBound_Mpre(bh, sigma, numPrePeriods; alpha=0.05)
            Mvec = collect(range(0.0, Mub; length=10))
        end
    end

    rows = NamedTuple[]
    if monotonicityDirection === nothing && biasDirection === nothing
        method === nothing && (method = "FLCI")
        Delta = "DeltaSD"
        for M in Mvec
            if method == "FLCI"
                temp = findOptimalFLCI(bh, sigma, M, numPrePeriods, numPostPeriods; l_vec=lv, alpha=alpha, seed=seed)
                push!(rows, (lb=temp.FLCI[1], ub=temp.FLCI[2], method="FLCI", Delta=Delta, M=M))
            else
                hf = method == "Conditional" ? "ARP" : (method == "C-F" ? "FLCI" : (method == "C-LF" ? "LF" : error("Method must equal one of: FLCI, Conditional, C-F or C-LF")))
                cs = computeConditionalCS_DeltaSD(bh, sigma, numPrePeriods, numPostPeriods;
                        l_vec=lv, alpha=alpha, M=M, hybrid_flag=hf, seed=seed)
                lb, ub = _ci_from_grid(cs)
                push!(rows, (lb=lb, ub=ub, method=method, Delta=Delta, M=M))
            end
        end
    elseif biasDirection !== nothing
        method === nothing && (method = "C-F")
        Delta = biasDirection == "positive" ? "DeltaSDPB" : "DeltaSDNB"
        for M in Mvec
            if method == "FLCI"
                @warn "You specified a sign restriction but method = FLCI. The FLCI does not use the sign restriction!"
                temp = findOptimalFLCI(bh, sigma, M, numPrePeriods, numPostPeriods; l_vec=lv, alpha=alpha, seed=seed)
                push!(rows, (lb=temp.FLCI[1], ub=temp.FLCI[2], method="FLCI", Delta=Delta, M=M))
            else
                hf = method == "Conditional" ? "ARP" : (method == "C-F" ? "FLCI" : (method == "C-LF" ? "LF" : error("Method must equal one of: FLCI, Conditional, C-F or C-LF")))
                cs = computeConditionalCS_DeltaSDB(bh, sigma, numPrePeriods, numPostPeriods;
                        l_vec=lv, alpha=alpha, M=M, biasDirection=biasDirection, hybrid_flag=hf, seed=seed)
                lb, ub = _ci_from_grid(cs)
                push!(rows, (lb=lb, ub=ub, method=method, Delta=Delta, M=M))
            end
        end
    else
        method === nothing && (method = "C-F")
        Delta = monotonicityDirection == "increasing" ? "DeltaSDI" : "DeltaSDD"
        for M in Mvec
            if method == "FLCI"
                @warn "You specified a shape restriction but method = FLCI. The FLCI does not use the shape restriction!"
                temp = findOptimalFLCI(bh, sigma, M, numPrePeriods, numPostPeriods; l_vec=lv, alpha=alpha, seed=seed)
                push!(rows, (lb=temp.FLCI[1], ub=temp.FLCI[2], method="FLCI", Delta=Delta, M=M))
            else
                hf = method == "Conditional" ? "ARP" : (method == "C-F" ? "FLCI" : (method == "C-LF" ? "LF" : error("Method must equal one of: FLCI, Conditional, C-F or C-LF")))
                cs = computeConditionalCS_DeltaSDM(bh, sigma, numPrePeriods, numPostPeriods;
                        l_vec=lv, alpha=alpha, M=M, monotonicityDirection=monotonicityDirection, hybrid_flag=hf, seed=seed)
                lb, ub = _ci_from_grid(cs)
                push!(rows, (lb=lb, ub=ub, method=method, Delta=Delta, M=M))
            end
        end
    end
    return _results_from_rows(rows)
end

"""
    createSensitivityResults_relativeMagnitudes(betahat, sigma, numPrePeriods, numPostPeriods; kwargs...)

Sensitivity analysis under the relative-magnitudes restriction Δ^RM(M̄)
(or Δ^SDRM(M̄)): report a robust confidence interval for each value of `M̄`, where
`M̄` bounds how large post-treatment violations of parallel trends can be relative
to the pre-treatment violations. `M̄ = 1` means "no larger than the worst
pre-period violation"; the breakdown `M̄` is the largest value for which the effect
stays significant.

# Keyword arguments
 * `Mbarvec = nothing`: vector of `M̄` values (default: 10 points on `[0, 2]`).
 * `bound = "deviation from parallel trends"`: use Δ^RM. Set to
  `"deviation from linear trend"` to instead bound the relative magnitude of the
  *non-linear* component (Δ^SDRM; requires `numPrePeriods >= 2`).
 * `method = "C-LF"`: `"C-LF"` (conditional least-favorable) or `"Conditional"`.
 * `l_vec = basisVector(1, numPostPeriods)`: the target `θ = l_vec' * τ_post`.
 * `monotonicityDirection` / `biasDirection`: optionally add a shape or sign
  restriction (mutually exclusive).
 * `alpha = 0.05`, `gridPoints = 1000`, `grid_lb`, `grid_ub`, `seed = 0`: as in
  [`computeConditionalCS_DeltaRM`](@ref).

# Returns
A [`SensitivityResults`](@ref) table with columns `lb`, `ub`, `method`, `Delta`,
`Mbar` - one row per value of `M̄`.

# Example
```julia
res  = createSensitivityResults_relativeMagnitudes(betahat, sigma, 3, 5; Mbarvec = 0.5:0.5:2)
orig = constructOriginalCS(betahat, sigma, 3, 5)
createSensitivityPlot_relativeMagnitudes(res, orig)
```
"""
function createSensitivityResults_relativeMagnitudes(betahat, sigma, numPrePeriods, numPostPeriods;
        bound="deviation from parallel trends", method="C-LF", Mbarvec=nothing,
        l_vec=basisVector(1, numPostPeriods), monotonicityDirection=nothing, biasDirection=nothing,
        alpha=0.05, gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
    warnIfNotSymmPSD(sigma)
    bh = _asvec(betahat); lv = _asvec(l_vec)

    Mbarvec === nothing && (Mbarvec = collect(range(0.0, 2.0; length=10)))
    (bound == "deviation from parallel trends" || bound == "deviation from linear trend") ||
        error("bound must equal either 'deviation from parallel trends' or 'deviation from linear trend'.")
    !(monotonicityDirection !== nothing && biasDirection !== nothing) ||
        error("Please select either a shape restriction or sign restriction (not both).")

    if method == "C-LF"
        hybrid_flag = "LF"; method_named = "C-LF"
    elseif method == "Conditional"
        hybrid_flag = "ARP"; method_named = "Conditional"
    else
        error("method must be either NULL, Conditional or C-LF.")
    end

    runcs(f, Delta) = begin
        rows = NamedTuple[]
        for Mbar in Mbarvec
            cs = f(Mbar)
            lb, ub = _ci_from_grid(cs)
            push!(rows, (lb=lb, ub=ub, method=method_named, Delta=Delta, Mbar=Mbar))
        end
        _results_from_rows(rows)
    end

    if bound == "deviation from parallel trends"
        if monotonicityDirection === nothing && biasDirection === nothing
            return runcs(Mbar -> computeConditionalCS_DeltaRM(bh, sigma, numPrePeriods, numPostPeriods;
                l_vec=lv, alpha=alpha, Mbar=Mbar, hybrid_flag=hybrid_flag,
                gridPoints=gridPoints, grid_ub=grid_ub, grid_lb=grid_lb, seed=seed), "DeltaRM")
        elseif monotonicityDirection !== nothing
            Delta = monotonicityDirection == "increasing" ? "DeltaRMI" :
                    (monotonicityDirection == "decreasing" ? "DeltaRMD" : error("monotonicityDirection must equal either increasing or decreasing."))
            return runcs(Mbar -> computeConditionalCS_DeltaRMM(bh, sigma, numPrePeriods, numPostPeriods;
                l_vec=lv, alpha=alpha, Mbar=Mbar, monotonicityDirection=monotonicityDirection,
                hybrid_flag=hybrid_flag, gridPoints=gridPoints, grid_ub=grid_ub, grid_lb=grid_lb, seed=seed), Delta)
        else
            Delta = biasDirection == "positive" ? "DeltaRMPB" :
                    (biasDirection == "negative" ? "DeltaRMNB" : error("biasDirection must equal either positive or negative."))
            return runcs(Mbar -> computeConditionalCS_DeltaRMB(bh, sigma, numPrePeriods, numPostPeriods;
                l_vec=lv, alpha=alpha, Mbar=Mbar, biasDirection=biasDirection,
                hybrid_flag=hybrid_flag, gridPoints=gridPoints, grid_ub=grid_ub, grid_lb=grid_lb, seed=seed), Delta)
        end
    else
        numPrePeriods == 1 && error("Error: not enough pre-periods for 'deviation from linear trend' (Delta^{SDRM} as base choice).")
        if monotonicityDirection === nothing && biasDirection === nothing
            return runcs(Mbar -> computeConditionalCS_DeltaSDRM(bh, sigma, numPrePeriods, numPostPeriods;
                l_vec=lv, alpha=alpha, Mbar=Mbar, hybrid_flag=hybrid_flag,
                gridPoints=gridPoints, grid_ub=grid_ub, grid_lb=grid_lb, seed=seed), "DeltaSDRM")
        elseif monotonicityDirection !== nothing
            Delta = monotonicityDirection == "increasing" ? "DeltaSDRMI" :
                    (monotonicityDirection == "decreasing" ? "DeltaSDRMD" : error("monotonicityDirection must equal either increasing or decreasing."))
            return runcs(Mbar -> computeConditionalCS_DeltaSDRMM(bh, sigma, numPrePeriods, numPostPeriods;
                l_vec=lv, alpha=alpha, Mbar=Mbar, monotonicityDirection=monotonicityDirection,
                hybrid_flag=hybrid_flag, gridPoints=gridPoints, grid_ub=grid_ub, grid_lb=grid_lb, seed=seed), Delta)
        else
            Delta = biasDirection == "positive" ? "DeltaSDRMPB" :
                    (biasDirection == "negative" ? "DeltaSDRMNB" : error("biasDirection must equal either positive or negative."))
            return runcs(Mbar -> computeConditionalCS_DeltaSDRMB(bh, sigma, numPrePeriods, numPostPeriods;
                l_vec=lv, alpha=alpha, Mbar=Mbar, biasDirection=biasDirection,
                hybrid_flag=hybrid_flag, gridPoints=gridPoints, grid_ub=grid_ub, grid_lb=grid_lb, seed=seed), Delta)
        end
    end
end
