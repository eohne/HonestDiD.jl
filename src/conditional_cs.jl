# computeConditionalCS_Delta* : confidence sets for each Delta restriction,
# ported from R/deltasd.R, deltasdb.R, deltasdm.R, deltarm.R, deltarmb.R,
# deltarmm.R, deltasdrm.R, deltasdrmb.R, deltasdrmm.R.
#
# Each returns a ConditionalCS(grid, accept) (or a Float64 length if returnLength).

_postidx(numPre, numPost) = (numPre+1):(numPre+numPost)

function _sdTheta(sigma, lv, numPre, numPost)
    post = _postidx(numPre, numPost)
    return sqrt(_scalar(lv' * sigma[post, post] * lv))
end

# vbar projection for the FLCI hybrid: min-norm minimizer of ||flci_l - A' vbar||^2,
# replacing R's CVXR quad-form solve (the minimizer is non-unique; min-norm is the
# canonical deterministic choice).
_project_vbar(A, flci_l) = pinv(Matrix(A')) * _asvec(flci_l)

# Delta^SD
"""
    computeConditionalCS_DeltaSD(betahat, sigma, numPrePeriods, numPostPeriods; kwargs...)

Confidence set for `l_vec' * τ_post` under the smoothness restriction
Δ^SD(M): the slope of the differential trend may change by at most `M` between
consecutive periods. Inference uses the conditional / hybrid moment-inequality
test of Andrews, Roth & Pakes (2019), inverted over a grid of candidate values of
the target parameter.

# Arguments
 * `betahat`: event-study coefficients, ordered pre-periods then post-periods
  (the reference period omitted).
 * `sigma`: covariance matrix of `betahat` (same ordering).
 * `numPrePeriods`, `numPostPeriods`: number of pre- and post-treatment coefficients.

# Keyword arguments
 * `l_vec = basisVector(1, numPostPeriods)`: defines `θ = l_vec' * τ_post`
  (default: the first post-period effect).
 * `M = 0`: smoothness bound; `M = 0` forces an exactly linear differential trend.
 * `alpha = 0.05`: confidence level is `1 - alpha`.
 * `hybrid_flag = "FLCI"`: first-stage hybrid - `"FLCI"` (conditional-FLCI),
  `"LF"` (conditional least-favorable) or `"ARP"` (pure conditional, no hybrid).
 * `hybrid_kappa = alpha/10`: size of the first-stage hybrid test.
 * `gridPoints = 1000`: number of grid points for the test inversion.
 * `grid_lb`, `grid_ub`: grid bounds (default: derived from the identified set / FLCI).
 * `postPeriodMomentsOnly = true`: drop moments involving only pre-period coefficients.
 * `returnLength = false`: if `true`, return the CI length (`Float64`) instead.
 * `seed = 0`: RNG seed for the least-favorable simulation (used only for `"LF"`).

# Returns
A [`ConditionalCS`](@ref) (`grid`, `accept`); use [`confidence_interval`](@ref)
to obtain `(lb, ub)`. With `returnLength = true`, returns the length as a `Float64`.

# Example
```julia
cs = computeConditionalCS_DeltaSD(betahat, sigma, 3, 5; M = 0.02, hybrid_flag = "FLCI")
confidence_interval(cs)
```

See also [`createSensitivityResults`](@ref) (which sweeps `M`),
[`computeConditionalCS_DeltaSDB`](@ref), [`computeConditionalCS_DeltaSDM`](@ref).
"""
function computeConditionalCS_DeltaSD(betahat, sigma, numPrePeriods, numPostPeriods;
        l_vec=basisVector(1, numPostPeriods), M=0, alpha=0.05, hybrid_flag="FLCI",
        hybrid_kappa=alpha / 10, returnLength=false, postPeriodMomentsOnly=true,
        gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    bh = _asvec(betahat); lv = _asvec(l_vec)
    A_SD = create_A_SD(numPrePeriods, numPostPeriods)
    d_SD = create_d_SD(numPrePeriods, numPostPeriods, M)

    if postPeriodMomentsOnly && numPostPeriods > 1
        postPeriodIndices = (numPrePeriods+1):size(A_SD, 2)
        rowsForARP = findall(r -> sum(A_SD[r, postPeriodIndices] .!= 0) > 0, 1:size(A_SD, 1))
    else
        rowsForARP = collect(1:size(A_SD, 1))
    end

    hybrid_list = Dict{Symbol,Any}(:hybrid_kappa => hybrid_kappa)

    if numPostPeriods == 1
        if hybrid_flag == "FLCI"
            flci = findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, lv; alpha=hybrid_kappa)
            hybrid_list[:flci_l] = flci.optimalVec
            hybrid_list[:flci_halflength] = flci.optimalHalfLength
            pe = dot(flci.optimalVec, bh)
            grid_ub === nothing && (grid_ub = pe + flci.optimalHalfLength)
            grid_lb === nothing && (grid_lb = pe - flci.optimalHalfLength)
        elseif hybrid_flag == "LF"
            hybrid_list[:lf_cv] = compute_least_favorable_cv(nothing, A_SD * sigma * A_SD', hybrid_kappa; seed=seed)
            if grid_ub === nothing && grid_lb === nothing
                idlb, idub = compute_IDset_DeltaSD(M, zeros(numPrePeriods + numPostPeriods), lv, numPrePeriods, numPostPeriods)
                sd = _sdTheta(sigma, lv, numPrePeriods, numPostPeriods)
                grid_ub = idub + 20 * sd; grid_lb = idlb - 20 * sd
            end
        elseif hybrid_flag == "ARP"
            if grid_ub === nothing && grid_lb === nothing
                idlb, idub = compute_IDset_DeltaSD(M, zeros(numPrePeriods + numPostPeriods), lv, numPrePeriods, numPostPeriods)
                sd = _sdTheta(sigma, lv, numPrePeriods, numPostPeriods)
                grid_ub = idub + 20 * sd; grid_lb = idlb - 20 * sd
            end
        else
            error("hybrid_flag must equal 'ARP' or 'FLCI' or 'LF'")
        end
        return APR_computeCI_NoNuis(bh, sigma, A_SD, d_SD, numPrePeriods, numPostPeriods, lv,
                                    alpha, returnLength, hybrid_flag, hybrid_list, grid_ub, grid_lb, gridPoints)
    else
        if hybrid_flag == "FLCI"
            flci = findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, lv; alpha=hybrid_kappa)
            hybrid_list[:flci_l] = flci.optimalVec
            hybrid_list[:vbar] = _project_vbar(A_SD, flci.optimalVec)
            hybrid_list[:flci_halflength] = flci.optimalHalfLength
            pe = dot(flci.optimalVec, bh)
            grid_ub === nothing && (grid_ub = pe + flci.optimalHalfLength)
            grid_lb === nothing && (grid_lb = pe - flci.optimalHalfLength)
        else
            idlb, idub = compute_IDset_DeltaSD(M, zeros(numPrePeriods + numPostPeriods), lv, numPrePeriods, numPostPeriods)
            sd = _sdTheta(sigma, lv, numPrePeriods, numPostPeriods)
            grid_ub === nothing && (grid_ub = idub + 20 * sd)
            grid_lb === nothing && (grid_lb = idlb - 20 * sd)
        end
        return ARP_computeCI(bh, sigma, numPrePeriods, numPostPeriods, A_SD, d_SD, lv, alpha,
                             hybrid_flag, hybrid_list, returnLength, grid_lb, grid_ub, gridPoints, rowsForARP)
    end
end

# Delta^SDB
"""
    computeConditionalCS_DeltaSDB(betahat, sigma, numPrePeriods, numPostPeriods; biasDirection="positive", kwargs...)

Confidence set under Δ^SDB(M) = the smoothness restriction Δ^SD(M)
intersected with a *sign* restriction on the post-period differential trend.
`biasDirection = "positive"` requires the post-period bias to be non-negative,
`"negative"` non-positive. All other arguments are as in
[`computeConditionalCS_DeltaSD`](@ref).
"""
function computeConditionalCS_DeltaSDB(betahat, sigma, numPrePeriods, numPostPeriods;
        M=0, l_vec=basisVector(1, numPostPeriods), alpha=0.05, hybrid_flag="FLCI",
        hybrid_kappa=alpha / 10, returnLength=false, biasDirection="positive",
        postPeriodMomentsOnly=true, gridPoints=10^3, grid_lb=nothing, grid_ub=nothing, seed=0)
    bh = _asvec(betahat); lv = _asvec(l_vec)
    A_SDB = create_A_SDB(numPrePeriods, numPostPeriods; biasDirection=biasDirection)
    d_SDB = create_d_SDB(numPrePeriods, numPostPeriods, M)

    if postPeriodMomentsOnly && numPostPeriods > 1
        postPeriodIndices = (numPrePeriods+1):size(A_SDB, 2)
        rowsForARP = findall(r -> sum(A_SDB[r, postPeriodIndices] .!= 0) > 0, 1:size(A_SDB, 1))
    else
        rowsForARP = collect(1:size(A_SDB, 1))
    end

    hybrid_list = Dict{Symbol,Any}(:hybrid_kappa => hybrid_kappa)

    if numPostPeriods == 1
        if hybrid_flag == "FLCI"
            flci = findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, lv; alpha=hybrid_kappa)
            hybrid_list[:flci_l] = flci.optimalVec
            hybrid_list[:flci_halflength] = flci.optimalHalfLength
            pe = dot(flci.optimalVec, bh)
            grid_ub === nothing && (grid_ub = pe + flci.optimalHalfLength)
            grid_lb === nothing && (grid_lb = pe - flci.optimalHalfLength)
        elseif hybrid_flag == "LF" || hybrid_flag == "ARP"
            if hybrid_flag == "LF"
                hybrid_list[:lf_cv] = compute_least_favorable_cv(nothing, A_SDB * sigma * A_SDB', hybrid_kappa; seed=seed)
            end
            if grid_ub === nothing && grid_lb === nothing
                idlb, idub = compute_IDset_DeltaSDB(M, zeros(numPrePeriods + numPostPeriods), lv, numPrePeriods, numPostPeriods, biasDirection)
                sd = _sdTheta(sigma, lv, numPrePeriods, numPostPeriods)
                grid_ub = idub + 20 * sd; grid_lb = idlb - 20 * sd
            end
        else
            error("hybrid_flag must equal 'ARP' or 'FLCI' or 'LF'")
        end
        return APR_computeCI_NoNuis(bh, sigma, A_SDB, d_SDB, numPrePeriods, numPostPeriods, lv,
                                    alpha, returnLength, hybrid_flag, hybrid_list, grid_ub, grid_lb, gridPoints)
    else
        if hybrid_flag == "FLCI"
            flci = findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, lv; alpha=hybrid_kappa)
            hybrid_list[:flci_l] = flci.optimalVec
            hybrid_list[:vbar] = _project_vbar(A_SDB, flci.optimalVec)
            hybrid_list[:flci_halflength] = flci.optimalHalfLength
            pe = dot(flci.optimalVec, bh)
            grid_ub === nothing && (grid_ub = pe + flci.optimalHalfLength)
            grid_lb === nothing && (grid_lb = pe - flci.optimalHalfLength)
        else
            idlb, idub = compute_IDset_DeltaSDB(M, zeros(numPrePeriods + numPostPeriods), lv, numPrePeriods, numPostPeriods, biasDirection)
            if biasDirection == "negative"
                idlb, idub = -idub, -idlb
            end
            sd = _sdTheta(sigma, lv, numPrePeriods, numPostPeriods)
            grid_ub === nothing && (grid_ub = idub + 20 * sd)
            grid_lb === nothing && (grid_lb = idlb - 20 * sd)
        end
        return ARP_computeCI(bh, sigma, numPrePeriods, numPostPeriods, A_SDB, d_SDB, lv, alpha,
                             hybrid_flag, hybrid_list, returnLength, grid_lb, grid_ub, gridPoints, rowsForARP)
    end
end

# Delta^SDM
"""
    computeConditionalCS_DeltaSDM(betahat, sigma, numPrePeriods, numPostPeriods; monotonicityDirection="increasing", kwargs...)

Confidence set under Δ^SDM(M) = the smoothness restriction Δ^SD(M)
intersected with a *monotonicity* restriction: the differential trend is
`"increasing"` or `"decreasing"`. All other arguments are as in
[`computeConditionalCS_DeltaSD`](@ref).
"""
function computeConditionalCS_DeltaSDM(betahat, sigma, numPrePeriods, numPostPeriods;
        M=0, l_vec=basisVector(1, numPostPeriods), alpha=0.05, monotonicityDirection="increasing",
        hybrid_flag="FLCI", hybrid_kappa=alpha / 10, returnLength=false,
        postPeriodMomentsOnly=true, gridPoints=10^3, grid_lb=nothing, grid_ub=nothing, seed=0)
    bh = _asvec(betahat); lv = _asvec(l_vec)
    A_SDM = create_A_SDM(numPrePeriods, numPostPeriods; monotonicityDirection=monotonicityDirection)
    d_SDM = create_d_SDM(numPrePeriods, numPostPeriods, M)

    if postPeriodMomentsOnly && numPostPeriods > 1
        postPeriodIndices = (numPrePeriods+1):size(A_SDM, 2)
        rowsForARP = findall(r -> sum(A_SDM[r, postPeriodIndices] .!= 0) > 0, 1:size(A_SDM, 1))
    else
        rowsForARP = collect(1:size(A_SDM, 1))
    end

    hybrid_list = Dict{Symbol,Any}(:hybrid_kappa => hybrid_kappa)

    if numPostPeriods == 1
        if hybrid_flag == "FLCI"
            flci = findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, lv; alpha=hybrid_kappa)
            hybrid_list[:flci_l] = flci.optimalVec
            hybrid_list[:flci_halflength] = flci.optimalHalfLength
            pe = dot(flci.optimalVec, bh)
            grid_ub === nothing && (grid_ub = pe + flci.optimalHalfLength)
            grid_lb === nothing && (grid_lb = pe - flci.optimalHalfLength)
        elseif hybrid_flag == "LF" || hybrid_flag == "ARP"
            if hybrid_flag == "LF"
                hybrid_list[:lf_cv] = compute_least_favorable_cv(nothing, A_SDM * sigma * A_SDM', hybrid_kappa; seed=seed)
            end
            if grid_ub === nothing && grid_lb === nothing
                idlb, idub = compute_IDset_DeltaSDM(M, zeros(numPrePeriods + numPostPeriods), lv, numPrePeriods, numPostPeriods, monotonicityDirection)
                sd = _sdTheta(sigma, lv, numPrePeriods, numPostPeriods)
                grid_ub = idub + 20 * sd; grid_lb = idlb - 20 * sd
            end
        else
            error("hybrid_flag must equal 'ARP' or 'FLCI' or 'LF'")
        end
        return APR_computeCI_NoNuis(bh, sigma, A_SDM, d_SDM, numPrePeriods, numPostPeriods, lv,
                                    alpha, returnLength, hybrid_flag, hybrid_list, grid_ub, grid_lb, gridPoints)
    else
        if hybrid_flag == "FLCI"
            flci = findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, lv; alpha=hybrid_kappa)
            hybrid_list[:flci_l] = flci.optimalVec
            hybrid_list[:vbar] = _project_vbar(A_SDM, flci.optimalVec)
            hybrid_list[:flci_halflength] = flci.optimalHalfLength
            # NOTE: R overwrites the grid here regardless of user input.
            pe = dot(flci.optimalVec, bh)
            grid_ub = pe + flci.optimalHalfLength
            grid_lb = pe - flci.optimalHalfLength
        else
            idlb, idub = compute_IDset_DeltaSDM(M, zeros(numPrePeriods + numPostPeriods), lv, numPrePeriods, numPostPeriods, monotonicityDirection)
            sd = _sdTheta(sigma, lv, numPrePeriods, numPostPeriods)
            grid_ub === nothing && (grid_ub = idub + 20 * sd)
            grid_lb === nothing && (grid_lb = idlb - 20 * sd)
        end
        return ARP_computeCI(bh, sigma, numPrePeriods, numPostPeriods, A_SDM, d_SDM, lv, alpha,
                             hybrid_flag, hybrid_list, returnLength, grid_lb, grid_ub, gridPoints, rowsForARP)
    end
end

# Relative-magnitude family (union over s, ±)
# One fixed-s, fixed-sign confidence band; returns the `accept` vector.
function _rm_fixedS_accept(A_s, d_s, bh, sigma, numPre, numPost, lv, alpha,
                           hybrid_flag, hybrid_kappa, postPeriodMomentsOnly, gridPoints, grid_ub, grid_lb, seed)
    (hybrid_flag == "LF" || hybrid_flag == "ARP") || error("hybrid_flag must equal 'ARP' or 'LF'")
    A = copy(A_s); d = copy(d_s)
    if postPeriodMomentsOnly
        if numPost > 1
            postPeriodIndices = (numPre+1):size(A, 2)
            rowsForARP = findall(r -> sum(A[r, postPeriodIndices] .!= 0) > 0, 1:size(A, 1))
        else
            postPeriodRows = findall(!=(0), A[:, end])
            A = A[postPeriodRows, :]; d = d[postPeriodRows]
            rowsForARP = collect(1:size(A, 1))
        end
    else
        rowsForARP = collect(1:size(A, 1))
    end

    hybrid_list = Dict{Symbol,Any}(:hybrid_kappa => hybrid_kappa)
    if numPost == 1
        if hybrid_flag == "LF"
            hybrid_list[:lf_cv] = compute_least_favorable_cv(nothing, A * sigma * A', hybrid_kappa; seed=seed)
        end
        ci = APR_computeCI_NoNuis(bh, sigma, A, d, numPre, numPost, lv, alpha, false,
                                  hybrid_flag, hybrid_list, grid_ub, grid_lb, gridPoints)
    else
        ci = ARP_computeCI(bh, sigma, numPre, numPost, A, d, lv, alpha, hybrid_flag, hybrid_list,
                           false, grid_lb, grid_ub, gridPoints, rowsForARP)
    end
    return ci.accept
end

# Generic union over s and ± sign. `build(s, max_positive)` returns (A, d).
function _rm_family_cs(build, min_s, bh, sigma, numPre, numPost, lv, Mbar, alpha,
                       hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly,
                       gridPoints, grid_ub, grid_lb, seed)
    s_indices = min_s:0
    sd = _sdTheta(sigma, lv, numPre, numPost)
    gub = grid_ub === nothing ? 20 * sd : grid_ub
    glb = grid_lb === nothing ? -20 * sd : grid_lb

    plus = zeros(gridPoints, length(s_indices))
    minus = zeros(gridPoints, length(s_indices))
    for (i, s) in enumerate(s_indices)
        Ap, dp = build(s, true)
        plus[:, i] = _rm_fixedS_accept(Ap, dp, bh, sigma, numPre, numPost, lv, alpha,
                                       hybrid_flag, hybrid_kappa, postPeriodMomentsOnly, gridPoints, gub, glb, seed)
        Am, dm = build(s, false)
        minus[:, i] = _rm_fixedS_accept(Am, dm, bh, sigma, numPre, numPost, lv, alpha,
                                        hybrid_flag, hybrid_kappa, postPeriodMomentsOnly, gridPoints, gub, glb, seed)
    end
    plus_max = vec(maximum(plus; dims=2))
    minus_max = vec(maximum(minus; dims=2))
    grid = collect(range(glb, gub; length=gridPoints))
    accept = max.(plus_max, minus_max)
    if returnLength
        gridLength = 0.5 .* (vcat(0.0, diff(grid)) .+ vcat(diff(grid), 0.0))
        return sum(accept .* gridLength)
    else
        return ConditionalCS(grid, accept)
    end
end

"""
    computeConditionalCS_DeltaRM(betahat, sigma, numPrePeriods, numPostPeriods; Mbar=0, kwargs...)

Confidence set under the relative-magnitudes restriction Δ^RM(M̄): the
largest post-treatment violation of parallel trends (between consecutive periods)
is at most `Mbar` times the largest pre-treatment violation. `Mbar = 1` means the
post-period violation is no larger than the worst pre-period violation.

Unlike the smoothness families, `hybrid_flag` defaults to `"LF"` and must be
`"LF"` or `"ARP"` (the FLCI hybrid is not available). The confidence set is the
union over the choice of which period attains the maximum pre-period violation
(and its sign). Other arguments are as in [`computeConditionalCS_DeltaSD`](@ref),
with `Mbar` in place of `M`.

See also [`createSensitivityResults_relativeMagnitudes`](@ref) (sweeps `Mbar`),
[`computeConditionalCS_DeltaSDRM`](@ref) (relative magnitudes of the *non-linear*
component).
"""
function computeConditionalCS_DeltaRM(betahat, sigma, numPrePeriods, numPostPeriods;
        l_vec=basisVector(1, numPostPeriods), Mbar=0, alpha=0.05, hybrid_flag="LF",
        hybrid_kappa=alpha / 10, returnLength=false, postPeriodMomentsOnly=true,
        gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    bh = _asvec(betahat); lv = _asvec(l_vec)
    build = (s, mp) -> (create_A_RM(numPrePeriods, numPostPeriods; Mbar=Mbar, s=s, max_positive=mp),
                        create_d_RM(numPrePeriods, numPostPeriods))
    return _rm_family_cs(build, -(numPrePeriods - 1), bh, sigma, numPrePeriods, numPostPeriods, lv,
                         Mbar, alpha, hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly,
                         gridPoints, grid_ub, grid_lb, seed)
end

"""
    computeConditionalCS_DeltaRMB(betahat, sigma, numPrePeriods, numPostPeriods; Mbar=0, biasDirection="positive", kwargs...)

Confidence set under Δ^RMB(M̄) = the relative-magnitudes restriction Δ^RM(M̄)
intersected with a sign restriction (`biasDirection = "positive"` / `"negative"`).
Other arguments as in [`computeConditionalCS_DeltaRM`](@ref).
"""
function computeConditionalCS_DeltaRMB(betahat, sigma, numPrePeriods, numPostPeriods;
        l_vec=basisVector(1, numPostPeriods), Mbar=0, alpha=0.05, hybrid_flag="LF",
        hybrid_kappa=alpha / 10, returnLength=false, biasDirection="positive",
        postPeriodMomentsOnly=true, gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    bh = _asvec(betahat); lv = _asvec(l_vec)
    build = (s, mp) -> (create_A_RMB(numPrePeriods, numPostPeriods; Mbar=Mbar, s=s, max_positive=mp, biasDirection=biasDirection),
                        create_d_RMB(numPrePeriods, numPostPeriods))
    return _rm_family_cs(build, -(numPrePeriods - 1), bh, sigma, numPrePeriods, numPostPeriods, lv,
                         Mbar, alpha, hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly,
                         gridPoints, grid_ub, grid_lb, seed)
end

"""
    computeConditionalCS_DeltaRMM(betahat, sigma, numPrePeriods, numPostPeriods; Mbar=0, monotonicityDirection="increasing", kwargs...)

Confidence set under Δ^RMM(M̄) = the relative-magnitudes restriction Δ^RM(M̄)
intersected with a monotonicity restriction
(`monotonicityDirection = "increasing"` / `"decreasing"`). Other arguments as in
[`computeConditionalCS_DeltaRM`](@ref).
"""
function computeConditionalCS_DeltaRMM(betahat, sigma, numPrePeriods, numPostPeriods;
        l_vec=basisVector(1, numPostPeriods), Mbar=0, alpha=0.05, hybrid_flag="LF",
        hybrid_kappa=alpha / 10, returnLength=false, monotonicityDirection="increasing",
        postPeriodMomentsOnly=true, gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    bh = _asvec(betahat); lv = _asvec(l_vec)
    build = (s, mp) -> (create_A_RMM(numPrePeriods, numPostPeriods; Mbar=Mbar, s=s, max_positive=mp, monotonicityDirection=monotonicityDirection),
                        create_d_RMM(numPrePeriods, numPostPeriods))
    return _rm_family_cs(build, -(numPrePeriods - 1), bh, sigma, numPrePeriods, numPostPeriods, lv,
                         Mbar, alpha, hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly,
                         gridPoints, grid_ub, grid_lb, seed)
end

"""
    computeConditionalCS_DeltaSDRM(betahat, sigma, numPrePeriods, numPostPeriods; Mbar=0, kwargs...)

Confidence set under Δ^SDRM(M̄) - relative magnitudes applied to the
*non-linearity* (second differences) of the differential trend: the largest
post-treatment deviation from a linear trend is at most `Mbar` times the largest
pre-treatment deviation. Requires `numPrePeriods >= 2` (a linear trend needs at
least two pre-periods to be identified). `hybrid_flag` is `"LF"` or `"ARP"`.
Other arguments as in [`computeConditionalCS_DeltaRM`](@ref).
"""
function computeConditionalCS_DeltaSDRM(betahat, sigma, numPrePeriods, numPostPeriods;
        l_vec=basisVector(1, numPostPeriods), Mbar=0, alpha=0.05, hybrid_flag="LF",
        hybrid_kappa=alpha / 10, returnLength=false, postPeriodMomentsOnly=true,
        gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    numPrePeriods == 1 && error("Error: not enough pre-periods (Delta^{SDRM} as base choice)!")
    bh = _asvec(betahat); lv = _asvec(l_vec)
    build = (s, mp) -> (create_A_SDRM(numPrePeriods, numPostPeriods; Mbar=Mbar, s=s, max_positive=mp),
                        create_d_SDRM(numPrePeriods, numPostPeriods))
    return _rm_family_cs(build, -(numPrePeriods - 2), bh, sigma, numPrePeriods, numPostPeriods, lv,
                         Mbar, alpha, hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly,
                         gridPoints, grid_ub, grid_lb, seed)
end

"""
    computeConditionalCS_DeltaSDRMB(betahat, sigma, numPrePeriods, numPostPeriods; Mbar=0, biasDirection="positive", kwargs...)

Confidence set under Δ^SDRMB(M̄) = Δ^SDRM(M̄) intersected with a sign
restriction (`biasDirection = "positive"` / `"negative"`). Requires
`numPrePeriods >= 2`. Other arguments as in [`computeConditionalCS_DeltaSDRM`](@ref).
"""
function computeConditionalCS_DeltaSDRMB(betahat, sigma, numPrePeriods, numPostPeriods;
        l_vec=basisVector(1, numPostPeriods), Mbar=0, alpha=0.05, hybrid_flag="LF",
        hybrid_kappa=alpha / 10, returnLength=false, biasDirection="positive",
        postPeriodMomentsOnly=true, gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    numPrePeriods == 1 && error("Error: not enough pre-periods (Delta^{SDRMB})!")
    bh = _asvec(betahat); lv = _asvec(l_vec)
    build = (s, mp) -> (create_A_SDRMB(numPrePeriods, numPostPeriods; Mbar=Mbar, s=s, max_positive=mp, biasDirection=biasDirection),
                        create_d_SDRMB(numPrePeriods, numPostPeriods))
    return _rm_family_cs(build, -(numPrePeriods - 2), bh, sigma, numPrePeriods, numPostPeriods, lv,
                         Mbar, alpha, hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly,
                         gridPoints, grid_ub, grid_lb, seed)
end

"""
    computeConditionalCS_DeltaSDRMM(betahat, sigma, numPrePeriods, numPostPeriods; Mbar=0, monotonicityDirection="increasing", kwargs...)

Confidence set under Δ^SDRMM(M̄) = Δ^SDRM(M̄) intersected with a monotonicity
restriction (`monotonicityDirection = "increasing"` / `"decreasing"`). Requires
`numPrePeriods >= 2`. Other arguments as in [`computeConditionalCS_DeltaSDRM`](@ref).
"""
function computeConditionalCS_DeltaSDRMM(betahat, sigma, numPrePeriods, numPostPeriods;
        l_vec=basisVector(1, numPostPeriods), Mbar=0, alpha=0.05, hybrid_flag="LF",
        hybrid_kappa=alpha / 10, returnLength=false, monotonicityDirection="increasing",
        postPeriodMomentsOnly=true, gridPoints=10^3, grid_ub=nothing, grid_lb=nothing, seed=0)
    numPrePeriods == 1 && error("Error: not enough pre-periods (Delta^{SDRMM})!")
    bh = _asvec(betahat); lv = _asvec(l_vec)
    build = (s, mp) -> (create_A_SDRMM(numPrePeriods, numPostPeriods; Mbar=Mbar, s=s, max_positive=mp, monotonicityDirection=monotonicityDirection),
                        create_d_SDRMM(numPrePeriods, numPostPeriods))
    return _rm_family_cs(build, -(numPrePeriods - 2), bh, sigma, numPrePeriods, numPostPeriods, lv,
                         Mbar, alpha, hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly,
                         gridPoints, grid_ub, grid_lb, seed)
end
