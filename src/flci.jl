# Fixed-length confidence interval (FLCI) construction, ported from R/flci.R.
# Replaces CVXR + ECOS with JuMP + Clarabel. The worst-case-bias problem is an
# SOCP (linear objective, one convex quadratic constraint); findLowestH is a QP.

# w <-> l space conversions (.wToLFn / .lToWFn).
function wToLFn(w::AbstractVector)
    n = length(w)
    if n == 1
        return [float(w[1])]
    end
    WtoLPreMat = Matrix{Float64}(I, n, n)
    for col in 1:(n-1)
        WtoLPreMat[col+1, col] = -1.0
    end
    return WtoLPreMat * w
end

function lToWFn(l_vec::AbstractVector)
    n = length(l_vec)
    if n == 1
        return [float(l_vec[1])]
    end
    lToWPostMat = Matrix{Float64}(I, n, n)
    for col in 1:(n-1)
        lToWPostMat[col+1, col] = 1.0
    end
    return lToWPostMat * l_vec
end

# .createMatricesForVarianceFromW
function createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)
    lv = _asvec(l_vec)
    pre = 1:numPrePeriods
    post = (numPrePeriods+1):size(sigma, 2)
    SigmaPre = sigma[pre, pre]
    SigmaPrePost = sigma[pre, post]
    SigmaPost = _scalar(lv' * sigma[post, post] * lv)
    WtoLPreMat = Matrix{Float64}(I, numPrePeriods, numPrePeriods)
    if numPrePeriods > 1
        for col in 1:(numPrePeriods-1)
            WtoLPreMat[col+1, col] = -1.0
        end
    end
    UstackWtoLPreMat = hcat(zeros(numPrePeriods, numPrePeriods), WtoLPreMat)
    A_quadratic_sd = UstackWtoLPreMat' * SigmaPre * UstackWtoLPreMat
    A_linear_sd = 2 .* (UstackWtoLPreMat' * SigmaPrePost * lv)
    A_constant_sd = SigmaPost
    return (A_quadratic_sd=Symmetric((A_quadratic_sd .+ A_quadratic_sd') ./ 2),
            A_linear_sd=A_linear_sd, A_constant_sd=A_constant_sd)
end

# Absolute-value linear constraint matrix (.createConstraints_AbsoluteValue).
function _A_absolutevalue(numPrePeriods)
    K = numPrePeriods
    L = [i >= j ? 1.0 : 0.0 for i in 1:K, j in 1:K]   # lower-tri ones incl diag
    return vcat(hcat(-Matrix{Float64}(I, K, K), L),
                hcat(-Matrix{Float64}(I, K, K), -L))
end

# Objective constant for the worst-case bias (.createObjectiveObjectForBias).
function _bias_objective_constant(numPostPeriods, l_vec)
    lv = _asvec(l_vec)
    s_terms = sum(abs(dot(1:s, lv[(numPostPeriods-s+1):numPostPeriods])) for s in 1:numPostPeriods)
    return s_terms - dot(1:numPostPeriods, lv)
end

# .findWorstCaseBiasGivenH
function findWorstCaseBiasGivenH(h, sigma, numPrePeriods, numPostPeriods, l_vec; M=1.0)
    lv = _asvec(l_vec)
    n = 2 * numPrePeriods
    Amat = createMatricesForVarianceFromW(sigma, numPrePeriods, lv)
    A_abs = _A_absolutevalue(numPrePeriods)
    threshold_sum = dot(1:numPostPeriods, lv)
    constant = _bias_objective_constant(numPostPeriods, lv)
    obj_lin = vcat(ones(numPrePeriods), zeros(numPrePeriods))

    model = _conic_model()
    @variable(model, x[1:n])
    @constraint(model, A_abs * x .<= 0)
    @constraint(model, sum(x[numPrePeriods+1:n]) == threshold_sum)
    @constraint(model, x' * Amat.A_quadratic_sd * x + dot(Amat.A_linear_sd, x) + Amat.A_constant_sd <= h^2)
    @objective(model, Min, constant + dot(obj_lin, x))
    optimize!(model)

    if _has_solution(model)
        status = _is_optimal(model) ? "optimal" : "optimal_inaccurate"
        objval = objective_value(model) * M
        xval = JuMP.value.(x)
    else
        status = _status_string(model)
        objval = Inf
        xval = fill(NaN, n)
    end
    optimal_w = xval[(numPrePeriods+1):n]
    optimal_l = any(isnan, optimal_w) ? fill(NaN, numPrePeriods) : wToLFn(optimal_w)
    return (status=status, value=objval, optimal_x=xval, optimal_w=optimal_w, optimal_l=optimal_l)
end

# .findLowestH (minimum-variance affine estimator).
function findLowestH(sigma, numPrePeriods, numPostPeriods, l_vec; sigmascale=10.0, maxscale=10)
    lv = _asvec(l_vec)
    n = 2 * numPrePeriods
    A_abs = _A_absolutevalue(numPrePeriods)
    threshold_sum = dot(1:numPostPeriods, lv)

    function solve_var(scale)
        Amat = createMatricesForVarianceFromW(scale .* sigma, numPrePeriods, lv)
        model = _conic_model()
        @variable(model, x[1:n])
        @constraint(model, A_abs * x .<= 0)
        @constraint(model, sum(x[numPrePeriods+1:n]) == threshold_sum)
        @objective(model, Min, x' * Amat.A_quadratic_sd * x + dot(Amat.A_linear_sd, x) + Amat.A_constant_sd)
        optimize!(model)
        ok = _has_solution(model)
        return (ok ? objective_value(model) : NaN), ok
    end

    val, ok = solve_var(1.0)
    if !ok
        iscale = 0
        while iscale < maxscale && !ok
            iscale += 1
            scaled = iscale > ceil(maxscale / 2) ? sigmascale^(ceil(Int, maxscale / 2) - iscale) : sigmascale^iscale
            v, ok = solve_var(scaled)
            val = v / (sigmascale^iscale)
        end
        ok || @warn "Error in optimization for h0 (tried rescaling)"
    end
    return sqrt(val)
end

# .computeSigmaLFromW
function computeSigmaLFromW(w, sigma, numPrePeriods, numPostPeriods, l_vec)
    Amat = createMatricesForVarianceFromW(sigma, numPrePeriods, l_vec)
    UstackW = vcat(zeros(length(w)), collect(float.(w)))
    return _scalar(UstackW' * Amat.A_quadratic_sd * UstackW + Amat.A_linear_sd' * UstackW + Amat.A_constant_sd)
end

# .findHForMinimumBias
function findHForMinimumBias(sigma, numPrePeriods, numPostPeriods, l_vec)
    lv = _asvec(l_vec)
    w = vcat(zeros(numPrePeriods - 1), dot(1:numPostPeriods, lv))
    hsquared = computeSigmaLFromW(w, sigma, numPrePeriods, numPostPeriods, lv)
    return sqrt(hsquared)
end

# .findOptimalCIDerivativeBisection
function findOptimalCIDerivativeBisection(a, b, M, numPoints, alpha, sigma,
                                          numPrePeriods, numPostPeriods, l_vec)
    f = function (h)
        biasDF = findWorstCaseBiasGivenH(h, sigma, numPrePeriods, numPostPeriods, l_vec)
        maxBias = M * biasDF.value
        if biasDF.value < Inf
            return qfoldednormal(1 - alpha, maxBias / h) * h
        else
            return NaN
        end
    end

    failtol = sqrt(eps(Float64))
    dif = min((b - a) / numPoints, abs(b) * (eps(Float64)^(1 / 3)))
    fa = f(a); fb = f(b)
    fpa = (f(a + dif) - fa) / dif
    fpb = (f(b - dif) - fb) / -dif
    iter = 1
    maxiter = 10 * ceil(Int, log(abs(b - a) / dif) / log(2))
    failed = false
    hstar = NaN

    if (fpa > fpb) || isnan(fa) || isnan(fb)
        failed = true
    elseif fpb < 0
        hstar = b
    elseif fpa > 0
        hstar = a
    else
        while !failed && abs(b - a) > dif
            iter += 1
            x = (a + b) / 2
            fpx = (f(x + dif) - f(x - dif)) / (2 * dif)
            failed = (fpx > fpb + failtol) || (fpx + failtol < fpa) || iter > maxiter
            if fpx > 0
                b = x
            else
                a = x
            end
        end
        hstar = (a + b) / 2
    end
    return failed ? NaN : hstar
end

# .findOptimalFLCI_helper
function findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, l_vec;
                                numPoints=100, alpha=0.05, seed=0)
    lv = _asvec(l_vec)
    h0 = findHForMinimumBias(sigma, numPrePeriods, numPostPeriods, lv)
    hMin = findLowestH(sigma, numPrePeriods, numPostPeriods, lv)
    hstar = findOptimalCIDerivativeBisection(hMin, h0, M, numPoints, alpha, sigma,
                                             numPrePeriods, numPostPeriods, lv)

    local optimal_l, halflength, status, mout
    if isnan(hstar)
        # Grid-search fallback.
        hGrid = range(hMin, h0; length=numPoints)
        best_hl = Inf
        optimal_l = fill(NaN, numPrePeriods); status = "failed"
        for h in hGrid
            bias = findWorstCaseBiasGivenH(h, sigma, numPrePeriods, numPostPeriods, lv)
            (bias.value == Inf) && continue
            (bias.status == "optimal" || bias.status == "optimal_inaccurate") || continue
            maxBias = bias.value * M
            hl = qfoldednormal(1 - alpha, maxBias / h) * h
            if hl < best_hl
                best_hl = hl; optimal_l = bias.optimal_l; status = bias.status
            end
        end
        halflength = best_hl; mout = M
    else
        bias = findWorstCaseBiasGivenH(hstar, sigma, numPrePeriods, numPostPeriods, lv)
        optimal_l = bias.optimal_l; status = bias.status; mout = M
        halflength = qfoldednormal(1 - alpha, (M * bias.value) / hstar) * hstar
    end

    optimalVec = vcat(optimal_l, lv)
    return (optimalVec=optimalVec, optimalPrePeriodVec=optimal_l,
            optimalHalfLength=halflength, M=mout, status=status)
end

"""
    findOptimalFLCI(betahat, sigma, M, numPrePeriods, numPostPeriods; l_vec=basisVector(1, numPostPeriods), numPoints=100, alpha=0.05, seed=0)

Compute the optimal fixed-length confidence interval (FLCI) for
`θ = l_vec' * τ_post` under the smoothness restriction Δ^SD(M). The FLCI is the
affine estimator `optimalVec' * betahat ± optimalHalfLength` whose half-length is
smallest among all affine estimators with valid `1 - alpha` coverage, optimizing
the bias–variance trade-off over the smoothness set.

# Arguments
 * `betahat`, `sigma`, `numPrePeriods`, `numPostPeriods`: the event study.
 * `M`: the smoothness bound.

# Keyword arguments
 * `l_vec = basisVector(1, numPostPeriods)`: the target parameter.
 * `alpha = 0.05`: confidence level is `1 - alpha`.
 * `numPoints = 100`: grid resolution for the internal bias/variance search.
 * `seed`: accepted for API compatibility (the half-length uses an analytic
  folded-normal quantile and is deterministic).

# Returns
A NamedTuple with fields:
 * `FLCI`: the interval `(lb, ub)`.
 * `optimalVec`: the affine weights on `betahat`.
 * `optimalHalfLength`: the half-length.
 * `M`, `status`: the bound used and the solver status.

[`createSensitivityResults`](@ref) with `method = "FLCI"` calls this across a grid
of `M`.
"""
function findOptimalFLCI(betahat, sigma, M, numPrePeriods, numPostPeriods;
                         l_vec=basisVector(1, numPostPeriods), numPoints=100, alpha=0.05, seed=0)
    bh = _asvec(betahat)
    res = findOptimalFLCI_helper(sigma, M, numPrePeriods, numPostPeriods, _asvec(l_vec);
                                 numPoints=numPoints, alpha=alpha, seed=seed)
    pe = dot(res.optimalVec, bh)
    FLCI = (pe - res.optimalHalfLength, pe + res.optimalHalfLength)
    return (FLCI=FLCI, optimalVec=res.optimalVec, optimalHalfLength=res.optimalHalfLength,
            M=res.M, status=res.status)
end
