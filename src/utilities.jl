# Utility functions ported from R/utilities.R, R/arp-nuisance.R (helpers),
# and R/flci.R (.qfoldednormal). These are internal helpers plus the exported
# `basisVector`.

# Coerce a possibly-matrix / possibly-scalar quantity to a plain Float64 scalar.
_scalar(x::Number) = float(x)
_scalar(x::AbstractArray) = float(only(x))

# Coerce betahat / l_vec inputs (vector or n×1 / 1×n matrix) to Vector{Float64}.
_asvec(x::AbstractVector) = collect(float.(x))
function _asvec(x::AbstractMatrix)
    @assert minimum(size(x)) == 1 "expected a vector but got a $(size(x,1))×$(size(x,2)) matrix"
    return collect(float.(vec(x)))
end

"""
    basisVector(index=1, size=1)

Return a length-`size` `Vector{Float64}` of zeros with a `1` at position `index`.
Mirrors R `basisVector`, which returns a column matrix; here we return a vector.
"""
function basisVector(index::Integer=1, size::Integer=1)
    v = zeros(Float64, size)
    v[index] = 1.0
    return v
end

# .selectionMat from R/utilities.R
# `selection` may be an integer or a vector of indices.
function selectionMat(selection, size::Integer; select::AbstractString="columns")
    sel = selection isa Integer ? [selection] : collect(selection)
    if select == "rows"
        m = zeros(Float64, length(sel), size)
        for (i, s) in enumerate(sel)
            m[i, s] = 1.0
        end
    else
        m = zeros(Float64, size, length(sel))
        for (i, s) in enumerate(sel)
            m[s, i] = 1.0
        end
    end
    return m
end

# .LeeCFN: c = Sigma*eta / (eta' Sigma eta)
function LeeCFN(eta::AbstractVecOrMat, Sigma::AbstractMatrix)
    e = vec(eta)
    return (Sigma * e) ./ _scalar(e' * Sigma * e)
end

# .VLoVUpFN: truncation limits for the conditional/Lee approach.
function VLoVUpFN(eta::AbstractVecOrMat, Sigma::AbstractMatrix,
                  A::AbstractMatrix, b::AbstractVector, z::AbstractVector)
    c = LeeCFN(eta, Sigma)
    Ac = A * c
    objective = (b .- A * z) ./ Ac
    negidx = findall(<(0.0), Ac)
    posidx = findall(>(0.0), Ac)
    VLo = isempty(negidx) ? -Inf : maximum(objective[negidx])
    VUp = isempty(posidx) ? Inf : minimum(objective[posidx])
    return (VLo, VUp)
end

# .warnIfNotSymmPSD
function warnIfNotSymmPSD(sigma::AbstractMatrix)
    asym = maximum(abs.(sigma .- sigma'))
    if asym > 0
        @warn "matrix sigma not exactly symmetric (largest asymmetry was $(asym))"
    end
    lambda = eigvals(Symmetric((sigma .+ sigma') ./ 2))
    if any(lambda .< 0)
        @warn "matrix sigma not numerically positive semi-definite (smallest eigenvalue was $(minimum(lambda)))"
    end
    return nothing
end

# .stopIfNotConformable
function stopIfNotConformable(betahat, sigma, numPrePeriods, numPostPeriods, l_vec)
    bvec = _asvec(betahat)
    betaL = length(bvec)
    sigmaR, sigmaC = size(sigma, 1), size(sigma, 2)
    sigmaR == sigmaC || error("expected a square matrix but sigma was $(sigmaR) by $(sigmaC)")
    sigmaR == betaL || error("betahat (length $(betaL)) and sigma ($(sigmaR) by $(sigmaC)) were non-conformable")
    numPeriods = numPrePeriods + numPostPeriods
    numPeriods == betaL || error("betahat (length $(betaL)) and pre + post periods ($(numPrePeriods) + $(numPostPeriods)) were non-conformable")
    length(_asvec(l_vec)) == numPostPeriods || error("l_vec (length $(length(_asvec(l_vec)))) and post periods ($(numPostPeriods)) were non-conformable")
    return nothing
end

# Truncated normal quantile (.norminvp / .norminvp_generalized)
# Quantile at probability `p` of a standard normal truncated to [l, u].
function norminvp(p::Real, l::Real, u::Real)
    if l == u
        return float(l)
    end
    d = truncated(Normal(0.0, 1.0), float(l), float(u))
    return quantile(d, float(p))
end

# .norminvp_generalized
function norminvp_generalized(p, l, u; mu=0.0, sd=1.0)
    μ = _scalar(mu); σ = _scalar(sd)
    ln = (_scalar(l) - μ) / σ
    un = (_scalar(u) - μ) / σ
    qn = norminvp(_scalar(p), ln, un)
    return μ + qn * σ
end

# Folded-normal quantile (.qfoldednormal)
# R does this by simulation. We invert the folded-normal CDF directly,
# F(x) = Phi((x-mu)/sd) - Phi((-x-mu)/sd) for x >= 0, by bisection - exact and
# deterministic. Vectorized over `mu`.
function _qfoldednormal_scalar(p::Real, mu::Real, sd::Real)
    p <= 0 && return 0.0
    Z = Normal(0.0, 1.0)
    F(x) = cdf(Z, (x - mu) / sd) - cdf(Z, (-x - mu) / sd)
    # Bracket the root: grow the upper bound until F(hi) >= p.
    lo = 0.0
    hi = abs(mu) + sd * 8.0 + 1.0
    iter = 0
    while F(hi) < p && iter < 200
        hi *= 2.0
        iter += 1
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        if F(mid) < p
            lo = mid
        else
            hi = mid
        end
        (hi - lo) < 1e-12 && break
    end
    return (lo + hi) / 2
end

qfoldednormal(p, mu::Real; sd::Real=1.0) = _qfoldednormal_scalar(p, mu, sd)
qfoldednormal(p, mu::AbstractVector; sd::Real=1.0) = [_qfoldednormal_scalar(p, m, sd) for m in mu]

# Reduced row echelon form + .construct_Gamma
function rref(A::AbstractMatrix; tol::Real=1e-10)
    M = Matrix{Float64}(A)
    nrows, ncols = size(M)
    lead = 1
    for r in 1:nrows
        lead > ncols && break
        i = r
        while abs(M[i, lead]) < tol
            i += 1
            if i > nrows
                i = r
                lead += 1
                lead > ncols && return M
            end
        end
        M[r, :], M[i, :] = M[i, :], M[r, :]
        M[r, :] ./= M[r, lead]
        for j in 1:nrows
            if j != r
                M[j, :] .-= M[j, lead] .* M[r, :]
            end
        end
        lead += 1
    end
    return M
end

leading_one(r::Integer, B::AbstractMatrix; tol::Real=1e-10) =
    findfirst(j -> abs(B[r, j]) > tol, 1:size(B, 2))

# .construct_Gamma: invertible matrix whose first row is l'.
function construct_Gamma(l::AbstractVector)
    barT = length(l)
    B = hcat(collect(float.(l)), Matrix{Float64}(I, barT, barT))
    rrefB = rref(B)
    leading_ones = [leading_one(r, rrefB) for r in 1:size(rrefB, 1)]
    Gamma = permutedims(B[:, leading_ones])
    abs(det(Gamma)) < 1e-12 && error("Something went wrong in RREF algorithm.")
    return Gamma
end
