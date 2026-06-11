# Builders for the (A, d) that define each Delta restriction set, i.e.
# {delta : A*delta <= d}. Ported from the R deltasd/deltarm/deltasdrm files and
# their sign/monotonicity variants.
#
# A few R quirks are deliberately kept so the matrices match: `a:b` counts down
# when a > b (so `1:0` is `c(1, 0)`), and assigning to row index 0 is a no-op.
# The `_rcolon` helper plus in-range guards reproduce that.

# R-style colon: counts down (guarded) when a > b
_rcolon(a::Integer, b::Integer) = a <= b ? collect(a:b) : collect(a:-1:b)

# Drop column `j` from a matrix.
_dropcol(A::AbstractMatrix, j::Integer) = A[:, setdiff(1:size(A, 2), j)]

# pracma::repmat(v, n, 1) for a 1×k row `v`: stack n identical rows.
_stackrows(v::AbstractMatrix, n::Integer) = n <= 0 ? zeros(0, size(v, 2)) : repeat(v, n, 1)

# Remove all-zero rows (||row||^2 <= 1e-10), matching R's zerorows filter.
function _dropzerorows(A::AbstractMatrix)
    keep = [sum(abs2, @view A[r, :]) > 1e-10 for r in 1:size(A, 1)]
    return A[keep, :]
end

# Delta^SD
function create_A_SD(numPrePeriods::Integer, numPostPeriods::Integer;
                     postPeriodMomentsOnly::Bool=false)
    P = numPrePeriods + numPostPeriods
    Atilde = zeros(P - 1, P + 1)
    for r in _rcolon(1, P - 1)
        (1 <= r <= P - 1) || continue
        Atilde[r, r:(r+2)] = [1.0, -2.0, 1.0]
    end
    Atilde = _dropcol(Atilde, numPrePeriods + 1)
    if postPeriodMomentsOnly
        postPeriodIndices = (numPrePeriods + 1):size(Atilde, 2)
        prePeriodOnlyRows = findall(r -> sum(Atilde[r, postPeriodIndices] .!= 0) == 0, 1:size(Atilde, 1))
        # R's `Atilde[-prePeriodOnlyRows, ]` keeps no rows when the drop set is
        # empty (the `-integer(0)` thing), so match that
        Atilde = isempty(prePeriodOnlyRows) ? Atilde[1:0, :] :
                 Atilde[setdiff(1:size(Atilde, 1), prePeriodOnlyRows), :]
    end
    return vcat(Atilde, -Atilde)
end

function create_d_SD(numPrePeriods::Integer, numPostPeriods::Integer, M::Real;
                     postPeriodMomentsOnly::Bool=false)
    A = create_A_SD(numPrePeriods, numPostPeriods; postPeriodMomentsOnly=postPeriodMomentsOnly)
    return fill(float(M), size(A, 1))
end

# monotonicity matrix A_M
function create_A_M(numPrePeriods::Integer, numPostPeriods::Integer,
                    monotonicityDirection::AbstractString, postPeriodMomentsOnly::Bool)
    P = numPrePeriods + numPostPeriods
    A_M = zeros(P, P)
    # On purpose NOT matching R here. R's `for (r in 1:(numPrePeriods-1))` becomes
    # `1:0 == c(1,0)` when numPrePeriods==1 and writes a stray A_M[1,2] = -1, which
    # weakens the restriction from delta_pre <= 0 to delta_pre <= delta_post1. The
    # native range is empty for numPre==1 and gives the intended delta_pre <= 0.
    # Identical to R for numPrePeriods >= 2.
    for r in 1:(numPrePeriods-1)
        A_M[r, r:(r+1)] = [1.0, -1.0]
    end
    A_M[numPrePeriods, numPrePeriods] = 1.0
    if numPostPeriods > 0
        A_M[numPrePeriods+1, numPrePeriods+1] = -1.0
        if numPostPeriods > 1
            for r in (numPrePeriods+2):(numPrePeriods+numPostPeriods)
                A_M[r, (r-1):r] = [1.0, -1.0]
            end
        end
    end
    if postPeriodMomentsOnly
        postPeriodIndices = (numPrePeriods + 1):size(A_M, 2)
        prePeriodOnlyRows = findall(r -> sum(A_M[r, postPeriodIndices] .!= 0) == 0, 1:size(A_M, 1))
        A_M = isempty(prePeriodOnlyRows) ? A_M[1:0, :] :
              A_M[setdiff(1:size(A_M, 1), prePeriodOnlyRows), :]
    end
    if monotonicityDirection == "decreasing"
        A_M = -A_M
    elseif monotonicityDirection != "increasing"
        error("direction must be 'increasing' or 'decreasing'")
    end
    return A_M
end

# Sign restriction (A_B)
function create_A_B(numPrePeriods::Integer, numPostPeriods::Integer, biasDirection::AbstractString)
    P = numPrePeriods + numPostPeriods
    A_B = -Matrix{Float64}(I, P, P)
    A_B = A_B[(numPrePeriods+1):(numPrePeriods+numPostPeriods), :]
    if biasDirection == "negative"
        A_B = -A_B
    elseif biasDirection != "positive"
        error("Input biasDirection must equal either `positive' or `negative'")
    end
    return A_B
end

# Delta^SDB
function create_A_SDB(numPrePeriods, numPostPeriods; biasDirection="positive", postPeriodMomentsOnly=false)
    A_SD = create_A_SD(numPrePeriods, numPostPeriods; postPeriodMomentsOnly=postPeriodMomentsOnly)
    A_B = create_A_B(numPrePeriods, numPostPeriods, biasDirection)
    return vcat(A_SD, A_B)
end

function create_d_SDB(numPrePeriods, numPostPeriods, M; postPeriodMomentsOnly=false)
    d_SD = create_d_SD(numPrePeriods, numPostPeriods, M; postPeriodMomentsOnly=postPeriodMomentsOnly)
    return vcat(d_SD, zeros(numPostPeriods))
end

# Delta^SDM
function create_A_SDM(numPrePeriods, numPostPeriods; monotonicityDirection="increasing", postPeriodMomentsOnly=false)
    A_M = create_A_M(numPrePeriods, numPostPeriods, monotonicityDirection, postPeriodMomentsOnly)
    A_SD = create_A_SD(numPrePeriods, numPostPeriods; postPeriodMomentsOnly=postPeriodMomentsOnly)
    return vcat(A_SD, A_M)
end

function create_d_SDM(numPrePeriods, numPostPeriods, M; postPeriodMomentsOnly=false)
    d_SD = create_d_SD(numPrePeriods, numPostPeriods, M; postPeriodMomentsOnly=postPeriodMomentsOnly)
    d_M = zeros(postPeriodMomentsOnly ? numPostPeriods : numPrePeriods + numPostPeriods)
    return vcat(d_SD, d_M)
end

# Delta^RM
# Builds the base RM matrix (first-difference relative-magnitude) for a given
# maximal period `s` and sign. `dropZero` removes the t=0 reference column.
function _A_RM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    P = numPrePeriods + numPostPeriods
    Atilde = zeros(P, P + 1)
    for r in 1:P
        Atilde[r, r:(r+1)] = [-1.0, 1.0]
    end
    v_max_dif = zeros(1, P + 1)
    v_max_dif[1, (numPrePeriods+s):(numPrePeriods+1+s)] = [-1.0, 1.0]
    if !max_positive
        v_max_dif = -v_max_dif
    end
    A_UB = vcat(_stackrows(v_max_dif, numPrePeriods), _stackrows(Mbar .* v_max_dif, numPostPeriods))
    A = vcat(Atilde .- A_UB, -Atilde .- A_UB)
    A = _dropzerorows(A)
    return A
end

function create_A_RM(numPrePeriods, numPostPeriods; Mbar=1, s=0, max_positive=true, dropZero=true)
    A = _A_RM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    return dropZero ? _dropcol(A, numPrePeriods + 1) : A
end

function create_d_RM(numPrePeriods, numPostPeriods; dropZero=true)
    A = create_A_RM(numPrePeriods, numPostPeriods; Mbar=0, s=0, dropZero=dropZero)
    return zeros(size(A, 1))
end

# Delta^RMB
function create_A_RMB(numPrePeriods, numPostPeriods; Mbar=1, s=0, max_positive=true, dropZero=true, biasDirection)
    A = _A_RM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    A_B = create_A_B(numPrePeriods, numPostPeriods, biasDirection)
    A = dropZero ? _dropcol(A, numPrePeriods + 1) : A
    return vcat(A, A_B)
end

function create_d_RMB(numPrePeriods, numPostPeriods; dropZero=true)
    A_RM = create_A_RM(numPrePeriods, numPostPeriods; Mbar=0, s=0, dropZero=dropZero)
    return zeros(size(A_RM, 1) + numPostPeriods)
end

# Delta^RMM
function create_A_RMM(numPrePeriods, numPostPeriods; Mbar=1, s=0, max_positive=true, dropZero=true, monotonicityDirection)
    A = _A_RM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    A_M = create_A_M(numPrePeriods, numPostPeriods, monotonicityDirection, false)
    A = dropZero ? _dropcol(A, numPrePeriods + 1) : A
    return vcat(A, A_M)
end

function create_d_RMM(numPrePeriods, numPostPeriods; dropZero=true)
    A_RM = create_A_RM(numPrePeriods, numPostPeriods; Mbar=0, s=0, dropZero=dropZero)
    return zeros(size(A_RM, 1) + numPrePeriods + numPostPeriods)
end

# Delta^SDRM
function _A_SDRM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    P = numPrePeriods + numPostPeriods
    Atilde = zeros(P - 1, P + 1)
    for r in 1:(P-1)
        Atilde[r, r:(r+2)] = [1.0, -2.0, 1.0]
    end
    v_max_dif = zeros(1, P + 1)
    v_max_dif[1, (numPrePeriods+1+s-2):(numPrePeriods+1+s)] = [1.0, -2.0, 1.0]
    if !max_positive
        v_max_dif = -v_max_dif
    end
    A_UB = vcat(_stackrows(v_max_dif, numPrePeriods - 1), _stackrows(Mbar .* v_max_dif, numPostPeriods))
    A = vcat(Atilde .- A_UB, -Atilde .- A_UB)
    A = _dropzerorows(A)
    return A
end

function create_A_SDRM(numPrePeriods, numPostPeriods; Mbar=1, s=0, max_positive=true, dropZero=true)
    A = _A_SDRM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    return dropZero ? _dropcol(A, numPrePeriods + 1) : A
end

function create_d_SDRM(numPrePeriods, numPostPeriods; dropZero=true)
    A = create_A_SDRM(numPrePeriods, numPostPeriods; Mbar=0, s=0, dropZero=dropZero)
    return zeros(size(A, 1))
end

# Delta^SDRMB
function create_A_SDRMB(numPrePeriods, numPostPeriods; Mbar=1, s=0, max_positive=true, dropZero=true, biasDirection)
    A = _A_SDRM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    A_B = create_A_B(numPrePeriods, numPostPeriods, biasDirection)
    A = dropZero ? _dropcol(A, numPrePeriods + 1) : A
    return vcat(A, A_B)
end

function create_d_SDRMB(numPrePeriods, numPostPeriods; dropZero=true)
    A_SDRM = create_A_SDRM(numPrePeriods, numPostPeriods; Mbar=0, s=0, dropZero=dropZero)
    return zeros(size(A_SDRM, 1) + numPostPeriods)
end

# Delta^SDRMM
function create_A_SDRMM(numPrePeriods, numPostPeriods; Mbar=1, s=0, max_positive=true, dropZero=true, monotonicityDirection)
    A = _A_SDRM_base(numPrePeriods, numPostPeriods, Mbar, s, max_positive)
    A_M = create_A_M(numPrePeriods, numPostPeriods, monotonicityDirection, false)
    A = dropZero ? _dropcol(A, numPrePeriods + 1) : A
    return vcat(A, A_M)
end

function create_d_SDRMM(numPrePeriods, numPostPeriods; dropZero=true)
    A_SDRM = create_A_SDRM(numPrePeriods, numPostPeriods; Mbar=0, s=0, dropZero=dropZero)
    return zeros(size(A_SDRM, 1) + numPrePeriods + numPostPeriods)
end
