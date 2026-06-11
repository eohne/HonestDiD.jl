# Lightweight, dependency-free result types. They implement the Tables.jl
# interface (so `DataFrame(r)`, `CSV.write(r)`, etc. all work without HonestDiD
# depending on DataFrames) and pretty-print to the REPL.

# Grid confidence set (computeConditionalCS_*)
"""
    ConditionalCS

Result of a `computeConditionalCS_Delta*` call: a test-inversion grid with fields

 * `grid::Vector{Float64}` - candidate values of the target parameter,
 * `accept::Vector{Float64}` - `1` if that value is in the confidence set, else `0`.

The confidence set is `{grid[i] : accept[i] == 1}`; use
[`confidence_interval`](@ref) for its `(lb, ub)`. Implements the Tables.jl
interface, so `DataFrame(cs)` / `Tables.rows(cs)` work if those packages are loaded.
"""
struct ConditionalCS
    grid::Vector{Float64}
    accept::Vector{Float64}
end

Tables.istable(::Type{ConditionalCS}) = true
Tables.columnaccess(::Type{ConditionalCS}) = true
Tables.columns(cs::ConditionalCS) = (grid=cs.grid, accept=cs.accept)
Tables.schema(::ConditionalCS) = Tables.Schema((:grid, :accept), (Float64, Float64))

"""
    confidence_interval(cs::ConditionalCS) -> (lb, ub)

Lower/upper endpoints of the accepted region (`(Inf, -Inf)` if empty).
"""
function confidence_interval(cs::ConditionalCS)
    acc = cs.grid[cs.accept .== 1]
    isempty(acc) ? (Inf, -Inf) : (minimum(acc), maximum(acc))
end

function Base.show(io::IO, ::MIME"text/plain", cs::ConditionalCS)
    lb, ub = confidence_interval(cs)
    n = length(cs.grid)
    if isempty(cs.grid) || lb > ub
        print(io, "ConditionalCS: $n grid points - no accepted values (empty CI)")
    else
        @printf(io, "ConditionalCS: %d grid points on [%.4g, %.4g]; accepted region = [%.6g, %.6g]",
                n, first(cs.grid), last(cs.grid), lb, ub)
    end
end

# Sensitivity table (createSensitivityResults*, constructOriginalCS)
"""
    SensitivityResults

A tidy results table returned by [`createSensitivityResults`](@ref),
[`createSensitivityResults_relativeMagnitudes`](@ref) and
[`constructOriginalCS`](@ref). Columns are accessed as properties - `res.lb`,
`res.ub`, `res.method`, `res.Delta`, and `res.M` or `res.Mbar`.

Pretty-prints as an aligned table and implements the Tables.jl interface, so it
interoperates with the data ecosystem without a DataFrames dependency:
`DataFrame(res)`, `CSV.write("f.csv", res)`, `Tables.rowtable(res)`, etc. all work
once the relevant package is loaded.
"""
struct SensitivityResults
    cols::NamedTuple
end

function Base.getproperty(r::SensitivityResults, s::Symbol)
    s === :cols && return getfield(r, :cols)
    return getproperty(getfield(r, :cols), s)
end
Base.propertynames(r::SensitivityResults) = (:cols, propertynames(getfield(r, :cols))...)

Tables.istable(::Type{SensitivityResults}) = true
Tables.columnaccess(::Type{SensitivityResults}) = true
Tables.columns(r::SensitivityResults) = getfield(r, :cols)

# Build from a vector of NamedTuple rows (column table), without DataFrames.
function _results_from_rows(rows::AbstractVector)
    isempty(rows) && return SensitivityResults(NamedTuple())
    ks = keys(rows[1])
    cols = NamedTuple{ks}(map(k -> [getproperty(r, k) for r in rows], ks))
    return SensitivityResults(cols)
end

_fmtcell(x::Real) = (isinteger(x) && abs(x) < 1e15) ? string(Int(x)) : @sprintf("%.5g", x)
_fmtcell(x::Missing) = ""
_fmtcell(x) = string(x)

function Base.show(io::IO, ::MIME"text/plain", r::SensitivityResults)
    cols = getfield(r, :cols)
    isempty(cols) && return print(io, "SensitivityResults (0 rows)")
    names = collect(keys(cols))
    n = length(cols[names[1]])
    body = [vcat(string(nm), [_fmtcell(cols[nm][i]) for i in 1:n]) for nm in names]
    w = [maximum(length, c) for c in body]
    println(io, "SensitivityResults ($n rows)")
    for row in 1:(n+1)
        print(io, "  ", join((rpad(body[c][row], w[c]) for c in eachindex(names)), "   "))
        println(io)
        row == 1 && println(io, "  ", join((repeat("─", w[c]) for c in eachindex(names)), "───"))
    end
end
