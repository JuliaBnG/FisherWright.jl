# Weighted sampling with replacement by Walker's alias method, built with
# Vose's algorithm (Vose 1991, IEEE Trans. Softw. Eng. 17: 972-975): O(n) to
# build, O(1) per draw, and no allocation once the table has reached the size
# of the largest set it is built for. `Distributions.sampler(Categorical(p))`
# allocates on every call and is therefore not used in the generation loop.

"""
    AliasTable(nmax = 0)

Preallocated alias table for weighted draws from up to `nmax` outcomes
without reallocation. Build it with `_alias_build!` and draw with
`_alias_draw`.

Outcome `k` is returned with probability
`(prob[k] + Σ_{j : alias[j] = k} (1 - prob[j])) / n`, which equals its
normalised weight up to rounding.
"""
struct AliasTable
    prob::Vector{Float64}
    alias::Vector{Int}
    small::Vector{Int}
    large::Vector{Int}
end

function AliasTable(nmax::Integer=0)
    v() = sizehint!(Int[], nmax)
    return AliasTable(sizehint!(Float64[], nmax), v(), v(), v())
end

Base.length(t::AliasTable) = length(t.prob)

"""
    _alias_build!(t::AliasTable, w, idx) -> AliasTable

Build `t` for `length(idx)` outcomes, outcome `k` having weight `w[idx[k]]`.
Weights must be finite and non-negative with a positive sum; an outcome of
weight zero is never drawn.
"""
function _alias_build!(t::AliasTable, w::AbstractVector{Float64}, idx::AbstractVector{Int})
    n = length(idx)
    n > 0 || throw(ArgumentError("no outcome to draw from"))
    total = 0.0
    @inbounds for i in idx
        x = w[i]
        (isfinite(x) && x >= 0) ||
            throw(ArgumentError("weights must be finite and non-negative; got $x"))
        total += x
    end
    (total > 0 && isfinite(total)) ||
        throw(ArgumentError("weights must have a positive, finite sum; got $total"))
    prob, alias, small, large = t.prob, t.alias, t.small, t.large
    resize!(prob, n)
    resize!(alias, n)
    empty!(small)
    empty!(large)
    scale = n / total
    @inbounds for k = 1:n
        p = w[idx[k]] * scale
        prob[k] = p
        alias[k] = k
        push!(p < 1 ? small : large, k)
    end
    # Pair a deficient outcome with one in excess, which donates the rest of
    # the cell; the donor is re-filed by what it has left.
    @inbounds while !isempty(small) && !isempty(large)
        s = pop!(small)
        l = pop!(large)
        alias[s] = l
        p = (prob[l] + prob[s]) - 1
        prob[l] = p
        push!(p < 1 ? small : large, l)
    end
    # Whatever is left has probability 1 up to rounding.
    @inbounds for k in small
        prob[k] = 1.0
    end
    @inbounds for k in large
        prob[k] = 1.0
    end
    empty!(small)
    empty!(large)
    return t
end

"""
    _alias_draw(t::AliasTable, rng) -> Int

One outcome `k` in `1:length(t)`: a uniform cell, then the cell's own outcome
or its alias. Two draws from `rng`.
"""
@inline function _alias_draw(t::AliasTable, rng)
    k = rand(rng, 1:length(t.prob))
    return @inbounds rand(rng) < t.prob[k] ? k : t.alias[k]
end

# Probability of each outcome implied by the table; for tests.
function _alias_probabilities(t::AliasTable)
    n = length(t)
    p = copy(t.prob)
    for j = 1:n
        t.alias[j] == j || (p[t.alias[j]] += 1 - t.prob[j])
    end
    return p ./ n
end
