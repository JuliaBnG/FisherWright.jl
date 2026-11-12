# Hermaphrodite populations (docs/dev/roadmap.md, "Mating schemes").

"""
    Monoecious(; selfing = nothing)

Mating in a population of hermaphrodites: every individual can be either
parent of an offspring, and an offspring of one parent is selfed. Each
offspring draws its first parent uniformly from the cohort, with replacement.
Its second parent is drawn by one of two models:

- `Monoecious()`: uniformly from the whole cohort, the first parent included,
  so offspring are selfed at random, with probability `1/N` (random union of
  gametes; the monoecious Wright–Fisher model of msprime's DTWF).
- `Monoecious(selfing = s)`: the first parent itself with probability `s`,
  and otherwise uniformly from the `N - 1` other individuals, so exactly a
  proportion `s` of the offspring are selfed in expectation, `0 ≤ s ≤ 1`.
  `selfing = 0` excludes selfing.

There are no sexes: the `males` record of a result is 0 for such a generation,
and the pedigree lists the first parent as `sire`, the second as `dam` and the
sex as `missing`. With a fitness model, both parents are drawn with
probability proportional to fitness, the second among the others when the
offspring is not selfed.

With a selfing rate `s` the equilibrium inbreeding coefficient is
`F = s/(2 - s)`, and the effective size is `N/(1 + F) = N(2 - s)/2`
(Pollak 1987; Nordborg and Donnelly 1997); with random selfing it is `N`.
[`effective_size`](@ref) gives the value realized by a run.

# Examples
```julia
using FisherWright

# A population of 50 that selfs half of its offspring: Ne ≈ 37.5
res = fisher_wright(50, 100, [1_000_000], 1.0; mating = Monoecious(selfing = 0.5), result = true)
effective_size(res)
```
"""
struct Monoecious <: AbstractMating
    selfing::Float64    # selfing rate; -1 for random selfing
    function Monoecious(; selfing::Union{Nothing,Real}=nothing)
        selfing === nothing && return new(-1.0)
        0 <= selfing <= 1 ||
            throw(ArgumentError("selfing is a proportion of the offspring and must lie in [0, 1]"))
        return new(Float64(selfing))
    end
end

_monoecious(::AbstractMating) = false
_monoecious(::Monoecious) = true

# Selfing excluded or below one needs a second individual; the census of
# every generation exceeds 1, so this always holds.
_check_mating(::Monoecious, ::AbstractDemography, ::Int, ::Int) = nothing

function Base.show(io::IO, m::Monoecious)
    print(io, "Monoecious(")
    m.selfing >= 0 && print(io, "selfing = ", m.selfing)
    print(io, ")")
end

# No sexes: st.sires and st.dams are left empty, so the record of males is
# 0. Per offspring the draws from `rng` are the first parent, then (for
# 0 < s < 1) a uniform for selfing, then the second parent unless selfed.
_mate!(m::Monoecious, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG) =
    _mono_mate!(st, m.selfing, f, no, rng, 0, 0)

# Kernel of `Monoecious` mating at selfing rate `s` (< 0: random selfing).
# With `level > 0` (AvoidRelatives), an outcrossed second parent related to
# the first at that level is redrawn among the others, up to `max_tries`
# draws in all; whether the offspring is selfed is decided before, as
# without avoidance. Returns the number of offspring of relatives accepted
# after `max_tries`.
function _mono_mate!(st::SimState, s::Float64, f::AbstractFitness, no::Int, rng::AbstractRNG,
    level::Int, max_tries::Int)
    n = st.n
    resize!(st.sex, n)
    empty!(st.sires)
    empty!(st.dams)
    pm, t = st.pm, st.tf
    neutral = f isa Neutral
    if !neutral
        w = st.w
        np = count(i -> w[i] > 0, 1:n)
        np > 0 || error("no parent with positive fitness: the population cannot reproduce")
        np > 1 || s == 1 || s < 0 || error(
            "only one parent has positive fitness, and selfing is below 1: the population cannot reproduce")
        idx = st.fidx
        resize!(idx, n)
        idx .= 1:n
        _alias_build!(t, w, idx)
    end
    forced = 0
    @inbounds for i = 1:no
        p = neutral ? rand(rng, 1:n) : _alias_draw(t, rng)
        q = s < 0 ? (neutral ? rand(rng, 1:n) : _alias_draw(t, rng)) :
            _selfed(s, rng) ? p : _draw_other(p, n, t, neutral, rng)
        if level > 0 && q != p
            k = 1
            while _related(st.rel, p, q, level)
                if k >= max_tries
                    forced += 1
                    break
                end
                q = _draw_other(p, n, t, neutral, rng)
                k += 1
            end
        end
        pm[i, 1] = p
        pm[i, 2] = q
    end
    return forced
end

# One of the individuals other than `p`: uniformly, or in proportion to
# fitness from the alias table `t`, by rejection.
@inline function _draw_other(p::Int, n::Int, t::AliasTable, neutral::Bool, rng)
    neutral && return _other(p, n, rng)
    q = _alias_draw(t, rng)
    while q == p
        q = _alias_draw(t, rng)
    end
    return q
end

# Whether an offspring is selfed at rate `s`; no draw for 0 or 1.
@inline _selfed(s::Float64, rng) = s == 1 || (s > 0 && rand(rng) < s)

# Uniformly one of the `n - 1` individuals other than `p`.
@inline function _other(p::Int, n::Int, rng)
    q = rand(rng, 1:n-1)
    return q >= p ? q + 1 : q
end
