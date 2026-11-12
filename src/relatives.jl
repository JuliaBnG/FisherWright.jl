# Avoidance of matings between relatives (docs/dev/roadmap.md, "Mating
# schemes").

const _LEVELS = (fullsib=1, halfsib=2, cousin=3)

"""
    AvoidRelatives(inner = RandomMating(); level = :halfsib, max_tries = 100)

Mating as by `inner`, a [`RandomMating`](@ref) or [`Monoecious`](@ref) scheme,
without matings between relatives. Each offspring draws its sire (first
parent) as `inner` does; a dam (second parent) related to it is rejected and
redrawn, as `inner` draws it, until an unrelated one is found. Relatives are
judged on the pedigree of the run, at one of three nested levels:

- `:fullsib`: the two share both parents;
- `:halfsib`: they share at least one parent (full and half sibs);
- `:cousin`: they share at least one parent or grandparent (sibs and first
  cousins, full or half).

After `max_tries` draws of the dam without an unrelated one the last pair is
kept: the result counts such offspring in `related_matings`, generation by
generation, so a population too small or too related to avoid relatives
shows in the record rather than stopping the run. With [`Monoecious`](@ref),
whether an offspring is selfed is decided as without avoidance, and only an
outcrossed second parent is redrawn, among the others.

The individuals of the starting population have unknown parents, and are
unrelated to one another; their offspring have unknown grandparents. A run
continued through `init` starts the pedigree anew in the same way.

Under selection the sire is drawn in proportion to fitness and so is every
draw of the dam: the dam of a sire is drawn in proportion to fitness among
the females unrelated to him. Rejection on the dam alone keeps the sires'
usage as under `inner`, while females with many relatives among the males
have fewer offspring.

# Examples
```julia
using FisherWright

# No full- or half-sib matings, 10 sires
m = AvoidRelatives(RandomMating(nmale = 10); level = :halfsib)
res = fisher_wright(50, 100, [1_000_000], 1.0; mating = m, result = true)
sum(res.related_matings)   # offspring of relatives that could not be avoided
```
"""
struct AvoidRelatives{M<:AbstractMating} <: AbstractMating
    inner::M
    level::Int
    max_tries::Int
    function AvoidRelatives(inner::M=RandomMating(); level::Symbol=:halfsib, max_tries::Integer=100) where {M}
        inner isa Union{RandomMating,Monoecious} || throw(ArgumentError(
            "AvoidRelatives wraps RandomMating or Monoecious; FamilySize assigns sires to dams as a whole"))
        haskey(_LEVELS, level) ||
            throw(ArgumentError("level must be :fullsib, :halfsib or :cousin"))
        max_tries > 0 || throw(ArgumentError("max_tries must be positive"))
        return new{M}(inner, _LEVELS[level], Int(max_tries))
    end
end

function Base.show(io::IO, m::AvoidRelatives)
    print(io, "AvoidRelatives(", m.inner, "; level = :", keys(_LEVELS)[m.level], ", max_tries = ", m.max_tries, ")")
end

_monoecious(m::AvoidRelatives) = _monoecious(m.inner)
_check_mating(m::AvoidRelatives, demog::AbstractDemography, n0::Int, nt::Int) =
    _check_mating(m.inner, demog, n0, nt)

# The parents and grandparents of the cohort, as indices into the cohorts
# before it (0 for unknown); only AvoidRelatives keeps them.
_track_relatives(::AbstractMating) = false
_track_relatives(::AvoidRelatives) = true

# Whether individuals `x` and `y` of the current cohort are related at
# `level`. Column `i` of `rel` holds the two parents (rows 1:2) and four
# grandparents (rows 3:6) of individual `i`; 0 is unknown and matches
# nothing. Parents are compared regardless of role, as a hermaphrodite can
# be either.
@inline function _related(rel::Matrix{Int}, x::Int, y::Int, level::Int)
    @inbounds begin
        a, b, c, d = rel[1, x], rel[2, x], rel[1, y], rel[2, y]
        if level == 1
            return a != 0 && b != 0 && ((a == c && b == d) || (a == d && b == c))
        end
        (a != 0 && (a == c || a == d)) && return true
        (b != 0 && (b == c || b == d)) && return true
        level == 2 && return false
        for k = 3:6
            g = rel[k, x]
            g == 0 && continue
            for l = 3:6
                g == rel[l, y] && return true
            end
        end
    end
    return false
end

# Parents and grandparents of the `no` offspring of the current matings,
# from those of the parents; `st.rel` and `st.rnew` are swapped.
function _update_relatives!(st::SimState, no::Int)
    rel, new, pm = st.rel, st.rnew, st.pm
    @inbounds for i = 1:no
        s, d = pm[i, 1], pm[i, 2]
        new[1, i], new[2, i] = s, d
        new[3, i], new[4, i] = rel[1, s], rel[2, s]
        new[5, i], new[6, i] = rel[1, d], rel[2, d]
    end
    st.rel, st.rnew = new, rel
    return nothing
end

_mate!(m::AvoidRelatives{Monoecious}, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG) =
    (st.forced[end] = _mono_mate!(st, m.inner.selfing, f, no, rng, m.level, m.max_tries); nothing)

# Sex as by the inner RandomMating; per offspring the sire, then the dam,
# then (if they are related) further dams. Without relatives the draws are
# those of RandomMating.
function _mate!(m::AvoidRelatives{RandomMating}, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG)
    resize!(st.sex, st.n)
    nm = _nmale(m.inner, st.n)
    if nm == 0
        _draw_sex!(st.sex, st.sires, st.dams, rng, typemax(Int))
    else
        _draw_fixed_sex!(st.sex, st.sires, st.dams, nm, rng)
    end
    sires, dams, pm, w = st.sires, st.dams, st.pm, st.w
    neutral = f isa Neutral
    if !neutral
        (any(i -> w[i] > 0, sires) && any(i -> w[i] > 0, dams)) || error(
            "no parent with positive fitness in one sex: the population cannot reproduce")
        _alias_build!(st.tm, w, sires)
        _alias_build!(st.tf, w, dams)
    end
    forced = 0
    @inbounds for i = 1:no
        s = neutral ? rand(rng, sires) : sires[_alias_draw(st.tm, rng)]
        d = neutral ? rand(rng, dams) : dams[_alias_draw(st.tf, rng)]
        k = 1
        while _related(st.rel, s, d, m.level)
            if k >= m.max_tries
                forced += 1
                break
            end
            d = neutral ? rand(rng, dams) : dams[_alias_draw(st.tf, rng)]
            k += 1
        end
        pm[i, 1] = s
        pm[i, 2] = d
    end
    st.forced[end] = forced
    return nothing
end
