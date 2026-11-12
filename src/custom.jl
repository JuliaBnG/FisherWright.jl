# User-defined mating (docs/dev/roadmap.md, "Mating schemes").

"""
    MatingContext

What a [`CustomMating`](@ref) callback sees of one generation, and where it
writes the matings. Parents are numbered `1:n` within their cohort.

- `generation::Int`: the generation being born, numbered as in the `census`
  record of the result (a run continued through `init` continues it).
- `n::Int`: the number of parents; `noffspring::Int`: the number of offspring.
- `pm`: the `noffspring × 2` mating table to fill: `pm[i, 1]` and `pm[i, 2]` are
  the sire and dam (first and second parent) of offspring `i`.
- `sex::Vector{Bool}`: the sex of each parent, `true` for male, to be set by a
  callback with `sexes = true`; ignored otherwise.
- `fitness`: the relative fitness of each parent under the fitness model of
  the run, as the built-in schemes use it (only ratios within the cohort
  matter; 1 for all under `Neutral`). The callback applies selection itself.
- `ids::UnitRange{Int}`: the pedigree IDs of the parents (see `first_id` of
  [`FisherWrightResult`](@ref)).
- `relatives`: a `6 × n` matrix: rows 1–2 are the parents and rows 3–6 the
  grandparents of each parent, as indices into the cohorts before (rows 1–2
  of a grandparent's children are those of its cohort); 0 is unknown, as for
  the starting population.
- `haplotypes`: the `2n` haplotypes of the parents (individual `i` is
  `2i - 1` and `2i`), neutral mutation positions, read-only.
- `genotypic`, `phenotypic`: the genotypic value and phenotype of each parent
  under [`TraitFitness`](@ref); 0 otherwise.
- `rng`: the random number generator of the run, to be used for every draw so
  that the run stays reproducible.
"""
struct MatingContext{PM,W,REL,H,G,R}
    generation::Int
    n::Int
    noffspring::Int
    pm::PM
    sex::Vector{Bool}
    fitness::W
    ids::UnitRange{Int}
    relatives::REL
    haplotypes::H
    genotypic::G
    phenotypic::G
    rng::R
end

"""
    CustomMating(f; sexes = false)

Mating by a user function: in every generation `f(ctx)` is called with a
[`MatingContext`](@ref) and fills `ctx.pm[i, 1]` and `ctx.pm[i, 2]`, the two
parents of each offspring `i = 1:ctx.noffspring`, as indices `1:ctx.n` into the
parental cohort. Its return value is ignored. Write it with `do` syntax:

```julia
assortative = CustomMating() do ctx
    # pair each offspring's parents among neighbours in phenotype
    order = sortperm(ctx.phenotypic)
    for i in 1:ctx.noffspring
        k = rand(ctx.rng, 1:ctx.n - 1)
        ctx.pm[i, 1], ctx.pm[i, 2] = order[k], order[k + 1]
    end
end
```

With `sexes = false` (the default) the parents are hermaphrodites, as under
[`Monoecious`](@ref): any parent can take either role, an offspring of one
parent in both is selfed, `males` is recorded as 0 and the pedigree sex as
`missing`. With `sexes = true` the callback also sets `ctx.sex[k]` for every
parent `k`, and every sire must be male and every dam female.

The table is checked after each call, and an index outside `1:ctx.n`, or a sire
or dam of the wrong sex, stops the run with an `ArgumentError`. Everything
else is the callback's responsibility: the distribution of family sizes, any
selection (`ctx.fitness` holds the fitness of the parents, which the run does
not otherwise use), and reproducibility, for which every random draw should
come from `ctx.rng`. The callback runs serially, once per generation; it
should not allocate in a long run if speed matters, and should not keep
references to the context, whose buffers are reused.

The record of the result (`paternal_pairs`, `maternal_pairs`, `selfing`) and
[`effective_size`](@ref) describe the realized matings as for the built-in
schemes.
"""
struct CustomMating{F} <: AbstractMating
    f::F
    sexes::Bool
end
CustomMating(f; sexes::Bool=false) = CustomMating{typeof(f)}(f, sexes)

Base.show(io::IO, m::CustomMating) = print(io, "CustomMating(", m.f, "; sexes = ", m.sexes, ")")

_monoecious(m::CustomMating) = !m.sexes
_track_relatives(::CustomMating) = true
_check_mating(::CustomMating, ::AbstractDemography, ::Int, ::Int) = nothing

function _mate!(m::CustomMating, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG)
    n = st.n
    resize!(st.sex, n)
    fill!(st.sex, false)
    w = view(st.w, 1:n)
    f isa Neutral && fill!(w, 1.0)
    ctx = MatingContext(st.ped.gen + 1, n, no, view(st.pm, 1:no, :), st.sex, w,
        st.ped.first:st.ped.first+n-1, view(st.rel, :, 1:n), view(st.prt, 1:2n),
        view(st.gval, 1:n), view(st.pval, 1:n), rng)
    m.f(ctx)
    length(st.sex) == n ||
        throw(ArgumentError("generation $(ctx.generation): the callback resized ctx.sex"))
    pm = st.pm
    @inbounds for i = 1:no, j = 1:2
        1 <= pm[i, j] <= n || throw(ArgumentError(
            "generation $(ctx.generation): parent $(pm[i, j]) of offspring $i is outside 1:$n"))
    end
    empty!(st.sires)
    empty!(st.dams)
    if m.sexes
        @inbounds for i = 1:no
            (st.sex[pm[i, 1]] && !st.sex[pm[i, 2]]) || throw(ArgumentError(
                "generation $(ctx.generation): offspring $i needs a male sire and a female dam"))
        end
        @inbounds for k = 1:n
            st.sex[k] ? push!(st.sires, k) : push!(st.dams, k)
        end
    end
    return nothing
end
