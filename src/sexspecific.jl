# Sex-specific fitness (docs/dev/roadmap.md, "Fitness models"): one fitness
# model for the males and one for the females of every generation. Fitness
# is evaluated for every parent under both models, before the sexes of the
# cohort are drawn, so that the draws of a run keep their order; sires are
# then drawn by the male fitness and dams by the female fitness.

"""
    SexSpecific(male, female)
    SexSpecific(; male = Neutral(), female = Neutral())

Fitness that depends on sex: a parent's fitness is that of the model `male`
if it is male and of `female` if it is female. Each is [`Neutral`](@ref),
[`SiteFitness`](@ref) or [`TraitFitness`](@ref). Pass it as the `fitness`
keyword of [`fisher_wright`](@ref).

Both models act on the same selected alleles, so they must describe the same
loci:

- two `SiteFitness` models need tables with the same positions and mutation
  rate; the selection coefficients and dominance may differ, for example
  `s > 0` in males and `s < 0` in females for sexually antagonistic
  selection;
- two `TraitFitness` models need the same `TraitLoci` table and the same
  environmental (co)variance, since every individual has one phenotype; the
  fitness functions may differ, for example `Stabilizing` around different
  optima. Traits expressed differently in the two sexes are two traits of
  one table (pleiotropy), the male fitness depending on one and the female
  fitness on the other;
- with `Neutral` on one side, selection acts in the other sex only.

Selection is soft and within each sex, as for the other models: every
offspring draws its sire among the males in proportion to the male fitness
and its dam among the females in proportion to the female fitness. The
mating scheme must draw separate sexes: [`RandomMating`](@ref),
[`FamilySize`](@ref) (with splits that can follow fitness) or
[`AvoidRelatives`](@ref) of `RandomMating`; [`Monoecious`](@ref) and
[`CustomMating`](@ref) are refused.

Fitness is evaluated for every parent under both models before the sexes of
the generation are drawn, so a run draws its random numbers in the same order
as a run with one model (each individual has one phenotype under a trait
model). The fitness record of the result (`mean_fitness`,
`mean_log_fitness`, `var_log_fitness`) is that of each parent under the model
of its own sex; the trait record is that of the trait model.

# Examples
```julia
using FisherWright

chr = [1_000_000, 1_000_000]
pos = fitness_loci(chr, 50; s = 0.0).pos
antagonistic = SexSpecific(SiteFitness(FitnessLoci(pos, 0.05; μ = 1e-5)),
                           SiteFitness(FitnessLoci(pos, -0.05; μ = 1e-5)))
res = fisher_wright(200, 100, chr, 1.0; fitness = antagonistic, result = true)
```
"""
struct SexSpecific{M<:AbstractFitness,F<:AbstractFitness} <: AbstractFitness
    male::M
    female::F
    function SexSpecific(male::AbstractFitness, female::AbstractFitness)
        for f in (male, female)
            f isa Union{Neutral,SiteFitness,TraitFitness} || throw(ArgumentError(
                "each sex takes Neutral, SiteFitness or TraitFitness; got $(typeof(f))"))
        end
        _same_loci(male, female)
        return new{typeof(male),typeof(female)}(male, female)
    end
end
SexSpecific(; male::AbstractFitness=Neutral(), female::AbstractFitness=Neutral()) = SexSpecific(male, female)

Base.show(io::IO, f::SexSpecific) = print(io, "SexSpecific(male = ", f.male, ", female = ", f.female, ")")

# The models of the two sexes must share the selected alleles.
_same_loci(::AbstractFitness, ::AbstractFitness) = nothing
_same_loci(::SiteFitness, ::TraitFitness) =
    throw(ArgumentError("SiteFitness and TraitFitness cannot be combined: they keep different loci"))
_same_loci(m::TraitFitness, f::SiteFitness) = _same_loci(f, m)
function _same_loci(m::SiteFitness, f::SiteFitness)
    (m.loci.pos == f.loci.pos && m.loci.μ == f.loci.μ) || throw(ArgumentError(
        "the fitness loci of the two sexes must have the same positions and mutation rate"))
    return nothing
end
function _same_loci(m::TraitFitness, f::TraitFitness)
    a, b = m.loci, f.loci
    (a.pos == b.pos && a.a == b.a && a.h == b.h && a.μ == b.μ) || throw(ArgumentError(
        "the two sexes must share one TraitLoci table: every individual has one genotypic value"))
    m.Σe == f.Σe || throw(ArgumentError(
        "the two sexes must have the same environmental (co)variance: every individual has one phenotype"))
    return nothing
end

# The model that carries the selected alleles: the male one if it has them.
_primary(f::SexSpecific) = _selected(f.male) ? f.male : f.female
_trait(f::SexSpecific) = f.male isa TraitFitness ? f.male : f.female isa TraitFitness ? f.female : nothing

_selected(f::SexSpecific) = _selected(_primary(f))
_sel_poisson(f::SexSpecific) = _sel_poisson(_primary(f))
_sel_mutation!(f::SexSpecific, st::SimState, p_sel::Poisson{Float64}) = _sel_mutation!(_primary(f), st, p_sel)
_check_fitness(f::SexSpecific, tbp::UInt32) = (_check_fitness(f.male, tbp); _check_fitness(f.female, tbp))
_initial_selected(init::FisherWrightResult, f::SexSpecific, nh::Int) = _initial_selected(init, _primary(f), nh)
_model_traits(f::SexSpecific) = (t = _trait(f); t === nothing ? nothing : ntraits(t))

# Absolute log fitness of the n parents under one model, into `lw`.
_sex_log_w!(lw, ::Neutral, st::SimState) = fill!(lw, 0.0)
_sex_log_w!(lw, f::TraitFitness, st::SimState) = _trait_log_w!(lw, f.w, st.fit.pval, st.n)
function _sex_log_w!(lw, f::SiteFitness, st::SimState)
    sel = st.fit.sel
    offset = 0.0
    for x in st.fit.ssubs
        offset += f.lhom[searchsortedfirst(f.loci.pos, x)]
    end
    Threads.@threads for i = 1:st.n
        @inbounds lw[i] = _log_fitness(f, sel[2i-1], sel[2i]) + offset
    end
    return lw
end

# Relative fitness, scaled so that the fittest has 1; all 0 if none is fit
# (mating then stops the run if that sex has no fit member).
function _scale_w!(w::Vector{Float64}, lw, n::Int)
    top = maximum(view(lw, 1:n))
    @inbounds for i = 1:n
        w[i] = top > -Inf ? exp(lw[i] - top) : 0.0
    end
    return w
end

function _fitness!(f::SexSpecific, st::SimState, rng::AbstractRNG)
    n = st.n
    # Sires are drawn by st.fit.w and dams by st.fit.wf: two vectors from here on.
    st.fit.wf === st.fit.w && (st.fit.wf = similar(st.fit.w))
    t = _trait(f)
    t === nothing || _trait_values!(t, st, rng)
    lm, lf = view(st.fit.lwraw, 1, 1:n), view(st.fit.lwraw, 2, 1:n)
    _sex_log_w!(lm, f.male, st)
    _sex_log_w!(lf, f.female, st)
    _scale_w!(st.fit.w, lm, n)
    _scale_w!(st.fit.wf, lf, n)
    return nothing
end

# The fitness record waits for the sexes of the parents.
_record_fitness!(::SexSpecific, ::SimState) = nothing
_record_trait!(f::SexSpecific, st::SimState) = _trait(f) === nothing ? _record_no_trait!(st) : nothing

# After mating: the fitness of each parent under the model of its sex.
_record_mated!(::AbstractFitness, ::SimState) = nothing
function _record_mated!(::SexSpecific, st::SimState)
    n, sex, lw = st.n, st.mate.sex, st.fit.lwraw
    m = s = s2 = 0.0
    @inbounds for i = 1:n
        x = lw[sex[i] ? 1 : 2, i]
        m += x
        s += exp(x)
    end
    m /= n
    @inbounds for i = 1:n
        s2 += (lw[sex[i] ? 1 : 2, i] - m)^2
    end
    push!(st.rec.mean_log_fitness, m)
    push!(st.rec.var_log_fitness, s2 / (n - 1))
    push!(st.rec.mean_fitness, s / n)
    return nothing
end

# Sexes must be drawn by the run.
_check_sexes(::AbstractMating, ::AbstractFitness) = nothing
function _check_sexes(m::AbstractMating, ::SexSpecific)
    (m isa CustomMating || _monoecious(m)) && throw(ArgumentError(
        "SexSpecific needs separate sexes drawn by the run (RandomMating, FamilySize or " *
        "AvoidRelatives of RandomMating), not $(m)"))
    return nothing
end
