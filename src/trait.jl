# Selection on a quantitative trait (docs/dev/roadmap.md, "Fitness
# models"): genotypic values at an explicit table of trait loci, a phenotype
# with environmental noise, and fitness as a function of the phenotype. The
# derived alleles at the trait loci live in the selected stream `st.sel`, as
# those of `SiteFitness`.

"""
    TraitLoci(pos, a; h = 0.5, μ = 0.0)

Table of the loci of a quantitative trait, fixed for a run.

- `pos`: genome-wide positions, strictly increasing.
- `a`: effect of the derived allele at each site (a number for all sites, or
  one per site): the genotypic value of the derived homozygote. Any finite
  value.
- `h`: dominance at each site (a number or one per site): the heterozygote
  has value `h * a`. With `h = 0.5` the site is additive.
- `μ`: mutation rate to the derived allele per site per haplotype per
  generation; `0` means the trait alleles of a run are those of its starting
  population.

The genotypic value of an individual is the sum over sites of `0`, `h*a` and
`a` for zero, one and two copies of the derived allele, so that an individual
without derived alleles has value 0. Sites are biallelic and there is no
back-mutation.

# Fields
`pos::Vector{UInt32}`, `a::Vector{Float64}`, `h::Vector{Float64}`, `μ::Float64`.

# See also
[`trait_loci`](@ref), [`TraitFitness`](@ref)
"""
struct TraitLoci
    pos::Vector{UInt32}
    a::Vector{Float64}
    h::Vector{Float64}
    μ::Float64
    function TraitLoci(pos::Vector{UInt32}, a::Vector{Float64}, h::Vector{Float64}, μ::Float64)
        K = length(pos)
        (length(a) == K && length(h) == K) ||
            throw(ArgumentError("pos, a and h must have the same length"))
        (isempty(pos) || pos[1] >= 1) || throw(ArgumentError("positions are 1-based"))
        all(k -> pos[k-1] < pos[k], 2:K) ||
            throw(ArgumentError("pos must be sorted and unique"))
        (all(isfinite, a) && all(isfinite, h)) || throw(ArgumentError("every a and h must be finite"))
        (isfinite(μ) && μ >= 0) || throw(ArgumentError("μ must be finite and non-negative"))
        return new(pos, a, h, μ)
    end
end

TraitLoci(pos::AbstractVector{<:Integer}, a; h=0.5, μ::Real=0.0) = TraitLoci(
    Vector{UInt32}(pos), _per_site(a, length(pos)), _per_site(h, length(pos)), Float64(μ))

Base.length(l::TraitLoci) = length(l.pos)

"""
    trait_loci(chr, nsites; a, h = 0.5, μ = 0.0, rng = Random.default_rng()) -> TraitLoci

Draw `nsites` distinct trait loci uniformly over a genome with chromosome
lengths `chr` and build their [`TraitLoci`](@ref) table.

- `a`: a number, used for every site, or a univariate distribution from
  Distributions.jl, from which one effect per site is drawn, for example
  `Normal(0, 0.1)`.
- `h`: a number, or a function `a -> h` applied to each effect.
- `μ`: mutation rate per site per haplotype per generation.

# Examples
```julia
using FisherWright, Distributions, Random

loci = trait_loci([1_000_000, 1_000_000], 500; a = Normal(0, 0.1), μ = 1e-5, rng = Xoshiro(1))
```
"""
function trait_loci(
    chr::AbstractVector{<:Integer},
    nsites::Integer;
    a,
    h=0.5,
    μ::Real=0.0,
    rng::AbstractRNG=Random.default_rng(),
)
    a isa Real || a isa UnivariateDistribution ||
        throw(ArgumentError("a must be a number or a univariate distribution"))
    # The sites are drawn as for fitness_loci, from the same rng.
    sites = fitness_loci(chr, nsites; s=0.0, rng=rng).pos
    av = a isa Real ? fill(Float64(a), Int(nsites)) : Float64.(rand(rng, a, Int(nsites)))
    hv = h isa Real ? fill(Float64(h), Int(nsites)) : Float64[h(x) for x in av]
    return TraitLoci(sites, av, hv, Float64(μ))
end

"""
    Stabilizing(θ, ω²)

Gaussian stabilizing selection on the phenotype `P`: fitness
`exp(-(P - θ)^2 / (2ω²))`, with optimum `θ` and width `ω² > 0`. For use in
[`TraitFitness`](@ref).
"""
struct Stabilizing
    θ::Float64
    ω²::Float64
    function Stabilizing(θ::Real, ω²::Real)
        (isfinite(θ) && isfinite(ω²) && ω² > 0) ||
            throw(ArgumentError("θ must be finite and ω² finite and positive"))
        return new(θ, ω²)
    end
end

"""
    Directional(β)

Directional selection on the phenotype `P`: fitness `exp(β P)`. For a normally
distributed phenotype with variance `σ²_P`, the selection differential is
`β σ²_P`. For use in [`TraitFitness`](@ref).
"""
struct Directional
    β::Float64
    Directional(β::Real) = isfinite(β) ? new(β) : throw(ArgumentError("β must be finite"))
end

# Log fitness of phenotype P.
_log_w(f::Stabilizing, P::Float64) = -(P - f.θ)^2 / (2f.ω²)
_log_w(f::Directional, P::Float64) = f.β * P
function _log_w(f, P::Float64)
    w = f(P)
    (w isa Real && w >= 0) || error("the fitness function returned $(repr(w)) for P = $P; a non-negative number is required")
    return log(Float64(w))
end

"""
    TraitFitness(loci::TraitLoci; σe = 0.0, fitness)

Selection on a quantitative trait. Each generation, every adult has the
genotypic value `G` given by `loci` (the sum over sites, including sites whose
derived allele was fixed and moved to the substitutions), the phenotype
`P = G + e` with `e ~ Normal(0, σe²)` drawn anew for every individual and
generation, and the fitness `fitness(P)`. `fitness` is [`Stabilizing`](@ref),
[`Directional`](@ref), or any function that returns a non-negative number.
Pass it as the `fitness` keyword of [`fisher_wright`](@ref).

Selection is soft and acts within each sex, as for [`SiteFitness`](@ref):
every offspring draws its sire among the males and its dam among the females
with probability proportional to fitness, with replacement.

The environmental deviations are drawn serially from the `rng` of the run,
before mating, so a run does not depend on the number of threads.

# Examples
```julia
using FisherWright, Distributions

chr = [1_000_000, 1_000_000]
loci = trait_loci(chr, 200; a = Normal(0, 0.2), μ = 1e-4)
f = TraitFitness(loci; σe = 1.0, fitness = Stabilizing(0.0, 5.0))
res = fisher_wright(100, 500, chr, 1.0; fitness = f, result = true)
res.mean_genotypic[end], res.var_genotypic[end]
```
"""
struct TraitFitness{F} <: AbstractFitness
    loci::TraitLoci
    σe::Float64
    w::F
    function TraitFitness(loci::TraitLoci; σe::Real=0.0, fitness)
        (isfinite(σe) && σe >= 0) || throw(ArgumentError("σe must be finite and non-negative"))
        return new{typeof(fitness)}(loci, Float64(σe), fitness)
    end
end

_selected(::TraitFitness) = true
_sel_poisson(f::TraitFitness) =
    f.loci.μ > 0 && length(f.loci) > 0 ? Poisson(length(f.loci) * f.loci.μ) : nothing
_sel_mutation!(f::TraitFitness, st::SimState, p_sel::Poisson{Float64}) = _sel_mutation!(f.loci.pos, st, p_sel)

# Genotypic value of an individual with selected haplotypes `x` and `y`,
# without the substituted sites: a merge-walk as in `_log_fitness`.
function _genotypic_value(loci::TraitLoci, x::Vector{UInt32}, y::Vector{UInt32})
    pos, a, h = loci.pos, loci.a, loci.h
    K, lx, ly = length(pos), length(x), length(y)
    g = 0.0
    i = j = k = 1
    @inbounds while i <= lx && j <= ly
        u, v = x[i], y[j]
        k = searchsortedfirst(pos, min(u, v), k, K, Base.Order.Forward)
        if u == v
            g += a[k]
            i += 1
            j += 1
        else
            g += h[k] * a[k]
            u < v ? (i += 1) : (j += 1)
        end
    end
    @inbounds while i <= lx
        k = searchsortedfirst(pos, x[i], k, K, Base.Order.Forward)
        g += h[k] * a[k]
        i += 1
    end
    @inbounds while j <= ly
        k = searchsortedfirst(pos, y[j], k, K, Base.Order.Forward)
        g += h[k] * a[k]
        j += 1
    end
    return g
end

# Genotypic value, phenotype and relative fitness of every parent, and the
# records of the generation. The value of the substituted sites, carried by
# every parent on both haplotypes, is added to G.
function _fitness!(f::TraitFitness, st::SimState, rng::AbstractRNG)
    n, sel, G, P, w = st.n, st.sel, st.gval, st.pval, st.w
    loci = f.loci
    offset = 0.0
    for x in st.ssubs
        offset += loci.a[searchsortedfirst(loci.pos, x)]
    end
    Threads.@threads for i = 1:n
        @inbounds G[i] = _genotypic_value(loci, sel[2i-1], sel[2i]) + offset
    end
    # Environmental deviations serially from rng: independent of threads.
    σe = f.σe
    @inbounds for i = 1:n
        P[i] = σe > 0 ? G[i] + σe * randn(rng) : G[i]
    end
    @inbounds for i = 1:n
        w[i] = _log_w(f.w, P[i])
    end
    g, p, lw = view(G, 1:n), view(P, 1:n), view(w, 1:n)
    push!(st.gmean, mean(g))
    push!(st.gvar, var(g))
    push!(st.pmean, mean(p))
    push!(st.pvar, var(p))
    push!(st.flmean, mean(lw))
    push!(st.flvar, var(lw))
    push!(st.fmean, sum(exp, lw) / n)
    top = maximum(lw)
    top > -Inf || error("every individual has zero fitness: the population cannot reproduce")
    @inbounds for i = 1:n
        w[i] = exp(w[i] - top)
    end
    return nothing
end

_record_fitness!(::TraitFitness, ::SimState) = nothing
_record_trait!(::TraitFitness, ::SimState) = nothing

# Errors that the first generation would otherwise hit midway.
function _check_fitness(f::TraitFitness, tbp::UInt32)
    (isempty(f.loci.pos) || f.loci.pos[end] <= tbp) ||
        throw(ArgumentError("a trait locus lies beyond the genome length $(Int(tbp))"))
    return nothing
end

"""
    to_haplotype(result::FisherWrightResult, loci::TraitLoci; include_fixed::Bool = true) -> Tuple{Haplotype, DataFrame}

Genotypes at the trait loci of a run with [`TraitFitness`](@ref), with their
true effects, as `to_haplotype(result, loci::FitnessLoci)`: one row per site of
the table, in its order, and a table with columns `chr`, `pos`, `a`, `h` and
`frq`. The genotypic value of individual `i` is
`sum(k -> (d = hap[k, 2i-1] + hap[k, 2i]; d == 0 ? 0.0 : d == 1 ? h[k] * a[k] : a[k]), 1:K)`
with `include_fixed = true`.
"""
function to_haplotype(result::FisherWrightResult, loci::TraitLoci; include_fixed::Bool=true)
    xy, cbp = _selected_matrix(result, loci, include_fixed)
    map = DataFrame(
        chr=Int8[searchsortedfirst(cbp, p) for p in loci.pos],
        pos=copy(loci.pos),
        a=copy(loci.a),
        h=copy(loci.h),
        frq=vec(sum(xy; dims=2)) ./ size(xy, 2),
    )
    return Haplotype(xy), map
end
