# Selection on quantitative traits (docs/dev/roadmap.md, "Fitness models"):
# genotypic values at an explicit table of trait loci, a phenotype with
# environmental noise, and fitness as a function of the phenotype. With
# several traits the loci are shared (pleiotropy): each site has an effect
# on every trait. The derived alleles at the trait loci live in the
# selected stream `st.fit.sel`, as those of `SiteFitness`.

"""
    TraitLoci(pos, a; h = 0.5, μ = 0.0)

Table of the loci of one or more quantitative traits, fixed for a run. With
`T > 1` traits every site has an effect on each of them (pleiotropy).

- `pos`: genome-wide positions, strictly increasing.
- `a`: effect of the derived allele, the genotypic value of the derived
  homozygote: a number (one trait, the same at every site), a vector (one
  trait, one effect per site), or a `K × T` matrix (`K` sites, `T` traits;
  row `k` is the effect of site `k` on each trait). Any finite value; a zero
  means the site does not affect that trait.
- `h`: dominance, so that the heterozygote has value `h * a`: a number, a
  vector with one value per site (the same for every trait), or a `K × T`
  matrix. With `h = 0.5` the site is additive.
- `μ`: mutation rate to the derived allele per site per haplotype per
  generation; `0` means the trait alleles of a run are those of its starting
  population.

The genotypic value of an individual for trait `t` is the sum over sites of
`0`, `h[k, t] * a[k, t]` and `a[k, t]` for zero, one and two copies of the
derived allele, so that an individual without derived alleles has value 0.
Sites are biallelic and there is no back-mutation.

# Fields
`pos::Vector{UInt32}`, `a::Matrix{Float64}` (`K × T`), `h::Matrix{Float64}`
(`K × T`), `μ::Float64`. The number of traits is [`ntraits`](@ref).

# See also
[`trait_loci`](@ref), [`TraitFitness`](@ref)
"""
struct TraitLoci
    pos::Vector{UInt32}
    a::Matrix{Float64}
    h::Matrix{Float64}
    μ::Float64
    at::Matrix{Float64}     # T × K: the effects by site, contiguous per site
    hat::Matrix{Float64}    # T × K: h .* a, the heterozygote values
    function TraitLoci(pos::Vector{UInt32}, a::Matrix{Float64}, h::Matrix{Float64}, μ::Float64)
        K = length(pos)
        size(a, 2) >= 1 || throw(ArgumentError("a must have at least one column (trait)"))
        (size(a, 1) == K && size(h) == size(a)) ||
            throw(ArgumentError("pos, a and h must have the same number of sites, and a and h the same traits"))
        (isempty(pos) || pos[1] >= 1) || throw(ArgumentError("positions are 1-based"))
        all(k -> pos[k-1] < pos[k], 2:K) ||
            throw(ArgumentError("pos must be sorted and unique"))
        (all(isfinite, a) && all(isfinite, h)) || throw(ArgumentError("every a and h must be finite"))
        (isfinite(μ) && μ >= 0) || throw(ArgumentError("μ must be finite and non-negative"))
        return new(pos, a, h, μ, permutedims(a), permutedims(h .* a))
    end
end

_effects(a::Real, K::Int) = fill(Float64(a), K, 1)
_effects(a::AbstractVector{<:Real}, ::Int) = reshape(Vector{Float64}(a), :, 1)
_effects(a::AbstractMatrix{<:Real}, ::Int) = Matrix{Float64}(a)
_effects(a, ::Int) = throw(ArgumentError("a must be a number, a vector or a matrix"))

_dominance(h::Real, K::Int, T::Int) = fill(Float64(h), K, T)
_dominance(h::AbstractVector{<:Real}, ::Int, T::Int) = repeat(Vector{Float64}(h), 1, T)
_dominance(h::AbstractMatrix{<:Real}, ::Int, ::Int) = Matrix{Float64}(h)
_dominance(h, ::Int, ::Int) = throw(ArgumentError("h must be a number, a vector or a matrix"))

function TraitLoci(pos::AbstractVector{<:Integer}, a; h=0.5, μ::Real=0.0)
    av = _effects(a, length(pos))
    return TraitLoci(Vector{UInt32}(pos), av, _dominance(h, length(pos), size(av, 2)), Float64(μ))
end

Base.length(l::TraitLoci) = length(l.pos)

"""
    ntraits(x) -> Int

The number of traits of a [`TraitLoci`](@ref) table, a [`TraitFitness`](@ref)
model, or the trait record of a [`FisherWrightResult`](@ref) (1 for a result
without a trait model).
"""
ntraits(l::TraitLoci) = size(l.a, 2)

"""
    trait_loci(chr, nsites; a, h = 0.5, μ = 0.0, rng = Random.default_rng()) -> TraitLoci

Draw `nsites` distinct trait loci uniformly over a genome with chromosome
lengths `chr` and build their [`TraitLoci`](@ref) table.

- `a`: a number, used for every site (one trait); a univariate distribution
  from Distributions.jl, from which one effect per site is drawn, for example
  `Normal(0, 0.1)` (one trait); or a multivariate distribution, from which
  one vector of effects per site is drawn, for `T = length(a)` traits with
  pleiotropic effects, for example `MvNormal(zeros(2), [1.0 0.5; 0.5 1.0])`,
  whose covariance matrix sets the correlation of the effects.
- `h`: a number, or a function `a -> h` applied to each effect (of each
  trait).
- `μ`: mutation rate per site per haplotype per generation.

The draws come from `rng` only: the sites first, then the effects.

# Examples
```julia
using FisherWright, Distributions, Random

loci = trait_loci([1_000_000, 1_000_000], 500; a = Normal(0, 0.1), μ = 1e-5, rng = Xoshiro(1))
two = trait_loci([1_000_000, 1_000_000], 500; a = MvNormal(zeros(2), [0.01 -0.005; -0.005 0.01]),
                 rng = Xoshiro(1))
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
    a isa Real || a isa UnivariateDistribution || a isa MultivariateDistribution ||
        throw(ArgumentError("a must be a number, a univariate or a multivariate distribution"))
    # The sites are drawn as for fitness_loci, from the same rng.
    sites = fitness_loci(chr, nsites; s=0.0, rng=rng).pos
    K = Int(nsites)
    av = a isa Real ? fill(Float64(a), K, 1) :
         a isa UnivariateDistribution ? reshape(Float64.(rand(rng, a, K)), :, 1) :
         Matrix{Float64}(permutedims(rand(rng, a, K)))
    hv = h isa Real ? fill(Float64(h), size(av)) : Float64[h(x) for x in av]
    return TraitLoci(sites, av, hv, Float64(μ))
end

"""
    Stabilizing(θ, ω²)
    Stabilizing(θ::AbstractVector, Ω::AbstractMatrix)
    Stabilizing(θ::AbstractVector, ω²::Real)

Gaussian stabilizing selection on the phenotype `P`. For one trait, fitness
`exp(-(P - θ)^2 / (2ω²))`, with optimum `θ` and width `ω² > 0`. For `T`
traits, fitness `exp(-(P - θ)' Ω⁻¹ (P - θ) / 2)` with the optimum `θ`
(length `T`) and a symmetric positive-definite `T × T` matrix `Ω`, whose
off-diagonal elements make selection correlational; a number `ω²` stands for
`ω² I`. For use in [`TraitFitness`](@ref).

# Fields
`θ::Vector{Float64}`, `Ω::Matrix{Float64}` (`[ω²;;]` for one trait).
"""
struct Stabilizing
    θ::Vector{Float64}
    Ω::Matrix{Float64}
    Ωi::Matrix{Float64}
    function Stabilizing(θ::AbstractVector{<:Real}, Ω::AbstractMatrix{<:Real})
        T = length(θ)
        (T >= 1 && size(Ω) == (T, T)) ||
            throw(ArgumentError("θ must have length T ≥ 1 and Ω size T × T"))
        (all(isfinite, θ) && all(isfinite, Ω)) || throw(ArgumentError("θ and Ω must be finite"))
        issymmetric(Ω) || throw(ArgumentError("Ω must be symmetric"))
        C = cholesky(Symmetric(Matrix{Float64}(Ω)); check=false)
        issuccess(C) || throw(ArgumentError("Ω must be positive definite"))
        return new(Vector{Float64}(θ), Matrix{Float64}(Ω), inv(C))
    end
end

function Stabilizing(θ::Real, ω²::Real)
    (isfinite(θ) && isfinite(ω²) && ω² > 0) ||
        throw(ArgumentError("θ must be finite and ω² finite and positive"))
    return Stabilizing([θ], fill(ω², 1, 1))
end
Stabilizing(θ::AbstractVector{<:Real}, ω²::Real) = Stabilizing(θ, Matrix{Float64}(ω² * I(length(θ))))

"""
    Directional(β)

Directional selection on the phenotype `P`: fitness `exp(β P)` for one trait,
or `exp(β'P)` for `T` traits with the selection gradient `β` (length `T`;
a zero element leaves that trait unselected directly). For a multivariate
normal phenotype with covariance matrix `P`, the selection differential is
`P β` (`β σ²_P` for one trait). For use in [`TraitFitness`](@ref).

# Fields
`β::Vector{Float64}`.
"""
struct Directional
    β::Vector{Float64}
    function Directional(β::AbstractVector{<:Real})
        (!isempty(β) && all(isfinite, β)) || throw(ArgumentError("β must be finite and non-empty"))
        return new(Vector{Float64}(β))
    end
end
Directional(β::Real) = Directional([β])

# The number of traits a fitness function of P is for; 0 for a user function.
_fitness_traits(f::Stabilizing) = length(f.θ)
_fitness_traits(f::Directional) = length(f.β)
_fitness_traits(_) = 0

# Log fitness of phenotype P: a number for one trait, a vector otherwise.
# One trait keeps the scalar formulae, so that its runs do not change.
_log_w(f::Stabilizing, P::Float64) = -(P - f.θ[1])^2 / (2f.Ω[1, 1])
_log_w(f::Directional, P::Float64) = f.β[1] * P
function _log_w(f::Stabilizing, P::AbstractVector{Float64})
    θ, Ωi, T = f.θ, f.Ωi, length(f.θ)
    q = 0.0
    @inbounds for j = 1:T, i = 1:T
        q += (P[i] - θ[i]) * Ωi[i, j] * (P[j] - θ[j])
    end
    return -q / 2
end
_log_w(f::Directional, P::AbstractVector{Float64}) = dot(f.β, P)
function _log_w(f, P)
    w = f(P)
    (w isa Real && w >= 0) || error("the fitness function returned $(repr(w)) for P = $P; a non-negative number is required")
    return log(Float64(w))
end

"""
    TraitFitness(loci::TraitLoci; σe = 0.0, Σe = nothing, fitness)

Selection on one or more quantitative traits. Each generation, every adult
has the genotypic values `G` given by `loci` (one per trait, the sum over
sites, including sites whose derived allele was fixed and moved to the
substitutions), the phenotype `P = G + e`, and the fitness `fitness(P)`.

The environmental deviations `e` are drawn anew for every individual and
generation: `e ~ Normal(0, σe²)` for one trait; for `T` traits,
`e ~ MvNormal(0, Σe)`, given either as `σe` (a number or a vector of `T`
standard deviations, independent across traits) or as the covariance matrix
`Σe` (symmetric positive semi-definite, `T × T`), not both. The environmental
covariance adds to the genetic covariance that the shared loci induce.

`fitness` is [`Stabilizing`](@ref), [`Directional`](@ref) (for as many traits
as `loci` has), or any function that returns a non-negative number; such a
function is called with the phenotype, a `Float64` for one trait and an
`AbstractVector{Float64}` of length `T` (not to be kept, its memory is
reused) otherwise. Pass the model as the `fitness` keyword of
[`fisher_wright`](@ref).

Selection is soft and acts within each sex, as for [`SiteFitness`](@ref):
every offspring draws its sire among the males and its dam among the females
with probability proportional to fitness, with replacement.

The environmental deviations are drawn serially from the `rng` of the run,
before mating, so a run does not depend on the number of threads. With `T`
traits each individual takes `T` standard normal draws, in trait order,
transformed by a fixed factor `L` of `Σe` (`L L' = Σe`).

# Examples
```julia
using FisherWright, Distributions

chr = [1_000_000, 1_000_000]
loci = trait_loci(chr, 200; a = Normal(0, 0.2), μ = 1e-4)
f = TraitFitness(loci; σe = 1.0, fitness = Stabilizing(0.0, 5.0))
res = fisher_wright(100, 500, chr, 1.0; fitness = f, result = true)
res.records.mean_genotypic[1, end], res.records.var_genotypic[1, 1, end]

# Two traits with negatively correlated effects; selection on the first.
two = trait_loci(chr, 200; a = MvNormal(zeros(2), [0.04 -0.02; -0.02 0.04]), μ = 1e-4)
f2 = TraitFitness(two; σe = [1.0, 1.0], fitness = Directional([0.5, 0.0]))
res2 = fisher_wright(100, 100, chr, 1.0; fitness = f2, result = true)
res2.records.mean_genotypic[:, end]          # the second trait responds by correlation
```
"""
struct TraitFitness{F} <: AbstractFitness
    loci::TraitLoci
    Σe::Matrix{Float64}
    L::Matrix{Float64}
    w::F
    function TraitFitness(loci::TraitLoci; σe=0.0, Σe=nothing, fitness)
        T = ntraits(loci)
        Σ, L = Σe === nothing ? _from_sd(σe, T) : _from_cov(σe, Σe, T)
        k = _fitness_traits(fitness)
        k == 0 || k == T ||
            throw(ArgumentError("the fitness model is for $k traits; the loci have $T"))
        return new{typeof(fitness)}(loci, Σ, L, fitness)
    end
end

function _from_sd(σe, T::Int)
    σ = σe isa Real ? fill(Float64(σe), T) :
        σe isa AbstractVector{<:Real} ? Vector{Float64}(σe) :
        throw(ArgumentError("σe must be a number or a vector"))
    length(σ) == T || throw(ArgumentError("σe must have one standard deviation per trait ($T)"))
    all(x -> isfinite(x) && x >= 0, σ) || throw(ArgumentError("σe must be finite and non-negative"))
    L = zeros(T, T)
    for t = 1:T
        L[t, t] = σ[t]
    end
    return L .^ 2, L
end

function _from_cov(σe, Σe, T::Int)
    σe isa Real && iszero(σe) || throw(ArgumentError("give the environmental (co)variance as σe or Σe, not both"))
    Σe isa AbstractMatrix{<:Real} && size(Σe) == (T, T) ||
        throw(ArgumentError("Σe must be a $T × $T matrix"))
    all(isfinite, Σe) && issymmetric(Σe) || throw(ArgumentError("Σe must be finite and symmetric"))
    Σ = Matrix{Float64}(Σe)
    C = cholesky(Symmetric(Σ); check=false)
    issuccess(C) && return Σ, Matrix(C.L)
    # Positive semi-definite: a square root from the eigendecomposition.
    λ, V = eigen(Symmetric(Σ))
    λ[1] >= -1e-10 * max(1.0, abs(λ[end])) || throw(ArgumentError("Σe must be positive semi-definite"))
    return Σ, V * Diagonal(sqrt.(max.(λ, 0.0)))
end

ntraits(f::TraitFitness) = ntraits(f.loci)
_model_traits(f::TraitFitness) = ntraits(f)

_selected(::TraitFitness) = true
_sel_poisson(f::TraitFitness) =
    f.loci.μ > 0 && length(f.loci) > 0 ? Poisson(length(f.loci) * f.loci.μ) : nothing
_sel_mutation!(f::TraitFitness, st::SimState, p_sel::Poisson{Float64}) = _sel_mutation!(f.loci.pos, st, p_sel)

# Genotypic values of an individual with selected haplotypes `x` and `y`,
# without the substituted sites, written to `g` (one per trait): a
# merge-walk as in `_log_fitness`.
function _genotypic_value!(g::AbstractVector{Float64}, loci::TraitLoci, x::Vector{UInt32}, y::Vector{UInt32})
    pos, at, hat = loci.pos, loci.at, loci.hat
    T = size(at, 1)
    K, lx, ly = length(pos), length(x), length(y)
    fill!(g, 0.0)
    i = j = k = 1
    @inbounds while i <= lx && j <= ly
        u, v = x[i], y[j]
        k = searchsortedfirst(pos, min(u, v), k, K, Base.Order.Forward)
        if u == v
            for t = 1:T
                g[t] += at[t, k]
            end
            i += 1
            j += 1
        else
            for t = 1:T
                g[t] += hat[t, k]
            end
            u < v ? (i += 1) : (j += 1)
        end
    end
    @inbounds while i <= lx
        k = searchsortedfirst(pos, x[i], k, K, Base.Order.Forward)
        for t = 1:T
            g[t] += hat[t, k]
        end
        i += 1
    end
    @inbounds while j <= ly
        k = searchsortedfirst(pos, y[j], k, K, Base.Order.Forward)
        for t = 1:T
            g[t] += hat[t, k]
        end
        j += 1
    end
    return g
end

_genotypic_value(loci::TraitLoci, x::Vector{UInt32}, y::Vector{UInt32}) =
    _genotypic_value!(zeros(ntraits(loci)), loci, x, y)

# Mean and (sample) covariance matrix of the columns of X (T × n), appended
# to the records. One trait uses mean and var of a contiguous vector, as
# before pleiotropy, so that its records do not change.
function _push_moments!(m::Vector{Float64}, v::Vector{Float64}, X::Matrix{Float64}, n::Int)
    T = size(X, 1)
    if T == 1
        x = view(vec(X), 1:n)
        push!(m, mean(x))
        push!(v, var(x))
    else
        Xn = view(X, :, 1:n)
        append!(m, vec(mean(Xn; dims=2)))
        append!(v, vec(cov(Xn; dims=2)))
    end
    return nothing
end

# Genotypic values and phenotypes of every parent, in `st.fit.gval` and
# `st.fit.pval`, and their records. The value of the substituted sites, carried
# by every parent on both haplotypes, is added to G.
function _trait_values!(f::TraitFitness, st::SimState, rng::AbstractRNG)
    n, sel, G, P = st.n, st.fit.sel, st.fit.gval, st.fit.pval
    loci = f.loci
    T = ntraits(loci)
    size(G, 1) == T || error("the state holds $(size(G, 1)) traits; the model has $T")
    offset = zeros(T)
    for x in st.fit.ssubs
        k = searchsortedfirst(loci.pos, x)
        for t = 1:T
            offset[t] += loci.at[t, k]
        end
    end
    Threads.@threads for i = 1:n
        g = view(G, :, i)
        _genotypic_value!(g, loci, sel[2i-1], sel[2i])
        @inbounds for t = 1:T
            g[t] += offset[t]
        end
    end
    # Environmental deviations serially from rng: independent of threads.
    L = f.L
    if T == 1
        σe = L[1, 1]
        @inbounds for i = 1:n
            P[1, i] = σe > 0 ? G[1, i] + σe * randn(rng) : G[1, i]
        end
    else
        copyto!(view(P, :, 1:n), view(G, :, 1:n))
        if any(!iszero, L)
            z = zeros(T)
            @inbounds for i = 1:n
                for t = 1:T
                    z[t] = randn(rng)
                end
                for t = 1:T, s = 1:T
                    P[t, i] += L[t, s] * z[s]
                end
            end
        end
    end
    _push_moments!(st.rec.mean_genotypic, st.rec.var_genotypic, G, n)
    _push_moments!(st.rec.mean_phenotypic, st.rec.var_phenotypic, P, n)
    return nothing
end

# Log fitness `fit(P)` of the phenotypes of the n parents, into `lw`.
function _trait_log_w!(lw::AbstractVector{Float64}, fit, P::Matrix{Float64}, n::Int)
    if size(P, 1) == 1
        @inbounds for i = 1:n
            lw[i] = _log_w(fit, P[1, i])
        end
    else
        @inbounds for i = 1:n
            lw[i] = _log_w(fit, view(P, :, i))
        end
    end
    return lw
end

# Genotypic values, phenotypes and relative fitness of every parent, and the
# records of the generation.
function _fitness!(f::TraitFitness, st::SimState, rng::AbstractRNG)
    n, w = st.n, st.fit.w
    _trait_values!(f, st, rng)
    _trait_log_w!(w, f.w, st.fit.pval, n)
    lw = view(w, 1:n)
    push!(st.rec.mean_log_fitness, mean(lw))
    push!(st.rec.var_log_fitness, var(lw))
    push!(st.rec.mean_fitness, sum(exp, lw) / n)
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
the table, in its order, and a table with columns `chr`, `pos`, the effects and
dominance, and `frq`. With one trait the effects and dominance are the columns
`a` and `h`; with `T` traits they are `a1, …, aT` and `h1, …, hT`. The genotypic
value of individual `i` for trait `t` is
`sum(k -> (d = hap[k, 2i-1] + hap[k, 2i]; d == 0 ? 0.0 : d == 1 ? h[k, t] * a[k, t] : a[k, t]), 1:K)`
with `include_fixed = true`, where `a = loci.a` and `h = loci.h`.
"""
function to_haplotype(result::FisherWrightResult, loci::TraitLoci; include_fixed::Bool=true)
    xy, cbp = _selected_matrix(result, loci, include_fixed)
    T = ntraits(loci)
    map = DataFrame(chr=Int8[searchsortedfirst(cbp, p) for p in loci.pos], pos=copy(loci.pos))
    if T == 1
        map.a = loci.a[:, 1]
        map.h = loci.h[:, 1]
    else
        for t = 1:T
            map[!, "a$t"] = loci.a[:, t]
        end
        for t = 1:T
            map[!, "h$t"] = loci.h[:, t]
        end
    end
    map.frq = vec(sum(xy; dims=2)) ./ size(xy, 2)
    return Haplotype(xy), map
end
