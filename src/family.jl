# Mating with controlled family sizes (docs/dev/roadmap.md, "Mating
# schemes"). The census size of a generation is fixed (constant or from
# `popsize`); its offspring are split among the dams, and every dam is
# mated to one sire, so her offspring are full sibs.

# How a total is split among parents: `:equal` (as equal as possible),
# `:random` (uniform multinomial) or a Gamma shape `k` (Dirichlet-multinomial
# with weights Gamma(k, 1)).
function _split_mode(x, what)
    x === :equal && return (:equal, 0.0)
    x === :random && return (:random, 0.0)
    x isa Real && isfinite(x) && x > 0 && return (:gamma, Float64(x))
    throw(ArgumentError("$what must be :equal, :random or a positive number (a Gamma shape); got $(repr(x))"))
end

"""
    FamilySize(; dams = :random, sires = :random, nsires = nothing, ndams = nothing,
               nmale = nothing, sexratio = nothing)

Mating with controlled family sizes. In every generation the `N_g` offspring
are split among the dams, and every dam is mated to one sire, so all her
offspring are full sibs; the dams are split among the sires. A sire's family
is the sum of the families of his dams. Sex is assigned to the parents as by
[`RandomMating`](@ref) with the same `nmale` and `sexratio`.

`dams` sets how the offspring are split among the dams, and `sires` how the
dams are split among the sires:

- `:equal`: as equal as possible; when the split is not exact, the parents
  that receive one more are drawn at random.
- `:random`: a uniform multinomial split, as each offspring choosing its dam
  (or each dam its sire) at random. Family sizes are then binomial, close to
  Poisson.
- a positive number `k`: a multinomial split with weights drawn from
  Gamma(`k`, 1) for every parent in every generation (Dirichlet-multinomial
  counts, close to negative binomial). Smaller `k` means more variable
  families; as `k → ∞` it approaches `:random`.

With `n_d` offspring per dam among `N_f` dams, two offspring share their dam
with probability `P_f = E[Σ n_d(n_d-1)] / (N(N-1))`; they share their sire
when they share their dam, or else when their dams share a sire. Two genes in
different offspring then come from the same parent with probability
`a = (P_f + P_m)/4`, and the effective size of the generation is about
`1/a + 1/2` (`1/a = N` for `N/2` males and `N/2` females under random mating),
which Hill's (1972) `8N/(V_m + V_f + 4)` approximates.

`nsires` and `ndams` restrict breeding to a subset of the parental cohort:
in every generation, `nsires` of the males and `ndams` of the females, each
a uniformly random subset, are the sires and dams, and the other individuals
have no offspring. The splits are then among these breeders, and the
coalescence probabilities above use `nsires` and `ndams` in place of the
numbers of males and females. With binomial sex, the sexes are redrawn until
there are at least `nsires` males and `ndams` females (exact rejection
sampling, as for the presence of both sexes); with `nmale` or `sexratio`,
every parental cohort must hold them, which is checked before the run.
[`Hierarchical`](@ref) is the nested design with equal splits.

**Selection.** With a fitness model ([`SiteFitness`](@ref),
[`TraitFitness`](@ref)), the weight of every parent in a split is multiplied
by its fitness. With `:random` splits each dam's expected number of
offspring, and each sire's expected number of dams, is proportional to its
fitness within its sex, and so is each sire's expected number of offspring;
with a Gamma shape the weights are Gamma(`k`) × fitness. `:equal` splits, and therefore [`Hierarchical`](@ref), cannot be
combined with selection, as equal families leave nothing to select on.
Breeder subsets (`nsires`, `ndams`) are drawn at random, and selection acts
among the breeders.

The census size is not affected. [`effective_size`](@ref) returns the effective size realized
by the family sizes of a run, `1/a` per generation.

# Examples
```julia
using FisherWright

# Equal dam families, 5 sires with equal numbers of dams
haps, _ = fisher_wright(50, 20, [100_000], 1.0; mating = FamilySize(dams = :equal, sires = :equal, nmale = 5))

# Highly variable dam families
haps, _ = fisher_wright(50, 20, [100_000], 1.0; mating = FamilySize(dams = 0.5))
```
"""
struct FamilySize <: AbstractMating
    sex::RandomMating
    dams::Symbol
    kd::Float64
    sires::Symbol
    km::Float64
    nsires::Int   # breeding sires per generation; 0 for all males
    ndams::Int    # breeding dams per generation; 0 for all females
end

function FamilySize(; dams=:random, sires=:random, nsires::Union{Nothing,Integer}=nothing,
    ndams::Union{Nothing,Integer}=nothing, nmale=nothing, sexratio=nothing)
    d, kd = _split_mode(dams, "dams")
    s, km = _split_mode(sires, "sires")
    nsires === nothing || nsires > 0 || throw(ArgumentError("nsires must be positive"))
    ndams === nothing || ndams > 0 || throw(ArgumentError("ndams must be positive"))
    return FamilySize(RandomMating(nmale=nmale, sexratio=sexratio), d, kd, s, km,
        something(nsires, 0), something(ndams, 0))
end

"""
    Hierarchical(sires, dams_per_sire; nmale = nothing, sexratio = nothing)

Nested mating design: in every generation `sires` males and
`sires * dams_per_sire` females, uniformly random subsets of the parental
cohort, are the breeders; each sire is mated to `dams_per_sire` of the dams,
and the offspring are split equally among the dams (the dams that receive
one more, when the split is not exact, are drawn at random). The other
individuals have no offspring.

The same as `FamilySize(dams = :equal, sires = :equal, nsires = sires,
ndams = sires * dams_per_sire; nmale, sexratio)`; see [`FamilySize`](@ref)
for sex assignment and the effective size.

# Examples
```julia
using FisherWright

# 4 sires with 5 dams each out of 100 individuals per generation
haps, _ = fisher_wright(100, 20, [100_000], 1.0; mating = Hierarchical(4, 5))
```
"""
function Hierarchical(sires::Integer, dams_per_sire::Integer; nmale=nothing, sexratio=nothing)
    (sires > 0 && dams_per_sire > 0) ||
        throw(ArgumentError("sires and dams_per_sire must be positive"))
    return FamilySize(dams=:equal, sires=:equal, nsires=sires, ndams=sires * dams_per_sire,
        nmale=nmale, sexratio=sexratio)
end

# Fixed sex numbers must hold the breeders in every parental cohort; with
# binomial sex, every parental cohort must be large enough for them.
function _check_mating(m::FamilySize, demog::AbstractDemography, n0::Int, nt::Int)
    _check_mating(m.sex, demog, n0, nt)
    (m.nsires == 0 && m.ndams == 0) && return nothing
    n = n0
    for g = 1:nt
        nm = _nmale(m.sex, n)
        if nm > 0
            m.nsires <= nm || throw(ArgumentError(
                "nsires = $(m.nsires) exceeds the $nm males among the parents of generation $g"))
            m.ndams <= n - nm || throw(ArgumentError(
                "ndams = $(m.ndams) exceeds the $(n - nm) females among the parents of generation $g"))
        else
            max(m.nsires, 1) + max(m.ndams, 1) <= n || throw(ArgumentError(
                "the $n parents of generation $g cannot hold $(m.nsires) sires and $(m.ndams) dams"))
        end
        n = _popsize(demog, g)
    end
    return nothing
end

# The first `k` of `v` become a uniformly random subset of size k, sorted
# (partial Fisher-Yates), copied into `out`; all of `v` when k is 0.
function _breeders!(out::Vector{Int}, v::Vector{Int}, k::Int, rng::AbstractRNG)
    resize!(out, length(v))
    copyto!(out, v)
    k == 0 && return out
    @inbounds for i = 1:k
        j = rand(rng, i:length(out))
        out[i], out[j] = out[j], out[i]
    end
    resize!(out, k)
    return sort!(out)
end

# Counts of a split of `total` among `n` parents into `cnt[1:n]`; `w` and `t`
# are scratch for the Gamma weights and their alias table.
function _split!(cnt::Vector{Int}, n::Int, total::Int, mode::Symbol, k::Float64, w::Vector{Float64},
    t::AliasTable, idx::Vector{Int}, rng::AbstractRNG)
    resize!(cnt, n)
    fill!(cnt, 0)
    if mode === :equal
        q, r = divrem(total, n)
        fill!(cnt, q)
        # r of the n parents, a uniformly random subset, get one more.
        need = r
        @inbounds for i = 1:n
            left = n - i + 1
            if need > 0 && (need == left || rand(rng) * left < need)
                cnt[i] += 1
                need -= 1
            end
        end
    elseif mode === :random
        @inbounds for _ = 1:total
            cnt[rand(rng, 1:n)] += 1
        end
    else
        g = Gamma(k)
        @inbounds for i = 1:n
            w[i] = rand(rng, g)
        end
        resize!(idx, n)
        @inbounds for i = 1:n
            idx[i] = i
        end
        # Weights of exactly zero have probability zero for k > 0, but
        # underflow can produce them: fall back to equal weights.
        any(i -> w[i] > 0, 1:n) || fill!(view(w, 1:n), 1.0)
        _alias_build!(t, w, idx)
        @inbounds for _ = 1:total
            cnt[_alias_draw(t, rng)] += 1
        end
    end
    return cnt
end

# As `_split!` among the parents `members`, with the weight of parent i
# multiplied by its fitness `fit[members[i]]`: `:random` splits in proportion
# to fitness, Gamma splits in proportion to Gamma(k) × fitness.
function _split_selected!(cnt::Vector{Int}, members::Vector{Int}, total::Int, mode::Symbol, k::Float64,
    fit::Vector{Float64}, sw::Vector{Float64}, t::AliasTable, idx::Vector{Int}, rng::AbstractRNG, what)
    n = length(members)
    resize!(cnt, n)
    fill!(cnt, 0)
    g = Gamma(mode === :gamma ? k : 1.0)
    @inbounds for i = 1:n
        sw[i] = (mode === :gamma ? rand(rng, g) : 1.0) * fit[members[i]]
    end
    any(i -> sw[i] > 0, 1:n) ||
        error("no $what with positive fitness: the population cannot reproduce")
    resize!(idx, n)
    @inbounds for i = 1:n
        idx[i] = i
    end
    _alias_build!(t, sw, idx)
    @inbounds for _ = 1:total
        cnt[_alias_draw(t, rng)] += 1
    end
    return cnt
end

_mate!(m::FamilySize, f::Neutral, st::SimState, no::Int, rng::AbstractRNG) = _family_mate!(m, f, st, no, rng)
_mate!(m::FamilySize, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG) =
    _family_mate!(m, f, st, no, rng)

function _family_mate!(m::FamilySize, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG)
    resize!(st.sex, st.n)
    nm = _nmale(m.sex, st.n)
    if nm == 0
        _draw_sex!(st.sex, st.sires, st.dams, rng, typemax(Int), max(m.nsires, 1), max(m.ndams, 1))
    else
        _draw_fixed_sex!(st.sex, st.sires, st.dams, nm, rng)
    end
    # The breeders: all males and females, or random subsets of them. The
    # record of males (st.sires) is left as it is.
    sires = _breeders!(st.bsires, st.sires, m.nsires, rng)
    dams = _breeders!(st.bdams, st.dams, m.ndams, rng)
    Nm, Nf = length(sires), length(dams)
    # Offspring per dam, then dams per sire; under selection in proportion
    # to fitness as well.
    f isa Neutral || _check_selection(m, f)
    if f isa Neutral
        nd = _split!(st.fcnt, Nf, no, m.dams, m.kd, st.w, st.tf, st.fidx, rng)
        ms = _split!(st.mcnt, Nm, Nf, m.sires, m.km, st.w, st.tm, st.fidx, rng)
    else
        nd = _split_selected!(st.fcnt, dams, no, m.dams, m.kd, st.w, st.swt, st.tf, st.fidx, rng, "dam")
        ms = _split_selected!(st.mcnt, sires, Nf, m.sires, m.km, st.w, st.swt, st.tm, st.fidx, rng, "sire")
    end
    # The sire of each dam: the slots of the sires, in a uniformly random
    # order (Fisher-Yates).
    sod = st.fidx
    resize!(sod, Nf)
    j = 0
    @inbounds for s = 1:Nm, _ = 1:ms[s]
        sod[j+=1] = s
    end
    @inbounds for i = Nf:-1:2
        k = rand(rng, 1:i)
        sod[i], sod[k] = sod[k], sod[i]
    end
    pm = st.pm
    i = 0
    @inbounds for d = 1:Nf, _ = 1:nd[d]
        i += 1
        pm[i, 1] = sires[sod[d]]
        pm[i, 2] = dams[d]
    end
    return nothing
end

# Selection needs a split that can follow fitness.
function _check_selection(m::FamilySize, f::AbstractFitness)
    (f isa Neutral || (m.dams !== :equal && m.sires !== :equal)) || throw(ArgumentError(
        "FamilySize with :equal splits (and Hierarchical) cannot be combined with selection: " *
        "equal families leave nothing to select on; use :random or a Gamma shape"))
    return nothing
end
