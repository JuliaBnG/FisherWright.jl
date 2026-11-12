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
    FamilySize(; dams = :random, sires = :random, nmale = nothing, sexratio = nothing)

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
`1/(2a) + 1/2`, which Hill's (1972) `8N/(V_m + V_f + 4)` approximates.

The census size is not affected. `FamilySize` cannot be combined with
selection yet. [`effective_size`](@ref) assumes random mating and does not
apply to runs with `FamilySize`.

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
end

function FamilySize(; dams=:random, sires=:random, nmale=nothing, sexratio=nothing)
    d, kd = _split_mode(dams, "dams")
    s, km = _split_mode(sires, "sires")
    return FamilySize(RandomMating(nmale=nmale, sexratio=sexratio), d, kd, s, km)
end

_check_mating(m::FamilySize, demog::AbstractDemography, n0::Int, nt::Int) = _check_mating(m.sex, demog, n0, nt)

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

function _mate!(m::FamilySize, ::Neutral, st::SimState, no::Int, rng::AbstractRNG)
    resize!(st.sex, st.n)
    nm = _nmale(m.sex, st.n)
    if nm == 0
        _draw_sex!(st.sex, st.sires, st.dams, rng, typemax(Int))
    else
        _draw_fixed_sex!(st.sex, st.sires, st.dams, nm, rng)
    end
    sires, dams = st.sires, st.dams
    Nm, Nf = length(sires), length(dams)
    # Offspring per dam, then dams per sire.
    nd = _split!(st.fcnt, Nf, no, m.dams, m.kd, st.w, st.tf, st.fidx, rng)
    ms = _split!(st.mcnt, Nm, Nf, m.sires, m.km, st.w, st.tm, st.fidx, rng)
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

_mate!(::FamilySize, ::AbstractFitness, ::SimState, ::Int, ::AbstractRNG) =
    throw(ArgumentError("FamilySize cannot be combined with selection yet"))
