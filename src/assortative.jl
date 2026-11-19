# Phenotypic assortative mating (docs/dev/roadmap.md, "Open after 0.8.0"):
# the dam of an offspring is drawn in proportion to her fitness times a
# Gaussian copula density of her phenotype and that of the sire.

"""
    Assortative(ρ; inner = RandomMating(), trait = 1)

Phenotypic assortative mating with correlation `ρ` between mates (`-1 < ρ < 1`; `ρ < 0` is disassortative
mating). Sexes are drawn as by `inner`, a [`RandomMating`](@ref) (binomial, or a fixed number or proportion of
males). Every offspring draws its sire among the males in proportion to fitness, as under random mating, and
then its dam among the females with probability proportional to

``w_d \\, c_ρ(z_s, z_d)``,

where `w_d` is the fitness of dam `d`, `z_s` and `z_d` are the normal scores of the phenotypes of trait `trait`
of the sire among the males and of the dam among the females (``Φ^{-1}((r - ½)/n)`` for rank `r` of `n`, tied
phenotypes sharing the average of their scores), and

``c_ρ(z_s, z_d) = \\exp\\left(\\frac{2ρ z_s z_d - ρ^2 (z_s^2 + z_d^2)}{2(1 - ρ^2)}\\right) \\Big/ \\sqrt{1 - ρ^2}``

is the density of the Gaussian copula with correlation `ρ`. Without selection the normal scores of the sire and
the dam of an offspring then have correlation `ρ` up to the discreteness of a finite cohort, whatever the
distribution of the phenotype; for a Gaussian phenotype that is the phenotypic correlation of mates. Selection
still acts within each sex through `w`. `ρ = 0` is random mating, with draws in a different order.

The phenotype comes from the fitness model, which must be a [`TraitFitness`](@ref) (or a [`SexSpecific`](@ref)
model with one); for assortative mating on a trait that is not under selection, use
`TraitFitness(loci; σe, fitness = P -> 1.0)`. The realized correlation of the phenotypes of sire and dam over the
matings of each generation is recorded in `records.mate_correlation`.

Draws from `rng` per generation: the sexes, as by `inner`; the sire of every offspring; then the dams, sire by
sire in the order of the sires in the cohort. Building the dam weights of one sire costs O(number of females),
so a generation costs O(distinct sires × females).

# Examples
```julia
using FisherWright

chr = fill(1_000_000, 10)
loci = trait_loci(chr, 1_000; a = 0.05)
res = fisher_wright(500, 50, chr, 1.0; result = true, mating = Assortative(0.5),
                    fitness = TraitFitness(loci; σe = 1.0, fitness = P -> 1.0))
res.records.mate_correlation[end]   # about 0.5
```
"""
struct Assortative <: AbstractMating
    ρ::Float64
    trait::Int
    inner::RandomMating
    function Assortative(ρ::Real; inner::RandomMating=RandomMating(), trait::Integer=1)
        -1 < ρ < 1 || throw(ArgumentError("ρ is a correlation and must lie strictly between -1 and 1"))
        trait >= 1 || throw(ArgumentError("trait must be positive"))
        return new(Float64(ρ), Int(trait), inner)
    end
end

Base.show(io::IO, m::Assortative) =
    print(io, "Assortative(", m.ρ, "; inner = ", m.inner, ", trait = ", m.trait, ")")

_check_mating(m::Assortative, demog::AbstractDemography, n0::Int, nt::Int) = _check_mating(m.inner, demog, n0, nt)

function _check_selection(m::Assortative, f::AbstractFitness)
    T = _model_traits(f)
    T === nothing && throw(ArgumentError(
        "Assortative mates on a phenotype and needs a TraitFitness (fitness = P -> 1.0 for a neutral trait)"))
    m.trait <= T || throw(ArgumentError("Assortative sorts on trait $(m.trait); the trait model has $T"))
    return nothing
end

# Normal scores of the phenotypes `P` of the individuals in `group`, into
# `z`: Φ⁻¹((r - ½)/n) for rank r, tied phenotypes sharing the mean score.
function _normal_scores!(z::Vector{Float64}, P::AbstractVector{Float64}, group::Vector{Int}, ord::Vector{Int})
    n = length(group)
    o = view(ord, 1:n)
    copyto!(o, group)
    sort!(o; by=i -> P[i])
    nd = Normal()
    k = 1
    @inbounds while k <= n
        j = k
        while j < n && P[o[j+1]] == P[o[k]]
            j += 1
        end
        s = 0.0
        for r = k:j
            s += quantile(nd, (r - 0.5) / n)
        end
        s /= j - k + 1
        for r = k:j
            z[o[r]] = s
        end
        k = j + 1
    end
    return z
end

function _mate!(m::Assortative, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG)
    mb = st.mate
    resize!(mb.sex, st.n)
    nm = _nmale(m.inner, st.n)
    if nm == 0
        _draw_sex!(mb.sex, mb.sires, mb.dams, rng, typemax(Int))
    else
        _draw_fixed_sex!(mb.sex, mb.sires, mb.dams, nm, rng)
    end
    sires, dams, pm, w, wf = mb.sires, mb.dams, mb.pm, st.fit.w, st.fit.wf
    (any(i -> w[i] > 0, sires) && any(i -> wf[i] > 0, dams)) || error(
        "no parent with positive fitness in one sex: the population cannot reproduce")
    _alias_build!(mb.tm, w, sires)
    @inbounds for i = 1:no
        pm[i, 1] = sires[_alias_draw(mb.tm, rng)]
    end

    P = view(st.fit.pval, m.trait, :)
    z, aw = mb.zsc, mb.aw
    _normal_scores!(z, P, sires, mb.ord)
    _normal_scores!(z, P, dams, mb.ord)

    # Offspring grouped by sire, in the order of the sires in the cohort.
    o = view(mb.ord, 1:no)
    copyto!(o, 1:no)
    sort!(o; by=i -> pm[i, 1], alg=MergeSort)
    ρ = m.ρ
    a, b = ρ / (1 - ρ^2), ρ^2 / (2(1 - ρ^2))
    k = 1
    @inbounds while k <= no
        s = pm[o[k], 1]
        j = k
        while j < no && pm[o[j+1], 1] == s
            j += 1
        end
        # Log weights of the dams for sire s; the terms in z_s alone and the
        # constant cancel in the normalization.
        zs = z[s]
        top = -Inf
        for d in dams
            lw = wf[d] > 0 ? log(wf[d]) + a * zs * z[d] - b * z[d]^2 : -Inf
            aw[d] = lw
            top = max(top, lw)
        end
        for d in dams
            aw[d] = exp(aw[d] - top)
        end
        _alias_build!(mb.tf, aw, dams)
        for r = k:j
            pm[o[r], 2] = dams[_alias_draw(mb.tf, rng)]
        end
        k = j + 1
    end
    st.rec.mate_correlation[end] = _mate_correlation(P, pm, no)
    return nothing
end

# Pearson correlation of the phenotypes of sire and dam over `no` matings;
# NaN if either has no variance.
function _mate_correlation(P::AbstractVector{Float64}, pm::Matrix{Int}, no::Int)
    ms = md = 0.0
    @inbounds for i = 1:no
        ms += P[pm[i, 1]]
        md += P[pm[i, 2]]
    end
    ms /= no
    md /= no
    sxy = sxx = syy = 0.0
    @inbounds for i = 1:no
        x, y = P[pm[i, 1]] - ms, P[pm[i, 2]] - md
        sxy += x * y
        sxx += x^2
        syy += y^2
    end
    return sxx > 0 && syy > 0 ? sxy / sqrt(sxx * syy) : NaN
end
