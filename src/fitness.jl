# Selection on an explicit table of fitness loci (docs/dev/roadmap.md,
# "Selection"). Selected alleles live in a second per-haplotype stream,
# `st.sel`, apart from the neutral mutations in `st.prt`: the positions a
# haplotype carries there are its derived alleles at fitness loci.

"""
    FitnessLoci(pos, s; h = 0.5, μ = 0.0)

Table of selected sites, fixed for a run.

- `pos`: genome-wide positions of the sites, strictly increasing.
- `s`: selection coefficient of the derived allele at each site (a number
  for all sites, or one per site), `s ≥ -1`.
- `h`: dominance at each site (a number or one per site), with
  `1 + h * s ≥ 0`.
- `μ`: mutation rate to the derived allele per site per haplotype per
  generation; `0` means the selected alleles of a run are those of its
  starting population.

A genotype at one site has fitness `1`, `1 + h*s` and `1 + s` with zero, one
and two copies of the derived allele; see [`SiteFitness`](@ref) for how sites
combine. Sites are biallelic and there is no back-mutation: a mutation at a
site the haplotype already carries has no effect.

# Fields
`pos::Vector{UInt32}`, `s::Vector{Float64}`, `h::Vector{Float64}`, `μ::Float64`.

# See also
[`fitness_loci`](@ref), [`SiteFitness`](@ref)
"""
struct FitnessLoci
    pos::Vector{UInt32}
    s::Vector{Float64}
    h::Vector{Float64}
    μ::Float64
    function FitnessLoci(pos::Vector{UInt32}, s::Vector{Float64}, h::Vector{Float64}, μ::Float64)
        K = length(pos)
        (length(s) == K && length(h) == K) ||
            throw(ArgumentError("pos, s and h must have the same length"))
        (isempty(pos) || pos[1] >= 1) || throw(ArgumentError("positions are 1-based"))
        all(k -> pos[k-1] < pos[k], 2:K) ||
            throw(ArgumentError("pos must be sorted and unique"))
        all(x -> isfinite(x) && x >= -1, s) ||
            throw(ArgumentError("every s must be finite and at least -1"))
        all(k -> isfinite(h[k]) && 1 + h[k] * s[k] >= 0, 1:K) ||
            throw(ArgumentError("every h must be finite with 1 + h*s ≥ 0"))
        (isfinite(μ) && μ >= 0) || throw(ArgumentError("μ must be finite and non-negative"))
        return new(pos, s, h, μ)
    end
end

_per_site(x::Real, K::Int) = fill(Float64(x), K)
_per_site(x::AbstractVector{<:Real}, ::Int) = Vector{Float64}(x)

FitnessLoci(pos::AbstractVector{<:Integer}, s; h=0.5, μ::Real=0.0) = FitnessLoci(
    Vector{UInt32}(pos), _per_site(s, length(pos)), _per_site(h, length(pos)), Float64(μ))

Base.length(l::FitnessLoci) = length(l.pos)

"""
    fitness_loci(chr, nsites; s, h = 0.5, μ = 0.0, rng = Random.default_rng()) -> FitnessLoci

Draw `nsites` distinct selected sites uniformly over a genome with chromosome
lengths `chr` and build their [`FitnessLoci`](@ref) table.

- `s`: a number, used for every site, or a univariate distribution from
  Distributions.jl, from which one coefficient per site is drawn (the
  distribution of fitness effects), for example `-Gamma(0.2, 0.05)`.
- `h`: a number, or a function `s -> h` applied to each coefficient.
- `μ`: mutation rate per site per haplotype per generation.

The draws come from `rng` only, so the table is reproducible and can be kept
with the simulated data as the true effects.

# Examples
```julia
using FisherWright, Distributions, Random

loci = fitness_loci([1_000_000, 1_000_000], 200; s = -Gamma(0.2, 0.05), μ = 1e-5, rng = Xoshiro(1))
```
"""
function fitness_loci(
    chr::AbstractVector{<:Integer},
    nsites::Integer;
    s,
    h=0.5,
    μ::Real=0.0,
    rng::AbstractRNG=Random.default_rng(),
)
    all(>(0), chr) || throw(ArgumentError("chromosome lengths must be positive"))
    L = sum(Int, chr)
    L < 2^32 || throw(ArgumentError("total genome length must be < 2^32 bp"))
    0 <= nsites <= L || throw(ArgumentError("nsites must lie in 0:$L"))
    # Distinct positions by rejection (few repeats while nsites ≪ L), or a
    # partial shuffle when the sites fill much of the genome.
    pos = if 2nsites <= L
        seen = Set{UInt32}()
        while length(seen) < nsites
            push!(seen, rand(rng, UInt32(1):UInt32(L)))
        end
        sort!(collect(seen))
    else
        sort!(UInt32.(randperm(rng, L)[1:nsites]))
    end
    sv = s isa Real ? fill(Float64(s), Int(nsites)) :
         s isa UnivariateDistribution ? Float64.(rand(rng, s, Int(nsites))) :
         throw(ArgumentError("s must be a number or a univariate distribution"))
    hv = h isa Real ? fill(Float64(h), Int(nsites)) : Float64[h(x) for x in sv]
    return FitnessLoci(pos, sv, hv, Float64(μ))
end

"""
    SiteFitness(loci::FitnessLoci)

Multiplicative fitness over the selected sites of `loci`: an individual's
fitness is the product over sites of `1`, `1 + h*s` or `1 + s` for zero, one
or two copies of the derived allele. Pass it as the `fitness` keyword of
[`fisher_wright`](@ref).

Selection is soft: fitness is evaluated on the adults of each generation, and
each offspring draws its sire and its dam with probability proportional to
fitness within each sex, with replacement. The census size is not affected,
and only fitness relative to the other parents of the same sex matters.
"""
struct SiteFitness <: AbstractFitness
    loci::FitnessLoci
    lhet::Vector{Float64}   # log(1 + h*s)
    lhom::Vector{Float64}   # log(1 + s)
    SiteFitness(loci::FitnessLoci) =
        new(loci, log1p.(loci.h .* loci.s), log1p.(loci.s))
end

# Whether a fitness model keeps selected alleles in `st.sel`.
_selected(::AbstractFitness) = false
_selected(::SiteFitness) = true

# Log fitness of an individual with selected haplotypes `a` and `b`: a
# merge-walk over both, a site on both adding log(1 + s) and a site on one
# log(1 + h*s). Sites are looked up in the table from the last hit onwards.
function _log_fitness(f::SiteFitness, a::Vector{UInt32}, b::Vector{UInt32})
    pos, lhet, lhom = f.loci.pos, f.lhet, f.lhom
    K, la, lb = length(pos), length(a), length(b)
    lw = 0.0
    i = j = k = 1
    @inbounds while i <= la && j <= lb
        x, y = a[i], b[j]
        k = searchsortedfirst(pos, min(x, y), k, K, Base.Order.Forward)
        if x == y
            lw += lhom[k]
            i += 1
            j += 1
        elseif x < y
            lw += lhet[k]
            i += 1
        else
            lw += lhet[k]
            j += 1
        end
    end
    # What is left of one haplotype is heterozygous.
    @inbounds while i <= la
        k = searchsortedfirst(pos, a[i], k, K, Base.Order.Forward)
        lw += lhet[k]
        i += 1
    end
    @inbounds while j <= lb
        k = searchsortedfirst(pos, b[j], k, K, Base.Order.Forward)
        lw += lhet[k]
        j += 1
    end
    return lw
end

# Relative fitness of every parent, scaled so that the fittest has w = 1,
# and the fitness record of the generation: log fitness is absolute, so the
# substituted sites, which every parent carries on both haplotypes, are
# added back.
function _fitness!(f::SiteFitness, st::SimState)
    w, sel, n = st.w, st.sel, st.n
    Threads.@threads for i = 1:n
        @inbounds w[i] = _log_fitness(f, sel[2i-1], sel[2i])
    end
    lw = view(w, 1:n)
    offset = 0.0
    for x in st.ssubs
        offset += f.lhom[searchsortedfirst(f.loci.pos, x)]
    end
    push!(st.flmean, mean(lw) + offset)
    push!(st.flvar, var(lw))
    push!(st.fmean, sum(x -> exp(x + offset), lw) / n)
    top = maximum(lw)
    top > -Inf || error("every individual has zero fitness: the population cannot reproduce")
    @inbounds for i = 1:n
        w[i] = exp(w[i] - top)
    end
    return nothing
end

_record_fitness!(::SiteFitness, ::SimState) = nothing

function _check_fitness(f::SiteFitness, tbp::UInt32)
    (isempty(f.loci.pos) || f.loci.pos[end] <= tbp) ||
        throw(ArgumentError("a fitness locus lies beyond the genome length $(Int(tbp))"))
    return nothing
end

# New mutations at the fitness loci on each haplotype of the current cohort
# (threaded). Haplotype i draws from `st.srng[i]`, after its neutral
# mutations: a Poisson number of hits on sites drawn uniformly from the
# table. A hit on a site the haplotype carries, or on a substituted site, is
# a no-op.
_sel_mutation!(f::SiteFitness, st::SimState, p_sel::Poisson{Float64}) = _sel_mutation!(f.loci.pos, st, p_sel)

function _sel_mutation!(pos::Vector{UInt32}, st::SimState, p_sel::Poisson{Float64})
    K = length(pos)
    sel, mbuf, nbuf, srng, subs = st.sel, st.smbuf, st.snbuf, st.srng, st.ssubs
    Threads.@threads for i = 1:2*st.n
        r = srng[i]
        nm = rand(r, p_sel)
        if nm > 0
            newm = nbuf[i]
            resize!(newm, nm)
            @inbounds for k = 1:nm
                newm[k] = pos[rand(r, 1:K)]
            end
            nm <= 32 ? sort!(newm; alg=InsertionSort) : sort!(newm)
            _drop_present!(newm, subs)
            if !isempty(newm)
                _insert_sorted!(mbuf[i], sel[i], newm)
                sel[i], mbuf[i] = mbuf[i], sel[i]
            end
        end
    end
    return nothing
end

# Number of new selected mutations per haplotype per generation, or nothing
# when the model has none.
_sel_poisson(::AbstractFitness) = nothing
_sel_poisson(f::SiteFitness) =
    f.loci.μ > 0 && length(f.loci) > 0 ? Poisson(length(f.loci) * f.loci.μ) : nothing

# Selected alleles must be sites of the table.
function _check_selected(v::Vector{UInt32}, loci, what)
    pos = loci.pos
    K = length(pos)
    k = 1
    for (m, x) in pairs(v)
        (m == 1 || v[m-1] < x) || throw(ArgumentError("$what: positions must be sorted and unique"))
        k = searchsortedfirst(pos, x, k, K, Base.Order.Forward)
        (k <= K && pos[k] == x) ||
            throw(ArgumentError("$what: position $(Int(x)) is not a site of the table of loci"))
    end
    return nothing
end

# Starting selected haplotypes and selected substitutions; copies of `init`.
_initial_selected(::Union{Nothing,Vector{Vector{UInt32}}}, ::AbstractFitness, nh::Int) =
    ([Vector{UInt32}() for _ = 1:nh], UInt32[])

function _initial_selected(init::FisherWrightResult, f::AbstractFitness, nh::Int)
    sel, subs = init.selected_haplotypes, init.selected_substitutions
    length(sel) == nh ||
        throw(ArgumentError("init has $(length(sel)) selected haplotypes; $nh are required"))
    if !_selected(f)
        (all(isempty, sel) && isempty(subs)) || throw(ArgumentError(
            "init carries selected alleles; pass the fitness model they belong to as `fitness`"))
        return ([Vector{UInt32}() for _ = 1:nh], UInt32[])
    end
    for (i, h) in pairs(sel)
        _check_selected(h, f.loci, "init selected haplotype $i")
    end
    _check_selected(subs, f.loci, "init.selected_substitutions")
    return ([_remove_sorted!(copy(h), subs) for h in sel], copy(subs))
end

"""
    to_haplotype(result::FisherWrightResult, loci::FitnessLoci; include_fixed::Bool = true) -> Tuple{Haplotype, DataFrame}

Genotypes at the fitness loci of a run with selection, with their true effects: the causal variants of the
simulated population.

Row `k` of the returned `BnGStructs.Haplotype` is site `loci.pos[k]`, for every site of the table whether or
not it segregates, and the columns are the haplotypes; an entry is `true` where the haplotype carries the
derived allele. Neutral mutations are not included, also not at the same coordinate.

# Keywords
- `include_fixed::Bool = true`: Sites in `result.selected_substitutions` are carried by every haplotype and are
  exported as such. With `false` they are exported as absent, like a lost allele. Note that the default differs
  from the other methods of `to_haplotype`: here a fixed site is a row of the table either way, and `true` is
  its genotype.

# Returns
- `Haplotype`: `length(loci)` loci × `2N` haplotypes.
- `DataFrame`: one row per site, with columns
  - `chr::Int8`: 1-based chromosome,
  - `pos::UInt32`: genome-wide position,
  - `s::Float64`, `h::Float64`: selection coefficient and dominance of the derived allele,
  - `frq::Float64`: frequency of the derived allele in the exported haplotypes.

# Errors
- Throws an `ArgumentError` if `loci` is empty, if a site lies beyond the genome of `result`, or if `result`
  carries a selected allele that is not a site of `loci`.

# Examples
```julia
using FisherWright

chr = [100_000, 100_000]
loci = fitness_loci(chr, 20; s = 0.05, μ = 1e-4)
res = fisher_wright(50, 200, chr, 1.0; fitness = SiteFitness(loci), result = true)
qtl, map = to_haplotype(res, loci)
dosage = [qtl[k, 2i-1] + qtl[k, 2i] for k in 1:length(loci), i in 1:50]   # sites × individuals
```
"""
# Genotypes at the sites of a table of loci, rows in table order.
function _selected_matrix(result::FisherWrightResult, loci, include_fixed::Bool)
    K = length(loci)
    K > 0 || throw(ArgumentError("the table of loci is empty"))
    cbp = result.chromosome_ends
    loci.pos[end] <= cbp[end] ||
        throw(ArgumentError("a locus of the table lies beyond the genome length $(Int(cbp[end]))"))
    sel, subs = result.selected_haplotypes, result.selected_substitutions
    for (i, h) in pairs(sel)
        _check_selected(h, loci, "selected haplotype $i")
    end
    _check_selected(subs, loci, "selected_substitutions")
    xy = extract_chip_bitarray(sel, loci.pos)
    if include_fixed
        for x in subs
            xy[searchsortedfirst(loci.pos, x), :] .= true
        end
    end
    return xy, cbp
end

function to_haplotype(result::FisherWrightResult, loci::FitnessLoci; include_fixed::Bool=true)
    xy, cbp = _selected_matrix(result, loci, include_fixed)
    map = DataFrame(
        chr=Int8[searchsortedfirst(cbp, p) for p in loci.pos],
        pos=copy(loci.pos),
        s=copy(loci.s),
        h=copy(loci.h),
        frq=vec(sum(xy; dims=2)) ./ size(xy, 2),
    )
    return Haplotype(xy), map
end
