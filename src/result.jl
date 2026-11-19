"""
    GenerationRecords

The per-generation record of a run, the `records` field of a
[`FisherWrightResult`](@ref). With `ngen` recorded generations, generation `g`
being the matings of the parents of census `g` and the birth of census
`g + 1`:

- `census::Vector{Int}` (`ngen + 1`): census sizes ``N_0, N_1, …, N_{ngen}``: `census[1]` is the starting
  population and `census[g + 1]` the number of individuals born in generation `g`.
- `males::Vector{Int}`: `males[g]` is the realized number of males among the `census[g]` parents of generation `g`;
  the remaining `census[g] - males[g]` parents are females. 0 for parents without sexes ([`Monoecious`](@ref)).
- `mean_fitness::Vector{Float64}`, `mean_log_fitness::Vector{Float64}`, `var_log_fitness::Vector{Float64}`:
  fitness of the parents of each generation, entry `g` describing the same `census[g]` adults as `males[g]`:
  the arithmetic mean of fitness, and the mean and the (sample) variance of its logarithm. Fitness is absolute,
  relative to an individual without any derived allele at the fitness loci, so substituted sites are included.
  In a neutral run these are `1`, `0` and `0`. If an adult has zero fitness (a lethal genotype), the mean of the
  logarithm is `-Inf` and its variance `NaN`. Under [`SexSpecific`](@ref), each parent under the model of its sex.
- `mean_genotypic::Matrix{Float64}`, `mean_phenotypic::Matrix{Float64}` (`T × ngen`), `var_genotypic::Array{Float64,3}`,
  `var_phenotypic::Array{Float64,3}` (`T × T × ngen`): under [`TraitFitness`](@ref) with `T` traits, the mean vector
  and the (sample) covariance matrix of the genotypic values and of the phenotypes of the parents of each generation,
  as for the fitness record: `mean_genotypic[t, g]` is the mean of trait `t` in generation `g` and
  `var_genotypic[:, :, g]` the genotypic covariance matrix (`var_genotypic[t, t, g]` the variance of trait `t`).
  `NaN` for generations without a trait model; `T = 1` ([`ntraits`](@ref)) for a run that never had one.
- `paternal_pairs::Vector{Float64}`, `maternal_pairs::Vector{Float64}`: the realized family sizes of each
  generation, as the number of ordered pairs of its offspring that share their sire, `Σ m(m - 1)` over the sires,
  and that share their dam, `Σ n(n - 1)` over the dams; [`effective_size`](@ref) uses them. `NaN` where not
  recorded. With hermaphrodites ([`Monoecious`](@ref)), `paternal_pairs` counts the ordered pairs `(i, j)` of
  distinct offspring in which the first parent of `i` is a parent of `j`, once for each role it has in `j`, and
  `maternal_pairs` the same for the second parent of `i`; with separate sexes these are the sums above.
- `selfing::Vector{Float64}`: the proportion of the offspring of each generation that are selfed (one parent in
  both roles); 0 with separate sexes.
- `related_matings::Vector{Int}`: under [`AvoidRelatives`](@ref), the number of offspring of each generation whose
  parents are relatives at the level avoided, kept because no unrelated mate was found within `max_tries` draws;
  0 under other schemes, which do not check relationships.

    GenerationRecords(n; ntraits = 1)
    GenerationRecords(; census, males, mean_fitness, ...)

The first form is the record of a population of `n` individuals without recorded generations. The second takes
the fields by keyword: `census` and `males` are required, and every other field defaults to its value for
generations without selection, trait model or family record (fitness `1`, `0`, `0`; traits `NaN` with
`ntraits` traits; family sizes `NaN`; selfing and related matings `0`). The lengths are checked.
"""
struct GenerationRecords
    census::Vector{Int}
    males::Vector{Int}
    mean_fitness::Vector{Float64}
    mean_log_fitness::Vector{Float64}
    var_log_fitness::Vector{Float64}
    mean_genotypic::Matrix{Float64}
    var_genotypic::Array{Float64,3}
    mean_phenotypic::Matrix{Float64}
    var_phenotypic::Array{Float64,3}
    paternal_pairs::Vector{Float64}
    maternal_pairs::Vector{Float64}
    selfing::Vector{Float64}
    related_matings::Vector{Int}
    function GenerationRecords(census, males, mean_fitness, mean_log_fitness, var_log_fitness,
        mean_genotypic, var_genotypic, mean_phenotypic, var_phenotypic, paternal_pairs, maternal_pairs, selfing,
        related_matings)
        ng = length(males)
        T = size(mean_genotypic, 1)
        length(census) == ng + 1 ||
            throw(ArgumentError("census must have one entry more than males (the starting population)"))
        all(v -> length(v) == ng, (mean_fitness, mean_log_fitness, var_log_fitness, paternal_pairs,
            maternal_pairs, selfing, related_matings)) &&
            size(mean_genotypic) == size(mean_phenotypic) == (T, ng) &&
            size(var_genotypic) == size(var_phenotypic) == (T, T, ng) && T >= 1 ||
            throw(ArgumentError("every record must have one entry per generation, " *
                "the trait records T × ngen and T × T × ngen"))
        return new(census, males, mean_fitness, mean_log_fitness, var_log_fitness, mean_genotypic, var_genotypic,
            mean_phenotypic, var_phenotypic, paternal_pairs, maternal_pairs, selfing, related_matings)
    end
end

function GenerationRecords(;
    census::AbstractVector{<:Integer},
    males::AbstractVector{<:Integer},
    ntraits::Integer=1,
    mean_fitness=ones(length(males)),
    mean_log_fitness=zeros(length(males)),
    var_log_fitness=zeros(length(males)),
    mean_genotypic=fill(NaN, ntraits, length(males)),
    var_genotypic=fill(NaN, ntraits, ntraits, length(males)),
    mean_phenotypic=fill(NaN, ntraits, length(males)),
    var_phenotypic=fill(NaN, ntraits, ntraits, length(males)),
    paternal_pairs=fill(NaN, length(males)),
    maternal_pairs=fill(NaN, length(males)),
    selfing=zeros(length(males)),
    related_matings=zeros(Int, length(males)),
)
    return GenerationRecords(Vector{Int}(census), Vector{Int}(males), Vector{Float64}(mean_fitness),
        Vector{Float64}(mean_log_fitness), Vector{Float64}(var_log_fitness), Matrix{Float64}(mean_genotypic),
        Array{Float64,3}(var_genotypic), Matrix{Float64}(mean_phenotypic), Array{Float64,3}(var_phenotypic),
        Vector{Float64}(paternal_pairs), Vector{Float64}(maternal_pairs), Vector{Float64}(selfing),
        Vector{Int}(related_matings))
end

GenerationRecords(n::Integer; ntraits::Integer=1) = GenerationRecords(; census=[Int(n)], males=Int[], ntraits=ntraits)

Base.:(==)(a::GenerationRecords, b::GenerationRecords) =
    all(f -> getfield(a, f) == getfield(b, f), fieldnames(GenerationRecords))
Base.isequal(a::GenerationRecords, b::GenerationRecords) =
    all(f -> isequal(getfield(a, f), getfield(b, f)), fieldnames(GenerationRecords))
Base.hash(r::GenerationRecords, h::UInt) = hash(map(f -> getfield(r, f), fieldnames(GenerationRecords)), h)

ngenerations(r::GenerationRecords) = length(r.males)

"""
    FisherWrightResult(active_haplotypes, chromosome_ends, substitutions;
                       selected_haplotypes, selected_substitutions, records, first_id = 1, pedigree)

The output of a Fisher-Wright simulation: the haplotypes of the last generation, the chromosome ends, the fixed
substitutions, and the record of the run.

# Fields
- `active_haplotypes::Vector{Vector{UInt32}}`: Active haplotype mutation positions (two per individual of the
  last generation). Each inner vector contains sorted, unique 1-based base pair coordinates of segregating mutations.
- `chromosome_ends::Vector{UInt32}`: Cumulative chromosome end positions in base pairs (1-based, strictly increasing).
- `substitutions::Vector{UInt32}`: Sorted, unique positions of mutations that reached fixation across all
  haplotypes and were removed from `active_haplotypes` during simulation.
- `selected_haplotypes::Vector{Vector{UInt32}}`: for each haplotype, the sorted positions of the derived alleles
  it carries at the fitness loci of a run with selection ([`SiteFitness`](@ref), [`TraitFitness`](@ref)); empty
  vectors in a neutral run. These alleles are kept apart from the neutral mutations in `active_haplotypes`.
- `selected_substitutions::Vector{UInt32}`: sorted positions of fitness loci whose derived allele reached fixation
  and was removed from `selected_haplotypes`.
- `records::GenerationRecords`: the per-generation record of the run (census, sexes, fitness, traits, family
  sizes, selfing, related matings); see [`GenerationRecords`](@ref).
- `first_id::Int`: the ID of individual 1 of `active_haplotypes`; individual `i` (haplotypes `2i - 1` and `2i`)
  has ID `first_id + i - 1`. Individuals are numbered in birth order across the run, the starting population
  first, and a run continued through `init` continues the numbering. 1 in a result built by hand.
- `pedigree::DataFrame`: with the `pedigree = k` keyword of [`fisher_wright`](@ref), the matings of the last `k`
  generations, one row per individual: `id`, `sire` and `dam` (IDs; 0 for the founders of the record, which are
  the parents of the oldest recorded generation, or the starting population; with hermaphrodites the first and
  the second parent, the same for a selfed individual), `sex` (`true` for male, the sex the individual had when
  its generation reproduced; `missing` for the last generation and for hermaphrodites) and `generation` (of
  birth, numbered as in `census`). Empty without the keyword.

The fields of the record were fields of the result itself before version 0.10. `result.census` and the other
record names still work, with a deprecation warning, and will be removed in 1.0: use `result.records.census`.

A run continued from a `FisherWrightResult` through the `init` keyword of [`fisher_wright`](@ref) appends its
generations to the record of that result, so generations are numbered from the first call; started from plain
haplotypes, the record begins anew. The constructor builds a result without recorded generations by default
(`records = GenerationRecords(length(active_haplotypes) ÷ 2)`); its keywords give, for example, the selected
alleles of a starting population, such as one copy of a new allele. A result is saved and read back with
[`write_result`](@ref) and [`read_result`](@ref).

# See also
[`fisher_wright`](@ref), [`effective_size`](@ref), [`to_haplotype`](@ref)
"""
struct FisherWrightResult
    active_haplotypes::Vector{Vector{UInt32}}
    chromosome_ends::Vector{UInt32}
    substitutions::Vector{UInt32}
    selected_haplotypes::Vector{Vector{UInt32}}
    selected_substitutions::Vector{UInt32}
    records::GenerationRecords
    first_id::Int
    pedigree::DataFrame
    function FisherWrightResult(active_haplotypes, chromosome_ends, substitutions, selected_haplotypes,
        selected_substitutions, records::GenerationRecords, first_id::Integer, pedigree::DataFrame)
        length(selected_haplotypes) == length(active_haplotypes) ||
            throw(ArgumentError("selected_haplotypes must have one entry per haplotype"))
        return new(active_haplotypes, chromosome_ends, substitutions, selected_haplotypes, selected_substitutions,
            records, first_id, pedigree)
    end
end

_empty_pedigree() = DataFrame(id=Int[], sire=Int[], dam=Int[], sex=Union{Missing,Bool}[], generation=Int[])

FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions;
    selected_haplotypes=[UInt32[] for _ in active_haplotypes],
    selected_substitutions=UInt32[],
    records::GenerationRecords=GenerationRecords(length(active_haplotypes) ÷ 2),
    first_id::Integer=1,
    pedigree::DataFrame=_empty_pedigree(),
) = FisherWrightResult(active_haplotypes, chromosome_ends, substitutions, selected_haplotypes,
    selected_substitutions, records, first_id, pedigree)

# The positional constructors of 0.9 (5 to 20 arguments, the record fields
# in the result), deprecated.
const _RECORD_FIELDS = fieldnames(GenerationRecords)
const _OLD_POSITIONAL = (:census, :males, :selected_haplotypes, :selected_substitutions, :mean_fitness,
    :mean_log_fitness, :var_log_fitness, :mean_genotypic, :var_genotypic, :mean_phenotypic, :var_phenotypic,
    :paternal_pairs, :maternal_pairs, :selfing, :related_matings, :first_id, :pedigree)
function FisherWrightResult(active_haplotypes, chromosome_ends, substitutions, census, males, rest...)
    Base.depwarn("the positional FisherWrightResult constructor with the record fields is deprecated; " *
        "use FisherWrightResult(haps, ends, subs; records = GenerationRecords(; census, males, ...), ...)",
        :FisherWrightResult)
    length(rest) <= 15 || throw(MethodError(FisherWrightResult,
        (active_haplotypes, chromosome_ends, substitutions, census, males, rest...)))
    kw = Dict{Symbol,Any}(zip(_OLD_POSITIONAL, (census, males, rest...)))
    sel = get(kw, :selected_haplotypes, [UInt32[] for _ in active_haplotypes])
    ssubs = get(kw, :selected_substitutions, UInt32[])
    first_id = get(kw, :first_id, 1)
    pedigree = get(kw, :pedigree, _empty_pedigree())
    rkw = (k => v for (k, v) in kw if k in _RECORD_FIELDS)
    T = haskey(kw, :mean_genotypic) ? size(kw[:mean_genotypic], 1) : 1
    return FisherWrightResult(active_haplotypes, chromosome_ends, substitutions, sel, ssubs,
        GenerationRecords(; ntraits=T, rkw...), first_id, pedigree)
end

# The record fields of 0.9, read through the result: deprecated.
@inline function Base.getproperty(r::FisherWrightResult, s::Symbol)
    if s in _RECORD_FIELDS
        Base.depwarn("result.$s is deprecated; use result.records.$s", :getproperty)
        return getfield(getfield(r, :records), s)
    end
    return getfield(r, s)
end

Base.:(==)(a::FisherWrightResult, b::FisherWrightResult) =
    all(f -> getfield(a, f) == getfield(b, f), fieldnames(FisherWrightResult))
Base.isequal(a::FisherWrightResult, b::FisherWrightResult) =
    all(f -> isequal(getfield(a, f), getfield(b, f)), fieldnames(FisherWrightResult))
Base.hash(r::FisherWrightResult, h::UInt) = hash(map(f -> getfield(r, f), fieldnames(FisherWrightResult)), h)

ntraits(r::GenerationRecords) = size(r.mean_genotypic, 1)
ntraits(r::FisherWrightResult) = ntraits(r.records)

"""
    effective_size(result::FisherWrightResult, generations = eachindex(result.records.males); from = :families) -> Float64
    effective_size(records::GenerationRecords, generations; from = :families) -> Float64

Effective population size realized in the generations of `result` that are given (one generation or a collection
of them; all generations of the run by default), from its record.

With `from = :families` (the default), from the realized family sizes. Two genes in different offspring of
generation `g` derive from the same parent with probability

``a_g = \\frac{\\sum_s m_s(m_s - 1) + \\sum_d n_d(n_d - 1)}{4 N_g (N_g - 1)}``,

where `m_s` and `n_d` are the numbers of offspring of sire `s` and dam `d` (`result.records.paternal_pairs[g]` and
`result.records.maternal_pairs[g]`) and `N_g = result.records.census[g + 1]`. The effective size of the generation is `1/a_g`, and
over several generations it is the harmonic mean, `1/ā`. This holds for every mating scheme
([`RandomMating`](@ref), [`FamilySize`](@ref), [`Hierarchical`](@ref)) and includes the variance of family size
that selection adds within a generation; it does not include the cumulative effect of heritable fitness over
generations (Robertson 1961). It is the value realized by the pedigree of the run, not an estimate from the
genetic data, and it varies from run to run as family sizes do.

With hermaphrodites ([`Monoecious`](@ref)), `a_g` counts a parent shared through either role, and a proportion
`σ_g = result.records.selfing[g]` of the offspring are selfed. Two genes in one individual then descend from one parent
with probability `σ_g`, so that two genes of different individuals, followed back, coalesce at the rate
`a_g/(2 - σ_g)` per generation instead of `a_g/2`, and the effective size of the generation is `(2 - σ_g)/(2a_g)`
(Nordborg and Donnelly 1997): `N/(1 + F)` with `F = σ/(2 - σ)` at a constant selfing rate, and about `N` under
random selfing. With separate sexes `σ_g = 0` and this is `1/a_g`; over several generations, the harmonic mean.

With `from = :sexes`, from the recorded numbers of male and female parents only, as in version 0.5 and 0.6:
`a_g = 1/(4N_m) + 1/(4N_f)`, the expectation under random mating, which gives Wright's (1931)
`4N_mN_f/(N_m + N_f)`. It does not apply to other mating schemes or to selection, and it throws for a generation of
hermaphrodites.

Under random mating the two agree on average. At constant numbers the asymptotic rate of loss of heterozygosity,
`1/(2N_e)`, corresponds to `N_e ≈ 1/a + 1/2`, because the two genes of one individual come from different parents.

# Errors
- Throws an `ArgumentError` if `generations` is empty, as for a result built without a census record, or if the
  family sizes of a generation are not recorded (`NaN`, a result built by hand) with `from = :families`; and a
  `BoundsError` for a generation outside the run.

# Examples
```julia
using FisherWright

res = fisher_wright(50, 20, [100_000], 1.0; mating = RandomMating(nmale = 5), result = true)
effective_size(res)                 # realized, about 18
effective_size(res; from = :sexes)  # 4 * 5 * 45 / 50 = 18.0
effective_size(res, 1:10)           # the first ten generations
```
"""
effective_size(result::FisherWrightResult, generations=eachindex(result.records.males); from::Symbol=:families) =
    effective_size(result.records, generations; from=from)

function effective_size(rec::GenerationRecords, generations=eachindex(rec.males); from::Symbol=:families)
    isempty(generations) && throw(ArgumentError("no generation to average over"))
    from in (:families, :sexes) || throw(ArgumentError("from must be :families or :sexes"))
    a = 0.0
    for g in generations
        if from === :sexes
            m = rec.males[g]
            m > 0 || throw(ArgumentError(
                "generation $g has no male parents (hermaphrodites); use from = :families"))
            f = rec.census[g] - m
            a += 1 / 4m + 1 / 4f
        else
            pp, mp = rec.paternal_pairs[g], rec.maternal_pairs[g]
            (isnan(pp) || isnan(mp)) &&
                throw(ArgumentError("the family sizes of generation $g are not recorded; use from = :sexes"))
            n = rec.census[g+1]
            a += (pp + mp) / (4n * (n - 1)) * 2 / (2 - rec.selfing[g])
        end
    end
    return length(generations) / a
end
"""
    _intersect_sorted!(dest::Vector{UInt32}, cand::AbstractVector{UInt32}, hap::AbstractVector{UInt32}) -> Vector{UInt32}

Compute the intersection `cand ∩ hap` and write the result into `dest`, returning `dest`.
Both `cand` and `hap` must be sorted in ascending order and unique. The output in `dest` is
sorted and unique.

# Details
- Uses binary search when `8 * length(cand) < length(hap)` with ``O(|\\text{cand}| \\log |\\text{hap}|)`` time.
- Uses linear two-pointer merge when lengths are comparable with ``O(|\\text{cand}| + |\\text{hap}|)`` time.
- Allocates no memory beyond resizing `dest`.
"""
function _intersect_sorted!(
    dest::Vector{UInt32},
    cand::AbstractVector{UInt32},
    hap::AbstractVector{UInt32},
)
    empty!(dest)
    m, n = length(cand), length(hap)
    (m == 0 || n == 0) && return dest
    if 8m < n # few candidates: probe the long vector
        lo = 1
        @inbounds for x in cand
            k = searchsortedfirst(hap, x, lo, n, Base.Order.Forward)
            k > n && break # every later candidate is larger still
            if hap[k] == x
                push!(dest, x)
                lo = k + 1
            else
                lo = k
            end
        end
    else # comparable lengths: linear merge
        i = j = 1
        @inbounds while i <= m && j <= n
            a, b = cand[i], hap[j]
            if a < b
                i += 1
            elseif a > b
                j += 1
            else
                push!(dest, a)
                i += 1
                j += 1
            end
        end
    end
    return dest
end

"""
    _remove_sorted!(hap::Vector{UInt32}, rm::AbstractVector{UInt32}) -> Vector{UInt32}

Delete every position present in `rm` from `hap` in place and return `hap`.
Both `hap` and `rm` must be sorted in ascending order and unique.
Runs in ``O(|\\text{hap}| + |\\text{rm}|)`` time with zero allocations.
"""
function _remove_sorted!(hap::Vector{UInt32}, rm::AbstractVector{UInt32})
    (isempty(hap) || isempty(rm)) && return hap
    n, m = length(hap), length(rm)
    w = 0
    j = 1
    @inbounds for i = 1:n
        x = hap[i]
        while j <= m && rm[j] < x
            j += 1
        end
        (j <= m && rm[j] == x) && continue
        w += 1
        hap[w] = x
    end
    resize!(hap, w)
    return hap
end

"""
    _drop_present!(v::Vector{UInt32}, set::AbstractVector{UInt32}) -> Vector{UInt32}

Delete from `v` in place every element also present in the sorted vector `set`, and return `v`.
Uses binary search per element, ``O(|v| \\log |set|)``, for a short `v` probed against a long `set`.
"""
function _drop_present!(v::Vector{UInt32}, set::AbstractVector{UInt32})
    isempty(set) && return v
    n = length(set)
    w = 0
    @inbounds for x in v
        k = searchsortedfirst(set, x)
        (k <= n && set[k] == x) && continue
        w += 1
        v[w] = x
    end
    resize!(v, w)
    return v
end

"""
    _fixed_mutations(haplotypes::Vector{Vector{UInt32}}) -> Vector{UInt32}

Identify and return sorted mutation positions carried by every haplotype in `haplotypes`.
Each inner vector in `haplotypes` must be sorted and unique.
"""
function _fixed_mutations(haplotypes::AbstractVector{Vector{UInt32}})
    isempty(haplotypes) && return UInt32[]
    fixed = copy(haplotypes[1])
    isempty(fixed) && return fixed
    scratch = UInt32[]
    @inbounds for k = 2:length(haplotypes)
        _intersect_sorted!(scratch, fixed, haplotypes[k])
        fixed, scratch = scratch, fixed
        isempty(fixed) && break
    end
    return fixed
end

"""
    _drop_mutations!(haplotypes::Vector{Vector{UInt32}}, mutations::Vector{UInt32}) -> Vector{Vector{UInt32}}

Remove `mutations` from every haplotype in `haplotypes` in place and return `haplotypes`.
Parallelized across haplotypes via `Threads.@threads`.
"""
function _drop_mutations!(haplotypes::AbstractVector{Vector{UInt32}}, mutations::Vector{UInt32})
    isempty(mutations) && return haplotypes
    Threads.@threads for i in eachindex(haplotypes)
        _remove_sorted!(haplotypes[i], mutations)
    end
    return haplotypes
end

"""
    _drop_mutations(haplotypes::Vector{Vector{UInt32}}, mutations::Vector{UInt32}) -> Vector{Vector{UInt32}}

Non-mutating version of [`_drop_mutations!`](@ref): returns a fresh vector of haplotypes with
`mutations` removed, leaving the input unmodified.
"""
function _drop_mutations(haplotypes::Vector{Vector{UInt32}}, mutations::Vector{UInt32})
    isempty(mutations) && return haplotypes
    cleaned = Vector{Vector{UInt32}}(undef, length(haplotypes))
    for (i, haplotype) in pairs(haplotypes)
        cleaned[i] = _remove_sorted!(copy(haplotype), mutations)
    end
    return cleaned
end

"""
    _fixation_step!(haplotypes::Vector{Vector{UInt32}}, substitutions::Vector{UInt32}) -> Vector{UInt32}

Extract positions carried by all `haplotypes`, remove them from `haplotypes` in place,
and return the updated, sorted substitution list merged with newly fixed mutations.
"""
function _fixation_step!(
    haplotypes::AbstractVector{Vector{UInt32}},
    substitutions::Vector{UInt32},
)
    fixed = _fixed_mutations(haplotypes)
    isempty(fixed) && return substitutions
    _drop_mutations!(haplotypes, fixed)
    return merge_sorted(substitutions, fixed)
end

"""
    _fixation_step(haplotypes::Vector{Vector{UInt32}}, substitutions::Vector{UInt32}) -> Tuple{Vector{Vector{UInt32}}, Vector{UInt32}}

Non-mutating version of [`_fixation_step!`](@ref), returning `(cleaned_haplotypes, updated_substitutions)`.
"""
function _fixation_step(
    haplotypes::Vector{Vector{UInt32}},
    substitutions::Vector{UInt32},
)
    fixed = _fixed_mutations(haplotypes)
    isempty(fixed) && return haplotypes, substitutions
    cleaned = _drop_mutations(haplotypes, fixed)
    updated = merge_sorted(substitutions, fixed)
    return cleaned, updated
end

"""
    _with_fixed(result::FisherWrightResult) -> Vector{Vector{UInt32}}

Return a copy of active haplotypes from `result` with all fixed `substitutions`
merged back into each haplotype in sorted order.
"""
function _with_fixed(result::FisherWrightResult)
    isempty(result.substitutions) && return result.active_haplotypes
    fixed = sort!(unique!(copy(result.substitutions)))
    merged = Vector{Vector{UInt32}}(undef, length(result.active_haplotypes))
    for (i, haplotype) in pairs(result.active_haplotypes)
        merged[i] = merge_sorted(haplotype, fixed)
    end
    return merged
end

# Haplotypes for dense export: the neutral mutations, with the substitutions
# when `include_fixed`, and merged with the selected alleles (and their
# substitutions) when `selected`. A coordinate carried in both streams
# appears once.
function _export_haplotypes(result::FisherWrightResult, include_fixed::Bool, selected::Bool)
    muts = include_fixed ? _with_fixed(result) : result.active_haplotypes
    selected || return muts
    sel = result.selected_haplotypes
    length(sel) == length(muts) ||
        throw(ArgumentError("selected_haplotypes must have one entry per haplotype"))
    subs = result.selected_substitutions
    return [merge_sorted(include_fixed ? merge_sorted(sel[i], subs) : sel[i], muts[i]) for i in eachindex(muts)]
end

"""
    _extract_fixed(result::FisherWrightResult) -> FisherWrightResult

Scan `result.active_haplotypes` for fixed mutations, remove them from active haplotypes,
and return a new `FisherWrightResult` with updated haplotypes and accumulated substitutions.
"""
function _extract_fixed(result::FisherWrightResult)
    cleaned, substitutions = _fixation_step(result.active_haplotypes, result.substitutions)
    return FisherWrightResult(cleaned, result.chromosome_ends, substitutions, result.selected_haplotypes,
        result.selected_substitutions, result.records, result.first_id, result.pedigree)
end


"""
    to_haplotype(result::FisherWrightResult; include_fixed = false, selected = false,
                 rng = Random.default_rng()) -> Tuple{Haplotype, VariantMap}

Convert a structured Fisher-Wright simulation result into a dense `BnGStructs.Haplotype` and the
`BnGStructs.VariantMap` of its loci.

# Arguments
- `result::FisherWrightResult`: Structured simulation result containing active haplotypes, chromosome ends, and substitutions.

# Keywords
- `include_fixed::Bool = false`: If `true`, merges fixed `substitutions` back into each haplotype before matrix construction.
  If `false`, only segregating (polymorphic) loci are exported.
- `selected::Bool = false`: If `true`, the derived alleles at the fitness loci of a run with selection
  (`selected_haplotypes`, and `selected_substitutions` with `include_fixed`) are exported together with the
  neutral mutations, as sequence data that contain the causal variants. A coordinate that carries both a neutral
  mutation and a selected allele is one locus, derived if either is. For the selected sites alone, with their
  effects, use `to_haplotype(result, loci::FitnessLoci)`.
- `rng = Random.default_rng()`: the generator of the placeholder nucleotides `ref` and `alt` (the model has
  none); give one for a reproducible map. The haplotypes do not depend on it.

# Returns
- `Haplotype`: A `BnGStructs.Haplotype` wrapping a `BitMatrix` of shape `(nlc, nhp)` (loci in rows, haplotypes in columns).
- `VariantMap`: the loci, with `chr` (1-based chromosome), `pos` (genome-wide base-pair coordinate, the
  BnGStructs convention), `ref` and `alt` (random distinct nucleotides) and `frq` (frequency of the derived
  allele among the haplotypes); `morgan` is empty, since a result does not keep the genetic map. Before version
  0.10 this was a `DataFrame` with the same columns; `DataFrame(chr = m.chr, pos = m.pos, ref = m.ref,
  alt = m.alt, frq = m.frq)` rebuilds it.

# Errors
- Throws an `ArgumentError` if `include_fixed = false` and no polymorphic loci exist across the haplotypes.

# Examples
```julia
using FisherWright

res = fisher_wright(20, 10, [100_000, 100_000], 1.0; result = true)
hap, vm = to_haplotype(res)
```
"""
function to_haplotype(result::FisherWrightResult; include_fixed::Bool=false, selected::Bool=false,
    rng::AbstractRNG=Random.default_rng())
    muts = _export_haplotypes(result, include_fixed, selected)
    xy, loci = muts2bitarray(muts, result.chromosome_ends; include_fixed=include_fixed, rng=rng)
    isempty(xy) && throw(ArgumentError("no polymorphic loci available for dense export"))
    vm = VariantMap(Vector{Int8}(loci.chr), Vector{UInt32}(loci.pos), Vector{Char}(loci.ref),
        Vector{Char}(loci.alt), Vector{Float32}(loci.frq), Float64[])
    return Haplotype(xy), vm
end

"""
    to_haplotype(result::FisherWrightResult, chip_positions::Vector{UInt32}; include_fixed::Bool = false, selected::Bool = false) -> Haplotype

Directly extract a dense `BnGStructs.Haplotype` at specified coordinate positions `chip_positions`
without materializing other simulated loci.

# Arguments
- `result::FisherWrightResult`: Structured simulation result.
- `chip_positions::Vector{UInt32}`: Sorted, unique 1-based base pair coordinates to extract.

# Keywords
- `include_fixed::Bool = false`: If `true`, includes fixed substitutions from `result` in the extraction.
- `selected::Bool = false`: If `true`, the selected alleles of a run with selection count as carried, as in the method above.

# Returns
- `Haplotype`: A `BnGStructs.Haplotype` wrapping a `BitMatrix` of shape `(length(chip_positions), 2*ne)`.
  Any requested coordinate not present in a haplotype is filled with 0 (false).

# Errors
- Throws an `ArgumentError` if `chip_positions` is not strictly sorted and unique.

# Examples
```julia
using FisherWright

res = fisher_wright(20, 10, [100_000], 1.0; result = true)
chip_pos = UInt32[10_000, 25_000, 50_000, 75_000]
hap = to_haplotype(res, chip_pos)
```
"""
function to_haplotype(result::FisherWrightResult, chip_positions::Vector{UInt32}; include_fixed::Bool=false, selected::Bool=false)
    muts = _export_haplotypes(result, include_fixed, selected)
    xy = extract_chip_bitarray(muts, chip_positions)
    return Haplotype(xy)
end

"""
    to_haplotype(result::FisherWrightResult, locus_set::BnGStructs.LocusSet, all_positions::Vector{UInt32}; include_fixed::Bool = false, selected::Bool = false) -> Haplotype

Directly extract a dense `BnGStructs.Haplotype` for coordinates selected by a `BnGStructs.LocusSet`.

The coordinates extracted are `all_positions[locus_set.loci]`.

# Arguments
- `result::FisherWrightResult`: Structured simulation result.
- `locus_set::BnGStructs.LocusSet`: Locus set selecting indices from `all_positions`.
- `all_positions::Vector{UInt32}`: Master vector of sorted, unique 1-based coordinates.

# Keywords
- `include_fixed::Bool = false`: If `true`, includes fixed substitutions in the extraction.
- `selected::Bool = false`: If `true`, the selected alleles of a run with selection count as carried.

# Returns
- `Haplotype`: A `BnGStructs.Haplotype` for the selected chip coordinates.
"""
function to_haplotype(result::FisherWrightResult, locus_set::BnGStructs.LocusSet, all_positions::Vector{UInt32}; include_fixed::Bool=false, selected::Bool=false)
    chip_positions = all_positions[locus_set.loci]
    return to_haplotype(result, chip_positions; include_fixed=include_fixed, selected=selected)
end
