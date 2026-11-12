"""
    FisherWrightResult(active_haplotypes, chromosome_ends, substitutions, census, males,
                       selected_haplotypes, selected_substitutions,
                       mean_fitness, mean_log_fitness, var_log_fitness,
                       mean_genotypic, var_genotypic, mean_phenotypic, var_phenotypic,
                       paternal_pairs, maternal_pairs, first_id, pedigree)
    FisherWrightResult(active_haplotypes, chromosome_ends, substitutions;
                       selected_haplotypes, selected_substitutions)

Structured container holding the output of a Fisher-Wright simulation, keeping
active haplotypes, chromosome breakpoints, and fixed substitutions together
with the census record of the run.

# Fields
- `active_haplotypes::Vector{Vector{UInt32}}`: Active haplotype mutation positions (two per individual of the
  last generation). Each inner vector contains sorted, unique 1-based base pair coordinates of segregating mutations.
- `chromosome_ends::Vector{UInt32}`: Cumulative chromosome end positions in base pairs (1-based, strictly increasing).
- `substitutions::Vector{UInt32}`: Sorted, unique positions of mutations that reached fixation across all
  haplotypes and were removed from `active_haplotypes` during simulation.
- `census::Vector{Int}`: Census sizes ``N_0, N_1, …, N_{nt}`` of the recorded generations: `census[1]` is
  the starting population and `census[g + 1]` the number of individuals born in generation `g`.
- `males::Vector{Int}`: `males[g]` is the realized number of males among the `census[g]` parents of generation `g`;
  the remaining `census[g] - males[g]` parents are females.

- `selected_haplotypes::Vector{Vector{UInt32}}`: for each haplotype, the sorted positions of the derived alleles
  it carries at the fitness loci of a run with selection ([`SiteFitness`](@ref)); empty vectors in a neutral run.
  These alleles are kept apart from the neutral mutations in `active_haplotypes`.
- `selected_substitutions::Vector{UInt32}`: sorted positions of fitness loci whose derived allele reached fixation
  and was removed from `selected_haplotypes`.

- `mean_fitness::Vector{Float64}`, `mean_log_fitness::Vector{Float64}`, `var_log_fitness::Vector{Float64}`:
  fitness of the parents of each generation, entry `g` describing the same `census[g]` adults as `males[g]`:
  the arithmetic mean of fitness, and the mean and the (sample) variance of its logarithm. Fitness is absolute,
  relative to an individual without any derived allele at the fitness loci, so substituted sites are included.
  In a neutral run these are `1`, `0` and `0`. If an adult has zero fitness (a lethal genotype), the mean of the
  logarithm is `-Inf` and its variance `NaN`.
- `mean_genotypic::Vector{Float64}`, `var_genotypic::Vector{Float64}`, `mean_phenotypic::Vector{Float64}`,
  `var_phenotypic::Vector{Float64}`: under [`TraitFitness`](@ref), the mean and (sample) variance of the genotypic
  value and of the phenotype of the parents of each generation, as for the fitness record; `NaN` for generations
  without a trait model.
- `paternal_pairs::Vector{Float64}`, `maternal_pairs::Vector{Float64}`: the realized family sizes of each
  generation, as the number of ordered pairs of its offspring that share their sire, `Σ m(m - 1)` over the sires,
  and that share their dam, `Σ n(n - 1)` over the dams; [`effective_size`](@ref) uses them. `NaN` in a result
  built by hand.
- `first_id::Int`: the ID of individual 1 of `active_haplotypes`; individual `i` (haplotypes `2i - 1` and `2i`)
  has ID `first_id + i - 1`. Individuals are numbered in birth order across the run, the starting population
  first, and a run continued through `init` continues the numbering. 1 in a result built by hand.
- `pedigree::DataFrame`: with the `pedigree = k` keyword of [`fisher_wright`](@ref), the matings of the last `k`
  generations, one row per individual: `id`, `sire` and `dam` (IDs; 0 for the founders of the record, which are
  the parents of the oldest recorded generation, or the starting population), `sex` (`true` for male, the sex
  the individual had when its generation reproduced; `missing` for the last generation) and `generation` (of
  birth, numbered as in `census`). Empty without the keyword.

A run continued from a `FisherWrightResult` through the `init` keyword of [`fisher_wright`](@ref) appends its
generations to the record of that result, so generations are numbered from the first call; started from plain
haplotypes, the record begins anew. The three-argument constructor builds a result without a record, with `census`
holding the number of individuals only and `males` empty; its keywords give the selected alleles of a starting
population, for example one copy of a new allele.

# See also
[`fisher_wright`](@ref), [`effective_size`](@ref), [`to_haplotype`](@ref)
"""
struct FisherWrightResult
    active_haplotypes::Vector{Vector{UInt32}}
    chromosome_ends::Vector{UInt32}
    substitutions::Vector{UInt32}
    census::Vector{Int}
    males::Vector{Int}
    selected_haplotypes::Vector{Vector{UInt32}}
    selected_substitutions::Vector{UInt32}
    mean_fitness::Vector{Float64}
    mean_log_fitness::Vector{Float64}
    var_log_fitness::Vector{Float64}
    mean_genotypic::Vector{Float64}
    var_genotypic::Vector{Float64}
    mean_phenotypic::Vector{Float64}
    var_phenotypic::Vector{Float64}
    paternal_pairs::Vector{Float64}
    maternal_pairs::Vector{Float64}
    first_id::Int
    pedigree::DataFrame
end

_empty_pedigree() = DataFrame(id=Int[], sire=Int[], dam=Int[], sex=Union{Missing,Bool}[], generation=Int[])

# Without a pedigree: individual i has ID i, and nothing is recorded.
FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions, mean_fitness, mean_log_fitness, var_log_fitness,
    mean_genotypic, var_genotypic, mean_phenotypic, var_phenotypic, paternal_pairs, maternal_pairs,
) = FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions, mean_fitness, mean_log_fitness, var_log_fitness,
    mean_genotypic, var_genotypic, mean_phenotypic, var_phenotypic, paternal_pairs, maternal_pairs,
    1, _empty_pedigree())

# Without a record of family sizes: NaN for every recorded generation.
FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions, mean_fitness, mean_log_fitness, var_log_fitness,
    mean_genotypic, var_genotypic, mean_phenotypic, var_phenotypic,
) = FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions, mean_fitness, mean_log_fitness, var_log_fitness,
    mean_genotypic, var_genotypic, mean_phenotypic, var_phenotypic,
    fill(NaN, length(males)), fill(NaN, length(males)))

# Without a trait record: NaN for every recorded generation.
FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions, mean_fitness, mean_log_fitness, var_log_fitness,
) = FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions, mean_fitness, mean_log_fitness, var_log_fitness,
    (fill(NaN, length(males)) for _ = 1:4)...)

# Without a fitness record: that of a neutral population, one entry per
# recorded generation.
FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions,
) = FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, census, males,
    selected_haplotypes, selected_substitutions,
    ones(length(males)), zeros(length(males)), zeros(length(males)))

FisherWrightResult(active_haplotypes, chromosome_ends, substitutions, census, males) =
    FisherWrightResult(active_haplotypes, chromosome_ends, substitutions, census, males,
        [UInt32[] for _ in active_haplotypes], UInt32[])

FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions;
    selected_haplotypes=[UInt32[] for _ in active_haplotypes],
    selected_substitutions=UInt32[],
) = FisherWrightResult(
    active_haplotypes, chromosome_ends, substitutions, [length(active_haplotypes) ÷ 2], Int[],
    selected_haplotypes, selected_substitutions)

"""
    effective_size(result::FisherWrightResult, generations = eachindex(result.males); from = :families) -> Float64

Effective population size realized in the generations of `result` that are given (one generation or a collection
of them; all generations of the run by default), from its record.

With `from = :families` (the default), from the realized family sizes. Two genes in different offspring of
generation `g` derive from the same parent with probability

``a_g = \\frac{\\sum_s m_s(m_s - 1) + \\sum_d n_d(n_d - 1)}{4 N_g (N_g - 1)}``,

where `m_s` and `n_d` are the numbers of offspring of sire `s` and dam `d` (`result.paternal_pairs[g]` and
`result.maternal_pairs[g]`) and `N_g = result.census[g + 1]`. The effective size of the generation is `1/a_g`, and
over several generations it is the harmonic mean, `1/ā`. This holds for every mating scheme
([`RandomMating`](@ref), [`FamilySize`](@ref), [`Hierarchical`](@ref)) and includes the variance of family size
that selection adds within a generation; it does not include the cumulative effect of heritable fitness over
generations (Robertson 1961). It is the value realized by the pedigree of the run, not an estimate from the
genetic data, and it varies from run to run as family sizes do.

With `from = :sexes`, from the recorded numbers of male and female parents only, as in version 0.5 and 0.6:
`a_g = 1/(4N_m) + 1/(4N_f)`, the expectation under random mating, which gives Wright's (1931)
`4N_mN_f/(N_m + N_f)`. It does not apply to other mating schemes or to selection.

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
function effective_size(result::FisherWrightResult, generations=eachindex(result.males); from::Symbol=:families)
    isempty(generations) && throw(ArgumentError("no generation to average over"))
    from in (:families, :sexes) || throw(ArgumentError("from must be :families or :sexes"))
    a = 0.0
    for g in generations
        if from === :sexes
            m = result.males[g]
            f = result.census[g] - m
            a += 1 / 4m + 1 / 4f
        else
            pp, mp = result.paternal_pairs[g], result.maternal_pairs[g]
            (isnan(pp) || isnan(mp)) &&
                throw(ArgumentError("the family sizes of generation $g are not recorded; use from = :sexes"))
            n = result.census[g+1]
            a += (pp + mp) / (4n * (n - 1))
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
    return FisherWrightResult(cleaned, result.chromosome_ends, substitutions, result.census, result.males,
        result.selected_haplotypes, result.selected_substitutions,
        result.mean_fitness, result.mean_log_fitness, result.var_log_fitness,
        result.mean_genotypic, result.var_genotypic, result.mean_phenotypic, result.var_phenotypic,
        result.paternal_pairs, result.maternal_pairs, result.first_id, result.pedigree)
end


"""
    to_haplotype(result::FisherWrightResult; include_fixed::Bool = false, selected::Bool = false) -> Tuple{Haplotype, DataFrame}

Convert a structured Fisher-Wright simulation result into a dense `BnGStructs.Haplotype`
matrix and an accompanying linkage map `DataFrame`.

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

# Returns
- `Haplotype`: A `BnGStructs.Haplotype` wrapping a `BitMatrix` of shape `(nlc, nhp)` (loci in rows, haplotypes in columns).
- `DataFrame`: Marker map with columns:
  - `chr::Int8`: 1-based chromosome identifier.
  - `pos::UInt32`: Genome-wide base pair coordinate.
  - `ref::Char`: Reference nucleotide ('A', 'C', 'G', or 'T').
  - `alt::Char`: Alternate nucleotide (distinct from `ref`).
  - `frq::Float32`: Alternate allele frequency in the population.

# Errors
- Throws an `ErrorException` if `include_fixed = false` and no polymorphic loci exist across the haplotypes.

# Examples
```julia
using FisherWright

res = fisher_wright(20, 10, [100_000, 100_000], 1.0; result = true)
hap, loci = to_haplotype(res)
```
"""
function to_haplotype(result::FisherWrightResult; include_fixed::Bool=false, selected::Bool=false)
    muts = _export_haplotypes(result, include_fixed, selected)
    xy, loci = muts2bitarray(muts, result.chromosome_ends; include_fixed=include_fixed)
    isempty(xy) && error("No polymorphic loci available for dense export")
    return Haplotype(xy), loci
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
