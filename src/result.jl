"""
    FisherWrightResult(active_haplotypes, chromosome_ends, substitutions)

Structured container holding the output of a Fisher-Wright simulation, keeping
active haplotypes, chromosome breakpoints, and fixed substitutions together.

# Fields
- `active_haplotypes::Vector{Vector{UInt32}}`: Active haplotype mutation positions (`2 * ne` vectors).
  Each inner vector contains sorted, unique 1-based base pair coordinates of segregating mutations.
- `chromosome_ends::Vector{UInt32}`: Cumulative chromosome end positions in base pairs (1-based, strictly increasing).
- `substitutions::Vector{UInt32}`: Sorted, unique positions of mutations that reached fixation across all
  haplotypes and were removed from `active_haplotypes` during simulation.

# See also
[`fisher_wright`](@ref), [`to_haplotype`](@ref)
"""
struct FisherWrightResult
    active_haplotypes::Vector{Vector{UInt32}}
    chromosome_ends::Vector{UInt32}
    substitutions::Vector{UInt32}
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
    _fixed_mutations(haplotypes::Vector{Vector{UInt32}}) -> Vector{UInt32}

Identify and return sorted mutation positions carried by every haplotype in `haplotypes`.
Each inner vector in `haplotypes` must be sorted and unique.
"""
function _fixed_mutations(haplotypes::Vector{Vector{UInt32}})
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
function _drop_mutations!(haplotypes::Vector{Vector{UInt32}}, mutations::Vector{UInt32})
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
    haplotypes::Vector{Vector{UInt32}},
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

"""
    _extract_fixed(result::FisherWrightResult) -> FisherWrightResult

Scan `result.active_haplotypes` for fixed mutations, remove them from active haplotypes,
and return a new `FisherWrightResult` with updated haplotypes and accumulated substitutions.
"""
function _extract_fixed(result::FisherWrightResult)
    cleaned, substitutions = _fixation_step(result.active_haplotypes, result.substitutions)
    return FisherWrightResult(cleaned, result.chromosome_ends, substitutions)
end


"""
    to_haplotype(result::FisherWrightResult; include_fixed::Bool = false) -> Tuple{Haplotype, DataFrame}

Convert a structured Fisher-Wright simulation result into a dense `BnGStructs.Haplotype`
matrix and an accompanying linkage map `DataFrame`.

# Arguments
- `result::FisherWrightResult`: Structured simulation result containing active haplotypes, chromosome ends, and substitutions.

# Keywords
- `include_fixed::Bool = false`: If `true`, merges fixed `substitutions` back into each haplotype before matrix construction.
  If `false`, only segregating (polymorphic) loci are exported.

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
function to_haplotype(result::FisherWrightResult; include_fixed::Bool=false)
    muts = include_fixed ? _with_fixed(result) : result.active_haplotypes
    xy, loci = muts2bitarray(muts, result.chromosome_ends; include_fixed=include_fixed)
    isempty(xy) && error("No polymorphic loci available for dense export")
    return Haplotype(xy), loci
end

"""
    to_haplotype(result::FisherWrightResult, chip_positions::Vector{UInt32}; include_fixed::Bool = false) -> Haplotype

Directly extract a dense `BnGStructs.Haplotype` at specified coordinate positions `chip_positions`
without materializing other simulated loci.

# Arguments
- `result::FisherWrightResult`: Structured simulation result.
- `chip_positions::Vector{UInt32}`: Sorted, unique 1-based base pair coordinates to extract.

# Keywords
- `include_fixed::Bool = false`: If `true`, includes fixed substitutions from `result` in the extraction.

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
function to_haplotype(result::FisherWrightResult, chip_positions::Vector{UInt32}; include_fixed::Bool=false)
    muts = include_fixed ? _with_fixed(result) : result.active_haplotypes
    xy = extract_chip_bitarray(muts, chip_positions)
    return Haplotype(xy)
end

"""
    to_haplotype(result::FisherWrightResult, locus_set::BnGStructs.LocusSet, all_positions::Vector{UInt32}; include_fixed::Bool = false) -> Haplotype

Directly extract a dense `BnGStructs.Haplotype` for coordinates selected by a `BnGStructs.LocusSet`.

The coordinates extracted are `all_positions[locus_set.loci]`.

# Arguments
- `result::FisherWrightResult`: Structured simulation result.
- `locus_set::BnGStructs.LocusSet`: Locus set selecting indices from `all_positions`.
- `all_positions::Vector{UInt32}`: Master vector of sorted, unique 1-based coordinates.

# Keywords
- `include_fixed::Bool = false`: If `true`, includes fixed substitutions in the extraction.

# Returns
- `Haplotype`: A `BnGStructs.Haplotype` for the selected chip coordinates.
"""
function to_haplotype(result::FisherWrightResult, locus_set::BnGStructs.LocusSet, all_positions::Vector{UInt32}; include_fixed::Bool=false)
    chip_positions = all_positions[locus_set.loci]
    return to_haplotype(result, chip_positions; include_fixed=include_fixed)
end
