"""
    FisherWrightResult

Structured Fisher-Wright simulation output that keeps the active haplotypes and
chromosome breakpoints together without changing the default tuple-based API.
"""
struct FisherWrightResult
    active_haplotypes::Vector{Vector{UInt32}}
    chromosome_ends::Vector{UInt32}
    substitutions::Vector{UInt32}
end

function _fixed_mutations(haplotypes::Vector{Vector{UInt32}})
    isempty(haplotypes) && return UInt32[]
    fixed = copy(haplotypes[1])
    for haplotype in haplotypes[2:end]
        fixed = intersect(fixed, haplotype)
        isempty(fixed) && break
    end
    sort!(unique!(fixed))
    return fixed
end

function _drop_mutations(haplotypes::Vector{Vector{UInt32}}, mutations::Vector{UInt32})
    isempty(mutations) && return haplotypes
    cleaned = Vector{Vector{UInt32}}(undef, length(haplotypes))
    for (i, haplotype) in pairs(haplotypes)
        cleaned[i] = isempty(haplotype) ? UInt32[] : setdiff(haplotype, mutations)
    end
    return cleaned
end

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

function _with_fixed(result::FisherWrightResult)
    isempty(result.substitutions) && return result.active_haplotypes
    fixed = sort!(unique!(copy(result.substitutions)))
    merged = Vector{Vector{UInt32}}(undef, length(result.active_haplotypes))
    for (i, haplotype) in pairs(result.active_haplotypes)
        merged[i] = merge_sorted(haplotype, fixed)
    end
    return merged
end

function _extract_fixed(result::FisherWrightResult)
    cleaned, substitutions = _fixation_step(result.active_haplotypes, result.substitutions)
    return FisherWrightResult(cleaned, result.chromosome_ends, substitutions)
end

"""
    to_haplotype(result::FisherWrightResult)

Convert a structured Fisher-Wright result to a dense `BnGStructs.Haplotype`
and the associated linkage map `DataFrame`.
"""
function to_haplotype(result::FisherWrightResult; include_fixed::Bool=false)
    muts = include_fixed ? _with_fixed(result) : result.active_haplotypes
    xy, loci = muts2bitarray(muts, result.chromosome_ends; include_fixed=include_fixed)
    isempty(xy) && error("No polymorphic loci available for dense export")
    return Haplotype(xy), loci
end