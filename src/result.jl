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

"""
    _intersect_sorted!(dest, cand, hap)

Write `cand ∩ hap` into `dest` and return it. Both inputs must be sorted and
unique; the result is sorted and unique.

Two strategies are used. When `cand` is much shorter than `hap` each candidate
is located by binary search, which costs O(|cand| log |hap|); otherwise a
linear merge costs O(|cand| + |hap|). Neither allocates beyond growing `dest`.
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
    _remove_sorted!(hap, rm)

Delete every position of `rm` from `hap` in place and return `hap`. Both must
be sorted and unique. Runs in O(|hap| + |rm|) with no allocation.
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
    _fixed_mutations(haplotypes)

Return the sorted positions carried by every haplotype. Haplotypes must be
sorted and unique.
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
    _drop_mutations!(haplotypes, mutations)

Remove `mutations` from every haplotype in place and return `haplotypes`.
"""
function _drop_mutations!(haplotypes::Vector{Vector{UInt32}}, mutations::Vector{UInt32})
    isempty(mutations) && return haplotypes
    Threads.@threads for i in eachindex(haplotypes)
        _remove_sorted!(haplotypes[i], mutations)
    end
    return haplotypes
end

"""
    _drop_mutations(haplotypes, mutations)

Non-mutating form of [`_drop_mutations!`](@ref): return fresh haplotypes with
`mutations` removed, leaving the input untouched.
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
    _fixation_step!(haplotypes, substitutions)

Move newly fixed positions out of `haplotypes` (in place) and return the
updated substitution list.
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
    _fixation_step(haplotypes, substitutions)

Non-mutating form of [`_fixation_step!`](@ref), returning
`(cleaned_haplotypes, updated_substitutions)`.
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
