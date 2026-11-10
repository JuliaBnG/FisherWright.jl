"""
    merge_sorted!(dest::Vector{T}, v1::AbstractVector{T}, v2::AbstractVector{T}) where {T} -> Vector{T}

Compute the sorted set-union of two sorted vectors `v1` and `v2` (with duplicate values removed),
writing the result into `dest` and returning `dest`.

`dest` is resized to hold the output, discarding any previous contents. By reusing the same `dest`
buffer across generations or iterations, allocation ceases once `dest` reaches steady-state capacity.
This is the allocation-free form of [`merge_sorted`](@ref).

# Arguments
- `dest::Vector{T}`: Destination buffer. Must not alias `v1` or `v2`.
- `v1::AbstractVector{T}`: First sorted input vector.
- `v2::AbstractVector{T}`: Second sorted input vector.

# Returns
- `Vector{T}`: The modified `dest` vector containing unique elements from `v1 ∪ v2` in ascending order.

# Preconditions
- `v1` and `v2` must each be sorted in non-decreasing order.
- `dest` must not alias `v1` or `v2`.

# Complexity
- Time: ``O(|\\text{v1}| + |\\text{v2}|)``.
- Space: zero additional allocations once `dest` has sufficient capacity.

# Errors
- Throws an `ArgumentError` if `dest === v1` or `dest === v2`.

# Examples
```julia
using FisherWright

buf = UInt32[]
FisherWright.merge_sorted!(buf, UInt32[1, 3, 5], UInt32[2, 3, 6])
buf == UInt32[1, 2, 3, 5, 6]
```
"""
function merge_sorted!(
    dest::Vector{T},
    v1::AbstractVector{T},
    v2::AbstractVector{T},
) where {T}
    (dest === v1 || dest === v2) &&
        throw(ArgumentError("merge_sorted! destination must not alias its inputs"))
    n1, n2 = length(v1), length(v2)
    resize!(dest, n1 + n2)
    i = j = k = 1
    has_last = false
    last = zero(T)  # ignored until has_last = true

    @inbounds while i <= n1 && j <= n2
        a = v1[i]
        b = v2[j]
        if a < b
            if !has_last || last != a
                dest[k] = a
                last = a
                has_last = true
                k += 1
            end
            i += 1
        elseif a > b
            if !has_last || last != b
                dest[k] = b
                last = b
                has_last = true
                k += 1
            end
            j += 1
        else
            if !has_last || last != a
                dest[k] = a
                last = a
                has_last = true
                k += 1
            end
            i += 1
            j += 1
        end
    end

    @inbounds while i <= n1
        a = v1[i]
        if !has_last || last != a
            dest[k] = a
            last = a
            has_last = true
            k += 1
        end
        i += 1
    end
    @inbounds while j <= n2
        b = v2[j]
        if !has_last || last != b
            dest[k] = b
            last = b
            has_last = true
            k += 1
        end
        j += 1
    end

    resize!(dest, k - 1)
    return dest
end

"""
    merge_sorted(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T} -> Vector{T}

Return the sorted set-union of two sorted vectors `v1` and `v2` as a newly allocated vector,
with duplicate values removed.

# Arguments
- `v1::AbstractVector{T}`: First sorted input vector.
- `v2::AbstractVector{T}`: Second sorted input vector.

# Returns
- `Vector{T}`: Fresh vector containing unique elements from `v1 ∪ v2` in ascending order.

# Preconditions
- `v1` and `v2` must each be sorted in non-decreasing order.
- Element type `T` must support `<` and `==`.

# Complexity
- Time: ``O(|\\text{v1}| + |\\text{v2}|)``.
- Space: allocates a new vector of length up to `length(v1) + length(v2)`.

# See also
[`merge_sorted!`](@ref)

# Examples
```julia
using FisherWright

FisherWright.merge_sorted([1, 4, 7], [2, 4, 6]) == [1, 2, 4, 6, 7]
```
"""
function merge_sorted(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T}
    return merge_sorted!(Vector{T}(undef, 0), v1, v2)
end

"""
    _isearchfirst(v::Vector{UInt32}, x::UInt32, lo::Int, hi::Int) -> Int

First index `k` in `lo:hi` with `v[k] >= x`, or `hi + 1` if there is none: the same as
`searchsortedfirst(v, x, lo, hi, Base.Order.Forward)` for sorted `v`.

Haplotype positions are close to uniform along the genome, so the index is first
guessed by linear interpolation between `v[lo]` and `v[hi]` and the bracket narrowed
with a second probe a few elements away. Up to four such rounds are made while the
bracket is wider than 32 elements, then a binary search finishes. On haplotypes too
large for the cache this needs a few memory reads instead of about ``\\log_2 n``,
each of which is a likely cache miss.
"""
function _isearchfirst(v::Vector{UInt32}, x::UInt32, lo::Int, hi::Int)
    @inbounds begin
        lo > hi && return lo
        v[lo] >= x && return lo
        v[hi] < x && return hi + 1
        # Invariant: v[lo] < x <= v[hi], so the answer is in lo+1:hi.
        rounds = 0
        while hi - lo > 32 && rounds < 4
            a, b = Int(v[lo]), Int(v[hi])
            g = lo + div((Int(x) - a) * (hi - lo), b - a)
            g = clamp(g, lo + 1, hi - 1)
            if v[g] < x
                lo = g
                t = min(g + 8, hi)
                v[t] >= x ? (hi = t) : (lo = t)
            else
                hi = g
                t = max(g - 8, lo)
                v[t] < x ? (lo = t) : (hi = t)
            end
            rounds += 1
        end
        return searchsortedfirst(v, x, lo + 1, hi, Base.Order.Forward)
    end
end
_isearchfirst(v::AbstractVector, x, lo::Int, hi::Int) =
    searchsortedfirst(v, x, lo, hi, Base.Order.Forward)

"""
    _insert_sorted!(dest::Vector{T}, v::Vector{T}, ins::AbstractVector{T}) where {T} -> Vector{T}

Sorted set-union of `v` and a short `ins`, written into `dest`: the same result as
[`merge_sorted!`](@ref)`(dest, v, ins)` when `v` is sorted and duplicate-free.

Each element of `ins` is located in `v` by [`_isearchfirst`](@ref), and the runs of `v` between
insertion points are block-copied, so the cost is ``O(k \\log n)`` comparisons plus one
memory copy of `v` (``k = |\\text{ins}|``, ``n = |\\text{v}|``), instead of a scalar
two-pointer pass over all of `v`. This is the per-haplotype mutation step of
[`fisher_wright`](@ref), where ``k \\approx 10`` and ``n`` is in the tens of thousands.

# Preconditions
- `v` is sorted and duplicate-free; `ins` is sorted (it may contain duplicates or
  elements of `v`, which are dropped).
- `dest` does not alias `v` or `ins`.
"""
function _insert_sorted!(dest::Vector{T}, v::Vector{T}, ins::AbstractVector{T}) where {T}
    (dest === v || dest === ins) &&
        throw(ArgumentError("_insert_sorted! destination must not alias its inputs"))
    n = length(v)
    resize!(dest, n + length(ins))
    k = 0 # last written index of dest
    i = 1 # next unread index of v
    @inbounds for x in ins
        p = _isearchfirst(v, x, i, n)
        c = p - i
        if c > 0
            copyto!(dest, k + 1, v, i, c)
            k += c
        end
        i = p
        # v[i:p-1] < x, so dest[k] == x only for a repeat within ins.
        ((p <= n && v[p] == x) || (k > 0 && dest[k] == x)) && continue
        k += 1
        dest[k] = x
    end
    c = n - i + 1
    if c > 0
        copyto!(dest, k + 1, v, i, c)
        k += c
    end
    resize!(dest, k)
    return dest
end
