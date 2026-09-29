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
