"""
    merge_sorted!(dest, v1, v2)

Write the sorted set-union of the sorted vectors `v1` and `v2` (duplicates
removed) into `dest`, and return `dest`.

`dest` is resized as needed and its previous contents are discarded. Because
the buffer is retained between calls, a `dest` that is reused across
generations stops allocating once it has grown to its steady-state size. This
is the allocation-free form of [`merge_sorted`](@ref).

Preconditions:
- `v1` and `v2` must each be sorted in nondecreasing order.
- `dest` must not alias `v1` or `v2`.

Complexity: O(length(v1) + length(v2)) time, no allocation once `dest` is large
enough.
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
    merge_sorted(v1, v2)

Return the sorted set-union of two sorted vectors `v1` and `v2` (duplicates removed).

Preconditions:
- `v1` and `v2` must each be sorted in nondecreasing order.
- Element type must support `<` and `==`.

Complexity: O(length(v1) + length(v2)) time, O(length(v1) + length(v2)) worst-case transient space.

Use [`merge_sorted!`](@ref) in hot loops to reuse a destination buffer instead
of allocating a fresh vector per call.
"""
function merge_sorted(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T}
    return merge_sorted!(Vector{T}(undef, 0), v1, v2)
end
