"""
    merge_sorted(v1, v2)

Return the sorted set-union of two sorted vectors `v1` and `v2` (duplicates removed).

Preconditions:
- `v1` and `v2` must each be sorted in nondecreasing order.
- Element type must support `<` and `==`.

Complexity: O(length(v1) + length(v2)) time, O(length(v1) + length(v2)) worst-case transient space.
"""
function merge_sorted(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T}
    n1, n2 = length(v1), length(v2)
    merged = Vector{T}(undef, n1 + n2)
    i = j = k = 1
    has_last = false
    last = zero(T)  # ignored until has_last = true

    @inbounds while i <= n1 && j <= n2
        a = v1[i];
        b = v2[j]
        if a < b
            if !has_last || last != a
                merged[k] = a;
                last = a;
                has_last = true;
                k += 1
            end
            i += 1
        elseif a > b
            if !has_last || last != b
                merged[k] = b;
                last = b;
                has_last = true;
                k += 1
            end
            j += 1
        else
            if !has_last || last != a
                merged[k] = a;
                last = a;
                has_last = true;
                k += 1
            end
            i += 1;
            j += 1
        end
    end

    @inbounds while i <= n1
        a = v1[i]
        if !has_last || last != a
            merged[k] = a;
            last = a;
            has_last = true;
            k += 1
        end
        i += 1
    end
    @inbounds while j <= n2
        b = v2[j]
        if !has_last || last != b
            merged[k] = b;
            last = b;
            has_last = true;
            k += 1
        end
        j += 1
    end

    resize!(merged, k - 1)
    return merged
end
