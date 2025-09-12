"""
    function merge_sorted(v1::Vector{UInt32}, v2::Vector{UInt32})
Merge two sorted vectors into one sorted vector without duplicates.
"""
function merge_sorted(v1::Vector{UInt32}, v2::Vector{UInt32})
    i, j, k = 1, 1, 1
    merged = Vector{UInt32}(undef, length(v1) + length(v2))

    # Merge while both vectors have elements
    while i <= length(v1) && j <= length(v2)
        if v1[i] < v2[j]
            if k == 1 || merged[k-1] != v1[i] # Avoid duplicates
                merged[k] = v1[i]
                k += 1
            end
            i += 1
        elseif v1[i] > v2[j]
            if k == 1 || merged[k-1] != v2[j] # Avoid duplicates
                merged[k] = v2[j]
                k += 1
            end
            j += 1
        else
            # Handle duplicates: add only one of the duplicate elements
            if k == 1 || merged[k-1] != v1[i]
                merged[k] = v1[i]
                k += 1
            end
            i += 1
            j += 1
        end
    end

    # Add remaining elements from v1
    while i <= length(v1)
        if k == 1 || merged[k-1] != v1[i] # Avoid duplicates
            merged[k] = v1[i]
            k += 1
        end
        i += 1
    end

    # Add remaining elements from v2
    while j <= length(v2)
        if k == 1 || merged[k-1] != v2[j] # Avoid duplicates
            merged[k] = v2[j]
            k += 1
        end
        j += 1
    end

    # Resize the merged vector to remove unused space
    return merged[1:(k-1)]
end
