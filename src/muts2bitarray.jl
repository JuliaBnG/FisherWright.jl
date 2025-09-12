"""
    muts2bitarray(muts, cbp; flip=false)

Convert per-haplotype mutation position vectors `muts` (UInt32, sorted, unique
per hap) and chromosome cumulative bp end positions `cbp` (UInt32) into:

- BitArray (loci × haplotypes) of presence (1) / absence (0)
- DataFrame with columns: chr (chromosome index), pos (bp), ref, alt, frq
  (allele frequency)

Options:
- flip::Bool: if true, randomly flip reference/alternate coding (invert bits) at
  randomly chosen loci (uniform, independent).

Fixed (monomorphic) loci are removed.
"""
function muts2bitarray(
    muts::Vector{Vector{UInt32}},
    cbp::Vector{UInt32};
    flip::Bool = false,
)
    nhp = length(muts)
    # Collect all mutations once
    total = 0
    @inbounds for t in muts
        total += length(t)
    end
    all_mts = Vector{UInt32}(undef, total)
    pos = 1
    @inbounds for t in muts
        lt = length(t)
        if lt > 0
            copyto!(all_mts, pos, t, 1, lt)
            pos += lt
        end
    end
    resize!(all_mts, pos-1)
    sort!(all_mts)
    # unique! in-place
    ulen = 0
    last = UInt32(0)
    @inbounds for i in eachindex(all_mts)
        v = all_mts[i]
        if i == 1 || v != last
            ulen += 1
            all_mts[ulen] = v
            last = v
        end
    end
    resize!(all_mts, ulen)
    nlc = length(all_mts)
    # Early exit
    if nlc == 0
        return BitArray(undef, 0, nhp),
        DataFrame(
            chr = Int8[],
            pos = UInt32[],
            ref = Char[],
            alt = Char[],
            frq = Float32[],
        )
    end
    # Build index (could alternatively binary search per mut)
    idx = Dict{UInt32,Int}()
    sizehint!(idx, nlc)
    @inbounds for (i, m) in enumerate(all_mts)
        idx[m] = i
    end
    xy = BitArray(undef, nlc, nhp)
    fill!(xy, false)
    Threads.@threads for i = 1:nhp
        hap = muts[i]
        @inbounds for mut in hap
            xy[idx[mut], i] = true
        end
    end
    if flip
        mask = rand(Bool, nlc)
        # Invert selected rows
        @inbounds for r = 1:nlc
            mask[r] && (xy[r, :] = .!view(xy, r, :))
        end
    end
    # Chromosome assignment via cumulative ends (cbp assumed sorted)
    chr = Vector{Int8}(undef, nlc)
    @inbounds for i = 1:nlc
        chr[i] = Int8(searchsortedfirst(cbp, all_mts[i]))
    end
    # Alleles: choose alt != ref per locus
    bases = Vector{Char}(['A', 'C', 'G', 'T'])
    ref = rand(bases, nlc)
    alt = Vector{Char}(undef, nlc)
    @inbounds for i = 1:nlc
        r = ref[i]
        # pick from the 3 remaining
        # simple rejection (fast enough)
        a = r
        while a == r
            a = bases[rand(1:4)]
        end
        alt[i] = a
    end
    counts = vec(sum(xy, dims = 2))                # Int counts
    polym = (counts .> 0) .& (counts .< nhp)     # polymorphic mask
    if !any(polym)
        return BitArray(undef, 0, nhp),
        DataFrame(
            chr = Int8[],
            pos = UInt32[],
            ref = Char[],
            alt = Char[],
            frq = Float32[],
        )
    end
    frq = Float32.(counts) ./ nhp                # derived allele frequency
    lmp = DataFrame(
        chr = chr[polym],
        pos = all_mts[polym],
        ref = ref[polym],
        alt = alt[polym],
        frq = frq[polym],
    )
    return xy[polym, :], lmp
end
