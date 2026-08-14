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
    flip::Bool=false,
    include_fixed::Bool=false,
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
            chr=Int8[],
            pos=UInt32[],
            ref=Char[],
            alt=Char[],
            frq=Float32[],
        )
    end

    xy = falses(nlc, nhp)
    # Fast parallel two-pointer / binary-search scan: both hap and all_mts are sorted!
    # Completely eliminates the Dict/hash table bottleneck
    Threads.@threads for i = 1:nhp
        hap = muts[i]
        nh = length(hap)
        if nh > 0
            ptr_all = 1
            ptr_hap = 1
            @inbounds while ptr_hap <= nh && ptr_all <= nlc
                m_hap = hap[ptr_hap]
                m_all = all_mts[ptr_all]
                if m_hap == m_all
                    xy[ptr_all, i] = true
                    ptr_hap += 1
                    ptr_all += 1
                elseif m_hap > m_all
                    ptr_all = searchsortedfirst(all_mts, m_hap, ptr_all + 1, nlc, Base.Order.Forward)
                else
                    ptr_hap += 1
                end
            end
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
        a = r
        while a == r
            a = bases[rand(1:4)]
        end
        alt[i] = a
    end
    counts = vec(sum(xy, dims=2))
    frq = Float32.(counts) ./ nhp
    if include_fixed
        lmp = DataFrame(chr=chr, pos=all_mts, ref=ref, alt=alt, frq=frq)
        return xy, lmp
    end
    polym = (counts .> 0) .& (counts .< nhp)
    if !any(polym)
        return BitArray(undef, 0, nhp),
        DataFrame(
            chr=Int8[],
            pos=UInt32[],
            ref=Char[],
            alt=Char[],
            frq=Float32[],
        )
    end
    lmp = DataFrame(
        chr=chr[polym],
        pos=all_mts[polym],
        ref=ref[polym],
        alt=alt[polym],
        frq=frq[polym],
    )
    return xy[polym, :], lmp
end

"""
    extract_chip_bitarray(muts::Vector{Vector{UInt32}}, chip_positions::Vector{UInt32}) -> BitArray

Directly extract a dense `BitArray` (k x nhp) at sorted, unique genomic
`chip_positions` from sparse per-haplotype `muts` vectors. Positions absent
from every haplotype produce all-false rows.
"""
function extract_chip_bitarray(muts::Vector{Vector{UInt32}}, chip_positions::Vector{UInt32})
    issorted(chip_positions) && allunique(chip_positions) ||
        throw(ArgumentError("chip_positions must be sorted and unique"))
    k = length(chip_positions)
    nhp = length(muts)
    xy = falses(k, nhp)

    Threads.@threads for i = 1:nhp
        hap = muts[i]
        nh = length(hap)
        if nh > 0
            ptr_chip = 1
            ptr_hap = 1
            @inbounds while ptr_hap <= nh && ptr_chip <= k
                m_hap = hap[ptr_hap]
                m_chip = chip_positions[ptr_chip]
                if m_hap == m_chip
                    xy[ptr_chip, i] = true
                    ptr_hap += 1
                    ptr_chip += 1
                elseif m_hap > m_chip
                    ptr_chip = searchsortedfirst(chip_positions, m_hap, ptr_chip + 1, k, Base.Order.Forward)
                else
                    ptr_hap = searchsortedfirst(hap, m_chip, ptr_hap + 1, nh, Base.Order.Forward)
                end
            end
        end
    end
    return xy
end
