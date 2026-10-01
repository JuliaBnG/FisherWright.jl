"""
    _atomic_set_bit!(chunks_ptr::Ptr{UInt64}, bit_idx::Int)

Atomically set the bit at 0-based linear index `bit_idx` to 1 in a `BitArray` chunks buffer.
Prevents data races when multiple threads write to adjacent columns that share a 64-bit chunk.
"""
@inline function _atomic_set_bit!(chunks_ptr::Ptr{UInt64}, bit_idx::Int)
    c_idx = bit_idx >>> 6
    bit = bit_idx & 63
    mask = UInt64(1) << bit
    Core.Intrinsics.atomic_pointermodify(chunks_ptr + c_idx * sizeof(UInt64), |, mask, :monotonic)
    return nothing
end

"""
    muts2bitarray(
        muts::Vector{Vector{UInt32}},
        cbp::Vector{UInt32};
        flip::Bool = false,
        include_fixed::Bool = false,
    ) -> Tuple{BitMatrix, DataFrame}

Convert sparse per-haplotype mutation coordinate vectors into a dense boolean presence/absence
matrix (`BitMatrix`) and an accompanying marker linkage map `DataFrame`.

# Arguments
- `muts::Vector{Vector{UInt32}}`: Vectors of 1-based mutation positions per haplotype (length `nhp = 2 * ne`).
  Each inner vector must be sorted and unique.
- `cbp::Vector{UInt32}`: Cumulative chromosome end positions in base pairs (1-based, strictly increasing).

# Keywords
- `flip::Bool = false`: If `true`, randomly inverts the reference/alternate allele assignment (flips bits)
  independently at each locus with probability 0.5.
- `include_fixed::Bool = false`: If `true`, includes monomorphic / fixed loci (mutations present in all haplotypes
  or absent from all haplotypes) in the output. If `false` (default), monomorphic loci are filtered out, returning
  only polymorphic loci.

# Returns
- `BitMatrix`: Dense matrix of dimensions `(nlc, nhp)` where rows correspond to loci and columns correspond to haplotypes.
  `true` (1) indicates the alternate allele, and `false` (0) indicates the reference allele.
- `DataFrame`: Marker map with columns:
  - `chr::Int8`: 1-based chromosome index determined by `cbp`.
  - `pos::UInt32`: Genome-wide base pair position.
  - `ref::Char`: Reference nucleotide ('A', 'C', 'G', or 'T').
  - `alt::Char`: Alternate nucleotide (sampled from the 3 nucleotides distinct from `ref`).
  - `frq::Float32`: Sample alternate allele frequency (`sum(xy, dims=2) / nhp`).

# Details
- Uses parallel two-pointer scanning / binary search across haplotypes, avoiding hash table overhead.
- If no qualifying loci remain (e.g. all mutations were fixed and `include_fixed = false`), returns an empty `0 × nhp`
  `BitMatrix` and an empty `DataFrame` with typed columns.

# Examples
```julia
using FisherWright

muts = [UInt32[100, 500], UInt32[500, 900], UInt32[100, 900]]
cbp = UInt32[1000]
xy, lmp = muts2bitarray(muts, cbp)
size(xy) # (3, 3)
```
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
    if nlc < 64 || nhp < 16
        for i = 1:nhp
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
    else
        p = pointer(xy.chunks)
        GC.@preserve xy begin
            Threads.@threads for i = 1:nhp
                hap = muts[i]
                nh = length(hap)
                if nh > 0
                    col_offset = (i - 1) * nlc
                    ptr_all = 1
                    ptr_hap = 1
                    @inbounds while ptr_hap <= nh && ptr_all <= nlc
                        m_hap = hap[ptr_hap]
                        m_all = all_mts[ptr_all]
                        if m_hap == m_all
                            _atomic_set_bit!(p, col_offset + ptr_all - 1)
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
    extract_chip_bitarray(
        muts::Vector{Vector{UInt32}},
        chip_positions::Vector{UInt32},
    ) -> BitMatrix

Directly extract a dense `BitMatrix` of shape `(k, nhp)` at targeted marker coordinates `chip_positions`
from sparse per-haplotype mutation vectors `muts`.

# Arguments
- `muts::Vector{Vector{UInt32}}`: Vectors of sorted, unique 1-based mutation positions per haplotype (`nhp = length(muts)`).
- `chip_positions::Vector{UInt32}`: Target marker coordinates in base pairs (`k = length(chip_positions)`).
  Must be strictly sorted in ascending order and unique.

# Returns
- `BitMatrix`: A `k × nhp` boolean matrix where entry `(r, c)` is `true` if haplotype `c`
  carries the mutation at `chip_positions[r]`, and `false` otherwise. Positions not observed in any
  haplotype produce all-false rows.

# Errors
- Throws an `ArgumentError` if `chip_positions` is not strictly sorted and unique.

# Details
- Uses a parallel two-pointer and binary search algorithm across haplotypes (`Threads.@threads`),
  completely bypassing the materialization and filtering of non-chip loci.

# Examples
```julia
using FisherWright

muts = [UInt32[100, 500], UInt32[200, 500]]
chip = UInt32[100, 200, 300, 500]
mat = extract_chip_bitarray(muts, chip)
size(mat) # (4, 2)
```
"""
function extract_chip_bitarray(muts::Vector{Vector{UInt32}}, chip_positions::Vector{UInt32})
    issorted(chip_positions) && allunique(chip_positions) ||
        throw(ArgumentError("chip_positions must be sorted and unique"))
    k = length(chip_positions)
    nhp = length(muts)
    xy = falses(k, nhp)

    if k < 64 || nhp < 16
        for i = 1:nhp
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
    else
        p = pointer(xy.chunks)
        GC.@preserve xy begin
            Threads.@threads for i = 1:nhp
                hap = muts[i]
                nh = length(hap)
                if nh > 0
                    col_offset = (i - 1) * k
                    ptr_chip = 1
                    ptr_hap = 1
                    @inbounds while ptr_hap <= nh && ptr_chip <= k
                        m_hap = hap[ptr_hap]
                        m_chip = chip_positions[ptr_chip]
                        if m_hap == m_chip
                            _atomic_set_bit!(p, col_offset + ptr_chip - 1)
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
        end
    end
    return xy
end
