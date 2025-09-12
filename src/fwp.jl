"""
    fisher_wright( # Fisher-Wright population
        ne::T1,
        nt::T2,
        chr::Vector{T3}, # chromosome lengths in bp,
        mr::Float64; # mutation rate per 1e8 base pair per meiosis
        M = 1e8,     # base pair per Morgan
    ) where {T1<:Integer, T2<:Integer, T3<:Integer}
Simulate a Fisher-Wright population of size `ne` by `nt` generations of random
mating. The chromosome lengths are given in `chr` (in Morgan). The mutation rate
is `mr` per 10⁸ bp per meiosis. The base pair per Morgan is `M` with default
value 10⁸. The mutations are stored in a `Set{UInt32}`. The maximum mutation
place in bp is hence 2³², which is 54.95 × 10⁸ bp. This is sufficient for most
of the genomes. An error occurred if the total length of the genome is larger
than this value.

Note:
1. parameter recombination rate was removed as it is confusing. Crossover rate
   is a more appropriate term. Yet, Crossovers is also defined in M, or,
   bp/Morgan, i.e., the chromosome length in M is also the expected number of
   crossovers of this chromosome per meiosis. The less bp per Morgan, the more
   crossovers.
2. mutation rate is asked for per 1e8 base pair, the typical number of bp per M,
   per meiosis. This is to avoid tedious numbers like 1e-8.
3. this function is **depricated**, use `fisher_wright(nid, chr, mr; M=1e8)`
   instead, where the `nid` is a vector of the number of individuals in each
   generation. This can simulate, e.g., bottlenecks, in history.

# The algorithm

1. Create two vectors of containers to store mutations in parents and offspring.
2. Generation of `M⋅m` new mutations ∈ ``[1, ~3×10⁹]`` for each haplotype
3. Insert these mutations into the parent generations
4. Randomly sample `Nₑ` pairs of ID in parent generation as sires and dams to
   offspring.
5. Splice sire and dam haplotypes into offspring's paternal haplotype
6. Swap parent and offspring storage
7. Repeat 1-6 for `nt` times.
"""
function fisher_wright(
    ne::T1,
    nt::T2,
    chr::Vector{T3},
    mr::Float64;
    M = 1e8,
) where {T1<:Integer,T2<:Integer,T3<:Integer}
    if !(ne > 1 && nt > 0 && all(chr .> 0) && 0.01 < mr < 20.0)
        throw(ArgumentError("Invalid parameter(s)"))
    end

    tg = sum(chr)
    tg < 2^32 || error("Total genome length must be < 2^32 bp for UInt32 storage")
    cbp = UInt32.(cumsum(chr))
    tbp = cbp[end]

    p_mut = Poisson(tg / M * mr)
    p_xo = Poisson.(chr / M)

    nh = 2 * ne
    prt = [Vector{UInt32}() for _ = 1:nh]
    off = [Vector{UInt32}() for _ = 1:nh]

    @info "Fisher-Wright population simulation start" ne nt total_bp=tg

    for g = 1:nt
        if g % 100 == 0
            print(
                '\r',
                ' '^8,
                "Generation $g / $nt, mean muts/haps: ",
                round(mean(length.(prt)); digits = 2),
            )
        end
        # Mutations (threaded)
        Threads.@threads for i = 1:nh
            nm = rand(p_mut)
            if nm > 0
                # Reuse a local buffer (allocate only when needed)
                newm = Vector{UInt32}(undef, nm)
                @inbounds for k = 1:nm
                    newm[k] = rand(UInt32(1):UInt32(tbp))
                end
                sort!(newm)
                prt[i] = merge_sorted(prt[i], newm)
            end
        end

        pm, _ = random_mate(ne, ne)

        Threads.@threads for i = 1:ne
            s = pm[i, 1];
            d = pm[i, 2]
            # Ensure recombine overwrites (or empty! target first if needed)
            rss = cobp(cbp, p_xo)
            recombine(prt[2s-1], prt[2s], off[2i-1], rss)
            rsd = cobp(cbp, p_xo)
            recombine(prt[2d-1], prt[2d], off[2i], rsd)
        end
        prt, off = off, prt
    end
    println()
    return prt, cbp
end
