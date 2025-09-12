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
function fisher_wright( # Fisher-Wright population
    ne::T1,
    nt::T2,
    chr::Vector{T3}, # chromosome lengths in bp,
    mr::Float64; # mutation rate per 1e8 base pair per meiosis
    M = 1e8,     # base pair per Morgan
) where {T1<:Integer,T2<:Integer,T3<:Integer}
    all([ne > 1, nt > 0, all(chr .> 0), 0.01 < mr < 20.0]) || error("Invalid parameter(s)")

    # genome parameters
    tg = sum(chr) # Total length of genome in bp
    tg < 2^32 || error("The total length of the genome is larger than 2^32 bp")
    pₘ = Poisson(tg / M * mr)  # Poisson distribution for mutations
    pᵣ = Poisson.(chr / M)     # Poisson distributions for crossovers
    cbp = UInt32.(cumsum(chr)) # Accumulated chromosome lengths in bp
    tbp = cbp[end] # Total length of genome in bp

    # Storage for mutations 
    nh = 2ne
    prt = [Vector{UInt32}() for _ = 1:nh]
    off = [Vector{UInt32}() for _ = 1:nh]

    # Random mating => a Fisher-Wright population
    @info "Fisher-Wright population:"
    for g = 1:nt
        if g % 100 == 0
            print('\r', ' '^8, "Generation: $g", " / ", nt)
            print(", mean mutations per haplotype: ", round(Int, mean(length.(prt))))
        end
        # create new mutations
        Threads.@threads for i = 1:nh
            nm = rand(pₘ)
            prt[i] = merge_sorted(prt[i], sort(UInt32.(rand(1:tbp, nm))))
        end

        pm, _ = random_mate(ne, ne) # pa and ma pairs

        # drop pa and ma's mutations to offspring
        Threads.@threads for i = 1:ne
            s, d = pm[i, :] # pa and ma for offspring i
            rss = cobp(cbp, pᵣ)
            recombine(prt[2s-1], prt[2s], off[2i-1], rss) # recombine pa haps
            rsd = cobp(cbp, pᵣ)
            recombine(prt[2d-1], prt[2d], off[2i], rsd) # recombine ma haps
        end
        prt, off = off, prt # swap
    end
    println()
    return prt, cbp
end

function fisher_wright(
    nid::Vector{T1}, # number of individuals in each generation
    chr::Vector{T2}, # chromosome lengths in bp,
    mr::Float64; # mutation rate per 1e8 base pair per meiosis
    M = 1e8,     # base pair per Morgan
) where {T1<:Integer,T2<:Integer}
    all([nid .> 1; chr .> 0; 0.01 < mr < 20.0]) || error("Invalid parameter(s)")
    cbp = UInt32.(cumsum(chr)) # Accumulated chromosome lengths in bp
    tcp = cbp[end] # Total length of genome in bp
    tcp < 2^32 || error("The total length of the genome can't be expressed in UInt32")
    pₘ = Poisson(tcp / M * mr)  # Poisson distribution for mutations
    pᵣ = Poisson.(chr / M)      # Poisson distributions for crossovers
    prt = [Vector{UInt32}() for _ = 1:2nid[1]] # generation 0 with 0 mutation

    @info "Fisher-Wright population:"
    for g = 1:length(nid)
        if g % 100 == 0
            print('\r', ' '^8, "Generation: $g", " / ", length(nid))
            print(", mean mutations per haplotype: ", round(Int, mean(length.(prt))))
        end
        # create new mutations
        nsd = length(prt) ÷ 2 # number of sires and dams
        Threads.@threads for i = 1:2nsd
            nm = rand(pₘ)
            prt[i] = merge_sorted(prt[i], sort(UInt32.(rand(1:tcp, nm))))
        end

        pm, _ = random_mate(nsd, nid[g]) # pa and ma pairs

        off = [Vector{UInt32}() for _ = 1:2nid[g]] # offspring haplotypes
        # drop pa and ma's mutations to offspring
        Threads.@threads for i = 1:nid[g]
            s, d = pm[i, :] # pa and ma for offspring i
            rss = cobp(cbp, pᵣ)
            recombine(prt[2s-1], prt[2s], off[2i-1], rss) # recombine pa haps
            rsd = cobp(cbp, pᵣ)
            recombine(prt[2d-1], prt[2d], off[2i], rsd) # recombine ma haps
        end
        prt = off # swap
    end
    println()
    return prt, cbp
end
