_mean_length(v) = isempty(v) ? 0.0 : sum(length, v) / length(v)

"""
    fisher_wright( # Fisher-Wright population
        ne::T1,
        nt::T2,
        chr::Vector{T3}, # chromosome lengths in bp,
        mr::Float64;     # mutation rate per `mut_base` base pairs per meiosis
        M = 1e8,         # base pair per Morgan
        mut_base = 1e8,  # base pairs per unit of `mr`
        result = false,
        fixation_interval = 1,
        verbose = false,
    ) where {T1<:Integer, T2<:Integer, T3<:Integer}
Simulate a Fisher-Wright population of size `ne` by `nt` generations of random
mating. The chromosome lengths are given in `chr` (in bp). The mutation rate is
`mr` per `mut_base` bp per meiosis, `mut_base` defaulting to 10⁸. The base pair
per Morgan is `M`, also with default value 10⁸. Mutations are stored as sorted
vectors of `UInt32` positions, so the maximum mutation place in bp is 2³², which
is 42.95 × 10⁸ bp. This is sufficient for most of the genomes. An error occurs
if the total length of the genome is larger than this value.

## Returns

1. The mutations
2. The accumulated base pairs on chrosomes

When `result = true` a [`FisherWrightResult`](@ref) is returned instead, which
additionally carries the positions that have fixed and been removed from the
active haplotypes.

## Note:
1. parameter recombination rate was removed as it is confusing. Crossover rate
   is a more appropriate term. Yet, Crossovers is also defined in M, or,
   bp/Morgan, i.e., the chromosome length in M is also the expected number of
   crossovers of this chromosome per meiosis. The less bp per Morgan, the more
   crossovers.
2. mutation rate is asked for per 1e8 base pair, the typical number of bp per M,
   per meiosis. This is to avoid tedious numbers like 1e-8. `M` and `mut_base`
   are separate parameters, so changing the recombination landscape through `M`
   leaves the mutation rate untouched.
3. `verbose = true` prints a progress line every 100 generations. The default is
   silent so that the function can be used inside libraries and benchmarks.
4. seeded runs reproduce only for a fixed `Threads.nthreads()`: the mutation and
   meiosis loops draw from task-local RNG streams whose seeds derive from how
   `Threads.@threads` splits the work.
5. with `result = false` fixed mutations are never extracted and stay in every
   haplotype for the whole run. Use `result = true` to reclaim them.
## The algorithm

1. Create two vectors of containers to store mutations in parents and offspring.
2. Generation of `M⋅m` new mutations ∈ ``[1, ~3×10⁹]`` for each haplotype
3. Insert these mutations into the parent generations
4. Randomly sample `Nₑ` pairs of ID in parent generation as sires and dams to
   offspring.
5. Splice sire and dam haplotypes into offspring's paternal haplotype
6. Swap parent and offspring storage
7. Repeat 1-6 for `nt` times.

All per-generation working storage (mutation buffers, merge targets, crossover
lists, mating tables) is allocated once and reused, so the steady-state
allocation rate is zero apart from haplotypes that grow.
"""
function fisher_wright(
    ne::T1,
    nt::T2,
    chr::Vector{T3},
    mr::Float64;
    M=1e8,
    mut_base=1e8,
    result::Bool=false,
    fixation_interval::Int=1,
    verbose::Bool=false,
) where {T1<:Integer,T2<:Integer,T3<:Integer}
    if !(ne > 1 && nt > 0 && all(chr .> 0) && 0.01 < mr < 20.0)
        throw(ArgumentError("Invalid parameter(s)"))
    end
    fixation_interval > 0 || throw(ArgumentError("fixation_interval must be positive"))
    mut_base > 0 || throw(ArgumentError("mut_base must be positive"))

    tg = sum(chr)
    tg < 2^32 || error("Total genome length must be < 2^32 bp for UInt32 storage")
    cbp = UInt32.(cumsum(chr))
    tbp = cbp[end]

    p_mut = Poisson(tg / mut_base * mr)
    recomb_map = uniform_recombination_map(chr; M=M)
    span = UInt32(1):UInt32(tbp)

    nh = 2 * ne
    prt = [Vector{UInt32}() for _ = 1:nh]
    off = [Vector{UInt32}() for _ = 1:nh]
    substitutions = UInt32[]

    # Reusable scratch: one merge target and one new-mutation buffer per
    # haplotype, one crossover buffer per mating, one mating table.
    mbuf = [Vector{UInt32}() for _ = 1:nh]
    nbuf = [Vector{UInt32}() for _ = 1:nh]
    cbuf = [Vector{UInt32}() for _ = 1:ne]
    pm = Matrix{Int}(undef, ne, 2)
    sex = Vector{Bool}(undef, ne)
    sires, dams = Int[], Int[]

    verbose &&
        @info "Fisher-Wright population simulation start" ne nt total_bp = Int(tg) threads =
            Threads.nthreads()

    for g = 1:nt
        if verbose && g % 100 == 0
            print(
                '\r',
                ' '^8,
                "Generation $g / $nt, mean muts/haps: ",
                round(_mean_length(prt); digits=2),
            )
        end
        # Mutations (threaded)
        Threads.@threads for i = 1:nh
            nm = rand(p_mut)
            if nm > 0
                newm = nbuf[i]
                resize!(newm, nm)
                @inbounds for k = 1:nm
                    newm[k] = rand(span)
                end
                nm <= 32 ? sort!(newm; alg=InsertionSort) : sort!(newm)
                # Merge into the scratch buffer, then swap it in: no allocation
                # once both buffers have reached their steady-state size.
                merge_sorted!(mbuf[i], prt[i], newm)
                prt[i], mbuf[i] = mbuf[i], prt[i]
            end
        end

        random_mate!(pm, sex, sires, dams)

        Threads.@threads for i = 1:ne
            s = pm[i, 1]
            d = pm[i, 2]
            # recombine empties its target, and cobp! its crossover buffer, so
            # both are safe to reuse across generations.
            co = cbuf[i]
            recombine(prt[2s-1], prt[2s], off[2i-1], cobp!(co, recomb_map))
            recombine(prt[2d-1], prt[2d], off[2i], cobp!(co, recomb_map))
        end
        prt, off = off, prt
        if result && (g % fixation_interval == 0 || g == nt)
            substitutions = _fixation_step!(prt, substitutions)
        end
    end
    verbose && println()
    return result ? FisherWrightResult(prt, cbp, substitutions) : (prt, cbp)
end

"""
    fisher_wright(
        sp::Species,
        nt::Integer,
        mr::Float64 = 1.0;
        M = sp.M,
        mut_base = 1e8,
        result::Bool = false,
        fixation_interval::Int = 1,
        verbose::Bool = false,
    )

Simulate a Fisher-Wright population using species parameters from a `BnGStructs.Species` object
(e.g., `Cattle(1000)`, `Pig(500)`, `Chicken(2000)`, or `GenericSpecies(...)`).

The population size `ne`, chromosome lengths `chr`, and base pairs per Morgan `M` are automatically
extracted from `sp`.
"""
function fisher_wright(
    sp::Species,
    nt::Integer,
    mr::Float64 = 1.0;
    M = sp.M,
    mut_base = 1e8,
    result::Bool = false,
    fixation_interval::Int = 1,
    verbose::Bool = false,
)
    return fisher_wright(
        Int(sp.nid),
        nt,
        sp.chromosome,
        mr;
        M = M,
        mut_base = mut_base,
        result = result,
        fixation_interval = fixation_interval,
        verbose = verbose,
    )
end
