_mean_length(v) = isempty(v) ? 0.0 : sum(length, v) / length(v)

"""
    fisher_wright(
        ne::Integer,
        nt::Integer,
        chr::Vector{<:Integer},
        mr::Float64;
        M = 1e8,
        mut_base = 1e8,
        result::Bool = false,
        fixation_interval::Int = 1,
        verbose::Bool = false,
    ) -> Tuple{Vector{Vector{UInt32}}, Vector{UInt32}} | FisherWrightResult

Simulate a diploid Fisher-Wright population of size `ne` over `nt` generations under
random mating, crossover recombination, and recurrent mutation.

Mutations are tracked as sorted, unique vectors of `UInt32` genomic coordinates, supporting
total genome lengths up to ``2^{32} \\approx 4.29 \\times 10^9`` base pairs (bp).

# Arguments
- `ne::Integer`: Effective diploid population size (number of individuals, `ne > 1`). Total haplotypes simulated is `2 * ne`.
- `nt::Integer`: Number of generations to simulate (`nt > 0`).
- `chr::Vector{<:Integer}`: Chromosome lengths in base pairs (all entries must be positive). Total genome length `sum(chr)` must be `< 2^32`.
- `mr::Float64`: Expected mutation rate per `mut_base` base pairs per meiosis (must satisfy `0.01 < mr < 20.0`).

# Keywords
- `M = 1e8`: Number of base pairs per Morgan (`M > 0`). A chromosome of length `L` experiences an expected `L / M` crossovers per meiosis.
- `mut_base = 1e8`: Base pairs per unit of mutation rate `mr` (`mut_base > 0`). Decoupled from `M` so adjusting recombination does not alter mutation intensity.
- `result::Bool = false`:
  - If `false` (default), returns `(prt, cbp)` where `prt` is the vector of `2 * ne` haplotypes and `cbp` is cumulative chromosome end coordinates. Fixed mutations remain in all haplotypes.
  - If `true`, returns a [`FisherWrightResult`](@ref), periodically extracting and tracking fixed substitutions to save memory and runtime.
- `fixation_interval::Int = 1`: Frequency (in generations) to scan for and extract fixed substitutions when `result = true`.
- `verbose::Bool = false`: If `true`, displays simulation progress every 100 generations. Default is `false`.

# Returns
- When `result = false`:
  - `prt::Vector{Vector{UInt32}}`: Active haplotypes (`2 * ne`), each containing sorted, unique 1-based genomic mutation positions.
  - `cbp::Vector{UInt32}`: Cumulative chromosome end positions in base pairs (`cumsum(chr)`).
- When `result = true`:
  - A [`FisherWrightResult`](@ref) containing `active_haplotypes`, `chromosome_ends`, and accumulated `substitutions`.

# Details
- **Recombination & Independent Assortment**: Recombination is modeled as a Poisson process along each chromosome with rate `chr[i] / M`. Independent assortment between chromosomes occurs with probability 0.5 at each chromosome boundary.
- **Multithreading & Reproducibility**: Meiosis and mutation generation are parallelized using `Threads.@threads`. Because task scheduling depends on `Threads.nthreads()`, runs with a fixed seed are strictly reproducible only when the thread count is held constant.
- **Allocation Efficiency**: Working buffers (mutations, crossovers, mating assignments) are allocated once at start and reused, achieving zero steady-state per-generation allocations.

# Examples
```julia
using FisherWright

# Simulate 50 individuals for 10 generations with two 100 kb chromosomes
haps, cbp = fisher_wright(50, 10, [100_000, 100_000], 1.0)
length(haps) # 100 haplotypes

# Return a structured result tracking fixed mutations
res = fisher_wright(50, 10, [100_000, 100_000], 1.0; result=true)
res.chromosome_ends
```
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
    ) -> Tuple{Vector{Vector{UInt32}}, Vector{UInt32}} | FisherWrightResult

Simulate a diploid Fisher-Wright population using species parameters from a
`BnGStructs.Species` instance (such as `Cattle`, `Pig`, `Chicken`, or `GenericSpecies`).

The population size `ne`, chromosome lengths `chr`, and default base pairs per Morgan `M`
are automatically extracted from `sp`:
- `ne = Int(sp.nid)`
- `chr = sp.chromosome`
- `M = sp.M`

# Arguments
- `sp::Species`: Species definition from `BnGStructs` specifying diploid population size (`nid`), chromosome lengths (`chromosome`), and recombination scale (`M`).
- `nt::Integer`: Number of generations to simulate (`nt > 0`).
- `mr::Float64 = 1.0`: Mutation rate per `mut_base` base pairs per meiosis (`0.01 < mr < 20.0`).

# Keywords
- `M = sp.M`: Number of base pairs per Morgan. Defaults to `sp.M`.
- `mut_base = 1e8`: Base pairs per unit of `mr`.
- `result::Bool = false`: If `true`, returns a [`FisherWrightResult`](@ref); otherwise returns `(prt, cbp)`.
- `fixation_interval::Int = 1`: Interval (in generations) to scan for and extract fixed substitutions when `result = true`.
- `verbose::Bool = false`: If `true`, prints progress every 100 generations.

# Examples
```julia
using BnGStructs, FisherWright

sp = GenericSpecies("Example", Int32(20), UInt32[100_000, 100_000], UInt32(50_000_000))
res = fisher_wright(sp, 10; result = true)
res.chromosome_ends
```
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
