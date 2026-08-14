"""
    RecombinationMap

Immutable recombination map with chromosome cumulative end positions and a
piecewise-constant crossover intensity per interval.

Fields `cbp`, `interval_ends` and `interval_rates` describe the map; `samplers`
holds the `Poisson` distribution for each interval, built once at construction
so that [`cobp!`](@ref) does not rebuild them on every meiosis.
"""
struct RecombinationMap
    cbp::Vector{UInt32}
    interval_ends::Vector{Vector{UInt32}}
    interval_rates::Vector{Vector{Float64}}
    samplers::Vector{Vector{Poisson{Float64}}}

    function RecombinationMap(
        cbp::Vector{UInt32},
        interval_ends::Vector{Vector{UInt32}},
        interval_rates::Vector{Vector{Float64}},
    )
        length(cbp) == length(interval_ends) == length(interval_rates) ||
            throw(ArgumentError("map dimensions do not match chromosome count"))
        prev_chr_end = UInt32(0)
        for i in eachindex(cbp)
            ends = interval_ends[i]
            rates = interval_rates[i]
            length(ends) == length(rates) ||
                throw(ArgumentError("interval ends and rates mismatch on chromosome $i"))
            isempty(ends) && throw(ArgumentError("chromosome $i has no intervals"))
            last(ends) == cbp[i] ||
                throw(ArgumentError("interval ends must terminate at chromosome end"))
            first_start = prev_chr_end + UInt32(1)
            ends[1] >= first_start ||
                throw(ArgumentError("invalid first interval end on chromosome $i"))
            for j in eachindex(rates)
                rates[j] >= 0 || throw(ArgumentError("interval rates must be nonnegative"))
                if j > 1
                    ends[j] > ends[j-1] ||
                        throw(ArgumentError("interval ends must be strictly increasing"))
                end
            end
            prev_chr_end = cbp[i]
        end
        samplers = [[Poisson(r) for r in rates] for rates in interval_rates]
        new(cbp, interval_ends, interval_rates, samplers)
    end
end

"""
    uniform_recombination_map(chr; M=1e8)

Create the default one-interval-per-chromosome recombination map used by the
existing Fisher-Wright model. `M` is the number of base pairs per Morgan, so
chromosome `i` gets an expected `chr[i] / M` crossovers per meiosis.
"""
function uniform_recombination_map(chr::Vector{T}; M=1e8) where {T<:Integer}
    all(chr .> 0) || throw(ArgumentError("chromosome lengths must be positive"))
    M > 0 || throw(ArgumentError("M must be positive"))
    cbp = UInt32.(cumsum(chr))
    ends = [UInt32[cbp[i]] for i in eachindex(cbp)]
    rates = [Float64[chr[i]/M] for i in eachindex(chr)]
    return RecombinationMap(cbp, ends, rates)
end

"""
    uniform_recombination_map(sp::Species; M=sp.M)

Create the default one-interval-per-chromosome recombination map for a `BnGStructs.Species` object.
"""
uniform_recombination_map(sp::Species; M=sp.M) = uniform_recombination_map(sp.chromosome; M=M)

"""
    recombine(h₁, h₂, hₒ, cross_overs; rng = Random.default_rng())

Recombine two parental haplotypes into an offspring haplotype.

`h₁` and `h₂` are the parent's two haplotypes, each a sorted vector of unique
mutation positions. `cross_overs` is a sorted vector of crossover positions.
Starting from a randomly chosen parental haplotype, the positions falling in
each successive segment are copied from alternating parents into `hₒ`, which is
emptied first and returned.

A position exactly equal to a crossover point belongs to the segment that
starts at that crossover.

`hₒ` is reused in place, so passing the same buffer across generations avoids
reallocation. The result is sorted and unique whenever the inputs are.

Preconditions: `h₁`, `h₂` and `cross_overs` are sorted in nondecreasing order.
"""
function recombine(
    h₁::Vector{UInt32},
    h₂::Vector{UInt32},
    hₒ::Vector{UInt32},
    cross_overs::Vector{UInt32};
    rng::AbstractRNG = Random.default_rng(),
)
    empty!(hₒ)

    m, n = length(h₁), length(h₂)
    i = j = 1
    # Randomly select which haplotype to start with
    o = rand(rng, Bool)

    @inbounds for co in cross_overs
        # First index at or beyond the crossover; the segment is [i, i₂-1].
        i₂ = searchsortedfirst(h₁, co, i, m, Base.Order.Forward)
        j₂ = searchsortedfirst(h₂, co, j, n, Base.Order.Forward)
        if o
            i₂ > i && append!(hₒ, view(h₁, i:(i₂-1)))
        else
            j₂ > j && append!(hₒ, view(h₂, j:(j₂-1)))
        end
        i, j = i₂, j₂
        o = !o # flip haplotype
    end

    # Trailing segment, from the last crossover to the end of the genome.
    if o
        i <= m && append!(hₒ, view(h₁, i:m))
    else
        j <= n && append!(hₒ, view(h₂, j:n))
    end
    return hₒ
end

"""
    cobp!(dest, map::RecombinationMap; rng = Random.default_rng())

Sample crossover positions from a piecewise-constant recombination map into
`dest`, which is emptied first and returned.

Within each interval the number of crossovers is Poisson with the interval's
rate and the positions are uniform. Each chromosome end is additionally emitted
as a crossover with probability 0.5, which makes chromosomes assort
independently.

This is the allocation-free form of [`cobp`](@ref): reusing `dest` across
meioses removes the per-meiosis allocations entirely.
"""
function cobp!(
    dest::Vector{UInt32},
    map::RecombinationMap;
    rng::AbstractRNG = Random.default_rng(),
)
    empty!(dest)
    prev_chr_end = UInt32(0)
    @inbounds for i in eachindex(map.cbp)
        chrom_start = prev_chr_end + UInt32(1)
        ends = map.interval_ends[i]
        samplers = map.samplers[i]
        for j in eachindex(ends)
            interval_start = j == 1 ? chrom_start : ends[j-1] + UInt32(1)
            interval_end = ends[j]
            nr = rand(rng, samplers[j])
            if nr > 0 && interval_start < interval_end
                base = length(dest)
                span = interval_start:(interval_end-UInt32(1))
                for _ = 1:nr
                    push!(dest, rand(rng, span))
                end
                if nr > 1
                    segment = view(dest, (base+1):length(dest))
                    nr <= 32 ? sort!(segment; alg=InsertionSort) : sort!(segment)
                end
            end
        end
        rand(rng) < 0.5 && push!(dest, map.cbp[i])
        prev_chr_end = map.cbp[i]
    end
    return dest
end

"""
    cobp(map::RecombinationMap; rng = Random.default_rng())
- Generate crossover points from a piecewise-constant recombination map.

Allocates a fresh vector; use [`cobp!`](@ref) with a reused buffer in hot loops.
"""
cobp(map::RecombinationMap; rng::AbstractRNG = Random.default_rng()) =
    cobp!(UInt32[], map; rng = rng)

"""
    function cobp(cbp::Vector{UInt32}, pᵣ::Vector{Poisson{Float64}})
- Generate crossover points for recombination.

The crossover points are generated based on the chromosome breakpoints and the
Poisson distribution of recombination events. The crossover points are generated
in the range of [1, bp-1] for each chromosome plus the accumulated breakpoints.
The crossover points are sorted and returned as a vector of UInt32.
"""
function cobp(
    cbp::Vector{UInt32},
    pᵣ::Vector{Poisson{Float64}};
    rng::AbstractRNG = Random.default_rng(),
)
    length(cbp) == length(pᵣ) || throw(ArgumentError("cbp and pᵣ must have the same length"))
    ends = [UInt32[cbp[i]] for i in eachindex(cbp)]
    rates = [Float64[mean(pᵣ[i])] for i in eachindex(pᵣ)]
    return cobp(RecombinationMap(cbp, ends, rates); rng = rng)
end
