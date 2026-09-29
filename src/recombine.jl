"""
    RecombinationMap(cbp, interval_ends, interval_rates)

Piecewise-constant crossover recombination map defined across chromosomes.

Stores cumulative chromosome ends, interval boundaries, and expected crossover rates
(Poisson intensities) per interval. Pre-computes Poisson samplers upon construction
to maximize performance in meiosis loops ([`cobp!`](@ref)).

# Fields
- `cbp::Vector{UInt32}`: Cumulative chromosome end positions in base pairs.
- `interval_ends::Vector{Vector{UInt32}}`: End coordinates for each recombination interval within each chromosome.
- `interval_rates::Vector{Vector{Float64}}`: Expected crossover counts (Poisson mean) per interval.
- `samplers::Vector{Vector{Poisson{Float64}}}`: Pre-constructed Poisson distributions for efficient sampling.

# Invariants
- The outer lengths of `cbp`, `interval_ends`, and `interval_rates` must match the number of chromosomes.
- Within each chromosome, `interval_ends[i]` and `interval_rates[i]` must have the same non-zero length.
- `interval_ends` must be strictly increasing and terminate at `cbp[i]`.
- `interval_rates` must be non-negative.

# See also
[`uniform_recombination_map`](@ref), [`cobp!`](@ref), [`cobp`](@ref)
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
    uniform_recombination_map(chr::Vector{<:Integer}; M = 1e8) -> RecombinationMap
    uniform_recombination_map(sp::Species; M = sp.M) -> RecombinationMap

Construct a uniform, single-interval-per-chromosome [`RecombinationMap`](@ref).

# Arguments
- `chr::Vector{<:Integer}`: Chromosome lengths in base pairs (all positive).
- `sp::Species`: A `BnGStructs.Species` object supplying chromosome lengths (`sp.chromosome`) and default `M` (`sp.M`).

# Keywords
- `M`: Base pairs per Morgan (default `1e8`, or `sp.M` for species). Chromosome `i` of length `L` receives an expected `L / M` crossovers per meiosis.

# Examples
```julia
using FisherWright

rmap = uniform_recombination_map([100_000, 200_000]; M = 1e8)
rmap.cbp == UInt32[100_000, 300_000]
```
"""
function uniform_recombination_map(chr::Vector{T}; M=1e8) where {T<:Integer}
    all(chr .> 0) || throw(ArgumentError("chromosome lengths must be positive"))
    M > 0 || throw(ArgumentError("M must be positive"))
    cbp = UInt32.(cumsum(chr))
    ends = [UInt32[cbp[i]] for i in eachindex(cbp)]
    rates = [Float64[chr[i]/M] for i in eachindex(chr)]
    return RecombinationMap(cbp, ends, rates)
end

uniform_recombination_map(sp::Species; M=sp.M) = uniform_recombination_map(sp.chromosome; M=M)

"""
    recombine(
        h₁::Vector{UInt32},
        h₂::Vector{UInt32},
        hₒ::Vector{UInt32},
        cross_overs::Vector{UInt32};
        rng::AbstractRNG = Random.default_rng(),
    ) -> Vector{UInt32}

Recombine two parental haplotypes `h₁` and `h₂` into an offspring haplotype `hₒ`
given a sorted list of `cross_overs`.

# Arguments
- `h₁::Vector{UInt32}`: First parental haplotype (sorted, unique mutation positions).
- `h₂::Vector{UInt32}`: Second parental haplotype (sorted, unique mutation positions).
- `hₒ::Vector{UInt32}`: Destination vector for the offspring haplotype. Emptied first and returned.
- `cross_overs::Vector{UInt32}`: Sorted crossover coordinates defining segments.

# Keywords
- `rng::AbstractRNG = Random.default_rng()`: Random number generator used to select which parental haplotype begins the sequence.

# Details
- Randomly chooses the starting parent (`h₁` or `h₂`) with equal probability (0.5).
- Swaps active parental templates at each coordinate in `cross_overs`.
- A mutation coordinate exactly matching a crossover point belongs to the segment starting at that crossover.
- Reuses `hₒ` in place without reallocation when previously sized, guaranteeing sorted and unique output.

# Preconditions
- `h₁`, `h₂`, and `cross_overs` must be sorted in non-decreasing order.
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
    cobp!(dest::Vector{UInt32}, map::RecombinationMap; rng::AbstractRNG = Random.default_rng()) -> Vector{UInt32}

Sample crossover positions from a [`RecombinationMap`](@ref) into `dest` in place,
emptying `dest` first and returning it.

# Details
- Within each interval of `map`, the number of crossover events is drawn from a Poisson distribution
  with the interval's rate, and coordinates are distributed uniformly at random.
- Each chromosome boundary (`map.cbp[i]`) is added as a crossover breakpoint with probability 0.5,
  simulating independent chromosome assortment during meiosis.
- Reusing `dest` across successive meioses eliminates allocations.
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
    cobp(map::RecombinationMap; rng::AbstractRNG = Random.default_rng()) -> Vector{UInt32}
    cobp(cbp::Vector{UInt32}, pᵣ::Vector{Poisson{Float64}}; rng::AbstractRNG = Random.default_rng()) -> Vector{UInt32}

Generate crossover breakpoints for recombination.

# Methods
- `cobp(map::RecombinationMap; rng)`: Allocates a new vector and samples crossover positions using [`cobp!`](@ref).
- `cobp(cbp, pᵣ; rng)`: Constructs a single-interval-per-chromosome [`RecombinationMap`](@ref) from
  cumulative chromosome boundaries `cbp` and per-chromosome Poisson distributions `pᵣ`, then samples crossover points.

# See also
[`cobp!`](@ref), [`RecombinationMap`](@ref)
"""
cobp(map::RecombinationMap; rng::AbstractRNG = Random.default_rng()) =
    cobp!(UInt32[], map; rng = rng)

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
