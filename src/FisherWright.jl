"""
    FisherWright

A high-performance forward-time Fisher-Wright population simulator for diploid individuals
under mutation, crossover recombination, and random mating.

Haplotypes are represented as sorted, unique vectors of `UInt32` mutation coordinates
(supporting genomes up to ``2^{32} \\approx 4.29 \\times 10^9`` base pairs), eliminating
dense matrix allocations during generation stepping.

# Main Features
- **Simulation**: [`fisher_wright`](@ref) supports both explicit parameter configurations and
  [`BnGStructs.Species`](https://juliabng.github.io/BnGStructs.jl/stable/) presets.
- **Fixed Mutation Tracking**: Automatically extracts fixed substitutions during simulation via
  [`FisherWrightResult`](@ref).
- **Dense Export**: Convert sparse mutation representations into dense bit arrays and linkage maps
  via [`muts2bitarray`](@ref) and [`to_haplotype`](@ref).
- **Targeted Chip Extraction**: Directly extract genotype arrays at predefined marker coordinates
  without full genome materialization via [`extract_chip_bitarray`](@ref) or [`to_haplotype`](@ref).
- **Recombination Landscapes**: Custom or uniform crossover recombination maps via
  [`RecombinationMap`](@ref) and [`uniform_recombination_map`](@ref).
- **Synthetic Data**: Generate fast synthetic diploid genotypes and haplotypes with controlled
  allele frequencies via [`quickGT`](@ref) and [`quickHap`](@ref).
"""
module FisherWright


using BnGStructs
using DataFrames
using Distributions
using Random
using Statistics

include("random-mate.jl")
include("merge-sorted.jl")
include("recombine.jl")
include("muts2bitarray.jl")
include("result.jl")
include("fwp.jl")
include("quick-genotypes.jl")

export fisher_wright, muts2bitarray, extract_chip_bitarray, quickGT, quickHap, FisherWrightResult, to_haplotype, RecombinationMap, uniform_recombination_map

end # module FisherWright
