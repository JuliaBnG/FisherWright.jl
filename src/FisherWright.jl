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

export fisher_wright, muts2bitarray, quickGT, quickHap, FisherWrightResult, to_haplotype, RecombinationMap, uniform_recombination_map

end # module FisherWright
