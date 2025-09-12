module FisherWright

using DataFrames
using Distributions
using Statistics

include("random-mate.jl")
include("merge-sorted.jl")
include("recombine.jl")
include("muts2bitarray.jl")
include("fwp.jl")
include("quick-genotypes.jl")

export fisher_wright, muts2bitarray, quickGT, quickHap

end # module FisherWright
