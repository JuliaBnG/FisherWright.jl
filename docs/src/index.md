# FisherWright.jl

`FisherWright.jl` simulates diploid Fisher-Wright populations under mutation,
recombination, and random mating. It is designed for generating population
genetic data while keeping each haplotype as a sorted, unique vector of
`UInt32` mutation positions.

## Installation

```julia
using Pkg
Pkg.add("FisherWright")
```

Then load the package:

```julia
using FisherWright
```

## Quick start

Simulate 100 diploid individuals for 100 generations across two 1 Mbp
chromosomes, then convert the polymorphic mutations to a dense allele matrix
and linkage map:

```julia
mutations, chromosome_ends = fisher_wright(
    100,
    100,
    [1_000_000, 1_000_000],
    1.0,
)

alleles, loci = muts2bitarray(mutations, chromosome_ends)
```

Rows of `alleles` are loci and columns are haplotypes. `loci` is a
`DataFrame` containing chromosome, base-pair position, reference allele,
alternate allele, and allele frequency.

See [Simulation](@ref) for rate parameters and reproducibility, and [Results
and export](@ref) for fixed substitutions and `BnGStructs.Haplotype` export.

## Contents

```@contents
Pages = [
    "manual/simulation.md",
    "manual/results.md",
    "manual/utilities.md",
    "api.md",
]
Depth = 2
```