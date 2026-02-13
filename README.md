# FisherWright.jl

FisherWright.jl is a Julia package for simulating mutation-drift equilibrium
Fisher-Wright populations.  It provides efficient tools for modeling population
genetics, generating haplotypes, and converting mutation data into bit arrays
and linkage maps. This package is suitable for researchers and students in
population genetics, evolutionary biology, and related fields.

[![Build Status](https://github.com/JuliaBnG/FisherWright.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/FisherWright.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaBnG/FisherWright.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaBnG/FisherWright.jl)

## Features
- Simulate Fisher-Wright populations of arbitrary population size, multiple auto
  chromosomes and over multiple generations
- Population merge and splits.
- Model recombination and mutation processes on chromosomes
- Convert mutation data to `BitArray` and linkage map (`DataFrame`)
- Efficient handling of large-scale genomic data

## Main Functions

### `fisher_wright`

Simulates a Fisher-Wright population of a given size and number of generations,
with specified chromosome lengths and mutation rates. Returns simulated
population data with mutations and recombination events.

### `muts2bitarray`

Converts a vector of mutation sets (haplotypes) and chromosome breakpoints into
a `BitArray` of haplotypes and a linkage map (`DataFrame`). Supports random
flipping of alleles and removal of fixed loci.

### `quickGT`, `quickHap`

Utility functions for fast genotype and haplotype generation from simulated
data.

## Installation

Add the package to your Julia environment:

```julia
using Pkg
Pkg.add(FisherWright)
```

## Usage Example
```julia
using FisherWright

# Simulate a population
muts, cbp = fisher_wright(100, 1000, [1_000_000, 1_000_000], 1.0)

# Convert mutations to bit array and linkage map
xy, lmp = muts2bitarray(muts, cbp; flip = true)
```

## License

MIT License. See LICENSE file for details.

## ToDo

1. [x] Add dependency `BnGStructs`, return mutations as a `Genotype`, or
   `Haplotype`.
2. [ ] Add split and merge sub-populations
