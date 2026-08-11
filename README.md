# FisherWright.jl

FisherWright.jl is a Julia package for simulating mutation-drift equilibrium
Fisher-Wright populations.  It provides efficient tools for modeling population
genetics, generating haplotypes, and converting mutation data into bit arrays
and linkage maps. This package is suitable for researchers and students in
population genetics, evolutionary biology, and related fields.

[![Build Status](https://github.com/JuliaBnG/FisherWright.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/FisherWright.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaBnG/FisherWright.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaBnG/FisherWright.jl)

The complete manual and API reference are available at
[xijiang.org/JuliaBnG/FisherWright](https://xijiang.org/JuliaBnG/FisherWright).

To build the static documentation for hosting:

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs --startup-file=no docs/make.jl
```

The generated site is in `docs/build/`.

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

Keyword arguments:

| Keyword | Default | Meaning |
|---|---|---|
| `M` | `1e8` | base pairs per Morgan (recombination only) |
| `mut_base` | `1e8` | base pairs per unit of `mr` (mutation only) |
| `result` | `false` | return a `FisherWrightResult` and extract fixed positions |
| `fixation_interval` | `1` | generations between fixation scans, when `result = true` |
| `verbose` | `false` | print a progress line every 100 generations |

Every haplotype is a sorted vector of unique `UInt32` positions, and that
invariant holds for the returned population.

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
Pkg.add("FisherWright")
```

## Usage Example
```julia
using FisherWright

# Simulate a population
muts, cbp = fisher_wright(100, 1000, [1_000_000, 1_000_000], 1.0)

# Convert mutations to bit array and linkage map
xy, lmp = muts2bitarray(muts, cbp; flip = true)
```

## Benchmark

Run the reproducible benchmark script from the package directory:

```bash
julia --project=. bench/benchmark_phase3.jl
```

Optional parameters:

```bash
julia --project=. bench/benchmark_phase3.jl 120 200 100000000,100000000 0.5 5
```

Arguments are interpreted as:

1. `ne`
2. `nt`
3. comma-separated chromosome lengths
4. `mr`
5. repetitions
6. seed

The script reports wall time, cumulative allocation, GC fraction and peak RSS,
and prints the thread count it ran with. Pin threads with `julia -t N` when
comparing against another simulator: cumulative allocation is not a memory
footprint, and wall time depends directly on the thread count.

## Changes in v0.3.0

**Bug fix.** `recombine` skipped the last position of each parental haplotype
inside the crossover loop, emitting it only through the trailing append. At
29 × 100 Mbp with 30 000 mutations per haplotype this made 46% of meioses
inherit the wrong positions and left 11% of offspring haplotypes unsorted,
breaking the sorted-and-unique invariant the rest of the package relies on.
Simulation output changes accordingly.

**Behaviour change.** `fisher_wright` is now silent by default; pass
`verbose = true` for the old progress output.

**API.** `M` no longer scales the mutation rate — use the new `mut_base`
keyword for that. Added allocation-free forms `merge_sorted!`, `cobp!` and
`random_mate!` for use in hot loops.

**Performance.** Per-generation working storage is allocated once and reused,
and the fixation scan uses sorted merges instead of hash sets. Measured on
29 × 100 Mbp, `mr = 1.0`, `result = true`:

| | before | after |
|---|---|---|
| ne=250, nt=200, 12 threads | 1.44 s, 3.12 GiB, 55% GC | 0.27 s, 0.07 GiB, 0.4% GC |
| ne=500, nt=200, 12 threads | 2.37 s, 6.04 GiB, 65% GC | 0.59 s, 0.15 GiB, 7.8% GC |
| ne=250, nt=200, 1 thread | 5.83 s, 3.08 GiB, 70% GC | 0.95 s, 0.07 GiB, 0.5% GC |

The meiosis inner loop (`cobp!` plus `recombine`) no longer allocates at all in
steady state.

## License

MIT License. See LICENSE file for details.

## ToDo

1. [x] Add dependency `BnGStructs`, return mutations as a `Genotype`, or
   `Haplotype`.
2. [ ] Add split and merge sub-populations
