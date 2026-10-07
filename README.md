# FisherWright.jl

FisherWright.jl simulates neutral, diploid Wright-Fisher populations forward in
time, to generate whole-genome founder populations for breeding and
quantitative-genetics simulations. Starting from no variation, it applies random
mating between two sexes, Poisson crossovers and recurrent mutation for a chosen
number of generations. Run long enough (see [Recommended Use](#recommended-use)),
the population reaches mutation-drift equilibrium. The segregating variants are
exported as a `BitMatrix` with a linkage map, as a targeted SNP-chip panel, or as
[BnGStructs.jl](https://github.com/JuliaBnG/BnGStructs.jl) haplotypes.

[![Build Status](https://github.com/JuliaBnG/FisherWright.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/FisherWright.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaBnG/FisherWright.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaBnG/FisherWright.jl)

The complete manual, API reference and population genetics validation are
available at
[juliabng.github.io/FisherWright](https://juliabng.github.io/FisherWright/).

To build the static documentation for hosting:

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs --startup-file=no docs/make.jl
```

The generated site is in `docs/build/`.

## Features
- Forward simulation of a constant-size, randomly mating diploid population with
  two sexes, any number of autosomes and genomes up to 2³² − 1 bp
- Poisson crossovers (uniform rate, `M` bp per Morgan) with free assortment
  between chromosomes, and recurrent neutral mutation
- Optional extraction of fixed mutations during the run, which keeps memory
  bounded over long burn-ins
- Export to a 1-bit-per-allele `BitArray` with a linkage map (`DataFrame`), to a
  targeted SNP-chip panel, or to `BnGStructs.Haplotype`
- Validated against neutral theory and msprime (see the manual)

Not modelled: selection, changes in population size, population structure or
migration, and non-uniform recombination maps in `fisher_wright`.

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
| `verbose` | `false` | report progress (mean muts/hap, substitutions, elapsed, ETA) to `stderr` |
| `progress_interval` | `100` | generations between progress reports, when `verbose = true` |

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

## Recommended Use

The population starts with no variation, so it is only at mutation-drift
equilibrium after enough generations. Diversity approaches equilibrium at rate
about `1/(2ne)` per generation, so `nt = 10ne` gets within about 1% of
equilibrium heterozygosity. Shorter runs give a population that is still
gaining diversity: `ne = 2000, nt = 2000`, for example, is far from equilibrium.

```julia
using FisherWright

ne = 200
chr = fill(10_000_000, 10)  # 10 chromosomes of 10 Mbp

# Burn in for 10ne generations and track fixed mutations separately
res = fisher_wright(ne, 10ne, chr, 1.0; result = true)

# Expected neutral values at equilibrium, for a quick check:
# θ = 4ne·μ·L with μ = mr / mut_base per bp per generation
θ = 4ne * (1.0 / 1e8) * sum(chr)                  # Σ2pq ≈ θ
S = θ * sum(1 / i for i in 1:2ne-1)               # segregating sites ≈ θ·aₙ
```

- Use `result = true` for long runs. Fixed mutations are then moved out of the
  haplotypes into `substitutions`, which keeps memory bounded.
- `mr` is the mutation rate per `mut_base` bp per generation, and `M` sets
  recombination only, so the two can be changed independently.
- Runs with a fixed seed are reproducible only with the same thread count
  (`julia -t N`).
- For whole-genome populations where only chip markers are needed, use
  `extract_chip_bitarray` or `to_haplotype(result, chip_positions)` rather
  than building the full `BitMatrix`.

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

## Changes in v0.3.3

Added direct, allocation-conscious extraction of selected chip coordinates:
`extract_chip_bitarray` returns a `BitMatrix`, while
`to_haplotype(result, chip_positions)` returns a `BnGStructs.Haplotype`.
Coordinates must be sorted and unique. A `LocusSet` overload supports
subsetting a shared coordinate vector.

## Changes in v0.3.5

**Bug fix.** `fisher_wright` added each generation's mutations to the parents
before mating, so the population it returned had gone through one round of
reproduction since its last mutations and was missing its newest singletons.
At equilibrium this gave about 7-8% fewer segregating sites than Watterson's
θ·aₙ (a loss of about 0.46θ) and than msprime's discrete-time Wright-Fisher
model; heterozygosity (Σ2pq) was barely affected. Each generation now mates
and recombines first, then mutates the offspring. Segregating sites now match
θ·aₙ within replicate error.

**Simulation output changes.** Populations now contain more rare variants, and
results for a given seed differ from v0.3.4. On the 10 × 100 Mbp benchmark
(`ne = 2000`, `nt = 2000`), the number of segregating sites went from 528,369
to 570,905, against 569,235 from msprime.

**Tests.** Added a check against neutral theory: after `20ne` generations,
segregating sites must match θ·aₙ and Σ2pq must match θ. The old code fails it.

## Changes in v0.3.7

**Bug fix (Thread-safety).** Resolved a silent data race in multithreaded `BitMatrix` construction in `muts2bitarray` and `extract_chip_bitarray`. Because Julia's `BitArray` packs elements column-major in 64-bit (`UInt64`) words, any locus count with `nlc % 64 != 0` shares the boundary word between column $i$ and column $i+1$. Multiple threads setting bits in adjacent columns under `Threads.@threads` previously performed non-atomic read-modify-writes, causing genotype bits to be overwritten and lost. Replaced with atomic bitwise OR (`Core.Intrinsics.atomic_pointermodify`) directly on chunks, guaranteeing deterministic, lossless exports with zero heap allocations and full thread scaling.

**Tests.** Added a regression testset verifying bit-for-bit equivalence between serial and threaded export across non-64-aligned locus counts ($nlc = 65$) over repeated iterations.

## Changes in v0.3.8

**Bug fix (chromosome-boundary assortment).** The independent-assortment
breakpoint between chromosomes is now placed at the first coordinate of the
next chromosome (`cbp[i] + 1`), rather than at the preceding chromosome's
terminal coordinate. This prevents a boundary marker from being assigned to
the wrong segment during recombination. A regression test verifies that the
terminal locus of one chromosome and the first locus of the next segregate
into complementary offspring haplotypes.

**Validation and reproducibility.** CI now exercises the threaded code path
with four Julia threads. The repository includes reproducible benchmark drivers
and captured validation output for the 1 Gb transient comparison and founder
simulation example.

## Contributing and support

See [CONTRIBUTING.md](CONTRIBUTING.md). Bug reports and questions go to the
[issue tracker](https://github.com/JuliaBnG/FisherWright.jl/issues).

## License

MIT License. See LICENSE file for details.
