# Simulation

## `fisher_wright`

```julia
fisher_wright(ne, nt, chr, mr; M=1e8, mut_base=1e8,
              result=false, fixation_interval=1, verbose=false)
```

Simulate `ne` diploid individuals for `nt` generations. `chr` is a vector of
positive chromosome lengths in base pairs, and `mr` is the expected number of
new mutations per `mut_base` base pairs per haplotype per meiosis.

| Keyword | Default | Meaning |
|---|---:|---|
| `M` | `1e8` | Base pairs per Morgan; controls recombination only. |
| `mut_base` | `1e8` | Base pairs per unit of `mr`; controls mutation only. |
| `result` | `false` | Return a `FisherWrightResult` and extract fixed mutations. |
| `fixation_interval` | `1` | Generations between fixation scans when `result=true`. |
| `verbose` | `false` | Print progress every 100 generations. |

With the default `result=false`, the function returns
`(haplotypes, chromosome_ends)`. There are `2ne` haplotypes, and
`chromosome_ends` contains cumulative chromosome endpoints. Mutation positions
are genome-wide coordinates, not per-chromosome coordinates.

```@example simulation
using FisherWright

haplotypes, chromosome_ends = fisher_wright(20, 10, [100_000, 100_000], 1.0)
(length(haplotypes), chromosome_ends)
```

## Species parameters

When using a `BnGStructs.Species`, pass the species object directly instead of
separately supplying its population size and chromosome lengths. The default
`M` is taken from the species, while all simulation keywords remain available:

```@example simulation
using BnGStructs

species = GenericSpecies(
    "Example",
    Int32(20),
    UInt32[100_000, 100_000],
    UInt32(50_000_000),
)
result = fisher_wright(species, 10; result=true)
result.chromosome_ends
```

## Recombination

`M` gives the number of base pairs per Morgan. A chromosome of length `L`
therefore has an expected `L / M` crossovers per meiosis. The default
recombination map has one uniform interval per chromosome:

```@example simulation
map = uniform_recombination_map([100_000, 250_000]; M=1e8)
map.cbp
```

The species overload uses the chromosome lengths and default `M` from the
species:

```@example simulation
uniform_recombination_map(species).cbp
```

Each chromosome end is treated as an independent assortment boundary. The
package does not expose a custom map through `fisher_wright`; `RecombinationMap`
and `cobp!` are available for advanced internal workflows.

## Reproducibility

Pass a seeded task-local random-number generator only to the utility functions
that accept `rng`. `fisher_wright` uses threaded loops, so seeded simulation
runs are reproducible only when the number of Julia threads is unchanged. Set
the thread count explicitly when comparing simulations:

```bash
julia -t 1 --project=. script.jl
```

## Parameter limits

`ne` must exceed one, `nt` must be positive, and chromosome lengths must be
positive. `mr` must be strictly between `0.01` and `20.0`. The total genome
length must be smaller than `2^32`, because positions use `UInt32`.
