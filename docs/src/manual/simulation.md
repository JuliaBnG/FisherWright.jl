# Simulation

## `fisher_wright`

```julia
fisher_wright(ne, nt, chr, mr; M=1e8, mut_base=1e8,
              result=false, fixation_interval=1, mutate=true,
              init=nothing, popsize=nothing, mating=RandomMating(),
              verbose=false, progress_interval=100)
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
| `mutate` | `true` | Generate new mutations; `false` leaves drift and recombination only. |
| `init` | `nothing` | Starting population: `2ne` haplotypes or a `FisherWrightResult`; `nothing` starts monomorphic. |
| `popsize` | `nothing` | Census size per generation: a vector of `nt` sizes or a function of the generation; `nothing` keeps `ne`. |
| `mating` | `RandomMating()` | Mating scheme; `RandomMating(nmale=m)` or `RandomMating(sexratio=p)` fixes the number or the proportion of males. |
| `verbose` | `false` | Report progress (mean mutations per haplotype, substitutions, elapsed time, ETA) to `stderr`. |
| `progress_interval` | `100` | Generations between progress reports when `verbose=true`. |

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

## Continuing a run

By default a run starts monomorphic. `init` supplies the starting
population instead, either as a vector of `2ne` haplotypes (each a
sorted, unique `Vector{UInt32}` of positions in `1:sum(chr)`) or as a
`FisherWrightResult` with the same chromosomes. Haplotypes `2i - 1` and
`2i` form individual `i`, as in the output, and these individuals are
the parents of generation 1. `init` is copied and left unmodified.

```@example simulation
burnin = fisher_wright(20, 200, [100_000, 100_000], 1.0; result=true)
replicate = fisher_wright(20, 50, [100_000, 100_000], 1.0;
                          init=burnin, result=true)
length(replicate.substitutions) >= length(burnin.substitutions)
```

The substitutions of a `FisherWrightResult` are carried over when
`result=true` and merged back into every haplotype when `result=false`,
so both return modes still describe the same population.

With `mutate=false` no new mutations arise, and the starting variants
change only by drift and recombination. This is the setting for checks
such as the decay of heterozygosity, ``H_t/H_0 = (1 - 1/(2N_e))^t``:

```@example simulation
drift = fisher_wright(20, 50, [100_000, 100_000], 1.0;
                      init=burnin, mutate=false, result=true)
issubset(reduce(union, drift.active_haplotypes),
         reduce(union, burnin.active_haplotypes))
```

The per-slot random streams are seeded anew from `rng` on every call.
Two calls of `nt₁` and `nt₂` generations are therefore a valid
realization of `nt₁ + nt₂` generations, but not the one a single call
gives for the same seed.

## Changing population size

By default the census size stays at `ne`. `popsize` gives the size of
every generation instead, either as a vector of `nt` integers or as a
function `g -> N_g`, which is called once for each `g = 1, …, nt` before
the run. Every size must exceed one. `ne` is then the starting size
``N_0``, and the run returns ``2N_{nt}`` haplotypes.

```@example simulation
sizes = [fill(10, 5); fill(40, 15)]   # bottleneck, then recovery
haplotypes, _ = fisher_wright(40, 20, [100_000, 100_000], 1.0;
                              popsize=sizes)
length(haplotypes)
```

```@example simulation
growth = g -> round(Int, 20 * 1.05^g)   # 5% growth per generation
haplotypes, _ = fisher_wright(20, 20, [100_000, 100_000], 1.0;
                              popsize=growth)
length(haplotypes) == 2growth(20)
```

Generation ``g`` draws the sire and dam of each of its ``N_g`` offspring
from the ``N_{g-1}`` individuals of the previous generation, whose sexes
are drawn as in the constant case. The drift of generation ``g``
therefore follows the number of parents, ``N_{g-1}``: its effective size
is that of the two-sex model for ``N_{g-1}`` individuals (see
[Effective size](@ref)), whatever the number of offspring. Mutation is
per haplotype, so the number of new mutations follows ``N_g``.

The population is not at mutation–drift equilibrium after a change in
size; combine `popsize` with `init` to burn in at constant size first.

Buffers and random streams are allocated once for the largest size of
the run, so the realization for a given seed depends on that maximum. A
trajectory that stays at `ne` reproduces the run without `popsize`.

## Sex ratio

Sex is not inherited. It is assigned anew to each cohort when it
reproduces, and every offspring then draws its sire uniformly from the
males and its dam uniformly from the females, with replacement. The
`mating` keyword takes a `RandomMating` scheme that sets how sex is
assigned to the ``N`` parents of a generation:

| Scheme | Males among ``N`` parents |
|---|---|
| `RandomMating()` | Binomial(``N``, ½), conditioned on both sexes being present (the default) |
| `RandomMating(nmale=m)` | exactly `m`, a uniformly random subset |
| `RandomMating(sexratio=p)` | exactly `round(Int, p * N)`, at least one of each sex |

```@example simulation
# 5 sires and 35 dams in every generation
haplotypes, _ = fisher_wright(40, 20, [100_000, 100_000], 1.0;
                              mating=RandomMating(nmale=5))
length(haplotypes)
```

With ``N_m`` males and ``N_f`` females the effective size of a
generation is ``N_e \approx 4N_mN_f/(N_m + N_f) + 1/2`` (the exact
recursion is under [Effective size](@ref), with
``a = 1/(4N_m) + 1/(4N_f)``): 18 in the example above, against a census
size of 40. `RandomMating(sexratio=0.5)` removes the variation in the
sex ratio of the default scheme, so that ``N_e \approx N + 1/2``
instead of ``N - 1/2`` for even ``N``.

`nmale` is a count and must be smaller than every parental cohort, which
matters when it is combined with `popsize`; `sexratio` is a proportion
of males and follows the census size. The fixed schemes draw sex
differently from `rng`, so for the same seed they share no realization
with the default scheme.

## Reproducibility

Pass an explicit generator, for example `rng = Xoshiro(2026)`.
`fisher_wright` seeds one stream per meiosis and haplotype slot from
`rng` and samples mating serially from it, so a run depends only on the
state of `rng`, not on the number of threads, the fixation interval, or
the return mode.

This holds for a fixed Julia version and package environment. The
realization also depends on the variate algorithms of Julia's `Random`
and of Distributions.jl, which can change between releases. For example,
`randn` differs between Julia 1.10 and 1.13, so Poisson draws with mean
≥ 6 can differ; in most genome-scale runs that includes the number of
new mutations per haplotype, and the same seed then gives a different
population. Keep the Julia version and `Manifest.toml` with simulated
data.

## Parameter limits

`ne` and every entry of `popsize` must exceed one, `nt` must be
positive, and chromosome lengths must be positive. `mr` must be strictly
between `0.01` and `20.0`. The total genome length must be smaller than
`2^32`, because positions use `UInt32`.
