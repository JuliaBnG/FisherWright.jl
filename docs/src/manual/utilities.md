# Utilities

## Synthetic haplotypes

`quickHap` generates phased haplotypes at independent loci as a
`BnGStructs.Haplotype` (`true` for the alternate allele). Allele frequencies
are sampled from a Beta distribution and restricted to the requested
minor-allele-frequency range. For genotypes with linkage and the allele
frequency spectrum of a population, simulate a small historical population
with `fisher_wright` instead; it is cheap.

```@example utilities
using FisherWright
using Random

haplotypes, frequencies = quickHap(100, 50; rng=MersenneTwister(2), return_p=true)
(typeof(haplotypes), size(haplotypes), extrema(frequencies))
```

It accepts:

| Keyword | Default | Meaning |
|---|---|---|
| `maf` | `0.2` | Strict lower minor-allele-frequency bound. |
| `qd` | `Beta(0.75, 0.75)` | Distribution used to sample allele frequencies. |
| `rng` | `Random.default_rng()` | Random-number generator. |
| `return_p` | `false` | Return sampled frequencies with the matrix. |

`maf` must be strictly between zero and `0.5`.

## Efficient internal helpers

For repeated low-level simulations, FisherWright provides in-place helpers:

- `merge_sorted!` writes a sorted set union into a reusable destination vector.
- `cobp!` samples crossover positions into a reusable destination vector.
- `random_mate!` writes mating pairs and a sex assignment into supplied buffers.

These functions are intentionally not part of the public export list. They are
documented in the [API reference](@ref) for advanced users, but their
interfaces may change more freely than the exported simulation and export API.
