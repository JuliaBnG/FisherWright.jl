# Utilities

## Synthetic genotype and haplotype matrices

`quickGT` generates diploid genotype counts (0, 1, or 2) and `quickHap`
generates haploid allele calls (0 or 1). Allele frequencies are sampled from a
Beta distribution and restricted to the requested minor-allele-frequency range.

```@example utilities
using FisherWright
using Random

genotypes, frequencies = quickGT(100, 50; rng=MersenneTwister(1), return_p=true)
(size(genotypes), extrema(frequencies))
```

```@example utilities
haplotypes = quickHap(100, 50; rng=MersenneTwister(2))
size(haplotypes)
```

Both functions accept:

| Keyword | Default | Meaning |
|---|---|---|
| `maf` | `0.1` for `quickGT`, `0.2` for `quickHap` | Strict lower minor-allele-frequency bound. |
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
