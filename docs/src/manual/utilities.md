# Utilities

Test data with linkage and the allele frequency spectrum of a population
come from a small historical population: `fisher_wright` with a few hundred
individuals and generations takes a fraction of a second, and
`to_haplotype` turns it into a `BnGStructs.Haplotype` and its
`VariantMap`. (`quickGT` and `quickHap`, which drew genotypes at
independent loci, were removed in 0.10.)

## Efficient internal helpers

For repeated low-level simulations, FisherWright provides in-place helpers:

- `merge_sorted!` writes a sorted set union into a reusable destination vector.
- `cobp!` samples crossover positions into a reusable destination vector.
- `random_mate!` writes mating pairs and a sex assignment into supplied buffers.

These functions are intentionally not part of the public export list. They are
documented in the [API reference](@ref) for advanced users, but their
interfaces may change more freely than the exported simulation and export API.
