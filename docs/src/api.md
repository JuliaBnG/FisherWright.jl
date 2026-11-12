# API reference

## Package

```@docs
FisherWright.FisherWright
```

## Simulation and results

Mating schemes, demographic trajectories and fitness models are on
[API reference: population models](@ref).

```@docs
fisher_wright
FisherWrightResult
effective_size
to_haplotype
```

## Dense export

```@docs
muts2bitarray
extract_chip_bitarray
```

## Recombination

```@docs
RecombinationMap
uniform_recombination_map
```

## Synthetic data

```@docs
quickGT
quickHap
```

## Advanced helpers

```@docs
FisherWright.merge_sorted
FisherWright.merge_sorted!
FisherWright.recombine
FisherWright._recombine!
FisherWright.cobp
FisherWright.cobp!
FisherWright.random_mate
FisherWright.random_mate!
```
