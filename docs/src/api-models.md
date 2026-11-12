# API reference: population models

The models passed to [`fisher_wright`](@ref) through its `mating`,
`popsize` and `fitness` keywords.

## Mating

```@docs
RandomMating
FamilySize
Hierarchical
Monoecious
AvoidRelatives
CustomMating
MatingContext
```

## Census size

```@docs
bottleneck
growth
```

## Selection at sites

```@docs
SiteFitness
FitnessLoci
fitness_loci
Neutral
```

## Selection on a trait

```@docs
TraitFitness
TraitLoci
trait_loci
Stabilizing
Directional
```
