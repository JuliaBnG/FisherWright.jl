# FisherWright.jl contributor guidance

## Commands

This is a Julia 1.10+ package. Bootstrap the package environment with:

```bash
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate()'
```

Run the complete test suite:

```bash
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

Tests are collected in `test/runtests.jl` and are not selector-enabled. For a
single assertion while iterating, load the package and execute the focused
check directly, for example:

```bash
julia --project=. --startup-file=no -e 'using Test, FisherWright; @test FisherWright.merge_sorted(UInt32[1, 3], UInt32[2, 3]) == UInt32[1, 2, 3]'
```

Build the documentation (the docs environment develops the local checkout):

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs --startup-file=no docs/make.jl
```

Run the reproducible performance benchmark with:

```bash
julia --project=. bench/benchmark_phase3.jl
```

There is no lint command configured in the repository.

## Architecture

`src/FisherWright.jl` is the module entry point and includes the implementation
files. `fwp.jl` owns the main `fisher_wright` loop: each generation adds
mutations, samples mating pairs, recombines each parent pair into offspring,
then swaps parent/offspring storage. The mutation and meiosis loops are
threaded; per-generation scratch buffers are allocated before the loop and
reused.

The simulation's canonical representation is a vector of haplotypes, where
each haplotype is a sorted, unique `Vector{UInt32}` of genome-wide mutation
positions. Chromosome locations are cumulative `UInt32` endpoints (`cbp`), not
per-chromosome coordinates. `recombine.jl` provides recombination-map
construction and crossover sampling; `random-mate.jl` supplies parent pairs;
`merge-sorted.jl` maintains mutation sets; and `result.jl` removes fixed
mutations into `FisherWrightResult.substitutions` when `result=true`.

`muts2bitarray.jl` and `to_haplotype` convert sparse simulation output to a
locus-by-haplotype `BitMatrix` plus linkage-map `DataFrame`. The default export
removes monomorphic loci; structured results require `include_fixed=true` to
put substitutions back into dense output. `quick-genotypes.jl` is independent
synthetic genotype/haplotype generation rather than part of the simulation
pipeline.

## Repository conventions

- Preserve the sorted-and-unique haplotype invariant through every mutation,
  recombination, fixation, and conversion path. In particular, crossover
  positions and inputs to `merge_sorted!`/`recombine` must be sorted.
- Prefer the in-place `!` helpers in simulation hot paths. Their destination
  buffers are emptied and reused; `merge_sorted!` must not receive a destination
  that aliases either input.
- Keep mutation positions and chromosome endpoints as `UInt32`, and reject
  genome sizes that cannot fit in that coordinate system.
- `M` controls recombination only; `mut_base` controls the mutation-rate basis
  only. Do not couple changes to these parameters.
- Preserve threaded behavior. Explicit RNGs are supported by utility helpers,
  while `fisher_wright` is reproducible with seeding only when
  `Threads.nthreads()` is fixed.
- Public API changes belong in the export list in `src/FisherWright.jl`, must
  have docstrings, and should be added to `docs/src/api.md`; docs run with
  `checkdocs = :exports`.
- Tests intentionally import internal helpers explicitly from `FisherWright`
  to cover low-level invariants. Use deterministic `MersenneTwister` instances
  for tests whose result depends on random draws; keep statistically sized
  simulation checks tolerant rather than exact.
