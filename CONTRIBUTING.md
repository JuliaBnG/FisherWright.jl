# Contributing to FisherWright.jl

Thank you for your interest in FisherWright.jl. The package is part of the
[JuliaBnG](https://github.com/JuliaBnG) packages for breeding and
quantitative-genetics simulation in Julia.

## Getting help and reporting problems

- **Questions and bug reports:** open an issue on the
  [issue tracker](https://github.com/JuliaBnG/FisherWright.jl/issues). For a
  bug, include the Julia and FisherWright versions (`using Pkg; Pkg.status()`),
  the number of threads, and a minimal example with a fixed seed. Results are
  reproducible for a given seed only with the same thread count.
- **Suspected statistical errors:** if simulated diversity, the site frequency
  spectrum, LD or crossover counts disagree with theory or with another
  simulator, say which statistic, which reference and how many replicates. The
  manual's "Population genetics validation" page describes the references and
  statistics already used.
- **Feature requests:** open an issue describing the use case. The package
  deliberately covers neutral, constant-size, randomly mating populations
  (see "Scope" below).

## Contributing code or documentation

1. Open an issue first for anything beyond a small fix, so the approach can be
   agreed before you spend time on it.
2. Fork the repository and create a branch from `main`.
3. Make the change, with tests in `test/runtests.jl`. Run them with
   `julia -t 4 --project -e 'using Pkg; Pkg.test()'`.
4. If you change the generation loop (`src/fwp.jl`), `recombine`, `cobp!`, the
   recombination map, `random_mate!` or `muts2bitarray`, rerun the full
   validation as described under "When to rerun" on the validation page, and
   include the new output files.
5. Update the docstrings, and the manual in `docs/src/` where relevant.
6. Open a pull request that describes what changed and why.

### Conventions

- Haplotypes are sorted vectors of unique `UInt32` positions. Every function
  that returns haplotypes must preserve this invariant.
- Hot loops (meiosis, mutation, merging) must not allocate in steady state.
  The test suite checks this for `cobp!`, `recombine` and `merge_sorted!`.
- New low-level random functions should accept an `rng` keyword, as `recombine`,
  `cobp!` and `random_mate!` do. `fisher_wright` itself uses the global RNG
  (`Random.seed!`).
- Statistical tests should be shown to fail for a known error, not only to pass
  for the correct model (see "Power of the CI tests" in the manual).

## Scope

FisherWright.jl generates neutral founder populations. Selection, changes in
population size, population structure and spatial models are out of scope.
SLiM, fwdpy11 and msprime cover these well. Proposals that extend the model
are welcome for discussion, but each addition must come with validation
against theory or an independent simulator.

## Governance and maintenance

The package is maintained by Xijiang Yu (Norwegian University of Life
Sciences). The maintainer reviews issues and pull requests, usually within two
weeks, and makes the final decision on changes and releases. Releases follow
semantic versioning and are registered in the Julia General Registry. User-visible
changes are listed in the README.

## Code of conduct

Be respectful and constructive. Harassment or personal attacks are not
tolerated. The maintainer may close discussions or block participants who
behave inappropriately.

## License

By contributing, you agree that your contributions are licensed under the MIT
License of this repository.
