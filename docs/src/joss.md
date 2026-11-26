# Changes since the JOSS paper

The JOSS paper (`paper/paper.md`) describes version 0.4: its validation
ran on v0.4.1 and its benchmark on v0.4.2. The paper is kept as
submitted, but the package has grown since, and a reader who installs it
now gets version 0.11. This page says which statements of the paper
still hold for the current release, which no longer do, what code
written for 0.4 needs to change, and how the paper's numbers can be
reproduced.

## What still holds

**The default model.** A call without the keywords added since 0.4,
such as `fisher_wright(ne, nt, chr, mr)`, runs the model of the paper:
constant census size, two sexes with binomial sex counts, sire and dam
drawn at random with replacement, Poisson crossovers with independent
assortment, neutral mutation.

**Reproducibility, and the realizations themselves.** Fixed seeds give
the realizations of v0.4.0 exactly. The test suite pins them with
golden-output digests (`test/golden.jl`, "Golden output (v0.4.0
realizations)"), in both return modes, and every release since 0.4.2
has passed it. A realization is independent of
the thread count, the fixation interval and whether substitutions are
returned separately, as the paper states, and, as it also states, it can
change with the Julia version where a Poisson mean reaches 6.

**Design.** Haplotypes are sorted `Vector{UInt32}` coordinates;
`recombine`, `merge_sorted!` and `cobp!` are the kernels of the paper,
and the tests still check that they allocate nothing after warm-up.
Mutation insertion has used block copies (`_insert_sorted!`) since 0.4.2,
with the same output. Fixed substitutions are extracted as described.

**Export.** `muts2bitarray` (a `BitMatrix` and a `DataFrame` map),
`extract_chip_bitarray` and `to_haplotype(result, positions)` are
unchanged.

**Validation.** The paper's validation, rerun with v0.11.0 for both of
its seeds, gives the same output as the logs of v0.4.1 that the paper
reports: every statistic, expectation and z-score, for FisherWright and
for the msprime reference, is identical to the last printed digit, and
all 24 checks pass with maximum ``|z|`` of 2.37 (seed 7) and 1.91 (seed
2026), as in the paper. Only the run times differ (another host). The
realizations are identical because the default path still reproduces
v0.4.0, and every Poisson mean in this setting is below 6, so they do
not depend on the Julia version either.

**Performance.** The default path has been benchmarked at every minor
release on the paper's run (N = 2,000, 2,000 generations, 10 × 100 Mb,
8 threads), each against the release before, with interleaved runs on
the same idle host:

| Releases | Host | Mean simulation time (s) |
|---|---|---|
| v0.4.2, v0.5.0 | Ryzen 9 5950X | 57.15, 57.11 |
| v0.5.0, v0.6.0 | Ryzen 9 3900X | 63.26, 63.25 |
| v0.6.0, v0.7.0 | Ryzen 9 3900X | 63.29, 63.28 |
| v0.7.1, v0.8.0 | Ryzen 9 5950X | 57.60, 57.72 |
| v0.4.2, v0.11.0 | Ryzen 9 5950X | 57.40, 57.42 |

The logs are `docs/src/validation/benchmark-*.txt`. The paper's 57.0 s
for v0.4.2 is on the Ryzen 9 5950X. Releases 0.9 to 0.11 left the
output of the default path unchanged (the golden digests), but 0.10.0
reorganized the internal state of a run, so v0.11.0 was compared with
v0.4.2 directly on the paper's host (2026-11-26, Julia 1.13.4): the two
differ by 0.03% in mean simulation time, well inside the spread of
their runs (57.08–57.89 s), and every run gives the same 571,686
polymorphic SNPs and 4 substitutions. The default path is as fast as
in the paper, after all the features added since.

## What no longer holds

**Scope.** The paper says that the package "does not implement
selection, population structure, or non-stationary demography" and that
mating is the random union of a sire and a dam. Since 0.5 it implements
the first and the last of these, and other mating schemes. Each is
opt-in, through a keyword whose default is the model of the paper:

| Release | Added | Keyword |
|---|---|---|
| 0.5.0 | census size per generation (bottleneck, growth, any trajectory); fixed number or proportion of males; start from a given population; runs without mutation | `popsize`, `mating = RandomMating(...)`, `init`, `mutate` |
| 0.6.0 | selection on a table of sites (`s`, `h`, distributions of effects) | `fitness = SiteFitness(...)` |
| 0.7.0 | selection on a quantitative trait (stabilizing, directional, any function); controlled family sizes and the nested sire–dam design; effective size realized by the family sizes | `fitness = TraitFitness(...)`, `mating = FamilySize(...)`, `Hierarchical(...)` |
| 0.8.0 | pedigree of the last k generations; hermaphrodites with selfing; avoidance of relatives; any mating rule as a function | `pedigree`, `Monoecious`, `AvoidRelatives`, `CustomMating` |
| 0.8.2 | a genetic length per chromosome | `morgan` |
| 0.9.0 | several traits on shared loci (pleiotropy) | `TraitLoci` with a matrix of effects |
| 0.9.1 | fitness that differs between the sexes | `SexSpecific` |
| 0.10.0 | a versioned file format for results | `write_result`, `read_result` |
| 0.11.0 | assortative mating by phenotype | `mating = Assortative(ρ)` |

Each model and mating scheme was validated before release, against
theory, against independent simulations of the same model, or against
matched SLiM or msprime models; the manual has a page for each group. Population structure and
migration are still not modelled, nor overlapping generations.

**The statement on `BnGsim.jl`.** The paper calls the breeding-scheme
simulator `BnGsim.jl` "in development". It has since been registered
and runs the breeding loop on `FisherWright.jl` founders.

## Code written for 0.4

Most code written for 0.4 runs unchanged. These changes can break it:

- **The record of a result** (0.10). Fields added from 0.5 to 0.9 that
  hold one entry per generation (`census`, `males`, the fitness and
  trait records, family sizes) moved into `result.records`. The fields
  of 0.4 (`active_haplotypes`, `chromosome_ends`, `substitutions`) are
  where they were. The positional constructors are deprecated; use
  `FisherWrightResult(haps, ends, subs)`, which 0.4 also had.
- **`to_haplotype(result)`** (0.10) returns a `BnGStructs.VariantMap`
  instead of a `DataFrame`, with the same columns; the haplotypes are
  the same.
- **`quickGT` and `quickHap`** (0.10) are removed: genotypes at
  independent loci are a poor stand-in for a population, and a small
  historical population simulated with `fisher_wright` takes a fraction
  of a second.
- **`fisher_wright(sp, nt[, mr])`**, the species method (0.8.1), takes
  the population size: `fisher_wright(sp, ne, nt[, mr])`, since
  `Species` no longer carries one (BnGStructs 0.4).
- **`effective_size(result)`** (0.7) returns the size realized by the
  family sizes; `from = :sexes` gives the value of 0.5 and 0.6.
- **BnGStructs** 0.3 or 0.4 is required (0.10).

The README lists every change by release.

## Reproducing the paper's numbers

Validation (about a minute per seed with 8 threads; needs `uv` for the
msprime reference):

```bash
julia -t 8 --project=. docs/src/validation/scripts/validate-popgen-theory.jl 50 100 2026
julia -t 8 --project=. docs/src/validation/scripts/validate-popgen-theory.jl 50 100 7
```

Benchmark (the paper's run; compare versions with interleaved runs on
an idle host):

```bash
julia --startup-file=no -t 8 --project=. bench/simulate-fisher-wright.jl 2000 2000 10 100000000 1.0 42
```

The paper's logs are `docs/src/validation/validate-popgen-seed2026.txt`,
`validate-popgen-seed7.txt` and
`benchmark-v0.5.0-vs-v0.4.2.txt`; the reruns with v0.11.0 are
`validate-popgen-v0.11.0-seed2026.txt`,
`validate-popgen-v0.11.0-seed7.txt` and
`benchmark-v0.11.0-vs-v0.4.2.txt`, in the same directory:

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-popgen-v0.11.0-seed2026.txt"), String)))
```
