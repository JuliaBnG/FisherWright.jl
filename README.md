# FisherWright.jl

FisherWright.jl simulates diploid Wright-Fisher populations forward in time, to
generate whole-genome founder populations for breeding and quantitative-genetics
simulations. Starting from no variation, it applies random mating between two
sexes, Poisson crossovers and recurrent mutation for a chosen number of
generations, neutrally or with selection at given sites. Run long enough (see [Recommended Use](#recommended-use)),
the population reaches mutation-drift equilibrium. The segregating variants are
exported as a `BitMatrix` with a linkage map, as a targeted SNP-chip panel, or as
[BnGStructs.jl](https://github.com/JuliaBnG/BnGStructs.jl) haplotypes.

[![Build Status](https://github.com/JuliaBnG/FisherWright.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/FisherWright.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaBnG/FisherWright.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaBnG/FisherWright.jl)

The complete manual, API reference and population genetics validation are
available at
[juliabng.github.io/FisherWright](https://juliabng.github.io/FisherWright/).

To build the static documentation for hosting:

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs --startup-file=no docs/make.jl
```

The generated site is in `docs/build/`.

## Features
- Forward simulation of a diploid population with two sexes or of
  hermaphrodites, any number of autosomes and genomes up to 2³² − 1 bp
- Random mating, or controlled family sizes with one sire per dam
  (equal, multinomial or overdispersed families; breeder subsets and
  the nested sire–dam design of livestock breeding); hermaphrodites with
  random or given selfing; avoidance of full-sib, half-sib or cousin
  matings; or any scheme written as a function of the cohort (IDs,
  relatives, fitness, phenotypes, haplotypes)
- Pedigree of the last k generations, with running IDs that continue
  across runs
- Constant census size or a given trajectory (bottleneck, growth, any
  vector or function of the generation), and a binomial or fixed sex
  ratio
- Runs can start from a given population and be continued, with or
  without new mutations
- Poisson crossovers, uniform within a chromosome, with a genetic length
  per chromosome (`morgan`, or `M` bp per Morgan) and free assortment
  between chromosomes, and recurrent neutral mutation
- Optional extraction of fixed mutations during the run, which keeps memory
  bounded over long burn-ins
- Selection at an explicit table of sites with given effects (`s`, `h`,
  or drawn from a distribution of fitness effects), or on a quantitative
  trait (stabilizing, directional or any function of the phenotype),
  or on several traits with pleiotropic effects at shared loci;
  the same or different (sex-specific) fitness in males and females;
  soft and within each sex; the causal variants and their effects can be
  exported
- Effective size realized by the family sizes of a run
- Export to a 1-bit-per-allele `BitArray` with a linkage map (`DataFrame`), to a
  targeted SNP-chip panel, or to `BnGStructs.Haplotype` with a
  `BnGStructs.VariantMap`
- A versioned binary file format for results (`write_result`,
  `read_result`), readable by every later release
- Validated against neutral theory, msprime, independent simulations
  of the same model under selection, and matched SLiM models of the
  mating schemes and of linked selection (see the manual)

Not modelled: hard selection and an emergent census size, overlapping
generations, population structure or migration, crossover interference,
and recombination rates that vary within a chromosome.

## Main Functions

### `fisher_wright`

Simulates a Fisher-Wright population of a given size and number of generations,
with specified chromosome lengths and mutation rates. Returns simulated
population data with mutations and recombination events.

Keyword arguments:

| Keyword | Default | Meaning |
|---|---|---|
| `M` | `1e8` | base pairs per Morgan (recombination only) |
| `morgan` | `nothing` | genetic length of each chromosome in Morgans; replaces `chr ./ M` |
| `mut_base` | `1e8` | base pairs per unit of `mr` (mutation only) |
| `result` | `false` | return a `FisherWrightResult` and extract fixed positions |
| `fixation_interval` | `1` | generations between fixation scans, when `result = true` |
| `mutate` | `true` | generate new mutations; `false` leaves drift and recombination only |
| `init` | `nothing` | starting population: `2ne` haplotypes or a `FisherWrightResult` to continue from |
| `popsize` | `nothing` | census size per generation: a vector of `nt` sizes or a function of the generation, for example `bottleneck(100, 10, 21, 5)` or `growth(50, 0.05)` |
| `mating` | `RandomMating()` | `RandomMating(nmale = m)` or `RandomMating(sexratio = p)` fixes the number or proportion of males; `FamilySize(...)` and `Hierarchical(s, d)` control family sizes; `Monoecious(selfing = s)` for hermaphrodites; `AvoidRelatives(inner; level)` avoids related mates; `CustomMating(f)` leaves the matings to a function |
| `fitness` | `Neutral()` | `SiteFitness(loci)` selects on the sites of a `FitnessLoci` table, `TraitFitness(loci; σe, Σe, fitness)` on one or more quantitative traits; needs `result = true` |
| `pedigree` | `0` | generations of matings to record in the result's `pedigree`; needs `result = true` |
| `verbose` | `false` | report progress (mean muts/hap, substitutions, elapsed, ETA) to `stderr` |
| `progress_interval` | `100` | generations between progress reports, when `verbose = true` |
| `rng` | `Random.default_rng()` | random number generator; the realization depends only on its state |

Every haplotype is a sorted vector of unique `UInt32` positions, and that
invariant holds for the returned population.

### `muts2bitarray`

Converts a vector of mutation sets (haplotypes) and chromosome breakpoints into
a `BitArray` of haplotypes and a linkage map (`DataFrame`). Supports random
flipping of alleles and removal of fixed loci.

### `quickHap`

Synthetic haplotypes at independent loci, as a `BnGStructs.Haplotype`.

### `write_result`, `read_result`

Save a `FisherWrightResult` to a file and read it back, for example a
burn-in that later sessions continue.

## Installation

Add the package to your Julia environment:

```julia
using Pkg
Pkg.add("FisherWright")
```

## Usage Example
```julia
using FisherWright

# Simulate a population
muts, cbp = fisher_wright(100, 1000, [1_000_000, 1_000_000], 1.0)

# Convert mutations to bit array and linkage map
xy, lmp = muts2bitarray(muts, cbp; flip = true)
```

## Recommended Use

The population starts with no variation, so it is only at mutation-drift
equilibrium after enough generations. Diversity approaches equilibrium at rate
about `1/(2ne)` per generation, so `nt = 10ne` gets within about 1% of
equilibrium heterozygosity. Shorter runs give a population that is still
gaining diversity: `ne = 2000, nt = 2000`, for example, is far from equilibrium.

```julia
using FisherWright

ne = 200
chr = fill(10_000_000, 10)  # 10 chromosomes of 10 Mbp

# Burn in for 10ne generations and track fixed mutations separately
res = fisher_wright(ne, 10ne, chr, 1.0; result = true)

# Expected neutral values at equilibrium, for a quick check:
# θ = 4ne·μ·L with μ = mr / mut_base per bp per generation (two sexes with
# binomial counts make Ne about ne - 0.5, a difference too small to matter here)
θ = 4ne * (1.0 / 1e8) * sum(chr)                  # Σ2pq ≈ θ
S = θ * sum(1 / i for i in 1:2ne-1)               # segregating sites ≈ θ·aₙ
```

- Use `result = true` for long runs. Fixed mutations are then moved out of the
  haplotypes into `substitutions`, which keeps memory bounded.
- `mr` is the mutation rate per `mut_base` bp per generation, and `M` sets
  recombination only, so the two can be changed independently.
- Pass `rng = Xoshiro(seed)` for a reproducible run. The realization
  depends only on the state of `rng`, not on the thread count, for a
  fixed Julia version and package environment: variate algorithms can
  change between Julia releases.
- For whole-genome populations where only chip markers are needed, use
  `extract_chip_bitarray` or `to_haplotype(result, chip_positions)` rather
  than building the full `BitMatrix`.

## Benchmark

Run the reproducible benchmark script from the package directory:

```bash
julia --project=. bench/benchmark_phase3.jl
```

Optional parameters:

```bash
julia --project=. bench/benchmark_phase3.jl 120 200 100000000,100000000 0.5 5
```

Arguments are interpreted as:

1. `ne`
2. `nt`
3. comma-separated chromosome lengths
4. `mr`
5. repetitions
6. seed

The script reports wall time, cumulative allocation, GC fraction and peak RSS,
and prints the thread count it ran with. Pin threads with `julia -t N` when
comparing against another simulator: cumulative allocation is not a memory
footprint, and wall time depends directly on the thread count.

## Changes in v0.3.0

**Bug fix.** `recombine` skipped the last position of each parental haplotype
inside the crossover loop, emitting it only through the trailing append. At
29 × 100 Mbp with 30 000 mutations per haplotype this made 46% of meioses
inherit the wrong positions and left 11% of offspring haplotypes unsorted,
breaking the sorted-and-unique invariant the rest of the package relies on.
Simulation output changes accordingly.

**Behaviour change.** `fisher_wright` is now silent by default; pass
`verbose = true` for the old progress output.

**API.** `M` no longer scales the mutation rate — use the new `mut_base`
keyword for that. Added allocation-free forms `merge_sorted!`, `cobp!` and
`random_mate!` for use in hot loops.

**Performance.** Per-generation working storage is allocated once and reused,
and the fixation scan uses sorted merges instead of hash sets. Measured on
29 × 100 Mbp, `mr = 1.0`, `result = true`:

| | before | after |
|---|---|---|
| ne=250, nt=200, 12 threads | 1.44 s, 3.12 GiB, 55% GC | 0.27 s, 0.07 GiB, 0.4% GC |
| ne=500, nt=200, 12 threads | 2.37 s, 6.04 GiB, 65% GC | 0.59 s, 0.15 GiB, 7.8% GC |
| ne=250, nt=200, 1 thread | 5.83 s, 3.08 GiB, 70% GC | 0.95 s, 0.07 GiB, 0.5% GC |

The meiosis inner loop (`cobp!` plus `recombine`) no longer allocates at all in
steady state.

## Changes in v0.3.3

Added direct, allocation-conscious extraction of selected chip coordinates:
`extract_chip_bitarray` returns a `BitMatrix`, while
`to_haplotype(result, chip_positions)` returns a `BnGStructs.Haplotype`.
Coordinates must be sorted and unique. A `LocusSet` overload supports
subsetting a shared coordinate vector.

## Changes in v0.3.5

**Bug fix.** `fisher_wright` added each generation's mutations to the parents
before mating, so the population it returned had gone through one round of
reproduction since its last mutations and was missing its newest singletons.
At equilibrium this gave about 7-8% fewer segregating sites than Watterson's
θ·aₙ (a loss of about 0.46θ) and than msprime's discrete-time Wright-Fisher
model; heterozygosity (Σ2pq) was barely affected. Each generation now mates
and recombines first, then mutates the offspring. Segregating sites now match
θ·aₙ within replicate error.

**Simulation output changes.** Populations now contain more rare variants, and
results for a given seed differ from v0.3.4. On the 10 × 100 Mbp benchmark
(`ne = 2000`, `nt = 2000`), the number of segregating sites went from 528,369
to 570,905, against 569,235 from msprime.

**Tests.** Added a check against neutral theory: after `20ne` generations,
segregating sites must match θ·aₙ and Σ2pq must match θ. The old code fails it.

## Changes in v0.3.7

**Bug fix (Thread-safety).** Resolved a silent data race in multithreaded `BitMatrix` construction in `muts2bitarray` and `extract_chip_bitarray`. Because Julia's `BitArray` packs elements column-major in 64-bit (`UInt64`) words, any locus count with `nlc % 64 != 0` shares the boundary word between column $i$ and column $i+1$. Multiple threads setting bits in adjacent columns under `Threads.@threads` previously performed non-atomic read-modify-writes, causing genotype bits to be overwritten and lost. Replaced with atomic bitwise OR (`Core.Intrinsics.atomic_pointermodify`) directly on chunks, guaranteeing deterministic, lossless exports with zero heap allocations and full thread scaling.

**Tests.** Added a regression testset verifying bit-for-bit equivalence between serial and threaded export across non-64-aligned locus counts ($nlc = 65$) over repeated iterations.

## Changes in v0.3.8

**Bug fix (chromosome-boundary assortment).** The independent-assortment
breakpoint between chromosomes is now placed at the first coordinate of the
next chromosome (`cbp[i] + 1`), rather than at the preceding chromosome's
terminal coordinate. This prevents a boundary marker from being assigned to
the wrong segment during recombination. A regression test verifies that the
terminal locus of one chromosome and the first locus of the next segregate
into complementary offspring haplotypes.

**Validation and reproducibility.** CI now exercises the threaded code path
with four Julia threads. The repository includes reproducible benchmark drivers
and captured validation output for the 1 Gb transient comparison and founder
simulation example.

## Changes in v0.3.9 and v0.3.10

**Robustness.** Property-based tests (Supposition.jl) check the
documented contracts against naive references. `muts2bitarray` throws
`ArgumentError` for more than 127 chromosomes and for positions beyond
the last chromosome end; `RecombinationMap` requires strictly
increasing, positive chromosome ends. Sexes are rejection-sampled
without a retry cap: with `ne = 2` the old 10-try cap threw about once
per 1,024 generations. Runs that completed before draw the same numbers.

## Changes in v0.4.0

**Bug fixes. Seeded output changes relative to v0.3.10.**

- `cobp!` draws the crossovers of an interval ending at `e` from
  `p+1:e`, where `p` is the previous interval end. It drew from `p:e-1`,
  so the gap before each interval's last position was never sampled and
  a crossover could fall on the chromosome start. The bias was about
  1/L.
- With `result = true`, a new mutation at a substituted (fixed and
  removed) site has no effect, as with `result = false`. Before, it
  re-entered the active haplotypes as a spurious segregating site.

**Reproducibility.** New `rng` keyword. Each meiosis and haplotype slot
owns a `Xoshiro` stream seeded once from `rng`, and mating is drawn
serially from `rng`. A realization therefore depends only on the state
of `rng`: it is the same for any thread count, for both return modes and
for any `fixation_interval`.

## Changes in v0.4.1

Documentation themes (JuliaBnG Sepia and Black). The simulation is
unchanged.

## Changes in v0.4.2

**Performance.** New mutations are inserted into each haplotype by block
copies between insertion points (`_insert_sorted!`) instead of an
element-by-element merge. The position searches in recombination and
mutation insertion use interpolation search (`_isearchfirst`), since
positions are close to uniform along the genome, and finish with a
binary search. On the 1 Gb benchmark (`ne = 2000`, `nt = 2000`, 10 × 100
Mbp, 8 threads, AMD Ryzen 9 5950X), simulation time went from 61.3 s to
56.9 s, 7% faster (means of three interleaved runs).

**Seeded output is unchanged.** New tests pin the v0.4.0 realizations
with golden-output digests, and the new kernels are tested against
`searchsortedfirst` and `merge_sorted!`.

**Documentation.** The reproducibility guarantee holds for a fixed Julia
version and package environment: `randn` differs between Julia 1.10 and
1.13, which changes Poisson draws with mean ≥ 6 and hence the
realization for a given seed.

**Validation.** Theoretical expectations use the effective size of the
two-sex model, Ne = 99.48 for N = 100, instead of N, and msprime's
diversity and SFS counts are scaled by Ne/N. In 1,000-replicate runs,
Σ2pq was 396.0 ± 0.6 against an expected 395.9. The validation and
benchmark logs were rerun at v0.4.1, and the power check of the CI tests
now also runs the chromosome-boundary test.

## Changes in v0.5.0

**Demography and mating.** The census size, the sex ratio and the
starting population can now be set. The generation is still discrete and
non-overlapping, mating is random between the sexes with sires and dams
drawn with replacement, and the model is neutral.

- `popsize`: census size per generation, as a vector of `nt` sizes or a
  function of the generation. Generation `g` draws the parents of its
  `N_g` offspring from the `N_{g-1}` individuals before it, so drift
  follows the number of parents. `bottleneck(n, nb, start, duration)`
  and `growth(n0, r; start, limit)` build common trajectories.
- `mating = RandomMating(nmale = m)` or `RandomMating(sexratio = p)`:
  exactly `m`, or `round(p * N)`, males in every parental cohort, drawn
  as a uniformly random subset. The default, `RandomMating()`, keeps
  binomial sex.
- `init`: start from given haplotypes or from a `FisherWrightResult`,
  for example to burn in once and continue in several scenarios.
- `mutate = false`: drift and recombination only, for checks such as
  the decay of heterozygosity from a population given by `init`.

**Results.** `FisherWrightResult` has two new fields: `census`, the
sizes N₀, …, N_nt, and `males`, the realized number of males among the
parents of each generation. A run continued from a result appends to its
record. `effective_size(result, generations)` returns Wright's
4N_mN_f/(N_m + N_f) for one generation or its harmonic mean over
several. The three-argument constructor still works; code that
constructs or destructures the struct by all its fields needs the two
new ones.

**Seeded output is unchanged** for calls that use none of the new
keywords: the golden-output digests of v0.4.0 pass. Internally the
generation loop moved into a typed engine (`src/simulate.jl`) with
buffers sized for the largest cohort of the run.

**Validation.** New CI tests compare heterozygosity with the exact
two-sex recursion: decay without mutation, through a bottleneck and an
expansion, and under fixed sex ratios (3♂/17♀, 10♂/10♀, and a
proportion of 0.2 males through a change in size). A new full
validation compares S, Σ2pq, the SFS and LD with msprime's DTWF under
the same demography, after a bottleneck and after an expansion: all 36
checks pass for two seeds at 200 replicates (manual, "Validation:
changing population size").

**Performance.** The default path is within run-to-run noise of v0.4.2.
On the 1 Gb benchmark (`ne = 2000`, `nt = 2000`, 10 × 100 Mbp, 8
threads, AMD Ryzen 9 5950X), simulation time was 57.15 s for v0.4.2 and
57.11 s for v0.5.0 (means of three interleaved runs; single runs ranged
over 0.1 s and 0.4 s), with the same output for the same seed. Peak
resident memory fell from 2,484 to 2,328 MiB. The log is
`docs/src/validation/benchmark-v0.5.0-vs-v0.4.2.txt`.

## Changes in v0.6.0

**Selection.** Selection on an explicit table of sites. The model stays
discrete-generation, dioecious and soft: fitness is evaluated on the
adults, and every offspring draws its sire among the males and its dam
among the females with probability proportional to fitness, with
replacement, so the census size is not affected.

- `FitnessLoci(pos, s; h, μ)`: the selected sites, with the selection
  coefficient and dominance of each derived allele (genotype fitness 1,
  1 + hs, 1 + s, multiplicative over sites) and a mutation rate per
  site. `fitness_loci(chr, nsites; s, h, μ, rng)` draws the sites and,
  if `s` is a distribution, their effects.
- `fitness = SiteFitness(loci)` in `fisher_wright`; `Neutral()` is the
  default. Selection requires `result = true`.
- The derived alleles at the fitness loci are kept apart from the
  neutral mutations, in the new fields `selected_haplotypes` and
  `selected_substitutions`, and pass through meiosis with the same
  crossovers. A starting population with selected alleles is given
  through `init`.
- `to_haplotype(result, loci)` exports the genotypes at the fitness loci
  with their effects; `selected = true` on the other `to_haplotype`
  methods adds the selected alleles to the neutral export.
- `FisherWrightResult` records, per generation, the mean fitness of the
  parents and the mean and variance of log fitness.

**Breaking changes.** `FisherWrightResult` has five new fields
(`selected_haplotypes`, `selected_substitutions`, `mean_fitness`,
`mean_log_fitness`, `var_log_fitness`); the three- and five-argument
constructors still work. New exported names: `FitnessLoci`,
`fitness_loci`, `SiteFitness`, `Neutral`.

**Seeded output is unchanged** for neutral runs: the golden-output
digests of v0.4.0 pass. Under selection, parents are drawn with a Vose
alias table.

**Relation to theory.** The fixation probability of a new allele in
this two-sex model is below that of an ideal Wright-Fisher population of
the same size by terms of order 1/N (5% at N = 50 and 2Ns = 10), from
the binomial sex ratio and from competition within the sex of the
carrier. With many selected sites, selection at unlinked sites lowers
the effective size at each site by about 1/(1 + 4C²), C² the variance of
relative fitness.

**Validation.** The engine is compared with independent naive
simulations of the same model: the fixation probability of one new copy
(eight cases, 200,000 replicates each, |z| ≤ 1.64), and the distribution
of the derived-allele count at selected sites under recurrent mutation,
for deleterious, neutral and beneficial sites (28 comparisons,
|z| ≤ 1.57). Two new CI tests use their reference values (manual,
"Validation: selection").

**Performance.** Neutral runs are as fast as in v0.5.0: on the 1 Gb
benchmark (`ne = 2000`, `nt = 2000`, 10 × 100 Mbp, 8 threads, idle AMD
Ryzen 9 3900X), 63.26 s for v0.5.0 and 63.25 s for v0.6.0 (means of
three interleaved runs; single runs within 0.2 s). Selection costs 4 to
15% more over 500 generations of the same genome, increasing with the
number of selected sites (10³ to 10⁵) and the rate of new selected
alleles (0.01 to 0.1 per haplotype per generation;
`bench/benchmark-selection.jl`). The log is
`docs/src/validation/benchmark-v0.6.0-vs-v0.5.0.txt`.

## Changes in v0.7.0

**Selection on a quantitative trait.** `fitness = TraitFitness(loci;
σe, fitness)` with a `TraitLoci(pos, a; h, μ)` table (or `trait_loci`):
the genotypic value G sums 0, h·a and a over the loci, including fixed
sites; the phenotype is P = G + e with e ~ N(0, σe²) drawn serially from
`rng`; fitness is `Stabilizing(θ, ω²)`, `Directional(β)` or any
non-negative function of P. Selection is soft and within each sex. The
result records the mean and variance of G and P of the parents of every
generation, and `to_haplotype(result, loci)` exports the trait loci with
their effects.

**Family sizes.** `mating = FamilySize(; dams, sires, nsires, ndams,
nmale, sexratio)` keeps the census size of each generation and splits
its offspring among the dams and the dams among the sires (`:equal`,
`:random`, or Dirichlet-multinomial with Gamma(k) weights); every dam
has one sire, so her offspring are full sibs. `nsires` and `ndams`
restrict breeding to random subsets of the cohort, and
`Hierarchical(s, d)` is the nested design of s sires with d dams each
and equal families. With selection, the split weights are multiplied by
fitness (`:equal` splits refuse selection).

**Realized effective size.** The result records the realized family
sizes of every generation (`paternal_pairs`, `maternal_pairs`), and
`effective_size(result)` returns 1/a with
a = (Σm(m − 1) + Σn(n − 1))/(4N(N − 1)), harmonic over generations:
correct for every mating scheme, and including the variance in family
size that selection adds.

**Breaking changes.**

- `effective_size(result)` now returns the effective size realized by
  the family sizes of the run. Under random mating it equals the
  previous value on average but varies between runs;
  `effective_size(result; from = :sexes)` gives the previous value,
  4NmNf/(Nm + Nf).
- `FisherWrightResult` has six new fields (`mean_genotypic`,
  `var_genotypic`, `mean_phenotypic`, `var_phenotypic`,
  `paternal_pairs`, `maternal_pairs`); the existing constructors still
  work.
- New exported names: `TraitLoci`, `trait_loci`, `TraitFitness`,
  `Stabilizing`, `Directional`, `FamilySize`, `Hierarchical`.

**Seeded output is unchanged** for runs with the default mating and no
selection: the golden-output digests of v0.4.0 pass.

**Validation.** Drift under every family-size scheme, against the exact
two-sex recursion with the coalescence probability of the scheme
(closed forms for all splits, breeder subsets included). One generation
of trait selection is tested to be unbiased given the parents, under
random mating and under family sizes. The neutral genetic variance of a
trait matches its exact expectation E[α²](2H_B − H_W) at every
generation (22 checks, |z| ≤ 1.93); the response to directional
selection is reported against the Gaussian prediction (manual,
"Validation: selection").

**Performance.** Runs with the default mating and no selection are as
fast as in v0.6.0: on the 1 Gb benchmark (`ne = 2000`, `nt = 2000`,
10 × 100 Mbp, 8 threads, idle AMD Ryzen 9 3900X), 63.29 s for v0.6.0
and 63.28 s for v0.7.0 (means of three interleaved runs), with the same
output. Over 500 generations of the same genome, selection at sites
costs 7 to 15% more, selection on a trait with 10,000 loci 16%, and
`FamilySize` costs nothing in itself (its runs are faster when families
are variable, because the effective size and so the diversity are
lower). The log is `docs/src/validation/benchmark-v0.7.0-vs-v0.6.0.txt`.

## Changes in v0.7.1

A maintenance release of v0.7.0, with no change to the source or to
seeded output. BnGStructs 0.3 is allowed (`BnGStructs = "0.2.5, 0.3"`),
so that packages on BnGStructs 0.3 can resolve with FisherWright. The
standard libraries have compat entries, the test suite runs Aqua, and
the repository has a `CITATION.cff`.

## Changes in v0.8.0

**Pedigree.** `pedigree = k` records the matings of the last k
generations in the result, as a `DataFrame` of `id`, `sire`, `dam`,
`sex` and `generation`, with the parents of the oldest recorded
generation as founders. Individuals are numbered in birth order across
the run, and a run continued through `init` continues the numbering;
`first_id` links the haplotypes of the result to the pedigree.
Recording does not change the run.

**Hermaphrodites.** `mating = Monoecious()` draws both parents from the
whole cohort, so that an individual selfs at random with probability
1/N (the discrete-time Wright–Fisher model); `Monoecious(selfing = s)`
selfs with probability s and otherwise mates with one of the others.
Under selection both parents are drawn in proportion to fitness. The
result records the proportion of selfed offspring per generation
(`selfing`), and `effective_size` allows for it.

**Avoidance of relatives.** `AvoidRelatives(inner; level, max_tries)`
wraps `RandomMating` or `Monoecious` and redraws the dam while she is a
full sib (`:fullsib`), shares a parent (`:halfsib`) or a parent or
grandparent (`:cousin`) with the sire. After `max_tries` draws the
related pair is kept, and the result counts such offspring per
generation (`related_matings`).

**User-defined mating.** `CustomMating(f; sexes)` calls `f(ctx)` every
generation to fill the mating table, with a `MatingContext` holding
the generation, the cohort's IDs, parents and grandparents, relative
fitness, genotypic values and phenotypes, haplotypes and the random
number generator; with `sexes = true` the function also sets the
sexes. The table is checked after each call.

**Breaking changes.**

- `FisherWrightResult` has four new fields (`selfing`,
  `related_matings`, `first_id`, `pedigree`); the existing constructors
  still work.
- New exported names: `Monoecious`, `AvoidRelatives`, `CustomMating`,
  `MatingContext`.

**Seeded output is unchanged** for the mating schemes and fitness models
of v0.7.0 (checked against v0.7.1 for random mating, fixed males,
`FamilySize`, `Hierarchical`, a bottleneck, site and trait selection),
and the golden-output digests of v0.4.0 pass.

**Validation.** Drift under `Monoecious`, with and without selfing,
against the exact recursion of the heterozygosity within and between
individuals (CI test). Every mating scheme against a matched SLiM model
written independently in Eidos, on SLiM's meiosis and pedigree (manual,
"Validation: mating schemes"): 2 × 19 comparisons, largest |z| 1.77.
Linked selection (background selection and sweeps at 2,000 sites)
against a matched SLiM model (manual, "Validation: linked
selection"): neutral sites, background selection and sweeps alone
agree within noise; with sweeps and background selection together,
FisherWright's π is 0.23% below SLiM's over 14,000 replicates (pooled
z = −2.79, not significant over the comparisons made).

**Performance.** Runs with the default mating and no selection are as
fast as in v0.7.1: on the 1 Gb benchmark (`ne = 2000`, `nt = 2000`,
10 × 100 Mbp, 8 threads, idle AMD Ryzen 9 5950X), 57.60 s for v0.7.1
and 57.72 s for v0.8.0 (means of three interleaved runs, which
overlap), with the same output. The new schemes and pedigree recording
cost at most 2% over random mating. The log is
`docs/src/validation/benchmark-v0.8.0-vs-v0.7.1.txt`.

## Changes in v0.8.1

**BnGStructs 0.4.** BnGStructs 0.4 removed the population size from
`Species` (`Species.nid`): a size belongs to a population, not to a
species. The species method now takes it explicitly,
`fisher_wright(sp, ne, nt[, mr]; kw...)`, and the compat allows
BnGStructs 0.2.5–0.4. The old form `fisher_wright(sp, nt[, mr])` is
deprecated: with BnGStructs 0.3 and a species that has a size it still
runs (with a deprecation warning), with 0.4 it throws an
`ArgumentError` that names the new form. Seeded output is unchanged.

## Changes in v0.8.2

**Per-chromosome genetic lengths.** `fisher_wright` and
`uniform_recombination_map` take `morgan`, the genetic length of each
chromosome in Morgans, in place of the genome-wide `chr ./ M`. The
species method reads `sp.morgan` (BnGStructs ≥ 0.3), so a species
built with its own genetic lengths now recombines with them; `M` or
`morgan` overrides it. `morgan = chr ./ M` gives the same realization
as `M`, and seeded output is unchanged for every species with default
genetic lengths, the named species included.

## Changes in v0.9.0

**Pleiotropy.** A `TraitLoci` table takes a `K × T` matrix of effects:
`T` traits on shared loci, each site with an effect on every trait,
and dominance per site or per site and trait. `trait_loci` draws the
effects from a multivariate distribution, for example
`MvNormal(0, Σa)`, whose covariance becomes the genetic covariance of
the traits. `TraitFitness(loci; σe, Σe, fitness)` adds environmental
deviations with standard deviations `σe` (independent) or covariance
matrix `Σe` (positive semi-definite), and selects on the phenotype
vector: `Stabilizing(θ, Ω)` (Gaussian, `Ω` positive definite, its
off-diagonal elements correlational selection), `Directional(β)`
(the selection gradient), or any function of the vector. New export:
`ntraits`.

**Breaking changes.**

- `mean_genotypic` and `mean_phenotypic` of `FisherWrightResult` are
  `T × ngen` matrices, and `var_genotypic` and `var_phenotypic`
  `T × T × ngen` arrays, the covariance matrix of each generation. A
  run with one trait (or none) has `T = 1`: `mean_genotypic[1, g]` and
  `var_genotypic[1, 1, g]`; linear indexing (`mean_genotypic[g]`)
  still gives generation `g`.
- `TraitLoci.a` and `TraitLoci.h` are `K × T` matrices; `Stabilizing`
  has fields `θ::Vector` and `Ω::Matrix` (was `θ`, `ω²`) and
  `Directional` has `β::Vector`. The constructors with numbers are
  unchanged.
- `TraitFitness` no longer has a field `σe`; it keeps `Σe` and a factor
  `L` with `L L' = Σe`.
- `MatingContext.genotypic` and `phenotypic` are `T × n` matrices
  (`view(ctx.phenotypic, 1, :)` for one trait).
- `to_haplotype(result, loci::TraitLoci)` names the effect columns
  `a1, …, aT` and `h1, …, hT` when `T > 1`; one trait keeps `a`, `h`.

Seeded output is unchanged for one trait: the scalar formulae and the
order of the environmental draws are kept, and runs with
`Directional`, `Stabilizing` and a user function, with and without
environmental noise, give the same populations and records as 0.8.2.
With `T` traits every individual takes `T` standard normal draws per
generation, in trait order.

## Changes in v0.9.1

**Sex-specific fitness.** `fitness = SexSpecific(male, female)` gives
males and females their own fitness model (`Neutral`, `SiteFitness` or
`TraitFitness`): sires are drawn in proportion to the male fitness and
dams to the female fitness. Two `SiteFitness` models share the
positions and mutation rate of their tables and may differ in `s` and
`h`, for sexually antagonistic selection; two `TraitFitness` models
share the trait loci and environmental (co)variance and may differ in
the fitness function; with `Neutral` on one side, selection acts in
one sex only. Sex-limited expression is two traits of one table, the
male fitness on one and the female fitness on the other. Fitness is
evaluated for every parent under both models before the sexes are
drawn, so a run draws its random numbers in the same order as with one
model. Mating schemes with separate sexes drawn by the run only
(`Monoecious` and `CustomMating` are refused). Seeded output of the
other models is unchanged. No breaking changes.

## Changes in v0.10.0

A structural release: the result, the internal state and the export
types are reorganized, without a new model. Seeded output is unchanged:
53 configurations (every fitness model with every mating scheme,
continuation, changing size, genetic lengths, no mutation, no result,
and the export functions) give the same haplotypes, records and
pedigrees as 0.9.1, at 1 and 4 threads.

**New.**

- `GenerationRecords`: the record of every generation (`census`,
  `males`, the fitness and trait records, family sizes, `selfing`,
  `related_matings`), with a keyword constructor that fills the fields
  not given and checks the lengths, and `effective_size` of a record.
- `write_result(path, result)` and `read_result(path)`: a binary result
  format defined byte by byte (little-endian, CRC-32C checksum), with
  a format version; every later release reads the files of every
  earlier version, and the tests read a file written by 0.10.0. It
  replaces `Serialization`, which depends on the Julia and package
  versions.
- `rng` keyword of `to_haplotype` and `muts2bitarray`, for the
  placeholder nucleotides `ref` and `alt`, which came from the global
  generator only; the default is unchanged.

**Breaking changes.**

- The per-generation record of `FisherWrightResult` is one field,
  `records::GenerationRecords`: `result.records.census`,
  `result.records.mean_fitness`, and so on. `result.census` and the
  other old names still work with a deprecation warning, and will be
  removed in 1.0.
- `FisherWrightResult(haps, ends, subs; selected_haplotypes,
  selected_substitutions, records, first_id, pedigree)` is the
  constructor. The positional constructors of 5 to 20 arguments are
  deprecated, and the constructor checks that `selected_haplotypes` has
  one entry per haplotype.
- `to_haplotype(result)` returns a `BnGStructs.VariantMap` instead of a
  `DataFrame`, with the same columns (`morgan` empty), and throws an
  `ArgumentError` (was `ErrorException`) when no locus is polymorphic.
  `muts2bitarray` and the loci tables of `to_haplotype(result, loci)`
  stay `DataFrame`s.
- `quickHap` returns a `BnGStructs.Haplotype` (was `Matrix{Int8}`), with
  the same alleles for the same `rng`; `Int8.(hp)` gives the old
  matrix. It throws an `ArgumentError` (was `ErrorException`) for `maf`
  outside `(0, 0.5)`.
- `quickGT` is removed. Genotypes at independent loci with Beta
  frequencies are a poor stand-in for a population; simulate a small
  historical population with `fisher_wright` (cheap), or sum the two
  haplotypes of each individual from `quickHap` for the same model.
- Requires BnGStructs 0.3 or 0.4 (0.2.5 is dropped, for `VariantMap`).
- Internal: `SimState` keeps the genome and the per-slot streams, and
  groups the rest into `mate` (mating table and the scratch of the
  mating schemes), `fit` (selected stream and fitness) and `rec` (the
  record as it grows); code that reached into `SimState`, such as a
  user fitness model writing `st.w`, uses `st.fit.w`.

## Contributing and support

See [CONTRIBUTING.md](CONTRIBUTING.md). Bug reports and questions go to the
[issue tracker](https://github.com/JuliaBnG/FisherWright.jl/issues).

## License

MIT License. See LICENSE file for details.
