# Simulation

## `fisher_wright`

```julia
fisher_wright(ne, nt, chr, mr; M=1e8, mut_base=1e8,
              result=false, fixation_interval=1, mutate=true,
              init=nothing, popsize=nothing, mating=RandomMating(),
              fitness=Neutral(), verbose=false, progress_interval=100)
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
| `fitness` | `Neutral()` | Fitness model; `SiteFitness(loci)` selects on the sites of a `FitnessLoci` table (needs `result=true`). |
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
sizes = g -> round(Int, 20 * 1.05^g)   # 5% growth per generation
haplotypes, _ = fisher_wright(20, 20, [100_000, 100_000], 1.0;
                              popsize=sizes)
length(haplotypes) == 2sizes(20)
```

`bottleneck` and `growth` build the two common trajectories as
functions of the generation. `bottleneck(n, nb, start, duration)` gives
`nb` individuals in generations `start` to `start + duration - 1` and
`n` otherwise (`recovery` sets a different size afterwards).
`growth(n0, r)` gives ``N_g = n_0(1 + r)^g``, rounded; `start` delays
it, `limit` is the ceiling of a growing and the floor of a declining
population, and a negative `r` is a decline.

```@example simulation
run = fisher_wright(40, 20, [100_000, 100_000], 1.0;
                    popsize=bottleneck(40, 10, 6, 5), result=true)
run.census[5:12]   # N_4, …, N_11
```

```@example simulation
run = fisher_wright(20, 40, [100_000, 100_000], 1.0;
                    popsize=growth(20, 0.05; limit=60), result=true)
run.census[[1, 2, 11, 41]]   # N_0, N_1, N_10, N_40
```

Both return plain functions, so a trajectory can also be built as a
vector and combined with others: `growth(20, 0.05).(1:40)`.

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

## Family sizes

`RandomMating` lets every offspring choose its sire and its dam
independently, so family sizes are binomial, close to Poisson.
`mating = FamilySize(; dams, sires, nmale, sexratio)` controls them
instead. The census size of each generation stays as given; its
offspring are split among the dams, and every dam is mated to one sire,
so all her offspring are full sibs. The dams are split among the sires.
Sex is assigned as by `RandomMating` with the same `nmale` or
`sexratio`.

| `dams`, `sires` | Split |
|---|---|
| `:equal` | as equal as possible; the parents that get one more are drawn at random |
| `:random` | uniform multinomial (the default) |
| `k > 0` | multinomial with Gamma(``k``, 1) weights per parent and generation: Dirichlet-multinomial, close to negative binomial; smaller ``k``, more variable |

```@example simulation
# 5 sires with 7 dams each, and equal dam families
equal = FamilySize(dams=:equal, sires=:equal, nmale=5)
haplotypes, _ = fisher_wright(40, 50, [100_000, 100_000], 1.0;
                              mating=equal)
length(haplotypes)
```

Two offspring share their dam with probability
``P_f = E[\sum n_d(n_d - 1)]/(N(N - 1))``. They share their sire if they
share their dam, or else if their dams share a sire, which happens with
probability ``\pi = E[\sum c_s(c_s - 1)]/(N_f(N_f - 1))`` for ``c_s``
dams per sire. Two genes in different offspring therefore come from the
same parent with probability ``a = (P_f + P_m)/4``, with
``P_m = P_f + (1 - P_f)\pi``, and the recursion of
[Effective size](@ref) applies with this ``a``; the effective size is
about ``1/a + 1/2`` (``1/a = N`` for ``N/2`` males and ``N/2`` females
under random mating). For the three splits of ``n`` items among ``m``
parents, the probability that two items share a parent is
``(r(q + 1)q + (m - r)q(q - 1))/(n(n - 1))`` with ``n = qm + r`` for
`:equal`, ``1/m`` for `:random`, and ``(k + 1)/(mk + 1)`` for Gamma
weights. In the example above this gives ``N_e \approx 21.8``,
against ``18.0`` for `RandomMating(nmale=5)`: equal dam families raise
it, and a few sires keep it low.

**Breeders.** `nsires` and `ndams` restrict breeding to random subsets
of the cohort: in every generation `nsires` of the males and `ndams` of
the females are the sires and dams, and the others have no offspring.
`Hierarchical(s, d)` is the nested design of livestock breeding, ``s``
sires each mated to ``d`` dams with equal families; it is
`FamilySize(dams=:equal, sires=:equal, nsires=s, ndams=s*d)`. The
coalescence probabilities above then use the numbers of breeders.

```@example simulation
# 4 sires × 5 dams out of 100 individuals per generation
design = Hierarchical(4, 5)
run = fisher_wright(100, 30, [100_000, 100_000], 1.0; mating=design,
                    result=true)
(run.census[end], run.males[end])
```

Here ``N_e \approx 14.6``, against ``99.5`` for random mating of all
100. With binomial sex the sexes are redrawn until there are enough
sires and dams, as they are redrawn until both sexes are present; with
`nmale` or `sexratio`, every parental cohort must hold the breeders,
which is checked before the run. `males` still records the number of
males among all parents.

`effective_size(result)` returns the effective size realized by the
family sizes of the run (see [Census record and effective size](@ref)),
which for `Hierarchical` equals the value above exactly. Neither
`FamilySize` nor `Hierarchical` can be combined with selection yet.

## Selection

By default the model is neutral. `fitness = SiteFitness(loci)` puts the
sites of a `FitnessLoci` table under selection. The table is fixed for
the run: each site has a position, a selection coefficient ``s`` of its
derived allele and a dominance ``h``, and the table has one mutation
rate ``\mu`` per site per haplotype per generation.

```@example simulation
loci = FitnessLoci([20_000, 60_000, 150_000], [0.05, -0.02, 0.1];
                   h=[0.5, 0.5, 0.0], μ=1e-4)
(loci.pos, loci.s, loci.h)
```

`fitness_loci` draws the sites uniformly over the genome, and the
coefficients from a distribution of fitness effects if `s` is a
distribution from Distributions.jl; `h` may be a function of ``s``:

```@example simulation
using Distributions, Random

dfe = fitness_loci([100_000, 100_000], 50; s=-Gamma(0.3, 0.1),
                   h=s -> 0.25, μ=1e-4, rng=Xoshiro(1))
(length(dfe), extrema(dfe.s))
```

**Model.** A genotype at one site has fitness ``1``, ``1 + hs`` or
``1 + s`` with zero, one or two copies of the derived allele, and the
fitness of an individual is the product over sites. Selection is soft
and acts within each sex: fitness is evaluated on the adults, and every
offspring draws its sire among the males and its dam among the females
with probability proportional to fitness, with replacement. The census
size does not change, and only fitness relative to the other parents of
the same sex matters. Sites are biallelic without back-mutation.

```@example simulation
run = fisher_wright(50, 200, [100_000, 100_000], 1.0;
                    fitness=SiteFitness(loci), result=true,
                    rng=Xoshiro(3))
(run.selected_substitutions, count(!isempty, run.selected_haplotypes))
```

The derived alleles at the fitness loci are kept apart from the neutral
mutations. They are returned in `selected_haplotypes`, one sorted vector
of positions per haplotype, and in `selected_substitutions` once fixed;
`active_haplotypes` holds neutral mutations only, and a neutral mutation
at the coordinate of a fitness locus has no effect on fitness. The two
recombine as one haplotype. Removing a fixed site is exact, because it
multiplies every fitness by the same factor. Selection therefore needs
`result=true`.

A selected allele that is present at the start is given through `init`,
as a result with `selected_haplotypes`; this is how a single new copy is
followed:

```@example simulation
chr = [100_000, 100_000]
one = FitnessLoci([20_000], 0.1)
sel = [UInt32[] for _ in 1:100]
push!(sel[1], 20_000)                    # one copy, on haplotype 1
start = FisherWrightResult([UInt32[] for _ in 1:100], UInt32.(cumsum(chr)),
                           UInt32[]; selected_haplotypes=sel)
later = fisher_wright(50, 20, chr, 1.0; init=start, mutate=false,
                      fitness=SiteFitness(one), result=true)
count(!isempty, later.selected_haplotypes)   # copies after 20 generations
```

**Relation to theory.** The fixation probability of a new allele is a
few percent below that of an ideal Wright–Fisher population of the same
size, by terms of order ``1/N``: the binomial sex ratio lowers the
variance effective size, and a rare allele competes only within the sex
of its carrier. See [Validation of selection](@ref).

If no member of one sex has positive fitness (lethal genotypes), the
population cannot reproduce and the run stops with an error. For the
same seed, a run with selection shares no realization with a neutral
run.

## Selection on a trait

`fitness = TraitFitness(loci; σe, fitness)` selects on a quantitative
trait instead of on the sites directly. `loci` is a `TraitLoci` table:
each site has a position, an effect ``a`` of its derived allele and a
dominance ``h``, and the table has a mutation rate ``\mu`` per site per
haplotype per generation. `trait_loci` draws the sites, and the effects
from a distribution if `a` is one.

**Model.** The genotypic value ``G`` of an individual is the sum over
sites of ``0``, ``ha`` and ``a`` for zero, one and two copies of the
derived allele (``h = ½`` is additive), including sites whose derived
allele was fixed and moved to `selected_substitutions`. The phenotype
is ``P = G + e``, with ``e \sim N(0, \sigma_e^2)`` drawn anew for every
adult in every generation, serially from `rng`, so that the run does
not depend on the number of threads. Fitness is a function of ``P``:

| Fitness | ``w(P)`` |
|---|---|
| `Directional(β)` | ``e^{\beta P}``; with a normal phenotype, the selection differential is ``\beta\sigma_P^2`` |
| `Stabilizing(θ, ω²)` | ``e^{-(P - \theta)^2/(2\omega^2)}`` |
| any function | its value, which must be non-negative; for example `P -> P > 0 ? 1.0 : 0.0` for truncation |

Selection is soft and within each sex, as with `SiteFitness`: every
offspring draws its sire among the males and its dam among the females
with probability proportional to ``w(P)``, with replacement.

```@example simulation
chr = [1_000_000, 1_000_000]
qtl = trait_loci(chr, 200; a=Normal(0, 0.2), μ=2e-4, rng=Xoshiro(1))
up = fisher_wright(100, 200, chr, 1.0; result=true, rng=Xoshiro(2),
                   fitness=TraitFitness(qtl; σe=1.0,
                                        fitness=Directional(0.5)))
g = [1, 50, 100, 150, 200]
(up.mean_genotypic[g], up.var_genotypic[g])
```

The result records, for the parents of every generation, the mean and
variance of ``G`` and of ``P`` (`mean_genotypic`, `var_genotypic`,
`mean_phenotypic`, `var_phenotypic`; `NaN` for generations without a
trait model), besides the fitness record. Continuing the population
under stabilizing selection around ``\theta = 2`` pulls the mean back:

```@example simulation
held = fisher_wright(100, 200, chr, 1.0; init=up, result=true,
                     rng=Xoshiro(3),
                     fitness=TraitFitness(qtl; σe=1.0,
                                          fitness=Stabilizing(2.0, 4.0)))
held.mean_genotypic[[200, 250, 300, 350, 400]]
```

`to_haplotype(result, qtl)` exports the genotypes at the trait loci
with their effects (columns `chr`, `pos`, `a`, `h`, `frq`), from which
the true genotypic values of the population follow. A fixed site adds
its ``a`` to every ``G``; it is kept as an offset when it is moved out
of the haplotypes, so the realization does not depend on
`fixation_interval`.

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
