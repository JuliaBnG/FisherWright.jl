# Results and export

## Default result

The default simulation result is a tuple:

```julia
haplotypes, chromosome_ends = fisher_wright(100, 100, [1_000_000], 1.0)
```

Every element of `haplotypes` is a sorted, unique `Vector{UInt32}`. A
mutation is present when its genome-wide position occurs in the haplotype.

Use `muts2bitarray` to obtain a dense representation:

```@example results
using FisherWright

mutations, chromosome_ends = fisher_wright(20, 10, [100_000], 1.0)
alleles, loci = muts2bitarray(mutations, chromosome_ends)
(size(alleles), propertynames(loci))
```

`alleles` is a `BitMatrix` with loci in rows and haplotypes in columns.
`loci` has these columns:

| Column | Meaning |
|---|---|
| `chr` | One-based chromosome index. |
| `pos` | Genome-wide base-pair position. |
| `ref` | Randomly generated reference base. |
| `alt` | Randomly generated alternate base, distinct from `ref`. |
| `frq` | Alternate-allele frequency across haplotypes. |

Fixed loci are removed by default. Pass `flip=true` to randomly invert the
reference/alternate coding at each locus.

## Structured result and substitutions

Set `result=true` to return a `FisherWrightResult`. At each
`fixation_interval`, mutations carried by every active haplotype are removed
from the active population and recorded as substitutions:

```@example results
result = FisherWrightResult(
    [UInt32[1, 3], UInt32[2], UInt32[], UInt32[1, 2, 3]],
    UInt32[3],
    UInt32[],
)
(length(result.active_haplotypes), result.chromosome_ends, result.substitutions)
```

This avoids retaining fixed mutations in every haplotype during long
simulations. The result fields are:

| Field | Meaning |
|---|---|
| `active_haplotypes` | Polymorphic mutation positions after fixed loci are removed. |
| `chromosome_ends` | Cumulative chromosome endpoints. |
| `substitutions` | Sorted fixed mutation positions. |
| `census` | Census sizes ``N_0, N_1, …, N_{nt}`` of the run. |
| `males` | Number of males among the parents of each generation. |
| `selected_haplotypes` | Derived alleles each haplotype carries at the fitness loci; empty in a neutral run. |
| `selected_substitutions` | Fitness loci whose derived allele was fixed and removed. |
| `mean_fitness`, `mean_log_fitness`, `var_log_fitness` | Fitness of the parents of each generation: mean, and mean and variance of its logarithm. |
| `mean_genotypic`, `var_genotypic`, `mean_phenotypic`, `var_phenotypic` | Under `TraitFitness`, the mean and variance of the genotypic value and of the phenotype of the parents of each generation; `NaN` otherwise. |

Convert a structured result to `BnGStructs.Haplotype` and a linkage map:

```@example results
haplotype, loci = to_haplotype(result)
(haplotype.nlc, haplotype.nhp, length(loci.pos))
```

Use `to_haplotype(result; include_fixed=true)` when fixed substitutions should
also appear in the dense output. If no polymorphic loci are available,
`to_haplotype` raises an error rather than returning an empty `Haplotype`.

## Census record and effective size

A result returned by `fisher_wright` records the census size of every
generation and the realized number of males among its parents.
`census[1]` is the starting population and `census[g + 1]` the number of
individuals born in generation `g`; `males[g]` of the `census[g]`
parents of generation `g` were males. A run continued from a
`FisherWrightResult` through `init` appends its generations to the
record of that result, so a burn-in and its continuation are numbered as
one run; started from plain haplotypes, the record begins anew. A result
built by hand, as above, has `census` holding the number of individuals
only and `males` empty.

```@example results
sizes = [fill(10, 5); fill(40, 15)]
run = fisher_wright(40, 20, [100_000], 1.0; popsize=sizes, result=true)
(run.census[1:8], run.males[1:8])
```

The record also holds the realized family sizes of every generation:
`paternal_pairs[g]` and `maternal_pairs[g]` are the numbers of ordered
pairs of offspring of generation `g` that share their sire,
``\sum_s m_s(m_s - 1)``, and that share their dam,
``\sum_d n_d(n_d - 1)``. Two genes in different offspring therefore
derive from the same parent with probability

```math
a_g = \frac{\sum_s m_s(m_s - 1) + \sum_d n_d(n_d - 1)}{4N_g(N_g - 1)},
```

``N_g`` the number of offspring, and `effective_size` returns
``1/a_g`` for a generation and the harmonic mean ``1/\bar{a}`` for
several:

```@example results
(effective_size(run), effective_size(run, 2:6), effective_size(run, 7:20))
```

This holds for every mating scheme and includes the variance of family
size that selection adds within a generation. It is the value realized
by the pedigree of the run, not an estimate from the genetic data, and
it varies from run to run as family sizes do. It depends on the
parents: generation 6 above has 40 offspring but the effective size of
its 10 parents.

`effective_size(run; from=:sexes)` uses the numbers of male and female
parents only, ``a_g = 1/(4N_m) + 1/(4N_f)``, which is the expectation of
the realized value under random mating and gives Wright's (1931)
``4N_mN_f/(N_m + N_f)``; it does not apply to other mating schemes or
to selection. Until version 0.6 this was the only form.

```@example results
(effective_size(run, 7:20), effective_size(run, 7:20; from=:sexes))
```

At constant numbers the asymptotic rate of loss of heterozygosity,
``1/(2N_e)``, corresponds to this value plus ½ (see
[Effective size](@ref)). Selection that is heritable lowers the
effective size further over generations (Robertson 1961), which a
record of one generation at a time does not include.

## Selected alleles

After a run with selection (see [Selection](@ref)) the derived alleles
at the fitness loci are in `selected_haplotypes` and
`selected_substitutions`, apart from the neutral mutations.
`to_haplotype(result, loci)` returns them as genotypes at the sites of
the table, with their true effects:

```@example results
using Random

chr = [100_000, 100_000]
loci = fitness_loci(chr, 200; s=0.02, μ=5e-4, rng=Xoshiro(1))
sel = fisher_wright(50, 300, chr, 1.0; mut_base=1e6,
                    fitness=SiteFitness(loci), result=true,
                    rng=Xoshiro(2))
qtl, effects = to_haplotype(sel, loci)
first(effects, 5)
```

```@example results
(segregating=count(x -> 0 < x < 1, effects.frq),
 fixed=count(==(1), effects.frq))
```

Row `k` of `qtl` is site `loci.pos[k]`, for every site of the table
whether or not it segregates, so the rows line up with the table and
with `effects`. A site whose derived allele was fixed and moved to
`selected_substitutions` is exported as carried by every haplotype
(`frq = 1`); pass `include_fixed=false` to export it as absent. The
dosage of individual `i` at site `k` is `qtl[k, 2i-1] + qtl[k, 2i]`.

The result also records, for the parents of every generation, the mean
fitness and the mean and variance of log fitness. Entry `g` describes
the same `census[g]` adults as `males[g]`. Fitness is absolute, relative
to an individual without any derived allele at the fitness loci, so
fixed sites count; a neutral run records `1`, `0` and `0`.

```@example results
g = [1, 100, 200, 300]
(sel.mean_fitness[g], sel.var_log_fitness[g])
```

Under soft selection the mean fitness itself has no effect on the
population; it is the quantity that mutation–selection balance and
the genetic load are stated in. A lethal genotype among the parents
makes the mean of log fitness `-Inf` and its variance `NaN`, while the
mean fitness stays finite.

To export the causal variants together with the neutral mutations, as
sequence data would contain them, add `selected=true` to the other
methods of `to_haplotype`:

```@example results
sequence, map = to_haplotype(sel; selected=true)
(size(sequence), count(in(loci.pos), map.pos))
```

A neutral mutation can fall on the coordinate of a fitness locus. In
the combined export such a coordinate is one locus, derived on a
haplotype that carries either; `to_haplotype(result, loci)` reports the
selected allele alone. `muts2bitarray(result.selected_haplotypes,
result.chromosome_ends)` gives the segregating selected sites as a
`BitMatrix`.

## Chip extraction

Extract a selected set of marker coordinates directly from a structured result
without first materializing every simulated locus:

```@example results
chip_positions = UInt32[10_000, 50_000, 90_000]
chip_haplotype = to_haplotype(result, chip_positions)
size(chip_haplotype)
```

`chip_positions` must be sorted and unique. Each requested position is
retained even when it is absent from the population, in which case its row is
all false. `extract_chip_bitarray(result.active_haplotypes, chip_positions)`
returns the corresponding `BitMatrix` when a `BnGStructs.Haplotype` wrapper is
not needed. A `BnGStructs.LocusSet` can also select indices from a shared,
sorted coordinate vector:

```julia
chip_haplotype = to_haplotype(result, locus_set, all_positions)
```
