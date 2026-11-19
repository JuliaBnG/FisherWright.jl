# Validation of the mating schemes

`docs/src/validation/scripts/validate-mating.sh` compares the mating
schemes of FisherWright with a matched model in SLiM 5
(`validate-mating.slim`), in which each scheme is written anew in
Eidos and runs on SLiM's own meiosis and pedigree. Where the expected
heterozygosity is known exactly, it is compared with that as well.
[Mutation–drift equilibrium](@ref) checks the diversity maintained by
the schemes under mutation against exact references.

**Model.** Both simulators run ``N = 40`` diploid individuals for 30
discrete generations without mutation, from a polymorphic start: each
of ``L = 1{,}000`` sites, on 20 chromosomes of 50 bp and 1 Morgan, is
carried by each haplotype with probability ½. In SLiM the model is
nonWF with explicit reproduction: every generation, the offspring of
each population are created from a table of sire and dam drawn by the
scheme, and the parents then die. SLiM has no sexes here: as in
FisherWright, the roles of sire and dam are assigned anew to the
parents of every generation (Bernoulli(½) per individual, redrawn
until both sexes, or the breeders a design needs, are present; or a
fixed number of males, a random subset). Relatives are judged on SLiM's
pedigree IDs (`pedigreeParentIDs`, `pedigreeGrandparentIDs`), with the
starting population of unknown parents, as in FisherWright. Crossovers
differ slightly: 49 gaps of 1/49 per chromosome in SLiM, a
Poisson(1) number at uniform positions in FisherWright; this changes
the linkage between sites, but not the expectations compared.

!!! note "Recombination map of the reported runs"
    The runs reported below used an earlier version of the map, with
    the gap of 0.5 one site too early: in SLiM a breakpoint at
    position ``x`` falls between ``x - 1`` and ``x``, and the interval
    of rate 0.5 ended at the last site of each chromosome rather than
    at the first site of the next. The last site of each chromosome
    was thus unlinked from it and linked at 1/49 to the next one. The
    statistics compared do not depend on linkage, so the results
    stand; the script has been corrected since (the 0.5 interval now
    ends at ``50(c + 1)``, the first site of chromosome ``c + 1``).

**Schemes.**

| Key | FisherWright | SLiM implementation |
|---|---|---|
| `rm` | `RandomMating()` | sire and dam uniform among the males and females |
| `rm5` | `RandomMating(nmale=5)` | the same with 5 males |
| `fsgamma` | `FamilySize(dams=1.0, sires=:random, nmale=10)` | dams of the offspring drawn with Gamma(1) weights; each dam's sire uniform among the males |
| `hier` | `Hierarchical(4, 3)` | 4 sires and 12 dams drawn from the parents; 3 dams per sire, equal dam families |
| `mono` | `Monoecious()` | both parents uniform among all; one parent twice is a selfed offspring (`addSelfed`) |
| `self` | `Monoecious(selfing=0.5)` | selfed with probability ½, otherwise one of the others |
| `avoid_half` | `AvoidRelatives(RandomMating(); level=:halfsib)` | as `rm`, the dam redrawn while she shares a parent with the sire, at most 100 draws |
| `avoid_cousin5` | `AvoidRelatives(RandomMating(nmale=5); level=:cousin)` | as `rm5`, the dam redrawn while she shares a parent or grandparent with the sire |
| `avoid_mono` | `AvoidRelatives(Monoecious(selfing=0); level=:fullsib)` | outcrossing only, the mate redrawn while a full sib |

**Statistics.** At the end, Σ2pq per site over the ``2N`` haplotypes
(with replacement), the observed heterozygosity ``H_{obs}``, and the
number of offspring of relatives that could not be avoided within 100
draws (`related_matings`, summed over the run). Each simulator runs
4,000 replicates per scheme (SLiM: 40 processes of 100 populations);
the table gives mean ± SE, the z-score of the difference between the
simulators, and for the first six schemes the exact expectation, from
the recursion of the heterozygosity within (``H_W``) and between
(``H_B``) individuals,

```math
H_B' = (1 - a)H_B + \tfrac{a}{2}H_W, \qquad
H_W' = (1 - \sigma)H_B + \tfrac{\sigma}{2}H_W,
```

with ``a`` the probability that two genes of different offspring come
from one parent (from the closed forms under [Family sizes](@ref)) and
``\sigma`` the selfing rate; ``\Sigma 2pq = H_W/(2N) + (1 - 1/N)H_B``
and ``H_{obs} = H_W``. Avoidance of relatives has no such reference,
which is why it is compared with an independent implementation.

**Resolution.** The standard errors are 0.04% to 0.17% of the means. A
difference between the simulators of 3.5 combined SE corresponds to
an error of about 0.5% (random mating) to 1% (the hierarchical
design) in the effective size over the 30 generations. The schemes
themselves are far apart: avoiding half-sib matings raises ``H_{obs}``
from 0.3455 to 0.3537, about 40 SE.

The run takes about 10 minutes on 16 cores:

```bash
docs/src/validation/scripts/validate-mating.sh path/to/slim outdir 4000 2026 16
```

To run some schemes only, list them in `SCHEMES`, as for the check of
cousin avoidance below:

```bash
SCHEMES="avoid_cousin5" docs/src/validation/scripts/validate-mating.sh path/to/slim outdir 16000 11 24
```

## Results, 4,000 replicates, seed 2026

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-mating-seed2026.txt"), String)))
```

## Results, 4,000 replicates, seed 7

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-mating-seed7.txt"), String)))
```

## Cousin avoidance, 16,000 replicates, seed 11

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-mating-cousin-seed11.txt"), String)))
```

## Reading the results

Across the 2 × 9 schemes and 2 × 19 non-degenerate comparisons
between the simulators, the largest ``|z|`` is 1.77 (`rm5`, seed 7,
both statistics). FisherWright is within ``|z| \le 1.66`` of the exact
expectation for the six schemes that have one. Neither simulator is
systematically above the other or above the expectation. Where no
offspring of relatives can be forced, both report none: with about
20 females, an unrelated dam is always found within 100 draws. Under cousin
avoidance with five sires, nearly every pair is related after a few
generations, and about 237 of the 1,200 offspring per run are forced
in both simulators.

One pattern prompted a further run. In both seeds, cousin avoidance
left FisherWright slightly more heterozygous than SLiM (z between 1.14
and 1.55, combined over the seeds about 1.9) with slightly fewer forced
matings. Several schemes and statistics were compared, so this is
unremarkable as a single result, but its consistency across seeds
called for a check. Rerun alone with four times the replicates (seed
11, above), the differences are −0.61, −0.38 and +0.17 SE, so the
pattern was noise.

Within the resolution given above, FisherWright's implementations of
the mating schemes agree with SLiM's pedigree and meiosis and, where
the expectation is known, with theory. Linkage was not compared: the
two crossover models differ (49 discrete gaps against continuous
positions).

## Scripts

```@eval
using Markdown, FisherWright
dir = joinpath(pkgdir(FisherWright), "docs", "src", "validation", "scripts")
Markdown.MD(Markdown.Code("text", read(joinpath(dir, "validate-mating.slim"), String)),
    Markdown.Code("julia", read(joinpath(dir, "validate-mating.jl"), String)),
    Markdown.Code("bash", read(joinpath(dir, "validate-mating.sh"), String)))
```

## Mutation–drift equilibrium

`docs/src/validation/scripts/validate-equilibrium.jl` checks the
diversity that the mating schemes maintain under mutation and drift.
The comparison above starts from a polymorphic population and has no
mutation, so it tests drift alone.

**Model.** ``N = 100`` diploid parents run for ``T = 2{,}000 = 20N``
generations from a monomorphic start, with 10 chromosomes of 10 Mb
(0.1 Morgan each) and neutral mutation at ``\mu = 10^{-8}`` per bp per
haplotype per generation. With ``N_e`` as low as 37, the run lasts at
least 27 ``N_e`` generations; the reference values below are for
generation ``T``, not the limit, so equilibrium is not assumed. The
schemes are:

- `RandomMating()` (binomial sex, the default);
- `RandomMating(nmale = 10)` (10♂/90♀, ``N_e \approx 36``);
- `RandomMating(sexratio = 0.5)`;
- `FamilySize(dams = 1.0, sires = :random, nmale = 20)`, overdispersed
  dam families with Gamma(1) weights;
- `FamilySize(dams = :equal, sires = :equal, nmale = 20)`;
- `Monoecious()`, `Monoecious(selfing = 0.2)` and
  `Monoecious(selfing = 0.5)`.

**References.** For infinite sites, the expected heterozygosity per bp
within an individual (``h_W``) and between two genes of different
individuals (``h_B``) follows exactly from

```math
h_B' = (1 - a)h_B + \tfrac{a}{2}h_W + 2\mu, \qquad
h_W' = (1 - \sigma)h_B + \tfrac{\sigma}{2}h_W + 2\mu,
```

with ``a`` and ``\sigma`` as above. For random mating with binomial
sex, ``a`` is averaged over the number of males, conditioned on both
sexes being present. The expectation of Σ2pq per bp over all ``2N``
haplotypes is ``h_W/(2N) + (1 - 1/N)h_B``, and that of ``H_{obs}`` is
``h_W``. Both are exact.

The statistic S is the number of segregating sites in a sample of one
random gene from each of ``n = 20`` random individuals. Its
expectation is
``\mu L\,\mathrm{E}[\sum_t k_t \mathbf{1}\{k_t \ge 2\}]``,
where ``k_t`` is the number of distinct genes in generation ``t`` that
are ancestral to the sample. Each such gene carries one transmission
from generation ``t - 1``; a mutation on the only lineage is fixed in
the sample. The expectation of the branch length comes from a backward
simulation of the sample's genealogy through a pedigree, drawn one
generation at a time from the documented rules of each scheme. That
code is written in the script, independently of the package's. The
backward simulation is checked against the exact ``h_B`` with
``n = 2`` (row `h_B/bp*`), with 10⁶ genealogies for each sample size
and scheme.

Kingman's approximation, ``h_B L \sum_{i<n} 1/i``, is printed in each
scheme's header line, but it is not a reference here. With 20 lineages
and ``a`` up to 0.028, several coalescences are expected in one
generation. In a first run with that approximation as the reference,
FisherWright's S was 0.2% (``N_e \approx 100``) to 1.4%
(``N_e \approx 37``) above it, at z up to 11, and the backward
simulation reproduced the excess.

The sites of FisherWright are biallelic. At ``\mu T = 2 \cdot 10^{-5}``,
the difference from infinite sites is far below the standard errors.
Each comparison is a z-score against the standard error of 4,000
forward replicates; for S, the error of the backward simulation is
added. The script exits with status 1 if any ``|z| > 3.5``. A seed
took about 3.3 hours on 32 threads, two seeds side by side on a 64-core
server:

```bash
julia -t 32 --project=. docs/src/validation/scripts/validate-equilibrium.jl 4000 2026
```

### Equilibrium results, seed 2026

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-equilibrium-seed2026.txt"), String)))
```

### Equilibrium results, seed 7

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-equilibrium-seed7.txt"), String)))
```

### Reading the equilibrium results

Both seeds pass. Of the 2 × 8 × 4 comparisons, the largest ``|z|`` is
2.26, and it is for the backward simulation itself (`h_B/bp*`, 10♂/90♀,
seed 7); FisherWright's own largest is 1.92 (`Hobs/bp`, 10♂/90♀, seed
2026, and `S_20`, selfing 0.2, seed 2026). Under the null hypothesis,
about three of the 64 would exceed 2; one does. Σ2pq and ``H_{obs}``
of one scheme are strongly correlated, so their z-scores are not
independent.

The comparisons resolve small differences in the effective size. The
standard errors of Σ2pq are about 0.1% of the mean (0.15% at
``N_e \approx 37``), so a bias of about 0.4% in ``N_e`` would be found
at ``|z| = 3.5``. The schemes differ by far more: ``N_e`` is 36.5,
41.3, 75.5, 87.6, 90.5, 99.5, 100 and 100.5. Between binomial sex and
exactly half males, the difference of one individual in ``N_e`` shows
in Σ2pq (3.9607 against 4.0009 × 10⁻⁶, seed 2026).

``H_{obs}`` under selfing tests the within-individual recursion
separately from Σ2pq. At selfing 0.5, ``H_{obs}/\Sigma 2pq`` is 0.68,
close to ``1 - F = 2/3`` with ``F = \sigma/(2 - \sigma)`` at
equilibrium, and both statistics match their exact expectations.

S is the one statistic here that depends on more than pairs of genes.
Its agreement with the backward simulation, within 0.2% (at most 1.9
combined SE), shows that FisherWright's genealogies of 20 genes have
the right total length under every scheme, including the multiple and
simultaneous coalescences that Kingman's approximation leaves out. The
approximation falls short of the backward reference by 0.1% to 1.5%,
rising with ``a``: the expected number of pairs
coalescing in one generation, ``\binom{20}{2} a``, is 1.9 at
``N_e \approx 100`` and 5.3 at ``N_e \approx 37``. Under these
schemes, Watterson's estimator of ``4N_e\mu`` from a sample of 20
genes is therefore biased upward by about as much. This is a property
of discrete generations at small ``N``, not of the simulator.

### Equilibrium script

```@eval
using Markdown, FisherWright
dir = joinpath(pkgdir(FisherWright), "docs", "src", "validation", "scripts")
Markdown.MD(Markdown.Code("julia", read(joinpath(dir, "validate-equilibrium.jl"), String)))
```

## Assortative mating

`docs/src/validation/scripts/validate-assortative.jl` checks
[`Assortative`](@ref) against Bulmer's (1980) recursion for the genetic
variance under assortative mating.

**Model.** ``N = 2{,}000`` diploids with two sexes, a neutral trait
(fitness 1) with ``K = 1{,}000`` additive loci, one per chromosome, so
that every pair of loci recombines freely; no mutation. The starting
population is drawn in Hardy–Weinberg and linkage equilibrium, with
allele frequencies from U(0.1, 0.9) and effects scaled to genic
variance ``V_a = 1``; ``V_E = 1``. Mating is assortative with
``ρ = 0.5``, 0.8 and −0.5 for ``T = 15`` generations, 200 replicates
each.

**Reference.** With free recombination and many loci of small effect,
the additive variance of the next parents is

```math
V_A' = c \, \frac{V_A (1 + m h^2)}{2} + \frac{V_a - d/K_{eff}}{2}, \qquad h^2 = V_A / V_P,
```

where ``m`` is the phenotypic correlation of mates, ``V_a`` the genic
variance (reduced by drift by ``1 - 1/2N`` per generation) and
``c = 1 - (1/N_m + 1/N_f)/2`` the population variance of the sires and
dams relative to the sample variance ``V_A`` (Bulmer 1980; Crow and
Kimura 1970). The disequilibrium ``d = V_A - V_a`` halves every
generation, and ``V_A`` approaches ``V_a/(1 - m h^2)``. The last term
is the share of ``d`` within loci: assortment also correlates the two
alleles of a parent at one locus, which lowers its heterozygosity and
the Mendelian sampling variance it passes on. The disequilibrium of
loci ``i`` and ``j`` is in proportion to ``x_i x_j``, ``x = α^2 pq``, so
that share is ``1/K_{eff}`` with ``K_{eff} = (\sum x)^2 / \sum x^2``,
about 317 here. It is 0.2 % of ``V_A`` at ``ρ = 0.8`` and is absent
from the infinitesimal model (``K \to \infty``). The script takes
``V_A``, ``V_P`` and ``m`` from the record of each run and tests the
one-step prediction (averaged over generations), the trajectory of the
recursion from the start at generations 2, 5 and 15, and the change of
the mean genotypic value, which is 0 in expectation.

### Results, 200 replicates, seed 2026

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-assortative-seed2026.txt"), String)))
```

### Results, 200 replicates, seed 7

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-assortative-seed7.txt"), String)))
```

### Reading the results

All 30 checks pass; the largest ``|z|`` is 2.55. The realized
correlation of mates is within 0.004 of ``ρ`` in every case. The
variance rises from 1 to 1.40 (``ρ = 0.5``) and 2.07–2.09 (``ρ = 0.8``,
still approaching its equilibrium of 2.21) and falls to 0.81
(``ρ = -0.5``), as the recursion gives, with standard errors of 0.1 %
to 0.4 % of ``V_A``.

The within-locus term was added after a first run without it. Then the
one-step residuals for ``ρ > 0`` were negative in both seeds (combined
z −2.1 at ``ρ = 0.5`` and −2.5 at ``ρ = 0.8``, about 0.1 % to 0.2 % of
``V_A``), and near zero for ``ρ < 0``: the sign and size this term
predicts. With it, the one-step residuals combined over the two seeds
are −1.2, −1.1 and −0.4. Whether the dams drawn have less variance than
the females, through the weights of the copula, was also checked: by
direct computation for 1,000 females, 0.02 % less, too little to
matter.

Two patterns remain, each within the limits. The changes of the mean
have ``|z|`` of 2.55, 2.33 and 2.05 under seed 2026 but below 1 under
seed 7, and the model has no mechanism for a change of the mean (sires
are drawn uniformly, the weights of the dams are symmetric in their
normal scores): noise. At ``ρ = -0.5``, ``V_A`` at generation 15 lies
below the recursion in both seeds (z −1.44 and −2.05, combined −2.5,
0.3 % to 0.5 % of ``V_A``), the largest combined deviation among the
checks. The recursion leaves out the departure of the regression of
``A`` on ``P`` from linearity, and with 30 checks a combined ``|z|`` of
2.5 is not unusual, so this is recorded here rather than explained.

### Assortative mating script

```@eval
using Markdown, FisherWright
dir = joinpath(pkgdir(FisherWright), "docs", "src", "validation", "scripts")
Markdown.MD(Markdown.Code("julia", read(joinpath(dir, "validate-assortative.jl"), String)))
```
