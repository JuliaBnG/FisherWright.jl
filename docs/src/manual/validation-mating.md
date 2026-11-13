# Validation of the mating schemes

`docs/src/validation/scripts/validate-mating.sh` compares the mating
schemes of FisherWright with a matched model in SLiM 5
(`validate-mating.slim`), in which each scheme is written anew in
Eidos and runs on SLiM's own meiosis and pedigree. Where the expected
heterozygosity is known exactly, it is compared with that as well.

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
