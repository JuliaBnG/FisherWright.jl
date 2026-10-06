---
title: 'FisherWright.jl: forward-time Wright–Fisher simulation of whole-genome founder populations in Julia'
tags:
  - Julia
  - population genetics
  - forward simulation
  - Wright–Fisher model
  - breeding simulation
  - genomic prediction
authors:
  - name: Xijiang Yu
    orcid: 0000-0000-0000-0000  # TODO: add ORCID
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Animal and Aquacultural Sciences, Faculty of Biosciences, Norwegian University of Life Sciences (NMBU), Ås, Norway
    index: 1
    ror: 04a1mvv97  # TODO: verify NMBU ROR ID
date: 6 October 2026
bibliography: paper.bib
---

# Summary

Stochastic simulations of breeding programs and complex traits start from a founder population, whose
allele frequencies and linkage disequilibrium (LD) determine marker–QTL associations, genomic prediction
accuracy and the response to selection. `FisherWright.jl` generates such founders by forward-time simulation
of a diploid Wright–Fisher population. The population has constant census size $N$, two sexes, random union
of a sire and a dam drawn with replacement, Poisson crossovers with free assortment between chromosomes, and
neutral mutation in a biallelic finite-sites model. Starting from no variation, it runs for a user-chosen
number of generations, moves fixed mutations out of the active haplotypes, and returns the segregating
variants. The variants are returned as a `BitMatrix` with a linkage map, as a targeted SNP-chip panel, or
directly as `BnGStructs.jl` haplotype containers for downstream breeding simulation. The package is written
entirely in Julia [@bezanson2017julia], is MIT-licensed, and is registered in the Julia General Registry.

# Statement of need

Quantitative geneticists simulating genomic selection need founders whose LD and site frequency spectrum
(SFS) arise from an explicit evolutionary process. Independent Bernoulli sampling gives no LD. Heuristic
approximations, such as using Sved's $E[r^2] = 1/(1 + 4N_e c)$ as a target, do not reproduce the LD decay of
a finite population under a MAF filter. Coalescent tools produce such founders efficiently, but a Julia
breeding pipeline must then call Python or C programs and exchange files with them. `FisherWright.jl` is for
users who want founders generated inside the Julia session, with the same data structures as the rest of
their simulation. It also suits teaching and method development, where the forward process itself should be
transparent and modifiable.

# State of the field

Coalescent simulators such as `msprime` [@kelleher2016efficient; @baumdicker2022efficient] and `MaCS`
[@chen2009fast] are the fastest way to obtain neutral equilibrium genomes. `MaCS` underlies founder
generation in `AlphaSimR` [@gaynor2021alphasimr]. General forward simulators such as `SLiM`
[@haller2019tree; @haller2026slim] and `fwdpy11` [@thornton2019polygenic] support selection, demography and
spatial structure, through the Eidos scripting language or Python with a C++ core. Breeding-oriented
simulators such as `QMSim` [@sargolzaei2009qmsim], `MoBPS` [@pook2020mobps] and `XSim`
[@cheng2015xsim; @chen2022xsim] focus on pedigrees and selection programs. `XSim` version 2 is implemented in
Julia.

We built a new package rather than contributing to these for three reasons. First, the `JuliaBnG` packages
needed a founder generator with no non-Julia dependencies, whose output is the `BnGStructs.jl` haplotype type
used throughout the pipeline. Second, it had to be small enough to validate exhaustively against theory.
Third, `FisherWright.jl` was refactored out of the author's breeding simulator `xyBnG.jl` [@yu2026xybng],
which is not publicly available, so that the founder step could be released, tested and cited on its own.
`FisherWright.jl` deliberately does not model selection, demography or population structure; users who need
these should use `SLiM` or `fwdpy11`. For neutral equilibrium founders where staying in Julia doesn't
matter, `msprime` is much faster.

# Software design

Each haplotype is a sorted `Vector{UInt32}` of derived-allele positions on a concatenated genome of up to
$2^{32} - 1$ bp. We chose this over a dense genotype matrix, which is infeasible at gigabase scale, and over
tree-sequence recording, which makes forward simulation efficient but would add a dependency and a second
data model. With sorted positions, meiosis is a single two-pointer pass over the two parental haplotypes
that switches strand at each sorted crossover point (`recombine`). New mutations are added by an in-place
sorted merge (`merge_sorted!`). Both kernels write into preallocated, thread-local buffers, and the test
suite checks that they do not allocate after warm-up. Mating pairs are processed in parallel with
`Threads.@threads`; results are reproducible for a fixed seed and thread count.

Fixed mutations would otherwise accumulate at rate $\mu L$ per generation and slow every meiosis. When
`result = true`, a scan at a chosen interval removes positions shared by all $2N$ haplotypes and records them
as substitutions in a `FisherWrightResult`. The active haplotypes then hold only segregating variants, and
memory stays bounded over long runs.

Export is separated from simulation. `muts2bitarray` packs the segregating variants into a
1-bit-per-allele `BitMatrix` with a `DataFrame` linkage map. `extract_chip_bitarray` and
`to_haplotype(result, positions)` take a sorted list of chip positions and build only those rows, without
materializing the whole-genome matrix.

The model is kept minimal: constant size, neutrality, random mating and a uniform recombination rate.
Every feature added to the simulator must also be validated, and the package is meant to be one step in a
larger pipeline, not a general simulator.

# Validation

The validation suite compares `FisherWright.jl` v0.3.8 with neutral theory and with `msprime`'s
discrete-time Wright–Fisher model under the same parameters. The setting is $N = 100$, $20N$ generations,
$10 \times 10$ Mb, $\mu = r = 10^{-8}$, and 50 replicates of each simulator. The 24 prespecified checks
cover crossover counts and dispersion, independent assortment, segregating sites against Watterson's
$\theta a_n$, $\sum 2pq$, the whole-population and sample SFS, and $\sigma^2_d$ in five distance bins
from 10 kb to 3 Mb. A check fails when $|z| > 3.5$. All checks pass for two seeds (maximum $|z|$ 2.49 and
2.24).

Faster CI tests run on every commit. They are shown to pass for five seeds, and to fail when $N_e$ is
halved, when the recombination rate is doubled, and when an earlier generation-loop bug is reintroduced.

The validation scripts, raw logs, and a comparison with an early release, which documents what changed and
why, are part of the online manual.

# Research impact statement

`FisherWright.jl` is the founder generator of the `JuliaBnG` packages, alongside `BnGStructs.jl`. An early
version (v0.1.x) generated the historical population of a stochastic evaluation of methane selection in the
Norwegian White Sheep breeding program [@yu2026effects]. That version predates the corrections in the
package changelog, and the effect of those corrections on the founders is documented in the manual.

The repository includes reproducible benchmark drivers. In one 2,000-generation run of a 1 Gb genome
($10 \times 100$ Mb) with $N = 2{,}000$, the simulation took 59.8 s on 8 threads. Packing 567,975
segregating sites for 4,000 haplotypes into a 271 MiB `BitMatrix` took a further 9.0 s. For comparison,
`SLiM` with tree-sequence recording took 51.6 s plus 11.2 s for export, and `msprime` DTWF took 9.6 s plus
7.6 s.

# AI usage disclosure

<!-- TODO: author to verify and complete; JOSS requires full disclosure. -->
Generative AI tools were used during development and writing:

- GitHub Copilot, for code suggestions and some commits (marked as co-authored in the Git history);
- GPT-5, for an early code review;
- Gemini, for drafts of earlier statistical tests, which were later replaced, and for manuscript editing;
- Google Antigravity, for manuscript wording;
- Claude (Anthropic), for code and manuscript review, validation and comparison scripts, documentation,
  and a first draft of this paper.

The author reviewed, edited and validated all AI-assisted output and made the design decisions. Every
reported number comes from logged runs in the repository.

# Acknowledgements

<!-- TODO -->

# References
