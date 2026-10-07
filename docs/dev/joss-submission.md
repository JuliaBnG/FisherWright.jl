# JOSS submission checklist

Status as of 2026-10-07 (v0.3.8).

## 1. Before submitting

- [ ] Commit and push the paper work and pending source changes
      (`README.md`, `docs/src/manual/simulation.md`, `paper/paper.md`,
      `paper/paper.bib`, `src/fwp.jl`, `test/runtests.jl`); add
      `.github/workflows/draft-pdf.yml` and the JOSS checklist and roadmap
      under `docs/dev/`. CI and the docs build must be green on the default
      branch.
- [ ] Build the draft PDF with the *Draft PDF* workflow (runs on changes
      under `paper/`, or manually from the Actions tab) and download the
      `paper` artifact. Check:
  - [ ] citations, math (Sved's $E[r^2]$, the line-broken
        `$L = 100$` / `Mb`) and bibliography render correctly;
  - [ ] ORCID and ROR links resolve.
- [x] Every bibliography entry has a DOI (`yu2026xybng`:
      `10.5281/zenodo.22515787`).
- [x] Update the front-matter `date` if submitting after 6 October 2026.
- [ ] Repository requirements:
  - [x] OSI license (MIT), `CONTRIBUTING.md`, documentation, automated
        tests, public history since 2025-09-12;
  - [x] README/docs have installation instructions, a usage example,
        and community guidelines (how to report issues and get support);
        link `CONTRIBUTING.md` from the README;
  - [ ] optional: `CITATION.cff`, code of conduct.
- [ ] Tag a release (e.g. v0.3.9) that contains the paper.

## 2. Submit

- [ ] Log in at <https://joss.theoj.org/papers/new> with ORCID.
- [ ] Fill in the form: repository
      `https://github.com/JuliaBnG/FisherWright.jl`, branch (only if not
      the default), version, topic area (e.g. Biological Sciences /
      population genetics), short description, AI-use and conflict of
      interest declarations consistent with the paper.
- [ ] In the pre-review issue (`openjournals/joss-reviews`), run:
  - [ ] `@editorialbot generate pdf`
  - [ ] `@editorialbot check references`
  - [ ] `@editorialbot check repository`

## 3. Pre-review and review

- [ ] Answer the editor's scope questions. Expect overlap questions about
      msprime, SLiM and XSim; the "build vs. contribute" paragraph in
      *State of the field* is the answer.
- [ ] Suggest reviewers if asked: Julia and forward population-genetics
      simulation background, no conflicts of interest (not NMBU or EAAP
      co-authors).
- [ ] In the `REVIEW` issue (usually 2 reviewers), expect them to run the
      tests and documentation examples and to check performance and
      validation claims against `docs/src/validation/`. Reported numbers
      must match the logs exactly.
- [ ] Address reviewer issues/PRs in the repository; reply in the review
      thread with links to the fixing commits. Re-run
      `@editorialbot generate pdf` after paper changes.

## 4. Acceptance

- [ ] Tag a final release that includes all review changes.
- [ ] Archive that release on Zenodo (GitHub–Zenodo integration or manual
      upload). Title and author list must match the paper, license must
      match the repository.
- [ ] Post the version and the archive DOI in the review thread; the
      editor runs `@editorialbot set <DOI> as archive` and
      `@editorialbot set <version> as version`, then recommends
      acceptance.
- [ ] After publication: add the JOSS badge and DOI to the README, add
      `CITATION.cff` or a "How to cite" section, and cite the JOSS paper
      in revised `xyBnG.jl` and methane papers where relevant.
