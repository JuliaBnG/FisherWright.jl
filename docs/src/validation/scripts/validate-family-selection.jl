# Drift under FamilySize with selection, against the pedigree.
#
#   julia -t 32 --project=. docs/src/validation/scripts/validate-family-selection.jl RUNS SEED
#
# At a neutral locus unlinked to every selected locus, the expected
# heterozygosity of generation t given the pedigree is
#
#     E[H_t | pedigree] = H* (1 − θ̄_t),
#
# H* the heterozygosity of the base from which the founders' alleles are
# drawn and θ̄_t the mean coancestry of generation t over all ordered pairs
# of individuals, self-coancestry (1 + F)/2 included. Mendelian segregation
# at an unlinked locus is independent of the genotypes selection acts on,
# so this holds exactly for any selection, including the cumulative effect
# of heritable fitness on drift (Robertson 1961), which acts through the
# pedigree. The check is that FisherWright's transmission of alleles and its
# pedigree record agree, generation after generation, under FamilySize with
# and without selection.
#
# Model: N = 100 diploids, FamilySize matings (offspring split among the
# dams, one sire per dam, dams split among the sires; binomial sex),
# 100 generations. Selection through TraitFitness: 400 additive loci
# (h = ½, effect ±0.5, sign per locus and run) on chromosomes 1–4, with
# mutation μ = 1e-3 per locus and generation from a monomorphic start,
# σe = 1. Neutral markers: 1,000 loci on chromosomes 5–24 (50 each), with
# founder alleles drawn independently at frequency ½ (H* = ½); no
# mutation at them (mut_base = 1e12). Every chromosome is 100 Mb and
# 1 Morgan, so the markers are unlinked to the trait loci.
#
# Scenarios:
#   1  FamilySize(dams = :random), no selection (control)
#   2  FamilySize(dams = :random), Directional(0.5)
#   3  FamilySize(dams = 1.0), Directional(0.5)   (Gamma(1) dam weights)
#   4  FamilySize(dams = 1.0), Stabilizing(0, 4)
#
# Per run and every 10 generations: H_t, the mean over markers of 2pq in
# generation t, and θ̄_t from the recorded pedigree. Test: the mean over
# runs of H_t − ½(1 − θ̄_t) is 0 (|z| ≤ 3.5, SE over runs). Reported for
# comparison: the drift effective size from θ̄ over generations 50–100,
# and `effective_size` (the within-generation value from the family
# sizes, which leaves out the cumulative effect).
using FisherWright, Random, Statistics, Printf

const N = 100
const T = 100
const STEP = 10
const NQ = 4                      # chromosomes with trait loci
const NM = 20                     # chromosomes with markers
const CHRLEN = 100_000_000
const KQ = 400
const KM = 50                     # markers per chromosome
const ALPHA = 0.5
const MUQ = 1e-3
const SCEN = [
    ("random, no selection", () -> FamilySize(dams = :random), nothing),
    ("random, Directional(0.5)", () -> FamilySize(dams = :random), Directional(0.5)),
    ("Gamma(1), Directional(0.5)", () -> FamilySize(dams = 1.0), Directional(0.5)),
    ("Gamma(1), Stabilizing(0, 4)", () -> FamilySize(dams = 1.0), Stabilizing(0.0, 4.0)),
]

chr = fill(CHRLEN, NQ + NM)
qpos = collect(UInt32, range(CHRLEN ÷ (2KQ ÷ NQ), step = NQ * CHRLEN ÷ KQ, length = KQ))
mpos = UInt32[(NQ + c - 1) * CHRLEN + CHRLEN ÷ (2KM) + (k - 1) * (CHRLEN ÷ KM) for c = 1:NM for k = 1:KM]

# Mean coancestry of the new generation from the previous one's matrix.
# `prev` maps an ID to its row in Θ.
function next_coancestry(Θ, prev, sire, dam)
    n = length(sire)
    s = [prev[x] for x in sire]
    d = [prev[x] for x in dam]
    Θn = Matrix{Float64}(undef, n, n)
    for j = 1:n, i = 1:j
        Θn[i, j] = Θn[j, i] = i == j ? (1 + Θ[s[i], d[i]]) / 2 :
                                  (Θ[s[i], s[j]] + Θ[s[i], d[j]] + Θ[d[i], s[j]] + Θ[d[i], d[j]]) / 4
    end
    return Θn
end

function heterozygosity(haps)
    c = zeros(Int, length(mpos))
    for h in haps, x in h
        k = searchsortedfirst(mpos, x)
        k <= length(mpos) && mpos[k] == x && (c[k] += 1)
    end
    n = length(haps)
    mean(2 * (x / n) * (1 - x / n) for x in c)
end

function one_run(scen, seed)
    name, mating, fit = SCEN[scen]
    rng = Xoshiro(seed)
    init = [UInt32[] for _ = 1:2N]
    for h = 1:2N, x in mpos
        rand(rng, Bool) && push!(init[h], x)
    end
    kw = (M = 1e8, mut_base = 1e12, result = true, mating = mating(), pedigree = STEP, rng = rng)
    if fit !== nothing
        α = ALPHA .* rand(rng, (-1.0, 1.0), KQ)
        loci = TraitLoci(qpos, 2 .* α; h = 0.5, μ = MUQ)
        kw = (; kw..., fitness = TraitFitness(loci; σe = 1.0, fitness = fit))
    end
    Θ = Matrix{Float64}(I_half(N))           # founders: unrelated, non-inbred
    ids = collect(1:N)
    out = zeros(T ÷ STEP + 1, 2)
    out[1, :] = [heterozygosity(init), mean(Θ)]
    res = nothing
    for s = 1:T÷STEP
        res = fisher_wright(N, STEP, chr, 0.02; kw..., init = s == 1 ? init : res)
        ped = res.pedigree
        for g in sort(unique(ped.generation[ped.sire .> 0]))
            rows = ped[ped.generation .== g, :]
            prev = Dict(x => i for (i, x) in enumerate(ids))
            Θ = next_coancestry(Θ, prev, rows.sire, rows.dam)
            ids = rows.id
        end
        ids == collect(res.first_id:res.first_id+N-1) || error("pedigree and haplotypes are out of step")
        out[s+1, :] = [heterozygosity(res.active_haplotypes), mean(Θ)]
    end
    return out, effective_size(res)
end

I_half(n) = [i == j ? 0.5 : 0.0 for i = 1:n, j = 1:n]

function main(runs, seed)
    t0 = time()
    println("# FisherWright ", pkgversion(FisherWright), ", Julia ", VERSION, ", threads ", Threads.nthreads(),
        ", ", gethostname(), "; $runs runs per scenario, seed $seed")
    println("N = $N, $T generations; trait: $KQ loci on $NQ chromosomes, |α| = $ALPHA, μ = $MUQ, σe = 1;")
    println("markers: $(NM * KM) on $NM other chromosomes (unlinked), founder alleles at ½, H* = ½")
    ok = true
    gens = 0:STEP:T
    for scen in eachindex(SCEN)
        tabs = Vector{Matrix{Float64}}(undef, runs)
        nes = zeros(runs)
        Threads.@threads for r = 1:runs
            tabs[r], nes[r] = one_run(scen, hash((seed, scen, r)))
        end
        println("\n", scen, "  ", SCEN[scen][1])
        @printf("%4s %18s %18s %18s %8s\n", "gen", "H (markers)", "½(1 − θ̄)", "difference", "z")
        for (k, g) in enumerate(gens)
            h = [x[k, 1] for x in tabs]
            e = [0.5 * (1 - x[k, 2]) for x in tabs]
            d = h .- e
            se = std(d) / sqrt(runs)
            z = se > 0 ? mean(d) / se : 0.0
            abs(z) <= 3.5 || (ok = false)
            @printf("%4d %9.5f ± %6.5f %9.5f ± %6.5f %9.5f ± %6.5f %8.2f\n", g, mean(h), std(h) / sqrt(runs),
                mean(e), std(e) / sqrt(runs), mean(d), se, z)
        end
        # drift effective size from θ̄ between generations 50 and 100
        i, j = findfirst(==(50), gens), findfirst(==(T), gens)
        ne_θ = [1 / (2 * (1 - ((1 - x[j, 2]) / (1 - x[i, 2]))^(1 / (T - 50)))) for x in tabs]
        @printf("drift Ne from θ̄ (generations 50–%d): %.2f ± %.2f; effective_size (family sizes): %.2f ± %.2f\n",
            T, mean(ne_θ), std(ne_θ) / sqrt(runs), mean(nes), std(nes) / sqrt(runs))
    end
    @printf("\n%s (%.0f s)\n", ok ? "all comparisons within |z| ≤ 3.5" : "FAILED: a comparison exceeds |z| = 3.5", time() - t0)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    runs = parse(Int, get(ARGS, 1, "16"))
    seed = parse(Int, get(ARGS, 2, "2026"))
    main(runs, seed) || exit(1)
end
