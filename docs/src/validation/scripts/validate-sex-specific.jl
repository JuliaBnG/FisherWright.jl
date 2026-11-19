# Validation of sex-specific fitness (SexSpecific) in FisherWright.jl
# against the deterministic two-sex recursion of single-locus selection.
#
# With FisherWright installed, run from this directory:
#
#   julia -t 8 validate-sex-specific.jl [reps = 200] [seed = 2026]
#
# One biallelic autosomal site under selection, N = 10,000 diploids, two
# sexes (Binomial(N, 1/2), conditioned on both), random mating, no mutation,
# a genome of one chromosome of 1 kb (recombination is irrelevant with one
# site). Genotypes AA, Aa, aa (A derived) have fitness 1 + s, 1 + hs and 1
# with (s, h) = (s_m, h_m) in males and (s_f, h_f) in females. Every
# individual of the starting population is a heterozygote.
#
# Reference: with soft selection within each sex, the frequency of A in the
# gametes of the selected males and females of a generation is
#     x' = [P_AA (1 + s_m) + P_Aa (1 + h_m s_m)/2] / w̄_m,
#     y' = [P_AA (1 + s_f) + P_Aa (1 + h_f s_f)/2] / w̄_f,
# where P are the genotype frequencies of the parents (the same in both
# sexes, since sex is independent of genotype), and the offspring have
# P_AA = x'y', P_Aa = x'(1 - y') + y'(1 - x') (Owen 1953; Kidwell et al.
# 1977). The expected frequency of A among the offspring is (x' + y')/2.
# This holds exactly in expectation for one generation given the parents
# (tested in CI, "Sex-specific fitness"). Over t generations the mean of a
# finite population lags the deterministic trajectory, since drift erodes
# pq by about t/(2N) and with it the change by selection: for additive
# selection of net strength s (s/2 per sex for one sex only) the lag is
# of order s pq t²/(8N), at most about 8e-4 here (cases 2 and 3 at
# generation 80), up to 0.7 standard errors of the mean at 200 replicates.
# More replicates make this bias visible; it is not an error of the
# engine.
#
# Scenarios:
#   1. sexually antagonistic, additive: s_m = 0.1, s_f = -0.1, h = 1/2;
#   2. selection in males only, s_m = 0.1, h = 1/2; and
#   3. the same in both sexes at half strength, s = 0.05: to first order in
#      s the two have the same trajectory, and the recursion gives the
#      difference at second order;
#   4. antagonistic with dominance reversal (the allele dominant where it
#      is favoured): s_m = 0.2, h_m = 0.8, s_f = -0.2, h_f = 0.2, a
#      protected polymorphism (Kidwell et al. 1977) with the equilibrium
#      given by the recursion.
#
# Each check is a z-score of the mean frequency over replicates against the
# recursion; a check fails when |z| > ZMAX, and the script exits with
# status 1.

using FisherWright, Random, Statistics, Printf

const ZMAX = 3.5
const N = 10_000
const CHR = [1_000]
const SITE = UInt32[500]

# Deterministic two-sex recursion from all heterozygotes: the expected
# offspring frequency after each of `T` generations.
function recursion(sm, hm, sf, hf, T)
    PAA, PAa = 0.0, 1.0
    out = zeros(T)
    for t = 1:T
        Paa = 1 - PAA - PAa
        wm = PAA * (1 + sm) + PAa * (1 + hm * sm) + Paa
        wf = PAA * (1 + sf) + PAa * (1 + hf * sf) + Paa
        x = (PAA * (1 + sm) + PAa * (1 + hm * sm) / 2) / wm
        y = (PAA * (1 + sf) + PAa * (1 + hf * sf) / 2) / wf
        out[t] = (x + y) / 2
        PAA, PAa = x * y, x * (1 - y) + y * (1 - x)
    end
    return out
end

model(s, h) = s == 0 ? Neutral() : SiteFitness(FitnessLoci(SITE, s; h=h))

frequency(r) = isempty(r.selected_substitutions) ? sum(length, r.selected_haplotypes) / 2N : 1.0

# Frequencies of A at the generations in `checks`, one row per replicate.
function simulate(fit, checks, reps, rng)
    start = FisherWrightResult([UInt32[] for _ = 1:2N], UInt32[CHR[1]], UInt32[];
        selected_haplotypes=[isodd(i) ? copy(SITE) : UInt32[] for i = 1:2N])
    Q = zeros(reps, length(checks))
    for r = 1:reps
        res, done = start, 0
        for (k, t) in pairs(checks)
            res = fisher_wright(N, t - done, CHR, 1.0; init=res, fitness=fit, result=true, mutate=false,
                rng=rng)
            done = t
            Q[r, k] = frequency(res)
        end
    end
    return Q
end

function main(reps=200, seed=2026)
    println("="^80)
    println(" FisherWright.jl validation of sex-specific fitness")
    println(" One site, N = $N, all heterozygous at the start, no mutation; reps = $reps, seed = $seed,")
    println(" threads = $(Threads.nthreads()). A check fails when |z| > $ZMAX")
    println("="^80)
    rng = Xoshiro(seed)
    ok = true
    cases = (
        ("1. antagonistic, additive", (0.1, 0.5, -0.1, 0.5), [10, 20, 40, 80]),
        ("2. males only", (0.1, 0.5, 0.0, 0.5), [10, 20, 40, 80]),
        ("3. both sexes at half strength", (0.05, 0.5, 0.05, 0.5), [10, 20, 40, 80]),
        ("4. antagonistic, dominance reversal", (0.2, 0.8, -0.2, 0.2), [10, 20, 40, 80, 160]),
    )
    for (label, (sm, hm, sf, hf), checks) in cases
        fit = SexSpecific(model(sm, hm), model(sf, hf))
        t = @elapsed Q = simulate(fit, checks, reps, rng)
        E = recursion(sm, hm, sf, hf, checks[end])
        @printf("\n  %s: s_m = %g, h_m = %g, s_f = %g, h_f = %g  (%.0f s)\n", label, sm, hm, sf, hf, t)
        @printf("  %-12s %21s %10s %8s\n", "generation", "mean frequency", "recursion", "z")
        for (k, g) in pairs(checks)
            m, se = mean(Q[:, k]), std(Q[:, k]) / sqrt(reps)
            z = (m - E[g]) / se
            pass = abs(z) <= ZMAX
            ok &= pass
            @printf("  %-12d %10.5f ± %-8.5f %10.5f %8.2f  %s\n", g, m, se, E[g], z, pass ? "pass" : "FAIL")
        end
        if label[1] == '4'
            @printf("  equilibrium of the recursion %.5f\n", recursion(sm, hm, sf, hf, 10_000)[end])
        end
    end
    E2, E3 = recursion(0.1, 0.5, 0.0, 0.5, 80), recursion(0.05, 0.5, 0.05, 0.5, 80)
    @printf("\n  Cases 2 and 3 by the recursion at generations 10, 40, 80: %s and %s\n",
        string(round.(E2[[10, 40, 80]]; digits=5)), string(round.(E3[[10, 40, 80]]; digits=5)))

    println("\n", "="^80)
    println(ok ? " All checks passed." : " FAILED")
    println("="^80)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(parse.(Int, ARGS)...) ? 0 : 1)
end
