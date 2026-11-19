# Validation of assortative mating (Assortative) in FisherWright.jl against
# Bulmer's (1980) recursion for the genetic variance under assortative
# mating.
#
# With FisherWright installed, run from this directory:
#
#   julia -t 8 validate-assortative.jl [reps = 200] [seed = 2026]
#
# Model. N diploids, two sexes (Binomial(N, 1/2), conditioned on both), a
# neutral trait (fitness 1) with K additive loci (h = 1/2: genotypic values
# 0, a/2, a), one locus per chromosome, so that every pair of loci
# recombines freely; no mutation. The starting population is drawn in
# Hardy-Weinberg and linkage equilibrium: allele frequencies p ~ U(0.1, 0.9),
# effects a ~ Normal, scaled so that the genic variance
# V_a = Σ 2pq (a/2)² is 1; environmental variance V_E = 1. The run then
# mates assortatively, with copula correlation ρ, for T generations.
#
# Reference. With free recombination and many loci of small effect
# (infinitesimal model), the additive variance of the parents of
# generation t + 1, given those of generation t, is
#     V_A' = (1/4) Var(A_s + A_d) + V_a / 2
#          = V_A (1 + m h²) / 2 + V_a / 2,   h² = V_A / V_P,
# where m is the phenotypic correlation of mates and the covariance of the
# breeding values of mates is m h² V_A for a linear regression of A on P
# (Bulmer 1980, ch. 9; Crow & Kimura 1970, §4.5). The disequilibrium
# component halves every generation, and the variance approaches
# V_a / (1 - m h²). Here V_A is the variance of the genotypic values of the
# parents, V_P that of their phenotypes, and m the realized correlation of
# the phenotypes of sire and dam over the matings, all from the record of
# the run (`var_genotypic`, `var_phenotypic`, `mate_correlation`); V_a is
# the genic variance of the start, reduced by drift by (1 - 1/(2N)) per
# generation. Sires and dams are drawn from the N_m males and N_f females
# of the cohort, whose breeding values have population variance
# V_A (1 - 1/N_m) and V_A (1 - 1/N_f) for a sample variance V_A, so the
# parental term carries the factor c = 1 - (1/N_m + 1/N_f)/2 (N_m from the
# record `males`). Finally, assortative mating also correlates the two
# alleles of a parent at the same locus, which lowers its heterozygosity
# and so the Mendelian sampling variance: the Mendelian term is
# (V_a - d_w)/2, where d_w, the within-locus part of the disequilibrium
# d = V_A - V_a, is d / K_eff with K_eff = (Σ x)² / Σ x², x = α² pq per
# locus (α = a/2), since the disequilibrium of loci i and j is in
# proportion to x_i x_j. K_eff is about 280 here; the term is 0.2 % of V_A
# at ρ = 0.8 and is left out of the infinitesimal model (K → ∞).
#
# Checks, each a z-score over replicates (|z| > ZMAX fails):
#   1. one-step: the mean over generations of V_A(t+1) minus its
#      prediction from generation t is 0;
#   2. trajectory: V_A at generations 2, 5 and T against the recursion
#      iterated from the start with the realized m and V_P;
#   3. the mean of the genotypic values does not change (assortative
#      mating changes no allele frequency in expectation).
# The recursion ignores the departure of A given P from linearity; the
# script prints the equilibrium of the recursion for comparison.

using FisherWright, Random, Statistics, Printf, Distributions

const ZMAX = 3.5
const N, K, T = 2_000, 1_000, 15
const CHR = fill(1_000, K)
const POS = UInt32.(cumsum(CHR) .- 500)          # one locus mid-chromosome
const VE = 1.0

# Effects scaled to genic variance 1 for frequencies p.
function loci_and_start(rng)
    p = rand(rng, Uniform(0.1, 0.9), K)
    a = randn(rng, K)
    a ./= sqrt(sum(2 .* p .* (1 .- p) .* (a ./ 2) .^ 2))
    loci = TraitLoci(POS, a; h=0.5, μ=0.0)
    sel = [POS[rand(rng, K) .< p] for _ = 1:2N]
    start = FisherWrightResult([UInt32[] for _ = 1:2N], UInt32.(cumsum(CHR)), UInt32[]; selected_haplotypes=sel)
    # genic variance of the starting sample
    q = [count(h -> insorted(POS[k], h), sel) / 2N for k = 1:K]
    x = (a ./ 2) .^ 2 .* q .* (1 .- q)
    Va = 2sum(x)
    return loci, start, Va, sum(x)^2 / sum(abs2, x)
end

function replicate(ρ, rng)
    loci, start, Va, Keff = loci_and_start(rng)
    f = TraitFitness(loci; σe=sqrt(VE), fitness=P -> 1.0)
    res = fisher_wright(N, T, CHR, 1.0; init=start, mutate=false, result=true, fitness=f,
        mating=Assortative(ρ), rng=rng)
    R = res.records
    VA = R.var_genotypic[1, 1, :]
    VP = R.var_phenotypic[1, 1, :]
    m = R.mate_correlation
    va(t) = Va * (1 - 1 / 2N)^(t - 1)
    c = [1 - (1 / R.males[t] + 1 / (R.census[t] - R.males[t])) / 2 for t = 1:T]
    mend(V, t) = (va(t) - (V - va(t)) / Keff) / 2
    pred1 = [c[t] * VA[t] * (1 + m[t] * VA[t] / VP[t]) / 2 + mend(VA[t], t) for t = 1:T-1]
    one = mean(VA[2:T] .- pred1)
    # iterated from the start, with the realized m and V_E
    traj = zeros(T)
    traj[1] = VA[1]
    for t = 1:T-1
        traj[t+1] = c[t] * traj[t] * (1 + m[t] * traj[t] / (traj[t] + VE)) / 2 + mend(traj[t], t)
    end
    return (; one, VA, traj, m, dmean=R.mean_genotypic[1, T] - R.mean_genotypic[1, 1], Va, Keff)
end

function check(label, x, expected, ok)
    m, s = mean(x), std(x) / sqrt(length(x))
    z = (m - expected) / s
    pass = abs(z) <= ZMAX
    @printf("  %-28s %10.5f ± %-8.5f %10.5f %8.2f  %s\n", label, m, s, expected, z, pass ? "pass" : "FAIL")
    return ok & pass
end

function main(reps=200, seed=2026)
    println("="^86)
    println(" FisherWright.jl validation of assortative mating against Bulmer's recursion")
    println(" N = $N, K = $K free loci, V_a = 1, V_E = $VE, T = $T, no mutation; reps = $reps, seed = $seed,")
    println(" threads = $(Threads.nthreads()). A check fails when |z| > $ZMAX")
    println("="^86)
    rng = Xoshiro(seed)
    ok = true
    for ρ in (0.5, 0.8, -0.5)
        t0 = time()
        rs = [replicate(ρ, rng) for _ = 1:reps]
        @printf("\n  ρ = %g  (%.0f s); realized mate correlation %.4f (generation 1) %.4f (generation %d); K_eff %.0f\n",
            ρ, time() - t0, mean(r.m[1] for r in rs), mean(r.m[T] for r in rs), T, mean(r.Keff for r in rs))
        @printf("  %-28s %21s %10s %8s\n", "", "simulated", "expected", "z")
        ok = check("one-step residual", [r.one for r in rs], 0.0, ok)
        for t in (2, 5, T)
            ok = check("V_A, generation $t - recursion", [r.VA[t] - r.traj[t] for r in rs], 0.0, ok)
        end
        ok = check("change of the mean", [r.dmean for r in rs], 0.0, ok)
        mT = mean(r.m[T] for r in rs)
        VaT = mean(r.Va for r in rs) * (1 - 1 / 2N)^(T - 1)
        eq = VaT                                  # V = Va / (1 - m V / (V + V_E)), by iteration
        for _ = 1:200
            eq = VaT / (1 - mT * eq / (eq + VE))
        end
        @printf("  V_A at generation %d: %.4f; recursion %.4f; its equilibrium for m = %.4f: %.4f (random mating: %.4f)\n",
            T, mean(r.VA[T] for r in rs), mean(r.traj[T] for r in rs), mT, eq, VaT)
    end
    println("\n", "="^86)
    println(ok ? " All checks passed." : " FAILED")
    println("="^86)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(parse.(Int, ARGS)...) ? 0 : 1)
end
