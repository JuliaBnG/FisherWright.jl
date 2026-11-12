# Validation of selection on a trait (TraitFitness) in FisherWright.jl
# against quantitative-genetic theory.
#
# With FisherWright installed, run from this directory:
#
#   julia -t 8 validate-trait.jl [reps = 1000] [seed = 2026]
#
# 1. Neutral genetic variance under mutation and drift, against its exact
#    expectation. N diploids, two sexes (Binomial(N, 1/2), conditioned on
#    both), random mating; K = 20,000 trait loci on 10 chromosomes of 1 Mb,
#    each additive (h = 1/2) with effect a ~ Normal(0, σa²) drawn anew for
#    every replicate, mutation μ = 1e-5 per site per haplotype per
#    generation (U = Kμ = 0.2), fitness constant, a population without
#    variation at the start. For allele indicators, the sample variance of
#    the dosage among the N individuals has expectation 2H_B - H_W, H_B and
#    H_W the heterozygosity of two genes in different individuals and in one
#    individual. With independent effects of zero mean the loci add, and the
#    allele substitution effect is α = a/2, so
#        E[var G] = E[α²] (2H_B - H_W),
#    where H_B and H_W follow the two-sex recursion with 2U new alleles per
#    pair of genes and generation (infinite sites; recurrent hits at a site
#    have probability about 4Neμ = 0.002 here). At equilibrium this is
#    2 Ne V_m with Ne = 1/a, V_m = 2U E[α²] (Lynch and Hill 1986). With
#    σe > 0 also E[var P] = E[var G] + σe².
# 2. For comparison only: the response to directional selection,
#    fitness exp(βP). For a Gaussian G and P = G + e, the fitness-weighted
#    mean of G is shifted by β V_G and its variance is unchanged (no Bulmer
#    effect), so the response per generation is about β V_G; it is exact in
#    expectation for the realized parents and their fitness (CI test
#    "Selection on a trait"), and the comparison shows how close the
#    finite-locus, finite-N process stays to the Gaussian prediction.
#    Within each sex the fitness-weighted mean of G over n parents is a
#    ratio estimator; for exp(βP) with Gaussian G and P its first-order
#    bias is -β V_G e^{β² V_P} / n, so the expected ratio is about
#    1 - e^{β² V_P} E[(1/N_m + 1/N_f)/2], which the script also prints.
#
# Each comparison in part 1 is a z-score over replicates; a check fails when
# |z| > ZMAX, and the script then exits with status 1.

using FisherWright, Random, Statistics, Printf, Distributions

const ZMAX = 3.5
const K, MU, SIGMA_A = 20_000, 1e-5, 0.2
const CHR = fill(1_000_000, 10)

function coalescence(N)
    w = [Float64(binomial(big(N), k) / big(2)^N) for k = 1:N-1]
    return sum(w[k] * (1 / 4k + 1 / 4(N - k)) for k = 1:N-1) / sum(w)
end

# E[2H_B - H_W] for the parents of generations 1..T from a population
# without variation.
function expected_dosage_var(N, T, U)
    a = coalescence(N)
    hB = hW = 0.0
    out = zeros(T)
    for g = 1:T
        out[g] = 2hB - hW
        hB, hW = (1 - a) * hB + a / 2 * hW + 2U, hB + 2U
    end
    return out
end

function neutral_case(N, T, σe, reps, rng)
    VG, VP = zeros(reps, T), zeros(reps, T)
    for r = 1:reps
        loci = trait_loci(CHR, K; a = Normal(0, SIGMA_A), μ = MU, rng = rng)
        f = TraitFitness(loci; σe = σe, fitness = P -> 1.0)
        res = fisher_wright(N, T, CHR, 0.02; mut_base = 1e12, fitness = f, result = true, rng = rng)
        VG[r, :] .= res.var_genotypic
        VP[r, :] .= res.var_phenotypic
    end
    return VG, VP
end

function main(reps = 1000, seed = 2026)
    println("="^80)
    println(" FisherWright.jl validation of selection on a trait")
    println(" K = $K additive loci, a ~ Normal(0, $(SIGMA_A)²), μ = $MU, reps = $reps, seed = $seed,")
    println(" threads = $(Threads.nthreads()). A check fails when |z| > $ZMAX")
    println("="^80)
    rng = Xoshiro(seed)
    ok = true
    Eα² = SIGMA_A^2 / 4
    println("\n[1/2] Neutral genetic variance against its exact expectation")
    for (N, T, σe) in ((20, 200, 0.0), (50, 500, 0.0), (50, 100, 1.0))
        t = @elapsed VG, VP = neutral_case(N, T, σe, reps, rng)
        E = Eα² .* expected_dosage_var(N, T, K * MU)
        @printf("\n  N = %d, T = %d, σe = %g  (%.0f s); equilibrium 2NeVm = %.4f with Ne = 1/a = %.2f\n",
            N, T, σe, t, 2 / coalescence(N) * 2K * MU * Eα², 1 / coalescence(N))
        @printf("  %-12s %21s %10s %8s\n", "generation", "mean var G", "expected", "z")
        for g in filter(<=(T), [2, 5, 10, 20, 50, 100, 200, 500])
            m, s = mean(VG[:, g]), std(VG[:, g]) / sqrt(reps)
            z = (m - E[g]) / s
            pass = abs(z) <= ZMAX
            ok &= pass
            @printf("  %-12d %10.5f ± %-8.5f %10.5f %8.2f  %s\n", g, m, s, E[g], z, pass ? "pass" : "FAIL")
        end
        if σe > 0
            d = VP[:, T] .- VG[:, T]
            z = (mean(d) - σe^2) / (std(d) / sqrt(reps))
            pass = abs(z) <= ZMAX
            ok &= pass
            @printf("  %-12s %10.5f ± %-8.5f %10.5f %8.2f  %s\n", "var P - var G", mean(d),
                std(d) / sqrt(reps), σe^2, z, pass ? "pass" : "FAIL")
        end
    end

    println("\n[2/2] Response to directional selection exp(βP), for comparison")
    N, burn, T, β, σe = 50, 500, 30, 0.5, 1.0
    ratios, small = Float64[], Float64[]
    w = [Float64(binomial(big(N), k)) for k = 1:N-1]
    Einv = sum(w[k] * (1 / k + 1 / (N - k)) / 2 for k = 1:N-1) / sum(w)
    for r = 1:min(reps, 200)
        loci = trait_loci(CHR, K; a = Normal(0, SIGMA_A), μ = MU, rng = rng)
        base = fisher_wright(N, burn, CHR, 0.02; mut_base = 1e12, fitness = TraitFitness(loci; fitness = P -> 1.0),
            result = true, rng = rng)
        sel = fisher_wright(N, T, CHR, 0.02; mut_base = 1e12, init = base, result = true, rng = rng,
            fitness = TraitFitness(loci; σe = σe, fitness = Directional(β)))
        g = burn+1:burn+T
        R = sel.mean_genotypic[g[end]] - sel.mean_genotypic[g[1]]
        VG = sel.var_genotypic[g[1:end-1]]
        push!(ratios, R / (β * sum(VG)))
        push!(small, 1 - sum(VG .* exp.(β^2 .* sel.var_phenotypic[g[1:end-1]])) / sum(VG) * Einv)
    end
    @printf("  N = %d, %d generations of selection after %d neutral, β = %g, σe = %g, %d replicates\n",
        N, T, burn, β, σe, length(ratios))
    @printf("  cumulative response / Σ β V_G = %.3f ± %.3f\n", mean(ratios), std(ratios) / sqrt(length(ratios)))
    @printf("  Gaussian prediction 1; with the finite-N bias of the weighted mean %.3f\n", mean(small))

    println("\n", "="^80)
    println(ok ? " All checks passed." : " FAILED")
    println("="^80)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(parse.(Int, ARGS)...) ? 0 : 1)
end
