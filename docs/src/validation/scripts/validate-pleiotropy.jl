# Validation of pleiotropy (TraitFitness with several traits on shared loci)
# in FisherWright.jl against quantitative-genetic theory.
#
# With FisherWright installed, run from this directory:
#
#   julia -t 8 validate-pleiotropy.jl [reps = 1000] [seed = 2026]
#
# 1. Neutral genetic covariance matrix under mutation and drift, against its
#    exact expectation. N diploids, two sexes (Binomial(N, 1/2), conditioned
#    on both), random mating; K = 20,000 trait loci on 10 chromosomes of
#    1 Mb, each additive (h = 1/2) for two traits, with effects
#    (a1, a2) ~ MvNormal(0, Σa) drawn anew for every replicate,
#    Σa = 0.2² [1 ρ; ρ 1], ρ = -0.5; mutation μ = 1e-5 per site per
#    haplotype per generation (U = Kμ = 0.2), fitness constant, a population
#    without variation at the start. As for one trait (validate-trait.jl),
#    the sample covariance of the dosages of one site among the N
#    individuals has expectation 2H_B - H_W, and with effects of zero mean
#    independent across sites the sites add, so for traits s and t
#        E[cov(G_s, G_t)] = E[α_s α_t] (2H_B - H_W),  E[α_s α_t] = Σa[s, t]/4,
#    the allele substitution effects being α = a/2: the pleiotropic
#    covariance of the effects becomes a genetic covariance. With an
#    environmental covariance Σe, also E[S_P - S_G] = Σe for the sample
#    covariance matrices S_P and S_G, since e is independent of G.
# 2. For comparison only: the direct and the correlated response to
#    directional selection on the first trait, fitness exp(β'P) with
#    β = (β1, 0). For Gaussian G and P = G + e, the fitness-weighted mean of
#    G is shifted by Gβ (Lande 1979), so trait 2 responds by β1 G21 per
#    generation, through the shared loci only. The ratio of the mean
#    cumulative response to the mean of Σ_g (G_g β) over replicates (a ratio
#    of means, since Σ_g G21 of one replicate can be near 0) is about 1, less the finite-N bias of the
#    weighted mean within each sex, first order -e^{β'Pβ} E[(1/N_m + 1/N_f)/2]
#    of the response for either trait, which the script also prints.
#
# Each comparison in part 1 is a z-score over replicates; a check fails when
# |z| > ZMAX, and the script then exits with status 1.

using FisherWright, Random, Statistics, Printf, Distributions

const ZMAX = 3.5
const K, MU, SIGMA_A, RHO = 20_000, 1e-5, 0.2, -0.5
const SIGMA_A_MAT = SIGMA_A^2 .* [1.0 RHO; RHO 1.0]
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

# Sample G and P covariance matrices of every generation, reps × 2 × 2 × T.
function neutral_case(N, T, Σe, reps, rng)
    SG, SP = zeros(reps, 2, 2, T), zeros(reps, 2, 2, T)
    for r = 1:reps
        loci = trait_loci(CHR, K; a = MvNormal(zeros(2), SIGMA_A_MAT), μ = MU, rng = rng)
        f = TraitFitness(loci; Σe = Σe, fitness = P -> 1.0)
        res = fisher_wright(N, T, CHR, 0.02; mut_base = 1e12, fitness = f, result = true, rng = rng)
        SG[r, :, :, :] .= res.records.var_genotypic
        SP[r, :, :, :] .= res.records.var_phenotypic
    end
    return SG, SP
end

function check(label, x, expected, ok)
    m, s = mean(x), std(x) / sqrt(length(x))
    z = (m - expected) / s
    pass = abs(z) <= ZMAX
    @printf("  %-16s %10.5f ± %-8.5f %10.5f %8.2f  %s\n", label, m, s, expected, z, pass ? "pass" : "FAIL")
    return ok & pass
end

function main(reps = 1000, seed = 2026)
    println("="^80)
    println(" FisherWright.jl validation of pleiotropy")
    println(" K = $K additive loci, (a1, a2) ~ MvNormal(0, $(SIGMA_A)² [1 ρ; ρ 1]), ρ = $RHO, μ = $MU,")
    println(" reps = $reps, seed = $seed, threads = $(Threads.nthreads()). A check fails when |z| > $ZMAX")
    println("="^80)
    rng = Xoshiro(seed)
    ok = true
    Eαα = SIGMA_A_MAT ./ 4
    zero2 = zeros(2, 2)
    println("\n[1/2] Neutral genetic covariance matrix against its exact expectation")
    for (N, T, Σe) in ((20, 200, zero2), (50, 500, zero2), (50, 100, [1.0 0.4; 0.4 0.5]))
        t = @elapsed SG, SP = neutral_case(N, T, Σe, reps, rng)
        E = expected_dosage_var(N, T, K * MU)
        @printf("\n  N = %d, T = %d, Σe = %s  (%.0f s)\n", N, T, Σe == zero2 ? "0" : string(Σe), t)
        @printf("  %-16s %21s %10s %8s\n", "generation", "mean", "expected", "z")
        for g in filter(<=(T), [2, 10, 50, 100, 200, 500])
            for (s, u, name) in ((1, 1, "G11"), (1, 2, "G12"), (2, 2, "G22"))
                ok = check("$g  $name", SG[:, s, u, g], Eαα[s, u] * E[g], ok)
            end
        end
        if Σe != zero2
            for (s, u, name) in ((1, 1, "P11 - G11"), (1, 2, "P12 - G12"), (2, 2, "P22 - G22"))
                ok = check("$T  $name", SP[:, s, u, T] .- SG[:, s, u, T], Σe[s, u], ok)
            end
        end
    end

    println("\n[2/2] Direct and correlated response to exp(β'P), β = (β1, 0), for comparison")
    N, burn, T, β1, σe = 50, 500, 30, 0.5, 1.0
    β = [β1, 0.0]
    R1, R2, D1, D2, small = Float64[], Float64[], Float64[], Float64[], Float64[]
    w = [Float64(binomial(big(N), k)) for k = 1:N-1]
    Einv = sum(w[k] * (1 / k + 1 / (N - k)) / 2 for k = 1:N-1) / sum(w)
    for r = 1:min(reps, 200)
        loci = trait_loci(CHR, K; a = MvNormal(zeros(2), SIGMA_A_MAT), μ = MU, rng = rng)
        base = fisher_wright(N, burn, CHR, 0.02; mut_base = 1e12, result = true, rng = rng,
            fitness = TraitFitness(loci; fitness = P -> 1.0))
        sel = fisher_wright(N, T, CHR, 0.02; mut_base = 1e12, init = base, result = true, rng = rng,
            fitness = TraitFitness(loci; σe = σe, fitness = Directional(β)))
        g = burn+1:burn+T
        R = sel.records.mean_genotypic[:, g[end]] .- sel.records.mean_genotypic[:, g[1]]
        Gβ = sum(sel.records.var_genotypic[:, :, k] * β for k in g[1:end-1])
        push!(R1, R[1]); push!(D1, Gβ[1])
        push!(R2, R[2]); push!(D2, Gβ[2])
        Pβ = [β' * sel.records.var_phenotypic[:, :, k] * β for k in g[1:end-1]]
        G11 = sel.records.var_genotypic[1, 1, g[1:end-1]]
        push!(small, 1 - sum(G11 .* exp.(Pβ)) / sum(G11) * Einv)
    end
    n = length(R1)
    # Ratio of means and its delta-method standard error.
    ratio(R, D) = (q = mean(R) / mean(D); (q, std(R .- q .* D) / sqrt(n) / abs(mean(D))))
    q1, s1 = ratio(R1, D1)
    q2, s2 = ratio(R2, D2)
    @printf("  N = %d, %d generations of selection after %d neutral, β = (%g, 0), σe = %g, %d replicates\n",
        N, T, burn, β1, σe, n)
    @printf("  trait 1 (direct):     mean response %8.4f, mean Σ (Gβ)_1 %8.4f, ratio %.3f ± %.3f\n",
        mean(R1), mean(D1), q1, s1)
    @printf("  trait 2 (correlated): mean response %8.4f, mean Σ (Gβ)_2 %8.4f, ratio %.3f ± %.3f\n",
        mean(R2), mean(D2), q2, s2)
    @printf("  Gaussian prediction 1; with the finite-N bias of the weighted mean %.3f\n", mean(small))

    println("\n", "="^80)
    println(ok ? " All checks passed." : " FAILED")
    println("="^80)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(parse.(Int, ARGS)...) ? 0 : 1)
end
