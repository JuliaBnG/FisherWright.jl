# Validation of selection in FisherWright.jl at many sites under recurrent
# mutation: the distribution of the derived-allele count at selected sites,
# including the fixed class (substitutions), after T generations from a
# population without variation. The package's engine is compared with an
# independent naive simulation of the same model at one site.
#
# With FisherWright installed, run from this directory:
#
#   julia -t 32 validate-selection-balance.jl [runs = 4000] [ref_reps = 400_000] [seed = 2026]
#
# Model. N diploids, two sexes (Binomial(N, 1/2) per generation, conditioned
# on both sexes), soft selection within each sex: each offspring draws its
# sire among the males and its dam among the females with probability
# proportional to fitness, with replacement. Genotype fitness 1, 1 + h*s,
# 1 + s at each site, multiplicative over sites. Each haplotype gains the
# derived allele at a site with probability 1 - exp(-μ) per generation (the
# engine draws Poisson(Kμ) hits uniformly over the K sites, which is
# Poisson(μ) per site); a hit on a site already carried has no effect.
# Order in a generation: fitness of the adults, mating, meiosis, mutation.
#
# Engine. `fisher_wright` with K sites, one per chromosome of 100 bp, so the
# sites recombine freely and nothing else does; neutral mutation is switched
# off in effect (mut_base = 1e12). Sites in one population share its
# pedigree, so the unit of replication is the run: class proportions are
# averaged over sites within a run, and standard errors are taken over runs.
#
# Reference. `naive_site` follows one site as genotype counts per individual,
# with no code shared with the package: sex by coin flips, parents by
# rejection sampling on fitness, one Bernoulli draw per transmitted gene and
# per mutation. It omits the other K - 1 sites, so it differs from the engine
# by their effect on this one: selection at unlinked sites adds variance in
# offspring number and lowers the effective size at every site (Robertson
# 1961), by a factor of about 1/(1 + 4C²) with C² the variance of relative
# fitness (Santiago and Caballero 1995, unlinked loci). The checked cases
# keep C² (printed as the variance of log fitness) at about 0.002 or less,
# where the effect is below 1%. The last case is not checked: with K = 50
# and C² ≈ 0.01 it shows the effect, against the reference at N and at
# N - 1, about 3% fewer individuals.
#
# The number of engine runs is `runs` × 50/K, so that each case has about
# the same number of sites. Each class is a z-score of the difference engine
# - reference; a checked case fails when |z| > ZMAX, and the script then
# exits with status 1. Deterministic values are printed for comparison only:
# the per-site derived frequency μ/(h|s|) and the mean fitness exp(-2Kμ)
# (Haldane-Muller), which hold for 4N h|s| ≫ 1 at mutation-selection
# balance; here 4N h|s| = 6 and 1.2.

using FisherWright, Random, Statistics, Printf

const ZMAX = 3.5

# Classes of the derived-allele count c in 0:2N; 2N is fixed (substituted).
classes(N) = [0:0, 1:1, 2:3, 4:9, 10:2N-1, 2N:2N]
label(r, N) = first(r) == last(r) ? (first(r) == 2N ? "fixed" : "$(first(r))") : "$(first(r))-$(last(r))"

# One site, `reps` replicates of T generations from c = 0; returns the final
# counts.
function naive_site(N, s, h, μ, T, reps, rng)
    w = (1.0, 1 + h * s, 1 + s)
    wmax = maximum(w)
    pm = 1 - exp(-μ)
    g, g2 = zeros(Int, N), zeros(Int, N)
    males, females = Int[], Int[]
    function parent(pool)
        while true
            i = rand(rng, pool)
            rand(rng) * wmax < w[g[i]+1] && return i
        end
    end
    out = zeros(Int, reps)
    for r = 1:reps
        fill!(g, 0)
        for _ = 1:T
            empty!(males)
            empty!(females)
            while isempty(males) || isempty(females)
                empty!(males)
                empty!(females)
                for i = 1:N
                    push!(rand(rng, Bool) ? males : females, i)
                end
            end
            for i = 1:N
                a, b = g[parent(males)], g[parent(females)]
                x = rand(rng) * 2 < a    # gene from the sire carries the derived allele
                y = rand(rng) * 2 < b
                x = x || rand(rng) < pm  # mutation of an ancestral gene
                y = y || rand(rng) < pm
                g2[i] = x + y
            end
            g, g2 = g2, g
        end
        out[r] = sum(g)
    end
    return out
end

function naive(N, s, h, μ, T, reps, seed)
    nt = Threads.nthreads()
    per = cld(reps, nt)
    tasks = [Threads.@spawn naive_site(N, s, h, μ, T, per, Xoshiro(hash((seed, k)))) for k = 1:nt]
    return reduce(vcat, fetch.(tasks))
end

# Derived-allele counts at the K sites of one engine run, and its mean
# fitness and variance of log fitness over the last generation's parents.
function engine_run(N, K, s, h, μ, T, rng)
    pos = collect(50:100:100K)
    loci = FitnessLoci(pos, s; h = h, μ = μ)
    res = fisher_wright(N, T, fill(100, K), 0.02; mut_base = 1e12, fitness = SiteFitness(loci),
        result = true, rng = rng)
    c = zeros(Int, K)
    for hap in res.selected_haplotypes, x in hap
        c[searchsortedfirst(loci.pos, x)] += 1
    end
    for x in res.selected_substitutions
        c[searchsortedfirst(loci.pos, x)] = 2N
    end
    return c, res.mean_fitness[end], res.var_log_fitness[end]
end

function engine(N, K, s, h, μ, T, runs, seed)
    nt = Threads.nthreads()
    per = cld(runs, nt)
    tasks = [Threads.@spawn begin
        rng = Xoshiro(hash((seed, k)))
        [engine_run(N, K, s, h, μ, T, rng) for _ = 1:per]
    end for k = 1:nt]
    return reduce(vcat, fetch.(tasks))
end

# (N, K, s, h, μ, T, checked)
const CASES = [
    (30, 10, -0.1, 0.5, 2e-3, 300, true),
    (30, 10, -0.1, 0.1, 1e-3, 300, true),
    (30, 50, 0.0, 0.5, 2e-3, 300, true),
    (30, 20, 0.05, 0.5, 5e-4, 300, true),
    (30, 50, -0.1, 0.5, 2e-3, 300, false),
]

# One table row: z-score of the engine's mean of `pe` (one value per run)
# against the reference's mean of `pr` (one value per replicate).
function row(name, pe, pr, checked)
    se_e, se_r = sqrt(var(pe) / length(pe)), sqrt(var(pr) / length(pr))
    z = (mean(pe) - mean(pr)) / sqrt(se_e^2 + se_r^2)
    pass = abs(z) <= ZMAX
    @printf("  %-10s %10.5f ± %-8.5f %10.5f ± %-8.5f %8.2f  %s\n", name, mean(pe), se_e, mean(pr), se_r, z,
        checked ? (pass ? "pass" : "FAIL") : "(not checked)")
    return pass || !checked
end

function balance_main(runs = 4000, ref_reps = 400_000, seed = 2026)
    println("="^88)
    println(" FisherWright.jl validation of selection under recurrent mutation")
    println(" derived-allele count per selected site after T generations; engine runs = $runs,")
    println(" naive reference reps = $ref_reps, seed = $seed, threads = $(Threads.nthreads())")
    println(" A check fails when |z| > $ZMAX")
    println("="^88)
    ok = true
    for (k, (N, K, s, h, μ, T, checked)) in enumerate(CASES)
        t = @elapsed begin
            eng = engine(N, K, s, h, μ, T, runs * 50 ÷ K, seed + k)
            ref = naive(N, s, h, μ, T, ref_reps, seed + 100k)
        end
        @printf("\n[%d/%d] N = %d, K = %d sites, s = %g, h = %g, μ = %g, T = %d, %d runs  (%.0f s)%s\n",
            k, length(CASES), N, K, s, h, μ, T, length(eng), t, checked ? "" : "  -- not checked")
        @printf("  %-10s %21s %21s %8s\n", "count", "engine", "naive reference", "z")
        for r in classes(N)
            ok &= row(label(r, N), [count(in(r), c) / K for (c, _, _) in eng], [x in r for x in ref], checked)
        end
        qe = [mean(c) / 2N for (c, _, _) in eng]
        ok &= row("mean q", qe, ref ./ 2N, checked)
        C² = mean(v for (_, _, v) in eng)
        @printf("  for comparison: var(log w) of the parents %.4f", C²)
        if s < 0
            @printf("; μ/(h|s|) = %.4f; mean fitness %.4f against exp(-2Kμ) = %.4f",
                μ / (h * abs(s)), mean(w for (_, w, _) in eng), exp(-2K * μ))
        end
        println()
        if !checked
            # The reference with one individual fewer, as a scale for the
            # loss of effective size predicted by 1/(1 + 4C²).
            ref1 = naive(N - 1, s, h, μ, T, ref_reps, seed + 200k)
            @printf("  predicted Ne/N = 1/(1 + 4C²) = %.3f; reference at N = %d: mean q %.5f ± %.5f, fixed %.5f ± %.5f\n",
                1 / (1 + 4C²), N - 1, mean(ref1 ./ 2(N - 1)), std(ref1 ./ 2(N - 1)) / sqrt(length(ref1)),
                mean(ref1 .== 2(N - 1)), std(ref1 .== 2(N - 1)) / sqrt(length(ref1)))
            @printf("  engine against the reference at N: mean q %+.1f%%, fixed %+.1f%%\n",
                100 * (mean(qe) / mean(ref ./ 2N) - 1),
                100 * (mean(count(==(2N), c) / K for (c, _, _) in eng) / mean(ref .== 2N) - 1))
        end
    end
    println("\n", "="^88)
    println(ok ? " All checks passed." : " FAILED")
    println("="^88)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(balance_main(parse.(Int, ARGS)...) ? 0 : 1)
end
