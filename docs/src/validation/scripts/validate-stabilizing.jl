# Genetic variance under mutation, stabilizing selection and drift:
# TraitFitness with Stabilizing against a matched SLiM model
# (validate-stabilizing.slim), with the diffusion of the stochastic house of
# cards for comparison. Driven by validate-stabilizing.sh:
#
#   julia -t 32 --project=. validate-stabilizing.jl alpha
#   julia -t 32 --project=. validate-stabilizing.jl engine OUTDIR RUNS SEED
#   julia --project=. validate-stabilizing.jl compare OUTDIR
#
# Model. N = 500 diploids; in every generation the roles of sire and dam
# are assigned anew (binomial sex, conditioned on both, the default of
# FisherWright); every offspring draws its sire among the males and its dam
# among the females with probability proportional to fitness, with
# replacement. Phenotype P = G + e, e ~ N(0, 1) drawn anew; fitness
# exp(-P²/(2ω²)), ω² = 49, optimum 0. K = 8,000 additive trait loci
# (h = 1/2), evenly spaced on 20 chromosomes of 100 Mb and 1 Morgan, in
# four interleaved classes of 2,000 with substitution effect α (each copy
# of the derived allele adds α to G) of fixed size and a random sign, drawn
# anew for each run; the sizes give B = 2Ne α²/V_s ≈ 1, 4, 16 and 64 at
# V_s = 51. Mutation to the derived allele at μ = 10^-3/(4N) per site per
# haplotype per generation, biallelic sites, no back mutation; a fixed site
# stays fixed. No neutral mutation.
#
# Protocol. From a population without variation, BURN = 20N generations,
# then a sample every GAP = 2N generations, NS = 21 samples in all (30,000
# generations per run). Per sample: the genic variance Σ 2pqα² and the
# number of fixed sites of each class, the mean genotypic value, and the
# mean over the preceding window (GAP generations; the last 2N of the
# burn-in for the first sample) of V_P, V_G and the variance of log fitness
# of the parents. The samples of a run are averaged, and standard errors
# are taken over runs, so the correlation between the samples of one run
# does not matter.
#
# Comparison. Each statistic, FisherWright against SLiM, as a z-score of
# the difference with the standard errors of both; the genic variance of
# a class per available (not fixed) site, in units of 4μV_s with the V_s
# pooled over both simulators. A comparison fails when |z| > 3.5, and the
# script then exits with status 1. The mean genotypic value is also
# compared with 0, its exact expectation (the effects have random signs).
#
# Theory, for comparison. The trait loci are biallelic with forward
# mutation only, so every allele is unique: this is the house-of-cards
# family of models (Turelli 1984; Bürger 2000), not the Gaussian (Kimura
# 1965; Lande 1975), which needs a continuum of alleles at each locus. With
# the mean at the optimum, a site with substitution effect α is under
# underdominant selection, Δx = -(α²/(2V_s)) x(1 - x)(1 - 2x), V_s = ω² +
# V_P (Robertson 1956), with drift variance x(1 - x)/(2Ne). For new
# mutations at rate 2Nμ per site, the stationary density is
#     φ(x) = θ (S(1) - S(x)) / (S(1) x (1 - x) ψ(x)),
#     ψ(x) = exp(B x(1 - x)),  S(x) = ∫_0^x ψ,  B = 2Ne α²/V_s,  θ = 4Ne μ,
# and the expected genic variance of a site, in units of 4μV_s, is B·I(B)
# with I = ∫_0^1 (S(1) - S(x)) / (S(1) ψ(x)) dx: B/2 (neutral) for B → 0
# and 1 (Turelli's house of cards) for B → ∞, above 1 in between. Bürger's
# stochastic house-of-cards formula is (B/2)/(1 + B/2) in these units. This
# single-locus diffusion leaves out the linkage disequilibrium that
# selection builds: carriers of an allele of effect α are found on
# backgrounds shifted by about -αk/(1 + k), k = V_G/V_s (selection shifts
# them by -αk each generation, unlinked loci halve the shift), which
# weakens the selection on the allele to about B(1 - k)². The diffusion
# with this correction is printed as well. Both leave out the drift of the
# mean about the optimum (variance about V_s/(2Ne); Lande 1976).

using FisherWright, Random, Statistics, Printf

const N = 500
const NCHR, CHRLEN = 20, 100_000_000
const KC = 2_000
const BNOM = [1.0, 4.0, 16.0, 64.0]
const NC = length(BNOM)
const K = NC * KC
const VS0 = 51.0                   # nominal V_s for the effect sizes
const OMEGA2, SIGMAE = 49.0, 1.0
const THETA = 1e-3                 # 4Nμ per site
const MU = THETA / 4N
const BURN, GAP, NS = 20N, 2N, 21
const ZMAX = 3.5
const COLS = ["va1", "va2", "va3", "va4", "fx1", "fx2", "fx3", "fx4", "gbar", "vp", "vg", "c2"]

function effective_size()
    w = [Float64(binomial(big(N), k)) for k = 1:N-1]
    a = sum(w[k] * (1 / 4k + 1 / 4(N - k)) for k = 1:N-1) / sum(w)
    return 1 / a + 0.5
end
const NE = effective_size()
const ALPHA = sqrt.(BNOM .* VS0 ./ 2NE)    # |α| per class

# B·I(B): the expected genic variance per site in units of 4μV_s, by the
# trapezoidal rule on a fine grid (the integrand is smooth).
function shc_diffusion(B; m = 200_000)
    x = range(0, 1, length = m + 1)
    ψ = exp.(B .* x .* (1 .- x))
    S = zeros(m + 1)
    for k = 2:m+1
        S[k] = S[k-1] + (ψ[k-1] + ψ[k]) / 2m
    end
    g = (S[end] .- S) ./ (S[end] .* ψ)
    I = (sum(g) - (g[1] + g[end]) / 2) / m
    return B * I
end

# The NS samples of one FisherWright run, one row each (the columns COLS).
function engine_run(seed)
    rng = Xoshiro(seed)
    pos = collect(range(CHRLEN ÷ (2K ÷ NCHR), step = NCHR * CHRLEN ÷ K, length = K))
    cls = [mod1(k, NC) for k = 1:K]
    α = ALPHA[cls] .* rand(rng, (-1.0, 1.0), K)
    loci = TraitLoci(pos, 2 .* α; h = 0.5, μ = MU)
    f = TraitFitness(loci; σe = SIGMAE, fitness = Stabilizing(0.0, OMEGA2))
    run(nt, init) = fisher_wright(N, nt, fill(CHRLEN, NCHR), 0.02; M = 1e8, mut_base = 1e12, fitness = f,
        result = true, init = init, rng = rng)
    rows = zeros(NS, length(COLS))
    res = run(BURN, nothing)
    for s = 1:NS
        s > 1 && (res = run(GAP, res))
        c = zeros(Int, K)
        for hap in res.selected_haplotypes, x in hap
            c[searchsortedfirst(pos, x)] += 1
        end
        for x in res.selected_substitutions
            k = searchsortedfirst(pos, x)
            c[k] = 2N
            rows[s, 4+cls[k]] += 1
        end
        for k = 1:K
            q = c[k] / 2N
            rows[s, cls[k]] += 2q * (1 - q) * α[k]^2
        end
        rows[s, 9] = sum(c[k] / 2N * 2α[k] for k = 1:K)
        w = (s == 1 ? 2N : GAP)
        g = length(res.var_phenotypic)-w+1:length(res.var_phenotypic)
        rows[s, 10] = mean(res.var_phenotypic[g])
        rows[s, 11] = mean(res.var_genotypic[g])
        rows[s, 12] = mean(res.var_log_fitness[g])
    end
    return rows
end

function engine(out, runs, seed)
    tabs = Vector{Matrix{Float64}}(undef, runs)
    Threads.@threads for r = 1:runs
        tabs[r] = engine_run(hash((seed, r)))
    end
    open(joinpath(out, "engine.tsv"), "w") do io
        println(io, "# FisherWright ", pkgversion(FisherWright), ", Julia ", VERSION, ", threads ",
            Threads.nthreads(), ", ", gethostname(), "; $runs runs, seed $seed")
        for r = 1:runs, s = 1:NS
            println(io, join(["STAB"; string(r); string(BURN + (s - 1) * GAP); string.(tabs[r][s, :])], '\t'))
        end
    end
end

# Per-run means of the columns, from rows "STAB  run  generation  COLS...".
function runmeans(file)
    acc = Dict{String,Vector{Vector{Float64}}}()
    for line in eachline(file)
        f = split(line, '\t')
        f[1] == "STAB" || continue
        push!(get!(acc, f[2], Vector{Float64}[]), parse.(Float64, f[4:end]))
    end
    return [mean(v) for v in values(acc)]
end

function compare(out)
    fw, sl = runmeans(joinpath(out, "engine.tsv")), runmeans(joinpath(out, "slim.tsv"))
    first = readline(joinpath(out, "engine.tsv"))
    println(startswith(first, "# ") ? first[3:end] :
            "engine.tsv has no provenance line (written by an earlier version of this script)")
    @printf("N = %d (Ne = %.2f), %d additive loci on %d × %d Mb (1 Morgan each), 4Nμ = %g per site; ",
        N, NE, K, NCHR, CHRLEN ÷ 1_000_000, THETA)
    @printf("fitness exp(-P²/(2·%g)), σe = %g\n", OMEGA2, SIGMAE)
    @printf("|α| = %s; burn-in %d generations, then %d samples %d generations apart\n",
        join((@sprintf("%.4f", a) for a in ALPHA), ", "), BURN, NS, GAP)
    @printf("runs: FisherWright %d, SLiM %d; per-run means of the samples, SE over runs\n\n", length(fw), length(sl))
    vs = OMEGA2 + mean([x[10] for x in [fw; sl]])
    k = mean([x[11] for x in [fw; sl]]) / vs
    failed = String[]
    function row(name, a, b; fmt = "%.4f", ratio = true)
        ma, sa = mean(a), std(a) / sqrt(length(a))
        mb, sb = mean(b), std(b) / sqrt(length(b))
        z = sa + sb > 0 ? (ma - mb) / sqrt(sa^2 + sb^2) : 0.0
        f = Printf.Format(fmt * " ± " * fmt)
        r = ratio && sa + sb > 0 ? @sprintf("%8.4f", ma / mb) : "       —"
        zs = sa + sb > 0 ? @sprintf("%7.2f", z) : "      —"
        @printf("%-22s %22s %22s %s %s\n", name, Printf.format(f, ma, sa), Printf.format(f, mb, sb), r, zs)
        abs(z) > ZMAX && push!(failed, name)
        return ma
    end
    @printf("%-22s %22s %22s %8s %7s\n", "statistic", "FisherWright", "SLiM", "ratio", "z")
    th = zeros(NC, 4)
    for j = 1:NC
        per(x) = x[j] / ((KC - x[4+j]) * 4MU * vs)
        m = row("V_a/(K_a 4μV_s), B=$(Int(BNOM[j]))", per.(fw), per.(sl))
        B = 2NE * ALPHA[j]^2 / vs
        th[j, :] = [B, shc_diffusion(B), shc_diffusion(B * (1 - k)^2) / (1 - k)^2, B / (2 + B)]
    end
    for j = 1:NC
        row("fixed per site, B=$(Int(BNOM[j]))", [x[4+j] / KC for x in fw], [x[4+j] / KC for x in sl]; fmt = "%.5f")
    end
    row("mean G", [x[9] for x in fw], [x[9] for x in sl]; ratio = false)
    row("V_P", [x[10] for x in fw], [x[10] for x in sl])
    row("V_G", [x[11] for x in fw], [x[11] for x in sl])
    row("var log w", [x[12] for x in fw], [x[12] for x in sl]; fmt = "%.6f")
    for (name, x) in (("FisherWright", fw), ("SLiM", sl))
        g = [y[9] for y in x]
        z = mean(g) / (std(g) / sqrt(length(g)))
        @printf("mean G against 0 (exact), %s: z = %.2f\n", name, z)
        abs(z) > ZMAX && push!(failed, "mean G of $name against 0")
    end
    @printf("\nTheory, for comparison (V_s = ω² + V_P = %.3f, k = V_G/V_s = %.4f):\n", vs, k)
    @printf("%-8s %8s %12s %12s %12s %12s\n", "class", "B", "diffusion", "LD-corr.", "SHC", "FW/LD-corr.")
    for j = 1:NC
        m = mean(x[j] / ((KC - x[4+j]) * 4MU * vs) for x in fw)
        @printf("%-8d %8.3f %12.4f %12.4f %12.4f %12.4f\n", j, th[j, 1], th[j, 2], th[j, 3], th[j, 4], m / th[j, 3])
    end
    isempty(failed) || (println("\nFAILED: ", join(failed, "; ")); exit(1))
end

if abspath(PROGRAM_FILE) == @__FILE__
    cmd = ARGS[1]
    if cmd == "alpha"
        println(join(ALPHA, ","))
    elseif cmd == "engine"
        engine(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]))
    elseif cmd == "compare"
        compare(ARGS[2])
    else
        error("unknown command $cmd")
    end
end
