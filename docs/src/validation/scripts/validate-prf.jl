# Site frequency spectrum at weakly selected sites against the Poisson
# random field (Sawyer and Hartl 1992); docs/dev/roadmap.md, validation 6.
#
#   julia -t 32 --project=. docs/src/validation/scripts/validate-prf.jl [runs = 1000] [seed = 2026]
#
# Model. N = 1,000 diploids, two sexes (binomial, the default), T = 20N =
# 20,000 generations from a population without variation; the slowest mode
# of a neutral polymorphism decays as exp(-t/2N), so e^-10 of the start is
# left. K = 500,000 selected sites, evenly spaced on 20 chromosomes of
# 100 Mb and 1 Morgan, in five interleaved classes of 100,000 sites with
# scaled selection α = 2Ne·s ∈ {-10, -2, 0, 2, 10}, h = 1/2 (fitness 1,
# 1 + s/2, 1 + s; selection within each sex), and derived-allele mutation
# at μ = 10^-4/(4N) per site per haplotype per generation. Neutral
# mutation is switched off in effect (mut_base = 10^12).
#
# Infinite sites. The sites of FisherWright are finite: a hit on a site
# already carried has no effect, and a fixed site stays fixed. With
# 4Nμ = 10^-4 per site, recurrent hits change the spectrum by about 10^-4,
# and the influx of new mutations is that of the sites not fixed by the end
# (K_a), which differs from that of the past 4N generations by under 0.1%
# (α = 10: about 0.5% of the sites fix in T generations). About 400 sites
# segregate at a time, so interference between them is negligible.
#
# Statistic. The expected spectrum ξ_i (i = 1, ..., n - 1) of a sample of
# n = 20 genes drawn without replacement from the 2N, per available site,
# computed exactly (hypergeometric) from the derived-allele counts of each
# run, so that it carries no sampling noise beyond the population's.
# Reported as ξ_i/θ in bins of i, with θ = 4Ne·μ; the unit of replication
# is the run, and standard errors are over runs.
#
# Reference. The Poisson random field,
#   ξ_i = θ ∫ C(n, i) x^i (1 - x)^(n-i) f(x) dx,
#   f(x) = (1 - e^{-α(1-x)}) / ((1 - e^{-α}) x (1 - x))      (θ/x for α = 0),
# the diffusion limit of the stationary spectrum under genic selection with
# Δx = (s/2) x (1 - x) and drift variance x(1 - x)/(2Ne). Ne = 1/a + 1/2
# with a = E[1/(4Nm) + 1/(4Nf)] over binomial sex (999.5 here), and
# s = α/(2Ne). The diffusion is a limit: at N = 1,000 its error is of order
# 1/N, about the resolution here. A bin fails when |z| > 3.5, and the
# script then exits with status 1.

using FisherWright, Random, Statistics, Printf

const N, T = 1_000, 20_000
const NCHR, CHRLEN = 20, 100_000_000
const K = 500_000
const ALPHA = [-10.0, -2.0, 0.0, 2.0, 10.0]
const THETA = 1e-4                 # 4Nμ per site
const MU = THETA / 4N
const NS = 20
const BINS = [1:1, 2:2, 3:4, 5:9, 10:19]
const ZMAX = 3.5

# Ne = 1/a + 1/2 under binomial sex, conditioned on both sexes.
function effective_size()
    w = [Float64(binomial(big(N), k)) for k = 1:N-1]
    a = sum(w[k] * (1 / 4k + 1 / 4(N - k)) for k = 1:N-1) / sum(w)
    return 1 / a + 0.5
end
const NE = effective_size()

# P(i | c): i derived genes in a sample of NS from 2N genes with c derived.
function hypergeometric()
    P = zeros(2N - 1, NS - 1)
    for c = 1:2N-1, i = 1:NS-1
        (i <= c && NS - i <= 2N - c) || continue
        p = Float64(binomial(NS, i))
        for j = 0:i-1
            p *= (c - j) / (2N - j)
        end
        for j = 0:NS-i-1
            p *= (2N - c - j) / (2N - i - j)
        end
        P[c, i] = p
    end
    return P
end
const HYP = hypergeometric()

# ξ_i/θ of the Poisson random field, by Simpson's rule on [0, 1]; the
# integrand C(n, i) x^(i-1) (1 - x)^(n-i-1) g(x) is a smooth function.
function prf(α, i; m = 20_000)
    g(x) = α == 0 ? 1 - x : expm1(-α * (1 - x)) / expm1(-α)
    F(x) = binomial(NS, i) * x^(i - 1) * (1 - x)^(NS - i - 1) * g(x)
    h = 1 / m
    s = F(0.0) + F(1.0)
    for k = 1:m-1
        s += (isodd(k) ? 4 : 2) * F(k * h)
    end
    return s * h / 3
end

# One run: ξ_i/θ per available site for every class, and the fixed fraction.
function run(loci, cls, seed)
    res = fisher_wright(N, T, fill(CHRLEN, NCHR), 0.02; M = 1e8, mut_base = 1e12,
        fitness = SiteFitness(loci), result = true, rng = Xoshiro(seed))
    c = zeros(Int, K)
    for hap in res.selected_haplotypes, x in hap
        c[searchsortedfirst(loci.pos, x)] += 1
    end
    nc = length(ALPHA)
    fixed = zeros(Int, nc)
    for x in res.selected_substitutions
        fixed[cls[searchsortedfirst(loci.pos, x)]] += 1
    end
    xi = zeros(nc, NS - 1)
    for k = 1:K
        0 < c[k] < 2N || continue
        @views xi[cls[k], :] .+= HYP[c[k], :]
    end
    avail = K ÷ nc .- fixed
    return xi ./ (avail .* 4NE * MU), fixed ./ (K ÷ nc)
end

function main(runs = 1000, seed = 2026)
    println("FisherWright ", pkgversion(FisherWright), ", Julia ", VERSION, ", threads ", Threads.nthreads())
    @printf("N = %d (Ne = %.2f), %d generations, %d sites on %d × %d Mb (1 Morgan each), 4Nμ = %g per site\n",
        N, NE, T, K, NCHR, CHRLEN ÷ 1_000_000, THETA)
    println("h = 1/2, α = 2Ne·s ∈ $(ALPHA); sample of $NS genes; $runs runs, seed $seed")
    println("ξ_i/θ per available site, summed over the bins of i; reference: Poisson random field")
    println()
    pos = collect(range(CHRLEN ÷ (2K ÷ NCHR), step = NCHR * CHRLEN ÷ K, length = K))
    cls = [mod1(k, length(ALPHA)) for k = 1:K]
    loci = FitnessLoci(pos, ALPHA[cls] ./ 2NE; h = 0.5, μ = MU)
    out = Vector{Any}(undef, runs)
    Threads.@threads for r = 1:runs
        out[r] = run(loci, cls, hash((seed, r)))
    end
    @printf("%-8s %-6s %22s %10s %8s %7s\n", "α", "i", "FisherWright", "PRF", "ratio", "z")
    failed = String[]
    for (j, α) in enumerate(ALPHA)
        ref = [prf(α, i) for i = 1:NS-1]
        for b in [BINS; [1:NS-1]]
            x = [sum(o[1][j, b]) for o in out]
            e = sum(ref[b])
            m, se = mean(x), std(x) / sqrt(runs)
            z = (m - e) / se
            lab = first(b) == last(b) ? "$(first(b))" : b == 1:NS-1 ? "all" : "$(first(b))-$(last(b))"
            @printf("%-8g %-6s %22s %10.4f %8.4f %7.2f\n", α, lab, @sprintf("%.4f ± %.4f", m, se), e, m / e, z)
            abs(z) > ZMAX && push!(failed, "α = $α, i = $lab")
        end
        @printf("%-8g fixed in %d generations: %.4f%% of the sites\n", α, T, 100 * mean(o[2][j] for o in out))
        flush(stdout)
    end
    isempty(failed) || (println("\nFAILED: ", join(failed, "; ")); exit(1))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(parse.(Int, ARGS)...)
end
