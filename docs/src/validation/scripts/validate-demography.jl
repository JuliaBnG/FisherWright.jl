# Validation of FisherWright.jl under a changing population size, against the
# exact two-sex recursion and against msprime's discrete-time Wright-Fisher
# model (DTWF) with the same demography.
#
# With FisherWright installed, run from this directory:
#
#   julia -t 8 validate-demography.jl [reps] [seed]
#
# Model. N = 100 diploids, two sexes (sex Binomial(N, 1/2) per generation,
# conditioned on both sexes), random mating, 10 chromosomes of 10 Mb,
# μ = r = 1e-8, free recombination between chromosomes. Each replicate runs
# 20N = 2000 generations from a population without variation (10 × 2N,
# mutation-drift equilibrium) and is then continued, through `init`, in two
# scenarios, each sampled as the whole final population:
#
#   bottleneck: 20 individuals for 20 generations, then 100 for 20 generations
#               (200 haplotypes);
#   expansion:  400 individuals for 100 generations (800 haplotypes).
#
# The two scenarios share their burn-in, so their results are correlated with
# each other; every check is within one scenario.
#
# References.
# 1. Σ2pq against its exact expectation. With a = E[1/(4N_m) + 1/(4N_f)] over
#    the sex counts of the N parents of a generation and U = μL new mutations
#    per haplotype (infinite sites),
#      H_B' = (1 - a)·H_B + (a/2)·H_W + 2U,   H_W' = H_B + 2U,
#    for two genes in different individuals and in one individual, and
#    Σ2pq = H_W/(2N) + (1 - 1/N)·H_B over the N individuals of a cohort. The
#    recursion is run from H = 0 through the same 2000 + nt generations, so
#    it assumes neither equilibrium nor an effective size.
# 2. S, Σ2pq, the SFS of the whole population and of a random sample of 20
#    haplotypes, and σ²_d per distance bin against msprime DTWF
#    (validate-popgen-msprime.py --epochs, via `uv`). msprime's DTWF picks
#    both parents of an individual from one pool of N' individuals, so a
#    cohort of N two-sex parents is represented by N' = Ne(N), the effective
#    size from the leading eigenvalue of the recursion above (`dioecious_ne`:
#    19.38, 99.48 and 399.50 for N = 20, 100 and 400). N' must be an integer:
#    each msprime chromosome draws floor(Ne) or ceil(Ne) for every epoch with
#    probabilities that make the mean Ne. The statistics are smooth in N', so
#    the mean over replicates is the value at Ne up to terms of order
#    (1/N')², below 0.1% here.
#    Generation g of FisherWright draws the parents of its N_g offspring from
#    the N_{g-1} individuals of generation g - 1; msprime draws the parents
#    of the individuals at time t - 1 from the size set for time t. Cohort g
#    is at time nt - g, so an epoch that starts with cohort g₀ in
#    FisherWright ends, backwards, at msprime time nt - g₀, and the earlier
#    size applies from time nt - g₀ + 1.
#
# Each comparison is a z-score, and a check fails when |z| > ZMAX. The script
# exits with status 1 if any check fails.

include(joinpath(@__DIR__, "validate-popgen-theory.jl"))
using FisherWright: bottleneck

const N0, NB, NX = 100, 20, 400
const TB, TR, TX = 20, 20, 100          # bottleneck, recovery, expansion (generations)

# Statistics of one population in the layout of validate-popgen-msprime.py.
function popstats(res)
    haps = res.active_haplotypes
    n = length(haps)
    cnt = Dict{UInt32,Int}()
    for h in haps, p in h
        cnt[p] = get(cnt, p, 0) + 1
    end
    k = [c for c in values(cnt) if 0 < c < n]
    f = k ./ n
    w = [count(c -> lo <= c <= hi, k) for (lo, hi) in WHOLE]
    sub = Dict{UInt32,Int}()
    for h in haps[randperm(n)[1:NSUB]], p in h
        sub[p] = get(sub, p, 0) + 1
    end
    ks = [c for c in values(sub) if c < NSUB]
    s = [count(in(g), ks) for g in SUB]
    bm, lmp = muts2bitarray(haps, res.chromosome_ends)
    num, den = ld_sums(bm, lmp.chr, lmp.pos)
    return vcat(length(k), sum(2 .* f .* (1 .- f)), w, s, num, den)
end

# One burn-in, continued in both scenarios.
function fw_replicate()
    chr = fill(CHRLEN, NCHR)
    burn = fisher_wright(N0, 20N0, chr, MU * 1e8; result = true)
    bot = fisher_wright(N0, TB + TR, chr, MU * 1e8; init = burn, result = true,
        popsize = bottleneck(N0, NB, 1, TB))
    exp = fisher_wright(N0, TX, chr, MU * 1e8; init = burn, result = true, popsize = g -> NX)
    @assert bot.census[2001:end] == [N0; fill(NB, TB); fill(N0, TR)]
    @assert exp.census[2001:end] == [N0; fill(NX, TX)]
    return popstats(bot), popstats(exp)
end

# E[1/(4N_m) + 1/(4N_f)] over the sex counts of N parents.
function coalescence(N)
    w = [Float64(binomial(big(N), k) / big(2)^N) for k = 1:N-1]
    return sum(w[k] * (1 / 4k + 1 / 4(N - k)) for k = 1:N-1) / sum(w)
end

# Exact E[Σ2pq] of the last cohort; sizes[g] is the number of offspring of
# generation g, the first parents are n0 individuals without variation, and
# U is the expected number of new mutations per haplotype per generation.
function expected_pi(n0, sizes, U)
    a = Dict(N => coalescence(N) for N in unique([n0; sizes]))
    hB = hW = 0.0
    N = n0
    for n in sizes
        hB, hW = (1 - a[N]) * hB + a[N] / 2 * hW + 2U, hB + 2U
        N = n
    end
    return hW / 2N + (1 - 1 / N) * hB
end

function msprime_scenario(epochs, samples, reps, seed)
    script = joinpath(@__DIR__, "validate-popgen-msprime.py")
    ep = ["$t:$x" for (t, x) in epochs]
    cmd = `uv run --with msprime --with numpy python $script --samples $samples --nchr $NCHR
           --chr-len $CHRLEN --mu $MU --r $MU --reps $reps --nsub $NSUB --maf $MAF --seed $seed --epochs $ep`
    rows = [parse.(Float64, split(l)) for l in eachline(cmd) if !isempty(strip(l))]
    return permutedims(reduce(hcat, rows))
end

function report(name, fw, ms, n, Eπ)
    col(M, j) = M[:, j]
    zdiff(a, b) = (mean(a) - mean(b)) / sqrt(se(a)^2 + se(b)^2)
    row(label, a, b, z) = @printf("  %-22s %9.1f ± %-6.1f %9.1f ± %-6.1f %8.2f  %s\n",
        label, mean(a), se(a), mean(b), se(b), z, check(label, z))
    hdr() = @printf("  %-22s %18s %18s %8s\n", "statistic", "FisherWright", "msprime DTWF", "z")

    section[] = "$name diversity"
    println("\n  Diversity")
    π = col(fw, 2)
    z = (mean(π) - Eπ) / se(π)
    @printf("  %-22s %9.1f ± %-6.1f %18.1f %8.2f  %s\n", "Σ2pq vs. recursion", mean(π), se(π), Eπ, z,
        check("Σ2pq vs. recursion", z))
    hdr()
    row("Σ2pq vs. msprime", π, col(ms, 2), zdiff(π, col(ms, 2)))
    row("S vs. msprime", col(fw, 1), col(ms, 1), zdiff(col(fw, 1), col(ms, 1)))

    section[] = "$name SFS whole population"
    println("  SFS, whole population ($n haplotypes)")
    hdr()
    for (j, (lo, hi)) in enumerate(WHOLE)
        hi = min(hi, n - 1)
        row(lo == hi ? "k = $lo" : "k = $lo-$hi", col(fw, 2 + j), col(ms, 2 + j), zdiff(col(fw, 2 + j), col(ms, 2 + j)))
    end
    section[] = "$name SFS sample"
    println("  SFS, random sample of $NSUB haplotypes")
    hdr()
    for (j, g) in enumerate(SUB)
        c = 2 + length(WHOLE) + j
        row(length(g) == 1 ? "k = $(first(g))" : "k = $(first(g))-$(last(g))", col(fw, c), col(ms, c), zdiff(col(fw, c), col(ms, c)))
    end

    section[] = "$name LD"
    println("  LD decay, pooled σ²_d (MAF ≥ $MAF)")
    @printf("  %-22s %18s %18s %8s\n", "distance", "FisherWright", "msprime DTWF", "z")
    o = 2 + length(WHOLE) + length(SUB)
    nb = length(BINS)
    for (j, (lo, hi)) in enumerate(BINS)
        r(M) = M[:, o+j] ./ M[:, o+nb+j]
        pooled(M) = sum(M[:, o+j]) / sum(M[:, o+nb+j])
        a, b = r(fw), r(ms)
        z = (pooled(fw) - pooled(ms)) / sqrt(se(a)^2 + se(b)^2)
        label = @sprintf("%g-%g Mb", lo / 1e6, hi / 1e6)
        @printf("  %-22s %9.4f ± %-6.4f %9.4f ± %-6.4f %8.2f  %s\n", label, pooled(fw), se(a), pooled(ms), se(b),
            z, check(label, z))
    end
end

function demography_main(reps = 100, seed = 2026)
    println("="^78)
    println(" FisherWright.jl validation under a changing population size")
    println(" N = $N0 for $(20N0) generations, then")
    println("   bottleneck: $NB for $TB generations, $N0 for $TR generations")
    println("   expansion:  $NX for $TX generations")
    println(" $NCHR × $(CHRLEN ÷ 10^6) Mb, μ = r = $MU, reps = $reps, threads = $(Threads.nthreads())")
    println(" A check fails when |z| > $ZMAX")
    println("="^78)
    Random.seed!(seed)

    println("\nSimulating $reps FisherWright replicates (one burn-in, two scenarios) ...")
    t = @elapsed runs = [fw_replicate() for _ = 1:reps]
    @printf("  done in %.1f s\n", t)
    fwb = permutedims(reduce(hcat, first.(runs)))
    fwx = permutedims(reduce(hcat, last.(runs)))

    ne = Dict(N => dioecious_ne(N) for N in (NB, N0, NX))
    @printf("\nEffective size of a two-sex cohort (msprime N'): %s\n",
        join([@sprintf("Ne(%d) = %.2f", N, ne[N]) for N in (NB, N0, NX)], ", "))
    # Backwards in time: the size at time t is that of the parents of the
    # cohort at time t - 1.
    epb = [(0, ne[N0]), (TR, ne[NB]), (TR + TB, ne[N0])]
    epx = [(0, ne[NX]), (TX, ne[N0])]
    println("Simulating $reps msprime DTWF replicates per scenario ...")
    t = @elapsed begin
        msb = msprime_scenario(epb, N0, reps, seed)
        msx = msprime_scenario(epx, NX, reps, seed + 1)
    end
    @printf("  done in %.1f s\n", t)

    U = MU * NCHR * CHRLEN
    burn = fill(N0, 20N0)
    println("\n[1/2] Bottleneck: $N0 → $NB ($TB generations) → $N0 ($TR generations)")
    report("bottleneck", fwb, msb, 2N0, expected_pi(N0, [burn; fill(NB, TB); fill(N0, TR)], U))
    println("\n[2/2] Expansion: $N0 → $NX ($TX generations)")
    report("expansion", fwx, msx, 2NX, expected_pi(N0, [burn; fill(NX, TX)], U))
    # What the Σ2pq check can resolve: one generation of the bottleneck or
    # of the expansion more, in units of the standard error of the mean.
    println("\n  Resolution of Σ2pq vs. recursion (change in E[Σ2pq] / SE)")
    eb(t) = expected_pi(N0, [burn; fill(NB, t); fill(N0, TR)], U)
    ex(t) = expected_pi(N0, [burn; fill(NX, t)], U)
    @printf("  %-44s %7.1f  (%5.1f SE)\n", "bottleneck of $(TB + 1) instead of $TB generations",
        eb(TB + 1) - eb(TB), (eb(TB + 1) - eb(TB)) / se(fwb[:, 2]))
    @printf("  %-44s %7.1f  (%5.1f SE)\n", "bottleneck of $(NB + 1) instead of $NB individuals",
        expected_pi(N0, [burn; fill(NB + 1, TB); fill(N0, TR)], U) - eb(TB),
        (expected_pi(N0, [burn; fill(NB + 1, TB); fill(N0, TR)], U) - eb(TB)) / se(fwb[:, 2]))
    @printf("  %-44s %7.1f  (%5.1f SE)\n", "expansion of $(TX + 10) instead of $TX generations",
        ex(TX + 10) - ex(TX), (ex(TX + 10) - ex(TX)) / se(fwx[:, 2]))
    @printf("  %-44s %7.1f\n", "E[Σ2pq] at the end of the burn-in", expected_pi(N0, burn, U))

    println("\n", "="^78)
    if isempty(failures)
        println(" All checks passed.")
    else
        println(" FAILED: ", join(failures, "; "))
    end
    println("="^78)
    return isempty(failures)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse.(Int, ARGS)
    exit(demography_main(args...) ? 0 : 1)
end
