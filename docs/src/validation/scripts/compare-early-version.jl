# Compare an early FisherWright release with the current one, at t = N and at
# equilibrium, on population statistics and on the spectrum of a founder sample.
#
# Run through compare-early-version.sh, which runs this script once per
# (version, t) arm. With a FisherWright environment active:
#
#   julia -t 8 compare-early-version.jl <tmult> [reps=30] [seed=2026]
#
# Model: N = 200 diploids, 10 chromosomes of 20 Mb, μ = 1e-8 per bp per
# meiosis, 1 cM/Mb, tmult·N generations from no variation.
#
# Population statistics (all 2N haplotypes):
#   S, Σ2pq, number of sites with MAF ≥ 0.2, and mean r² among MAF ≥ 0.2 sites
#   on the same chromosome in three distance bins.
# Founder statistics (NF individuals sampled from the final generation; QTL
# candidates are the sites segregating in the sample, as in a breeding
# simulation that samples QTL from its founders):
#   number of segregating sites, fraction with MAF ≥ 0.2, fraction with
#   MAF < 0.05, mean 2pq per site, and share of Σ2pq from MAF < 0.05 sites.
#
# Versions before v0.3.0 return `(haplotypes, cbp)` and can leave haplotypes
# unsorted or with repeated positions (see the v0.3.0 notes in the README).
# Positions are counted once per haplotype, and the number of affected
# haplotypes is reported.

using FisherWright, Random, Statistics, Printf

const tmult = parse(Int, ARGS[1])
const reps = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 30
const seed = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 2026
const N, NF, NCHR, CLEN = 200, 75, 10, 20_000_000
const BINS = [(0, 500_000), (500_000, 2_000_000), (2_000_000, 5_000_000)]
const old_api = pkgversion(FisherWright) < v"0.3.0"

# Early releases print progress to stdout.
simulate() = old_api ?
    first(redirect_stdout(() -> fisher_wright(N, tmult * N, fill(CLEN, NCHR), 1.0), devnull)) :
    fisher_wright(N, tmult * N, fill(CLEN, NCHR), 1.0; result = true).active_haplotypes

function counts(haps)
    cnt = Dict{UInt32,Int}()
    for h in haps, p in unique(h)
        cnt[p] = get(cnt, p, 0) + 1
    end
    cnt
end

function population_stats(haps)
    n = length(haps)
    S, H, com = 0, 0.0, UInt32[]
    for (p, c) in counts(haps)
        0 < c < n || continue
        q = c / n
        S += 1
        H += 2q * (1 - q)
        0.2 ≤ q ≤ 0.8 && push!(com, p)
    end
    sort!(com)
    idx = Dict(p => i for (i, p) in enumerate(com))
    X = zeros(length(com), n)
    for (j, h) in enumerate(haps), p in h
        i = get(idx, p, 0)
        i > 0 && (X[i, j] = 1.0)
    end
    pf = vec(mean(X, dims = 2))
    r2 = [Float64[] for _ in BINS]
    for i in eachindex(com), k = i+1:length(com)
        d = Int(com[k]) - Int(com[i])
        d ≥ BINS[end][2] && break
        div(Int(com[i]) - 1, CLEN) == div(Int(com[k]) - 1, CLEN) || continue
        D = mean(view(X, i, :) .* view(X, k, :)) - pf[i] * pf[k]
        v = D^2 / (pf[i] * (1 - pf[i]) * pf[k] * (1 - pf[k]))
        for (b, (lo, hi)) in enumerate(BINS)
            lo ≤ d < hi && push!(r2[b], v)
        end
    end
    (S = S, H = H, ncom = length(com), r2 = mean.(r2))
end

function founder_stats(haps)
    ids = randperm(N)[1:NF]
    fh = [haps[h] for i in ids for h in (2i - 1, 2i)]
    n = length(fh)
    q = [c / n for c in values(counts(fh)) if 0 < c < n]
    maf = min.(q, 1 .- q)
    h2 = 2 .* q .* (1 .- q)
    (seg = length(q), fcom = mean(maf .≥ 0.2), frare = mean(maf .< 0.05),
     h2 = mean(h2), rareshare = sum(h2[maf .< 0.05]) / sum(h2))
end

bad(haps) = (count(!issorted, haps), count(!allunique, haps), sum(h -> length(h) - length(unique(h)), haps))

Random.seed!(seed)
res = map(1:reps) do _
    haps = simulate()
    (population_stats(haps), founder_stats(haps), bad(haps))
end

ms(v) = (v = filter(!isnan, v); (mean(v), std(v) / sqrt(length(v))))
pop, fnd = first.(res), getindex.(res, 2)
@printf("FisherWright v%s, N = %d, t = %dN, %d × %d Mb, μ = 1e-8, 1 cM/Mb, reps = %d, seed = %d, threads = %d\n",
    pkgversion(FisherWright), N, tmult, NCHR, CLEN ÷ 1_000_000, reps, seed, Threads.nthreads())
@printf("  final generation, per replicate: %.1f haplotypes unsorted, %.1f with repeated positions (%.1f repeats in total)\n",
    (mean(getindex.(last.(res), i)) for i = 1:3)...)
println("  population (2N = $(2N) haplotypes)")
@printf("    S                     %10.1f ± %6.1f\n", ms([r.S for r in pop])...)
@printf("    Σ2pq                  %10.1f ± %6.1f\n", ms([r.H for r in pop])...)
@printf("    sites with MAF ≥ 0.2  %10.1f ± %6.1f\n", ms([r.ncom for r in pop])...)
for (b, (lo, hi)) in enumerate(BINS)
    @printf("    r², MAF ≥ 0.2, %3.1f–%3.1f Mb %8.4f ± %6.4f\n", lo / 1e6, hi / 1e6, ms([r.r2[b] for r in pop])...)
end
println("  founder sample ($NF individuals, $(2NF) haplotypes; sites segregating in the sample)")
@printf("    segregating sites     %10.1f ± %6.1f\n", ms([r.seg for r in fnd])...)
@printf("    fraction MAF ≥ 0.2    %10.4f ± %6.4f\n", ms([r.fcom for r in fnd])...)
@printf("    fraction MAF < 0.05   %10.4f ± %6.4f\n", ms([r.frare for r in fnd])...)
@printf("    mean 2pq per site     %10.4f ± %6.4f\n", ms([r.h2 for r in fnd])...)
@printf("    Σ2pq share, MAF < 0.05 %9.4f ± %6.4f\n", ms([r.rareshare for r in fnd])...)
println()
