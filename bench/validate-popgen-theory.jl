# Population genetics validation of FisherWright.jl against theory and msprime.
#
#   julia -t 8 --project=. bench/validate-popgen-theory.jl [reps] [ne] [seed]
#
# 1. Crossovers per meiosis: Poisson counts per chromosome, 0.5 assortment.
# 2. Equilibrium diversity: S vs. Watterson's θ·aₙ, Σ2pq vs. θ.
# 3. Site frequency spectrum: whole population vs. msprime DTWF, and a random
#    sample of `nsub` ≪ 2N haplotypes vs. the coalescent θ/k.
# 4. LD decay: pooled σ²_d = ΣD²/Σp₁q₁p₂q₂ per distance bin vs. msprime DTWF.
#    Ohta-Kimura and Sved curves are printed for reference only.
#
# FisherWright runs `ne` diploids for 20ne generations (10×2N, equilibrium) on
# 10 chromosomes of 10 Mb, μ = r = 1e-8, and `reps` replicates; msprime DTWF
# (via `uv`, bench/validate-popgen-msprime.py) simulates the same model. Each
# comparison is a z-score, and a check fails when |z| > ZMAX. The script exits
# with status 1 if any check fails.

using FisherWright: uniform_recombination_map, cobp!, fisher_wright, muts2bitarray
using Random, Statistics, Printf

const ZMAX = 3.5
const NCHR, CHRLEN, MU = 10, 10_000_000, 1e-8
const NSUB, MAF = 20, 0.1
const WHOLE = [(1, 1), (2, 2), (3, 5), (6, 20), (21, typemax(Int))]
const SUB = [1:1, 2:2, 3:4, 5:9, 10:19]
const BINS = [(10_000, 30_000), (30_000, 100_000), (100_000, 300_000), (300_000, 1_000_000), (1_000_000, 3_000_000)]

const failures = String[]
const section = Ref("")

function check(name, z)
    ok = abs(z) <= ZMAX
    ok || push!(failures, section[] * ": " * strip(name))
    return ok ? "pass" : "FAIL"
end

se(x) = std(x) / sqrt(length(x))

function ld_sums(bm, chr, pos)
    n = size(bm, 2)
    p = vec(sum(bm; dims = 2)) ./ n
    num, den = zeros(length(BINS)), zeros(length(BINS))
    for c in unique(chr)
        k = findall(i -> chr[i] == c && MAF <= p[i] <= 1 - MAF, eachindex(p))
        G = Float64.(bm[k, :])
        q = p[k]
        D = G * G' ./ n .- q * q'
        V = (q .* (1 .- q)) * (q .* (1 .- q))'
        d = Int.(pos[k])' .- Int.(pos[k])
        for (j, (lo, hi)) in enumerate(BINS)
            m = (lo .<= d) .& (d .< hi)
            num[j] += sum(abs2, D[m])
            den[j] += sum(V[m])
        end
    end
    return num, den
end

# Per-replicate statistics in the same layout as validate-popgen-msprime.py.
function fw_replicate(ne)
    res = fisher_wright(ne, 20ne, fill(CHRLEN, NCHR), MU * 1e8; result = true)
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

function msprime_replicates(ne, reps, seed)
    script = joinpath(@__DIR__, "validate-popgen-msprime.py")
    cmd = `uv run --with msprime --with numpy python $script --ne $ne --nchr $NCHR
           --chr-len $CHRLEN --mu $MU --r $MU --reps $reps --nsub $NSUB --maf $MAF --seed $seed`
    rows = [parse.(Float64, split(l)) for l in eachline(cmd) if !isempty(strip(l))]
    return permutedims(reduce(hcat, rows))
end

function validate_crossovers(; K = 100_000)
    section[] = "crossovers"
    println("\n[1/4] Crossovers per meiosis (K = $K, chromosomes of 50 and 100 Mb, M = 1e8)")
    rmap = uniform_recombination_map([50_000_000, 100_000_000]; M = 1e8)
    buf = UInt32[]
    c1, c2, as = zeros(Int, K), zeros(Int, K), zeros(Int, K)
    for k = 1:K
        cobp!(buf, rmap)
        c1[k] = count(<(50_000_000), buf)
        c2[k] = count(x -> 50_000_000 < x < 150_000_000, buf)
        as[k] = count(==(50_000_000), buf)
    end
    @printf("  %-28s %9s %9s %8s  %s\n", "statistic", "observed", "expected", "z", "")
    for (name, x, λ) in (("chr 1 mean", c1, 0.5), ("chr 2 mean", c2, 1.0))
        z = (mean(x) - λ) / sqrt(λ / K)
        @printf("  %-28s %9.4f %9.4f %8.2f  %s\n", name, mean(x), λ, z, check(name, z))
        d = var(x) / mean(x)
        z = (d - 1) / sqrt(2 / (K - 1))
        @printf("  %-28s %9.4f %9.4f %8.2f  %s\n", replace(name, "mean" => "variance/mean"), d, 1.0, z,
            check(replace(name, "mean" => "variance/mean"), z))
    end
    z = (mean(as) - 0.5) / sqrt(0.25 / K)
    @printf("  %-28s %9.4f %9.4f %8.2f  %s\n", "assortment P(switch)", mean(as), 0.5, z, check("assortment", z))
end

function main(reps = 50, ne = 100, seed = 2026)
    println("="^78)
    println(" FisherWright.jl population genetics validation")
    println(" ne = $ne, nt = $(20ne), $NCHR × $(CHRLEN ÷ 10^6) Mb, μ = r = $MU, reps = $reps, threads = $(Threads.nthreads())")
    println(" A check fails when |z| > $ZMAX")
    println("="^78)
    Random.seed!(seed)
    validate_crossovers()

    println("\nSimulating $reps FisherWright replicates ...")
    t = @elapsed fw = permutedims(reduce(hcat, [fw_replicate(ne) for _ = 1:reps]))
    @printf("  done in %.1f s\n", t)
    println("Simulating $reps msprime DTWF replicates ...")
    t = @elapsed ms = msprime_replicates(ne, reps, seed)
    @printf("  done in %.1f s\n", t)

    n = 2ne
    θ = 4ne * MU * NCHR * CHRLEN
    aₙ = sum(1 / i for i = 1:n-1)
    col(M, j) = M[:, j]
    zdiff(a, b) = (mean(a) - mean(b)) / sqrt(se(a)^2 + se(b)^2)
    zth(a, e) = (mean(a) - e) / se(a)
    row(name, a, b, e, z) = @printf("  %-22s %9.1f ± %-6.1f %9.1f ± %-6.1f %9.1f %8.2f  %s\n",
        name, mean(a), se(a), mean(b), se(b), e, z, check(name, z))
    hdr() = @printf("  %-22s %18s %18s %9s %8s\n", "statistic", "FisherWright", "msprime DTWF", "theory", "z")

    # Σ2pq over the whole population uses population frequencies, so its
    # expectation is θ(1 - 1/n) rather than θ.
    section[] = "diversity"
    println("\n[2/4] Equilibrium diversity (z: FisherWright vs. theory, then vs. msprime DTWF)")
    hdr()
    row("S (θ·aₙ)", col(fw, 1), col(ms, 1), θ * aₙ, zth(col(fw, 1), θ * aₙ))
    row("S vs. msprime", col(fw, 1), col(ms, 1), θ * aₙ, zdiff(col(fw, 1), col(ms, 1)))
    row("Σ2pq (θ(1-1/n))", col(fw, 2), col(ms, 2), θ * (1 - 1 / n), zth(col(fw, 2), θ * (1 - 1 / n)))
    row("Σ2pq vs. msprime", col(fw, 2), col(ms, 2), θ * (1 - 1 / n), zdiff(col(fw, 2), col(ms, 2)))
    @printf("  %-22s %9.3f %28.3f\n", "S / θ·aₙ", mean(col(fw, 1)) / (θ * aₙ), mean(col(ms, 1)) / (θ * aₙ))

    println("\n[3/4] Site frequency spectrum")
    section[] = "SFS whole population"
    println("  Whole population ($n haplotypes), z: FisherWright vs. msprime DTWF")
    hdr()
    for (j, (lo, hi)) in enumerate(WHOLE)
        hi = min(hi, n - 1)
        e = θ * sum(1 / k for k = lo:hi)
        row(lo == hi ? "k = $lo" : "k = $lo-$hi", col(fw, 2 + j), col(ms, 2 + j), e, zdiff(col(fw, 2 + j), col(ms, 2 + j)))
    end
    section[] = "SFS sample"
    println("  Random sample of $NSUB haplotypes, z: FisherWright vs. θ/k")
    hdr()
    for (j, g) in enumerate(SUB)
        e = θ * sum(1 / k for k in g)
        c = 2 + length(WHOLE) + j
        row(length(g) == 1 ? "k = $(first(g))" : "k = $(first(g))-$(last(g))", col(fw, c), col(ms, c), e, zth(col(fw, c), e))
    end

    section[] = "LD"
    println("\n[4/4] LD decay, pooled σ²_d (MAF ≥ $MAF), z: FisherWright vs. msprime DTWF")
    @printf("  %-14s %18s %18s %9s %8s %8s\n", "distance", "FisherWright", "msprime DTWF", "z", "Ohta-K.", "Sved")
    o = 2 + length(WHOLE) + length(SUB)
    nb = length(BINS)
    for (j, (lo, hi)) in enumerate(BINS)
        r(M) = M[:, o+j] ./ M[:, o+nb+j]
        pooled(M) = sum(M[:, o+j]) / sum(M[:, o+nb+j])
        a, b = r(fw), r(ms)
        z = (pooled(fw) - pooled(ms)) / sqrt(se(a)^2 + se(b)^2)
        ρ = 4ne * sqrt(lo * hi) * 1e-8
        name = @sprintf("%g-%g Mb", lo / 1e6, hi / 1e6)
        @printf("  %-14s %9.4f ± %-6.4f %9.4f ± %-6.4f %9.2f %8.4f %8.4f  %s\n", name, pooled(fw), se(a),
            pooled(ms), se(b), z, (10 + ρ) / (22 + 13ρ + ρ^2), 1 / (1 + ρ), check(name, z))
    end

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
    exit(main(args...) ? 0 : 1)
end
