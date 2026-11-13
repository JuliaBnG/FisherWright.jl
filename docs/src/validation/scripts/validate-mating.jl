# Mating schemes of FisherWright against a matched SLiM model
# (validate-mating.slim; docs/dev/roadmap.md, validation 4).
#
# Both simulators start from the same kind of polymorphic population (each
# of L = 1,000 sites on 20 chromosomes of 50 bp and 1 Morgan carried by each
# haplotype with probability 1/2), run NT = 30 generations of N = 40
# individuals without mutation under one mating scheme, and report Σ2pq/L,
# the observed heterozygosity H_obs and the number of offspring of
# relatives that could not be avoided. The schemes are implemented
# independently in Eidos, on SLiM's meiosis and pedigree. For each scheme
# and statistic the table gives the mean ± SE over R replicates of each
# simulator, the z-score of their difference, and, where it is known, the
# exact expectation (recursion of the heterozygosity within and between
# individuals) with the z-score of FisherWright against it.
#
# Usage: julia --project -t 8 validate-mating.jl SLIMDIR R SEED [SCHEMES]
# where SLIMDIR holds the SLiM output, one file <scheme>.tsv per scheme
# (columns scheme, replicate, Σ2pq/L, H_obs, forced), and SCHEMES is an
# optional space-separated list of scheme keys (all by default).

using FisherWright, Random, Statistics, Printf

const N, NT, L = 40, 30, 1_000
const CHR = fill(50, 20)

const SCHEMES = [
    ("rm", "RandomMating()", RandomMating()),
    ("rm5", "RandomMating(nmale=5)", RandomMating(nmale=5)),
    ("fsgamma", "FamilySize(dams=1.0, sires=:random, nmale=10)", FamilySize(dams=1.0, sires=:random, nmale=10)),
    ("hier", "Hierarchical(4, 3)", Hierarchical(4, 3)),
    ("mono", "Monoecious()", Monoecious()),
    ("self", "Monoecious(selfing=0.5)", Monoecious(selfing=0.5)),
    ("avoid_half", "AvoidRelatives(RandomMating(); level=:halfsib)", AvoidRelatives(RandomMating(); level=:halfsib)),
    ("avoid_cousin5", "AvoidRelatives(RandomMating(nmale=5); level=:cousin)",
        AvoidRelatives(RandomMating(nmale=5); level=:cousin)),
    ("avoid_mono", "AvoidRelatives(Monoecious(selfing=0); level=:fullsib)",
        AvoidRelatives(Monoecious(selfing=0.0); level=:fullsib)),
]

# Σ2pq/L (with replacement, over the 2N haplotypes) and H_obs.
function stats(haps)
    c = zeros(Int, L)
    for h in haps, p in h
        c[p] += 1
    end
    n = length(haps)
    H = sum(x -> 2 * (x / n) * (1 - x / n), c) / L
    O = mean(length(symdiff(haps[2i-1], haps[2i])) for i = 1:n÷2) / L
    return H, O
end

function run_fw(m, R, rng)
    out = zeros(R, 3)
    for r = 1:R
        start = [UInt32.(findall(rand(rng, Bool, L))) for _ = 1:2N]
        res = fisher_wright(N, NT, CHR, 1.0; M=50, init=start, mutate=false, mating=m, result=true, rng=rng)
        out[r, 1], out[r, 2] = stats(FisherWright._with_fixed(res))
        out[r, 3] = sum(res.related_matings)
    end
    return out
end

function read_slim(path)
    rows = [split(l, '\t') for l in eachline(path) if count(==('\t'), l) == 4]
    return [parse(Float64, r[k]) for r in rows, k = 3:5]
end

# Exact expectation of (Σ2pq/L, H_obs) after NT generations, from the
# probability a that two genes of different offspring come from one parent
# and the proportion σ of selfed offspring (both constant).
function expected(a, σ)
    hB = hW = 0.5
    for _ = 1:NT
        hB, hW = (1 - a) * hB + a / 2 * hW, (1 - σ) * hB + σ * hW / 2
    end
    return hW / 2N + (1 - 1 / N) * hB, hW
end

# Probability that two of n items split among m parents share one.
share(n, m, split) = split === :random ? 1 / m :
                     split === :equal ? (r = n % m; q = n ÷ m; (r * (q + 1) * q + (m - r) * q * (q - 1)) / (n * (n - 1))) :
                     (split + 1) / (m * split + 1)     # Gamma(split)
famA(n, Nm, Nf, ds, ss) = (Pf = share(n, Nf, ds); (Pf + Pf + (1 - Pf) * share(Nf, Nm, ss)) / 4)

function theory(key)
    key == "rm" || key == "rm5" || key == "fsgamma" || key == "hier" || key == "mono" || key == "self" || return nothing
    key == "mono" && return expected(1 / N, 1 / N)
    key == "self" && return expected(1 / N, 0.5)
    key == "rm5" && return expected(1 / 20 + 1 / 140, 0.0)
    key == "fsgamma" && return expected(famA(N, 10, 30, 1.0, :random), 0.0)
    key == "hier" && return expected(famA(N, 4, 12, :equal, :equal), 0.0)
    # Binomial sex, conditioned on both sexes: E[a] over the number of males.
    p = [binomial(big(N), k) for k = 1:N-1]
    a = sum(Float64(p[k]) * (1 / 4k + 1 / 4(N - k)) for k = 1:N-1) / Float64(sum(p))
    return expected(a, 0.0)
end

function main(slimdir, R, seed, only=nothing)
    println("FisherWright ", pkgversion(FisherWright), ", Julia ", VERSION, ", threads ", Threads.nthreads())
    println("N = $N, $NT generations, L = $L sites on 20 chromosomes of 1 Morgan, R = $R replicates per simulator, seed $seed")
    println()
    @printf("%-14s %-8s %20s %20s %7s %10s %7s\n", "scheme", "stat", "FisherWright", "SLiM", "z", "theory", "z(FW)")
    rng = Xoshiro(seed)
    for (key, label, m) in SCHEMES
        only === nothing || key in only || continue
        fw = run_fw(m, R, rng)
        sl = read_slim(joinpath(slimdir, key * ".tsv"))
        th = theory(key)
        println("# ", label, "  (SLiM replicates: ", size(sl, 1), ")")
        for (k, name) in enumerate(("Σ2pq/L", "H_obs", "forced"))
            a, b = fw[:, k], sl[:, k]
            ma, mb = mean(a), mean(b)
            sa, sb = std(a) / sqrt(length(a)), std(b) / sqrt(length(b))
            if sa == 0 && sb == 0
                @printf("%-14s %-8s %20s %20s\n", key, name, @sprintf("%.4g", ma), @sprintf("%.4g", mb))
                continue
            end
            z = (ma - mb) / sqrt(sa^2 + sb^2)
            t = th === nothing || k == 3 ? "" : @sprintf("%.5f", th[k])
            zt = th === nothing || k == 3 ? "" : @sprintf("%.2f", (ma - th[k]) / sa)
            @printf("%-14s %-8s %20s %20s %7.2f %10s %7s\n", key, name,
                @sprintf("%.5f ± %.5f", ma, sa), @sprintf("%.5f ± %.5f", mb, sb), z, t, zt)
        end
    end
end

main(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]), length(ARGS) > 3 ? split(ARGS[4]) : nothing)
