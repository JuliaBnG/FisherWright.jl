# Mutation–drift diversity under the mating schemes of FisherWright.jl
# (docs/dev/roadmap.md, validation 2: equilibrium π and S).
#
#   julia -t 16 --project=. docs/src/validation/scripts/validate-equilibrium.jl [reps = 1000] [seed = 2026] [breps = 10^6]
#
# reps forward replicates and breps backward genealogies per scheme.
#
# Model: N = 100 diploid parents, T = 2,000 generations (20N) from a
# monomorphic start, 10 chromosomes of 10 Mb (0.1 Morgan each), neutral
# mutation μ = 10^-8 per bp per haplotype per generation.
#
# Schemes: random mating with binomial sex (the default), 10 males and 90
# females, exactly half males; FamilySize with overdispersed dam families
# (Gamma(1) weights, 20 sires drawn at random per dam) and with equal
# families (20 sires, 4 dams each, 80 dams with one or two offspring);
# Monoecious with random selfing (1/N), and with selfing 0.2 and 0.5.
#
# Reference. For every scheme the expected heterozygosity per bp between
# two genes of different individuals (h_B) and within one individual (h_W)
# follows exactly, for infinite sites, from
#   h_B' = (1 - a) h_B + (a/2) h_W + 2μ,   h_W' = (1 - σ) h_B + (σ/2) h_W + 2μ,
# with a the probability that two genes of different offspring come from
# one parent and σ the selfing rate (0 with separate sexes). a is
# 1/(4Nm) + 1/(4Nf) under random mating (averaged over the binomial number
# of males, conditioned on both sexes), (P_f + P_m)/4 from the closed forms
# of the splits under FamilySize, and 1/N for hermaphrodites. Compared:
#   Σ2pq per bp over all 2N haplotypes, with expectation
#       h_W/(2N) + (1 - 1/N) h_B           (exact);
#   observed heterozygosity per bp, with expectation h_W   (exact);
#   S in a sample of one random gene from each of n = 20 individuals,
#       μ · L · E[Σ_t k_t 1{k_t ≥ 2}]       (exact, by simulation),
# where k_t is the number of distinct genes at generation t ancestral to
# the sample (each carries one transmission, and so μ mutations per bp,
# from generation t - 1; a mutation on the only lineage is fixed in the
# sample). E[·] comes from an independent backward simulation of the
# genealogy of the sample through the pedigree of each scheme, drawn
# generation by generation from the documented rules of the scheme (code
# below, not that of the package). The backward simulation is checked
# against the exact h_B with n = 2 genes. Kingman's approximation,
# h_B · L · Σ_{i<n} 1/i, is printed for comparison: with n = 20 genes and
# a up to 0.028 several coalescences are expected per generation, and the
# approximation falls short by up to 1.4% here.
# The sites of FisherWright are biallelic; at μT = 2·10^-5 the difference
# from infinite sites is far below the standard errors here.
#
# Each comparison is a z-score against the replicate standard error, with
# that of the backward simulation added for S; the script exits with
# status 1 if a comparison has |z| > 3.5.

using FisherWright, Random, Statistics, Printf

const N, T = 100, 2_000
const NCHR, CHRLEN, MU = 10, 10_000_000, 1e-8
const L = NCHR * CHRLEN
const NS = 20
const ZMAX = 3.5

# Probability that two distinct items of `total` fall to the same one of `n`
# parents under the splits of FamilySize: E[Σ c(c - 1)] / (total(total - 1)).
function same_parent(mode, k, total, n)
    total < 2 && return 0.0
    if mode === :equal
        q, r = divrem(total, n)
        return (r * (q + 1) * q + (n - r) * q * (q - 1)) / (total * (total - 1))
    end
    return mode === :random ? 1 / n : (k + 1) / (n * k + 1)
end

# a = (P_f + P_m)/4 of FamilySize for N parents with Nm males and N offspring.
function family_a(m::FamilySize, Nm)
    Nf = N - Nm
    Pf = same_parent(m.dams, m.kd, N, Nf)
    ps = same_parent(m.sires, m.km, Nf, Nm)
    return (Pf + (Pf + (1 - Pf) * ps)) / 4
end

# a under binomial sex, conditioned on both sexes being present.
function binomial_a()
    w = [Float64(binomial(big(N), k)) for k = 1:N-1]
    return sum(w[k] * (1 / 4k + 1 / 4(N - k)) for k = 1:N-1) / sum(w)
end

# Independent generators of the parents of one generation of N offspring
# from N parents: gene 1 of offspring i comes from parent pm[i, 1] (the
# sire), gene 2 from pm[i, 2] (the dam). Sexes are a random subset of the
# parents, redrawn every generation.
function sexes(rng, nm)
    p = randperm(rng, N)
    return p[1:nm], p[nm+1:N]
end

function binomial_males(rng)
    while true
        k = count(_ -> rand(rng, Bool), 1:N)
        0 < k < N && return k
    end
end

function random_pm!(pm, rng, nm)
    m, f = sexes(rng, nm == 0 ? binomial_males(rng) : nm)
    for i = 1:N
        pm[i, 1], pm[i, 2] = rand(rng, m), rand(rng, f)
    end
end

# `total` items among n parents: as equal as possible (the r that get one
# more at random), or multinomial with Exp(1) = Gamma(1) weights.
function split(rng, total, n, mode)
    c = zeros(Int, n)
    if mode === :equal
        q, r = divrem(total, n)
        c .= q
        c[randperm(rng, n)[1:r]] .+= 1
    else
        w = cumsum(randexp(rng, n))
        for _ = 1:total
            c[searchsortedfirst(w, rand(rng) * w[end])] += 1
        end
    end
    return c
end

# Offspring split among the dams, dams among the sires (in random order).
function family_pm!(pm, rng, nm, dams, sires)
    m, f = sexes(rng, nm)
    nd = split(rng, N, N - nm, dams)
    sod = sires === :random ? rand(rng, 1:nm, N - nm) :
          shuffle!(rng, reduce(vcat, fill.(1:nm, split(rng, N - nm, nm, :equal))))
    i = 0
    for d = 1:N-nm, _ = 1:nd[d]
        i += 1
        pm[i, 1], pm[i, 2] = m[sod[d]], f[d]
    end
end

# First parent uniform; the second uniform (s < 0), or the first with
# probability s and otherwise uniform among the others.
function mono_pm!(pm, rng, s)
    for i = 1:N
        p = rand(rng, 1:N)
        q = s < 0 ? rand(rng, 1:N) : rand(rng) < s ? p : mod1(p + rand(rng, 1:N-1), N)
        pm[i, 1], pm[i, 2] = p, q
    end
end

const SCHEMES = [
    ("RandomMating()", RandomMating(), binomial_a(), 0.0, (pm, r) -> random_pm!(pm, r, 0)),
    ("RandomMating(nmale=10)", RandomMating(nmale = 10), 1 / 40 + 1 / 360, 0.0,
        (pm, r) -> random_pm!(pm, r, 10)),
    ("RandomMating(sexratio=0.5)", RandomMating(sexratio = 0.5), 1 / 200 + 1 / 200, 0.0,
        (pm, r) -> random_pm!(pm, r, 50)),
    ("FamilySize(dams=1.0, nmale=20)", FamilySize(dams = 1.0, sires = :random, nmale = 20),
        family_a(FamilySize(dams = 1.0, sires = :random, nmale = 20), 20), 0.0,
        (pm, r) -> family_pm!(pm, r, 20, :gamma, :random)),
    ("FamilySize(equal, nmale=20)", FamilySize(dams = :equal, sires = :equal, nmale = 20),
        family_a(FamilySize(dams = :equal, sires = :equal, nmale = 20), 20), 0.0,
        (pm, r) -> family_pm!(pm, r, 20, :equal, :equal)),
    ("Monoecious()", Monoecious(), 1 / N, 1 / N, (pm, r) -> mono_pm!(pm, r, -1.0)),
    ("Monoecious(selfing=0.2)", Monoecious(selfing = 0.2), 1 / N, 0.2, (pm, r) -> mono_pm!(pm, r, 0.2)),
    ("Monoecious(selfing=0.5)", Monoecious(selfing = 0.5), 1 / N, 0.5, (pm, r) -> mono_pm!(pm, r, 0.5)),
]

# Σ_t k_t 1{k_t ≥ 2} over the T generations for one random gene of each of
# n random individuals: every gene picks its parent from the pedigree and
# one of the parent's two genes at random (one meiosis per transmission).
function branch_length(pm!, n, rng)
    genes = [(i, rand(rng, 1:2)) for i in randperm(rng, N)[1:n]]
    pm = Matrix{Int}(undef, N, 2)
    len = 0
    for _ = 1:T
        length(genes) < 2 && break
        len += length(genes)
        pm!(pm, rng)
        for j in eachindex(genes)
            i, s = genes[j]
            genes[j] = (pm[i, s], rand(rng, 1:2))
        end
        unique!(genes)
    end
    return len
end

# h_B and h_W per bp after T generations from a monomorphic start.
function expected(a, σ)
    hB = hW = 0.0
    for _ = 1:T
        hB, hW = (1 - a) * hB + a / 2 * hW + 2MU, (1 - σ) * hB + σ / 2 * hW + 2MU
    end
    return hB, hW
end

function replicate(m, seed)
    haps, _ = fisher_wright(N, T, fill(CHRLEN, NCHR), MU * 1e8; M = 1e8, mating = m, rng = Xoshiro(seed))
    nh = length(haps)
    cnt = Dict{UInt32,Int}()
    for h in haps, p in h
        cnt[p] = get(cnt, p, 0) + 1
    end
    H = sum(c -> 2 * (c / nh) * (1 - c / nh), values(cnt); init = 0.0) / L
    O = mean(length(symdiff(haps[2i-1], haps[2i])) for i = 1:N) / L
    # One random gene of each of NS random individuals.
    rng = Xoshiro(hash((seed, :sample)))
    s = Dict{UInt32,Int}()
    for i in randperm(rng, N)[1:NS], p in haps[2i-rand(rng, 0:1)]
        s[p] = get(s, p, 0) + 1
    end
    S = count(<(NS), values(s))
    return H, O, S
end

function main(reps = 1000, seed = 2026, breps = 1_000_000)
    println("FisherWright ", pkgversion(FisherWright), ", Julia ", VERSION, ", threads ", Threads.nthreads())
    println("N = $N, $T generations, $NCHR × $(CHRLEN ÷ 1_000_000) Mb, μ = r = $MU per bp, ",
        "$reps replicates per scheme, seed $seed; S in a sample of $NS genes from distinct individuals; ",
        "$breps backward genealogies per scheme and sample size")
    println()
    @printf("%-31s %-8s %24s %22s %7s\n", "scheme", "stat", "FisherWright", "expected", "z")
    hsum = sum(1 / i for i = 1:NS-1)
    failed = String[]
    for (k, (name, m, a, σ, pm!)) in enumerate(SCHEMES)
        out = zeros(reps, 3)
        Threads.@threads for r = 1:reps
            out[r, :] .= replicate(m, hash((seed, k, r)))
        end
        back = zeros(breps, 2)
        Threads.@threads for r = 1:breps
            rng = Xoshiro(hash((seed, :genealogy, k, r)))
            back[r, 1] = branch_length(pm!, 2, rng)
            back[r, 2] = branch_length(pm!, NS, rng)
        end
        hB, hW = expected(a, σ)
        @printf("# %s: a = %.6f, σ = %.3f, Ne = h_B/(4μ) = %.2f, Kingman S_%d = %.2f\n",
            name, a, σ, hB / 4MU, NS, hB * L * hsum)
        rows = (("Σ2pq/bp", out[:, 1], hW / 2N + (1 - 1 / N) * hB, 0.0),
            ("Hobs/bp", out[:, 2], hW, 0.0),
            # The backward simulation (n = 2) against the exact h_B.
            ("h_B/bp*", MU * back[:, 1], hB, 0.0),
            ("S_$NS", out[:, 3], MU * L * mean(back[:, 2]), MU * L * std(back[:, 2]) / sqrt(breps)))
        for (stat, x, e, ese) in rows
            mx, se = mean(x), std(x) / sqrt(length(x))
            z = (mx - e) / sqrt(se^2 + ese^2)
            if startswith(stat, "S")
                fmt, ef = @sprintf("%.2f ± %.2f", mx, se), @sprintf("%.2f ± %.2f", e, ese)
            else
                fmt, ef = @sprintf("%.4e ± %.2e", mx, se), @sprintf("%.4e", e)
            end
            @printf("%-31s %-8s %24s %22s %7.2f\n", name, stat, fmt, ef, z)
            abs(z) > ZMAX && push!(failed, "$name $stat")
        end
        flush(stdout)
    end
    isempty(failed) || (println("\nFAILED: ", join(failed, "; ")); exit(1))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(parse.(Int, ARGS)...)
end
