# Validation of an admixture pulse built with FisherWright.jl through `init`:
# two populations split from one ancestral population, evolve apart, and are
# merged once; the admixed population is then continued.
#
# With FisherWright installed, run from this directory:
#
#   julia -t 8 validate-admixture.jl [reps] [seed]
#
# Model. N = 500 diploids in every population, two sexes (sex Binomial(N, 1/2)
# per generation, conditioned on both sexes), random mating, 2 chromosomes of
# 50 Mb at 1e8 bp per Morgan (0.5 M each, Poisson crossovers, so Haldane's map
# holds exactly), free recombination between chromosomes, μ = 1e-8 per bp.
#
#   1. Ancestral population: 10N = 5000 generations from no variation.
#   2. Split: populations A and B each start from the ancestral haplotypes
#      (`init`) and evolve TSPLIT = 100 generations, with mutation, in separate
#      calls (`result = false`, so sites fixed in one population stay in its
#      haplotypes and are polymorphic after the merge).
#   3. Merge: NA = 150 individuals drawn at random from A and N - NA = 350 from
#      B, as whole individuals (haplotypes 2i - 1 and 2i of a source individual
#      stay together), so α = NA/N = 0.3. A second variant shuffles the 2N
#      merged haplotypes among individuals before continuing.
#   4. Continuation: 20 generations at size N with `mutate = false`, one call
#      per generation, so the same sites are followed throughout.
#
# Sites: polymorphic in the merged population with minor-allele frequency at
# least 0.05, thinned evenly to at most 1000 per chromosome. Every statistic
# is computed on the whole population (2N haplotypes), not a sample.
#
# References.
# 1. Ancestry proportion. With p_A, p_B the allele frequencies of the merged
#    individuals from A and B and δ = p_A - p_B, the least-squares estimate
#      α̂_t = Σ (p_t - p_B) δ / Σ δ²
#    equals α at t = 0 exactly, and E[α̂_t] = α for t > 0, since drift without
#    mutation and selection leaves E[p_t] = p_0.
# 2. Decay of linkage disequilibrium. Let D_t be the LD of a pair of sites at
#    recombination fraction r (Haldane, r = (1 - e^{-2d})/2; r = 1/2 between
#    chromosomes), and C_0 the covariance of the two sites between the two
#    haplotypes of one individual at the merge,
#      C_0 = Σ_i (x_{i1} y_{i2} + x_{i2} y_{i1}) / (2N) - p_x p_y.
#    A gamete is parental with probability 1 - r and joins the alleles of the
#    parent's two haplotypes otherwise, so E[D_1] = λ[(1 - r) D_0 + r C_0]
#    and, from generation 1 on, when individuals are random unions of
#    gametes, E[D_{t+1}] = λ (1 - r) E[D_t], with the drift factor
#    λ = 1 - 1/(2Ne) and Ne from the two-sex recursion (`dioecious_ne`):
#      E[D_t] = λ^t (1 - r)^{t-1} [(1 - r) D_0 + r C_0],   t ≥ 1.
#    A merge of whole individuals has C_0 ≈ α(1 - α) δ_x δ_y, the
#    between-source part of D_0, so that part of the LD does not decay in the
#    first generation (the offspring of the merged individuals are F1, not yet
#    recombinant between the sources); after a shuffle C_0 ≈ 0 and the usual
#    λ^t (1 - r)^t D_0 holds. For each distance bin the check compares
#      S_t = Σ D_t D_0   with   P_t = Σ E[D_t] D_0
#    over the pairs of the bin; E[S_t | merge] = P_t, so S_t - P_t has mean
#    zero over replicates. The ratio S/P is shown for reading; where E[D_t] is
#    a small fraction of D_0 (large r, late t) it is noisy. The prediction
#    neglects terms of order r/N in the drift factor (below 0.1% here).
#
# Each comparison is a z-score over independent replicates (each with its own
# ancestral history), and a check fails when |z| > ZMAX. The two variants share
# their histories and merged individuals, so they are correlated with each
# other; every check is within one variant. The resolution section reports the
# z-score of the prediction that ignores C_0, λ^t (1 - r)^t D_0, for the merge
# of whole individuals. The script exits with status 1 if any check fails.

include(joinpath(@__DIR__, "validate-popgen-theory.jl"))
using Random

const AN, ANC, TSPLIT, NA = 500, 5000, 100, 150
const ACHR, AMR, AM = [50_000_000, 50_000_000], 1.0, 1e8
const AMAF, MAXSITES = 0.05, 1000
const TIMES = [1, 2, 3, 5, 10, 20]
const ABINS = [(0.0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5)]
const BINNAMES = [["$(a)-$(b) M" for (a, b) in ABINS]; "unlinked"]

# 0/1 matrix (2N × m) of the alleles at the sorted sites `pos`.
function allele_matrix(haps, pos)
    X = zeros(length(haps), length(pos))
    for (i, h) in enumerate(haps)
        j = 1
        for x in h
            while j <= length(pos) && pos[j] < x
                j += 1
            end
            j > length(pos) && break
            pos[j] == x && (X[i, j] = 1.0)
        end
    end
    return X
end

# Sites polymorphic in `haps` with MAF ≥ AMAF, at most MAXSITES per chromosome.
function choose_sites(haps, cbp)
    count = Dict{UInt32,Int}()
    for h in haps, x in h
        count[x] = get(count, x, 0) + 1
    end
    n = length(haps)
    all = sort!([x for (x, c) in count if AMAF <= c / n <= 1 - AMAF])
    sites = UInt32[]
    lo = UInt32(0)
    for e in cbp
        s = filter(x -> lo < x <= e, all)
        k = length(s) <= MAXSITES ? eachindex(s) : round.(Int, range(1, length(s), MAXSITES))
        append!(sites, s[k])
        lo = e
    end
    return sites
end

ld(X, p) = X' * X / size(X, 1) - p * p'

function within_cov(X, p)
    X1, X2 = X[1:2:end, :], X[2:2:end, :]
    return (X1' * X2 + X2' * X1) / size(X, 1) - p * p'
end

# Bin of each pair (i < j), 0 for none; recombination fraction of each pair.
function pair_bins(sites, cbp)
    m = length(sites)
    chr = [searchsortedfirst(cbp, x) for x in sites]
    bin, r = zeros(Int, m, m), fill(0.5, m, m)
    for j = 1:m, i = 1:j-1
        if chr[i] != chr[j]
            bin[i, j] = length(ABINS) + 1
        else
            d = abs(Float64(sites[j]) - Float64(sites[i])) / AM
            r[i, j] = (1 - exp(-2d)) / 2
            b = findfirst(((a, c),) -> a < d <= c, ABINS)
            bin[i, j] = b === nothing ? 0 : b
        end
    end
    return bin, r
end

# One replicate: per variant, α̂ at TIMES and (S, P, naive P) per bin and time.
function admixture_replicate(seed, λ)
    rng = Xoshiro(seed)
    α = NA / AN
    anc, cbp = fisher_wright(AN, ANC, ACHR, AMR; M=AM, rng)
    a, _ = fisher_wright(AN, TSPLIT, ACHR, AMR; M=AM, init=anc, rng)
    b, _ = fisher_wright(AN, TSPLIT, ACHR, AMR; M=AM, init=anc, rng)
    pick(pop, ids) = [pop[2i-k] for i in ids for k in (1, 0)]
    whole = vcat(pick(a, randperm(rng, AN)[1:NA]), pick(b, randperm(rng, AN)[1:AN-NA]))
    perm = randperm(rng, 2AN)
    shuffled = whole[perm]
    sites = choose_sites(whole, cbp)
    bin, r = pair_bins(sites, cbp)
    nb = length(ABINS) + 1
    out = map(((whole, 1:2AN .<= 2NA), (shuffled, perm .<= 2NA))) do (pop, isA)
        X0 = allele_matrix(pop, sites)
        p0 = vec(mean(X0; dims=1))
        pA, pB = vec(mean(X0[isA, :]; dims=1)), vec(mean(X0[.!isA, :]; dims=1))
        δ = pA - pB
        α̂(p) = sum((p - pB) .* δ) / sum(δ .^ 2)
        D0, C0 = ld(X0, p0), within_cov(X0, p0)
        α0 = α̂(p0)
        αt = Float64[]
        S, P, Q = zeros(nb, length(TIMES)), zeros(nb, length(TIMES)), zeros(nb, length(TIMES))
        cur = pop
        for t = 1:last(TIMES)
            cur, _ = fisher_wright(AN, 1, ACHR, AMR; M=AM, init=cur, mutate=false, rng)
            k = findfirst(==(t), TIMES)
            k === nothing && continue
            X = allele_matrix(cur, sites)
            p = vec(mean(X; dims=1))
            push!(αt, α̂(p))
            D = ld(X, p)
            for j in axes(bin, 2), i = 1:j-1
                bi = bin[i, j]
                bi == 0 && continue
                ρ = r[i, j]
                S[bi, k] += D[i, j] * D0[i, j]
                P[bi, k] += λ^t * (1 - ρ)^(t - 1) * ((1 - ρ) * D0[i, j] + ρ * C0[i, j]) * D0[i, j]
                Q[bi, k] += λ^t * (1 - ρ)^t * D0[i, j]^2
            end
        end
        u = bin .== nb
        X = α * (1 - α) * δ * δ'
        (; α0, αt, S, P, Q, nsites=length(sites),
            c0=sum(C0[u] .* D0[u]) / sum(D0[u] .^ 2), x0=sum(X[u] .* D0[u]) / sum(D0[u] .^ 2))
    end
    return out
end

function admixture_main(reps=40, seed=2026)
    λ = 1 - 1 / (2 * dioecious_ne(AN))
    α = NA / AN
    println("="^78)
    println(" FisherWright admixture pulse: $reps replicates, seed $seed, $(Threads.nthreads()) threads")
    println(" N = $AN, ancestral $ANC generations, split $TSPLIT generations, α = $α")
    println(" A check fails when |z| > $ZMAX")
    println("="^78)
    res = Vector{Any}(undef, reps)
    t0 = time()
    Threads.@threads for i = 1:reps
        res[i] = admixture_replicate(seed + i, λ)
    end
    @printf(" %d replicates in %.0f s; sites per replicate %d-%d\n", reps, time() - t0,
        extrema(r[1].nsites for r in res)...)
    for (v, name) in enumerate(["merge of whole individuals", "merge of shuffled haplotypes"])
        section[] = name
        R = [r[v] for r in res]
        println("\n", name)
        @printf("  α̂ at the merge: max |α̂₀ - α| = %.1e\n", maximum(abs(x.α0 - α) for x in R))
        @printf("  unlinked pairs, slope on D₀: C₀ %.3f, between-source part α(1 - α)δδ' %.3f\n",
            mean(x.c0 for x in R), mean(x.x0 for x in R))
        println("  ancestry proportion α̂_t (expected $α)")
        println("  t          α̂_t ± SE           z")
        for (k, t) in enumerate(TIMES)
            x = [y.αt[k] for y in R]
            z = (mean(x) - α) / se(x)
            @printf("  %-4d %9.4f ± %.4f   %7.2f  %s\n", t, mean(x), se(x), z, check("α̂ at t = $t", z))
        end
        println("  LD decay: S = Σ D_t D_0 against P = Σ E[D_t] D_0, per distance bin")
        println("  t   distance        S/P      z(S - P)")
        for (k, t) in enumerate(TIMES), b in eachindex(BINNAMES)
            d = [y.S[b, k] - y.P[b, k] for y in R]
            z = mean(d) / se(d)
            ratio = mean(y.S[b, k] for y in R) / mean(y.P[b, k] for y in R)
            @printf("  %-3d %-13s %7.3f   %7.2f  %s\n", t, BINNAMES[b], ratio, z,
                check("LD t = $t, $(BINNAMES[b])", z))
        end
        if v == 1
            println("\n  Resolution: z of the prediction that ignores C₀, λ^t (1 - r)^t D₀")
            println("  t   distance      S/P_naive   z(S - P_naive)")
            for (k, t) in enumerate(TIMES[1:3]), b in eachindex(BINNAMES)
                d = [y.S[b, k] - y.Q[b, k] for y in R]
                @printf("  %-3d %-13s %7.3f   %7.2f\n", t, BINNAMES[b],
                    mean(y.S[b, k] for y in R) / mean(y.Q[b, k] for y in R), mean(d) / se(d))
            end
        end
    end
    println("\n", "="^78)
    println(isempty(failures) ? " All checks passed." : " FAILED: " * join(failures, "; "))
    println("="^78)
    return isempty(failures)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse.(Int, ARGS)
    exit(admixture_main(args...) ? 0 : 1)
end
