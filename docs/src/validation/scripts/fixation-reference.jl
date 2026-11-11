# Reference fixation probabilities for selection in FisherWright.jl.
#
#   julia -t 8 fixation-reference.jl [reps = 10_000_000] [seed = 2026]
#
# FisherWright's selection is soft and acts within each sex: every offspring
# draws its sire among the males and its dam among the females with
# probability proportional to fitness, with replacement, and sex is
# Binomial(N, 1/2) per generation, conditioned on both sexes (or a fixed
# number of males). This differs from the ideal Wright-Fisher population
# behind Kimura's diffusion formula by terms of order 1/N, so the reference
# for the implementation is an independent simulation of the same model:
#
# `naive_fixation` follows one locus as genotype counts 0, 1, 2 per
# individual, deliberately naive and sharing no code with the package: sex by
# coin flips, parents by rejection sampling on fitness, one Bernoulli per
# transmitted gene. It starts from one copy in a random individual and runs
# to loss or fixation.
#
# Part 1 gives the reference values used by the CI test "Fixation
# probability under selection" (test/runtests.jl), with their standard
# errors. Part 2 shows how the model relates to theory: the exact fixation
# probability of the ideal Wright-Fisher chain with 2N genes (`wf_exact`),
# Kimura's diffusion value, and the two-sex model with binomial and with
# equal sexes, for 2Ns = 10 at N = 50 and N = 500.

using Random, Printf, LinearAlgebra, Distributions

# Fitness of genotypes 0, 1, 2: 1, 1 + h*s, 1 + s. `nmale = 0` draws sex as
# Bernoulli(1/2) per individual until both sexes are present.
function naive_fixation(N, s, h, reps, rng; nmale = 0)
    w = (1.0, 1 + h * s, 1 + s)
    wmax = maximum(w)
    g, g2 = zeros(Int, N), zeros(Int, N)
    males, females = Int[], Int[]
    function parent(pool)
        while true
            i = rand(rng, pool)
            rand(rng) * wmax < w[g[i]+1] && return i
        end
    end
    fixed = 0
    for _ = 1:reps
        fill!(g, 0)
        g[rand(rng, 1:N)] = 1
        while true
            empty!(males)
            empty!(females)
            if nmale == 0
                while isempty(males) || isempty(females)
                    empty!(males)
                    empty!(females)
                    for i = 1:N
                        push!(rand(rng, Bool) ? males : females, i)
                    end
                end
            else
                perm = randperm(rng, N)
                append!(males, view(perm, 1:nmale))
                append!(females, view(perm, nmale+1:N))
            end
            # A sex without any fit individual cannot occur for s > -1.
            for i = 1:N
                a, b = g[parent(males)], g[parent(females)]
                g2[i] = (rand(rng) * 2 < a) + (rand(rng) * 2 < b)
            end
            g, g2 = g2, g
            c = sum(g)
            c == 0 && break
            if c == 2N
                fixed += 1
                break
            end
        end
    end
    return fixed
end

# `reps` replicates split over the threads, each with its own stream.
function naive(N, s, h, reps, seed; nmale = 0)
    nt = Threads.nthreads()
    per = cld(reps, nt)
    tasks = [Threads.@spawn naive_fixation(N, s, h, per, Xoshiro(hash((seed, k))); nmale = nmale) for k = 1:nt]
    R = per * nt
    p = sum(fetch, tasks) / R
    return p, sqrt(p * (1 - p) / R)
end

# Exact fixation probability of one copy in an ideal Wright-Fisher
# population of M = 2N genes: random union of gametes, genotype fitness
# 1, 1 + h*s, 1 + s, binomial sampling of M genes.
function wf_exact(M, s, h)
    P = zeros(M + 1, M + 1)
    for i = 0:M
        p = i / M
        wbar = (1 - p)^2 + 2p * (1 - p) * (1 + h * s) + p^2 * (1 + s)
        q = (p^2 * (1 + s) + p * (1 - p) * (1 + h * s)) / wbar
        P[i+1, :] = pdf.(Binomial(M, q), 0:M)
    end
    return ((I-P[2:M, 2:M])\P[2:M, M+1])[1]
end

# Kimura (1962) for h = 1/2, one copy among 2N genes, Ne = N.
kimura(N, s) = s == 0 ? 1 / 2N : (1 - exp(-s)) / (1 - exp(-2N * s))

# Cases of the CI test: (N, s, h).
const CI_CASES = [(30, 0.2, 0.5), (30, 0.2, 0.0), (30, 0.2, 1.0), (30, 0.0, 0.5), (30, -0.1, 0.5)]

function main(reps = 10_000_000, seed = 2026)
    println("="^78)
    println(" Fixation probability of one new copy: naive two-sex reference")
    println(" reps = $reps per case, seed = $seed, threads = $(Threads.nthreads())")
    println("="^78)
    println("\n[1/2] Reference values for the CI test (binomial sex)")
    @printf("  %4s %7s %5s %12s %10s %12s\n", "N", "s", "h", "u", "SE", "ideal WF")
    for (N, s, h) in CI_CASES
        p, e = naive(N, s, h, reps, seed)
        @printf("  %4d %7.3f %5.2f %12.6f %10.6f %12.6f\n", N, s, h, p, e, wf_exact(2N, s, h))
    end
    println("\n[2/2] The two-sex model against the ideal Wright-Fisher population, 2Ns = 10, h = 1/2")
    @printf("  %5s %7s %10s %10s %22s %22s\n", "N", "s", "Kimura", "ideal WF", "two sexes, binomial", "two sexes, N/2 males")
    for (N, s) in ((50, 0.1), (500, 0.01))
        pb, eb = naive(N, s, 0.5, reps ÷ (N ÷ 50), seed + 1)
        pf, ef = naive(N, s, 0.5, reps ÷ (N ÷ 50), seed + 2; nmale = N ÷ 2)
        @printf("  %5d %7.3f %10.5f %10.5f %12.5f ± %-7.5f %12.5f ± %-7.5f\n", N, s, kimura(N, s),
            wf_exact(2N, s, 0.5), pb, eb, pf, ef)
    end
    println("="^78)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(parse.(Int, ARGS)...)
end
