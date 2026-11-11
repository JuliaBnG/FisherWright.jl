# Validation of selection in FisherWright.jl: fixation probability of one new
# copy of a selected allele, simulated with the package's engine, against an
# independent naive simulation of the same model (fixation-reference.jl).
#
# With FisherWright installed, run from this directory:
#
#   julia -t 8 validate-selection.jl [reps = 200_000] [ref_reps = 2_000_000] [seed = 2026]
#
# Model. N diploids, two sexes, one selected site with genotype fitness 1,
# 1 + h*s, 1 + s. Each generation every offspring draws its sire among the
# males and its dam among the females with probability proportional to
# fitness, with replacement (soft selection within each sex). Sex is
# Binomial(N, 1/2) per generation conditioned on both sexes, or exactly N/2
# males. A replicate starts with one copy on a random haplotype of an
# otherwise monomorphic population, without further mutation, and runs to
# loss or fixation.
#
# Each case is a z-score of the difference of two proportions, engine against
# naive reference; a check fails when |z| > ZMAX, and the script then exits
# with status 1. The ideal Wright-Fisher value is printed for comparison
# only: the two-sex model is below it by terms of order 1/N.

include(joinpath(@__DIR__, "fixation-reference.jl"))
using FisherWright: FisherWright, SimState, _simulate!, ConstantSize, RandomMating, SiteFitness, FitnessLoci,
    uniform_recombination_map
using Statistics

const ZMAX = 3.5

# One copy introduced on a random haplotype; the engine advances one
# generation at a time until the allele is lost or fixed.
function engine_fixation(N, s, h, reps, rng; mating = RandomMating())
    rmap = uniform_recombination_map([100]; M = 1e8)
    f = SiteFitness(FitnessLoci([50], s; h = h))
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], N, rng)
    fixed = 0
    for _ = 1:reps
        foreach(empty!, st.sel)
        push!(st.sel[rand(rng, 1:2N)], UInt32(50))
        while true
            _simulate!(st, ConstantSize(N), mating, f, 1, rmap, Poisson(1.0), UInt32(1):UInt32(100);
                mutate = false, track = false, fixation_interval = 1, verbose = false, progress_interval = 100,
                rng = rng)
            c = count(!isempty, view(st.sel, 1:2N))
            c == 0 && break
            if c == 2N
                fixed += 1
                break
            end
        end
        resize!(st.census, 1)
        empty!(st.males)
    end
    return fixed / reps
end

# (N, s, h, equal sexes)
const CASES = [
    (30, 0.2, 0.5, false), (30, 0.2, 0.0, false), (30, 0.2, 1.0, false), (30, 0.0, 0.5, false),
    (30, -0.1, 0.5, false), (50, 0.1, 0.5, false), (50, 0.1, 0.5, true), (50, 0.02, 0.5, false),
]

function selection_main(reps = 200_000, ref_reps = 2_000_000, seed = 2026)
    println("="^86)
    println(" FisherWright.jl validation of selection: fixation probability of one new copy")
    println(" engine reps = $reps, naive reference reps = $ref_reps, seed = $seed, threads = $(Threads.nthreads())")
    println(" A check fails when |z| > $ZMAX")
    println("="^86)
    @printf("  %4s %6s %5s %-9s %21s %21s %7s  %-5s %9s\n", "N", "s", "h", "sexes", "engine", "naive reference", "z", "",
        "ideal WF")
    rng = Xoshiro(seed)
    ok = true
    for (k, (N, s, h, equal)) in enumerate(CASES)
        mating = equal ? RandomMating(sexratio = 0.5) : RandomMating()
        u = engine_fixation(N, s, h, reps, rng; mating = mating)
        r, er = naive(N, s, h, ref_reps, seed + k; nmale = equal ? N ÷ 2 : 0)
        eu = sqrt(r * (1 - r) / reps)
        z = (u - r) / sqrt(eu^2 + er^2)
        pass = abs(z) <= ZMAX
        ok &= pass
        @printf("  %4d %6.2f %5.2f %-9s %10.5f ± %-8.5f %10.5f ± %-8.5f %7.2f  %-5s %9.5f\n", N, s, h,
            equal ? "N/2 males" : "binomial", u, eu, r, er, z, pass ? "pass" : "FAIL", wf_exact(2N, s, h))
    end
    println("="^86)
    println(ok ? " All checks passed." : " FAILED")
    println("="^86)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(selection_main(parse.(Int, ARGS)...) ? 0 : 1)
end
