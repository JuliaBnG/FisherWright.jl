# Cost of selection in fisher_wright as a function of the number of
# selected sites K and the rate of new selected alleles per haplotype.
#
#   julia -t 8 --project=. bench/benchmark-selection.jl [ne = 2000] [nt = 500] [reps = 2] [seed = 42]
#
# Genome and neutral process as in bench/simulate-fisher-wright.jl: 10 × 100
# Mb, μ = r = 1e-8 per bp per generation. Selection: `fitness_loci` with K
# sites drawn uniformly, s = -0.01, h = 0.5, and a per-site rate μ_sel such
# that K·μ_sel new selected alleles arise per haplotype per generation. Each
# configuration is timed `reps` times, interleaved with the neutral run, and
# the table reports the median wall time and its ratio to the neutral median.

using FisherWright, Random, Printf, Statistics

function run(ne, nt, fitness, seed)
    GC.gc()
    return @elapsed fisher_wright(ne, nt, fill(100_000_000, 10), 1.0; M = 1e8, mut_base = 1e8,
        result = true, fixation_interval = 10, fitness = fitness, rng = Xoshiro(seed))
end

function main(ne = 2000, nt = 500, reps = 2, seed = 42)
    chr = fill(100_000_000, 10)
    grid = [(K, u) for K in (1_000, 10_000, 100_000) for u in (0.01, 0.1)]
    models = Pair{String,Any}["neutral" => Neutral()]
    for (K, u) in grid
        loci = fitness_loci(chr, K; s = -0.01, h = 0.5, μ = u / K, rng = Xoshiro(seed))
        push!(models, @sprintf("K = %6d, Kμ = %4.2f", K, u) => SiteFitness(loci))
    end
    # Compile every method on a small run first.
    for (_, f) in models
        fisher_wright(10, 2, [10_000], 1.0; result = true, fitness = f isa SiteFitness ?
            SiteFitness(FitnessLoci([5_000], -0.01; μ = 1e-3)) : f)
    end
    println("="^72)
    println(" Selection overhead: ne = $ne, nt = $nt, 10 × 100 Mb, threads = $(Threads.nthreads()), reps = $reps")
    println("="^72)
    t = Dict(name => Float64[] for (name, _) in models)
    for r = 1:reps, (name, f) in models
        push!(t[name], run(ne, nt, f, seed + r))
    end
    base = median(t["neutral"])
    @printf("  %-24s %12s %10s   %s\n", "model", "median (s)", "ratio", "runs (s)")
    for (name, _) in models
        @printf("  %-24s %12.2f %10.3f   %s\n", name, median(t[name]), median(t[name]) / base,
            join((@sprintf("%.2f", x) for x in t[name]), ", "))
    end
    println("="^72)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(parse.(Int, ARGS)...)
end
