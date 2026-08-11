using FisherWright
using Random
using Statistics
using Printf

function parse_chr(s::String)
    parts = split(s, ',')
    return parse.(Int, parts)
end

"""
Peak resident set size of this process in KiB, or -1 where /proc is absent.

`@timed` reports *cumulative* bytes allocated, which is not a memory footprint;
this is the number to compare against `/usr/bin/time -f '%M'`.
"""
function peak_rss_kb()
    Sys.islinux() || return -1
    for line in eachline("/proc/self/status")
        startswith(line, "VmHWM:") && return parse(Int, split(line)[2])
    end
    return -1
end

function one_run(ne::Int, nt::Int, chr::Vector{Int}, mr::Float64, seed::Int)
    GC.gc()
    trial = @timed begin
        Random.seed!(seed)
        fisher_wright(ne, nt, chr, mr; result=true, fixation_interval=1)
    end
    return trial.time, trial.bytes, trial.gctime
end

function main()
    ne = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    nt = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    chr = length(ARGS) >= 3 ? parse_chr(ARGS[3]) : [100_000_000, 100_000_000, 100_000_000]
    mr = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    reps = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 3
    seed = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 1234

    println("FisherWright phase-3 benchmark")
    println("julia_version=$(VERSION)")
    # Thread count is part of the result: seeded runs only reproduce for a
    # fixed nthreads, and wall time depends on it directly.
    println("threads=$(Threads.nthreads())")
    println("ne=$ne nt=$nt chr=$(join(chr, ',')) mr=$mr reps=$reps seed=$seed")

    # Warm-up on a tiny problem: it compiles the same methods as the real run
    # without paying for a second full-size simulation.
    one_run(10, 5, chr, mr, seed)

    times = Float64[]
    bytes = Int[]
    gctimes = Float64[]
    for _ = 1:reps
        t, b, gct = one_run(ne, nt, chr, mr, seed)
        push!(times, t)
        push!(bytes, b)
        push!(gctimes, gct)
    end

    @printf("time_s median=%.6f mean=%.6f min=%.6f max=%.6f\n",
        median(times), mean(times), minimum(times), maximum(times))
    @printf("bytes  median=%d mean=%.1f min=%d max=%d\n",
        Int(median(bytes)), mean(bytes), minimum(bytes), maximum(bytes))
    @printf("gc_pct median=%.2f\n", 100 * median(gctimes ./ times))
    @printf("peak_rss_kb=%d\n", peak_rss_kb())
end

main()
