using FisherWright
using Random
using Statistics
using Printf

function parse_chr(s::String)
    parts = split(s, ',')
    return parse.(Int, parts)
end

function one_run(ne::Int, nt::Int, chr::Vector{Int}, mr::Float64)
    GC.gc()
    trial = @timed begin
        Random.seed!(1234)
        fisher_wright(ne, nt, chr, mr; result=true, fixation_interval=1)
    end
    return trial.time, trial.bytes
end

function main()
    ne = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    nt = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    chr = length(ARGS) >= 3 ? parse_chr(ARGS[3]) : [100_000_000, 100_000_000, 100_000_000]
    mr = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    reps = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 3

    println("FisherWright phase-3 benchmark")
    println("julia_version=$(VERSION)")
    println("ne=$ne nt=$nt chr=$(join(chr, ',')) mr=$mr reps=$reps")

    # Warm-up to reduce one-time compilation impact in repeated runs.
    one_run(ne, nt, chr, mr)

    times = Float64[]
    bytes = Int[]
    for _ = 1:reps
        t, b = one_run(ne, nt, chr, mr)
        push!(times, t)
        push!(bytes, b)
    end

    @printf("time_s median=%.6f mean=%.6f min=%.6f max=%.6f\n", median(times), mean(times), minimum(times), maximum(times))
    @printf("bytes  median=%d mean=%.1f min=%d max=%d\n", Int(median(bytes)), mean(bytes), minimum(bytes), maximum(bytes))
end

main()