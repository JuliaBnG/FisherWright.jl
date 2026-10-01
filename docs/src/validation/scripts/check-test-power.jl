# Power check for the population genetics tests in test/runtests.jl.
#
# Run from the repository root:
#
#   julia -t 4 --project=. docs/src/validation/scripts/check-test-power.jl
#
# A test is only useful if it passes for the correct model and fails for a
# wrong one. This script takes the block between the "BEGIN/END population
# genetics validation" markers in test/runtests.jl and runs it unchanged for
# several seeds, then with one known error introduced at a time:
#
#   - simulated Ne halved (the theory still assumes Ne = 100);
#   - recombination rate doubled (M halved), in the simulation and in the
#     crossover test;
#   - the pre-v0.3.5 generation loop (src/fwp.jl from commit 4a24282~1), which
#     mutated parents before mating.
#
# For each case it prints pass/total per testset and the evaluated expression
# of every failed test.

using Test, Printf

const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))
const PREFIX_COMMIT = "4a24282~1"
const SIM = "fisher_wright(ne, 20ne, chr, 1.0; result = true)"
const XMAP = "rmap = uniform_recombination_map(chr; M=M)"

# Records results without printing or throwing, so every case runs to the end.
mutable struct Collect <: Test.AbstractTestSet
    description::String
    results::Vector{Any}
end
Collect(desc; kw...) = Collect(desc, Any[])
Test.record(ts::Collect, r) = (push!(ts.results, r); r)
function Test.finish(ts::Collect)
    Test.get_testset_depth() > 0 && Test.record(Test.get_testset(), ts)
    return ts
end

"The validation block of test/runtests.jl, split into setup code and top-level testsets."
function validation_block()
    src = read(joinpath(ROOT, "test", "runtests.jl"), String)
    a = findfirst("# BEGIN population genetics validation", src)
    b = findfirst("# END population genetics validation", src)
    (a === nothing || b === nothing) && error("validation markers not found in test/runtests.jl")
    code = src[first(a):first(b)-1]
    i = first(findfirst("\n@testset ", code))
    return code[1:i], split(strip(code[i:end]), r"\n(?=@testset )")
end

"Source of the package with src/fwp.jl replaced by the pre-fix version."
function prefix_source()
    dir = mktempdir()
    cp(joinpath(ROOT, "src"), joinpath(dir, "src"))
    write(joinpath(dir, "src", "fwp.jl"), read(`git -C $ROOT show $PREFIX_COMMIT:src/fwp.jl`, String))
    return joinpath(dir, "src", "FisherWright.jl")
end

function run_case(setup, sets; edits = (), package = nothing)
    m = Module()
    Core.eval(m, :(using Test, Random, Statistics))
    Core.eval(m, :(const Collect = $Collect))
    if package === nothing
        Core.eval(m, :(using FisherWright: fisher_wright, muts2bitarray, uniform_recombination_map, cobp!))
    else
        Base.include(m, package)
        Core.eval(m, :(using .FisherWright: fisher_wright, muts2bitarray, uniform_recombination_map, cobp!))
    end
    Base.include_string(m, setup)
    return map(sets) do s
        s = replace(s, "@testset \"" => "@testset Collect \"")
        for e in edits
            s = replace(s, e)
        end
        Base.include_string(m, s)
    end
end

function report(label, results)
    println(label)
    function show_set(ts, path)
        tests = filter(r -> !(r isa Collect), ts.results)
        if !isempty(tests)
            npass = count(r -> r isa Test.Pass, tests)
            @printf("  %-72s %2d/%-2d %s\n", path, npass, length(tests), npass == length(tests) ? "pass" : "FAIL")
            for r in tests
                r isa Test.Fail && println("      failed: ", r.data)
                r isa Test.Error && println("      error: ", first(split(sprint(show, r), '\n')))
            end
        end
        foreach(r -> r isa Collect && show_set(r, path * " / " * r.description), ts.results)
    end
    foreach(ts -> show_set(ts, ts.description), results)
    println()
end

function main()
    setup, sets = validation_block()
    any(s -> occursin(SIM, s), sets) || error("simulation call not found in the validation block")
    any(s -> occursin(XMAP, s), sets) || error("recombination map call not found in the validation block")

    println("Power check of the population genetics tests in test/runtests.jl")
    println("Julia $(VERSION), threads = $(Threads.nthreads())\n")
    report("Correct model, CI seeds (2026 for the theory tests, 42 for the crossover test)", run_case(setup, sets))
    for s in (1, 7, 42, 99)
        seeds = ("Random.seed!(2026)" => "Random.seed!($s)", "Random.seed!(42)" => "Random.seed!($s)")
        report("Correct model, seed $s", run_case(setup, sets; edits = seeds))
    end
    report("Known error: simulated Ne halved (50); theory assumes 100",
        run_case(setup, sets; edits = (SIM => replace(SIM, "fisher_wright(ne," => "fisher_wright(ne ÷ 2,"),)))
    report("Known error: recombination rate doubled (M halved)",
        run_case(setup, sets; edits = (SIM => replace(SIM, "result = true" => "M = 5e7, result = true"),
            XMAP => replace(XMAP, "M=M)" => "M=M / 2)"))))
    report("Known error: pre-v0.3.5 generation loop (src/fwp.jl from $PREFIX_COMMIT)",
        run_case(setup, sets; package = prefix_source()))
end

main()
