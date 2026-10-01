using FisherWright
using Random
using Printf

function peak_rss_mb()
    if Sys.islinux()
        for line in eachline("/proc/self/status")
            if startswith(line, "VmHWM:")
                kb = parse(Int, split(line)[2])
                return kb / 1024.0
            end
        end
    end
    return -1.0
end

function main()
    ne = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
    nt = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
    nchr = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10
    chr_len = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 100_000_000 # 100 Mb per autosome
    mr = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1.0 # 1.0 per 1e8 bp = 1e-8 / bp
    seed = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 42

    total_len = nchr * chr_len
    chr_vec = fill(chr_len, nchr)

    println("=================================================================")
    println(" FisherWright.jl Simulation")
    println("=================================================================")
    println("Diploid Individuals (ne): $ne (haplotypes: $(2 * ne))")
    println("Generations (nt)        : $nt")
    println("Autosomes (nchr)        : $nchr autosomes × $(chr_len / 1e6) Mb = $(total_len / 1e6) Mb total")
    println("Mutation rate           : $(mr * 1e-8) / bp / generation")
    println("Recombination rate      : 1e-8 / bp / generation (1 cM/Mb)")
    println("Independent assortment  : Yes (0.5 crossover probability at autosome ends)")
    println("Threads                 : $(Threads.nthreads())")
    println("Random seed             : $seed")
    println("-----------------------------------------------------------------")

    Random.seed!(seed)
    GC.gc()

    # Warm-up JIT compilation on a minimal population to measure pure execution
    let warmup_res = fisher_wright(10, 2, [10_000], 1.0; result=true, fixation_interval=1, verbose=false)
        muts2bitarray(warmup_res.active_haplotypes, warmup_res.chromosome_ends)
    end
    GC.gc()

    # 1. Forward-in-time Wright-Fisher simulation
    t_sim = @elapsed res = fisher_wright(
        ne,
        nt,
        chr_vec,
        mr;
        M = 1e8,
        mut_base = 1e8,
        result = true,
        fixation_interval = 10,
        verbose = false,
    )

    # 2. Conversion to BitMatrix
    t_bit = @elapsed bitmat, locus_map = muts2bitarray(res.active_haplotypes, res.chromosome_ends)

    n_loci, n_haps = size(bitmat)
    bitmat_bytes = Base.summarysize(bitmat)
    bitmat_mb = bitmat_bytes / (1024 * 1024)
    peak_rss = peak_rss_mb()

    @printf("Simulation runtime      : %.3f s\n", t_sim)
    @printf("BitMatrix conversion    : %.3f s\n", t_bit)
    @printf("Total runtime           : %.3f s\n", t_sim + t_bit)
    @printf("Polymorphic SNPs found  : %d\n", n_loci)
    @printf("Fixed substitutions     : %d\n", length(res.substitutions))
    @printf("BitMatrix dimension     : %d loci × %d haplotypes\n", n_loci, n_haps)
    @printf("BitMatrix memory size   : %.2f MiB (%.2f bits/element)\n", bitmat_mb, (bitmat_bytes * 8) / (n_loci * n_haps))
    @printf("Process peak RSS        : %.2f MiB\n", peak_rss)
    println("=================================================================")

    return bitmat, locus_map
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
