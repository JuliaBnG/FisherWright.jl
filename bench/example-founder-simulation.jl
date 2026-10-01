using FisherWright
using BnGStructs
using StatsBase
using Random
using Printf

function main()
    println("=================================================================")
    println(" FisherWright.jl Application Example: Founder Simulation & Chip")
    println("=================================================================")
    chr_lengths = fill(100_000_000, 10)
    N = 1000
    Nt = 10 * N
    println("Species architecture   : 10 autosomes × 100 Mb (1 Gb total)")
    println("Population size (N)     : $N diploids (2,000 haplotypes)")
    println("Burn-in generations (Nt): $Nt (10N)")
    println("Threads                 : $(Threads.nthreads())")
    println("-----------------------------------------------------------------")

    # JIT warm-up
    fisher_wright(10, 2, [10_000], 1.0; result=true, verbose=false)

    # 1. Forward simulation
    t0 = time()
    result = fisher_wright(N, Nt, chr_lengths, 1.0; result=true, fixation_interval=50, verbose=false)
    t_sim = time() - t0
    @printf("Simulation runtime      : %.1f s\n", t_sim)

    # 2. WGS ascertainment
    t0 = time()
    _, loci = to_haplotype(result)
    poly_candidates = loci[0.05 .<= loci.frq .<= 0.95, :pos]
    chip_markers = sort(sample(poly_candidates, min(50_000, length(poly_candidates)); replace=false))
    t_asc = time() - t0
    @printf("Ascertainment runtime   : %.2f s (found %d MAF candidates)\n", t_asc, length(poly_candidates))

    # 3. Targeted chip extraction
    t0 = time()
    chip_haplotypes = to_haplotype(result, chip_markers)
    t_ext = time() - t0
    @printf("Chip extraction runtime : %.2f s\n", t_ext)
    println("Extracted $(chip_haplotypes.nlc) chip markers across $(chip_haplotypes.nhp) haplotypes.")
    println("=================================================================")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
