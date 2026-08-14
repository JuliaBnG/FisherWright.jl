# Results and export

## Default result

The default simulation result is a tuple:

```julia
haplotypes, chromosome_ends = fisher_wright(100, 100, [1_000_000], 1.0)
```

Every element of `haplotypes` is a sorted, unique `Vector{UInt32}`. A
mutation is present when its genome-wide position occurs in the haplotype.

Use `muts2bitarray` to obtain a dense representation:

```@example results
using FisherWright

mutations, chromosome_ends = fisher_wright(20, 10, [100_000], 1.0)
alleles, loci = muts2bitarray(mutations, chromosome_ends)
(size(alleles), propertynames(loci))
```

`alleles` is a `BitMatrix` with loci in rows and haplotypes in columns.
`loci` has these columns:

| Column | Meaning |
|---|---|
| `chr` | One-based chromosome index. |
| `pos` | Genome-wide base-pair position. |
| `ref` | Randomly generated reference base. |
| `alt` | Randomly generated alternate base, distinct from `ref`. |
| `frq` | Alternate-allele frequency across haplotypes. |

Fixed loci are removed by default. Pass `flip=true` to randomly invert the
reference/alternate coding at each locus.

## Structured result and substitutions

Set `result=true` to return a `FisherWrightResult`. At each
`fixation_interval`, mutations carried by every active haplotype are removed
from the active population and recorded as substitutions:

```@example results
result = FisherWrightResult(
    [UInt32[1, 3], UInt32[2], UInt32[], UInt32[1, 2, 3]],
    UInt32[3],
    UInt32[],
)
(length(result.active_haplotypes), result.chromosome_ends, result.substitutions)
```

This avoids retaining fixed mutations in every haplotype during long
simulations. The result fields are:

| Field | Meaning |
|---|---|
| `active_haplotypes` | Polymorphic mutation positions after fixed loci are removed. |
| `chromosome_ends` | Cumulative chromosome endpoints. |
| `substitutions` | Sorted fixed mutation positions. |

Convert a structured result to `BnGStructs.Haplotype` and a linkage map:

```@example results
haplotype, loci = to_haplotype(result)
(haplotype.nlc, haplotype.nhp, length(loci.pos))
```

Use `to_haplotype(result; include_fixed=true)` when fixed substitutions should
also appear in the dense output. If no polymorphic loci are available,
`to_haplotype` raises an error rather than returning an empty `Haplotype`.

## Chip extraction

Extract a selected set of marker coordinates directly from a structured result
without first materializing every simulated locus:

```@example results
chip_positions = UInt32[10_000, 50_000, 90_000]
chip_haplotype = to_haplotype(result, chip_positions)
size(chip_haplotype)
```

`chip_positions` must be sorted and unique. Each requested position is
retained even when it is absent from the population, in which case its row is
all false. `extract_chip_bitarray(result.active_haplotypes, chip_positions)`
returns the corresponding `BitMatrix` when a `BnGStructs.Haplotype` wrapper is
not needed. A `BnGStructs.LocusSet` can also select indices from a shared,
sorted coordinate vector:

```julia
chip_haplotype = to_haplotype(result, locus_set, all_positions)
```
