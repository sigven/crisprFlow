# Genes that lost many of their sgRNAs

Lists genes where at least a given fraction of sgRNAs were not detected
at a stage. Hits for such genes rest on few sgRNAs and should be
interpreted with care.

## Usage

``` r
qc_genes_with_lost_guides(rep_data, stage, min_frac = 0.4)
```

## Arguments

- rep_data:

  Output of [`qc_guide_representation`](qc_guide_representation.md).

- stage:

  Name of the stage (one of the names of `stages` given to
  [`qc_guide_representation`](qc_guide_representation.md)).

- min_frac:

  Minimum fraction of lost sgRNAs (default is 0.4, e.g. 2 of 5).

## Value

A tibble from `rep_data$gene_summary`, one row per gene, sorted by
fraction and number of lost sgRNAs.

## Examples

``` r
rep_data <- qc_guide_representation(
  exampleCounts,
  stages = list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3")))
qc_genes_with_lost_guides(rep_data, stage = "T0")
#> # A tibble: 5 × 6
#>   stage Gene  n_guides n_detected n_lost frac_lost
#>   <fct> <chr>    <int>      <int>  <int>     <dbl>
#> 1 T0    AKAP1        5          3      2       0.4
#> 2 T0    EDNRB        5          3      2       0.4
#> 3 T0    FTSJ3        5          3      2       0.4
#> 4 T0    MDM4         5          3      2       0.4
#> 5 T0    MERTK        5          3      2       0.4
```
