# Detect sgRNAs at one or more reference stages

Classifies each targeting sgRNA as detected or not at each stage of the
screen (e.g. plasmid, T0), and summarises sgRNA loss per gene. Use it to
see how much of the library is lost before selection starts.

## Usage

``` r
qc_guide_representation(
  counts,
  stages,
  ctrl_pattern = "NEG_CTRL",
  min_reads = 30
)
```

## Arguments

- counts:

  Count matrix: data frame with columns `sgRNA`, `Gene` and one column
  of raw counts per sample (see section 'Count matrix').

- stages:

  Named list; each element gives the sample IDs (columns of `counts`)
  pooled for that stage, e.g.
  `list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3"))`. The order
  of the list sets the order of stages.

- ctrl_pattern:

  Regular expression matched against the `Gene` column to identify
  control sgRNAs (default is "NEG_CTRL", which matches the controls of
  all libraries shipped with crisprFlow).

- min_reads:

  An sgRNA is detected at a stage if its mean count per sample in that
  stage is at least `min_reads` (default is 30).

## Value

A list with

- guide_summary:

  one row per sgRNA and stage: `sgRNA`, `Gene`, `stage`, `mean_count`,
  `detected`

- gene_summary:

  one row per gene and stage: `n_guides`, `n_detected`, `n_lost`,
  `frac_lost`

- threshold:

  the `min_reads` used

Pass it to [`qc_plot_guide_detection`](qc_plot_guide_detection.md) and
[`qc_genes_with_lost_guides`](qc_genes_with_lost_guides.md).

## Examples

``` r
rep_data <- qc_guide_representation(
  exampleCounts,
  stages = list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3")))
head(rep_data$gene_summary)
#> # A tibble: 6 × 6
#>   stage   Gene   n_guides n_detected n_lost frac_lost
#>   <fct>   <chr>     <int>      <int>  <int>     <dbl>
#> 1 Plasmid AAK1          5          5      0         0
#> 2 Plasmid ABL2          5          5      0         0
#> 3 Plasmid ACTB          5          5      0         0
#> 4 Plasmid ACTR5         5          5      0         0
#> 5 Plasmid ACVR1C        5          5      0         0
#> 6 Plasmid ADH1B         5          5      0         0
```
