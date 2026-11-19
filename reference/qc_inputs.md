# Inputs shared by the QC functions

The `qc_*` functions check the quality of pooled CRISPR screen count
data. Rather than parsing sample names, they are driven by two tables: a
count matrix (`counts`) and a sample table (`meta`). The same functions
therefore work for screens with different designs.

## Arguments

- counts:

  Count matrix: data frame with columns `sgRNA`, `Gene` and one column
  of raw counts per sample (see section 'Count matrix').

- meta:

  Sample table: data frame with one row per sample and columns
  `sample_id`, `label`, `group` and `replicate` (see section 'Sample
  table'). Only samples in `meta` are used.

- ctrl_pattern:

  Regular expression matched against the `Gene` column to identify
  control sgRNAs (default is "NEG_CTRL", which matches the controls of
  all libraries shipped with crisprFlow).

- summary_df:

  Count summary: data frame with columns `sample_id`, `total_reads` and
  `mapped` (see section 'Count summary').

- colour_by:

  Name of the column in `meta` used to colour samples.

- colours:

  Optional named character vector of colours, with names matching the
  values of the `colour_by` column (e.g.
  `c(Untreated = "#0073C2", Drug = "#CD534C")`). If NULL (default), the
  RColorBrewer "Dark2" palette is used.

- facet_by:

  Optional name of a column in `meta` used to split the plot into panels
  (default is NULL, no panels).

- interactive:

  Logical; if TRUE (default), return an interactive plotly object,
  otherwise a ggplot object.

- base_size:

  Base font size of the plot.

## Count matrix (`counts`)

A data frame with one row per sgRNA and the columns

- sgRNA:

  sgRNA identifier

- Gene:

  gene targeted by the sgRNA (control sgRNAs carry a gene name matching
  `ctrl_pattern`, e.g. `NONE_TARGETING_NEG_CTRL` or `CTRLA_NEG_CTRL`)

- \<sample_id\>:

  one column of raw (non-normalized) read counts per sample

This is the layout written by MAGeCK and guide-counter, see
[`qc_read_counts`](qc_read_counts.md) and
[`qc_read_count_dir`](qc_read_count_dir.md).

## Sample table (`meta`)

A data frame with one row per sample. The row order sets the order in
which samples are shown in plots. Required columns:

- sample_id:

  name of the sample's column in `counts`

- label:

  short, unique display label for the sample

- group:

  biological group / condition (preferably a factor, whose levels then
  set the order of groups)

- replicate:

  replicate identifier (character)

Any further columns (e.g. `cell_line`, `day`) can be used for colouring,
faceting or point shapes through the `colour_by`, `facet_by` and
`shape_by` arguments. Use [`qc_check_meta`](qc_check_meta.md) to check a
sample table against a count matrix.

## Count summary (`summary_df`)

A data frame with one row per sample, with at least the columns
`sample_id`, `total_reads` (sequenced reads) and `mapped` (reads
assigned to a library sgRNA), as returned by
[`qc_read_count_summary`](qc_read_count_summary.md) and
[`qc_read_summary_dir`](qc_read_summary_dir.md).

## Font

Plots use the font family set by the option `crisprFlow.font` (default
`"sans"`), e.g. `options(crisprFlow.font = "Arial")`.

## Examples

``` r
## Simulated screen shipped with crisprFlow
head(exampleCounts)
#> # A tibble: 6 × 12
#>   sgRNA  Gene  Plasmid  T0_1  T0_2  T0_3 Untreated_1 Untreated_2 Untreated_3
#>   <chr>  <chr>   <int> <int> <int> <int>       <int>       <int>       <int>
#> 1 AAK1_1 AAK1      227   192   236   201         193         227         412
#> 2 AAK1_2 AAK1      299   251   261   336         220         353         270
#> 3 AAK1_3 AAK1      224   125   146   185         164         194         219
#> 4 AAK1_4 AAK1      548   479   430   703         557         544         376
#> 5 AAK1_5 AAK1      183   201   158   171         221         267         152
#> 6 ABL2_1 ABL2      658   434   481   551         786         724         773
#> # ℹ 3 more variables: Drug_1 <int>, Drug_2 <int>, Drug_3 <int>
exampleSampleMeta
#> # A tibble: 10 × 5
#>    sample_id   label       group     replicate   day
#>    <chr>       <chr>       <fct>     <chr>     <int>
#>  1 Plasmid     Plasmid     Plasmid   1             0
#>  2 T0_1        T0_1        T0        1             2
#>  3 T0_2        T0_2        T0        2             2
#>  4 T0_3        T0_3        T0        3             2
#>  5 Untreated_1 Untreated_1 Untreated 1             9
#>  6 Untreated_2 Untreated_2 Untreated 2             9
#>  7 Untreated_3 Untreated_3 Untreated 3             9
#>  8 Drug_1      Drug_1      Drug      1             9
#>  9 Drug_2      Drug_2      Drug      2             9
#> 10 Drug_3      Drug_3      Drug      3             9
qc_check_meta(exampleSampleMeta, exampleCounts)
```
