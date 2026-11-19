# Plot mapped and unmapped reads per sample

Horizontal stacked bars of mapped and unmapped reads per sample, with
the mapping rate printed on each bar. Samples below `min_pct_mapped` are
marked "(low)".

## Usage

``` r
qc_plot_read_mapping(
  summary_df,
  meta,
  colour_by = "group",
  colours = NULL,
  facet_by = NULL,
  min_pct_mapped = 0.65,
  interactive = TRUE,
  base_size = 11
)
```

## Arguments

- summary_df:

  Count summary: data frame with columns `sample_id`, `total_reads` and
  `mapped` (see section 'Count summary').

- meta:

  Sample table: data frame with one row per sample and columns
  `sample_id`, `label`, `group` and `replicate` (see section 'Sample
  table'). Only samples in `meta` are used.

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

- min_pct_mapped:

  Minimum acceptable fraction of mapped reads, from 0 to 1 (default is
  0.65, the minimum used by MAGeCK-VISPR).

- interactive:

  Logical; if TRUE (default), return an interactive plotly object,
  otherwise a ggplot object.

- base_size:

  Base font size of the plot.

## Value

A plotly object (`interactive = TRUE`) or a ggplot object.

## Examples

``` r
qc_plot_read_mapping(exampleCountSummary, exampleSampleMeta,
                     interactive = FALSE)

```
