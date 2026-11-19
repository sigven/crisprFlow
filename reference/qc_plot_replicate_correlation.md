# Plot sample-to-sample correlations

Heatmap of pairwise Pearson correlations between samples, computed on
log2(CPM + 1) of targeting sgRNAs (see [`qc_log2cpm`](qc_log2cpm.md)).
Replicates should correlate more strongly with each other than with
other groups.

## Usage

``` r
qc_plot_replicate_correlation(
  counts,
  meta,
  ctrl_pattern = "NEG_CTRL",
  limits = NULL,
  show_values = TRUE,
  interactive = TRUE,
  base_size = 10
)
```

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

- limits:

  Optional colour scale limits, `c(min, max)`. Default (NULL) runs from
  the lowest correlation (rounded down) to 1.

- show_values:

  Logical; print correlation values in the tiles (default is TRUE).

- interactive:

  Logical; if TRUE (default), return an interactive plotly object,
  otherwise a ggplot object.

- base_size:

  Base font size of the plot (default is 10).

## Value

A plotly object (`interactive = TRUE`) or a ggplot object.

## Examples

``` r
qc_plot_replicate_correlation(exampleCounts, exampleSampleMeta,
                              interactive = FALSE)

```
