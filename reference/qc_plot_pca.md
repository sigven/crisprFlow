# Plot a PCA of samples

Principal component analysis of log2(CPM + 1) of targeting sgRNAs (see
[`qc_log2cpm`](qc_log2cpm.md)); sgRNAs without variance are removed.
Replicates should cluster together.

## Usage

``` r
qc_plot_pca(
  counts,
  meta,
  ctrl_pattern = "NEG_CTRL",
  colour_by = "group",
  colours = NULL,
  shape_by = NULL,
  n_top = NULL,
  pcs = c(1, 2),
  label_points = TRUE,
  interactive = TRUE,
  base_size = 12
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

- colour_by:

  Name of the column in `meta` used to colour samples.

- colours:

  Optional named character vector of colours, with names matching the
  values of the `colour_by` column (e.g.
  `c(Untreated = "#0073C2", Drug = "#CD534C")`). If NULL (default), the
  RColorBrewer "Dark2" palette is used.

- shape_by:

  Optional name of a column in `meta` used for point shapes (default is
  NULL).

- n_top:

  Optional number of most variable sgRNAs to use (default is NULL, all
  sgRNAs).

- pcs:

  The two principal components to plot (default is `c(1, 2)`).

- label_points:

  Logical; label points with `meta$label` (default is TRUE; static plots
  only - interactive plots show labels on hover).

- interactive:

  Logical; if TRUE (default), return an interactive plotly object,
  otherwise a ggplot object.

- base_size:

  Base font size of the plot (default is 12).

## Value

A plotly object (`interactive = TRUE`) or a ggplot object.

## Examples

``` r
qc_plot_pca(exampleCounts, exampleSampleMeta, interactive = FALSE)

```
