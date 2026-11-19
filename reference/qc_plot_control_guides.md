# Plot control sgRNA abundance across samples

Heatmap of the abundance (log10 CPM + 1) of each control sgRNA in each
sample. Controls should be phenotypically neutral, so their abundance
should be stable; missing or outlier controls stand out as off-colour
tiles.

## Usage

``` r
qc_plot_control_guides(
  counts,
  meta,
  ctrl_pattern = "NEG_CTRL",
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

- interactive:

  Logical; if TRUE (default), return an interactive plotly object,
  otherwise a ggplot object.

- base_size:

  Base font size of the plot (default is 10).

## Value

A plotly object (`interactive = TRUE`) or a ggplot object; NULL (with a
message) if no sgRNA matches `ctrl_pattern`.

## Examples

``` r
qc_plot_control_guides(exampleCounts, exampleSampleMeta,
                       interactive = FALSE)

```
