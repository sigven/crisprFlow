# Within-group replicate scatter plots

For each group in `meta`, plots all pairwise combinations of replicates
as log2(CPM + 1) scatter plots, with a diagonal reference line and
Pearson r annotated per panel. The result is a faceted grid with groups
as rows and replicate pairs as columns — suitable for QC reports and
supplementary figures.

## Usage

``` r
qc_plot_replicate_scatter(
  counts,
  meta,
  ctrl_pattern = "NEG_CTRL",
  by = "group",
  point_colour = "#444444",
  point_alpha = 0.25,
  point_size = 0.6,
  axis_limits = NULL,
  interactive = TRUE,
  base_size = 11
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

- by:

  Column in `meta` defining replicate groups (default `"group"`). Groups
  with only one sample are skipped.

- point_colour:

  Colour for scatter points (default `"#444444"`).

- point_alpha:

  Alpha for points (default `0.25`).

- point_size:

  Point size (default `0.6`).

- axis_limits:

  Optional numeric vector `c(min, max)` applied to both axes. If `NULL`
  (default) the range is set automatically from the data.

- interactive:

  Logical; if TRUE (default), return an interactive plotly object,
  otherwise a ggplot object.

- base_size:

  Base font size of the plot.

## Value

A ggplot object (faceted grid). Plotly is not supported for this plot
due to the number of panels; `interactive` is accepted for API
consistency but always returns ggplot.

## Details

Only targeting sgRNAs (guides not matching `ctrl_pattern`) are included
in both the scatter points and the correlation calculation, matching
[`qc_replicate_correlation_table`](qc_replicate_correlation_table.md).

## Examples

``` r
qc_plot_replicate_scatter(exampleCounts, exampleSampleMeta,
                          interactive = FALSE)

```
