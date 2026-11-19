# Plot sgRNA count distributions

Per-sample density of read counts of targeting sgRNAs (log10 scale),
with each sample's median (dashed) and the low-count threshold (dotted).
Narrow, overlapping distributions indicate even library representation.

## Usage

``` r
qc_plot_count_distribution(
  counts,
  meta,
  ctrl_pattern = "NEG_CTRL",
  colour_by = "replicate",
  colours = NULL,
  facet_by = "group",
  low_threshold = 30,
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

- colour_by:

  Name of the column in `meta` used to colour samples (default is
  "replicate").

- colours:

  Optional named character vector of colours, with names matching the
  values of the `colour_by` column (e.g.
  `c(Untreated = "#0073C2", Drug = "#CD534C")`). If NULL (default), the
  RColorBrewer "Dark2" palette is used.

- facet_by:

  Optional name of a column in `meta` used to split the plot into panels
  (default is "group").

- low_threshold:

  Read count drawn as a dotted reference line (default is 30).

- base_size:

  Base font size of the plot.

## Value

A ggplot object.

## Examples

``` r
qc_plot_count_distribution(exampleCounts, exampleSampleMeta)

```
