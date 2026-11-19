# Plot the fraction of reads mapping to each negative-control group

For each control type defined in `ctrl_patterns`, computes the
percentage of total raw reads per sample that map to that control group,
and plots the values as a dot plot faceted by control type. Stable
proportions across samples confirm that control guides are evenly
represented and that library composition has not shifted.

## Usage

``` r
qc_plot_control_read_proportion(
  counts,
  meta,
  ctrl_patterns = c(`Negative controls` = "NEG_CTRL"),
  colour_by = "group",
  colours = NULL,
  y_limits = NULL,
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

- ctrl_patterns:

  Named character vector mapping display labels to regex patterns used
  to identify control guides in the `Gene` column of `counts`. For
  example:
  `c("Non-targeting" = "NONE_TARGETING", "Safe-targeting" = "SAFE_TARGETING")`.
  Each pattern is passed to
  [`str_detect`](https://stringr.tidyverse.org/reference/str_detect.html).

- colour_by:

  Column in `meta` used for point colour (default `"group"`).

- colours:

  Optional named character vector of colours for the `colour_by` levels.

- y_limits:

  Optional numeric vector of length 2 (`c(min, max)`) passed to
  `ylim()`. If `NULL` (default), limits are set automatically.

- interactive:

  Logical; if TRUE (default), return an interactive plotly object,
  otherwise a ggplot object.

- base_size:

  Base font size of the plot.

## Value

A plotly object (`interactive = TRUE`) or a ggplot object; NULL (with a
message) if no guides match any pattern in `ctrl_patterns`.

## Examples

``` r
qc_plot_control_read_proportion(
  exampleCounts, exampleSampleMeta,
  ctrl_patterns = c("Negative controls" = "NEG_CTRL"),
  interactive = FALSE
)

```
