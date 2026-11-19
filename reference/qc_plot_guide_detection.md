# Plot detected and lost sgRNAs per stage

Stacked bars of detected vs. absent/low targeting sgRNAs at each stage.

## Usage

``` r
qc_plot_guide_detection(rep_data, base_size = 12)
```

## Arguments

- rep_data:

  Output of [`qc_guide_representation`](qc_guide_representation.md).

- base_size:

  Base font size of the plot (default is 12).

## Value

A ggplot object.

## Examples

``` r
rep_data <- qc_guide_representation(
  exampleCounts,
  stages = list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3")))
qc_plot_guide_detection(rep_data)

```
