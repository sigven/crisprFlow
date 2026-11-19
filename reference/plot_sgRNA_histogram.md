# Plot sgRNA count distribution histogram

This function generates histograms of sgRNA count distributions for
different treatments and timepoints

## Usage

``` r
plot_sgRNA_histogram(
  counts = NULL,
  ggsci_style = "JCO",
  treatment_name = "FGF401",
  control_name = "DMSO",
  num_bins = 100,
  title_x_top_margin = 12,
  title_y_right_margin = 12,
  plot_fontsize = 11
)
```

## Arguments

- counts:

  A data frame containing sgRNA counts with sgRNA IDs as row names

- ggsci_style:

  The style for the color palette from the ggsci package (default is
  "JCO")

- treatment_name:

  The name of the treatment used in the experiment (default is "
  FGF401")

- control_name:

  The name of the control used in the experiment (default is " DMSO")

- num_bins:

  The number of bins to use in the histogram (default is 100)

- title_x_top_margin:

  Top margin for the x-axis title (default is 12)

- title_y_right_margin:

  Right margin for the y-axis title (default is

- plot_fontsize:

  Font size for the plot text elements (default is 11

## Value

A list of ggplot2 objects representing the sgRNA count distribution
histograms
