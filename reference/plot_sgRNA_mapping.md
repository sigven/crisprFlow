# Plot sgRNA mapping statistics bar plot

This function generates a bar plot showing the mapping statistics of
sgRNA reads

## Usage

``` r
plot_sgRNA_mapping(
  sgRNA_summary,
  ggsci_style = "JCO",
  treatment_name = "FGF401",
  control_name = "DMSO",
  title_x_top_margin = 12,
  title_y_right_margin = 12,
  plot_fontsize = 12
)
```

## Arguments

- sgRNA_summary:

  A data frame containing sgRNA summary statistics including read counts
  and mapping percentages

- ggsci_style:

  The style for the color palette from the ggsci package (default is
  "JCO")

- treatment_name:

  The name of the treatment used in the experiment (default is "
  FGF401")

- control_name:

  The name of the control used in the experiment (default is " DMSO")

- title_x_top_margin:

  Top margin for the x-axis title (default is 12)

- title_y_right_margin:

  Right margin for the y-axis title (default is 12)

- plot_fontsize:

  Font size for the plot text elements (default is 12)

## Value

A ggplot2 object representing the sgRNA mapping statistics bar plot
