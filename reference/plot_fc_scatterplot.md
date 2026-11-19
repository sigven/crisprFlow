# Plot fold change scatterplot, drug vs. T0, control vs. T0

This function generates a scatterplot comparing fold changes between
drug treatment and control conditions

## Usage

``` r
plot_fc_scatterplot(
  sgRNA_fc_df = NULL,
  x_column = "T4_control_vs_T0_FC",
  y_column = "T4_drug_vs_T0_FC",
  x_lab = "Control vs T0",
  y_lab = "Drug vs T0",
  title = "T4 - sgRNA abundance LFC",
  title_y_right_margin = 12,
  title_x_top_margin = 12,
  ggsci_style = "JCO",
  treatment_name = "FGF401",
  per_sgRNA = TRUE,
  highlight_genes = TRUE,
  plot_fontsize = 13
)
```

## Arguments

- sgRNA_fc_df:

  A data frame containing sgRNA fold change data

- x_column:

  The name of the column to use for the x-axis (default is
  "T4_control_vs_T0_FC")

- y_column:

  The name of the column to use for the y-axis (default is
  "T4_drug_vs_T0_FC")

- x_lab:

  The label for the x-axis (default is "T4 - Control vs T0 - fold change
  (normalized sgRNA count)")

- y_lab:

  The label for the y-axis (default is "T4 - Drug vs T0 - fold change
  (normalized sgRNA count)")

- title:

  The title of the scatterplot (default is "T4 - sgRNA fold change")

- title_y_right_margin:

  Right margin for the y-axis title (default is 12)

- title_x_top_margin:

  Top margin for the x-axis title (default is 12)

- ggsci_style:

  The style for the color palette from the ggsci package (default is
  "JCO")

- treatment_name:

  The name of the treatment used in the experiment (default is "
  FGF401")

- per_sgRNA:

  Logical indicating whether to plot at the sgRNA level (TRUE) or gene
  level (FALSE) (default is TRUE)

- highlight_genes:

  Logical indicating whether to highlight specific genes on the
  scatterplot (default is TRUE)

- plot_fontsize:

  Font size for the plot text elements (default is 11

## Value

A ggplot2 object representing the fold change scatterplot
