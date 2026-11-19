# Plot sgRNA Gini index bar plot

This function generates a bar plot of the Gini index for sgRNA
distributions

## Usage

``` r
plot_sgRNA_gini_index(sgRNA_summary, ggsci_style = "JCO", plot_fontsize = 11)
```

## Arguments

- sgRNA_summary:

  A data frame containing sgRNA summary statistics including Gini index

- ggsci_style:

  The style for the color palette from the ggsci package (default is
  "JCO")

- plot_fontsize:

  Font size for the plot text elements (default is 11

## Value

A ggplot2 object representing the Gini index bar plot
