# Plot PCA of sgRNA count data

This function performs PCA on sgRNA count data and generates a scatter
plot of the first two principal components.

## Usage

``` r
plot_sgRNA_pca(
  title = "Title",
  sgRNA_counts,
  ggsci_style = "JCO",
  treatment_name = "FGF401",
  plot_fontsize = 11
)
```

## Arguments

- title:

  The title of the PCA plot (default is "Title")

- sgRNA_counts:

  A data frame containing sgRNA counts with sgRNA IDs as row names

- ggsci_style:

  The style for the color palette from the ggsci package (default is
  "JCO")

- treatment_name:

  The name of the treatment used in the experiment (default is "
  FGF401")

- plot_fontsize:

  Font size for the plot text elements (default is 11

## Value

A ggplot2 object representing the PCA scatter plot
