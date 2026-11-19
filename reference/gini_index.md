# Compute the Gini coefficient/index for a numeric vector x.

This function calculates the Gini index for a numeric vector. The Gini
index is a measure of statistical dispersion that represents the
inequality among values of a frequency distribution. It ranges from 0
(perfect equality) to 1 (perfect inequality).

## Usage

``` r
gini_index(x, na.rm = TRUE, normalize_counts = TRUE, log_transform = TRUE)
```

## Arguments

- x:

  A numeric vector for which to calculate the Gini index.

- na.rm:

  A logical indicating whether to remove NA values from the input vector
  (default is TRUE).

- normalize_counts:

  A logical indicating whether to normalize the counts to counts per
  million and log2-transform them (default is TRUE).

- log_transform:

  A logical indicating whether to log2-transform

## Value

A numeric value representing the Gini index of the input vector.
