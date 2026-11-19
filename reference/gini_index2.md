# Compute the Gini coefficient/index for a numeric vector x.

Calculate Gini index - MAGeCK alike (same results as MAGeCK, no count
normalization - only log2 transform). The Gini index is a measure of
statistical dispersion that represents the inequality among values of a
frequency distribution. It ranges from 0 (perfect equality) to 1
(perfect inequality).

## Usage

``` r
gini_index2(
  x,
  normalize_counts = FALSE,
  log_transform = TRUE,
  na.rm = TRUE,
  allow.negative = FALSE,
  digits = 4
)
```

## Arguments

- x:

  numeric vector (sgRNA counts). Non-negative values expected.

- normalize_counts:

  A logical indicating whether to normalize

- log_transform:

  A logical indicating whether to log2-transform the values in x
  (default is TRUE).

- na.rm:

  logical; remove NA pairs (x and w) before computing. Default TRUE.

- allow.negative:

  logical; if FALSE (default) stops when x has negatives. If TRUE, it
  will shift x by adding -min(x) (so smallest becomes 0). (Shifting
  changes the Gini interpretation — use with care.)

- digits:

  numeric; round to the specified number of decimal places.

## Value

numeric scalar in \[0,1\] (or NaN if undefined)
