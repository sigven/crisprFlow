# Read sgRNA count data from a MAGeCK/guide-counter count file

This function reads sgRNA count data from a specified
MAGeCK/guide-counter count file, filters out non-targeting controls, and
returns a data frame

## Usage

``` r
read_sgRNA_counts(sgRNA_count_file)
```

## Arguments

- sgRNA_count_file:

  The path to the MAGeCK/guide-counter count file

## Value

A data frame containing sgRNA counts with sgRNA IDs as row names
