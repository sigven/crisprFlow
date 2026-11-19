# Read sgRNA count data from a MAGeCK/guide-counter count file and separate targeted sgRNAs and controls

Read sgRNA count data from a MAGeCK/guide-counter count file and
separate targeted sgRNAs and controls

## Usage

``` r
read_sgRNA_counts2(sgRNA_count_file = NULL, normalize = T)
```

## Arguments

- sgRNA_count_file:

  The path to the MAGeCK/guide-counter count file

- normalize:

  Logical indicating whether to normalize sgRNA counts to counts per
  million (default is TRUE)

## Value

A data frame containing sgRNA counts with sgRNA IDs as row names
