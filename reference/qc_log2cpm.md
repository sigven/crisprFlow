# log2 CPM matrix of targeting sgRNAs

Normalizes counts of targeting sgRNAs (controls removed) to counts per
million and log2-transforms them, log2(CPM + 1).

## Usage

``` r
qc_log2cpm(counts, sample_ids, ctrl_pattern = "NEG_CTRL")
```

## Arguments

- counts:

  Count matrix: data frame with columns `sgRNA`, `Gene` and one column
  of raw counts per sample (see section 'Count matrix').

- sample_ids:

  Sample columns of `counts` to include.

- ctrl_pattern:

  Regular expression matched against the `Gene` column to identify
  control sgRNAs (default is "NEG_CTRL", which matches the controls of
  all libraries shipped with crisprFlow).

## Value

A numeric matrix with one row per targeting sgRNA (row names are sgRNA
IDs) and one column per sample.

## Examples

``` r
m <- qc_log2cpm(exampleCounts, exampleSampleMeta$sample_id)
m[1:3, 1:4]
#>         Plasmid     T0_1     T0_2     T0_3
#> AAK1_1 9.213399 8.991235 9.304064 9.049755
#> AAK1_2 9.610267 9.377150 9.449108 9.789926
#> AAK1_3 9.194238 8.373576 8.612652 8.930320
```
