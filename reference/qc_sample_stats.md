# Per-sample sgRNA count statistics

Computes per-sample metrics from the count matrix alone, complementing
the count summary (total reads, mapping rate): number of zero and
low-count targeting sgRNAs, median count and the share of reads on
control sgRNAs.

## Usage

``` r
qc_sample_stats(counts, meta, ctrl_pattern = "NEG_CTRL", low_threshold = 30)
```

## Arguments

- counts:

  Count matrix: data frame with columns `sgRNA`, `Gene` and one column
  of raw counts per sample (see section 'Count matrix').

- meta:

  Sample table: data frame with one row per sample and columns
  `sample_id`, `label`, `group` and `replicate` (see section 'Sample
  table'). Only samples in `meta` are used.

- ctrl_pattern:

  Regular expression matched against the `Gene` column to identify
  control sgRNAs (default is "NEG_CTRL", which matches the controls of
  all libraries shipped with crisprFlow).

- low_threshold:

  Targeting sgRNAs with fewer reads than this are counted as low-count
  (default is 30).

## Value

A tibble with one row per sample: `sample_id`, `assigned` (total
counts), `n_zero` and `n_low` (number of targeting sgRNAs with zero /
fewer than `low_threshold` reads), `median_count` (median count of
targeting sgRNAs) and `pct_ctrl` (percent of counts on control sgRNAs).

## Examples

``` r
qc_sample_stats(exampleCounts, exampleSampleMeta)
#> # A tibble: 10 × 6
#>    sample_id   assigned n_zero n_low median_count pct_ctrl
#>    <chr>          <int>  <int> <int>        <dbl>    <dbl>
#>  1 Plasmid       402626      0     0         344.     4.86
#>  2 T0_1          398709      0    36         333      5.19
#>  3 T0_2          393451      0    36         332.     4.96
#>  4 T0_3          401478      0    36         338.     5.35
#>  5 Untreated_1   396544      1    46         337      5.27
#>  6 Untreated_2   400334      0    43         336.     5.43
#>  7 Untreated_3   400652      0    43         341      5.35
#>  8 Drug_1        392594      2    53         320      5.02
#>  9 Drug_2        402617      0    48         319      5.48
#> 10 Drug_3        395738      2    71         244      5.22
```
