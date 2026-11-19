# Within-group replicate correlations

Summarises pairwise Pearson correlations (log2 CPM of targeting sgRNAs)
between replicates of the same group. Groups with a single sample are
left out.

## Usage

``` r
qc_replicate_correlation_table(
  counts,
  meta,
  ctrl_pattern = "NEG_CTRL",
  by = "group"
)
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

- by:

  Name of the column in `meta` defining replicate groups (default is
  "group").

## Value

A tibble with one row per group: the group, `n_reps` (number of
replicates), `min_r` and `mean_r` (lowest and mean pairwise
correlation).

## Examples

``` r
qc_replicate_correlation_table(exampleCounts, exampleSampleMeta)
#> # A tibble: 3 × 4
#>   group     n_reps min_r mean_r
#>   <fct>      <int> <dbl>  <dbl>
#> 1 T0             3 0.919  0.920
#> 2 Untreated      3 0.939  0.941
#> 3 Drug           3 0.830  0.868
```
