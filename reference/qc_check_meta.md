# Check a sample table (and its match to a count matrix)

Checks that the sample table has the required columns, unique sample IDs
and labels, and - if a count matrix is given - that every sample in the
sample table has a count column. Called by all QC functions that take
`meta`, so problems are reported up front with a clear message.

## Usage

``` r
qc_check_meta(meta, counts = NULL)
```

## Arguments

- meta:

  Sample table: data frame with one row per sample and columns
  `sample_id`, `label`, `group` and `replicate` (see section 'Sample
  table'). Only samples in `meta` are used.

- counts:

  Optional count matrix to check `meta` against (see
  [`qc_inputs`](qc_inputs.md)).

## Value

`meta`, invisibly, if all checks pass (otherwise an error).

## Examples

``` r
qc_check_meta(exampleSampleMeta, exampleCounts)
```
