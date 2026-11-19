# Read a count summary file

Reads a MAGeCK / guide-counter count summary (e.g.
`<prefix>.countsummary.txt`) and harmonises the column names, so that it
can be used as `summary_df` in the QC functions.

## Usage

``` r
qc_read_count_summary(summary_file)
```

## Arguments

- summary_file:

  Path to the tab-separated count summary file. It must contain the
  columns `Label`, `Reads` and `Mapped`.

## Value

A count summary with one row per sample and the columns `sample_id`
(from `Label`), `total_reads`, `mapped`, and - when present -
`percentage`, `total_sgrnas`, `zero_counts` and `gini_index`.
