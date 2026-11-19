# Read and combine per-sample count summaries

Reads all per-sample count summary files in a directory (e.g. the
`<sample_label>.countsummary.txt` files written by
[`count_sgRNA_guidecounter_fastq`](count_sgRNA_guidecounter_fastq.md))
and binds them into one table.

## Usage

``` r
qc_read_summary_dir(count_dir, pattern = "\\.countsummary\\.txt$")
```

## Arguments

- count_dir:

  Directory with the per-sample summary files.

- pattern:

  Regular expression matching the summary file names (default matches
  files ending in `.countsummary.txt`).

## Value

A count summary, see
[`qc_read_count_summary`](qc_read_count_summary.md).
