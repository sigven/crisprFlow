# Simulated CRISPR screen: count summary

Per-sample count summary for [`exampleCounts`](exampleCounts.md), in the
layout returned by [`qc_read_count_summary`](qc_read_count_summary.md).

## Usage

``` r
exampleCountSummary
```

## Format

A tibble with 10 rows (samples) and 7 columns:

- sample_id:

  column name of the sample in `exampleCounts`

- total_reads:

  number of sequenced reads

- mapped:

  number of reads assigned to a library sgRNA

- percentage:

  fraction of reads mapped

- total_sgrnas:

  number of sgRNAs in the library

- zero_counts:

  number of sgRNAs with zero reads

- gini_index:

  Gini index of the sgRNA count distribution

## Source

Simulated, see `data-raw/qc_example_data.R`.
