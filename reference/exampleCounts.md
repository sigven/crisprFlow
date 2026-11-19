# Simulated CRISPR screen: sgRNA counts

Simulated raw sgRNA counts for a small screen with the CHIP1 library
(200 genes with 5 sgRNAs each, plus 50 control sgRNAs), used in the
examples of the QC functions ([`qc_inputs`](qc_inputs.md)). Samples are
a plasmid sample, T0 and day-9 untreated / drug-treated samples, in
triplicate. The data contain deliberate QC issues: some sgRNAs are lost
between plasmid and T0, and sample `Drug_3` is noisy and has a low
mapping rate.

## Usage

``` r
exampleCounts
```

## Format

A tibble with 1,050 rows and 12 columns: `sgRNA`, `Gene`, and one column
of raw counts per sample (see
[`exampleSampleMeta`](exampleSampleMeta.md)).

## Source

Simulated, see `data-raw/qc_example_data.R`.
