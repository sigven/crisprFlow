# Simulated CRISPR screen: sample table

Sample table for [`exampleCounts`](exampleCounts.md), in the layout
expected by the QC functions (see [`qc_inputs`](qc_inputs.md)).

## Usage

``` r
exampleSampleMeta
```

## Format

A tibble with 10 rows (samples) and 5 columns:

- sample_id:

  column name of the sample in `exampleCounts`

- label:

  display label

- group:

  Plasmid, T0, Untreated or Drug (factor)

- replicate:

  replicate number (character)

- day:

  day of sampling (0 = plasmid)

## Source

Simulated, see `data-raw/qc_example_data.R`.
