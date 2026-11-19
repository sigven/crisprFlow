# Read a merged count matrix

Reads a tab-separated count matrix with columns `sgRNA`, `Gene` and one
column per sample (e.g. `<prefix>.count.txt` from MAGeCK count, or a
merged guide-counter count file).

## Usage

``` r
qc_read_counts(counts_file)
```

## Arguments

- counts_file:

  Path to the tab-separated count file.

## Value

A count matrix (see [`qc_inputs`](qc_inputs.md)).
