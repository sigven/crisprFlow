# Read and merge per-sample count files

Reads all per-sample count files in a directory (e.g. the
`<sample_label>.count.txt` files written by
[`count_sgRNA_guidecounter_fastq`](count_sgRNA_guidecounter_fastq.md))
and merges them into one count matrix. sgRNAs missing from a sample get
a count of zero.

## Usage

``` r
qc_read_count_dir(count_dir, pattern = "\\.count\\.txt$")
```

## Arguments

- count_dir:

  Directory with the per-sample count files.

- pattern:

  Regular expression matching the count file names (default matches
  files ending in `.count.txt`).

## Value

A count matrix (see [`qc_inputs`](qc_inputs.md)), with one column per
file.
