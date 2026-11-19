# Estimate log2 fold change of sgRNA counts between arbitrary conditions

Normalizes raw sgRNA counts to CPM and computes log2 fold changes for
any set of condition pairs specified via a comparisons table. Replicate
columns are resolved through the samplesheet and collapsed either by
averaging CPM before computing FC, or by computing per-replicate FC and
then averaging.

## Usage

``` r
estimate_fold_change(
  counts = NULL,
  samplesheet = NULL,
  comparisons = NULL,
  replicate_method = c("mean_cpm_then_fc", "per_replicate_fc"),
  pseudo_count = 1,
  output_dir = NULL,
  output_prefix = "crisprFlow",
  overwrite = FALSE,
  control_pattern = "NEG_CTRL"
)
```

## Arguments

- counts:

  A data frame with columns `sgRNA`, `Gene`, and one column per sample
  whose names match those derived from the samplesheet.

- samplesheet:

  A data frame with at least columns `sample_label` and `condition`. An
  optional `replicate` column drives the mapping from condition to count
  columns.

- comparisons:

  A data frame with columns `numerator`, `denominator`, and `label`
  defining which condition pairs to compare and the prefix used for the
  output FC columns.

- replicate_method:

  How to handle multiple replicates per condition. `"mean_cpm_then_fc"`
  (default) averages CPM across replicates first, then computes a single
  log2FC. `"per_replicate_fc"` computes log2FC for each replicate pair
  (matched when counts are equal, all-pairwise otherwise) and reports
  the mean; a `_log2FC_sd` column is also added when more than one pair
  exists.

- pseudo_count:

  Pseudo-count added before taking log2 ratios to avoid log(0) (default
  1).

- output_dir:

  Directory where output TSV files are written.

- output_prefix:

  Prefix for the output file names (default `"crisprFlow"`).

- overwrite:

  If `FALSE` (default) and output files already exist, they are read and
  returned without recomputing.

- control_pattern:

  Pattern identifying negative-control sgRNAs; not used in the
  computation but retained for consistency with the rest of the package
  API.

## Value

A named list with elements `per_sgRNA` and `per_gene`. Each is a data
frame containing mean CPM columns for every referenced condition, log2FC
columns (and optionally SD columns) for every comparison, and — for
`per_gene` — values averaged across the sgRNAs targeting each gene.

## Details

The samplesheet must contain at least `sample_label` and `condition`
columns. If a `replicate` column is present, count column names are
expected to follow the pattern `<sample_label>_<replicate>` (e.g.
`Bottom_B1`); otherwise just `<sample_label>` is used.

When only a single replicate exists for either condition in a
comparison, `replicate_method` is ignored and the single-column CPM
value is used directly.
