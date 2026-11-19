# Run MAGeCK MLE for hit detection

This function runs MAGeCK MLE for hit detection using a provided sample
sheet and count table. It constructs a design matrix based on the sample
information and executes the MAGeCK MLE command with specified
parameters. The results are saved in the specified output directory.

## Usage

``` r
detect_hits_mageck_mle(
  samplesheet_csv = NULL,
  count_table_txt = NULL,
  output_dir = NULL,
  norm_method = c("median", "total", "control", "none"),
  mageck_bin = "/Users/sigven/miniconda3/bin/mageck",
  threads = 1,
  permutation_round = 4,
  late_time_start = 10,
  adjust_method = "fdr",
  output_prefix = "mageck_mle_hits",
  control_pattern = "NEG_CTRL"
)
```

## Arguments

- samplesheet_csv:

  Path to the sample sheet CSV file

- count_table_txt:

  Path to the count table text file (from MAGeCK count or guidecounter)

- output_dir:

  Directory where output files will be saved

- norm_method:

  Normalization method to use ("median", "total", "control" or "none";
  default is "median"). With "control", counts are normalized to the
  control genes in the count table (genes matching `control_pattern`)

- mageck_bin:

  Path to the MAGeCK binary

- threads:

  Number of threads to use (default is 1)

- permutation_round:

  Number of permutation rounds (default is 2, recommended 10 (longer
  time))

- late_time_start:

  Time point to start considering as "late" (default is 10)

- adjust_method:

  Method for p-value adjustment (default is "fdr")

- output_prefix:

  Prefix for output files (default is "crisprFlowRun")

- control_pattern:

  Pattern to identify control genes (default is "NEG_CTRL")

## Value

A data frame containing the gene summary results from MAGeCK MLE
