# Count sgRNA from FASTQ files using MAGECK

This function counts sgRNA sequences from FASTQ files using the MAGECK
tool. It requires a sample sheet CSV file containing sample information,
the sgRNA library to use, and parameters for output directory and file
naming.

## Usage

``` r
count_sgRNA_mageck(
  samplesheet_csv = NULL,
  sgRNA_library = "DTKP",
  mageck_bin = "/Users/sigven/miniconda3/bin/mageck",
  norm_method = "median",
  count_N = FALSE,
  output_dir = NULL,
  output_prefix = "FGFR4",
  overwrite = FALSE
)
```

## Arguments

- samplesheet_csv:

  A character string specifying the path to a CSV file containing sample
  information with mandatory columns: fastq, sample_label, condition.

- sgRNA_library:

  A character string specifying the sgRNA library to use. Must be one of
  'ACOC','DTKP','GEEX','MEPR','PROT','TMMO','CHIP1'. Default is "DTKP".

- mageck_bin:

  A character string specifying the path to the MAGECK binary.

- norm_method:

  Method for normalization ('none','median','total') (default is
  "median").

- count_N:

  A logical indicating whether to count 'N' bases in sgRNA sequences
  (default is FALSE).

- output_dir:

  A character string specifying the directory where the count results
  will be saved.

- output_prefix:

  A character string to prefix output files. Default is "FGFR4".

- overwrite:

  A logical indicating whether to overwrite existing output files
  (default is FALSE).

## Value

None. The function performs counting and saves results to files.
