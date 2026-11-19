# Count sgRNA from a single FASTQ file using MAGeCK

This function counts sgRNA sequences from a single FASTQ file using the
MAGeCK tool (`mageck count`). It is the MAGeCK counterpart of
`count_sgRNA_guidecounter_fastq`, taking one FASTQ file and a sample
label, and producing output files named by the sample label.

## Usage

``` r
count_sgRNA_mageck_fastq(
  fname_fastq = NULL,
  sample_label = NULL,
  sgRNA_library = "DTKP",
  mageck_bin = "/Users/sigven/miniconda3/bin/mageck",
  norm_method = "none",
  count_N = FALSE,
  output_dir = NULL,
  overwrite = FALSE
)
```

## Arguments

- fname_fastq:

  A character string specifying the path to the FASTQ file to be
  counted.

- sample_label:

  A character string specifying the label for the sample being counted.

- sgRNA_library:

  A character string specifying the sgRNA library to use. Must be one of
  'ACOC','DTKP','GEEX','MEPR','PROT','TMMO','CHIP1'. Default is "DTKP".

- mageck_bin:

  A character string specifying the path to the MAGeCK binary.

- norm_method:

  Method for normalization ('none','median','total') (default is
  "none"). Only affects MAGeCK's normalized count file, which is not
  retained - the returned counts are raw counts.

- count_N:

  A logical indicating whether to count 'N' bases in sgRNA sequences
  (default is FALSE).

- output_dir:

  A character string specifying the directory where the count results
  will be saved.

- overwrite:

  A logical indicating whether to overwrite existing output files
  (default is FALSE).

## Value

A list containing two data frames: 'stats' with summary statistics and
'counts' with sgRNA counts.
