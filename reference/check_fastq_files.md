# Check for expected FASTQ files in a directory

This function checks if all expected FASTQ files are present in the
specified directory based on provided timepoints and treatments.

## Usage

``` r
check_fastq_files(
  timepoints = c("T0", "T4", "T10", "T14"),
  treatments = c("DMSO", "FGF401"),
  fastq_dir = NULL
)
```

## Arguments

- timepoints:

  A character vector of timepoints to check (default is
  c("T0","T4","T10","T14")).

- treatments:

  A character vector of treatments to check (default is
  c("DMSO","FGF401")).

- fastq_dir:

  A character string specifying the directory containing the FASTQ
  files.

## Value

None. The function stops execution if any expected FASTQ files are
missing.
