# Count sgRNA from FASTQ files using guide-counter

This function counts sgRNA sequences from FASTQ files using the
guide-counter tool. It requires the path to the guide-counter binary,
the directory containing the FASTQ files, and parameters for timepoints,
treatments and sgRNA library.

## Usage

``` r
count_sgRNA_guidecounter(
  samplesheet_csv = NULL,
  sgRNA_library = "DTKP",
  min_sgRNA_length = 17,
  max_sgRNA_length = 25,
  guide_counter_bin = "/Users/sigven/miniconda3/bin/guide-counter",
  count_exact = TRUE,
  offset_sample_size = "100000",
  offset_min_fraction = 0.0025,
  output_dir = NULL,
  output_prefix = "FGFR4",
  overwrite = FALSE
)
```

## Arguments

- samplesheet_csv:

  A character string specifying the path to a CSV file

- sgRNA_library:

  A character string specifying the sgRNA library to use. Must be one of
  'ACOC','DTKP','GEEX','MEPR','PROT','TMMO','CHIP1'. Default is "DTKP".

- min_sgRNA_length:

  An integer specifying the minimum sgRNA length to consider. Default is
  17.

- max_sgRNA_length:

  An integer specifying the maximum sgRNA length to consider (default is
  25).

- guide_counter_bin:

  A character string specifying the path to the guide-counter binary.

- count_exact:

  A logical indicating whether to use exact matching in guide-counter
  (default is TRUE).

- offset_sample_size:

  An integer specifying the number of reads to be examined for offset
  estimation (default is 100K).

- offset_min_fraction:

  A numeric value specifying the minimum fraction of reads supporting an
  offset to consider it valid (default is 0.0025).

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
