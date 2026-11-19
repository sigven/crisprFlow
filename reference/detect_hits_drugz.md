# Run DrugZ for hit detection

This function is a placeholder for running DrugZ for hit detection.

## Usage

``` r
detect_hits_drugz(
  samplesheet_csv = NULL,
  count_table_txt = NULL,
  pseudo_count = 5,
  output_dir = NULL,
  output_prefix = "FGF401",
  drugz_python = "/Users/sigven/miniconda3/envs/drugz/bin/python"
)
```

## Arguments

- samplesheet_csv:

  Path to the sample sheet CSV file

- count_table_txt:

  Path to the count table text file (from MAGeCK/guidecounter)

- pseudo_count:

  Pseudo count to add (default is 5)

- output_dir:

  Directory where output files will be saved

- output_prefix:

  Prefix for output files (default is "FGF401")

- drugz_python:

  Path to the DrugZ python interpreter (Conda environment)

## Value

A list containing DrugZ hit detection results
