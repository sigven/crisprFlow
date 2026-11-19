# Read sgRNA summary statistics from a MAGeCK/guide-counter summary file

This function reads sgRNA summary statistics from a specified
MAGeCK/guide-counter summary file, cleans the column names, and returns
a data frame

## Usage

``` r
read_sgRNA_summary(
  sgRNA_count_summary_file = NULL,
  timecourse_drug_screen = TRUE
)
```

## Arguments

- sgRNA_count_summary_file:

  The path to the MAGeCK/guide-counter count summary file

- timecourse_drug_screen:

  logical indicating whether the data is from a timecourse drug screen
  (default is TRUE)

## Value

A data frame containing cleaned sgRNA summary statistics
