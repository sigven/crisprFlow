# Function to clean treatment and timepoint labels from a data frame This function separates the 'label' column into 'treatment' and 'timepoint' columns, removes any existing 'treatment' and 'timepoint' columns, and arranges the data frame by timepoint

Function to clean treatment and timepoint labels from a data frame This
function separates the 'label' column into 'treatment' and 'timepoint'
columns, removes any existing 'treatment' and 'timepoint' columns, and
arranges the data frame by timepoint

## Usage

``` r
clean_treatment_timepoint_label(df = NULL, timecourse_drug_screen = TRUE)
```

## Arguments

- df:

  A data frame containing a 'label' column with treatment and timepoint
  information

- timecourse_drug_screen:

  Logical indicating whether the data is from a timecourse drug screen
  (default is TRUE)
