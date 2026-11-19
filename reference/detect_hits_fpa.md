# Run FPA for hit detection

This function runs the FPA (Fitness Perturbation Analysis) method for
hit detection using a provided count table. It processes the input data,
executes the FPA analysis with specified parameters, and returns the
results.

## Usage

``` r
detect_hits_fpa(
  data = NULL,
  output_dir = NULL,
  output_prefix = "fpa_hits",
  fpa_t0_cutoff = 0,
  fpa_differential = FALSE,
  fpa_control_condition = "DMSO",
  fpa_smoothing = "Laplace",
  fpa_fitness_computation = "Regression",
  fpa_statistics_method = "Density",
  fpa_statistic = "Sum",
  fpa_iterations = 10000,
  control_pattern = "NEG_CTRL"
)
```

## Arguments

- data:

  Data frame containing sgRNA count data (from MAGeCK count or
  guidecounter)

- output_dir:

  Directory where output files will be saved

- output_prefix:

  Prefix for output files (default is "crisprFlowRun")

- fpa_t0_cutoff:

  Minimum count threshold for T0 (default is 0)

- fpa_differential:

  Boolean indicating if differential analysis should be performed
  (default is FALSE)

- fpa_control_condition:

  Condition to be used as control (default is "DMSO")

- fpa_smoothing:

  Smoothing method to use ("Gradient", "Laplace", or "Both")

- fpa_fitness_computation:

  Fitness computation method to use ("Regression", "SVD", or "Both")

- fpa_statistics_method:

  Statistics method to use ("Density" or "Rank")

- fpa_statistic:

  Statistic to use ("Median", "Sum", or "Both")

- fpa_iterations:

  Number of iterations for the analysis (default is 10000)

- control_pattern:

  Pattern to identify control genes (default is "NEG_CTRL")

## Value

A list containing the FPA analysis results
