# Check that the counted sgRNAs match the intended library

Compares sgRNA IDs in the count matrix with those in the library. A
mismatch (e.g. counts generated against a different library) invalidates
every downstream step, so this is worth running first.

## Usage

``` r
qc_check_library_concordance(counts, library, min_overlap = 0.95)
```

## Arguments

- counts:

  Count matrix: data frame with columns `sgRNA`, `Gene` and one column
  of raw counts per sample (see section 'Count matrix').

- library:

  sgRNA library, data frame with an `sgRNA` column (see
  [`qc_read_library`](qc_read_library.md)).

- min_overlap:

  Minimum fraction of shared sgRNAs - both of the library and of the
  count matrix - for the two to be called concordant (default is 0.95).

## Value

A one-row tibble with the number of sgRNAs in the library, in the count
matrix and in both, the shared fractions, and `concordant` (TRUE/FALSE).

## Examples

``` r
## Not concordant: the example data use only a subset of CHIP1
qc_check_library_concordance(exampleCounts, qc_read_library("CHIP1"))
#> # A tibble: 1 × 6
#>   n_library_guides n_count_guides n_shared frac_library_in_counts
#>              <int>          <int>    <int>                  <dbl>
#> 1            11985           1050     1050                 0.0876
#> # ℹ 2 more variables: frac_counts_in_library <dbl>, concordant <lgl>
```
