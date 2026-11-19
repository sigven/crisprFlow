# Read an sgRNA library

Reads an sgRNA library, either one shipped with crisprFlow (by name) or
a library file (MAGeCK layout: sgRNA ID, sequence, gene; comma- or
tab-separated, header optional).

## Usage

``` r
qc_read_library(sgRNA_library = NULL, library_file = NULL)
```

## Arguments

- sgRNA_library:

  Name of a library shipped with crisprFlow
  ('ACOC','DTKP','GEEX','MEPR','PROT','TMMO','CHIP1'). Ignored if
  `library_file` is given.

- library_file:

  Optional path to a library file.

## Value

A data frame with columns `sgRNA`, `seq` and `Gene`.

## Examples

``` r
lib <- qc_read_library("CHIP1")
head(lib)
#> # A tibble: 6 × 3
#>   sgRNA  seq                  Gene 
#>   <chr>  <chr>                <chr>
#> 1 A1CF_1 GCAGCATCCCAACCAGGTGG A1CF 
#> 2 A1CF_2 ATGACTCTCATACTCCACGA A1CF 
#> 3 A1CF_3 GAGTCATCGAGCAGCTGCCA A1CF 
#> 4 A1CF_4 AGTCACCCTAGCAAAACCAG A1CF 
#> 5 A1CF_5 GATCCCACCACAACCTACCT A1CF 
#> 6 AAK1_1 GAAATGGTCAACCTGTACAG AAK1 
```
