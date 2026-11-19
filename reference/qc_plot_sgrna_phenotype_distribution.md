# sgRNA-level phenotype distribution with per-gene strip chart

Produces a two-panel figure. The top panel shows the kernel density of
`fc_col` log2 fold-changes for all real-gene sgRNAs (guide-counter exact
counts). The bottom panel shows jittered dot strips for the top
`n_promoters` and top `n_suppressors` hit genes ranked by Chronos
`tvb_effect`, coloured by hit classification. Guide consistency across
sgRNAs is the key QC signal: a robust hit should shift most or all of
its guides in the same direction.

## Usage

``` r
qc_plot_sgrna_phenotype_distribution(
  sgrna_fc,
  hit_df,
  class_colours,
  fc_col = "Top_vs_Bottom_log2FC",
  n_promoters = 6L,
  n_suppressors = 6L,
  ctrl_pattern = "NEG_CTRL|SAFE_TARGETING",
  restrict_genes = NULL,
  x_limits = c(-4.5, 4.5),
  interactive = FALSE,
  base_size = 11
)
```

## Arguments

- sgrna_fc:

  Data frame of per-sgRNA fold changes, as produced by guide-counter.
  Must contain columns `sgRNA`, `Gene`, and the column named by
  `fc_col`.

- hit_df:

  Gene-level hit table. Must contain columns `Gene`, `ctrl_type`
  (`"Real gene"` for real genes), `classification`, and `tvb_effect`
  (Top-vs-Bottom Chronos effect).

- class_colours:

  Named character vector mapping classification labels to hex colours.

- fc_col:

  Name of the fold-change column in `sgrna_fc` to plot (default
  `"Top_vs_Bottom_log2FC"`).

- n_promoters:

  Number of top invasion promoters to show (default 6).

- n_suppressors:

  Number of top invasion suppressors to show (default 6).

- ctrl_pattern:

  Regex matching control sgRNA gene names to exclude from the density
  and strips (default `"NEG_CTRL|SAFE_TARGETING"`).

- restrict_genes:

  Optional character vector of gene names to restrict candidate hits to
  (e.g. genes confirmed by both Chronos and MAGeCK). When `NULL`
  (default) all classified hit genes are eligible.

- x_limits:

  Numeric vector of length 2 giving x-axis limits (default
  `c(-4.5, 4.5)`).

- interactive:

  Logical. `FALSE` (default) returns a combined patchwork static figure
  (density + strips). `TRUE` returns a plotly strip chart only (with
  hover tooltips showing sgRNA, gene, and LFC).

- base_size:

  Base font size passed to `.theme_qc()` (default 11).

## Value

A patchwork ggplot object (`interactive = FALSE`) or a plotly object
(`interactive = TRUE`).
