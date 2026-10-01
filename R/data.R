#' Simulated CRISPR screen: sgRNA counts
#'
#' Simulated raw sgRNA counts for a small screen with the CHIP1 library
#' (200 genes with 5 sgRNAs each, plus 50 control sgRNAs), used in the
#' examples of the QC functions (\code{\link{qc_inputs}}). Samples are a
#' plasmid sample, T0 and day-9 untreated / drug-treated samples, in
#' triplicate. The data contain deliberate QC issues: some sgRNAs are lost
#' between plasmid and T0, and sample \code{Drug_3} is noisy and has a low
#' mapping rate.
#'
#' @format A tibble with 1,050 rows and 12 columns: \code{sgRNA},
#' \code{Gene}, and one column of raw counts per sample (see
#' \code{\link{exampleSampleMeta}}).
#' @source Simulated, see \code{data-raw/qc_example_data.R}.
"exampleCounts"

#' Simulated CRISPR screen: sample table
#'
#' Sample table for \code{\link{exampleCounts}}, in the layout expected by
#' the QC functions (see \code{\link{qc_inputs}}).
#'
#' @format A tibble with 10 rows (samples) and 5 columns:
#' \describe{
#'   \item{sample_id}{column name of the sample in \code{exampleCounts}}
#'   \item{label}{display label}
#'   \item{group}{Plasmid, T0, Untreated or Drug (factor)}
#'   \item{replicate}{replicate number (character)}
#'   \item{day}{day of sampling (0 = plasmid)}
#' }
#' @source Simulated, see \code{data-raw/qc_example_data.R}.
"exampleSampleMeta"

#' Simulated CRISPR screen: count summary
#'
#' Per-sample count summary for \code{\link{exampleCounts}}, in the layout
#' returned by \code{\link{qc_read_count_summary}}.
#'
#' @format A tibble with 10 rows (samples) and 7 columns:
#' \describe{
#'   \item{sample_id}{column name of the sample in \code{exampleCounts}}
#'   \item{total_reads}{number of sequenced reads}
#'   \item{mapped}{number of reads assigned to a library sgRNA}
#'   \item{percentage}{fraction of reads mapped}
#'   \item{total_sgrnas}{number of sgRNAs in the library}
#'   \item{zero_counts}{number of sgRNAs with zero reads}
#'   \item{gini_index}{Gini index of the sgRNA count distribution}
#' }
#' @source Simulated, see \code{data-raw/qc_example_data.R}.
"exampleCountSummary"
