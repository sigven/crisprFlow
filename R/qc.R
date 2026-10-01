#' Inputs shared by the QC functions
#'
#' The \code{qc_*} functions check the quality of pooled CRISPR screen count
#' data. Rather than parsing sample names, they are driven by two tables:
#' a count matrix (\code{counts}) and a sample table (\code{meta}). The same
#' functions therefore work for screens with different designs.
#'
#' @section Count matrix (\code{counts}):
#' A data frame with one row per sgRNA and the columns
#' \describe{
#'   \item{sgRNA}{sgRNA identifier}
#'   \item{Gene}{gene targeted by the sgRNA (control sgRNAs carry a gene name
#'   matching \code{ctrl_pattern}, e.g. \code{NONE_TARGETING_NEG_CTRL} or
#'   \code{CTRLA_NEG_CTRL})}
#'   \item{<sample_id>}{one column of raw (non-normalized) read counts per
#'   sample}
#' }
#' This is the layout written by MAGeCK and guide-counter, see
#' \code{\link{qc_read_counts}} and \code{\link{qc_read_count_dir}}.
#'
#' @section Sample table (\code{meta}):
#' A data frame with one row per sample. The row order sets the order in
#' which samples are shown in plots. Required columns:
#' \describe{
#'   \item{sample_id}{name of the sample's column in \code{counts}}
#'   \item{label}{short, unique display label for the sample}
#'   \item{group}{biological group / condition (preferably a factor, whose
#'   levels then set the order of groups)}
#'   \item{replicate}{replicate identifier (character)}
#' }
#' Any further columns (e.g. \code{cell_line}, \code{day}) can be used for
#' colouring, faceting or point shapes through the \code{colour_by},
#' \code{facet_by} and \code{shape_by} arguments. Use
#' \code{\link{qc_check_meta}} to check a sample table against a count
#' matrix.
#'
#' @section Count summary (\code{summary_df}):
#' A data frame with one row per sample, with at least the columns
#' \code{sample_id}, \code{total_reads} (sequenced reads) and \code{mapped}
#' (reads assigned to a library sgRNA), as returned by
#' \code{\link{qc_read_count_summary}} and \code{\link{qc_read_summary_dir}}.
#'
#' @section Font:
#' Plots use the font family set by the option \code{crisprFlow.font}
#' (default \code{"sans"}), e.g. \code{options(crisprFlow.font = "Arial")}.
#'
#' @param counts Count matrix: data frame with columns \code{sgRNA},
#' \code{Gene} and one column of raw counts per sample (see section
#' 'Count matrix').
#' @param meta Sample table: data frame with one row per sample and columns
#' \code{sample_id}, \code{label}, \code{group} and \code{replicate} (see
#' section 'Sample table'). Only samples in \code{meta} are used.
#' @param ctrl_pattern Regular expression matched against the \code{Gene}
#' column to identify control sgRNAs (default is "NEG_CTRL", which matches
#' the controls of all libraries shipped with crisprFlow).
#' @param summary_df Count summary: data frame with columns
#' \code{sample_id}, \code{total_reads} and \code{mapped} (see section
#' 'Count summary').
#' @param colour_by Name of the column in \code{meta} used to colour samples.
#' @param colours Optional named character vector of colours, with names
#' matching the values of the \code{colour_by} column (e.g.
#' \code{c(Untreated = "#0073C2", Drug = "#CD534C")}). If NULL (default),
#' the RColorBrewer "Dark2" palette is used.
#' @param facet_by Optional name of a column in \code{meta} used to split
#' the plot into panels (default is NULL, no panels).
#' @param interactive Logical; if TRUE (default), return an interactive
#' plotly object, otherwise a ggplot object.
#' @param base_size Base font size of the plot.
#'
#' @name qc_inputs
#' @examples
#' ## Simulated screen shipped with crisprFlow
#' head(exampleCounts)
#' exampleSampleMeta
#' qc_check_meta(exampleSampleMeta, exampleCounts)
NULL


## ── Internal helpers ──────────────────────────────────────────────────────────

.qc_font <- function() {
  getOption("crisprFlow.font", "sans")
}

.theme_qc <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = .qc_font()),
      strip.text = ggplot2::element_text(
        size = base_size, face = "bold", family = .qc_font()),
      strip.background = ggplot2::element_rect(
        fill = "grey95", colour = "grey80"),
      legend.position = "bottom",
      panel.spacing = ggplot2::unit(0.8, "lines")
    )
}

## Wrap a ggplot in plotly, with a horizontal legend below the plot
.to_plotly <- function(p, tooltip = "text") {
  plotly::ggplotly(p, tooltip = tooltip) |>
    plotly::layout(
      font = list(family = .qc_font()),
      legend = list(orientation = "h", x = 0, y = -0.15)
    )
}

## Named colours for the values of a metadata column
.qc_colours <- function(values, colours = NULL) {
  values <- if (is.factor(values)) levels(droplevels(values))
            else unique(as.character(values))
  if (!is.null(colours)) {
    missing <- setdiff(values, names(colours))
    if (length(missing) > 0) {
      stop("No colour given for: ", paste(missing, collapse = ", "))
    }
    return(colours)
  }
  pal <- scales::brewer_pal(palette = "Dark2")(8)
  stats::setNames(rep_len(pal, length(values)), values)
}

## Attach metadata to a long table keyed by sample_id, keeping meta row order
.join_meta <- function(df, meta) {
  df |>
    dplyr::inner_join(meta, by = "sample_id") |>
    dplyr::mutate(label = factor(.data$label, levels = meta$label))
}

.check_meta_column <- function(meta, column, arg) {
  if (!is.null(column) && !column %in% colnames(meta)) {
    stop("'", arg, "' column '", column, "' not found in meta")
  }
}

.is_ctrl <- function(counts, ctrl_pattern) {
  stringr::str_detect(counts$Gene, ctrl_pattern)
}


## ── Input checks ──────────────────────────────────────────────────────────────

#' Check a sample table (and its match to a count matrix)
#'
#' Checks that the sample table has the required columns, unique sample IDs
#' and labels, and - if a count matrix is given - that every sample in the
#' sample table has a count column. Called by all QC functions that take
#' \code{meta}, so problems are reported up front with a clear message.
#'
#' @inheritParams qc_inputs
#' @param counts Optional count matrix to check \code{meta} against (see
#' \code{\link{qc_inputs}}).
#'
#' @return \code{meta}, invisibly, if all checks pass (otherwise an error).
#' @export
#'
#' @examples
#' qc_check_meta(exampleSampleMeta, exampleCounts)
#'
qc_check_meta <- function(meta, counts = NULL) {
  assertthat::assert_that(
    is.data.frame(meta),
    msg = "meta must be a data frame"
  )
  required <- c("sample_id", "label", "group", "replicate")
  missing <- setdiff(required, colnames(meta))
  assertthat::assert_that(
    length(missing) == 0,
    msg = paste0("meta is missing required column(s): ",
                 paste(missing, collapse = ", "))
  )
  assertthat::assert_that(
    !anyDuplicated(meta$sample_id),
    !anyDuplicated(meta$label),
    msg = "meta$sample_id and meta$label must be unique"
  )
  if (!is.null(counts)) {
    .check_counts(counts)
    missing <- setdiff(meta$sample_id, colnames(counts))
    assertthat::assert_that(
      length(missing) == 0,
      msg = paste0("Sample(s) in meta without a column in counts: ",
                   paste(missing, collapse = ", "))
    )
  }
  invisible(meta)
}

.check_counts <- function(counts) {
  assertthat::assert_that(
    is.data.frame(counts),
    all(c("sgRNA", "Gene") %in% colnames(counts)),
    msg = "counts must be a data frame with columns 'sgRNA' and 'Gene'"
  )
}


## ── Data loading ──────────────────────────────────────────────────────────────

#' Read a merged count matrix
#'
#' Reads a tab-separated count matrix with columns \code{sgRNA},
#' \code{Gene} and one column per sample (e.g. \code{<prefix>.count.txt}
#' from MAGeCK count, or a merged guide-counter count file).
#'
#' @param counts_file Path to the tab-separated count file.
#'
#' @return A count matrix (see \code{\link{qc_inputs}}).
#' @export
#'
qc_read_counts <- function(counts_file) {
  assertthat::assert_that(
    file.exists(counts_file),
    msg = paste0("counts_file not found: ", counts_file)
  )
  counts <- readr::read_tsv(counts_file, show_col_types = FALSE)
  .check_counts(counts)
  counts
}

#' Read and merge per-sample count files
#'
#' Reads all per-sample count files in a directory (e.g. the
#' \code{<sample_label>.count.txt} files written by
#' \code{\link{count_sgRNA_guidecounter_fastq}}) and merges them into one
#' count matrix. sgRNAs missing from a sample get a count of zero.
#'
#' @param count_dir Directory with the per-sample count files.
#' @param pattern Regular expression matching the count file names
#' (default matches files ending in \code{.count.txt}).
#'
#' @return A count matrix (see \code{\link{qc_inputs}}), with one column
#' per file.
#' @export
#'
qc_read_count_dir <- function(count_dir, pattern = "\\.count\\.txt$") {
  files <- list.files(count_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No count files matching '", pattern, "' in ", count_dir)
  }
  files |>
    purrr::map(\(f) readr::read_tsv(f, show_col_types = FALSE)) |>
    purrr::reduce(\(a, b) dplyr::full_join(a, b, by = c("sgRNA", "Gene"))) |>
    dplyr::mutate(dplyr::across(
      -dplyr::all_of(c("sgRNA", "Gene")), \(x) dplyr::coalesce(x, 0)))
}

## Harmonise column names of MAGeCK / guide-counter count summaries
.standardise_summary <- function(df) {
  names(df) <- tolower(names(df))
  rename_map <- c(
    reads = "total_reads", totalsgrnas = "total_sgrnas",
    zerocounts = "zero_counts", giniindex = "gini_index",
    percentmapped = "percentage"
  )
  hit <- names(df) %in% names(rename_map)
  names(df)[hit] <- rename_map[names(df)[hit]]
  dplyr::rename(df, sample_id = "label")
}

#' Read a count summary file
#'
#' Reads a MAGeCK / guide-counter count summary (e.g.
#' \code{<prefix>.countsummary.txt}) and harmonises the column names, so that
#' it can be used as \code{summary_df} in the QC functions.
#'
#' @param summary_file Path to the tab-separated count summary file. It must
#' contain the columns \code{Label}, \code{Reads} and \code{Mapped}.
#'
#' @return A count summary with one row per sample and the columns
#' \code{sample_id} (from \code{Label}), \code{total_reads}, \code{mapped},
#' and - when present - \code{percentage}, \code{total_sgrnas},
#' \code{zero_counts} and \code{gini_index}.
#' @export
#'
qc_read_count_summary <- function(summary_file) {
  assertthat::assert_that(
    file.exists(summary_file),
    msg = paste0("summary_file not found: ", summary_file)
  )
  readr::read_tsv(summary_file, show_col_types = FALSE) |>
    .standardise_summary()
}

#' Read and combine per-sample count summaries
#'
#' Reads all per-sample count summary files in a directory (e.g. the
#' \code{<sample_label>.countsummary.txt} files written by
#' \code{\link{count_sgRNA_guidecounter_fastq}}) and binds them into one
#' table.
#'
#' @param count_dir Directory with the per-sample summary files.
#' @param pattern Regular expression matching the summary file names
#' (default matches files ending in \code{.countsummary.txt}).
#'
#' @return A count summary, see \code{\link{qc_read_count_summary}}.
#' @export
#'
qc_read_summary_dir <- function(count_dir,
                                pattern = "\\.countsummary\\.txt$") {
  files <- list.files(count_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No summary files matching '", pattern, "' in ", count_dir)
  }
  files |>
    purrr::map(\(f) readr::read_tsv(f, show_col_types = FALSE)) |>
    purrr::list_rbind() |>
    .standardise_summary()
}

#' Read an sgRNA library
#'
#' Reads an sgRNA library, either one shipped with crisprFlow (by name) or
#' a library file (MAGeCK layout: sgRNA ID, sequence, gene; comma- or
#' tab-separated, header optional).
#'
#' @param sgRNA_library Name of a library shipped with crisprFlow
#' ('ACOC','DTKP','GEEX','MEPR','PROT','TMMO','CHIP1'). Ignored if
#' \code{library_file} is given.
#' @param library_file Optional path to a library file.
#'
#' @return A data frame with columns \code{sgRNA}, \code{seq} and
#' \code{Gene}.
#' @export
#'
#' @examples
#' lib <- qc_read_library("CHIP1")
#' head(lib)
#'
qc_read_library <- function(sgRNA_library = NULL, library_file = NULL) {
  if (is.null(library_file)) {
    check_sgRNA_library(sgRNA_library)
    library_file <- system.file(
      "extdata", "sgRNA_library", sgRNA_library,
      paste0("sgRNA_library_", sgRNA_library, ".csv"),
      package = "crisprFlow")
  }
  assertthat::assert_that(
    file.exists(library_file),
    msg = paste0("library_file not found: ", library_file)
  )
  first <- readLines(library_file, n = 1)
  delim <- if (grepl("\t", first)) "\t" else ","
  has_header <- !grepl("^[ACGTN]+$",
                       stringr::str_split_1(first, delim)[2],
                       ignore.case = TRUE)
  readr::read_delim(library_file, delim = delim, col_names = has_header,
                    show_col_types = FALSE) |>
    dplyr::select(1:3) |>
    stats::setNames(c("sgRNA", "seq", "Gene"))
}


## ── Library concordance ───────────────────────────────────────────────────────

#' Check that the counted sgRNAs match the intended library
#'
#' Compares sgRNA IDs in the count matrix with those in the library. A
#' mismatch (e.g. counts generated against a different library) invalidates
#' every downstream step, so this is worth running first.
#'
#' @inheritParams qc_inputs
#' @param library sgRNA library, data frame with an \code{sgRNA} column (see
#' \code{\link{qc_read_library}}).
#' @param min_overlap Minimum fraction of shared sgRNAs - both of the
#' library and of the count matrix - for the two to be called concordant
#' (default is 0.95).
#'
#' @return A one-row tibble with the number of sgRNAs in the library, in
#' the count matrix and in both, the shared fractions, and \code{concordant}
#' (TRUE/FALSE).
#' @export
#'
#' @examples
#' ## Not concordant: the example data use only a subset of CHIP1
#' qc_check_library_concordance(exampleCounts, qc_read_library("CHIP1"))
#'
qc_check_library_concordance <- function(counts, library,
                                         min_overlap = 0.95) {
  .check_counts(counts)
  n_lib <- nrow(library)
  n_counts <- nrow(counts)
  n_shared <- length(intersect(counts$sgRNA, library$sgRNA))
  tibble::tibble(
    n_library_guides = n_lib,
    n_count_guides = n_counts,
    n_shared = n_shared,
    frac_library_in_counts = n_shared / n_lib,
    frac_counts_in_library = n_shared / n_counts,
    concordant = n_shared / n_lib >= min_overlap &&
      n_shared / n_counts >= min_overlap
  )
}


## ── Per-sample statistics ─────────────────────────────────────────────────────

#' Per-sample sgRNA count statistics
#'
#' Computes per-sample metrics from the count matrix alone, complementing
#' the count summary (total reads, mapping rate): number of zero and
#' low-count targeting sgRNAs, median count and the share of reads on
#' control sgRNAs.
#'
#' @inheritParams qc_inputs
#' @param low_threshold Targeting sgRNAs with fewer reads than this are
#' counted as low-count (default is 30).
#'
#' @return A tibble with one row per sample: \code{sample_id},
#' \code{assigned} (total counts), \code{n_zero} and \code{n_low} (number of
#' targeting sgRNAs with zero / fewer than \code{low_threshold} reads),
#' \code{median_count} (median count of targeting sgRNAs) and
#' \code{pct_ctrl} (percent of counts on control sgRNAs).
#' @export
#'
#' @examples
#' qc_sample_stats(exampleCounts, exampleSampleMeta)
#'
qc_sample_stats <- function(counts, meta, ctrl_pattern = "NEG_CTRL",
                            low_threshold = 30) {
  qc_check_meta(meta, counts)
  is_ctrl <- .is_ctrl(counts, ctrl_pattern)
  purrr::map(meta$sample_id, \(s) {
    x <- counts[[s]]
    tibble::tibble(
      sample_id = s,
      assigned = sum(x),
      n_zero = sum(x[!is_ctrl] == 0),
      n_low = sum(x[!is_ctrl] < low_threshold),
      median_count = stats::median(x[!is_ctrl]),
      pct_ctrl = 100 * sum(x[is_ctrl]) / sum(x)
    )
  }) |>
    purrr::list_rbind()
}


## ── Normalisation ─────────────────────────────────────────────────────────────

#' log2 CPM matrix of targeting sgRNAs
#'
#' Normalizes counts of targeting sgRNAs (controls removed) to counts per
#' million and log2-transforms them, log2(CPM + 1).
#'
#' @inheritParams qc_inputs
#' @param sample_ids Sample columns of \code{counts} to include.
#'
#' @return A numeric matrix with one row per targeting sgRNA (row names are
#' sgRNA IDs) and one column per sample.
#' @export
#'
#' @examples
#' m <- qc_log2cpm(exampleCounts, exampleSampleMeta$sample_id)
#' m[1:3, 1:4]
#'
qc_log2cpm <- function(counts, sample_ids, ctrl_pattern = "NEG_CTRL") {
  .check_counts(counts)
  m <- counts[!.is_ctrl(counts, ctrl_pattern), ] |>
    dplyr::select(dplyr::all_of(c("sgRNA", sample_ids))) |>
    tibble::column_to_rownames("sgRNA") |>
    as.matrix()
  log2(sweep(m, 2, colSums(m), "/") * 1e6 + 1)
}


## ── Plot: read mapping ────────────────────────────────────────────────────────

#' Plot mapped and unmapped reads per sample
#'
#' Horizontal stacked bars of mapped and unmapped reads per sample, with the
#' mapping rate printed on each bar. Samples below \code{min_pct_mapped} are
#' marked "(low)".
#'
#' @inheritParams qc_inputs
#' @param min_pct_mapped Minimum acceptable fraction of mapped reads, from 0
#' to 1 (default is 0.65, the minimum used by MAGeCK-VISPR).
#'
#' @return A plotly object (\code{interactive = TRUE}) or a ggplot object.
#' @export
#'
#' @examples
#' qc_plot_read_mapping(exampleCountSummary, exampleSampleMeta,
#'                      interactive = FALSE)
#'
qc_plot_read_mapping <- function(summary_df, meta, colour_by = "group",
                                 colours = NULL, facet_by = NULL,
                                 min_pct_mapped = 0.65,
                                 interactive = TRUE, base_size = 11) {
  qc_check_meta(meta)
  .check_meta_column(meta, colour_by, "colour_by")
  .check_meta_column(meta, facet_by, "facet_by")
  fill_values <- c(.qc_colours(meta[[colour_by]], colours),
                   "Unmapped reads" = "#CCCCCC")

  df <- summary_df |>
    dplyr::select(dplyr::all_of(c("sample_id", "total_reads", "mapped"))) |>
    .join_meta(meta) |>
    dplyr::mutate(
      pct_mapped = .data$mapped / .data$total_reads,
      unmapped = .data$total_reads - .data$mapped,
      fill_key = as.character(.data[[colour_by]]),
      flag = dplyr::if_else(.data$pct_mapped < min_pct_mapped, " (low)", "")
    ) |>
    tidyr::pivot_longer(c("mapped", "unmapped"),
                        names_to = "status", values_to = "reads") |>
    dplyr::mutate(
      fill_key = dplyr::if_else(
        .data$status == "unmapped", "Unmapped reads", .data$fill_key),
      tooltip = paste0(
        "<b>", .data$label, "</b><br>",
        dplyr::if_else(.data$status == "mapped", "Mapped", "Unmapped"),
        " reads: ", scales::comma(.data$reads),
        "<br>% mapped: ", scales::percent(.data$pct_mapped, 0.1)),
      label = factor(.data$label, levels = rev(meta$label)),
      ## mapped reads first (left), unmapped last, in every bar
      fill_key = factor(.data$fill_key, levels = names(fill_values)),
      status = factor(.data$status, levels = c("unmapped", "mapped"))
    )

  p <- ggplot2::ggplot(
    df, ggplot2::aes(y = .data$label, x = .data$reads / 1e6,
                     fill = .data$fill_key, group = .data$status,
                     text = .data$tooltip)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      data = dplyr::filter(df, .data$status == "mapped"),
      ggplot2::aes(label = paste0(
        scales::percent(.data$pct_mapped, 0.1), .data$flag)),
      colour = "grey15", size = 3.2, hjust = 0, x = 0.1,
      family = .qc_font()) +
    ggplot2::scale_fill_manual(values = fill_values, name = NULL) +
    ggplot2::labs(y = NULL, x = "Total reads (millions)") +
    .theme_qc(base_size)

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(.data[[facet_by]]), scales = "free_y")
  }

  if (interactive) .to_plotly(p) else p
}


## ── Plot: generic per-sample metric ───────────────────────────────────────────

#' Plot a per-sample metric
#'
#' Bar plot of any per-sample metric, e.g. Gini index, zero-count sgRNAs or
#' percent reads on controls (see \code{\link{qc_sample_stats}}), with
#' optional reference lines.
#'
#' @inheritParams qc_inputs
#' @param df Data frame with a \code{sample_id} column and the metric
#' column, e.g. output of \code{\link{qc_sample_stats}} or a count summary.
#' @param value Name of the metric column in \code{df}.
#' @param y_label Y-axis title (default is \code{value}).
#' @param thresholds Optional numeric vector of y-values drawn as dashed
#' reference lines (e.g. \code{c(early = 0.1, late = 0.2)} for Gini index).
#' @param value_format Function formatting metric values in tooltips
#' (default shows three decimals).
#'
#' @return A plotly object (\code{interactive = TRUE}) or a ggplot object.
#' @export
#'
#' @examples
#' stats <- qc_sample_stats(exampleCounts, exampleSampleMeta)
#' qc_plot_sample_metric(stats, exampleSampleMeta, value = "pct_ctrl",
#'                       y_label = "% reads on control sgRNAs",
#'                       interactive = FALSE)
#' qc_plot_sample_metric(exampleCountSummary, exampleSampleMeta,
#'                       value = "gini_index", y_label = "Gini index",
#'                       thresholds = c(0.1, 0.2), interactive = FALSE)
#'
qc_plot_sample_metric <- function(
    df, meta, value, y_label = value,
    colour_by = "group", colours = NULL,
    facet_by = NULL, thresholds = NULL,
    value_format = scales::label_number(accuracy = 0.001),
    interactive = TRUE, base_size = 11) {
  qc_check_meta(meta)
  .check_meta_column(meta, colour_by, "colour_by")
  .check_meta_column(meta, facet_by, "facet_by")
  assertthat::assert_that(
    value %in% colnames(df),
    msg = paste0("Column '", value, "' not found in df")
  )
  d <- df |>
    dplyr::select(dplyr::all_of(c("sample_id", value))) |>
    .join_meta(meta) |>
    dplyr::mutate(tooltip = paste0(
      "<b>", .data$label, "</b><br>", y_label, ": ",
      value_format(.data[[value]])))

  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$label, y = .data[[value]],
                    fill = .data[[colour_by]], text = .data$tooltip)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(
      values = .qc_colours(meta[[colour_by]], colours), name = NULL) +
    ggplot2::labs(x = NULL, y = y_label) +
    .theme_qc(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 60, hjust = 1, size = 8))

  if (!is.null(thresholds)) {
    p <- p + ggplot2::geom_hline(
      yintercept = thresholds, linetype = "dashed",
      colour = "grey40", linewidth = 0.5)
  }
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_grid(
      cols = ggplot2::vars(.data[[facet_by]]),
      scales = "free_x", space = "free_x")
  }

  if (interactive) .to_plotly(p) else p
}


## ── Plot: count distributions ─────────────────────────────────────────────────

#' Plot sgRNA count distributions
#'
#' Per-sample density of read counts of targeting sgRNAs (log10 scale),
#' with each sample's median (dashed) and the low-count threshold (dotted).
#' Narrow, overlapping distributions indicate even library representation.
#'
#' @inheritParams qc_inputs
#' @param colour_by Name of the column in \code{meta} used to colour samples
#' (default is "replicate").
#' @param facet_by Optional name of a column in \code{meta} used to split
#' the plot into panels (default is "group").
#' @param low_threshold Read count drawn as a dotted reference line
#' (default is 30).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' qc_plot_count_distribution(exampleCounts, exampleSampleMeta)
#'
qc_plot_count_distribution <- function(counts, meta,
                                       ctrl_pattern = "NEG_CTRL",
                                       colour_by = "replicate",
                                       colours = NULL,
                                       facet_by = "group",
                                       low_threshold = 30,
                                       base_size = 11) {
  qc_check_meta(meta, counts)
  .check_meta_column(meta, colour_by, "colour_by")
  .check_meta_column(meta, facet_by, "facet_by")
  d <- counts[!.is_ctrl(counts, ctrl_pattern), ] |>
    dplyr::select(dplyr::all_of(meta$sample_id)) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "sample_id", values_to = "count") |>
    .join_meta(meta)

  meds <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      unique(c("label", facet_by, colour_by))))) |>
    dplyr::summarise(med = stats::median(.data$count) + 1,
                     .groups = "drop")

  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$count + 1, colour = .data[[colour_by]],
                    group = .data$label)) +
    ggplot2::geom_density(linewidth = 0.6) +
    ggplot2::geom_vline(
      data = meds,
      ggplot2::aes(xintercept = .data$med, colour = .data[[colour_by]]),
      linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = low_threshold, colour = "grey50",
                        linetype = "dotted") +
    ggplot2::scale_x_log10(breaks = 10^(0:5),
                           labels = scales::label_comma()) +
    ggplot2::annotation_logticks(sides = "b") +
    ggplot2::scale_colour_manual(
      values = .qc_colours(meta[[colour_by]], colours), name = colour_by) +
    ggplot2::labs(x = "Read count per sgRNA (log10, +1 pseudocount)",
                  y = "Density") +
    .theme_qc(base_size)

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_by]]))
  }
  p
}


## ── Replicate correlation ─────────────────────────────────────────────────────

#' Plot sample-to-sample correlations
#'
#' Heatmap of pairwise Pearson correlations between samples, computed on
#' log2(CPM + 1) of targeting sgRNAs (see \code{\link{qc_log2cpm}}).
#' Replicates should correlate more strongly with each other than with
#' other groups.
#'
#' @inheritParams qc_inputs
#' @param limits Optional colour scale limits, \code{c(min, max)}. Default
#' (NULL) runs from the lowest correlation (rounded down) to 1.
#' @param show_values Logical; print correlation values in the tiles
#' (default is TRUE).
#' @param base_size Base font size of the plot (default is 10).
#'
#' @return A plotly object (\code{interactive = TRUE}) or a ggplot object.
#' @export
#'
#' @examples
#' qc_plot_replicate_correlation(exampleCounts, exampleSampleMeta,
#'                               interactive = FALSE)
#'
qc_plot_replicate_correlation <- function(counts, meta,
                                          ctrl_pattern = "NEG_CTRL",
                                          limits = NULL, show_values = TRUE,
                                          interactive = TRUE,
                                          base_size = 10) {
  qc_check_meta(meta, counts)
  m <- qc_log2cpm(counts, meta$sample_id, ctrl_pattern)
  cm <- stats::cor(m, method = "pearson")
  lab <- stats::setNames(meta$label, meta$sample_id)

  d <- tibble::as_tibble(cm, rownames = "x") |>
    tidyr::pivot_longer(-"x", names_to = "y", values_to = "r") |>
    dplyr::mutate(
      x = factor(lab[.data$x], levels = meta$label),
      y = factor(lab[.data$y], levels = rev(meta$label)),
      tooltip = paste0(.data$x, " vs ", .data$y,
                       "<br>Pearson r = ", round(.data$r, 4))
    )
  if (is.null(limits)) limits <- c(floor(min(d$r) * 100) / 100, 1)

  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$r,
                    text = .data$tooltip)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_distiller(
      palette = "Blues", direction = 1, limits = limits,
      oob = scales::squish, name = "Pearson r") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = .qc_font()),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right")
  if (show_values) {
    ## white text on the darker tiles
    d$text_colour <- dplyr::if_else(
      d$r > mean(limits), "white", "black")
    p <- p + ggplot2::geom_text(
      data = d,
      ggplot2::aes(label = sprintf("%.2f", .data$r),
                   colour = .data$text_colour),
      size = 2.2, show.legend = FALSE) +
      ggplot2::scale_colour_identity()
  }

  if (interactive) .to_plotly(p) else p
}

#' Within-group replicate correlations
#'
#' Summarises pairwise Pearson correlations (log2 CPM of targeting
#' sgRNAs) between replicates of the same group. Groups with a single
#' sample are left out.
#'
#' @inheritParams qc_inputs
#' @param by Name of the column in \code{meta} defining replicate groups
#' (default is "group").
#'
#' @return A tibble with one row per group: the group, \code{n_reps}
#' (number of replicates), \code{min_r} and \code{mean_r} (lowest and mean
#' pairwise correlation).
#' @export
#'
#' @examples
#' qc_replicate_correlation_table(exampleCounts, exampleSampleMeta)
#'
qc_replicate_correlation_table <- function(counts, meta,
                                           ctrl_pattern = "NEG_CTRL",
                                           by = "group") {
  qc_check_meta(meta, counts)
  .check_meta_column(meta, by, "by")
  cm <- stats::cor(qc_log2cpm(counts, meta$sample_id, ctrl_pattern))
  meta |>
    dplyr::group_by(.data[[by]]) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::summarise(
      n_reps = dplyr::n(),
      r = list(cm[.data$sample_id, .data$sample_id][
        upper.tri(diag(dplyr::n()))]),
      .groups = "drop"
    ) |>
    dplyr::mutate(min_r = purrr::map_dbl(.data$r, min),
                  mean_r = purrr::map_dbl(.data$r, mean)) |>
    dplyr::select(-"r")
}


## ── Plot: PCA ─────────────────────────────────────────────────────────────────

#' Plot a PCA of samples
#'
#' Principal component analysis of log2(CPM + 1) of targeting sgRNAs (see
#' \code{\link{qc_log2cpm}}); sgRNAs without variance are removed.
#' Replicates should cluster together.
#'
#' @inheritParams qc_inputs
#' @param shape_by Optional name of a column in \code{meta} used for point
#' shapes (default is NULL).
#' @param n_top Optional number of most variable sgRNAs to use (default is
#' NULL, all sgRNAs).
#' @param pcs The two principal components to plot (default is
#' \code{c(1, 2)}).
#' @param label_points Logical; label points with \code{meta$label}
#' (default is TRUE; static plots only - interactive plots show labels on
#' hover).
#' @param base_size Base font size of the plot (default is 12).
#'
#' @return A plotly object (\code{interactive = TRUE}) or a ggplot object.
#' @export
#'
#' @examples
#' qc_plot_pca(exampleCounts, exampleSampleMeta, interactive = FALSE)
#'
qc_plot_pca <- function(counts, meta, ctrl_pattern = "NEG_CTRL",
                        colour_by = "group", colours = NULL,
                        shape_by = NULL, n_top = NULL, pcs = c(1, 2),
                        label_points = TRUE, interactive = TRUE,
                        base_size = 12) {
  qc_check_meta(meta, counts)
  .check_meta_column(meta, colour_by, "colour_by")
  .check_meta_column(meta, shape_by, "shape_by")
  m <- qc_log2cpm(counts, meta$sample_id, ctrl_pattern)
  v <- apply(m, 1, stats::var)
  m <- m[v > 0, , drop = FALSE]
  if (!is.null(n_top)) {
    v <- v[v > 0]
    m <- m[order(-v)[seq_len(min(n_top, nrow(m)))], , drop = FALSE]
  }

  pca <- stats::prcomp(t(m), center = TRUE, scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pc <- paste0("PC", pcs)

  d <- tibble::as_tibble(pca$x[, pc], rownames = "sample_id") |>
    stats::setNames(c("sample_id", "PCa", "PCb")) |>
    .join_meta(meta) |>
    dplyr::mutate(tooltip = paste0(
      "<b>", .data$label, "</b><br>",
      pc[1], ": ", round(.data$PCa, 2), "<br>",
      pc[2], ": ", round(.data$PCb, 2)))

  aes_shape <- if (!is.null(shape_by)) {
    ggplot2::aes(shape = .data[[shape_by]])
  } else {
    ggplot2::aes()
  }
  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$PCa, y = .data$PCb,
                    colour = .data[[colour_by]], text = .data$tooltip)) +
    ggplot2::geom_point(aes_shape, size = 3.5, alpha = 0.9) +
    ggplot2::scale_colour_manual(
      values = .qc_colours(meta[[colour_by]], colours), name = NULL) +
    ggplot2::labs(x = paste0(pc[1], ": ", pct[pcs[1]], "% variance"),
                  y = paste0(pc[2], ": ", pct[pcs[2]], "% variance"),
                  shape = NULL) +
    .theme_qc(base_size)
  if (label_points && !interactive) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$label), size = 2.8,
      show.legend = FALSE, family = .qc_font(),
      max.overlaps = 30, seed = 42)
  }

  if (interactive) .to_plotly(p) else p
}


## ── Guide representation ──────────────────────────────────────────────────────

#' Detect sgRNAs at one or more reference stages
#'
#' Classifies each targeting sgRNA as detected or not at each stage of the
#' screen (e.g. plasmid, T0), and summarises sgRNA loss per gene. Use it to
#' see how much of the library is lost before selection starts.
#'
#' @inheritParams qc_inputs
#' @param stages Named list; each element gives the sample IDs (columns of
#' \code{counts}) pooled for that stage, e.g.
#' \code{list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3"))}. The
#' order of the list sets the order of stages.
#' @param min_reads An sgRNA is detected at a stage if its mean count per
#' sample in that stage is at least \code{min_reads} (default is 30).
#'
#' @return A list with
#' \describe{
#'   \item{guide_summary}{one row per sgRNA and stage: \code{sgRNA},
#'   \code{Gene}, \code{stage}, \code{mean_count}, \code{detected}}
#'   \item{gene_summary}{one row per gene and stage: \code{n_guides},
#'   \code{n_detected}, \code{n_lost}, \code{frac_lost}}
#'   \item{threshold}{the \code{min_reads} used}
#' }
#' Pass it to \code{\link{qc_plot_guide_detection}} and
#' \code{\link{qc_genes_with_lost_guides}}.
#' @export
#'
#' @examples
#' rep_data <- qc_guide_representation(
#'   exampleCounts,
#'   stages = list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3")))
#' head(rep_data$gene_summary)
#'
qc_guide_representation <- function(counts, stages,
                                    ctrl_pattern = "NEG_CTRL",
                                    min_reads = 30) {
  .check_counts(counts)
  assertthat::assert_that(
    is.list(stages), !is.null(names(stages)), all(names(stages) != ""),
    msg = "stages must be a named list of sample IDs"
  )
  missing <- setdiff(unlist(stages), colnames(counts))
  assertthat::assert_that(
    length(missing) == 0,
    msg = paste0("Sample(s) in stages not found in counts: ",
                 paste(missing, collapse = ", "))
  )
  targeting <- counts[!.is_ctrl(counts, ctrl_pattern), ]

  per_stage <- purrr::imap(stages, \(ids, stage) {
    mean_count <- rowMeans(as.matrix(targeting[, ids, drop = FALSE]))
    tibble::tibble(sgRNA = targeting$sgRNA, Gene = targeting$Gene,
                   stage = stage, mean_count = mean_count,
                   detected = mean_count >= min_reads)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(stage = factor(.data$stage, levels = names(stages)))

  gene_summary <- per_stage |>
    dplyr::group_by(.data$stage, .data$Gene) |>
    dplyr::summarise(n_guides = dplyr::n(),
                     n_detected = sum(.data$detected), .groups = "drop") |>
    dplyr::mutate(n_lost = .data$n_guides - .data$n_detected,
                  frac_lost = .data$n_lost / .data$n_guides)

  list(guide_summary = per_stage, gene_summary = gene_summary,
       threshold = min_reads)
}

#' Plot detected and lost sgRNAs per stage
#'
#' Stacked bars of detected vs. absent/low targeting sgRNAs at each stage.
#'
#' @param rep_data Output of \code{\link{qc_guide_representation}}.
#' @param base_size Base font size of the plot (default is 12).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' rep_data <- qc_guide_representation(
#'   exampleCounts,
#'   stages = list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3")))
#' qc_plot_guide_detection(rep_data)
#'
qc_plot_guide_detection <- function(rep_data, base_size = 12) {
  d <- rep_data$guide_summary |>
    dplyr::count(.data$stage, .data$detected) |>
    dplyr::group_by(.data$stage) |>
    dplyr::mutate(pct = .data$n / sum(.data$n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(status = factor(
      dplyr::if_else(.data$detected, "Detected", "Absent/low"),
      levels = c("Absent/low", "Detected")))

  ggplot2::ggplot(
    d, ggplot2::aes(x = .data$stage, y = .data$n, fill = .data$status)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(scales::comma(.data$n), "\n(",
                                  scales::percent(.data$pct, 0.1), ")")),
      position = ggplot2::position_stack(vjust = 0.5), size = 3.2,
      colour = "white", fontface = "bold", family = .qc_font()) +
    ggplot2::scale_fill_manual(
      values = c(Detected = "#4393C3", `Absent/low` = "#D6604D"),
      name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = NULL, y = "Number of targeting sgRNAs",
                  subtitle = paste0("Detected = mean >= ",
                                    rep_data$threshold,
                                    " reads per sample in stage")) +
    .theme_qc(base_size) +
    ggplot2::theme(legend.position = "right")
}

#' Genes that lost many of their sgRNAs
#'
#' Lists genes where at least a given fraction of sgRNAs were not detected
#' at a stage. Hits for such genes rest on few sgRNAs and should be
#' interpreted with care.
#'
#' @param rep_data Output of \code{\link{qc_guide_representation}}.
#' @param stage Name of the stage (one of the names of \code{stages} given
#' to \code{\link{qc_guide_representation}}).
#' @param min_frac Minimum fraction of lost sgRNAs (default is 0.4, e.g. 2
#' of 5).
#'
#' @return A tibble from \code{rep_data$gene_summary}, one row per gene,
#' sorted by fraction and number of lost sgRNAs.
#' @export
#'
#' @examples
#' rep_data <- qc_guide_representation(
#'   exampleCounts,
#'   stages = list(Plasmid = "Plasmid", T0 = c("T0_1", "T0_2", "T0_3")))
#' qc_genes_with_lost_guides(rep_data, stage = "T0")
#'
qc_genes_with_lost_guides <- function(rep_data, stage, min_frac = 0.4) {
  assertthat::assert_that(
    stage %in% levels(rep_data$gene_summary$stage),
    msg = paste0("stage must be one of: ",
                 paste(levels(rep_data$gene_summary$stage), collapse = ", "))
  )
  rep_data$gene_summary |>
    dplyr::filter(.data$stage == !!stage, .data$frac_lost >= min_frac) |>
    dplyr::arrange(dplyr::desc(.data$frac_lost), dplyr::desc(.data$n_lost))
}


## ── Plot: control guides ──────────────────────────────────────────────────────

#' Plot control sgRNA abundance across samples
#'
#' Heatmap of the abundance (log10 CPM + 1) of each control sgRNA in each
#' sample. Controls should be phenotypically neutral, so their abundance
#' should be stable; missing or outlier controls stand out as off-colour
#' tiles.
#'
#' @inheritParams qc_inputs
#' @param base_size Base font size of the plot (default is 10).
#'
#' @return A plotly object (\code{interactive = TRUE}) or a ggplot object;
#' NULL (with a message) if no sgRNA matches \code{ctrl_pattern}.
#' @export
#'
#' @examples
#' qc_plot_control_guides(exampleCounts, exampleSampleMeta,
#'                        interactive = FALSE)
#'
qc_plot_control_guides <- function(counts, meta, ctrl_pattern = "NEG_CTRL",
                                   interactive = TRUE, base_size = 10) {
  qc_check_meta(meta, counts)
  is_ctrl <- .is_ctrl(counts, ctrl_pattern)
  if (!any(is_ctrl)) {
    message("No control guides match '", ctrl_pattern, "'")
    return(NULL)
  }
  tot <- colSums(counts[meta$sample_id])
  d <- counts[is_ctrl, ] |>
    dplyr::select(dplyr::all_of(c("sgRNA", "Gene", meta$sample_id))) |>
    tidyr::pivot_longer(-c("sgRNA", "Gene"),
                        names_to = "sample_id", values_to = "count") |>
    dplyr::mutate(cpm = .data$count / tot[.data$sample_id] * 1e6) |>
    .join_meta(meta) |>
    dplyr::mutate(tooltip = paste0(
      "<b>", .data$sgRNA, "</b><br>", .data$label,
      "<br>count: ", scales::comma(.data$count),
      "<br>CPM: ", round(.data$cpm, 1)))

  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$label, y = .data$sgRNA,
                    fill = log10(.data$cpm + 1), text = .data$tooltip)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "log10(CPM+1)") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = .qc_font()),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 6),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right")

  if (interactive) .to_plotly(p) else p
}
