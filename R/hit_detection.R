#' Run MAGeCK MLE for hit detection
#'
#' This function runs MAGeCK MLE for hit detection using a provided sample sheet
#' and count table. It constructs a design matrix based on the sample information
#' and executes the MAGeCK MLE command with specified parameters.
#' The results are saved in the specified output directory.
#'
#' @param samplesheet_csv Path to the sample sheet CSV file
#' @param count_table_txt Path to the count table text file (from MAGeCK
#' count or guidecounter)
#' @param output_dir Directory where output files will be saved
#' @param norm_method Normalization method to use ("median", "total",
#' "control" or "none"; default is "median"). With "control", counts are
#' normalized to the control genes in the count table (genes matching
#' \code{control_pattern})
#' @param mageck_bin Path to the MAGeCK binary
#' @param threads Number of threads to use (default is 1)
#' @param permutation_round Number of permutation rounds (default is 2, recommended 10 (longer time))
#' @param late_time_start Time point to start considering as "late" (default is
#' 10)
#' @param adjust_method Method for p-value adjustment (default is "fdr")
#' @param output_prefix Prefix for output files (default is "crisprFlowRun")
#' @param control_pattern Pattern to identify control genes (default is "NEG_CTRL")
#' @return A data frame containing the gene summary results from MAGeCK MLE
#'
#' @export
#'
detect_hits_mageck_mle <- function(
    samplesheet_csv = NULL,
    count_table_txt = NULL,
    output_dir = NULL,
    norm_method = c("median", "total", "control", "none"),
    mageck_bin = "/Users/sigven/miniconda3/bin/mageck",
    threads = 1,
    permutation_round = 4,
    late_time_start = 10,
    adjust_method = "fdr",
    output_prefix = "mageck_mle_hits",
    control_pattern = "NEG_CTRL"){

  norm_method <- match.arg(norm_method)
  assertthat::assert_that(
    !is.null(mageck_bin),
    msg = "Please provide path to mageck binary"
  )
  assertthat::assert_that(
    file.exists(mageck_bin),
    msg = glue::glue("mageck binary not found at: ", mageck_bin)
  )
  if(is.null(count_table_txt) | !file.exists(count_table_txt)){
    stop("Please provide path to '<project>.counts.txt' file (from mageck count or guidecounter)")
  }
  if(is.null(samplesheet_csv) | !file.exists(samplesheet_csv)){
    stop("Please provide path to sample_sheet CSV file")
  }
  sample_data <-
    readr::read_csv(file = samplesheet_csv,
                    show_col_types = F, na = c(".","", "NA"))
  if(!all(c('fastq','sample_label','condition2',
            'condition','T0','time_point') %in%
          colnames(sample_data))){
    stop(
      paste0("sample_sheet must contain columns: ",
             "fastq, sample_label, condition, T0, time_point"))
  }

  assertthat::assert_that(
    dir.exists(output_dir),
    msg = glue::glue("output_dir not found at: '{output_dir}'")
  )

  ## Control genes (for --control-gene) are taken from the count table,
  ## so that they match the sgRNA library used for counting
  control_genes_fname <- NULL
  if(norm_method == "control"){
    control_genes <- readr::read_tsv(
      count_table_txt, show_col_types = F) |>
      dplyr::filter(
        stringr::str_detect(.data$Gene, pattern = control_pattern)) |>
      dplyr::pull(.data$Gene) |>
      unique() |>
      sort()
    assertthat::assert_that(
      length(control_genes) > 0,
      msg = glue::glue(
        "No control genes matching '{control_pattern}' found in ",
        "count table - cannot use norm_method = 'control'")
    )
    control_genes_fname <- file.path(
      output_dir, paste0(output_prefix, "_control_genes.txt"))
    writeLines(control_genes, control_genes_fname)
  }
  design_matrix_fname <- file.path(
    output_dir, paste0(output_prefix, "_design_matrix.txt"))


  ## https://sourceforge.net/p/mageck/wiki/advanced_tutorial/
  ## h-tutorial-4-make-full-use-of-mageck-mle-for-more-complicated-experimental-design-eg-paired-samples-time-series
  design_matrix <-
    sample_data |>
    dplyr::mutate(
      sample = sample_label,
      baseline = 1, ## apparently 1 for all samples
      common = ifelse(T0 != 1, 1, 0),
      late_vs_early = ifelse(time_point >= late_time_start, 1, 0),
      drug_treated = ifelse(condition2 == "drug", 1, 0)
    ) |>
    dplyr::select(
      sample,
      baseline,
      common,
      late_vs_early,
      drug_treated
    ) |>
    dplyr::distinct() |>
    dplyr::arrange(
      sample
    )

  readr::write_tsv(
    design_matrix,
    file = design_matrix_fname,
    col_names = TRUE, quote = "none")


  args <- c(
    "mle",
    "--count-table",glue::glue("{count_table_txt}"),
    "--design-matrix",glue::glue("{design_matrix_fname}"),
    "--adjust-method",glue::glue("{adjust_method}"),
    "--norm-method",glue::glue("{norm_method}"),
    "--permutation-round",glue::glue("{permutation_round}"),
    "-n",glue::glue("{file.path(output_dir, output_prefix)}")
  )

  if(!is.null(control_genes_fname)){
    args <- c(args, "--control-gene",glue::glue("{control_genes_fname}"))
  }

  message("Running 'mageck mle' with command:\n")
  message(glue::glue("{mageck_bin} {paste(args, collapse=' ')}\n"))

  p <- processx::process$new(
    mageck_bin, args,
    stdout = file.path(
      output_dir, glue::glue("{output_prefix}_stdout.log")),
    stderr = file.path(
      output_dir, glue::glue("{output_prefix}_stderr.log"))
  )
  p$wait()
  if(p$get_exit_status() != 0){
    stop("Error running mageck mle, please check log files")
  }
  message("'mageck mle' finished successfully.\n")

  gene_output_fname <-
    file.path(
      output_dir,
      glue::glue("{output_prefix}.gene_summary.txt")
    )

  hits <- data.frame()

  if(!file.exists(gene_output_fname)){
    stop(glue::glue(
      "Expected output file not found: {gene_output_fname} Please check log files for errors."))
  } else {
    hits <- as.data.frame(readr::read_tsv(
      gene_output_fname,
      show_col_types = F,
      na = c(".","", "NA"))) |>
      janitor::clean_names() |>
      dplyr::rename(
        sgRNA = sg_rna,
        Gene = gene,
      ) |>
      dplyr::filter(
        !stringr::str_detect(
          .data$Gene, pattern = control_pattern)
      ) |>
      dplyr::arrange(
        .data$common_fdr
      )
    message(glue::glue("Output file written: {gene_output_fname}\n"))

  }

  return(hits)

}

#' Run FPA for hit detection
#'
#' This function runs the FPA (Fitness Perturbation Analysis) method for hit detection
#' using a provided count table. It processes the input data, executes the FPA analysis
#' with specified parameters, and returns the results.
#'
#' @param data Data frame containing sgRNA count data (from MAGeCK count or guidecounter)
#' @param fpa_t0_cutoff Minimum count threshold for T0 (default is 0)
#' @param fpa_differential Boolean indicating if differential analysis should be performed
#' (default is FALSE)
#' @param fpa_control_condition Condition to be used as control (default is "DMSO")
#' @param fpa_smoothing Smoothing method to use ("Gradient", "Laplace", or "Both")
#' @param fpa_fitness_computation Fitness computation method to use ("Regression", "SVD", or "Both")
#' @param fpa_statistics_method Statistics method to use ("Density" or "Rank")
#' @param fpa_statistic Statistic to use ("Median", "Sum", or "Both")
#' @param fpa_iterations Number of iterations for the analysis (default is 10000)
#' @param output_dir Directory where output files will be saved
#' @param output_prefix Prefix for output files (default is "crisprFlowRun")
#' @param control_pattern Pattern to identify control genes (default is "NEG_CTRL")
#' @return A list containing the FPA analysis results
#'
#' @export
#'
detect_hits_fpa <- function(
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
    control_pattern = "NEG_CTRL"){

  if(file.exists(
    file.path(output_dir, glue::glue("{output_prefix}.rds")))){
    message(glue::glue("FPA output file already exists: ",
                       "{file.path(output_dir, glue::glue('{output_prefix}.rds'))}"))
    hits <- readRDS(
      file.path(output_dir, glue::glue("{output_prefix}.rds")))
    return(hits)

  }

  assertthat::assert_that(
    !is.null(data),
    msg = "Please provide data frame with sgRNA counts (from MAGeCK count or guidecounter)"
  )
  assertthat::assert_that(
    is.data.frame(data),
    msg = "data must be a data frame"
  )
  required_cols <- c('sgRNA','Gene')
  assertable::assert_colnames(
    data,
    required_cols,
    only_colnames = F,
    quiet = T
  )

  if(is.null(output_dir) | !dir.exists(output_dir)){
    stop("Please provide path to output_dir")
  }

  ##smoothing can be either "Gradient, "Laplace", or "Both"
  assertthat::assert_that(
    fpa_smoothing %in% c("Gradient","Laplace","Both"),
    msg = "fpa_smoothing must be one of: 'Gradient', 'Laplace', or 'Both'"
  )

  ## fitness computation can be either "Regression", "SVD", or "Both"
  assertthat::assert_that(
    fpa_fitness_computation %in% c("Regression","SVD","Both"),
    msg = "fpa_fitness_computation must be one of: 'Regression', 'SVD', or 'Both'"
  )
  ## statistics method can be either "Density", or "Rank"
  assertthat::assert_that(
    fpa_statistics_method %in% c("Density","Rank"),
    msg = "fpa_statistics_method must be one of: 'Density', or 'Rank'"
  )
  ## statistic can be either "Median", "Sum", "Both"
  assertthat::assert_that(
    fpa_statistic %in% c("Median","Sum","Both"),
    msg = "fpa_statistic must be one of: 'Median', 'Sum', or 'Both'"
  )

  hits <-
    CRISPRFitnessPerturbationAnalysis::CRISPR.screen.analysis(
      data = data,
      Differential = fpa_differential,
      T0.cutoff = fpa_t0_cutoff,
      Control.condition = fpa_control_condition,
      Smoothing = fpa_smoothing,
      Fitness.computation = fpa_fitness_computation,
      sgRNA.control.id = "NEG_CTRL",
      Species = "Human",
      Statistics.method = fpa_statistics_method,
      Statistic = fpa_statistic,
      Iterations = fpa_iterations
    )

  saveRDS(
    hits,
    file = file.path(
      output_dir, glue::glue("{output_prefix}.rds"))
  )
  message(glue::glue("FPA output file written: ",
                     "{file.path(output_dir, glue::glue('{output_prefix}.rds'))}\n"))

  return(hits)

}

#' Run DrugZ for hit detection
#'
#' This function is a placeholder for running DrugZ for hit detection.
#'
#' @param samplesheet_csv Path to the sample sheet CSV file
#' @param count_table_txt Path to the count table text file (from MAGeCK/guidecounter)
#' @param pseudo_count Pseudo count to add (default is 5)
#' @param output_dir Directory where output files will be saved
#' @param output_prefix Prefix for output files (default is "FGF401")
#' @param drugz_python Path to the DrugZ python interpreter (Conda environment)
#'
#' @return A list containing DrugZ hit detection results
#'
#' @export
#'
detect_hits_drugz <- function(
    samplesheet_csv = NULL,
    count_table_txt = NULL,
    pseudo_count = 5,
    output_dir = NULL,
    output_prefix = "FGF401",
    drugz_python = "/Users/sigven/miniconda3/envs/drugz/bin/python") {

  assertthat::assert_that(
    !is.null(drugz_python),
    msg = "Please provide path to drugZ script"
  )
  assertthat::assert_that(
    file.exists(drugz_python),
    msg = glue::glue("drugZ script not found at: ", drugz_python)
  )
  if(is.null(count_table_txt) | !file.exists(count_table_txt)){
    stop("Please provide path to '<project>.counts.txt' file (from mageck count or guidecounter)")
  }
  if(is.null(samplesheet_csv) | !file.exists(samplesheet_csv)){
    stop("Please provide path to sample_sheet CSV file")
  }

  assertthat::assert_that(
    is.numeric(pseudo_count) & pseudo_count > 0,
    msg = "pseudo_count must be numeric and greater than 0"
  )

  assertthat::assert_that(
    !is.null(output_dir) & dir.exists(output_dir),
    msg = "Please provide path to output_dir"
  )

  assertthat::assert_that(
    is.character(output_prefix),
    msg = "output_prefix must be a character string"
  )

  sample_data <-
    readr::read_csv(file = samplesheet_csv,
                    show_col_types = F, na = c(".","", "NA"))
  if(!all(c('fastq','sample_label','condition2',
            'condition','T0','time_point') %in%
          colnames(sample_data))){
    stop(
      paste0("sample_sheet must contain columns: ",
             "fastq, sample_label, condition, T0, time_point"))
  }

  control_samples <- sample_data |>
    dplyr::filter(condition2 == "control" & time_point != 0) |>
    dplyr::pull(sample_label) |>
    paste0(collapse = ",")

  treatment_samples <- sample_data |>
    dplyr::filter(condition2 == "drug" & time_point != 0) |>
    dplyr::pull(sample_label) |>
    paste0(collapse = ",")

  assertthat::assert_that(
    dir.exists(output_dir),
    msg = glue::glue("output_dir not found at: '{output_dir}'")
  )

  drugz_script <-
    file.path(
      system.file(
        'extdata', package='crisprFlow'), 'drugz.py')

  assertthat::assert_that(
    file.exists(drugz_script),
    msg = glue::glue("drugZ script not found at: ", drugz_script)
  )

  hits_output_fname <-
    file.path(
      output_dir, glue::glue("{output_prefix}.txt")
    )
  hits_output_fname_fc <-
    file.path(
      output_dir, glue::glue("{output_prefix}_fc.txt")
    )

  args <- c(
    drugz_script,
    "-i",glue::glue("{count_table_txt}"),
    "-c",glue::glue("{control_samples}"),
    "-x",glue::glue("{treatment_samples}"),
    "-p",glue::glue("{pseudo_count}"),
    "-f",glue::glue("{hits_output_fname_fc}"),
    "-o",glue::glue("{hits_output_fname}")
  )

  message("Running 'drugZ' with command:\n")
  message(glue::glue("{drugz_python} {drugz_script} {paste(args, collapse=' ')}\n"))

  p <- processx::process$new(
    drugz_python, args,
    stdout = file.path(
      output_dir, glue::glue("{output_prefix}_stdout.log")),
    stderr = file.path(
      output_dir, glue::glue("{output_prefix}_stderr.log"))
  )
  p$wait()
  if(p$get_exit_status() != 0){
    stop("Error running drugZ, please check log files")
  }
  message("'drugZ' finished successfully.\n")


  hits <- list()

  hits[['df']] <- data.frame()
  hits[['fc']] <- data.frame()

  if(!file.exists(hits_output_fname)){
    stop(glue::glue(
      "Expected output file not found: {hits_output_fname}",
      "Please check log files for errors."))
  } else {
    hits[['df']] <- as.data.frame(readr::read_tsv(
      hits_output_fname,
      show_col_types = F,
      na = c(".","", "NA"))) |>
      janitor::clean_names()
    message(
      glue::glue("Output file written: {hits_output_fname}\n"))
  }

  if(!file.exists(hits_output_fname_fc)){
    stop(glue::glue(
      "Expected output file not found: ",
      "{hits_output_fname_fc} ",
      "Please check log files for errors."))
  } else {
    hits[['fc']] <- as.data.frame(readr::read_tsv(
      hits_output_fname_fc,
      show_col_types = F,
      na = c(".","", "NA"))) |>
      janitor::clean_names()
    message(glue::glue(
      "Output file written: ",
      "{hits_output_fname_fc}\n"))
  }

  return(hits)

}



#' Estimate log2 fold change of sgRNA counts between arbitrary conditions
#'
#' Normalizes raw sgRNA counts to CPM and computes log2 fold changes for any
#' set of condition pairs specified via a comparisons table. Replicate columns
#' are resolved through the samplesheet and collapsed either by averaging CPM
#' before computing FC, or by computing per-replicate FC and then averaging.
#'
#' The samplesheet must contain at least \code{sample_label} and
#' \code{condition} columns. If a \code{replicate} column is present, count
#' column names are expected to follow the pattern
#' \code{<sample_label>_<replicate>} (e.g. \code{Bottom_B1}); otherwise just
#' \code{<sample_label>} is used.
#'
#' When only a single replicate exists for either condition in a comparison,
#' \code{replicate_method} is ignored and the single-column CPM value is used
#' directly.
#'
#' @param counts A data frame with columns \code{sgRNA}, \code{Gene}, and one
#'   column per sample whose names match those derived from the samplesheet.
#' @param samplesheet A data frame with at least columns \code{sample_label}
#'   and \code{condition}. An optional \code{replicate} column drives the
#'   mapping from condition to count columns.
#' @param comparisons A data frame with columns \code{numerator},
#'   \code{denominator}, and \code{label} defining which condition pairs to
#'   compare and the prefix used for the output FC columns.
#' @param replicate_method How to handle multiple replicates per condition.
#'   \code{"mean_cpm_then_fc"} (default) averages CPM across replicates first,
#'   then computes a single log2FC. \code{"per_replicate_fc"} computes log2FC
#'   for each replicate pair (matched when counts are equal, all-pairwise
#'   otherwise) and reports the mean; a \code{_log2FC_sd} column is also added
#'   when more than one pair exists.
#' @param pseudo_count Pseudo-count added before taking log2 ratios to avoid
#'   log(0) (default 1).
#' @param output_dir Directory where output TSV files are written.
#' @param output_prefix Prefix for the output file names (default
#'   \code{"crisprFlow"}).
#' @param overwrite If \code{FALSE} (default) and output files already exist,
#'   they are read and returned without recomputing.
#' @param control_pattern Pattern identifying negative-control sgRNAs; not used
#'   in the computation but retained for consistency with the rest of the
#'   package API.
#'
#' @return A named list with elements \code{per_sgRNA} and \code{per_gene}.
#'   Each is a data frame containing mean CPM columns for every referenced
#'   condition, log2FC columns (and optionally SD columns) for every
#'   comparison, and — for \code{per_gene} — values averaged across the sgRNAs
#'   targeting each gene.
#'
#' @export
#'
estimate_fold_change <- function(
    counts = NULL,
    samplesheet = NULL,
    comparisons = NULL,
    replicate_method = c("mean_cpm_then_fc",
                         "per_replicate_fc"),
    pseudo_count = 1,
    output_dir = NULL,
    output_prefix = "crisprFlow",
    overwrite = FALSE,
    control_pattern = "NEG_CTRL") {

  replicate_method <- match.arg(replicate_method)

  fold_change_results <- list(
    per_sgRNA = data.frame(),
    per_gene  = data.frame()
  )

  assertthat::assert_that(
    !is.null(output_dir) & dir.exists(output_dir),
    msg = "Please provide a valid path to output_dir"
  )

  output_fname_fc_sgrna <- file.path(
    output_dir, glue::glue("{output_prefix}_fc_sgRNA.tsv"))
  output_fname_fc_gene <- file.path(
    output_dir, glue::glue("{output_prefix}_fc_gene.tsv"))

  if (!overwrite &
      file.exists(output_fname_fc_sgrna) &
      file.exists(output_fname_fc_gene)) {
    message(glue::glue(
      "Fold change output files already exist:\n",
      "{output_fname_fc_sgrna}\n",
      "{output_fname_fc_gene}\n"))
    fold_change_results[["per_sgRNA"]] <- as.data.frame(
      readr::read_tsv(output_fname_fc_sgrna,
                      show_col_types = FALSE, na = c(".", "", "NA")))
    fold_change_results[["per_gene"]] <- as.data.frame(
      readr::read_tsv(output_fname_fc_gene,
                      show_col_types = FALSE, na = c(".", "", "NA")))
    return(fold_change_results)
  }

  ## --- input validation ---------------------------------------------------

  assertthat::assert_that(
    !is.null(counts) & is.data.frame(counts),
    msg = "counts must be a data frame"
  )
  assertable::assert_colnames(
    counts, c("sgRNA", "Gene"), only_colnames = FALSE, quiet = TRUE)

  assertthat::assert_that(
    !is.null(samplesheet) & is.data.frame(samplesheet),
    msg = "samplesheet must be a data frame"
  )
  assertable::assert_colnames(
    samplesheet, c("sample_label", "condition"),
    only_colnames = FALSE, quiet = TRUE)

  assertthat::assert_that(
    !is.null(comparisons) & is.data.frame(comparisons),
    msg = "comparisons must be a data frame"
  )
  assertable::assert_colnames(
    comparisons, c("numerator", "denominator", "label"),
    only_colnames = FALSE, quiet = TRUE)

  assertthat::assert_that(
    is.numeric(pseudo_count) & pseudo_count >= 0,
    msg = "pseudo_count must be a non-negative number"
  )

  ## --- build condition -> count-column mapping ----------------------------

  has_replicate_col <-
    "replicate" %in% colnames(samplesheet) &&
    !all(is.na(samplesheet[["replicate"]]))

  sample_map <- samplesheet |>
    dplyr::select(dplyr::any_of(
      c("sample_label", "condition", "replicate"))) |>
    dplyr::distinct() |>
    dplyr::mutate(
      count_col = if (has_replicate_col) {
        paste0(.data$sample_label, "_", .data$replicate)
      } else {
        .data$sample_label
      }
    ) |>
    dplyr::filter(.data$count_col %in% colnames(counts))

  if (nrow(sample_map) == 0) {
    stop(paste0(
      "No count columns could be matched from the samplesheet. ",
      "Check that sample_label",
      if (has_replicate_col) " + replicate" else "",
      " columns correspond to count table column names."))
  }

  ## --- CPM normalisation --------------------------------------------------

  ## Normalise every matched sample column independently so that columns
  ## with different sequencing depths are comparable.
  sample_cols <- unique(sample_map$count_col)
  norm_counts <- counts
  norm_counts[, sample_cols] <- lapply(
    counts[, sample_cols, drop = FALSE],
    function(x) x / sum(x, na.rm = TRUE) * 1e6
  )

  ## --- build result data frame --------------------------------------------

  result_df <- counts[, c("sgRNA", "Gene")]

  ## Add mean CPM per condition for every condition referenced in comparisons
  all_conditions <- unique(
    c(comparisons$numerator, comparisons$denominator))

  for (cond in all_conditions) {
    cond_cols <- sample_map |>
      dplyr::filter(.data$condition == cond) |>
      dplyr::pull(.data$count_col)

    if (length(cond_cols) == 0) {
      stop(glue::glue("No count columns found for condition: '{cond}'"))
    }

    cpm_col <- paste0(cond, "_mean_CPM")
    result_df[[cpm_col]] <- round(
      if (length(cond_cols) == 1) {
        norm_counts[[cond_cols]]
      } else {
        rowMeans(norm_counts[, cond_cols, drop = FALSE], na.rm = TRUE)
      },
      digits = 4)
  }

  ## --- compute fold changes -----------------------------------------------

  for (i in seq_len(nrow(comparisons))) {
    num_cond   <- comparisons$numerator[i]
    den_cond   <- comparisons$denominator[i]
    comp_label <- comparisons$label[i]

    num_cols <- sample_map |>
      dplyr::filter(.data$condition == num_cond) |>
      dplyr::pull(.data$count_col)
    den_cols <- sample_map |>
      dplyr::filter(.data$condition == den_cond) |>
      dplyr::pull(.data$count_col)

    use_mean_first <-
      replicate_method == "mean_cpm_then_fc" ||
      length(num_cols) <= 1 ||
      length(den_cols) <= 1

    fc_col <- paste0(comp_label, "_log2FC")

    if (use_mean_first) {
      avg_num <- if (length(num_cols) == 1) {
        norm_counts[[num_cols]]
      } else {
        rowMeans(norm_counts[, num_cols, drop = FALSE], na.rm = TRUE)
      }
      avg_den <- if (length(den_cols) == 1) {
        norm_counts[[den_cols]]
      } else {
        rowMeans(norm_counts[, den_cols, drop = FALSE], na.rm = TRUE)
      }
      result_df[[fc_col]] <- round(
        log2((avg_num + pseudo_count) / (avg_den + pseudo_count)),
        digits = 4)

    } else {
      ## per-replicate FC: matched when equal counts, all-pairwise otherwise
      if (length(num_cols) == length(den_cols)) {
        pairs <- data.frame(
          num = num_cols, den = den_cols, stringsAsFactors = FALSE)
      } else {
        pairs <- expand.grid(
          num = num_cols, den = den_cols, stringsAsFactors = FALSE)
      }

      fc_matrix <- vapply(seq_len(nrow(pairs)), function(p) {
        log2((norm_counts[[pairs$num[p]]] + pseudo_count) /
               (norm_counts[[pairs$den[p]]] + pseudo_count))
      }, numeric(nrow(norm_counts)))

      if (!is.matrix(fc_matrix)) {
        fc_matrix <- matrix(fc_matrix, ncol = 1)
      }

      result_df[[fc_col]] <- round(rowMeans(fc_matrix), digits = 4)

      if (ncol(fc_matrix) > 1) {
        sd_col <- paste0(comp_label, "_log2FC_sd")
        result_df[[sd_col]] <- round(
          apply(fc_matrix, 1, sd), digits = 4)
      }
    }
  }

  fold_change_results[["per_sgRNA"]] <- result_df

  ## per-gene summary: mean across sgRNAs for all numeric columns
  fold_change_results[["per_gene"]] <- result_df |>
    dplyr::group_by(.data$Gene) |>
    dplyr::reframe(
      dplyr::across(where(is.numeric), ~ round(mean(.x, na.rm = TRUE),
                                               digits = 4))
    ) |>
    dplyr::distinct()

  ## --- write output -------------------------------------------------------

  readr::write_tsv(
    fold_change_results[["per_sgRNA"]],
    file = output_fname_fc_sgrna,
    col_names = TRUE, quote = "none")
  message(glue::glue(
    "Fold change per sgRNA output file written: {output_fname_fc_sgrna}\n"))

  readr::write_tsv(
    fold_change_results[["per_gene"]],
    file = output_fname_fc_gene,
    col_names = TRUE, quote = "none")
  message(glue::glue(
    "Fold change per gene output file written: {output_fname_fc_gene}\n"))

  return(fold_change_results)

}
