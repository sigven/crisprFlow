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
#' @param norm_method Normalization method to use ("median", "total", or "
#' control")
#' @param mageck_bin Path to the MAGeCK binary
#' @param threads Number of threads to use (default is 1)
#' @param permutation_round Number of permutation rounds (default is 2, recommended 10 (longer time))
#' @param late_time_start Time point to start considering as "late" (default is
#' 10)
#' @param adjust_method Method for p-value adjustment (default is "fdr")
#' @param output_prefix Prefix for output files (default is "crisprFlowRun")
#' @param control_pattern Pattern to identify control genes (default is "TARGETING_NEG_CTRL")
#' @return A data frame containing the gene summary results from MAGeCK MLE
#'
#' @export
#'
detect_hits_mageck_mle <- function(
    samplesheet_csv = NULL,
    count_table_txt = NULL,
    output_dir = NULL,
    norm_method = c("median", "total", "control"),
    mageck_bin = "/Users/sigven/miniconda3/bin/mageck",
    threads = 1,
    permutation_round = 4,
    late_time_start = 10,
    adjust_method = "fdr",
    output_prefix = "mageck_mle_hits",
    control_pattern = "TARGETING_NEG_CTRL"){

  assertthat::assert_that(
    norm_method %in% c("median","none","total")
  )
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

  control_genes_fname <-
    file.path(system.file(
      'extdata', 'control_genes',
      package='crisprFlow'), 'control_genes.txt')
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

  if(!is.null(control_genes_fname) & norm_method == "control"){
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
#' @param control_pattern Pattern to identify control genes (default is "TARGETING_NEG_CTRL")
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
    control_pattern = "TARGETING_NEG_CTRL"){

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



#' Estimate log2 fold change of sgRNA counts between drug and control
#'
#' This function estimates the log2 fold change of sgRNA counts
#' between drug and control conditions at various timepoints.
#' It normalizes the counts to counts per million (CPM)
#' and computes the log2 fold change for each timepoint.
#'
#' @param counts A data frame containing raw sgRNA counts with columns:
#' 'sgRNA', 'Gene', and count columns for different conditions
#' (e.g., 'T4_DMSO', 'T4_FGF401', etc.).
#' @param control_name A character string specifying the control condition name
#' (default is "DMSO").
#' @param output_dir Directory where output files will be saved (not used in current implementation).
#' @param output_prefix Prefix for output files (not used in current implementation).
#' @param overwrite Boolean indicating whether to overwrite existing files
#' @param drug_name A character string specifying the drug condition name
#' (default is "FGF401").
#' @return A data frame with normalized counts and log2 fold change columns
#' for each timepoint.
#' @export
#'
estimate_fold_change <- function(
    counts = NULL,
    control_name = "DMSO",
    output_dir = NULL,
    output_prefix = "FGF401",
    overwrite = FALSE,
    drug_name = "FGF401"){

  ## initialize results list
  fold_change_results <- list()
  for(cat in c('per_sgRNA','per_gene')){
    fold_change_results[[cat]] <- data.frame()
  }

  assertthat::assert_that(
    !is.null(output_dir) & dir.exists(output_dir),
    msg = "Please provide path to output_dir"
  )

  output_fname_fc_sgrna <-
    file.path(
      output_dir,
      glue::glue("{output_prefix}_fc_sgRNA.tsv"))
  output_fname_fc_gene <-
    file.path(
      output_dir,
      glue::glue("{output_prefix}_fc_gene.tsv"))
  if(!overwrite &
     file.exists(output_fname_fc_sgrna) &
     file.exists(output_fname_fc_gene)){
    message(glue::glue("Fold change output files already exist:\n",
                       "{output_fname_fc_sgrna}\n",
                       "{output_fname_fc_gene}\n"))
    fold_change_results[['per_sgRNA']] <-
      as.data.frame(readr::read_tsv(
        output_fname_fc_sgrna,
        show_col_types = F,
        na = c(".","", "NA")))
    fold_change_results[['per_gene']] <-
      as.data.frame(readr::read_tsv(
        output_fname_fc_gene,
        show_col_types = F,
        na = c(".","", "NA")))
    return(fold_change_results)
  }


  assertthat::assert_that(
    !is.null(counts),
    msg = "Please provide counts data frame"
  )
  assertthat::assert_that(
    is.data.frame(counts),
    msg = "counts must be a data frame"
  )
  assertable::assert_colnames(
    counts,
    c("sgRNA","Gene","T0"),
    only_colnames = F,
    quiet = T
  )

  df <- counts

  ## Raw counts differ by sequencing depth, so first normalize —
  ## typically to counts per million (CPM)
  norm_counts <- df |>
    dplyr::mutate(
      dplyr::across(-c(sgRNA, Gene), ~ .x / sum(.x) * 1e6))

  drug_cols <- grep(
    drug_name, colnames(norm_counts), value = TRUE)
  control_cols <- grep(
    control_name, colnames(norm_counts), value = TRUE)
  assertthat::assert_that(
    length(drug_cols) > 0,
    msg = glue::glue("No columns found for drug: {drug_name}")
  )
  assertthat::assert_that(
    length(control_cols) > 0,
    msg = glue::glue("No columns found for control: {control_name}")
  )
  # Extract timepoint suffixes shared between drug and control columns
  timepoints <- intersect(
    stringr::str_remove(
      drug_cols, paste0("_",drug_name,"$")),
    stringr::str_remove(
      control_cols, paste0("_",control_name,"$"))
  )

  if(length(timepoints) == 0){
    stop("No matching timepoints found between drug and control columns")
  }

  logfc_df <- norm_counts

  ## Calculate log2 fold change for each timepoint
  ## drug vs control and control vs T0
  for (tp in timepoints) {
    drug_col <- paste0(tp,"_",drug_name)
    ctrl_col <- paste0(tp,"_",control_name)
    t0_col <- "T0"
    if(!(drug_col %in% colnames(logfc_df)) |
       !(ctrl_col %in% colnames(logfc_df))){
      next
    }
    fc_col <- paste0(tp,"_drug_vs_control_FC")

    ## add 1 to avoid log2(0)

    ## Drug vs Control
    logfc_df[[fc_col]] <- round(
      log2((logfc_df[[drug_col]] + 1) /
             (logfc_df[[ctrl_col]] + 1)), digits = 4)

    ## Control vs T0
    fc_ctrl_col <- paste0(tp,"_control_vs_T0_FC")
    logfc_df[[fc_ctrl_col]] <- round(
      log2((logfc_df[[ctrl_col]] + 1) /
             (logfc_df[[t0_col]] + 1)), digits = 4)

    ## Drug vs T0
    fc_drug_col <- paste0(tp,"_drug_vs_T0_FC")
    logfc_df[[fc_drug_col]] <- round(
      log2((logfc_df[[drug_col]] + 1) /
             (logfc_df[[t0_col]] + 1)), digits = 4)

  }

  ## round normalized sgRNA counts to 4 decimal places
  for(col in colnames(logfc_df)){
    if(endsWith(col,"_FC") |
       col == "sgRNA" |
       col == "Gene"){
      next
    }
    if(is.numeric(logfc_df[[col]])){
      logfc_df[[col]] <-
        round(logfc_df[[col]], digits = 4)
    }
  }

  fold_change_results[['per_sgRNA']] <- logfc_df

  cols_to_summarise <-
    c("sgRNA","Gene",
      grep(colnames(logfc_df),
           pattern = "_FC$",
           value = TRUE))

  fold_change_results[['per_gene']] <- logfc_df |>
    dplyr::group_by(Gene) |>
    dplyr::reframe(
      dplyr::across(
        -c(sgRNA), ~ round(mean(.x), digits = 4)
      )
    ) |>
    dplyr::distinct()

  ## write output files
  readr::write_tsv(
    fold_change_results[['per_sgRNA']],
    file = output_fname_fc_sgrna,
    col_names = TRUE, quote = "none")
  message(glue::glue("Fold change per sgRNA output file written: ",
                     "{output_fname_fc_sgrna}\n"))
  readr::write_tsv(
    fold_change_results[['per_gene']],
    file = output_fname_fc_gene,
    col_names = TRUE, quote = "none")
  message(glue::glue("Fold change per gene output file written: ",
                     "{output_fname_fc_gene}\n"))

  return(fold_change_results)

}


