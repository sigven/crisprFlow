#' Count sgRNA from FASTQ files using MAGECK
#'
#' This function counts sgRNA sequences from FASTQ files
#' using the MAGECK tool. It requires a sample sheet CSV file
#' containing sample information, the sgRNA library to use,
#' and parameters for output directory and file naming.
#'
#' @param samplesheet_csv A character string specifying the path to a CSV file
#' containing sample information with mandatory columns: fastq, sample_label, condition.
#' @param sgRNA_library A character string specifying the sgRNA library to use.
#' Must be one of 'ACOC','DTKP','GEEX','MEPR','PROT','TMMO'. Default is "DTKP".
#' @param mageck_bin A character string specifying the path
#' to the MAGECK binary.
#' @param norm_method Method for normalization ('none','median','total')
#' (default is "median").
#' @param count_N A logical indicating whether to count 'N' bases
#' in sgRNA sequences (default is FALSE).
#' @param output_dir A character string specifying the directory
#' where the count results will be saved.
#' @param output_prefix A character string to prefix output files.
#'                      Default is "FGFR4".
#' @param overwrite A logical indicating whether to overwrite existing output
#' files (default is FALSE).
#' @return None. The function performs counting and saves results to files.
#' @export
#'
count_sgRNA_mageck <- function(
    samplesheet_csv = NULL,
    sgRNA_library = "DTKP",
    mageck_bin = "/Users/sigven/miniconda3/bin/mageck",
    norm_method = "median",
    count_N = FALSE,
    output_dir = NULL,
    output_prefix = "FGFR4",
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(output_dir),
    msg = "Please provide output_dir"
  )
  assertthat::assert_that(
    dir.exists(output_dir),
    msg = glue::glue("output_dir not found at: '{output_dir}'")
  )

  if(overwrite == FALSE &
     file.exists(file.path(output_dir,
                           glue::glue("{output_prefix}.count.txt"))) &
     file.exists(file.path(output_dir,
                           glue::glue("{output_prefix}.countsummary.txt")))){
    message("Count file and summary file already exist - returning existing count data")
    count_results <- list()
    count_results[['counts']] <- data.frame()
    count_results[['stats']] <- data.frame()

    counts <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{output_prefix}.count.txt")),
      show_col_types = F))
    stats <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{output_prefix}.countsummary.txt")),
      show_col_types = F))
    count_results <- list()
    count_results$stats <- stats
    count_results$counts <- counts
    return(count_results)

  }

  assertthat::assert_that(
    is.character(sgRNA_library)
  )
  if(!sgRNA_library %in% c('ACOC','DTKP','GEEX','MEPR',
                           'PROT','TMMO')){
    msg = paste0(
      "sgRNA_library must be one of ",
      "'ACOC','DTKP','GEEX','MEPR',",
      "'PROT','TMMO'")
  }
  assertthat::assert_that(
    norm_method %in% c("median","none","total")
  )
  assertthat::assert_that(
    !is.null(mageck_bin),
    msg = "Please provide path to mageck binary"
  )
  assertthat::assert_that(
    file.exists(mageck_bin),
    msg = glue::glue("mageck binary not found at: ", mageck_bin
    )
  )
  if(is.null(samplesheet_csv) | !file.exists(samplesheet_csv)){
    stop("Please provide path to sample_sheet CSV file")
  }
  sample_data <-
    readr::read_csv(file = samplesheet_csv,
                    show_col_types = F, na = c(".","", "NA"))
  if(!all(c('fastq','sample_label',
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

  sgRNA_libfile_csv <- file.path(
    system.file(
      "extdata", "sgRNA_library", package = "crisprFlow"),
    sgRNA_library,
    glue::glue("sgRNA_library_{sgRNA_library}.csv")
  )

  assertthat::assert_that(
    file.exists(sgRNA_libfile_csv),
    msg = glue::glue("sgRNA library file not found at: '{sgRNA_libfile_csv}'")
  )

  T0_label <-
    sample_data[sample_data$T0 == 1,]$sample_label
  fastq_files <-
    paste(sample_data$fastq, collapse=" ")
  sample_labels <-
    paste(sample_data$sample_label, collapse=",")
  timepoints <- as.character(
    glue::glue('T{unique(sample_data$time_point)}'))
  treatments <- unique(sample_data$condition)
  treatments <- treatments[!is.na(treatments)]

  if(!"T0" %in% timepoints){
    msg = "T0 must be included in timepoint vector"
  }

  assertthat::assert_that(
    length(timepoints) >= 1,
    length(treatments) >= 1,
    length(sgRNA_library) == 1
  )

  mageck_count_bin <- glue::glue("{mageck_bin} count")
  if(count_N){
    mageck_count_bin <- glue::glue("{mageck_bin} count --count-n")
  }

  mageck_count_command <-
    glue::glue(
      "{mageck_count_bin} -l {sgRNA_libfile_csv} ",
      "--day0-label {T0_label} ",
      "--norm-method {norm_method} ",
      "--sample-label {sample_labels} ",
      "-n {output_dir}/{output_prefix} ",
      "--fastq {fastq_files} > {output_dir}/{output_prefix}.mageck_count.log 2>&1"
    )
  system(mageck_count_command)

  ## delete Rmarkdown templates from MAGeCK - not used
  ## delete redunant log file
  system(glue::glue(
    "rm -f {output_dir}/{output_prefix}_countsummary.R*"))
  system(glue::glue(
    "rm -f {output_dir}/{output_prefix}.count_report.Rmd"))
  system(glue::glue(
    "rm -f {output_dir}/{output_prefix}.mageck_count.log"))

  count_results <- list()
  count_results[['counts']] <- data.frame()
  count_results[['stats']] <- data.frame()

  if(!file.exists(file.path(
    output_dir,
    glue::glue("{output_prefix}.count.txt")))){
    stop("Counting failed - output count file not found.")
  } else {
    message("Counting completed.")
    counts <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{output_prefix}.count.txt")),
      show_col_types = F))
    stats <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{output_prefix}.countsummary.txt")),
      show_col_types = F))
    count_results <- list()
    count_results$stats <- stats
    count_results$counts <- counts

  }
  return(count_results)


}


#' Count sgRNA from FASTQ files using guide-counter
#'
#' This function counts sgRNA sequences from FASTQ files
#' using the guide-counter tool. It requires the path to
#' the guide-counter binary, the directory containing
#' the FASTQ files, and parameters for timepoints, treatments
#' and sgRNA library.
#'
#' @param samplesheet_csv A character string specifying the path to a CSV file
#' @param sgRNA_library A character string specifying the sgRNA library to use.
#' Must be one of 'ACOC','DTKP','GEEX','MEPR','PROT','TMMO'. Default is "DTKP".
#' @param min_sgRNA_length An integer specifying the minimum sgRNA length to consider.
#' Default is 17.
#' @param max_sgRNA_length An integer specifying the maximum sgRNA length to consider
#' (default is 25).
#' @param guide_counter_bin A character string specifying the path
#' to the guide-counter binary.
#' @param count_exact A logical indicating whether to use
#' exact matching in guide-counter (default is TRUE).
#' @param offset_sample_size An integer specifying the number of
#' reads to be examined for offset estimation (default is 100K).
#' @param offset_min_fraction A numeric value specifying the minimum fraction
#' of reads supporting an offset to consider it valid (default is 0.0025).
#' @param output_dir A character string specifying the directory
#' where the count results will be saved.
#' @param output_prefix A character string to prefix output files.
#'                      Default is "FGFR4".
#' @param overwrite A logical indicating whether to overwrite existing output
#' files (default is FALSE).
#' @return None. The function performs counting and saves results to files.
#' @export
#'
count_sgRNA_guidecounter <- function(
    samplesheet_csv = NULL,
    sgRNA_library = "DTKP",
    min_sgRNA_length = 17,
    max_sgRNA_length = 25,
    guide_counter_bin =
      "/Users/sigven/miniconda3/bin/guide-counter",
    count_exact = TRUE,
    offset_sample_size = "100000",
    offset_min_fraction = 0.0025,
    output_dir = NULL,
    output_prefix = "FGFR4",
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(output_dir),
    msg = "Please provide output_dir"
  )
  assertthat::assert_that(
    dir.exists(output_dir),
    msg = glue::glue("output_dir not found at: '{output_dir}'")
  )

  if(overwrite == FALSE &
     file.exists(
       file.path(
         output_dir,
         glue::glue("{output_prefix}.count.txt"))) &
     file.exists(
       file.path(
         output_dir,
         glue::glue("{output_prefix}.countsummary.txt")))){
    message("Count file and summary file already exist - returning existing count data")
    count_results <- list()
    count_results[['counts']] <- data.frame()
    count_results[['stats']] <- data.frame()

    counts <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{output_prefix}.count.txt")),
      show_col_types = F))
    stats <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{output_prefix}.countsummary.txt")),
      show_col_types = F))
    count_results <- list()
    count_results$stats <- stats
    count_results$counts <- counts
    return(count_results)

  }


  if(is.null(samplesheet_csv) | !file.exists(samplesheet_csv)){
    stop("Please provide path to sample_sheet CSV file")
  }
  sample_data <-
    readr::read_csv(
      file = samplesheet_csv,
      show_col_types = F, na = c(".","", "NA"))
  if(!all(c('fastq','sample_label',
            'condition','T0','time_point') %in%
          colnames(sample_data))){
    stop(
      paste0("sample_sheet must contain columns: ",
             "fastq, sample_label, condition, T0, time_point"))
  }

  T0_label <-
    unique(sample_data[sample_data$T0 == 1,]$sample_label)
  fastq_files <-
    paste(sample_data$fastq, collapse=" ")
  sample_labels <-
    paste(sample_data$sample_label, collapse=",")
  sample_labels <- sample_data$sample_label
  if("replicate" %in% colnames(sample_data)){
    sample_labels <- paste(
      sample_labels, sample_data$replicate, sep="_")
  }
  sample_data$sample_label2 <- sample_labels

  timepoints <- as.character(
    glue::glue('T{unique(sample_data$time_point)}'))
  treatments <- unique(sample_data$condition)
  treatments <- treatments[!is.na(treatments)]

  if(!"T0" %in% timepoints){
    msg = "T0 must be included in timepoint vector"
  }

  fastq_dict <- list()
  for(s in sample_data$sample_label2){
    fastq_dict[[s]] <-
      sample_data[sample_data$sample_label2 == s,]$fastq
  }

  assertthat::assert_that(
    is.character(timepoints),
    is.character(treatments),
    is.character(sgRNA_library)
  )
  assertthat::assert_that(
    length(timepoints) >= 1,
    length(treatments) >= 1,
    length(sgRNA_library) == 1
  )
  assertthat::assert_that(
    !is.null(guide_counter_bin),
    msg = "Please provide path to guide_counter binary"

  )

  if(!sgRNA_library %in% c('ACOC','DTKP','GEEX','MEPR',
                           'PROT','TMMO')){
    msg = paste0(
      "sgRNA_library must be one of ",
      "'ACOC','DTKP','GEEX','MEPR',",
      "'PROT','TMMO'")
  }

  if(!"T0" %in% timepoints){
    msg = "T0 must be included in timepoint vector"
  }

  guide_counter_command_basic <-
    glue::glue(
      "{guide_counter_bin} count --control-pattern _CTRL ",
      "--offset-min-fraction {offset_min_fraction} ",
      "--offset-sample-size {offset_sample_size} "
    )
  if(count_exact){
    guide_counter_command_basic <-
      glue::glue(
        "{guide_counter_command_basic} --exact-match "
      )
  }

  # for reproducibility
  set.seed(12345)

  # Pool of characters: uppercase, lowercase, digits
  chars <- c(letters, LETTERS, 0:9)

  # Sample 9 characters and collapse into one string
  rand_str <- paste0(sample(chars, 9, replace = TRUE), collapse = "")
  message(glue::glue("Random string for output files: {rand_str}"))

  ## Construct commands for each timepoint and treatment
  ## Iterate over all sgRNA lengths in library (from 17 to 25)
  commands <- list()
  i <- min_sgRNA_length
  while(i <= max_sgRNA_length){
    ## sgRNA library file for this length
    sgRNA_libfile_csv <- file.path(
      system.file(
        "extdata", "sgRNA_library", package = "crisprFlow"),
      sgRNA_library,
      glue::glue("sgRNA_library_{sgRNA_library}_{i}.csv")
    )
    if(file.exists(sgRNA_libfile_csv)){

      # commands[['T0']] <-
      #   glue::glue(
      #     "{guide_counter_command_basic}",
      #     " --input {fastq_dict[['T0']]}",
      #     " --library {sgRNA_libfile_csv}",
      #     " --output {output_dir}/{output_prefix}_{rand_str}_T0_{i}",
      #     " > /dev/null 2>&1"
      #   )
      # system(commands[['T0']])

      for(label in sample_labels){
      #for(ti in setdiff(timepoints,"T0")){
      #  for(tr in treatments){
          #label <- paste0(ti,'_',tr)
          if(!label %in% names(fastq_dict)){
            next
          }
          cmd = glue::glue(
            "{guide_counter_command_basic}",
            " --input {fastq_dict[[label]]}",
            " --library {sgRNA_libfile_csv}",
            " --output {output_dir}/{output_prefix}_{rand_str}_{label}_{i}",
            " > /dev/null 2>&1")
          system(cmd)
        #}
      }
    } else {
      message(glue::glue("sgRNA library file not found at: {sgRNA_libfile_csv}"))
    }
    i <- i + 1
  }


  ## Execute commands
  stats_all_samples <- data.frame()
  counts_all_samples <- data.frame()
  t0_seen <- FALSE

  for(label in sample_labels){
  #for(ti in timepoints){
    ## don't do T0 two times
    #for(tr in treatments){
      #label <- glue::glue("{ti}_{tr}")
      #if(ti == 'T0'){
      #  label <- 'T0'
      #}

      #if(ti == 'T0' & t0_seen){
      #  next
      #}
      #t0_seen <- TRUE

      stats <- data.frame()
      counts <- data.frame()
      i <- min_sgRNA_length
      while(i <= max_sgRNA_length){
        count_file_i <-
          glue::glue(
            "{output_dir}/{output_prefix}_{rand_str}_{label}_{i}.counts.txt")
        if(file.exists(count_file_i)){
          counts_i <- readr::read_tsv(
            count_file_i, show_col_types = F) |>
            as.data.frame()
          if(NROW(counts_i) > 0){
              colnames(counts_i) <- c('sgRNA','Gene',label)
              counts <- counts |>
                dplyr::bind_rows(counts_i)

              stats_file_i <-
                glue::glue(
                  "{output_dir}/{output_prefix}_{rand_str}_{label}_{i}.stats.txt")
              stats_i <- readr::read_tsv(
                stats_file_i, show_col_types = F)
              stats_i$label <- label
              stats <- stats |>
                dplyr::bind_rows(stats_i)
          }
        }

        i <- i + 1
      }
      #cat(label, "\n")
      stats_summarised <- as.data.frame(
        stats |>
          dplyr::select(
            -dplyr::any_of(
              c("mean_reads_per_guide",
               "mean_reads_essential",
               "mean_reads_nonessential",
               "mean_reads_control",
               "mean_reads_other")
            )
          ) |>
          dplyr::rename(
            File = file,
            Label = label,
            Reads = total_reads,
          ) |>
          dplyr::group_by(
            .data$File,
            .data$Label,
            .data$Reads
          ) |>
          dplyr::reframe(
            Mapped = sum(.data$mapped_reads),
            Percentage = round(
              sum(.data$mapped_reads) / .data$Reads, digits=4),
            TotalsgRNAs = sum(.data$total_guides),
            Zerocounts = sum(.data$zero_read_guides)
          ) |>
          dplyr::distinct())

      stats_all_samples <- stats_all_samples |>
        dplyr::bind_rows(stats_summarised) |>
        dplyr::distinct()

      if(NROW(counts_all_samples) == 0){
        counts_all_samples <- counts
      }else{
        counts_all_samples <- counts_all_samples |>
          dplyr::left_join(counts, by = c("sgRNA","Gene")) |>
          dplyr::distinct()
      }

    #}
  }

  gini_indices <- data.frame()
  for(label in sample_labels)
  #for(label in sample_data$sample_label){
    if(label %in% colnames(counts_all_samples)){
      gini_index <- data.frame(
        Label = label,
        GiniIndex = crisprFlow::gini_index2(
          x = counts_all_samples[[label]],
          log_transform = TRUE,
          na.rm = TRUE)
      )
      gini_indices <- dplyr::bind_rows(
        gini_indices, gini_index)
    #}
  }

  stats_all_samples <-
    stats_all_samples |>
    dplyr::left_join(gini_indices, by = "Label") |>
    dplyr::distinct()


  readr::write_tsv(
    as.data.frame(stats_all_samples),
    file = file.path(
      output_dir,
      glue::glue("{output_prefix}.countsummary.txt")),
    col_names = TRUE, quote = "none")

  readr::write_tsv(
    as.data.frame(counts_all_samples),
    file = file.path(
      output_dir,
      glue::glue("{output_prefix}.count.txt")),
    col_names = TRUE, quote = "none")

  system(glue::glue("rm -f {output_dir}/{output_prefix}_{rand_str}*.txt"))
  message("Counting completed.")

  count_results <- list()
  count_results$stats <- as.data.frame(stats_all_samples)
  count_results$counts <- as.data.frame(counts_all_samples)
  return(count_results)
}


#' Count sgRNA from FASTQ files using guide-counter (alternative version)
#'
#' This function counts sgRNA sequences from FASTQ files using the
#' guide-counter tool. It requires the path to the guide-counter binary,
#' the directory containing the FASTQ files, and other parameters for counting.
#'
#' @param fname_fastq A character string specifying the path to the FASTQ file to be counted.
#' @param sample_label A character string specifying the label for the sample being counted.
#' @param sgRNA_library A character string specifying the sgRNA library to use. Must
#' be one of 'ACOC','DTKP','GEEX','MEPR','PROT','TMMO'. Default is "DTKP".
#' @param min_sgRNA_length An integer specifying the minimum sgRNA length to consider.
#' Default is 17.
#' @param max_sgRNA_length An integer specifying the maximum sgRNA length to
#' consider (default is 25).
#' @param guide_counter_bin A character string specifying the path to the
#' guide-counter binary.
#' @param count_exact A logical indicating whether to use exact matching in
#' guide-counter (default is TRUE).
#' @param offset_sample_size An integer specifying the number of reads to be
#' examined for offset estimation (default is 100K).
#' @param offset_min_fraction A numeric value specifying the minimum fraction
#' of reads supporting an offset to consider it valid (default is 0.0025).
#' @param output_dir A character string specifying the directory where the
#' count results will be saved.
#' @param output_prefix A character string to prefix output files. Default is "FGFR4".
#' @param overwrite A logical indicating whether to overwrite existing
#' output files (default is FALSE).
#'
#' @return A list containing two data frames: 'stats' with summary
#' statistics and 'counts' with sgRNA counts.
#'
#' @export
#'
count_sgRNA_guidecounter_fastq <- function(
    fname_fastq = NULL,
    sample_label = NULL,
    sgRNA_library = "DTKP",
    min_sgRNA_length = 17,
    max_sgRNA_length = 25,
    guide_counter_bin =
      "/Users/sigven/miniconda3/bin/guide-counter",
    count_exact = TRUE,
    offset_sample_size = "100000",
    offset_min_fraction = 0.0025,
    output_dir = NULL,
    overwrite = FALSE){

  assertthat::assert_that(
    !is.null(output_dir),
    msg = "Please provide output_dir"
  )
  assertthat::assert_that(
    dir.exists(output_dir),
    msg = glue::glue("output_dir not found at: '{output_dir}'")
  )

  if(overwrite == FALSE &
     file.exists(
       file.path(
         output_dir,
         glue::glue("{sample_label}.count.txt"))) &
     file.exists(
       file.path(
         output_dir,
         glue::glue("{sample_label}.countsummary.txt")))){
    message("Count file and summary file already exist - returning existing count data")
    count_results <- list()
    count_results[['counts']] <- data.frame()
    count_results[['stats']] <- data.frame()

    counts <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{sample_label}.count.txt")),
      show_col_types = F))
    stats <- as.data.frame(readr::read_tsv(
      file.path(
        output_dir,
        glue::glue("{sample_label}.countsummary.txt")),
      show_col_types = F))
    count_results <- list()
    count_results$stats <- stats
    count_results$counts <- counts
    return(count_results)

  }

  if(is.null(fname_fastq) | !file.exists(fname_fastq)){
    stop("Please provide path to FASTQ file")
  }
  if(is.null(sample_label)){
    stop("Please provide sample_label for this FASTQ file")
  }
  assertthat::assert_that(
    !is.null(guide_counter_bin),
    msg = "Please provide path to guide_counter binary"

  )

  if(!sgRNA_library %in% c('ACOC','DTKP','GEEX','MEPR',
                           'PROT','TMMO')){
    msg = paste0(
      "sgRNA_library must be one of ",
      "'ACOC','DTKP','GEEX','MEPR',",
      "'PROT','TMMO'")
  }

  guide_counter_command_basic <-
    glue::glue(
      "{guide_counter_bin} count --control-pattern _CTRL ",
      "--offset-min-fraction {offset_min_fraction} ",
      "--offset-sample-size {offset_sample_size} "
    )
  if(count_exact){
    guide_counter_command_basic <-
      glue::glue(
        "{guide_counter_command_basic} --exact-match "
      )
  }

  # for reproducibility
  set.seed(12345)

  # Pool of characters: uppercase, lowercase, digits
  chars <- c(letters, LETTERS, 0:9)

  # Sample 9 characters and collapse into one string
  rand_str <- paste0(sample(chars, 9, replace = TRUE), collapse = "")
  message(glue::glue("Random string for output files: {rand_str}"))

  ## Construct commands for each timepoint and treatment
  ## Iterate over all sgRNA lengths in library (from 17 to 25)
  commands <- list()
  i <- min_sgRNA_length
  while(i <= max_sgRNA_length){
    ## sgRNA library file for this length
    sgRNA_libfile_csv <- file.path(
      system.file(
        "extdata", "sgRNA_library", package = "crisprFlow"),
      sgRNA_library,
      glue::glue("sgRNA_library_{sgRNA_library}_{i}.csv")
    )
    if(file.exists(sgRNA_libfile_csv)){
      cmd = glue::glue(
        "{guide_counter_command_basic}",
        " --input {fname_fastq}",
        " --library {sgRNA_libfile_csv}",
        " --output {output_dir}/{sample_label}_{rand_str}_{i}",
        " > /dev/null 2>&1")
      system(cmd)
    } else {
      message(
        glue::glue("sgRNA library file not found at: {sgRNA_libfile_csv}"))
    }
    i <- i + 1
  }


  stats <- data.frame()
  counts <- data.frame()
  i <- min_sgRNA_length
  while(i <= max_sgRNA_length){
    count_file_i <-
      glue::glue(
        "{output_dir}/{sample_label}_{rand_str}_{i}.counts.txt")
    if(file.exists(count_file_i)){
      counts_i <- readr::read_tsv(
        count_file_i, show_col_types = F) |>
        as.data.frame()
      if(NROW(counts_i) > 0){
        colnames(counts_i) <- c('sgRNA','Gene',sample_label)
        counts <- counts |>
          dplyr::bind_rows(counts_i)

        stats_file_i <-
          glue::glue(
            "{output_dir}/{sample_label}_{rand_str}_{i}.stats.txt")
        stats_i <- readr::read_tsv(
          stats_file_i, show_col_types = F)
        stats_i$label <- sample_label
        stats <- stats |>
          dplyr::bind_rows(stats_i)
      }
    }

    i <- i + 1
  }

  stats_summarised <- as.data.frame(
    stats |>
      dplyr::select(
        -dplyr::any_of(
          c("mean_reads_per_guide",
            "mean_reads_essential",
            "mean_reads_nonessential",
            "mean_reads_control",
            "mean_reads_other")
        )
      ) |>
      dplyr::rename(
        File = file,
        Label = label,
        Reads = total_reads,
      ) |>
      dplyr::group_by(
        .data$File,
        .data$Label,
        .data$Reads
      ) |>
      dplyr::reframe(
        Mapped = sum(.data$mapped_reads),
        Percentage = round(
          sum(.data$mapped_reads) / .data$Reads, digits=4),
        TotalsgRNAs = sum(.data$total_guides),
        Zerocounts = sum(.data$zero_read_guides)
      ) |>
      dplyr::distinct())

  gini_index <- data.frame(
    Label = sample_label,
    GiniIndex = crisprFlow::gini_index2(
      x = counts[[sample_label]],
      log_transform = TRUE,
      na.rm = TRUE)
  )

  stats_summarised <- stats_summarised |>
    dplyr::left_join(
      gini_index, by = c("Label")) |>
    dplyr::distinct()



  readr::write_tsv(
    as.data.frame(stats_summarised),
    file = file.path(
      output_dir,
      glue::glue("{sample_label}.countsummary.txt")),
    col_names = TRUE, quote = "none")

  readr::write_tsv(
    as.data.frame(counts),
    file = file.path(
      output_dir,
      glue::glue("{sample_label}.count.txt")),
    col_names = TRUE, quote = "none")

  system(glue::glue("rm -f {output_dir}/{sample_label}_{rand_str}*.txt"))
  message("Counting completed.")

  count_results <- list()
  count_results$stats <- as.data.frame(stats_summarised)
  count_results$counts <- as.data.frame(counts)
  return(count_results)
}


#' Check for expected FASTQ files in a directory
#'
#' This function checks if all expected FASTQ files are present
#' in the specified directory based on provided timepoints and treatments.
#'
#' @param timepoints A character vector of timepoints to check
#' (default is c("T0","T4","T10","T14")).
#' @param treatments A character vector of treatments to check
#' (default is c("DMSO","FGF401")).
#' @param fastq_dir A character string specifying the directory
#' containing the FASTQ files.
#' @return None. The function stops execution if any expected
#' FASTQ files are missing.
#' @export
#'
check_fastq_files <- function(
    timepoints = c("T0","T4","T10","T14"),
    treatments = c("DMSO","FGF401"),
    fastq_dir = NULL){

  assertthat::assert_that(
    is.character(timepoints),
    is.character(treatments)
  )
  assertthat::assert_that(
    dir.exists(fastq_dir),
    msg = glue::glue("fastq_dir not found at: ", fastq_dir)
  )

  expected_fastq_files <- c()
  for(tp in timepoints){
    for(tr in treatments){
      if(tp == "T0"){
        expected_fastq_files <- c(
          expected_fastq_files, "T0.fastq.gz")
      } else {
        expected_fastq_files <- c(
          expected_fastq_files,
          glue::glue("{tp}_{tr}.fastq.gz")
        )
      }
    }
  }

  expected_fastq_files <-
    unique(expected_fastq_files)
  actual_fastq_files <-
    list.files(fastq_dir)
  missing_fastq_files <-
    setdiff(expected_fastq_files, actual_fastq_files)

  if(length(missing_fastq_files) > 0){
    stop(paste("Missing FASTQ files:",
               paste(missing_fastq_files, collapse = ", ")))
  } else {
    message("All expected FASTQ files are present.")
  }

}

#' Compute the Gini coefficient/index for a numeric vector x.
#'
#' This function calculates the Gini index for a numeric vector.
#' The Gini index is a measure of statistical dispersion that
#' represents the inequality among values of a frequency distribution.
#' It ranges from 0 (perfect equality) to 1 (perfect inequality).
#'
#' @param x A numeric vector for which to calculate the Gini index.
#' @param na.rm A logical indicating whether to remove NA values
#' from the input vector (default is TRUE).
#' @param normalize_counts A logical indicating whether to normalize
#' the counts to counts per million and log2-transform them
#' (default is TRUE).
#' @param log_transform A logical indicating whether to log2-transform
#' @return A numeric value representing the Gini index of the input vector.
#' @export
gini_index <- function(x,
                       na.rm = TRUE,
                       normalize_counts = TRUE,
                       log_transform = TRUE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  x <- as.numeric(x)
  if(normalize_counts){
    x <- (x / sum(x)) * 1e6
  }
  if(log_transform){
    x <- log2(x + 1)
  }
  n <- length(x)
  mu <- mean(x)
  return(2 * cov(x, rank(x)) / (mu * n))
}


#' Compute the Gini coefficient/index for a numeric vector x.
#'
#' Calculate Gini index - MAGeCK alike (same results as MAGeCK,
#' no count normalization - only log2 transform). The Gini index is
#' a measure of statistical dispersion that represents the inequality
#' among values of a frequency distribution. It ranges from 0 (perfect
#' equality) to 1 (perfect inequality).
#'
#' @param x numeric vector (sgRNA counts). Non-negative values expected.
#' @param normalize_counts A logical indicating whether to normalize
#' @param log_transform A logical indicating whether to log2-transform
#'        the values in x (default is TRUE).
#' @param na.rm logical; remove NA pairs (x and w) before computing. Default TRUE.
#' @param allow.negative logical; if FALSE (default) stops when x has negatives.
#'        If TRUE, it will shift x by adding -min(x) (so smallest becomes 0).
#'        (Shifting changes the Gini interpretation — use with care.)
#' @param digits numeric; round to the specified number of decimal places.
#' @return numeric scalar in [0,1] (or NaN if undefined)
#' @export
#'
gini_index2 <- function(
    x,
    normalize_counts = FALSE,
    log_transform = TRUE,
    na.rm = TRUE,
    allow.negative = FALSE,
    digits = 4){

  ## Check inputs
  if (!is.numeric(x)) stop("x must be numeric.")
  if(na.rm){
    x <- x[!is.na(x)]
  }
  x <- as.numeric(x)

  ## Preprocessing of sgRNA counts
  if(normalize_counts){
    x <- (x / sum(x)) * 1e6
  }
  if(log_transform){
    x <- log2(x + 1)
  }

  n <- length(x)
  if (n == 0) return(NA_real_)

  ## Handle negatives
  if (any(x < 0)) {
    if (!allow.negative) stop("x contains negative values. Set allow.negative=TRUE to shift them (changes interpretation).")
    shift <- -min(x)
    x <- x + shift
  }

  ## Main computation
  # - If all zero or all identical => Gini = 0
  if (all(x == 0)) return(0)
  # Sort ascending
  o <- order(x)
  xs <- x[o]
  n <- length(xs)
  mu <- mean(xs)
  if (mu == 0) return(0) # all zeros -> Gini 0
  # Use the formula: G = (2 * sum_i i * x_i) / (n * sum x) - (n+1)/n
  i <- seq_len(n)
  gini_index <- (2 * sum(i * xs)) / (n * sum(xs)) - (n + 1) / n
  # numerical guard: clamp to [0,1]
  gini_index <- max(0, min(1, gini_index))

  gini_index <- round(gini_index, digits = digits)
  return(as.numeric(gini_index))

}
