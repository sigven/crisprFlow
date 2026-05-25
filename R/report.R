
#' Function to clean treatment and timepoint labels from a data frame
#' This function separates the 'label' column into 'treatment' and 'timepoint'
#' columns, removes any existing 'treatment' and 'timepoint' columns,
#' and arranges the data frame by timepoint
#'
#' @param df A data frame containing a 'label' column with treatment and timepoint information
#' @param timecourse_drug_screen Logical indicating whether
#' the data is from a timecourse drug screen (default is TRUE)
#'
#' @export
clean_treatment_timepoint_label <- function(
    df = NULL,
    timecourse_drug_screen = TRUE){

  if(timecourse_drug_screen == FALSE){
    return(df)
  }

  if("label" %in% colnames(df)){

    if("treatment" %in% colnames(df)){
      df$treatment <- NULL
    }
    if("timepoint" %in% colnames(df)){
      df$timepoint <- NULL
    }

    df <- df |>
      tidyr::separate(
        .data$label,
        into = c("timepoint","treatment"),
        sep = "_",
        remove = F,
        fill = "right") |>
      dplyr::distinct() |>
      dplyr::mutate(
        timepoint_num = as.numeric(
          stringr::str_replace(
            .data$timepoint,"T","")
        )
      ) |>
      dplyr::arrange(.data$timepoint_num) |>
      dplyr::select(-c("timepoint_num")) |>
      dplyr::mutate(
        treatment = dplyr::if_else(
          label == "T0",
          as.character("T0"),
          as.character(.data$treatment)
        )
      ) |>
      dplyr::mutate(
        timepoint = factor(
          .data$timepoint,
          levels = unique(.data$timepoint)),
        treatment = factor(
          .data$treatment, levels =
            c(setdiff(
              unique(.data$treatment),
              "T0"), "T0")))

  }

  return(df)

}

#' Read sgRNA count data from a MAGeCK/guide-counter count file
#'
#' This function reads sgRNA count data from a specified MAGeCK/guide-counter
#' count file, filters out non-targeting controls, and returns a data frame
#'
#' @param sgRNA_count_file The path to the MAGeCK/guide-counter count file
#'
#' @return A data frame containing sgRNA counts with sgRNA IDs as row names
#'
#' @export
#'
read_sgRNA_counts <- function(
    sgRNA_count_file){

  sgRNA_count_data <- as.data.frame(
    readr::read_tsv(
      file = sgRNA_count_file,
      show_col_types = F) |>
      dplyr::filter(
        !stringr::str_detect(Gene, "TARGETING")) |>
      dplyr::select(-c("Gene")))

  rownames(sgRNA_count_data) <- sgRNA_count_data$sgRNA
  sgRNA_count_data$sgRNA <- NULL

  return(sgRNA_count_data)


}

#' Read sgRNA count data from a MAGeCK/guide-counter count file
#' and separate targeted sgRNAs and controls
#'
#' @param sgRNA_count_file The path to the MAGeCK/guide-counter count file
#' @param normalize Logical indicating whether to normalize sgRNA counts
#' to counts per million (default is TRUE)
#' @return A data frame containing sgRNA counts with sgRNA IDs as row names
#'
#' @export
#'
read_sgRNA_counts2 <- function(
    sgRNA_count_file = NULL,
    normalize = T){

  sgRNA_data <- list()
  sgRNA_data[['all']] <- data.frame()
  sgRNA_data[['targeted']] <- data.frame()
  sgRNA_data[['controls']] <- data.frame()

  if(!is.null(sgRNA_count_file) &
     file.exists(sgRNA_count_file)){
    sgRNA_count_data <- as.data.frame(
      readr::read_tsv(
        file = sgRNA_count_file,
        show_col_types = F))

    if(normalize == T){
      # Normalize the sgRNA counts
      for(c in colnames(sgRNA_count_data)){
        if(c %in% c("sgRNA","Gene")){
          next
        }
        sgRNA_count_data[[c]] <-
          sgRNA_count_data[[c]] /
          sum(sgRNA_count_data[[c]], na.rm = T) * 1e6
      }
    }

    sgRNA_count_controls <-
      sgRNA_count_data |>
      dplyr::filter(
        stringr::str_detect(Gene, "NEG_CTRL"))

    sgRNA_count_targeted <-
      sgRNA_count_data |>
      dplyr::filter(
        !stringr::str_detect(Gene, "NEG_CTRL"))

    sgRNA_data[['all']] <- sgRNA_count_data
    sgRNA_data[['targeted']] <- sgRNA_count_targeted
    sgRNA_data[['controls']] <- sgRNA_count_controls

  }

  return(sgRNA_data)

}


#' Read sgRNA summary statistics from a MAGeCK/guide-counter summary file
#'
#' This function reads sgRNA summary statistics from a specified
#' MAGeCK/guide-counter summary file, cleans the column names,
#' and returns a data frame
#'
#' @param sgRNA_count_summary_file The path to the MAGeCK/guide-counter count
#' summary file
#' @param timecourse_drug_screen logical indicating whether the data is from a timecourse drug screen (default is TRUE)
#' @return A data frame containing cleaned sgRNA summary statistics
#' @export
#'
read_sgRNA_summary <- function(
    sgRNA_count_summary_file = NULL,
    timecourse_drug_screen = TRUE){

  assertthat::assert_that(
    !is.null(sgRNA_count_summary_file),
    msg = "sgRNA_count_summary_file parameter must be provided"
  )
  assertthat::assert_that(
    file.exists(sgRNA_count_summary_file),
    msg = "sgRNA_count_summary_file does not exist at the specified path"
  )

  sgRNA_count_summary <- as.data.frame(
    readr::read_tsv(
      file = sgRNA_count_summary_file,
      show_col_types = F)) |>
    janitor::clean_names() |>
    dplyr::rename(total_reads = reads,
                  total_sgrnas = totalsg_rn_as) |>
    crisprFlow::clean_treatment_timepoint_label(
      timecourse_drug_screen = timecourse_drug_screen
    )

  return(sgRNA_count_summary)

}

#' Plot PCA of sgRNA count data
#'
#' This function performs PCA on sgRNA count data and generates a scatter plot of the first two principal components.
#'
#' @param title The title of the PCA plot (default is "Title")
#' @param sgRNA_counts A data frame containing sgRNA counts with sgRNA IDs as row names
#' @param ggsci_style The style for the color palette from the ggsci
#' package (default is "JCO")
#' @param treatment_name The name of the treatment used in the experiment (default is "
#' FGF401")
#' @param plot_fontsize Font size for the plot text elements (default is 11
#'
#' @return A ggplot2 object representing the PCA scatter plot
#' @export
plot_sgRNA_pca <- function(
    title = "Title",
    sgRNA_counts,
    ggsci_style = "JCO",
    treatment_name = "FGF401",
    plot_fontsize = 11){

  assertthat::assert_that(
    !is.null(sgRNA_counts),
    msg = "sgRNA_counts parameter must be provided")
  assertthat::assert_that(
    is.data.frame(sgRNA_counts),
    msg = "sgRNA_counts must be a data frame")
  assertable::assert_colnames(
    sgRNA_counts,
    c("Gene","sgRNA"),
    only_colnames = F, quiet = TRUE
  )

  sgRNA_counts_mat <-
    sgRNA_counts |>
    dplyr::select(-c("Gene"))

  rownames(sgRNA_counts_mat) <-
    sgRNA_counts_mat$sgRNA
  sgRNA_counts_mat$sgRNA <- NULL



  ## Remove zero-variance sgRNAs
  counts_transposed <- t(sgRNA_counts_mat)
  x <- caret::nearZeroVar(
    as.matrix(counts_transposed), saveMetrics = TRUE)
  counts_transposed_noZeroVar <-
    counts_transposed[,x$zeroVar == F]

  ## Perform PCA
  pca_result <- stats::prcomp(
    counts_transposed_noZeroVar, scale. = TRUE)  # scale = TRUE standardizes variables

  pca_df <- as.data.frame(pca_result$x)  # Extract scores
  pca_df$label <- rownames(pca_df)
  rownames(pca_df) <- NULL

  ## Clean treatment and timepoint labels
  pca_df <- crisprFlow::clean_treatment_timepoint_label(
    df = pca_df)

  ## Calculate variance explained for PC1 and PC2
  pca_1_var_explained <-
    round(summary(pca_result)$importance[2,][1] *100, digits = 0)
  pca_2_var_explained <-
    round(summary(pca_result)$importance[2,][2] *100, digits = 0)

  ## Generate PCA plot
  p <- ggplot2::ggplot(
    pca_df, ggplot2::aes(
      x = PC1, y = PC2)) +
    ggplot2::geom_point(
      size = 2.7, ggplot2::aes(shape = timepoint,
                               colour = treatment)) +
    ggplot2::xlab(paste0("PC1: ",pca_1_var_explained,"%"))+
    ggplot2::ylab(paste0("PC2: ",pca_2_var_explained,"%")) +
    ggplot2::theme_classic() +
    ggsci::scale_colour_jco() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = plot_fontsize, vjust = 0.5),
      axis.text.y = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize - 1),
      axis.title.x = ggplot2::element_text(
        family = "Helvetica",
        size = plot_fontsize,
        margin = ggplot2::margin(
          b = 12)),
      axis.title.y = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize),
      legend.text = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize)
    )

  plots <- list()
  plots[['ggplot']] <- p
  plots[['plotly']] <- plotly::ggplotly(p)

  ## Adjust legend in plotly plot - make it horizontal
  ## - remove legend title
  plots[['plotly']]$x$layout$legend$title <- NULL
  plots[['plotly']] <- plotly::layout(
    plots[['plotly']],
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = -0.2
    ),
    margin = list(
      l = 80,
      r = 20,
      t = 50,
      b = 50)
  )

  return(plots)
}


plot_FC_volcano <- function(
    fc_data = NULL,
    per_gene = TRUE,
    mageck_mle_hits = NULL,
    fpa_hits = NULL){

  # Placeholder for future implementation
  assertthat::assert_that(
    !is.null(fc_df),
    msg = "fc_df parameter must be provided"
  )
  ## check that fc_data is a list containing two data.frames
  ## 'per_sgRNA' and 'per_gene'
  assertthat::assert_that(
    is.list(fc_data),
    msg = "fc_data must be a list"
  )
  assertthat::assert_that(
    all(c("per_sgRNA","per_gene") %in% names(fc_data)),
    msg = "fc_data must contain 'per_sgRNA' and 'per_gene'"
  )
  assertthat::assert_that(
    is.data.frame(fc_data$per_sgRNA),
    msg = "fc_data$per_sgRNA must be a data.frame"
  )
  assertthat::assert_that(
    is.data.frame(fc_data$per_gene),
    msg = "fc_data$per_gene must be a data.frame"
  )

  assertthat::assert_that(
    per_gene == TRUE,
    msg = "Currently only per_gene = TRUE is supported"
  )

  fc_df <- fc_data$per_gene
  assertable::assert_colnames(
    fc_df,
    "Gene",
    only_colnames = FALSE,
    quiet = TRUE
  )

  fc_df_long <- fc_df |>
    tidyr::pivot_longer(
      cols = -c(Gene),
      names_to = c("Time", "Target", "Reference"),
      names_pattern = "(T\\d+)_(.*)_vs_(.*)_FC"
    ) |>
    dplyr::mutate(
      Comparison = paste0(Target, " vs ", Reference),
      Time = factor(Time, levels = c("T4", "T10", "T14"))
    ) |>
    dplyr::rename(
      log2FC = value
    )


}


#' Plot sgRNA zero counts bar plot
#'
#' This function generates a bar plot showing the number of sgRNAs with zero counts
#'
#' @param sgRNA_summary A data frame containing sgRNA summary statistics
#' including zero counts
#' @param ggsci_style The style for the color palette from the ggsci
#' package (default is "JCO")
#' @param plot_fontsize Font size for the plot text elements (default is 12
#'
#' @return A ggplot2 object representing the sgRNA zero counts bar plot
#' @export
plot_sgRNA_missing <- function(
    sgRNA_summary,
    ggsci_style = "JCO",
    plot_fontsize = 12){

  assertthat::assert_that(
    is.data.frame(sgRNA_summary),
    msg = "sgRNA_summary must be a data frame"
  )
  assertable::assert_colnames(
    sgRNA_summary,
    c("treatment","timepoint","zerocounts"),
    only_colnames = F, quiet = TRUE
  )

  sgRNA_summary <- sgRNA_summary |>
    dplyr::mutate(
      pct_missing_guides =
        round((zerocounts / total_sgrnas) * 100, digits = 2)
    )

  p <- ggplot2::ggplot(
    data = sgRNA_summary,
    ggplot2::aes(
      x = timepoint,
      y = pct_missing_guides,
      fill = treatment,
      text =
        paste("timepoint:", timepoint,
              "<br>treatment: ", treatment,
              "<br>percent missing guides: ",pct_missing_guides))) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_dodge()) +
    ggplot2::theme_classic() +
    ggsci::scale_fill_jco() +
    ggplot2::scale_y_continuous(
      breaks = c(2, 4, 6, 8, 10, 12)) +
    ggplot2::ylab("Percent missing - zero sgRNA count") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = plot_fontsize, vjust = 0.5),
      axis.text.y = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize - 1),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize)
    )


  plots <- list()
  plots[['ggplot']] <- p
  plots[['plotly']] <- plotly::ggplotly(p, tooltip = "text")

  ## Adjust legend in plotly plot - make it horizontal
  ## - remove legend title
  plots[['plotly']]$x$layout$legend$title <- NULL
  plots[['plotly']] <- plotly::layout(
    plots[['plotly']],
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = -0.2
    )
  )


  return(plots)

}


#' Plot sgRNA mapping statistics bar plot
#'
#' This function generates a bar plot showing the mapping statistics of sgRNA reads
#'
#' @param sgRNA_summary A data frame containing sgRNA summary statistics
#' including read counts and mapping percentages
#' @param ggsci_style The style for the color palette from the ggsci
#' package (default is "JCO")
#' @param treatment_name The name of the treatment used in the experiment (default is "
#' FGF401")
#' @param control_name The name of the control used in the experiment (default is "
#' DMSO")
#' @param title_x_top_margin Top margin for the x-axis title (default is 12)
#' @param title_y_right_margin Right margin for the y-axis title (default is 12)
#' @param plot_fontsize Font size for the plot text elements (default is 12)
#'
#' @return A ggplot2 object representing the sgRNA mapping statistics bar plot
#' @export
#'
#'
plot_sgRNA_mapping <- function(
    sgRNA_summary,
    ggsci_style = "JCO",
    treatment_name = "FGF401",
    control_name = "DMSO",
    title_x_top_margin = 12,
    title_y_right_margin = 12,
    plot_fontsize = 12){

  assertthat::assert_that(
    is.data.frame(sgRNA_summary),
    msg = "sgRNA_summary must be a data frame"
  )
  assertable::assert_colnames(
    sgRNA_summary,
    c("treatment","timepoint","total_reads",
      "mapped","percentage"),
    only_colnames = F, quiet = TRUE
  )

  mapped_stats <-
    dplyr::select(
      sgRNA_summary,
      c("treatment","timepoint",
        "total_reads","mapped",
        "percentage")
    ) |>
    dplyr::rename(reads = mapped) |>
    dplyr::mutate(mapping_status = "mapped") |>
    dplyr::distinct()

  unmapped_stats <-
    mapped_stats |>
    dplyr::mutate(mapping_status = "unmapped") |>
    dplyr::mutate(reads = total_reads - reads) |>
    dplyr::mutate(percentage = 1 - percentage) |>
    dplyr::mutate(
      percentage = paste0(percentage * 100,"%"))

  mapped_stats$percentage <-
    paste0(mapped_stats$percentage * 100,"%")

  all_mapped_stats <-
    dplyr::bind_rows(
      mapped_stats, unmapped_stats
    ) |>
    dplyr::mutate(
      mapping_status = factor(
        mapping_status,
        levels = c("unmapped","mapped"))
    ) |>
    dplyr::mutate(mapping_rank = dplyr::if_else(
      mapping_status == "unmapped",0,1
    )) |>
    dplyr::mutate(
      color_label =
        paste0(treatment," / ", mapping_status)
    ) |>
    dplyr::mutate(
      color_label = factor(
        color_label,
        levels =
          c(paste0(treatment_name," / unmapped"),
            paste0(treatment_name," / mapped"),
            paste0(control_name," / unmapped"),
            paste0(control_name," / mapped"),
            "T0 / unmapped",
            "T0 / mapped"))
    )

  colors_unmapped <-
    ggsci::pal_jco(palette = c("default"), alpha = 0.6)(3)
  colors_mapped <-
    ggsci::pal_jco(palette = c("default"), alpha = 1)(3)

  p <- ggplot2::ggplot(
    data = all_mapped_stats,
    ggplot2::aes(x = timepoint,
                 y = reads,
                 fill = color_label)) +
    ggplot2::geom_bar(
      stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::ylab("sgRNA read count") +
    ggplot2::scale_fill_manual(
      values = c(colors_unmapped[2],colors_mapped[2],
                 colors_unmapped[1],colors_mapped[1],
                 colors_unmapped[3],colors_mapped[3]
      ), name = "treatment / read mapping status") +
    ggplot2::geom_text(
      ggplot2::aes(
        x = timepoint,
        label = percentage),
      colour = "white",
      size = 3,
      position = ggplot2::position_stack(vjust=0.6)) +
    ggplot2::theme(
      text = ggplot2::element_text(
        size = plot_fontsize - 9),
      strip.text.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(
        size = plot_fontsize, vjust = 0.5),
      axis.text.y = ggplot2::element_text(
        family = "Helvetica",
        size = plot_fontsize - 1),
      axis.title.x = ggplot2::element_blank(),
      #axis.title.x = ggplot2::element_text(
      #  family = "Helvetica",
      #  size = plot_fontsize,
      #  margin = ggplot2::margin(t = title_x_top_margin)),
      axis.title.y = ggplot2::element_text(
        family = "Helvetica",
        size = plot_fontsize,
        margin = ggplot2::margin(
          r = title_y_right_margin),
        vjust = 0.5),
      legend.text = ggplot2::element_text(
        family = "Helvetica",
        size = plot_fontsize)
    ) +
    ## use "free_x" to allow each facet
    ## to have its own x-scale (so no shared empty bars)
    ggplot2::facet_wrap(
      ~treatment, scales = "free_x", drop = TRUE)

  plots <- list()
  plots[['ggplot']] <- p
  plots[['plotly']] <- plotly::ggplotly(p)

  ## Adjust legend in plotly plot - make it horizontal
  ## - remove legend title
  plots[['plotly']]$x$layout$legend$title <- NULL
  plots[['plotly']] <- plotly::layout(
    plots[['plotly']],
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = -0.2
    ),
    margin = list(
      l = 120,
      r = 20,
      t = 50,
      b = 50)
  )

  return(plots)

}

#' Plot sgRNA Gini index bar plot
#'
#' This function generates a bar plot of the Gini index for sgRNA distributions
#'
#' @param sgRNA_summary A data frame containing sgRNA summary statistics including Gini index
#' @param ggsci_style The style for the color palette from the ggsci
#' package (default is "JCO")
#' @param plot_fontsize Font size for the plot text elements (default is 11
#'
#' @return A ggplot2 object representing the Gini index bar plot
#'
#' @export
#'
plot_sgRNA_gini_index <- function(
    sgRNA_summary,
    ggsci_style = "JCO",
    plot_fontsize = 11){

  p <- ggplot2::ggplot(
    data = sgRNA_summary,
    ggplot2::aes(x = timepoint,
                 y = gini_index,
                 fill = treatment)) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_dodge2()) +
    ggplot2::theme_classic() +
    ggsci::scale_fill_jco() +
    ggplot2::ylab("Gini index") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = plot_fontsize, vjust = 0.5),
      axis.text.y = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize - 1),
      #axis.title.x = ggplot2::element_text(
      #  family = "Helvetica", size = plot_fontsize),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize),
      legend.text = ggplot2::element_text(
        family = "Helvetica", size = plot_fontsize)
    )

  plots <- list()
  plots[['ggplot']] <- p
  plots[['plotly']] <- plotly::ggplotly(p)

  ## Adjust legend in plotly plot - make it horizontal
  ## - remove legend title
  plots[['plotly']]$x$layout$legend$title <- NULL
  plots[['plotly']] <- plotly::layout(
    plots[['plotly']],
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = -0.2
    ),
    margin = list(
      l = 80,
      r = 20,
      t = 50,
      b = 50)
  )

  return(plots)

}

#' Plot sgRNA count distribution histogram
#'
#' This function generates histograms of sgRNA count distributions
#' for different treatments and timepoints
#'
#' @param counts A data frame containing sgRNA counts with sgRNA IDs as row names
#' @param ggsci_style The style for the color palette from the ggsci
#' package (default is "JCO")
#' @param treatment_name The name of the treatment used in the experiment (default is "
#' FGF401")
#' @param control_name The name of the control used in the experiment (default is "
#' DMSO")
#' @param num_bins The number of bins to use in the histogram (default is 100)
#' @param title_x_top_margin Top margin for the x-axis title (default is 12)
#' @param title_y_right_margin Right margin for the y-axis title (default is
#' @param plot_fontsize Font size for the plot text elements (default is 11
#'
#' @return A list of ggplot2 objects representing the sgRNA count
#' distribution histograms
#' @export
#'
plot_sgRNA_histogram <- function(
    counts = NULL,
    ggsci_style = "JCO",
    treatment_name = "FGF401",
    control_name = "DMSO",
    num_bins = 100,
    title_x_top_margin = 12,
    title_y_right_margin = 12,
    plot_fontsize = 11){

  assertthat::assert_that(
    !is.null(counts),
    msg = "sgRNA_counts parameter must be provided")

  assertthat::assert_that(
    is.data.frame(counts),
    msg = "sgRNA_counts must be a data frame")

  #assertthat::assert_that(
  #  !is.null(rownames(counts)),
  #  msg = "sgRNA_counts must have row names representing sgRNA IDs")

  palette <- ggsci::pal_jco(
    palette = c("default"), alpha = 1)(3)

  sgRNA_df <- counts
  #sgRNA_df$sgRNA_ID <- rownames(sgRNA_df)
  sgRNA_long <- sgRNA_df |>
    tidyr::pivot_longer(
      cols = -c("sgRNA","Gene"),
      names_to = "label",
      values_to = "count") |>
    crisprFlow::clean_treatment_timepoint_label() |>
    dplyr::mutate(
      treatment2 = dplyr::case_when(
        treatment == treatment_name ~ "treatment",
        treatment == control_name ~ "control",
        TRUE ~ "T0"
      )
    ) |>
    dplyr::mutate(
      count = count + 1
    ) |>
    dplyr::filter(!is.na(count))

  readcounts <- list()
  for(tr in c('treatment','control','T0')){
    readcounts[[tr]] <- list()
    for(tp in unique(sgRNA_long$timepoint)){
      if(tp == "T0" & tr != "T0"){
        next
      }
      if(tp != "T0" & tr == "T0"){
        next
      }
      readcounts[[tr]][[tp]] <- list()
      readcounts[[tr]][[tp]][['plot']] <- NULL
      readcounts[[tr]][[tp]][['data']] <-
        sgRNA_long |>
        dplyr::filter(
          timepoint == tp &
            treatment2 == tr
        )

    }
  }

  for(tp in unique(sgRNA_long$timepoint)){
    for(tr in c("treatment","control",'T0')){
      df <- readcounts[[tr]][[tp]][['data']]
      if(NROW(df) == 0){
        next
      }

      ggplot_title <- "T0"
      if(tr != "T0"){
        ggplot_title <-
          paste0(
            unique(df$treatment),
            " - ", tp)
      }

      fill_color <- NULL
      if(tr == "T0"){
        fill_color <- palette[3]
      }else if(tr == "control"){
        fill_color <- palette[1]
      }else{
        fill_color <- palette[2]
      }

      median_readcount_val <-
        median(df$count, na.rm = TRUE)

      suppressWarnings({
        p <-
          ggplot2::ggplot(
            df,
            ggplot2::aes(x = count)) +
          ggplot2::geom_histogram(
            ggplot2::aes(
              text = ggplot2::after_stat(
                paste("frequency:", count,
                      "<br>bin range (read count): ",
                      round(xmin, 1), "–", round(xmax, 1),
                      "<br>bin midpoint (read count):", round(x,1)))),
            bins = num_bins,
            fill = fill_color,
            color = "black",
            alpha = 0.9) +
          ggplot2::scale_x_log10(
            breaks = c(1, 10, 100, 1000, 10000),
            labels = scales::label_number()) +
          ggplot2::geom_vline(
            xintercept = median_readcount_val,
            color = "#A73030FF",
            linetype = "dashed",
            linewidth = 0.7) +
          ggplot2::theme_minimal(base_size = 10) +
          ggplot2::ggtitle(ggplot_title) +
          ggplot2::labs(
            x = "sgRNA read count",
            y = "frequency") +
          ggplot2::annotation_logticks(sides = "b") +
          ggplot2::theme(
            text = ggplot2::element_text(size = plot_fontsize),
            axis.title.x = ggplot2::element_text(
              margin = ggplot2::margin(t = title_x_top_margin)),
            axis.title.y = ggplot2::element_text(
              margin = ggplot2::margin(r = title_y_right_margin)))
      })
      readcounts[[tr]][[tp]][['plot']] <- list()
      readcounts[[tr]][[tp]][['plot']][['ggplot']] <- p
      suppressWarnings(
        readcounts[[tr]][[tp]][['plot']][['plotly']] <-
          plotly::ggplotly(p, tooltip = "text")
      )
    }
  }

  return(readcounts)

}

#' Plot fold change scatterplot, drug vs. T0, control vs. T0
#'
#' This function generates a scatterplot comparing fold changes
#' between drug treatment and control conditions
#'
#' @param sgRNA_fc_df A data frame containing sgRNA fold change data
#' @param x_column The name of the column to use for the x-axis (default
#' is "T4_control_vs_T0_FC")
#' @param y_column The name of the column to use for the y-axis (default
#' is "T4_drug_vs_T0_FC")
#' @param x_lab The label for the x-axis (default is "T4 -
#' Control vs T0 - fold change (normalized sgRNA count)")
#' @param y_lab The label for the y-axis (default is "T4 -
#' Drug vs T0 - fold change (normalized sgRNA count)")
#' @param title The title of the scatterplot (default is "T4 - sgRNA fold change")
#' @param title_y_right_margin Right margin for the y-axis title (default is 12)
#' @param title_x_top_margin Top margin for the x-axis title (default is 12)
#' @param ggsci_style The style for the color palette from the ggsci
#' package (default is "JCO")
#' @param treatment_name The name of the treatment used in the experiment (default is "
#' FGF401")
#' @param per_sgRNA Logical indicating whether to plot at the sgRNA level
#' (TRUE) or gene level (FALSE) (default is TRUE)
#' @param highlight_genes Logical indicating whether to highlight specific genes
#' on the scatterplot (default is TRUE)
#' @param plot_fontsize Font size for the plot text elements (default is 11
#' @return A ggplot2 object representing the fold change scatterplot
#'
#' @export
#'
plot_fc_scatterplot <- function(
  sgRNA_fc_df = NULL,
  x_column = "T4_control_vs_T0_FC",
  y_column = "T4_drug_vs_T0_FC",
  x_lab = "Control vs T0",
  y_lab = "Drug vs T0",
  title = "T4 - sgRNA abundance LFC",
  title_y_right_margin = 12,
  title_x_top_margin = 12,
  ggsci_style = "JCO",
  treatment_name = "FGF401",
  per_sgRNA = TRUE, # TRUE = sgRNA-level, FALSE = gene-level
  highlight_genes = TRUE,
  plot_fontsize = 13){

  assertthat::assert_that(
    !is.null(sgRNA_fc_df),
    msg = "sgRNA_fc_df parameter must be provided"
  )
  assertthat::assert_that(
    is.data.frame(sgRNA_fc_df),
    msg = "sgRNA_fc_df must be a data frame"
  )
  assertable::assert_colnames(
    sgRNA_fc_df,
    c(x_column, y_column,"Gene"),
    only_colnames = F, quiet = TRUE
  )

  plots <- list()

  plots[['ggplot']] <- sgRNA_fc_df |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = !!rlang::sym(x_column),
        y = !!rlang::sym(y_column),
        z = Gene))
  if(per_sgRNA == TRUE){
    plots[['ggplot']] <- sgRNA_fc_df |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = !!rlang::sym(x_column),
          y = !!rlang::sym(y_column)))
  }
  plots[['ggplot']] <-
    plots[['ggplot']] +
    ggplot2::geom_point(
      colour = "#A73030FF") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::ggtitle(title) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab) +
    ggplot2::theme(
      legend.position = "none",
      text = ggplot2::element_text(
        family = "Helvetica",
        size = plot_fontsize),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = title_x_top_margin)),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = title_y_right_margin)))

  plots[['plotly']] <-
    plotly::ggplotly(plots[['ggplot']])

  ## - remove legend title
  plots[['plotly']]$x$layout$legend$title <- NULL
  return(plots)

}



