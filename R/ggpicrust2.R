#' This function integrates pathway name/description annotations, ten of the most advanced differential abundance (DA) methods, and visualization of DA results.
#'
#' @param file A character string representing the file path of the input file containing KO abundance data in picrust2 export format. The input file should have KO identifiers in the first column and sample identifiers in the first row. The remaining cells should contain the abundance values for each KO-sample pair.
#' @param data An optional data.frame containing KO abundance data in the same format as the input file. If provided, the function will use this data instead of reading from the file. By default, this parameter is set to NULL.
#' @param metadata A tibble, consisting of sample information
#' @param group A character, name of the group
#' @param pathway A character, consisting of "EC", "KO", "MetaCyc"
#' @param daa_method a character specifying the method for differential abundance analysis, default is "ALDEx2", choices are:
#' - "ALDEx2": ANOVA-Like Differential Expression tool for high throughput sequencing data
#' - "DESeq2": Differential expression analysis based on the negative binomial distribution using DESeq2
#' - "edgeR": Exact test for differences between two groups of negative-binomially distributed counts using edgeR
#' - "limma voom": Limma-voom framework for the analysis of RNA-seq data
#' - "metagenomeSeq": Fit logistic regression models to test for differential abundance between groups using metagenomeSeq
#' - "LinDA": Linear models for differential abundance analysis of microbiome compositional data
#' - "Maaslin2": Multivariate Association with Linear Models (MaAsLin2) for differential abundance analysis
#' @param ko_to_kegg A character to control the conversion of KO abundance to KEGG abundance
#' @param filter_for_prokaryotes Logical. If TRUE (default), filters out KEGG pathways
#'   that are specific to eukaryotes (e.g., human diseases, organismal systems) when
#'   ko_to_kegg = TRUE. Set to FALSE to include all KEGG pathways.
#' @param p.adjust a character specifying the method for p-value adjustment, default is "BH", choices are:
#'- "BH": Benjamini-Hochberg correction
#'- "holm": Holm's correction
#'- "bonferroni": Bonferroni correction
#'- "hochberg": Hochberg's correction
#'- "fdr": False discovery rate correction
#'- "none": No p-value adjustment.
#' @param order A character to control the order of the main plot rows
#' @param p_values_bar A character to control if the main plot has the p_values bar
#' @param x_lab A character to control the x-axis label name, you can choose from "feature","pathway_name" and "description"
#' @param select A vector consisting of pathway names to be selected
#' @param reference A character, a reference group level for several DA methods
#' @param colors A vector consisting of colors number
#' @param p_values_threshold A numeric value specifying the threshold for statistical
#'   significance of differential abundance. Pathways with adjusted p-values below this
#'   threshold will be displayed in the plot. Default is 0.05.
#'
#' @return A list containing:
#' \itemize{
#'   \item Numbered elements (1, 2, ...): Sub-lists for each DA method, each containing:
#'     \itemize{
#'       \item \code{plot}: A ggplot2 error bar plot visualizing the differential abundance results
#'       \item \code{results}: A data frame of differential abundance results for that method
#'     }
#'   \item \code{abundance}: The processed abundance data (KEGG pathway or original) for downstream analysis
#'   \item \code{metadata}: The metadata data frame
#'   \item \code{group}: The group variable name used in the analysis
#'   \item \code{daa_results_df}: The complete annotated DAA results data frame
#'   \item \code{ko_to_kegg}: Logical indicating whether KO to KEGG conversion was performed
#' }
#' These additional fields allow seamless integration with \code{\link{pathway_pca}} and
#' \code{\link{pathway_heatmap}} for further visualization without re-preparing data.
#' @export
#' @examples
#' \dontrun{
#' # Load necessary data: abundance data and metadata
#' abundance_file <- "path/to/your/abundance_file.tsv"
#' metadata <- read.csv("path/to/your/metadata.csv")
#'
#' # Run ggpicrust2 with input file path
#' results_file_input <- ggpicrust2(file = abundance_file,
#'                                  metadata = metadata,
#'                                  group = "your_group_column",
#'                                  pathway = "KO",
#'                                  daa_method = "LinDA",
#'                                  ko_to_kegg = "TRUE",
#'                                  order = "pathway_class",
#'                                  p_values_bar = TRUE,
#'                                  x_lab = "pathway_name")
#'
#' # Run ggpicrust2 with imported data.frame
#' abundance_data <- read_delim(abundance_file, delim="\t", col_names=TRUE, trim_ws=TRUE)
#'
#' # Run ggpicrust2 with input data
#' results_data_input <- ggpicrust2(data = abundance_data,
#'                                  metadata = metadata,
#'                                  group = "your_group_column",
#'                                  pathway = "KO",
#'                                  daa_method = "LinDA",
#'                                  ko_to_kegg = "TRUE",
#'                                  order = "pathway_class",
#'                                  p_values_bar = TRUE,
#'                                  x_lab = "pathway_name")
#'
#' # Access the plot and results dataframe for the first DA method
#' example_plot <- results_file_input[[1]]$plot
#' example_results <- results_file_input[[1]]$results
#'
#' # Use the example data in ggpicrust2 package
#' data(ko_abundance)
#' data(metadata)
#' results_file_input <- ggpicrust2(data = ko_abundance,
#'                                  metadata = metadata,
#'                                  group = "Environment",
#'                                  pathway = "KO",
#'                                  daa_method = "LinDA",
#'                                  ko_to_kegg = TRUE,
#'                                  order = "pathway_class",
#'                                  p_values_bar = TRUE,
#'                                  x_lab = "pathway_name")
#' # Analyze the EC or MetaCyc pathway
#' data(metacyc_abundance)
#' results_file_input <- ggpicrust2(data = metacyc_abundance,
#'                                  metadata = metadata,
#'                                  group = "Environment",
#'                                  pathway = "MetaCyc",
#'                                  daa_method = "LinDA",
#'                                  ko_to_kegg = FALSE,
#'                                  order = "group",
#'                                  p_values_bar = TRUE,
#'                                  x_lab = "description")
#'
#' # Use the returned data for PCA analysis (no need to re-prepare data)
#' pca_plot <- pathway_pca(
#'   abundance = results_file_input$abundance,
#'   metadata = results_file_input$metadata,
#'   group = results_file_input$group
#' )
#'
#' # Use the returned data for heatmap (filter significant pathways first)
#' sig_features <- results_file_input$daa_results_df %>%
#'   dplyr::filter(p_adjust < 0.05) %>%
#'   dplyr::pull(feature)
#' if (length(sig_features) > 0) {
#'   heatmap_plot <- pathway_heatmap(
#'     abundance = results_file_input$abundance[sig_features, , drop = FALSE],
#'     metadata = results_file_input$metadata,
#'     group = results_file_input$group
#'   )
#' }
#'}
ggpicrust2 <- function(file = NULL,
                       data = NULL,
                       metadata,
                       group,
                       pathway,
                       daa_method = "ALDEx2",
                       ko_to_kegg = FALSE,
                       filter_for_prokaryotes = TRUE,
                       p.adjust = "BH",
                       order = "group",
                       p_values_bar = TRUE,
                       x_lab = NULL,
                       select = NULL,
                       reference = NULL,
                       colors = NULL,
                       p_values_threshold = 0.05) {
  # Input validation
  if (is.null(file) && is.null(data)) {
    stop("Error: Please provide either a 'file' path (character string) or a 'data' data frame.")
  }

  if (!is.null(file) && !is.null(data)) {
    warning("Both 'file' and 'data' provided. Using 'data' and ignoring 'file'.")
    file <- NULL  # Set file to NULL to use data
  }

  # Validate file parameter if provided
  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1) {
      stop("Error: 'file' parameter must be a character string representing a file path. If you have already loaded your data into R, please use the 'data' parameter instead.")
    }
    if (!file.exists(file)) {
      stop("Error: File does not exist: ", file)
    }
  }

  # Validate data parameter if provided
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("Error: 'data' parameter must be a data frame. If you want to provide a file path, please use the 'file' parameter instead.")
    }
  }

  # Validate daa_method early

  if (daa_method == "Lefse") {
    stop("The 'Lefse' method is not suitable for ggpicrust2() as it does not output p-values.")
  }


  # Step 1: Load abundance data
  if (ko_to_kegg) {
    message("Converting KO to KEGG...\n")
    abundance <- if (!is.null(file)) {
      ko2kegg_abundance(file, filter_for_prokaryotes = filter_for_prokaryotes)
    } else {
      ko2kegg_abundance(data = data, filter_for_prokaryotes = filter_for_prokaryotes)
    }
  } else {
    message("Reading input data...\n")
    abundance <- if (!is.null(file)) {
      read_abundance_file(file)
    } else {
      data
    }
    abundance <- as.data.frame(abundance)
    abundance <- clean_ko_abundance(abundance)
    rownames(abundance) <- abundance[, 1]
    abundance <- abundance[, -1]
  }

  # Align abundance and metadata once for consistent downstream behavior.
  # This prevents Group/order mismatches in pathway_errorbar() when metadata
  # row order differs from abundance column order.
  aligned <- align_samples(abundance, metadata, verbose = FALSE)
  abundance <- aligned$abundance
  metadata <- aligned$metadata

  # Step 2: Differential abundance analysis
  message("Performing pathway differential abundance analysis...\n")
  daa_results_df <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = group,
    daa_method = daa_method,
    select = select,
    p.adjust = p.adjust,
    reference = reference
  )

  # Check for significant biomarkers
  num_significant <- sum(daa_results_df$p_adjust < p_values_threshold, na.rm = TRUE)
  if (num_significant == 0) {
    warning(sprintf("No statistically significant biomarkers found (p_adjust < %g). ", p_values_threshold),
            "Analysis will continue for visualization purposes.", call. = FALSE)
  }

  # Step 3: Pathway annotation
  message("Annotating pathways...\n")
  daa_results_df <- pathway_annotation(
    pathway = pathway,
    ko_to_kegg = ko_to_kegg,
    daa_results_df = daa_results_df,
    p_adjust_threshold = p_values_threshold
  )

  # Step 4: Create plots for each method
  message("Creating pathway error bar plots...\n")
  plot_result_list <- list()
  Group_vec <- metadata[[group]]
  names(Group_vec) <- colnames(abundance)

  for (i in seq_along(unique(daa_results_df$method))) {
    method_name <- unique(daa_results_df$method)[i]
    daa_sub_method_results_df <- daa_results_df[daa_results_df$method == method_name, ]

    combination_bar_plot <- pathway_errorbar(
      abundance = abundance,
      daa_results_df = daa_sub_method_results_df,
      Group = Group_vec,
      ko_to_kegg = ko_to_kegg,
      p_value_bar = p_values_bar,
      order = order,
      colors = colors,
      select = select,
      x_lab = x_lab,
      p_values_threshold = p_values_threshold
    )

    if (is.null(combination_bar_plot)) {
      message(sprintf("Plot %d skipped (no data for method: %s)\n", i, method_name))
    } else {
      message(sprintf("Plot %d created.\n", i))
    }

    plot_result_list[[i]] <- list(plot = combination_bar_plot, results = daa_sub_method_results_df)
  }

  message("ggpicrust2 analysis completed.\n")

  # Step 5: Add data for downstream analysis

  plot_result_list$abundance <- abundance
  plot_result_list$metadata <- metadata
  plot_result_list$group <- group
  plot_result_list$daa_results_df <- daa_results_df
  plot_result_list$ko_to_kegg <- ko_to_kegg

  return(plot_result_list)
}
