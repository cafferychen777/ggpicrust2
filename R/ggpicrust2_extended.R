#' Integrated analysis with ggpicrust2 including GSEA
#'
#' This function extends the ggpicrust2 functionality to include Gene Set Enrichment Analysis (GSEA).
#'
#' @param ... Parameters passed to ggpicrust2()
#' @param run_gsea Logical value indicating whether to perform GSEA analysis
#' @param gsea_params List of parameters to pass to pathway_gsea()
#'
#' @return A list containing ggpicrust2 and GSEA results
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#'
#' # Run integrated analysis with camera method (recommended)
#' integrated_results <- ggpicrust2_extended(
#'   data = ko_abundance,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway = "KO",
#'   daa_method = "LinDA",
#'   ko_to_kegg = TRUE,
#'   run_gsea = TRUE,
#'   gsea_params = list(
#'     method = "camera"
#'   )
#' )
#'
#' # Run with covariate adjustment
#' integrated_results_adj <- ggpicrust2_extended(
#'   data = ko_abundance,
#'   metadata = metadata,
#'   group = "Disease",
#'   pathway = "KO",
#'   daa_method = "LinDA",
#'   ko_to_kegg = TRUE,
#'   run_gsea = TRUE,
#'   gsea_params = list(
#'     method = "camera",
#'     covariates = c("age", "sex")
#'   )
#' )
#'
#' # Access DAA results
#' daa_results <- integrated_results$daa_results
#'
#' # Access GSEA results
#' gsea_results <- integrated_results$gsea_results
#'
#' # Access plots
#' daa_plot <- integrated_results$daa_plot
#' gsea_plot <- integrated_results$gsea_plot
#' }
ggpicrust2_extended <- function(..., 
                               run_gsea = FALSE, 
                               gsea_params = list()) {
  
  # Run standard ggpicrust2 analysis
  ggpicrust2_results <- ggpicrust2(...)
  
  # Extract arguments for GSEA
  args <- list(...)
  
  # Initialize results list
  results <- list(
    daa_results = ggpicrust2_results
  )
  
  # Run GSEA if requested
  if (run_gsea) {
    # Check if required packages are available
    if (!requireNamespace("limma", quietly = TRUE)) {
      warning("Package 'limma' is required for GSEA analysis. Skipping GSEA. ",
              "Install with: BiocManager::install('limma')")
    } else {
      # Extract necessary data for GSEA
      if (!is.null(args$data)) {
        abundance <- args$data
      } else if (!is.null(args$file)) {
        abundance <- read_abundance_file(args$file)
      } else {
        stop("No abundance data provided for GSEA analysis")
      }
      
      # Prepare abundance data
      if (is.data.frame(abundance)) {
        # Convert tibble to data.frame and set row names
        abundance <- as.data.frame(abundance)
        # Assuming first column contains feature IDs
        rownames(abundance) <- abundance[, 1]
        abundance <- abundance[, -1, drop = FALSE]
      }
      
      # Extract metadata and group
      metadata <- args$metadata
      group <- args$group
      
      # Determine pathway type
      pathway_type <- "KEGG"
      if (!is.null(args$pathway)) {
        if (args$pathway %in% c("KO", "KEGG")) {
          pathway_type <- "KEGG"
        } else if (args$pathway == "MetaCyc") {
          pathway_type <- "MetaCyc"
        } else if (args$pathway == "EC") {
          pathway_type <- "EC"
        }
      }
      
      # Set default GSEA parameters
      default_params <- list(
        abundance = abundance,
        metadata = metadata,
        group = group,
        pathway_type = pathway_type,
        method = "camera",
        min_size = 5,
        max_size = 500,
        p.adjust = "BH"
      )
      
      # Override defaults with user-provided parameters
      gsea_params <- utils::modifyList(default_params, gsea_params)
      
      # Run GSEA analysis
      message("Performing Gene Set Enrichment Analysis...")
      gsea_results <- do.call(pathway_gsea, gsea_params)
      
      # Annotate GSEA results
      message("Annotating GSEA results...")
      annotated_gsea_results <- gsea_pathway_annotation(
        gsea_results = gsea_results,
        pathway_type = pathway_type
      )
      
      # Create visualization
      message("Creating GSEA visualization...")
      gsea_plot <- visualize_gsea(
        gsea_results = annotated_gsea_results,
        plot_type = "barplot",
        n_pathways = 20,
        sort_by = "p.adjust"
      )
      
      # Compare GSEA and DAA results
      message("Comparing GSEA and DAA results...")
      daa_results <- ggpicrust2_results[[1]]$results
      comparison <- compare_gsea_daa(
        gsea_results = annotated_gsea_results,
        daa_results = daa_results,
        plot_type = "venn"
      )
      
      # Add GSEA results to output
      results$gsea_results <- annotated_gsea_results
      results$gsea_plot <- gsea_plot
      results$comparison <- comparison
    }
  }
  
  return(results)
}
