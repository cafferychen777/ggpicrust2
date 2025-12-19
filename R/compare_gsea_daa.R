#' Compare GSEA and DAA results
#'
#' This function compares the results from Gene Set Enrichment Analysis (GSEA) and
#' Differential Abundance Analysis (DAA) to identify similarities and differences.
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param daa_results A data frame containing DAA results from the pathway_daa function
#' @param plot_type A character string specifying the visualization type: "venn", "upset", "scatter", or "heatmap"
#' @param p_threshold A numeric value specifying the significance threshold
#'
#' @return A ggplot2 object or a list containing the plot and comparison results
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#'
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#'
#' # Run GSEA analysis (using camera method - recommended)
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "camera"
#' )
#'
#' # Run DAA analysis
#' daa_results <- pathway_daa(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment"
#' )
#'
#' # Compare results
#' comparison <- compare_gsea_daa(
#'   gsea_results = gsea_results,
#'   daa_results = daa_results,
#'   plot_type = "venn"
#' )
#' }
compare_gsea_daa <- function(gsea_results,
                            daa_results,
                            plot_type = "venn",
                            p_threshold = 0.05) {

  # Input validation
  if (!is.data.frame(gsea_results)) {
    stop("'gsea_results' must be a data frame")
  }

  if (!is.data.frame(daa_results)) {
    stop("'daa_results' must be a data frame")
  }

  if (length(plot_type) != 1 || !plot_type %in% c("venn", "upset", "scatter", "heatmap")) {
    stop("plot_type must be one of 'venn', 'upset', 'scatter', or 'heatmap'")
  }

  # Check if required columns exist
  if (!all(c("pathway_id", "p.adjust") %in% colnames(gsea_results))) {
    stop("GSEA results missing required columns: pathway_id, p.adjust")
  }

  if (!all(c("feature", "p_adjust") %in% colnames(daa_results))) {
    stop("DAA results missing required columns: feature, p_adjust")
  }

  # Extract significant pathways from each analysis
  # Convert to character to ensure consistent types for set operations and ggVennDiagram
  sig_gsea <- as.character(gsea_results$pathway_id[gsea_results$p.adjust < p_threshold])
  sig_daa <- as.character(daa_results$feature[daa_results$p_adjust < p_threshold])

  # Find overlapping and unique pathways
  overlap <- intersect(sig_gsea, sig_daa)
  gsea_only <- setdiff(sig_gsea, sig_daa)
  daa_only <- setdiff(sig_daa, sig_gsea)

  # Create comparison results
  comparison_results <- list(
    overlap = overlap,
    gsea_only = gsea_only,
    daa_only = daa_only,
    n_overlap = length(overlap),
    n_gsea_only = length(gsea_only),
    n_daa_only = length(daa_only),
    n_gsea_total = length(sig_gsea),
    n_daa_total = length(sig_daa)
  )

  # Create visualization based on plot_type
  if (plot_type == "venn") {
    # Check if required package is available
    if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
      warning("Package 'ggVennDiagram' is required for Venn diagrams. Using a basic plot instead.")

      # Create a basic representation
      counts <- c(comparison_results$n_gsea_only,
                 comparison_results$n_overlap,
                 comparison_results$n_daa_only)

      labels <- c("GSEA only", "Overlap", "DAA only")

      p <- ggplot2::ggplot(data.frame(counts = counts, labels = labels),
                         ggplot2::aes(x = labels, y = counts, fill = labels)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          title = "Comparison of Significant Pathways",
          x = "",
          y = "Number of Pathways",
          fill = ""
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5)
        )
    } else {
      # Create a proper Venn diagram
      venn_list <- list(
        GSEA = sig_gsea,
        DAA = sig_daa
      )

      p <- ggVennDiagram::ggVennDiagram(venn_list) +
        ggplot2::labs(title = "Overlap of Significant Pathways") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }

  } else if (plot_type == "upset") {
    # Check if required package is available
    if (!requireNamespace("UpSetR", quietly = TRUE)) {
      warning("Package 'UpSetR' is required for UpSet plots. Using a basic plot instead.")

      # Create a basic representation
      counts <- c(comparison_results$n_gsea_only,
                 comparison_results$n_overlap,
                 comparison_results$n_daa_only)

      labels <- c("GSEA only", "Overlap", "DAA only")

      p <- ggplot2::ggplot(data.frame(counts = counts, labels = labels),
                         ggplot2::aes(x = labels, y = counts, fill = labels)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          title = "Comparison of Significant Pathways",
          x = "",
          y = "Number of Pathways",
          fill = ""
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5)
        )
    } else {
      # Create a proper UpSet plot
      # First, create a binary matrix for the UpSetR input
      all_pathways <- unique(c(sig_gsea, sig_daa))
      upset_matrix <- matrix(0, nrow = length(all_pathways), ncol = 2)
      rownames(upset_matrix) <- all_pathways
      colnames(upset_matrix) <- c("GSEA", "DAA")

      upset_matrix[sig_gsea, "GSEA"] <- 1
      upset_matrix[sig_daa, "DAA"] <- 1

      # Convert to data frame for UpSetR
      upset_df <- as.data.frame(upset_matrix)

      # Create UpSet plot
      # We can't easily convert UpSetR plots to ggplot objects
      # So we'll create a basic ggplot representation instead
      p <- UpSetR::upset(upset_df, nsets = 2, order.by = "freq")

      # Create a simple ggplot representation
      intersect_size <- sum(upset_matrix[, "GSEA"] & upset_matrix[, "DAA"])
      gsea_only <- sum(upset_matrix[, "GSEA"]) - intersect_size
      daa_only <- sum(upset_matrix[, "DAA"]) - intersect_size

      plot_data <- data.frame(
        set = c("GSEA only", "Intersection", "DAA only"),
        count = c(gsea_only, intersect_size, daa_only)
      )

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$set, y = .data$count, fill = .data$set)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          title = "Set Intersection",
          x = "",
          y = "Count",
          fill = "Set"
        ) +
        ggplot2::theme_minimal()
    }

  } else if (plot_type == "scatter") {
    # Create a scatter plot comparing p-values or effect sizes

    # Merge the results
    merged_results <- merge(
      gsea_results[, c("pathway_id", "NES", "p.adjust")],
      daa_results[, c("feature", "log_2_fold_change", "p_adjust")],
      by.x = "pathway_id",
      by.y = "feature",
      all = FALSE
    )

    # If no overlapping pathways, return a message
    if (nrow(merged_results) == 0) {
      warning("No overlapping pathways found for scatter plot")
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No overlapping pathways found") +
        ggplot2::theme_void()
    } else {
      # Create scatter plot
      p <- ggplot2::ggplot(merged_results,
                         ggplot2::aes(x = .data$NES, y = .data$log_2_fold_change,
                                    color = -log10(.data$p.adjust))) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "blue", high = "red") +
        ggplot2::labs(
          title = "Comparison of GSEA and DAA Results",
          x = "Normalized Enrichment Score (GSEA)",
          y = "Log2 Fold Change (DAA)",
          color = "-log10(p.adjust)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray")
    }

  } else if (plot_type == "heatmap") {
    # Create a heatmap
    # This would require more complex implementation
    # For now, we'll return a placeholder
    warning("Heatmap plot not yet implemented")
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "Heatmap plot not yet implemented") +
      ggplot2::theme_void()
  }

  # Return both the plot and the comparison results
  return(list(
    plot = p,
    results = comparison_results
  ))
}
