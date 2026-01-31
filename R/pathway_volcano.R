#' Volcano Plot for Pathway Differential Abundance Analysis
#'
#' Creates a volcano plot to visualize the results of differential abundance analysis,
#' showing both statistical significance (-log10 p-value) and effect size (log2 fold change).
#'
#' @param daa_results A data frame containing differential abundance analysis results,
#'   typically from \code{\link{pathway_daa}}. Must contain columns for fold change,
#'   p-values, and optionally pathway names.
#' @param fc_col Character string specifying the column name for log2 fold change values.
#'   Default is "log2_fold_change". Legacy name "log2FoldChange" is also accepted.
#' @param p_col Character string specifying the column name for adjusted p-values.
#'   Default is "p_adjust".
#' @param label_col Character string specifying the column name for pathway labels.
#'   Default is "pathway_name". If NULL, no labels will be shown.
#' @param fc_threshold Numeric. Absolute fold change threshold for significance.
#'   Default is 1 (2-fold change). Pathways with |log2FC| > fc_threshold are considered
#'   biologically significant.
#' @param p_threshold Numeric. P-value threshold for statistical significance.
#'   Default is 0.05.
#' @param label_top_n Integer. Number of top significant pathways to label.
#'   Default is 10. Set to 0 to disable labels.
#' @param point_size Numeric. Size of points in the plot. Default is 2.
#' @param point_alpha Numeric. Transparency of points (0-1). Default is 0.6.
#' @param colors Named character vector with colors for "Down", "Not Significant", and "Up".
#'   Default uses blue for down-regulated, grey for non-significant, and red for up-regulated.
#' @param show_threshold_lines Logical. Whether to show dashed lines for fold change
#'   and p-value thresholds. Default is TRUE.
#' @param title Character string for plot title. Default is
#'   "Volcano Plot: Pathway Differential Abundance".
#' @param x_lab Character string for x-axis label. Default is "log2 Fold Change".
#' @param y_lab Character string for y-axis label. Default is "-log10(Adjusted P-value)".
#'
#' @return A ggplot2 object that can be further customized or saved.
#'
#' @details
#' The volcano plot is a scatter plot that shows statistical significance (y-axis)

#' versus fold change (x-axis). Points are colored by significance:
#' \itemize{
#'   \item \strong{Up (red)}: Pathways with p < p_threshold AND log2FC > fc_threshold
#'   \item \strong{Down (blue)}: Pathways with p < p_threshold AND log2FC < -fc_threshold
#'   \item \strong{Not Significant (grey)}: All other pathways
#' }
#'
#' The function automatically labels the top N most significant pathways using
#' \code{ggrepel::geom_text_repel()} if the ggrepel package is installed.
#'
#' @examples
#' \dontrun{
#' library(ggpicrust2)
#' library(tibble)
#'
#' # Load example data
#' data("ko_abundance")
#' data("metadata")
#'
#' # Convert KO to KEGG abundance
#' kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
#'
#' # Run differential abundance analysis
#' daa_results <- pathway_daa(
#'   abundance = kegg_abundance,
#'   metadata = metadata,
#'   group = "Environment",
#'   daa_method = "LinDA"
#' )
#'
#' # Annotate results
#' daa_annotated <- pathway_annotation(
#'   pathway = "KO",
#'   ko_to_kegg = TRUE,
#'   daa_results_df = daa_results
#' )
#'
#' # Create volcano plot
#' volcano_plot <- pathway_volcano(daa_annotated)
#' print(volcano_plot)
#'
#' # Customize the plot
#' volcano_plot <- pathway_volcano(
#'   daa_annotated,
#'   fc_threshold = 0.5,
#'   p_threshold = 0.01,
#'   label_top_n = 15,
#'   colors = c("Down" = "darkblue", "Not Significant" = "lightgrey", "Up" = "darkred")
#' )
#' }
#'
#' @seealso \code{\link{pathway_daa}}, \code{\link{pathway_annotation}},
#'   \code{\link{pathway_errorbar}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_text
#'   scale_color_manual labs theme_minimal theme element_text
#' @importFrom dplyr filter mutate case_when slice_max
#' @importFrom stats setNames
#'
#' @export
pathway_volcano <- function(daa_results,
                            fc_col = "log2_fold_change",
                            p_col = "p_adjust",
                            label_col = "pathway_name",
                            fc_threshold = 1,
                            p_threshold = 0.05,
                            label_top_n = 10,
                            point_size = 2,
                            point_alpha = 0.6,
                            colors = c("Down" = "#3182bd",
                                       "Not Significant" = "grey60",
                                       "Up" = "#de2d26"),
                            show_threshold_lines = TRUE,
                            title = "Volcano Plot: Pathway Differential Abundance",
                            x_lab = "log2 Fold Change",
                            y_lab = "-log10(Adjusted P-value)") {

 # Input validation
 if (!is.data.frame(daa_results)) {
   stop("'daa_results' must be a data frame.")
 }

 # Backward compatibility: accept legacy column name
 if (fc_col == "log2_fold_change" && !fc_col %in% colnames(daa_results) &&
     "log2FoldChange" %in% colnames(daa_results)) {
   fc_col <- "log2FoldChange"
 }

 require_column(daa_results, fc_col, "daa_results")
 require_column(daa_results, p_col, "daa_results")

 if (!is.null(label_col) && !label_col %in% colnames(daa_results)) {
   warning(paste0("Column '", label_col, "' not found. Labels will not be shown."))
   label_col <- NULL
 }

 # Prepare data
 df <- daa_results

 # Filter out rows with missing values in key columns
 df <- df[!is.na(df[[fc_col]]) & !is.na(df[[p_col]]), ]

 if (nrow(df) == 0) {
   stop("No valid data rows after removing NA values.")
 }

 # Add significance categories and -log10(p)
 df$neg_log10_p <- -log10(df[[p_col]])
 # Handle infinite values (when p-value is 0)
 if (any(is.infinite(df$neg_log10_p))) {
   max_finite <- max(df$neg_log10_p[is.finite(df$neg_log10_p)], na.rm = TRUE)
   df$neg_log10_p[is.infinite(df$neg_log10_p)] <- max_finite * 1.1
 }
 df$significance <- ifelse(
   df[[p_col]] < p_threshold & df[[fc_col]] > fc_threshold, "Up",
   ifelse(df[[p_col]] < p_threshold & df[[fc_col]] < -fc_threshold, "Down",
          "Not Significant")
 )

 # Create base plot
 p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[fc_col]], y = .data$neg_log10_p,
                                        color = .data$significance)) +
   ggplot2::geom_point(alpha = point_alpha, size = point_size) +
   ggplot2::scale_color_manual(values = colors) +
   ggplot2::labs(
     x = x_lab,
     y = y_lab,
     color = "Significance",
     title = title
   ) +
   ggplot2::theme_minimal() +
   ggplot2::theme(
     legend.position = "top",
     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
     axis.text = ggplot2::element_text(color = "black"),
     panel.grid.minor = ggplot2::element_blank()
   )

 # Add threshold lines if requested
 if (show_threshold_lines) {
   p <- p +
     ggplot2::geom_hline(yintercept = -log10(p_threshold),
                         linetype = "dashed", alpha = 0.5, color = "grey40") +
     ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                         linetype = "dashed", alpha = 0.5, color = "grey40")
 }

 # Add labels for top significant pathways
 if (!is.null(label_col) && label_top_n > 0) {
   # Get top pathways by significance (excluding non-significant and NA labels)
   top_pathways <- df[df$significance != "Not Significant", ]
   # Filter out rows with NA or empty labels
   top_pathways <- top_pathways[!is.na(top_pathways[[label_col]]) &
                                 top_pathways[[label_col]] != "", ]
   if (nrow(top_pathways) > 0) {
     top_pathways <- top_pathways[order(-top_pathways$neg_log10_p), ]
     top_pathways <- head(top_pathways, label_top_n)

     # Use ggrepel if available for better label placement
     if (requireNamespace("ggrepel", quietly = TRUE)) {
       p <- p + ggrepel::geom_text_repel(
         data = top_pathways,
         ggplot2::aes(label = .data[[label_col]]),
         size = 2.5,
         max.overlaps = 20,
         segment.alpha = 0.5,
         segment.size = 0.3,
         show.legend = FALSE
       )
     } else {
       # Fallback without ggrepel
       p <- p + ggplot2::geom_text(
         data = top_pathways,
         ggplot2::aes(label = .data[[label_col]]),
         size = 2.5,
         vjust = -0.8,
         hjust = 0.5,
         check_overlap = TRUE,
         show.legend = FALSE
       )
     }
   }
 }

 # Add count summary in subtitle
 n_up <- sum(df$significance == "Up")
 n_down <- sum(df$significance == "Down")
 n_ns <- sum(df$significance == "Not Significant")

 p <- p + ggplot2::labs(
   subtitle = paste0("Up: ", n_up, " | Down: ", n_down, " | Not Sig: ", n_ns)
 )

 return(p)
}
