#' Ridge Plot for GSEA Results
#'
#' Creates a ridge plot (joy plot) to visualize the distribution of gene/KO
#' abundances or fold changes for enriched pathways from GSEA analysis.
#' This helps interpret whether pathways are predominantly up- or down-regulated.
#'
#' @param gsea_results A data frame containing GSEA results from \code{\link{pathway_gsea}}.
#'   Must contain pathway_id column and either NES or direction column.
#' @param abundance A data frame or matrix containing the original abundance data
#'   (genes/KOs as rows, samples as columns) used in the GSEA analysis.
#' @param metadata A data frame containing sample metadata with group information.
#' @param group Character string specifying the column name in metadata for grouping.
#' @param pathway_reference A data frame containing pathway-to-gene mappings.
#'   Must have columns: pathway_id (or go_id for GO) and a column containing
#'   gene/KO members (semicolon-separated). If NULL, attempts to use built-in
#'   KEGG or GO reference data.
#' @param pathway_type Character string specifying the pathway type: "KEGG", "GO", or "MetaCyc".
#'   Default is "KEGG".
#' @param n_pathways Integer specifying the number of top pathways to display.
#'   Default is 10.
#' @param sort_by Character string specifying how to sort pathways:
#'   "NES" (Normalized Enrichment Score), "pvalue", or "p.adjust".
#'   Default is "p.adjust".
#' @param show_direction Logical. If TRUE, colors ridges by enrichment direction.
#'   Default is TRUE.
#' @param colors Named character vector with colors for "Up" and "Down" directions.
#'   Default is blue for down-regulated and red for up-regulated.
#' @param title Character string for plot title.
#' @param x_lab Character string for x-axis label.
#' @param scale_height Numeric value controlling the overlap of ridges.
#'   Default is 0.9. Higher values create more overlap.
#' @param alpha Numeric value for ridge transparency (0-1). Default is 0.7.
#'
#' @return A ggplot2 object that can be further customized or saved.
#'
#' @details
#' The ridge plot displays the distribution of gene abundances (or fold changes)
#' for genes within each enriched pathway. This visualization helps to:
#' \itemize{
#'   \item Understand the overall direction of change for each pathway
#'   \item Identify pathways with consistent vs. heterogeneous gene expression
#'   \item Compare the magnitude of changes across pathways
#' }
#'
#' The plot requires the \code{ggridges} package to be installed.
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
#' # Run GSEA (using camera method - recommended)
#' gsea_results <- pathway_gsea(
#'   abundance = ko_abundance %>% column_to_rownames("#NAME"),
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "camera"
#' )
#'
#' # Create ridge plot
#' ridge_plot <- pathway_ridgeplot(
#'   gsea_results = gsea_results,
#'   abundance = ko_abundance %>% column_to_rownames("#NAME"),
#'   metadata = metadata,
#'   group = "Environment",
#'   n_pathways = 10
#' )
#' print(ridge_plot)
#' }
#'
#' @seealso \code{\link{pathway_gsea}}, \code{\link{visualize_gsea}},
#'   \code{\link{pathway_volcano}}
#'
#' @importFrom ggplot2 ggplot aes labs theme_minimal theme element_text
#'   scale_fill_manual element_blank
#' @importFrom dplyr filter mutate arrange slice_head left_join
#' @importFrom tidyr pivot_longer
#' @importFrom stats setNames
#'
#' @export
pathway_ridgeplot <- function(gsea_results,
                               abundance,
                               metadata,
                               group,
                               pathway_reference = NULL,
                               pathway_type = "KEGG",
                               n_pathways = 10,
                               sort_by = "p.adjust",
                               show_direction = TRUE,
                               colors = c("Down" = "#3182bd", "Up" = "#de2d26"),
                               title = "Ridge Plot: Gene Distribution in Enriched Pathways",
                               x_lab = "log2 Fold Change",
                               scale_height = 0.9,
                               alpha = 0.7) {

  # Check for ggridges package
 if (!requireNamespace("ggridges", quietly = TRUE)) {
   stop("Package 'ggridges' is required for ridge plots. ",
        "Please install it with: install.packages('ggridges')")
 }

 # Input validation
 if (!is.data.frame(gsea_results)) {
   stop("'gsea_results' must be a data frame.")
 }

 if (!is.data.frame(abundance) && !is.matrix(abundance)) {
   stop("'abundance' must be a data frame or matrix.")
 }

 if (!group %in% colnames(metadata)) {
   stop(paste0("Column '", group, "' not found in metadata."))
 }

 # Convert abundance to matrix if needed
 if (is.data.frame(abundance)) {
   abundance <- as.matrix(abundance)
 }

 # Determine pathway ID column
 pathway_id_col <- if ("pathway_id" %in% colnames(gsea_results)) {
   "pathway_id"
 } else if ("go_id" %in% colnames(gsea_results)) {
   "go_id"
 } else if ("ID" %in% colnames(gsea_results)) {
   "ID"
 } else {
   stop("Cannot find pathway ID column in gsea_results. Expected: pathway_id, go_id, or ID")
 }

 # Determine pathway name column
 pathway_name_col <- if ("pathway_name" %in% colnames(gsea_results)) {
   "pathway_name"
 } else if ("go_name" %in% colnames(gsea_results)) {
   "go_name"
 } else if ("Description" %in% colnames(gsea_results)) {
   "Description"
 } else {
   pathway_id_col  # Use ID as name if no name column
 }

 # Determine direction/NES column
 has_nes <- "NES" %in% colnames(gsea_results)
 has_direction <- "direction" %in% colnames(gsea_results) ||
                   "Direction" %in% colnames(gsea_results)

 direction_col <- if ("direction" %in% colnames(gsea_results)) {
   "direction"
 } else if ("Direction" %in% colnames(gsea_results)) {
   "Direction"
 } else {
   NULL
 }

 # Sort and select top pathways
 sort_col <- if (sort_by %in% colnames(gsea_results)) {
   sort_by
 } else if ("pvalue" %in% colnames(gsea_results)) {
   "pvalue"
 } else {
   stop("Cannot find sorting column in gsea_results.")
 }

 # Filter and sort
 df <- gsea_results
 df <- df[!is.na(df[[sort_col]]), ]
 df <- df[order(df[[sort_col]]), ]
 df <- head(df, n_pathways)

 if (nrow(df) == 0) {
   stop("No pathways to display after filtering.")
 }

 # Get pathway-gene mappings
 if (is.null(pathway_reference)) {
   # Load reference data using unified loader
   if (pathway_type == "KEGG") {
     kegg_ref <- load_reference_data("ko_to_kegg")
     # Aggregate KO IDs by pathway_id
     pathway_reference <- stats::aggregate(
       ko_id ~ pathway_id + pathway_name,
       data = kegg_ref,
       FUN = function(x) paste(unique(x), collapse = ";")
     )
     colnames(pathway_reference)[colnames(pathway_reference) == "ko_id"] <- "ko_members"
   } else if (pathway_type == "GO") {
     pathway_reference <- load_reference_data("ko_to_go")
   } else {
     stop("Please provide pathway_reference parameter for ", pathway_type, " pathways.")
   }
 }

 # Calculate fold changes between groups
 group_vec <- metadata[[group]]
 group_levels <- unique(group_vec)

 if (length(group_levels) != 2) {
   warning("Ridge plot works best with 2 groups. Using first two groups.")
   group_levels <- group_levels[1:2]
 }

 # Calculate mean abundance per group
 samples_g1 <- which(group_vec == group_levels[1])
 samples_g2 <- which(group_vec == group_levels[2])

 mean_g1 <- rowMeans(abundance[, samples_g1, drop = FALSE], na.rm = TRUE)
 mean_g2 <- rowMeans(abundance[, samples_g2, drop = FALSE], na.rm = TRUE)

 # Calculate log2 fold change using unified function
 pseudocount <- calculate_pseudocount(as.vector(abundance))
 log2fc <- calculate_log2_fold_change(mean_g1, mean_g2, pseudocount = pseudocount)
 names(log2fc) <- rownames(abundance)

 # Determine gene member column in pathway_reference
 # Support both wide format (semicolon-separated) and long format (one gene per row)
 gene_col_candidates <- c("ko_members", "KO", "ko_id", "genes", "gene", "ec_numbers")
 gene_col <- NULL
 for (col in gene_col_candidates) {
   if (col %in% colnames(pathway_reference)) {
     gene_col <- col
     break
   }
 }
 if (is.null(gene_col)) {
   # Fallback: try to find any column with semicolon-separated values
   for (col in colnames(pathway_reference)) {
     if (any(grepl(";", pathway_reference[[col]], fixed = TRUE))) {
       gene_col <- col
       break
     }
   }
 }
 if (is.null(gene_col)) {
   stop("Cannot find gene member column in pathway_reference. ",
        "Expected columns: ", paste(gene_col_candidates, collapse = ", "))
 }

 # Build data for ridge plot
 # Using list + do.call(rbind, ...) for O(n) vs O(n²) performance
 ridge_data_list <- vector("list", nrow(df))

 # Determine ref_id_col once outside the loop
 ref_id_col <- if ("pathway_id" %in% colnames(pathway_reference)) {
   "pathway_id"
 } else if ("go_id" %in% colnames(pathway_reference)) {
   "go_id"
 } else {
   colnames(pathway_reference)[1]
 }

 for (i in seq_len(nrow(df))) {
   pid <- df[[pathway_id_col]][i]
   pname <- df[[pathway_name_col]][i]

   # Truncate long names
   if (nchar(pname) > 50) {
     pname <- paste0(substr(pname, 1, 47), "...")
   }

   # Get direction if available
   direction <- if (!is.null(direction_col)) {
     as.character(df[[direction_col]][i])
   } else if (has_nes) {
     ifelse(df$NES[i] > 0, "Up", "Down")
   } else {
     "Unknown"
   }

   ref_rows <- pathway_reference[pathway_reference[[ref_id_col]] == pid, ]

   if (nrow(ref_rows) > 0) {
     # Handle both long format (multiple rows) and wide format (semicolon-separated)
     if (nrow(ref_rows) > 1) {
       # Long format: each row is one gene
       genes <- unique(as.character(ref_rows[[gene_col]]))
     } else {
       # Wide format: genes are semicolon-separated in one row
       genes_str <- ref_rows[[gene_col]][1]
       genes <- unlist(strsplit(as.character(genes_str), ";"))
     }
     genes <- trimws(genes)

     # Get fold changes for these genes
     fc_values <- log2fc[names(log2fc) %in% genes]

     if (length(fc_values) > 0) {
       ridge_data_list[[i]] <- data.frame(
         pathway = pname,
         pathway_id = pid,
         direction = direction,
         log2fc = fc_values,
         stringsAsFactors = FALSE
       )
     }
   }
 }

 # Combine all data frames at once (O(n) instead of O(n²))
 non_null_list <- ridge_data_list[!sapply(ridge_data_list, is.null)]
 if (length(non_null_list) == 0) {
   stop("No gene data found for the selected pathways. ",
        "Check that pathway_reference matches your abundance data.")
 }
 ridge_data <- do.call(rbind, non_null_list)

 if (is.null(ridge_data) || nrow(ridge_data) == 0) {
   stop("No gene data found for the selected pathways. ",
        "Check that pathway_reference matches your abundance data.")
 }

 # Order pathways by their appearance in the sorted results
 pathway_order <- unique(df[[pathway_name_col]])
 pathway_order <- sapply(pathway_order, function(x) {
   if (nchar(x) > 50) paste0(substr(x, 1, 47), "...") else x
 })
 ridge_data$pathway <- factor(ridge_data$pathway, levels = rev(pathway_order))

 # Create ridge plot
 p <- ggplot2::ggplot(ridge_data, ggplot2::aes(x = .data$log2fc, y = .data$pathway))

 if (show_direction && "direction" %in% colnames(ridge_data)) {
   p <- p + ggridges::geom_density_ridges(
     ggplot2::aes(fill = .data$direction),
     alpha = alpha,
     scale = scale_height
   ) +
     ggplot2::scale_fill_manual(values = colors, name = "Direction")
 } else {
   p <- p + ggridges::geom_density_ridges(
     fill = "#756bb1",
     alpha = alpha,
     scale = scale_height
   )
 }

 p <- p +
   ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", alpha = 0.7) +
   ggplot2::labs(
     x = x_lab,
     y = NULL,
     title = title
   ) +
   ggplot2::theme_minimal() +
   ggplot2::theme(
     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
     axis.text.y = ggplot2::element_text(size = 9),
     legend.position = "top",
     panel.grid.minor = ggplot2::element_blank()
   )

 return(p)
}
