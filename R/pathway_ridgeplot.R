#' Ridge Plot for GSEA Results
#'
#' Creates a ridge plot (joy plot) to visualize the distribution of gene/KO
#' abundances or fold changes for enriched pathways from GSEA analysis.
#' This helps interpret whether pathways are predominantly up- or down-regulated.
#'
#' @param gsea_results A data frame containing GSEA results from \code{\link{pathway_gsea}}.
#'   Must contain pathway_id column and either NES or direction column.
#' @param abundance A data frame or matrix containing the original abundance data
#'   (genes/KOs as rows, samples as columns) used in the GSEA analysis. Data
#'   frames may also provide a leading non-numeric feature ID column (for
#'   example \code{#NAME}, \code{feature}, or \code{pathway}); it is converted
#'   to row names before sample alignment.
#' @param metadata A data frame containing sample metadata with group information.
#' @param group Character string specifying the column name in metadata for grouping.
#' @param comparison Optional character vector of length 2 specifying
#'   \code{c(group1, group2)} for fold-change calculation. The ridge plot
#'   displays \code{log2(group2 / group1)}. If NULL, the data must contain
#'   exactly two group levels and their factor-level order is used. For
#'   multi-group metadata, specify \code{comparison} explicitly so the plot
#'   matches the GSEA contrast being interpreted.
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
#' Fold changes are calculated as \code{log2(group2 / group1)}, where
#' \code{group1} and \code{group2} come from \code{comparison} or, for a
#' two-group dataset, from the factor-level order of \code{group}.
#' GSEA \code{pathway_id} values must be non-empty and unique. Missing or
#' empty \code{pathway_name} values are displayed as their corresponding
#' \code{pathway_id}.
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
                               comparison = NULL,
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

  abundance <- normalize_abundance_feature_ids(
    abundance,
    context = "pathway_ridgeplot() abundance"
  )
  validate_metadata(metadata)
  validate_abundance(abundance, min_samples = 2)
  validate_count_parameter(n_pathways, "n_pathways")
  validate_choice(sort_by, c("NES", "pvalue", "p.adjust"), "sort_by")
  valid_pathway_types <- c("KEGG", "GO", "MetaCyc")
  if (!is.character(pathway_type) || length(pathway_type) != 1 ||
      is.na(pathway_type) || !pathway_type %in% valid_pathway_types) {
    stop("pathway_type must be one of: ",
         paste(valid_pathway_types, collapse = ", "),
         call. = FALSE)
  }
  show_direction <- normalize_logical_flag(show_direction, "show_direction")
  validate_probability_threshold(alpha, "alpha", allow_zero = TRUE)
  if (!is.numeric(scale_height) || length(scale_height) != 1 ||
      is.na(scale_height) || !is.finite(scale_height) || scale_height <= 0) {
    stop("'scale_height' must be a single positive finite numeric value.",
         call. = FALSE)
  }

  # Convert abundance to matrix if needed
  if (is.data.frame(abundance)) {
    abundance <- as.matrix(abundance)
  }
  if (is.null(rownames(abundance))) {
    stop("'abundance' must have row names matching the gene/KO identifiers in pathway_reference.",
         call. = FALSE)
  }
  validate_numeric_matrix(abundance)
  validate_finite_numeric_values(as.vector(abundance), "abundance", "abundance",
                                 allow_na = FALSE)

  # Align samples before calculating group means. The fold-change displayed in
  # the ridge plot is group2/group1 across abundance columns, so metadata rows
  # must be in the same sample order as the abundance matrix columns.
  aligned <- align_samples(abundance, metadata, verbose = FALSE)
  abundance <- as.matrix(aligned$abundance)
  metadata <- aligned$metadata
  validate_group(metadata, group, min_groups = 2)
  if (any(is.na(metadata[[group]]))) {
    stop("Group column '", group, "' contains NA values after sample alignment.",
         call. = FALSE)
  }

  # Validate required columns (pathway_gsea() output format)
  require_column(gsea_results, "pathway_id", "gsea_results")
  gsea_results$pathway_id <- validate_nonempty_character_column(
    gsea_results$pathway_id,
    "pathway_id",
    "gsea_results"
  )
  if (anyDuplicated(gsea_results$pathway_id)) {
    duplicated_pathways <- unique(gsea_results$pathway_id[duplicated(gsea_results$pathway_id)])
    stop("gsea_results contains duplicate pathway_id values: ",
         paste(utils::head(duplicated_pathways, 5), collapse = ", "),
         ". Subset to one result row per pathway before plotting.",
         call. = FALSE)
  }
  require_column(gsea_results, sort_by, "gsea_results")
  if ("pvalue" %in% colnames(gsea_results)) {
    validate_probability_values(gsea_results$pvalue, "pvalue", "gsea_results")
  }
  if ("p.adjust" %in% colnames(gsea_results)) {
    validate_probability_values(gsea_results$p.adjust, "p.adjust", "gsea_results")
  }

  # pathway_name is optional, fallback row-wise to pathway_id
  if ("pathway_name" %in% colnames(gsea_results)) {
    pathway_labels <- as.character(gsea_results$pathway_name)
    missing_labels <- is.na(pathway_labels) | !nzchar(trimws(pathway_labels))
    pathway_labels[missing_labels] <- gsea_results$pathway_id[missing_labels]
    gsea_results$.ridge_pathway_label <- pathway_labels
  } else {
    gsea_results$.ridge_pathway_label <- gsea_results$pathway_id
  }
  pathway_name_col <- ".ridge_pathway_label"

  # direction is optional, can be derived from NES
  has_nes <- "NES" %in% colnames(gsea_results)
  if (has_nes) {
    validate_finite_numeric_values(gsea_results$NES, "NES", "gsea_results")
  }
  direction_col <- if ("direction" %in% colnames(gsea_results)) "direction" else NULL

  # Filter and sort. P-values/FDRs sort ascending; NES sorts by absolute effect
  # size so the top list can include both enriched and depleted pathways.
  df <- gsea_results
  df <- df[!is.na(df[[sort_by]]), ]
  if (sort_by == "NES") {
    df <- df[order(abs(df[[sort_by]]), decreasing = TRUE), ]
  } else {
    df <- df[order(df[[sort_by]]), ]
  }
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

 # Calculate fold changes between groups. The comparison must be tied to
 # group labels, not to the first groups encountered in sample order; otherwise
 # simply reordering abundance columns can flip the x-axis interpretation.
 group_factor <- droplevels(factor(metadata[[group]]))
 observed_levels <- levels(group_factor)
 if (is.null(comparison)) {
   if (length(observed_levels) != 2) {
     stop(
       "pathway_ridgeplot() requires 'comparison = c(group1, group2)' ",
       "when metadata contains ", length(observed_levels), " group levels. ",
       "The ridge plot displays log2(group2 / group1), so the comparison ",
       "must match the GSEA contrast being interpreted.",
       call. = FALSE
     )
   }
   group_levels <- observed_levels
 } else {
   comparison <- validate_nonempty_character_column(
     comparison,
     "comparison",
     "pathway_ridgeplot()"
   )
   if (length(comparison) != 2) {
     stop("'comparison' must be NULL or a character vector of length 2: c(group1, group2).",
          call. = FALSE)
   }
   if (anyDuplicated(comparison)) {
     stop("'comparison' must contain two distinct group labels.",
          call. = FALSE)
   }
   missing_groups <- setdiff(comparison, observed_levels)
   if (length(missing_groups) > 0) {
     stop(
       "'comparison' group(s) not found in metadata after sample alignment: ",
       paste(missing_groups, collapse = ", "),
       ". Available groups: ", paste(observed_levels, collapse = ", "),
       call. = FALSE
     )
   }
   group_levels <- comparison
 }
 group_vec <- as.character(group_factor)

 # Calculate mean abundance per group
 samples_g1 <- which(as.character(group_vec) == group_levels[1])
 samples_g2 <- which(as.character(group_vec) == group_levels[2])

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
     if (any(grepl(";", as.character(pathway_reference[[col]]), fixed = TRUE))) {
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
   pid <- df[["pathway_id"]][i]
   pname <- df[[pathway_name_col]][i]

   # Truncate long names
   if (nchar(pname) > 50) {
     pname <- paste0(substr(pname, 1, 47), "...")
   }

   # Get direction if available
   direction <- if (!is.null(direction_col)) {
     d <- as.character(df[[direction_col]][i])
     if (is.na(d)) "Unknown" else d
   } else if (has_nes && !is.na(df$NES[i])) {
     if (df$NES[i] > 0) "Up" else "Down"
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
         gene_id = names(fc_values),
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

 if (show_direction) {
   if (!is.character(colors) || is.null(names(colors))) {
     stop("'colors' must be a named character vector when show_direction = TRUE.",
          call. = FALSE)
   }
   invalid_colors <- vapply(colors, function(x) {
     !tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(e) FALSE)
   }, logical(1))
   if (any(invalid_colors)) {
     stop("Invalid color value(s) in 'colors': ",
          paste(names(colors)[invalid_colors], collapse = ", "),
          call. = FALSE)
   }
   missing_directions <- setdiff(unique(as.character(ridge_data$direction)),
                                 names(colors))
   if (length(missing_directions) > 0) {
     if (all(missing_directions == "Unknown")) {
       colors <- c(colors, Unknown = "#737373")
     } else {
       stop("'colors' is missing value(s) for direction(s): ",
            paste(missing_directions, collapse = ", "),
            call. = FALSE)
     }
   }
 }

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
