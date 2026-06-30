#' Compare GSEA and DAA results
#'
#' This function compares the results from Gene Set Enrichment Analysis (GSEA) and
#' Differential Abundance Analysis (DAA) to identify similarities and differences.
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param daa_results A data frame containing DAA results from the pathway_daa function
#' @param plot_type A single character string specifying the visualization
#'   type: "venn", "upset", or "scatter"
#' @param p_threshold A numeric value specifying the significance threshold.
#'   Must be in the range (0, 1].
#'
#' @details
#' Venn and UpSet plots compare significant pathway sets with unique pathway
#' IDs. The scatter plot compares effect sizes and therefore requires one row
#' per pathway in each input. If a result table contains multiple methods,
#' contrasts, or group pairs for the same pathway, filter it to a single
#' effect-size context before using \code{plot_type = "scatter"}.
#' For scatter plots, GSEA and DAA directions must be explicit. Preranked
#' GSEA positive NES values represent \code{gsea_results$group1} versus
#' \code{gsea_results$group2}, while DAA \code{log2_fold_change} values
#' represent \code{daa_results$group2 / daa_results$group1}. The DAA effect
#' size is therefore aligned to the GSEA-positive direction before plotting.
#' For Venn and UpSet plots, if both inputs include \code{group1} and
#' \code{group2}, each input must represent one comparable group pair. This
#' prevents pathways found in different biological contrasts from being counted
#' as method agreement.
#'
#' @return A list with two elements: \code{plot} (a ggplot2 object, or an
#'   UpSetR object when \code{plot_type = "upset"} and UpSetR is installed)
#'   and \code{results} (a named list with the overlap, GSEA-only, and
#'   DAA-only pathway sets plus their counts). For \code{plot_type = "scatter"},
#'   \code{results$scatter_data} contains the merged effect-size table,
#'   including \code{daa_log2_fold_change_aligned}.
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

  validate_choice(plot_type, c("venn", "upset", "scatter"), "plot_type")
  validate_probability_threshold(p_threshold, "p_threshold")

  # Check if required columns exist. Baseline (for set-membership plots)
  # is just the ID + adjusted p. The scatter variant additionally plots
  # effect sizes (NES vs log2_fold_change), so it needs those columns;
  # checking here rather than only inside the scatter branch means the
  # user gets an up-front, single error naming exactly what's missing
  # instead of a cryptic "undefined columns selected" from merge().
  if (!all(c("pathway_id", "p.adjust") %in% colnames(gsea_results))) {
    stop("GSEA results missing required columns: pathway_id, p.adjust")
  }

  if (!all(c("feature", "p_adjust") %in% colnames(daa_results))) {
    stop("DAA results missing required columns: feature, p_adjust")
  }
  validate_probability_values(gsea_results[["p.adjust"]], "p.adjust", "gsea_results")
  validate_probability_values(daa_results[["p_adjust"]], "p_adjust", "daa_results")
  gsea_pathways <- validate_nonempty_character_column(gsea_results[["pathway_id"]],
                                                      "pathway_id",
                                                      "gsea_results")
  daa_features <- validate_nonempty_character_column(daa_results[["feature"]],
                                                    "feature",
                                                    "daa_results")

  if (plot_type == "scatter") {
    missing_gsea <- setdiff("NES", colnames(gsea_results))
    missing_daa <- setdiff("log2_fold_change", colnames(daa_results))
    if (length(missing_gsea) > 0 || length(missing_daa) > 0) {
      stop(sprintf(
        "plot_type = 'scatter' requires effect-size columns: %s%s%s. Use plot_type = 'venn' or 'upset' if effect sizes are unavailable.",
        if (length(missing_gsea)) paste0("gsea_results needs '", missing_gsea, "'") else "",
        if (length(missing_gsea) && length(missing_daa)) "; " else "",
        if (length(missing_daa)) paste0("daa_results needs '", missing_daa, "'") else ""
      ), call. = FALSE)
    }
    validate_finite_numeric_values(gsea_results[["NES"]], "NES", "gsea_results")
    validate_finite_numeric_values(daa_results[["log2_fold_change"]],
                                   "log2_fold_change", "daa_results")

    missing_gsea_direction <- setdiff(c("group1", "group2"), colnames(gsea_results))
    missing_daa_direction <- setdiff(c("group1", "group2"), colnames(daa_results))
    if (length(missing_gsea_direction) > 0 || length(missing_daa_direction) > 0) {
      stop(
        "plot_type = 'scatter' requires explicit direction columns for effect-size comparison. ",
        if (length(missing_gsea_direction) > 0) {
          paste0("gsea_results missing: ",
                 paste(missing_gsea_direction, collapse = ", "), ". ")
        } else "",
        if (length(missing_daa_direction) > 0) {
          paste0("daa_results missing: ",
                 paste(missing_daa_direction, collapse = ", "), ". ")
        } else "",
        "Re-run pathway_gsea() with an explicit 'comparison' and use DAA results with group1/group2 columns.",
        call. = FALSE
      )
    }
    gsea_group1 <- validate_nonempty_character_column(gsea_results[["group1"]],
                                                      "group1",
                                                      "gsea_results")
    gsea_group2 <- validate_nonempty_character_column(gsea_results[["group2"]],
                                                      "group2",
                                                      "gsea_results")
    daa_group1 <- validate_nonempty_character_column(daa_results[["group1"]],
                                                    "group1",
                                                    "daa_results")
    daa_group2 <- validate_nonempty_character_column(daa_results[["group2"]],
                                                    "group2",
                                                    "daa_results")

    gsea_pairs <- unique(data.frame(group1 = gsea_group1,
                                    group2 = gsea_group2,
                                    stringsAsFactors = FALSE))
    daa_pairs <- unique(data.frame(group1 = daa_group1,
                                  group2 = daa_group2,
                                  stringsAsFactors = FALSE))
    if (nrow(gsea_pairs) > 1 || nrow(daa_pairs) > 1) {
      stop(
        "plot_type = 'scatter' requires one GSEA comparison and one DAA group pair. ",
        if (nrow(gsea_pairs) > 1) {
          paste0("gsea_results contains ", nrow(gsea_pairs), " comparisons. ")
        } else "",
        if (nrow(daa_pairs) > 1) {
          paste0("daa_results contains ", nrow(daa_pairs), " group pairs. ")
        } else "",
        "Filter each input to a single effect-size direction before plotting.",
        call. = FALSE
      )
    }

    duplicated_gsea <- unique(gsea_pathways[duplicated(gsea_pathways)])
    duplicated_daa <- unique(daa_features[duplicated(daa_features)])
    if (length(duplicated_gsea) > 0 || length(duplicated_daa) > 0) {
      stop(
        "plot_type = 'scatter' requires one effect-size row per pathway in each input. ",
        if (length(duplicated_gsea) > 0) {
          paste0("Duplicate pathway_id values in gsea_results: ",
                 paste(utils::head(duplicated_gsea, 5), collapse = ", "), ". ")
        } else "",
        if (length(duplicated_daa) > 0) {
          paste0("Duplicate feature values in daa_results: ",
                 paste(utils::head(duplicated_daa, 5), collapse = ", "), ". ")
        } else "",
        "Filter to one GSEA method/contrast and one DAA method/group pair before plotting effect sizes.",
        call. = FALSE
      )
    }
  } else {
    validate_compare_gsea_daa_set_direction(
      gsea_results,
      daa_results,
      plot_type
    )
  }

  # Extract significant pathways from each analysis
  # Convert to character to ensure consistent types for set operations and ggVennDiagram
  sig_gsea <- unique(gsea_pathways[which(gsea_results$p.adjust < p_threshold)])
  sig_daa <- unique(daa_features[which(daa_results$p_adjust < p_threshold)])

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
      # Create UpSet plot
      all_pathways <- unique(c(sig_gsea, sig_daa))
      upset_matrix <- matrix(0, nrow = length(all_pathways), ncol = 2)
      rownames(upset_matrix) <- all_pathways
      colnames(upset_matrix) <- c("GSEA", "DAA")

      upset_matrix[sig_gsea, "GSEA"] <- 1
      upset_matrix[sig_daa, "DAA"] <- 1

      upset_df <- as.data.frame(upset_matrix)
      p <- UpSetR::upset(upset_df, nsets = 2, order.by = "freq")
    }

  } else if (plot_type == "scatter") {
    # Create a scatter plot comparing p-values or effect sizes

    gsea_for_merge <- data.frame(
      pathway_id = gsea_pathways,
      NES = gsea_results[["NES"]],
      p.adjust = gsea_results[["p.adjust"]],
      gsea_group1 = gsea_group1,
      gsea_group2 = gsea_group2,
      stringsAsFactors = FALSE
    )
    daa_for_merge <- data.frame(
      feature = daa_features,
      log2_fold_change = daa_results[["log2_fold_change"]],
      p_adjust = daa_results[["p_adjust"]],
      daa_group1 = daa_group1,
      daa_group2 = daa_group2,
      stringsAsFactors = FALSE
    )

    # Merge the results
    merged_results <- merge(
      gsea_for_merge,
      daa_for_merge,
      by.x = "pathway_id",
      by.y = "feature",
      all = FALSE
    )

    # If no overlapping pathways, return a message
    if (nrow(merged_results) == 0) {
      warning("No overlapping pathways found for scatter plot")
      comparison_results$scatter_data <- merged_results
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No overlapping pathways found") +
        ggplot2::theme_void()
    } else {
      daa_aligned_to_gsea <- merged_results$gsea_group1 == merged_results$daa_group2 &
        merged_results$gsea_group2 == merged_results$daa_group1
      daa_opposite_to_gsea <- merged_results$gsea_group1 == merged_results$daa_group1 &
        merged_results$gsea_group2 == merged_results$daa_group2
      incompatible_direction <- !(daa_aligned_to_gsea | daa_opposite_to_gsea)
      if (any(incompatible_direction)) {
        incompatible_rows <- merged_results$pathway_id[incompatible_direction]
        stop(
          "GSEA and DAA group directions are not comparable for overlapping pathway(s): ",
          paste(utils::head(incompatible_rows, 5), collapse = ", "),
          ". GSEA uses group1 versus group2, while DAA log2_fold_change uses group2/group1.",
          call. = FALSE
        )
      }
      merged_results$daa_log2_fold_change_aligned <- ifelse(
        daa_aligned_to_gsea,
        merged_results$log2_fold_change,
        -merged_results$log2_fold_change
      )
      merged_results$gsea_neg_log10_p_adjust <-
        -log10(pmax(merged_results$p.adjust, .Machine$double.xmin))
      comparison_results$scatter_data <- merged_results

      # Create scatter plot
      p <- ggplot2::ggplot(merged_results,
                         ggplot2::aes(x = .data$NES, y = .data$daa_log2_fold_change_aligned,
                                    color = .data$gsea_neg_log10_p_adjust)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "blue", high = "red") +
        ggplot2::labs(
          title = "Comparison of GSEA and DAA Results",
          x = "Normalized Enrichment Score (GSEA)",
          y = "Log2 Fold Change (DAA, aligned to GSEA direction)",
          color = "-log10(p.adjust)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray")
    }
  }

  # Return both the plot and the comparison results
  return(list(
    plot = p,
    results = comparison_results
  ))
}

validate_compare_gsea_daa_set_direction <- function(gsea_results,
                                                    daa_results,
                                                    plot_type) {
  direction_cols <- c("group1", "group2")
  gsea_has_direction <- all(direction_cols %in% colnames(gsea_results))
  daa_has_direction <- all(direction_cols %in% colnames(daa_results))
  if (!gsea_has_direction || !daa_has_direction) {
    return(invisible(NULL))
  }

  gsea_group1 <- validate_nonempty_character_column(gsea_results[["group1"]],
                                                    "group1",
                                                    "gsea_results")
  gsea_group2 <- validate_nonempty_character_column(gsea_results[["group2"]],
                                                    "group2",
                                                    "gsea_results")
  daa_group1 <- validate_nonempty_character_column(daa_results[["group1"]],
                                                  "group1",
                                                  "daa_results")
  daa_group2 <- validate_nonempty_character_column(daa_results[["group2"]],
                                                  "group2",
                                                  "daa_results")

  gsea_pairs <- unique(data.frame(group1 = gsea_group1,
                                  group2 = gsea_group2,
                                  stringsAsFactors = FALSE))
  daa_pairs <- unique(data.frame(group1 = daa_group1,
                                group2 = daa_group2,
                                stringsAsFactors = FALSE))
  if (nrow(gsea_pairs) > 1 || nrow(daa_pairs) > 1) {
    stop(
      "plot_type = '", plot_type,
      "' compares significant pathway sets and requires one comparable group pair per input. ",
      if (nrow(gsea_pairs) > 1) {
        paste0("gsea_results contains ", nrow(gsea_pairs), " group pairs. ")
      } else "",
      if (nrow(daa_pairs) > 1) {
        paste0("daa_results contains ", nrow(daa_pairs), " group pairs. ")
      } else "",
      "Filter each input to a single biological comparison before comparing sets.",
      call. = FALSE
    )
  }

  same_pair <- identical(gsea_pairs$group1, daa_pairs$group1) &&
    identical(gsea_pairs$group2, daa_pairs$group2)
  reversed_pair <- identical(gsea_pairs$group1, daa_pairs$group2) &&
    identical(gsea_pairs$group2, daa_pairs$group1)
  if (!same_pair && !reversed_pair) {
    stop(
      "GSEA and DAA group pairs are not comparable for plot_type = '",
      plot_type, "': GSEA uses ",
      gsea_pairs$group1, " vs ", gsea_pairs$group2,
      ", while DAA uses ", daa_pairs$group1, " vs ", daa_pairs$group2,
      ". Filter to matching biological comparisons before comparing sets.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
