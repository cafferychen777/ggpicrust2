#' Generate Abundance Statistics Table for Pathway Analysis
#'
#' This function generates a table containing mean relative abundance, 
#' standard deviation, and log2 fold change statistics for pathways, 
#' similar to the data used in pathway_errorbar plots but returned as a 
#' data frame instead of a plot.
#'
#' @param abundance A data frame or matrix containing predicted functional 
#'        pathway abundance, with pathways/features as rows and samples as 
#'        columns. The column names should match the sample names in metadata.
#' @param daa_results_df A data frame containing differential abundance 
#'        analysis results from pathway_daa function. Must contain columns: 
#'        feature, group1, group2, p_adjust.
#' @param Group A vector containing group assignments for each sample in the
#'        same order as the columns in abundance matrix. Alternatively, if
#'        metadata is provided, this should match the order of samples in metadata.
#' @param ko_to_kegg Logical value indicating whether to use KO to KEGG 
#'        conversion. Default is FALSE.
#' @param p_values_threshold Numeric value for p-value threshold to filter 
#'        significant features. Default is 0.05.
#' @param select Character vector of specific features to include. If NULL, 
#'        all significant features are included.
#' @param max_features Maximum number of features to include in the table.
#'        Default is 30.
#' @param metadata Optional data frame containing sample metadata. If provided,
#'        the Group vector will be reordered to match the abundance column order.
#' @param sample_col Character string specifying the column name in metadata
#'        that contains sample identifiers. Default is "sample_name".
#'
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item \code{feature}: Feature/pathway identifier
#'   \item \code{group1}: Reference group name
#'   \item \code{group2}: Comparison group name  
#'   \item \code{mean_rel_abundance_group1}: Mean relative abundance for group1
#'   \item \code{sd_rel_abundance_group1}: Standard deviation of relative 
#'         abundance for group1
#'   \item \code{mean_rel_abundance_group2}: Mean relative abundance for group2
#'   \item \code{sd_rel_abundance_group2}: Standard deviation of relative 
#'         abundance for group2
#'   \item \code{log2_fold_change}: Log2 fold change (group2/group1)
#'   \item \code{p_adjust}: Adjusted p-value from differential analysis
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("ko_abundance")
#' data("metadata")
#' 
#' # Convert KO abundance to KEGG pathways
#' kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
#' 
#' # Perform differential abundance analysis
#' daa_results_df <- pathway_daa(
#'   abundance = kegg_abundance,
#'   metadata = metadata,
#'   group = "Environment",
#'   daa_method = "ALDEx2"
#' )
#' 
#' # Filter for specific method
#' daa_sub_method_results_df <- daa_results_df[
#'   daa_results_df$method == "ALDEx2_Welch's t test", 
#' ]
#' 
#' # Annotate results
#' daa_annotated_sub_method_results_df <- pathway_annotation(
#'   pathway = "KO",
#'   daa_results_df = daa_sub_method_results_df,
#'   ko_to_kegg = TRUE
#' )
#' 
#' # Generate abundance statistics table
#' abundance_stats_table <- pathway_errorbar_table(
#'   abundance = kegg_abundance,
#'   daa_results_df = daa_annotated_sub_method_results_df,
#'   Group = metadata$Environment,
#'   ko_to_kegg = TRUE,
#'   p_values_threshold = 0.05
#' )
#' 
#' # View the results
#' head(abundance_stats_table)
#' }
#'
#' @export
pathway_errorbar_table <- function(abundance,
                                  daa_results_df,
                                  Group,
                                  ko_to_kegg = FALSE,
                                  p_values_threshold = 0.05,
                                  select = NULL,
                                  max_features = 30,
                                  metadata = NULL,
                                  sample_col = "sample_name") {
  
  # Input validation
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("'abundance' must be a matrix or data frame")
  }

  if (!is.data.frame(daa_results_df)) {
    stop("'daa_results_df' must be a data frame")
  }

  # Handle Group vector ordering
  if (!is.null(metadata)) {
    # If metadata is provided, reorder Group to match abundance column order
    if (!sample_col %in% colnames(metadata)) {
      stop("Sample column '", sample_col, "' not found in metadata")
    }

    # Create a mapping from sample to group
    sample_to_group <- setNames(Group, metadata[[sample_col]])

    # Reorder Group to match abundance column order
    abundance_samples <- colnames(abundance)
    Group <- sample_to_group[abundance_samples]

    # Check for missing samples
    if (any(is.na(Group))) {
      missing_samples <- abundance_samples[is.na(Group)]
      stop("Some samples in abundance data not found in metadata: ",
           paste(missing_samples, collapse = ", "))
    }
  }
  
  required_cols <- c("feature", "group1", "group2", "p_adjust")
  missing_cols <- setdiff(required_cols, colnames(daa_results_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in daa_results_df: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  if (length(Group) != ncol(abundance)) {
    stop("Length of Group must match number of columns in abundance matrix")
  }
  
  # Check for single method
  if (nlevels(factor(daa_results_df$method)) != 1) {
    stop("The 'method' column in the 'daa_results_df' data frame contains ",
         "more than one method. Please filter it to contain only one method.")
  }
  
  if (nlevels(factor(daa_results_df$group1)) != 1 || 
      nlevels(factor(daa_results_df$group2)) != 1) {
    stop("The 'group1' or 'group2' column in the 'daa_results_df' data frame ",
         "contains more than one group. Please filter each to contain only ",
         "one group.")
  }
  
  # Filter for significant features
  daa_results_filtered_df <- daa_results_df[
    daa_results_df$p_adjust < p_values_threshold, 
  ]
  
  if (!is.null(select)) {
    daa_results_filtered_sub_df <- daa_results_filtered_df[
      daa_results_filtered_df$feature %in% select, 
    ]
  } else {
    daa_results_filtered_sub_df <- daa_results_filtered_df
  }
  
  if (nrow(daa_results_filtered_sub_df) > max_features) {
    warning("The number of features with statistical significance exceeds ", 
            max_features, ". Consider using 'select' parameter or increasing ",
            "'max_features' to include all features.")
  }
  
  if (nrow(daa_results_filtered_sub_df) == 0) {
    stop("No features with statistical significance found. ",
         "Consider adjusting p_values_threshold.")
  }
  
  # Get group names
  group1_name <- unique(daa_results_filtered_sub_df$group1)[1]
  group2_name <- unique(daa_results_filtered_sub_df$group2)[1]
  
  # Calculate abundance statistics using the helper function
  # Create metadata that matches the abundance column order
  # The Group vector should be in the same order as colnames(abundance)
  abundance_metadata <- data.frame(
    sample = colnames(abundance),
    group_col = Group,
    stringsAsFactors = FALSE
  )

  # Ensure the Group vector length matches the number of samples
  if (length(Group) != ncol(abundance)) {
    stop("Length of Group vector (", length(Group), ") does not match number of samples in abundance data (", ncol(abundance), ")")
  }

  abundance_stats <- calculate_abundance_stats(
    abundance = abundance,
    metadata = abundance_metadata,
    group = "group_col",
    features = daa_results_filtered_sub_df$feature,
    group1 = group1_name,
    group2 = group2_name
  )
  
  # Merge with DAA results to include p-values and other information
  result_table <- merge(
    abundance_stats,
    daa_results_filtered_sub_df[, c("feature", "p_adjust", "group1", "group2")],
    by = "feature",
    all.x = TRUE
  )
  
  # Reorder columns for better readability
  column_order <- c(
    "feature", "group1", "group2", 
    "mean_rel_abundance_group1", "sd_rel_abundance_group1",
    "mean_rel_abundance_group2", "sd_rel_abundance_group2", 
    "log2_fold_change", "p_adjust"
  )
  
  result_table <- result_table[, column_order]
  
  # Order by p-value (most significant first)
  result_table <- result_table[order(result_table$p_adjust), ]
  
  # Reset row names
  rownames(result_table) <- NULL
  
  return(result_table)
}
