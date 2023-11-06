#' Compare the Consistency of Statistically Significant Features
#'
#' This function compares the consistency and inconsistency of statistically significant features obtained
#' using different methods in `pathway_daa` from the `ggpicrust2` package. It creates a report showing the number of common and
#' different features identified by each method, and the features themselves.
#' @name compare_daa_results
#' @param daa_results_list A list of data frames containing statistically significant features obtained using different methods.
#' @param method_names A character vector of names for each method used.
#' @param p_values_threshold A numeric value representing the threshold for the p-values. Features with p-values less than this
#' threshold are considered statistically significant. Default is 0.05.
#'
#' @return A data frame with the comparison results. The data frame has the following columns:
#' \itemize{
#'   \item \code{method}: The name of the method.
#'   \item \code{num_features}: The total number of statistically significant features obtained by the method.
#'   \item \code{num_common_features}: The number of features that are common to other methods.
#'   \item \code{num_diff_features}: The number of features that are different from other methods.
#'   \item \code{diff_features}: The names of the features that are different from other methods.
#' }
#'
#' @examples
#' \donttest{
#' library(magrittr)
#' library(ggpicrust2)
#' library(tibble)
#' data("metacyc_abundance")
#' data("metadata")
#'
#' # Run pathway_daa function for multiple methods
#' methods <- c("DESeq2", "edgeR","Maaslin2")
#' daa_results_list <- lapply(methods, function(method) {
#' pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"),
#' metadata = metadata, group = "Environment", daa_method = method)
#' })
#'
#' names(daa_results_list) <- methods
#' # Correct Maaslin2 feature names by replacing dots with hyphens.
#' # Note: When using Maaslin2 as the differential abundance analysis method,
#' # it modifies the original feature names by replacing hyphens (-) with dots (.).
#' # This replacement can cause inconsistencies when trying to compare results from Maaslin2
#' # with those from other methods that do not modify feature names.
#' # Therefore, this line of code reverses that replacement, converting the dots back into
#' # hyphens for accurate and consistent comparisons across different methods.
#' daa_results_list[["Maaslin2"]]$feature <- gsub("\\.", "-", daa_results_list[["Maaslin2"]]$feature)
#'
#' # Compare results across different methods
#' comparison_results <- compare_daa_results(daa_results_list = daa_results_list,
#' method_names = c("DESeq2", "edgeR", "Maaslin2"))
#' }
#' @export
utils::globalVariables(c("group1","group2"))
compare_daa_results <- function(daa_results_list, method_names, p_values_threshold = 0.05) {
  # Compare the consistency and inconsistency of statistically significant features obtained using different methods in pathway_daa.

  # Arguments:
  # daa_results_list -- A list containing statistically significant features obtained using different methods.
  # method_names -- A character vector representing the names of each method.
  # p_values_threshold -- A numeric value representing the threshold for the p-values. Features with p-values less than this threshold are considered statistically significant.

  # Check input parameters
  stopifnot(is.list(daa_results_list))
  stopifnot(is.character(method_names))
  stopifnot(is.numeric(p_values_threshold))

  # Initialize a list to store the features obtained by each method
  features <- list()

  # Get the features obtained by each method
  for (i in seq_along(daa_results_list)) {
    if (all(c("group1", "group2") %in% colnames(daa_results_list[[i]]))) {
      # This is a pairwise comparison method
      group_combinations <- unique(daa_results_list[[i]][, c("group1", "group2")])
      features[[i]] <- lapply(seq_len(nrow(group_combinations)), function(j) {
        subset(daa_results_list[[i]], group1 == group_combinations[j, "group1"] & group2 == group_combinations[j, "group2"] & p_adjust < p_values_threshold)$feature
      })
    } else if (all(c("group1", "group2", "group3") %in% colnames(daa_results_list[[i]]))) {
      # This is a multi-group comparison method
      features[[i]] <- daa_results_list[[i]][daa_results_list[[i]]$p_adjust < p_values_threshold, "feature"]
    } else {
      stop("Unknown comparison type in daa_results_list[[", i, "]].")
    }
  }

  # Calculate the intersection and union of the features obtained by each method
  intersect_features <- Reduce(intersect, unlist(features, recursive = FALSE))
  union_features <- unique(unlist(features))

  # Calculate the differences in the features obtained by each method
  diff_features <- lapply(features, function(x) setdiff(x, intersect_features))

  # Compare features across different groups using the statistical principle
  for (i in seq_along(diff_features)) {
    for (j in seq_along(diff_features)) {
      if (i != j) {
        common_features <- intersect(diff_features[[i]], diff_features[[j]])
        if (length(common_features) > 0) {
          diff_features[[i]] <- setdiff(diff_features[[i]], common_features)
          diff_features[[j]] <- setdiff(diff_features[[j]], common_features)
        }
      }
    }
  }

  # Initialize a data frame to store the comparison results
  comparison_results <- data.frame(
    method = character(),
    num_features = integer(),
    num_common_features = integer(),
    num_diff_features = integer(),
    diff_features = character(),
    stringsAsFactors = FALSE
  )

  # Output the comparison results and store them in the data frame
  cat("Comparing", length(daa_results_list), "methods:\n\n")
  for (i in seq_along(daa_results_list)) {
    num_features <- length(unlist(features[[i]]))
    num_common_features <- length(intersect_features)
    num_diff_features <- length(unlist(diff_features[[i]]))
    diff_features_names <- paste(unlist(diff_features[[i]]), collapse = ", ")

    cat("The", method_names[i], "method obtained", num_features, "statistically significant features.\n")
    cat("The number of features that are common to other methods is", num_common_features, "\n")
    cat("The number of features that are different from other methods is", num_diff_features, "\n")
    cat("The names of the features that are different from other methods are", diff_features_names, "\n\n")

    comparison_results <- rbind(comparison_results, data.frame(
      method = method_names[i],
      num_features = num_features,
      num_common_features = num_common_features,
      num_diff_features = num_diff_features,
      diff_features = diff_features_names,
      stringsAsFactors = FALSE
    ))
  }

  cat("The number of features that are common to all methods is", length(intersect_features), "\n")
  cat("The names of the features that are common to all methods are", paste(intersect_features, collapse = ", "), "\n")
  cat("The number of features that are obtained by any of the methods is", length(union_features), "\n")
  cat("The names of the features that are obtained by any of the methods are", paste(union_features, collapse = ", "), "\n")

  return(comparison_results)
}
