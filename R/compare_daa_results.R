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
#'   \item \code{common_features}: The names of the features that are common to all methods.
#'   \item \code{diff_features}: The names of the features that are different from other methods.
#' }
#'
#' @examples
#' # Minimal DAA-like results from three methods (no external dependencies required)
#' deseq2_df <- data.frame(
#'   feature = c("ko00010", "ko00020", "ko00564"),
#'   group1 = c("A", "A", "A"),
#'   group2 = c("B", "B", "B"),
#'   p_adjust = c(0.01, 0.20, 0.03),
#'   stringsAsFactors = FALSE
#' )
#'
#' edgeR_df <- data.frame(
#'   feature = c("ko00010", "ko00680", "ko00564"),
#'   group1 = c("A", "A", "A"),
#'   group2 = c("B", "B", "B"),
#'   p_adjust = c(0.02, 0.04, 0.01),
#'   stringsAsFactors = FALSE
#' )
#'
#' maaslin2_df <- data.frame(
#'   feature = c("ko00010", "ko03030", "ko00564"),
#'   group1 = c("A", "A", "A"),
#'   group2 = c("B", "B", "B"),
#'   p_adjust = c(0.03, 0.02, 0.04),
#'   stringsAsFactors = FALSE
#' )
#'
#' daa_results_list <- list(DESeq2 = deseq2_df, edgeR = edgeR_df, Maaslin2 = maaslin2_df)
#' comparison_results <- compare_daa_results(
#'   daa_results_list = daa_results_list,
#'   method_names = c("DESeq2", "edgeR", "Maaslin2"),
#'   p_values_threshold = 0.05
#' )
#' comparison_results
#' @export
utils::globalVariables(c("group1","group2"))
compare_daa_results <- function(daa_results_list, method_names, p_values_threshold = 0.05) {
  # Compare the consistency and inconsistency of statistically significant features obtained using different methods in pathway_daa.

  # Arguments:
  # daa_results_list -- A list containing statistically significant features obtained using different methods.
  # method_names -- A character vector representing the names of each method.
  # p_values_threshold -- A numeric value representing the threshold for the p-values. Features with p-values less than this threshold are considered statistically significant.


  # Input validation with clear error messages
  if (!is.list(daa_results_list)) {
    stop("'daa_results_list' must be a list of DAA result data frames")
  }
  if (length(daa_results_list) == 0) {
    stop("'daa_results_list' cannot be empty")
  }
  if (!is.character(method_names)) {
    stop("'method_names' must be a character vector")
  }
  if (length(method_names) != length(daa_results_list)) {
    stop(sprintf("'method_names' length (%d) must match 'daa_results_list' length (%d)",
                 length(method_names), length(daa_results_list)))
  }
  if (!is.numeric(p_values_threshold) || length(p_values_threshold) != 1) {
    stop("'p_values_threshold' must be a single numeric value")
  }
  if (p_values_threshold <= 0 || p_values_threshold > 1) {
    stop("'p_values_threshold' must be between 0 and 1")
  }

  # Initialize a list to store the features obtained by each method
  features <- list()

  # Get the features obtained by each method
  for (i in seq_along(daa_results_list)) {
    if (!all(c("group1", "group2") %in% colnames(daa_results_list[[i]]))) {
      stop("DAA results in daa_results_list[[", i, "]] must contain 'group1' and 'group2' columns")
    }
    group_combinations <- unique(daa_results_list[[i]][, c("group1", "group2")])
    features[[i]] <- lapply(seq_len(nrow(group_combinations)), function(j) {
      subset(daa_results_list[[i]], group1 == group_combinations[j, "group1"] & group2 == group_combinations[j, "group2"] & p_adjust < p_values_threshold)$feature
    })
  }

  # Flatten nested list structure (from group combinations) to simple vectors
  features_flat <- lapply(features, function(x) unique(unlist(x)))

  # Calculate the intersection and union of the features obtained by each method
  intersect_features <- Reduce(intersect, features_flat)
  union_features <- unique(unlist(features_flat))

  # Calculate the differences in the features obtained by each method (exclude common)
  diff_features <- lapply(features_flat, function(x) setdiff(x, intersect_features))

  # Find features unique to each method (not shared with any other method)
  # Using duplicated() for correct handling of features appearing in 3+ methods
  all_diff <- unlist(diff_features)
  is_duplicated <- duplicated(all_diff) | duplicated(all_diff, fromLast = TRUE)
  unique_only_features <- unique(all_diff[!is_duplicated])
  diff_features <- lapply(diff_features, function(x) x[x %in% unique_only_features])

  # Initialize a data frame to store the comparison results
  comparison_results <- data.frame(
    method = character(),
    num_features = integer(),
    num_common_features = integer(),
    num_diff_features = integer(),
    common_features = character(),
    diff_features = character(),
    stringsAsFactors = FALSE
  )

  # Output the comparison results and store them in the data frame
  message("Comparing ", length(daa_results_list), " methods:\n")
  for (i in seq_along(daa_results_list)) {
    num_features <- length(features_flat[[i]])
    num_common_features <- length(intersect_features)
    num_diff_features <- length(diff_features[[i]])
    common_features_names <- paste(intersect_features, collapse = ", ")
    diff_features_names <- paste(diff_features[[i]], collapse = ", ")

    message("The ", method_names[i], " method obtained ", num_features, " statistically significant features.")
    message("The number of features that are common to other methods is ", num_common_features)
    message("The number of features that are different from other methods is ", num_diff_features)
    message("The names of the features that are different from other methods are ", diff_features_names, "\n")

    comparison_results <- rbind(comparison_results, data.frame(
      method = method_names[i],
      num_features = num_features,
      num_common_features = num_common_features,
      num_diff_features = num_diff_features,
      common_features = common_features_names,
      diff_features = diff_features_names,
      stringsAsFactors = FALSE
    ))
  }

  message("The number of features that are common to all methods is ", length(intersect_features))
  message("The names of the features that are common to all methods are ", paste(intersect_features, collapse = ", "))
  message("The number of features that are obtained by any of the methods is ", length(union_features))
  message("The names of the features that are obtained by any of the methods are ", paste(union_features, collapse = ", "))

  return(comparison_results)
}
