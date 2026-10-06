#' Compare the Consistency of Statistically Significant Features
#'
#' This function compares the consistency and inconsistency of statistically significant features obtained
#' using different methods in `pathway_daa` from the `ggpicrust2` package. It creates a report showing the number of common and
#' different features identified by each method, and the features themselves.
#' @name compare_daa_results
#' @param daa_results_list A list of data frames containing statistically significant features obtained using different methods.
#' @param method_names A character vector of names for each method used.
#' @param p_values_threshold A numeric value representing the threshold for the p-values. Features with p-values less than this
#' threshold are considered statistically significant. Default is 0.05. Must be in the range (0, 1].
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
#' @details
#' Each list element is one discovery set. Split multi-test output such as
#' ALDEx2 by its \code{method} column before comparison; otherwise the tests'
#' discoveries are pooled within that element. This function does not rerun
#' testing or adjust p-values again.
#'
#' For multi-group DAA output, each discovery is compared as a
#' \code{feature + group-pair} unit. The group pair is treated as unordered for
#' this set-level comparison, because the function compares whether methods
#' identified a significant difference, not the effect-size direction. This
#' prevents the same feature from being counted as method-consistent when
#' different methods found it in different pairwise contrasts, while still
#' treating \code{A vs B} and \code{B vs A} as the same biological comparison.
#' Rows whose two group labels are identical are invalid self-comparisons and
#' are rejected.
#' If all significant discoveries share one group pair, the printed feature
#' lists use feature IDs only for backward-readable output; otherwise feature
#' lists include the canonical contrast as \code{feature [group1 vs group2]}.
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

canonical_daa_group_pair <- function(group1, group2) {
  group1 <- as.character(group1)
  group2 <- as.character(group2)
  data.frame(
    group1 = pmin(group1, group2),
    group2 = pmax(group1, group2),
    stringsAsFactors = FALSE
  )
}

format_daa_comparison_units <- function(keys, lookup, include_contrast) {
  if (length(keys) == 0) {
    return("")
  }
  rows <- lookup[match(keys, lookup$key), , drop = FALSE]
  labels <- if (include_contrast) {
    sprintf("%s [%s vs %s]", rows$feature, rows$group1, rows$group2)
  } else {
    rows$feature
  }
  paste(labels, collapse = ", ")
}

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
  validate_probability_threshold(p_values_threshold, "p_values_threshold")
  method_names <- canonicalize_daa_method_names(method_names)
  method_names_chr <- as.character(method_names)
  if (any(is.na(method_names_chr) | !nzchar(trimws(method_names_chr)))) {
    stop("'method_names' must contain non-empty method labels without NA.",
         call. = FALSE)
  }
  if (anyDuplicated(method_names_chr)) {
    duplicated_methods <- unique(method_names_chr[duplicated(method_names_chr)])
    stop("'method_names' contains duplicated labels: ",
         paste(utils::head(duplicated_methods, 5), collapse = ", "),
         ". Use distinct labels so comparison rows are unambiguous.",
         call. = FALSE)
  }

  # Initialize a list to store the feature/contrast discoveries obtained by
  # each method. A DAA row tests a specific pairwise group pair; collapsing
  # across different pairs would overstate cross-method agreement in
  # multi-group analyses. The pair is canonicalized as unordered because this
  # function compares significant discoveries, not effect-size direction.
  feature_units <- vector("list", length(daa_results_list))

  # Get the features obtained by each method
  for (i in seq_along(daa_results_list)) {
    result_i <- daa_results_list[[i]]
    context_i <- paste0("daa_results_list[[", i, "]]")
    validate_dataframe(
      result_i,
      required_cols = c("feature", "group1", "group2", "p_adjust"),
      param_name = context_i
    )
    validate_probability_values(
      result_i$p_adjust,
      "p_adjust",
      context_i
    )
    feature_chr <- validate_nonempty_character_column(result_i$feature,
                                                      "feature", context_i)
    group1_chr <- validate_nonempty_character_column(result_i$group1,
                                                     "group1", context_i)
    group2_chr <- validate_nonempty_character_column(result_i$group2,
                                                     "group2", context_i)
    self_comparison <- group1_chr == group2_chr
    if (any(self_comparison)) {
      stop(
        context_i,
        " contains invalid self-comparison rows where 'group1' equals ",
        "'group2': ",
        paste(utils::head(unique(group1_chr[self_comparison]), 5),
              collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    sig <- !is.na(result_i$p_adjust) & result_i$p_adjust < p_values_threshold
    canonical_pair <- canonical_daa_group_pair(group1_chr[sig],
                                               group2_chr[sig])
    units_i <- data.frame(
      feature = feature_chr[sig],
      group1 = canonical_pair$group1,
      group2 = canonical_pair$group2,
      stringsAsFactors = FALSE
    )
    feature_units[[i]] <- units_i[!duplicated(units_i), , drop = FALSE]
  }

  # Assign opaque integer IDs after deduplicating the exact structured keys.
  # Serializing feature/group columns with a separator can collide when user
  # labels contain that separator; data-frame joins preserve key boundaries.
  unit_lookup <- unique(do.call(rbind, feature_units))
  if (is.null(unit_lookup) || nrow(unit_lookup) == 0) {
    unit_lookup <- data.frame(
      feature = character(0),
      group1 = character(0),
      group2 = character(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(unit_lookup) <- NULL
  unit_lookup$key <- as.character(seq_len(nrow(unit_lookup)))
  unit_lookup <- unit_lookup[, c("key", "feature", "group1", "group2"),
                             drop = FALSE]
  features_flat <- lapply(feature_units, function(units) {
    if (nrow(units) == 0) return(character(0))
    dplyr::left_join(
      units,
      unit_lookup,
      by = c("feature", "group1", "group2")
    )$key
  })
  include_contrast <- nrow(unique(unit_lookup[, c("group1", "group2"), drop = FALSE])) > 1

  # Calculate the intersection and union of the features obtained by each method
  intersect_features <- Reduce(intersect, features_flat)
  union_features <- unique(unlist(features_flat))

  # "diff_features" for method i = features method i reports whose
  # significance is NOT endorsed by every other method. Equivalently:
  # method_i_features \ intersection_of_all_methods. A feature found by
  # 2 of 3 methods is a disagreement from the perspective of each of
  # those 2 methods, so it MUST appear in both of their diff lists.
  #
  # Previously there was a second filter that kept only features
  # appearing in a single method, which silently collapsed "found by
  # 2/3" into "unanimous" and zeroed out the disagreement signal that
  # this function exists to surface.
  diff_features <- lapply(features_flat, function(x) setdiff(x, intersect_features))

  common_features_names <- format_daa_comparison_units(
    intersect_features,
    unit_lookup,
    include_contrast
  )
  comparison_results <- data.frame(
    method = method_names_chr,
    num_features = lengths(features_flat),
    num_common_features = rep(length(intersect_features),
                              length(features_flat)),
    num_diff_features = lengths(diff_features),
    common_features = rep(common_features_names, length(features_flat)),
    diff_features = vapply(
      diff_features,
      format_daa_comparison_units,
      character(1),
      lookup = unit_lookup,
      include_contrast = include_contrast
    ),
    stringsAsFactors = FALSE
  )

  # Report the comparison summary.
  message("Comparing ", length(daa_results_list), " methods:\n")
  unit_noun <- if (include_contrast) "feature/contrast units" else "features"
  for (i in seq_along(daa_results_list)) {
    message("The ", method_names_chr[i], " method obtained ",
            comparison_results$num_features[i],
            " statistically significant ", unit_noun, ".")
    message("The number of ", unit_noun,
            " that are common to other methods is ",
            comparison_results$num_common_features[i])
    message("The number of ", unit_noun,
            " that are different from other methods is ",
            comparison_results$num_diff_features[i])
    message("The names of the ", unit_noun,
            " that are different from other methods are ",
            comparison_results$diff_features[i], "\n")
  }

  message("The number of ", unit_noun,
          " that are common to all methods is ", length(intersect_features))
  message("The names of the ", unit_noun,
          " that are common to all methods are ",
          format_daa_comparison_units(intersect_features, unit_lookup, include_contrast))
  message("The number of ", unit_noun,
          " that are obtained by any of the methods is ", length(union_features))
  message("The names of the ", unit_noun,
          " that are obtained by any of the methods are ",
          format_daa_comparison_units(union_features, unit_lookup, include_contrast))

  comparison_results
}
