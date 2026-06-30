#' Import Differential Abundance Analysis (DAA) results from MicrobiomeAnalyst
#'
#' This function imports DAA results from an external platform such as MicrobiomeAnalyst. It can be used to compare the results obtained from different platforms.
#' @name import_MicrobiomeAnalyst_daa_results
#' @param file_path a character string specifying the path to the CSV file containing the DAA results from MicrobiomeAnalyst. If this parameter is NULL and no data frame is provided, an error will be thrown. Default is NULL.
#' @param data a data frame containing the DAA results from MicrobiomeAnalyst. Feature identifiers can be stored in a feature/name column or in non-default row names. P-value and adjusted p-value columns are detected by semantic names such as \code{Pvalues} and \code{FDR}; \code{Statistics} and fold-change columns are optional. If this parameter is NULL and no file path is provided, an error will be thrown. Default is NULL.
#' @param method a single non-empty character string specifying the method used for the DAA. This will be added as a new column in the returned data frame. Default is "MicrobiomeAnalyst".
#' @param group_levels a character vector specifying at least two unique group levels for the DAA. These values will be added as new columns in the returned data frame. Default is c("control", "treatment").
#'
#' @return a data frame containing the DAA results from MicrobiomeAnalyst with
#'   validated \code{feature}, \code{p_values}, and \code{p_adjust} columns,
#'   optional \code{Statistics} and \code{log2_fold_change} columns when present
#'   in the imported result, plus additional columns for the method and group
#'   levels.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a CSV file named "DAA_results.csv" in your current directory
#' daa_results <- import_MicrobiomeAnalyst_daa_results(file_path = "DAA_results.csv")
#' }
#'
#' @export
import_MicrobiomeAnalyst_daa_results <- function(file_path = NULL,
                                                data = NULL,
                                                method = "MicrobiomeAnalyst",
                                                group_levels = c("control", "treatment")) {
  if (!is.null(data) && !is.null(file_path)) {
    warning("Both file_path and data were provided. Using data and ignoring file_path.",
            call. = FALSE)
  }

  if (is.null(data)) {
    if (is.null(file_path)) {
      stop("Please provide either a file_path or a data frame.")
    }
    if (!file.exists(file_path)) {
      stop("file_path does not exist: ", file_path, call. = FALSE)
    }
    data <- utils::read.csv(file_path,
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
  }

  validate_dataframe(data, param_name = "data")

  if (!is.character(method) || length(method) != 1 ||
      is.na(method) || !nzchar(method)) {
    stop("method must be a single non-empty character value.",
         call. = FALSE)
  }

  if (!is.character(group_levels) || length(group_levels) < 2 ||
      anyNA(group_levels) || any(!nzchar(group_levels)) ||
      anyDuplicated(group_levels)) {
    stop("group_levels must contain at least two unique, non-empty character values without NA values.",
         call. = FALSE)
  }

  data <- standardize_microbiomeanalyst_columns(data)

  if (anyNA(data$feature) || any(!nzchar(as.character(data$feature)))) {
    stop("Column 'feature' must contain non-empty feature identifiers.",
         call. = FALSE)
  }
  if (anyDuplicated(as.character(data$feature))) {
    dup_features <- unique(as.character(data$feature[duplicated(as.character(data$feature))]))
    stop("Column 'feature' contains duplicated identifiers: ",
         paste(utils::head(dup_features, 5), collapse = ", "),
         call. = FALSE)
  }

  numeric_cols <- intersect(
    c("p_values", "p_adjust", "Statistics", "log2_fold_change"),
    colnames(data)
  )
  for (col in numeric_cols) {
    if (!is.numeric(data[[col]])) {
      parsed <- suppressWarnings(as.numeric(data[[col]]))
      non_missing_input <- !is.na(data[[col]])
      introduced_na <- is.na(parsed) & non_missing_input
      if (any(introduced_na)) {
        stop("Column '", col, "' must be numeric or numeric-like.",
             call. = FALSE)
      }
      data[[col]] <- parsed
    }
  }

  validate_probability_values(data$p_values, "p_values", "data")
  validate_probability_values(data$p_adjust, "p_adjust", "data")
  if ("Statistics" %in% colnames(data)) {
    validate_finite_numeric_values(data$Statistics, "Statistics", "data")
  }
  if ("log2_fold_change" %in% colnames(data)) {
    validate_finite_numeric_values(data$log2_fold_change,
                                   "log2_fold_change",
                                   "data")
  }

  # Create a new column for method
  data$method <- method

  # Create new columns for group levels
  for (i in seq_along(group_levels)) {
    data[paste0("group", i)] <- group_levels[i]
  }

  return(data)
}

normalize_microbiomeanalyst_column_name <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", x))
}

resolve_microbiomeanalyst_column <- function(data,
                                             aliases,
                                             target,
                                             required = TRUE,
                                             allow_blank_first = FALSE) {
  column_names <- names(data)
  normalized_names <- normalize_microbiomeanalyst_column_name(column_names)
  normalized_aliases <- normalize_microbiomeanalyst_column_name(aliases)
  matches <- which(normalized_names %in% normalized_aliases)

  if (allow_blank_first && length(matches) == 0 && ncol(data) > 0 &&
      (!nzchar(column_names[1]) || is.na(column_names[1]))) {
    matches <- 1L
  }

  if (length(matches) == 0) {
    if (required) {
      stop("MicrobiomeAnalyst DAA results must contain a column for '",
           target, "'. Recognized names include: ",
           paste(aliases, collapse = ", "),
           call. = FALSE)
    }
    return(NA_integer_)
  }

  if (length(matches) > 1) {
    stop("Ambiguous MicrobiomeAnalyst DAA columns for '", target, "': ",
         paste(column_names[matches], collapse = ", "),
         ". Keep only one semantic column before importing.",
         call. = FALSE)
  }

  matches
}

has_non_default_rownames <- function(data) {
  rn <- rownames(data)
  !is.null(rn) &&
    length(rn) == nrow(data) &&
    !identical(rn, as.character(seq_len(nrow(data)))) &&
    !anyNA(rn) &&
    all(nzchar(rn))
}

standardize_microbiomeanalyst_columns <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  names(data) <- as.character(names(data))

  feature_aliases <- c(
    "feature", "features", "id", "ids", "name", "names",
    "rowname", "rownames", "otu", "otunames", "otu names",
    "ko", "kegg", "pathway", "pathway_id"
  )
  p_value_aliases <- c(
    "p_values", "pvalues", "pvalue", "p-value", "p.value",
    "pval", "pvals", "raw_p", "rawp", "p"
  )
  p_adjust_aliases <- c(
    "p_adjust", "padjust", "p.adjust", "fdr", "adjpvalues",
    "adjpvalue", "adj.p.values", "adj.p.value", "padj",
    "qvalue", "qvalues", "adjustedpvalue", "adjustedpvalues"
  )
  statistic_aliases <- c(
    "Statistics", "Statistic", "statistics", "statistic",
    "stats", "stat", "tstat", "tstatistic"
  )
  log2fc_aliases <- c(
    "log2_fold_change", "log2foldchange", "log2FC", "logFC",
    "log2fc", "logfc", "logFoldChange", "log_fold_change"
  )

  feature_col <- resolve_microbiomeanalyst_column(
    data,
    feature_aliases,
    "feature",
    required = FALSE,
    allow_blank_first = TRUE
  )
  if (is.na(feature_col)) {
    if (!has_non_default_rownames(data)) {
      stop("MicrobiomeAnalyst DAA results must contain feature identifiers ",
           "either in a feature/name column or in non-default row names.",
           call. = FALSE)
    }
    data <- data.frame(
      feature = rownames(data),
      data,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  semantic_columns <- list(
    feature = resolve_microbiomeanalyst_column(
      data, feature_aliases, "feature",
      required = TRUE,
      allow_blank_first = TRUE
    ),
    p_values = resolve_microbiomeanalyst_column(
      data, p_value_aliases, "p_values", required = TRUE
    ),
    p_adjust = resolve_microbiomeanalyst_column(
      data, p_adjust_aliases, "p_adjust", required = TRUE
    ),
    Statistics = resolve_microbiomeanalyst_column(
      data, statistic_aliases, "Statistics", required = FALSE
    ),
    log2_fold_change = resolve_microbiomeanalyst_column(
      data, log2fc_aliases, "log2_fold_change", required = FALSE
    )
  )

  for (target in names(semantic_columns)) {
    idx <- semantic_columns[[target]]
    if (!is.na(idx)) {
      names(data)[idx] <- target
    }
  }

  data
}
