#' Data Utilities for ggpicrust2
#'
#' Internal utility functions for data preprocessing, sample matching,
#' and format standardization.
#'
#' @name data_utils
#' @keywords internal
NULL

# =============================================================================
# File I/O
# =============================================================================

#' Read abundance data from a delimited file
#'
#' Reads PICRUSt2-style abundance data from .tsv, .txt, or .csv files.
#' Delimiter is auto-detected from the file extension (csv -> comma, others -> tab).
#'
#' @param file_path Path to the input file
#' @return A data frame (tibble) with the parsed abundance data
#' @noRd
read_abundance_file <- function(file_path) {
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("'file_path' must be a single file path string")
  }
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  file_ext <- tolower(tools::file_ext(file_path))
  if (!file_ext %in% c("txt", "tsv", "csv")) {
    stop(
      "Unsupported file format '.", file_ext, "'. ",
      "Accepted formats: .tsv, .txt, .csv"
    )
  }

  delimiter <- if (file_ext == "csv") "," else "\t"

  abundance <- readr::read_delim(
    file_path,
    delim = delimiter,
    show_col_types = FALSE,
    progress = FALSE
  )

  if (ncol(abundance) < 2) {
    stop("Input file must contain at least 2 columns (IDs + samples)")
  }

  abundance
}

# =============================================================================
# KO Abundance Cleaning
# =============================================================================

#' Clean KO Abundance Data
#'
#' Standardizes KO ID format by removing common prefixes (e.g., "ko:" from
#' PICRUSt 2.6.2 format). This function is idempotent - calling it multiple
#' times produces the same result.
#'
#' @param abundance Data frame with KO IDs in first column
#' @param verbose Logical, whether to print messages (default: TRUE)
#' @return Data frame with cleaned KO IDs
#' @noRd
clean_ko_abundance <- function(abundance, verbose = TRUE) {
  if (!is.data.frame(abundance) || ncol(abundance) < 2) {
    return(abundance)
  }

  ko_ids <- abundance[[1]]
  original_ids <- ko_ids

  # Remove "ko:" prefix (PICRUSt 2.6.2 format)
  ko_ids <- gsub("^ko:", "", ko_ids)

  # Count changes
  n_cleaned <- sum(original_ids != ko_ids)

  if (n_cleaned > 0 && verbose) {
    message(sprintf("Standardized %d KO IDs (removed 'ko:' prefix)", n_cleaned))
  }

  abundance[[1]] <- ko_ids
  abundance
}

# =============================================================================
# Sample Matching Utilities
# =============================================================================
#' Check if Two Sample Vectors Have Sufficient Overlap
#'
#' @param vec1 First vector of sample identifiers
#' @param vec2 Second vector of sample identifiers
#' @param threshold Minimum proportion of overlap required (default: 0.5)
#' @return Logical indicating whether vectors have sufficient overlap
#' @noRd
samples_match <- function(vec1, vec2, threshold = 0.5, require_unique = FALSE) {
  vec1 <- as.character(vec1)
  vec2 <- as.character(vec2)
  n_common <- length(intersect(vec1, vec2))
  min_length <- min(length(vec1), length(vec2))

  if (min_length == 0) return(FALSE)

  # A true sample identifier column has one row per sample, so its values
  # must be unique. Without this guard, a low-cardinality categorical column
  # (e.g. "group", "batch") whose levels coincidentally include a few sample
  # names can pass the overlap threshold.
  if (require_unique && anyDuplicated(vec1)) return(FALSE)

  n_common / min_length >= threshold
}

#' Find Sample Column in Metadata
#'
#' Identifies which column in metadata contains sample identifiers that match
#' the abundance data columns. Searches in priority order: standard names,
#' then all columns, then rownames.
#'
#' @param metadata Data frame containing sample metadata
#' @param abundance_samples Character vector of sample names from abundance data
#' @return Column name if found, or ".rownames" if rownames match, or NULL
#' @noRd
find_sample_column <- function(metadata, abundance_samples) {
  abundance_samples <- as.character(abundance_samples)

  # Priority 1: Standard column names (case variations)
  standard_names <- c(
    "sample", "Sample", "SAMPLE",
    "sample_name", "Sample_Name", "SampleName", "samplename",
    "sample_id", "Sample_ID", "SampleID", "sampleid",
    "samples", "Samples"
  )

  # Priority 1 uses a more lenient overlap threshold (50%) because a
  # standard column name ("sample", "sample_id", etc.) is itself strong
  # evidence that this is the sample identifier -- partial overlap is
  # usually just metadata filtering upstream. But uniqueness is not
  # optional: a "sample" column with duplicate values is not a sample ID
  # no matter what it's called, so enforce require_unique = TRUE here
  # too to prevent mis-detection on metadata whose standard-named column
  # happens to be a subject_id / replicate_id with repeats.
  for (col in standard_names) {
    if (col %in% colnames(metadata)) {
      if (samples_match(metadata[[col]], abundance_samples,
                        require_unique = TRUE)) {
        return(col)
      }
    }
  }

  # Priority 2 scans every column without a name signal, so we require much
  # stronger evidence: (a) column values must be unique within metadata (a
  # real sample ID is one-per-row), and (b) near-complete overlap with the
  # abundance columns (>= 90%). A loose 50% threshold here can silently
  # pick up a subject_id or batch column whose levels happen to share a
  # few strings with the sample names.
  for (col in colnames(metadata)) {
    if (samples_match(metadata[[col]], abundance_samples,
                      threshold = 0.9, require_unique = TRUE)) {
      return(col)
    }
  }

  # Priority 3: rownames -- also require uniqueness and high overlap, since
  # default integer rownames ("1", "2", ...) should not be mistaken for
  # sample identifiers.
  if (samples_match(rownames(metadata), abundance_samples,
                    threshold = 0.9, require_unique = TRUE)) {
    return(".rownames")
  }

  NULL
}

#' Align Abundance Data and Metadata by Samples
#'
#' Ensures abundance columns and metadata rows are in the same sample order.
#' Handles sample column detection, validation, and reordering in one step.
#'
#' @param abundance Data frame or matrix with samples as columns (and features as rows)
#' @param metadata Data frame with sample metadata
#' @param sample_col Column name containing sample IDs (auto-detected if NULL)
#' @param verbose Logical, whether to print messages (default: TRUE)
#' @return List containing:
#'   \item{abundance}{Aligned abundance data}
#'   \item{metadata}{Aligned metadata}
#'   \item{sample_col}{Name of the sample column used}
#'   \item{n_samples}{Number of aligned samples}
#' @noRd
align_samples <- function(abundance, metadata, sample_col = NULL, verbose = TRUE) {
  # Get abundance sample names
  abundance_samples <- colnames(abundance)

  if (length(abundance_samples) == 0) {
    stop("Abundance data has no column names (sample identifiers)")
  }

  # Find sample column if not provided
  if (is.null(sample_col)) {
    sample_col <- find_sample_column(metadata, abundance_samples)

    if (is.null(sample_col)) {
      stop(
        "Cannot find matching sample identifiers between abundance and metadata.\n",
        "  Abundance samples (first 5): ",
        paste(head(abundance_samples, 5), collapse = ", "),
        if (length(abundance_samples) > 5) "..." else "", "\n",
        "  Metadata columns: ", paste(colnames(metadata), collapse = ", "), "\n",
        "  Metadata rownames (first 5): ",
        paste(head(rownames(metadata), 5), collapse = ", "),
        if (nrow(metadata) > 5) "..." else ""
      )
    }

    if (verbose) {
      if (sample_col == ".rownames") {
        message("Using metadata rownames as sample identifiers")
      } else {
        message(sprintf("Using column '%s' as sample identifier", sample_col))
      }
    }
  }

  # Handle rownames case
  if (sample_col == ".rownames") {
    metadata <- as.data.frame(metadata)
    metadata$.sample_id <- rownames(metadata)
    sample_col <- ".sample_id"
  }

  # Validate sample_col exists
 if (!sample_col %in% colnames(metadata)) {
    stop(sprintf("Sample column '%s' not found in metadata", sample_col))
  }

  metadata_samples <- as.character(metadata[[sample_col]])

  # A sample identifier column is a primary key: one row per sample.
  # find_sample_column() already enforces uniqueness when auto-detecting,
  # but users who pass `sample_col` explicitly bypass that path. Without
  # this guard, duplicated IDs silently collapse in match() below
  # (match() returns only the first hit), so a metadata row is dropped
  # without warning and the caller gets stats on the wrong sample count.
  dup_ids <- unique(metadata_samples[duplicated(metadata_samples)])
  if (length(dup_ids) > 0) {
    stop(sprintf(
      "Sample column '%s' contains duplicated IDs: %s. Each metadata row must correspond to a unique sample.",
      sample_col,
      paste(head(dup_ids, 5), collapse = ", ")
    ), call. = FALSE)
  }

  # Find common samples
  common_samples <- intersect(abundance_samples, metadata_samples)

  if (length(common_samples) == 0) {
    stop(
      "No common samples found between abundance and metadata.\n",
      "  Abundance samples: ", paste(head(abundance_samples, 3), collapse = ", "), "...\n",
      "  Metadata samples: ", paste(head(metadata_samples, 3), collapse = ", "), "..."
    )
  }

  # Report mismatches if any
  missing_in_abundance <- setdiff(metadata_samples, abundance_samples)
  missing_in_metadata <- setdiff(abundance_samples, metadata_samples)

  if (verbose && (length(missing_in_abundance) > 0 || length(missing_in_metadata) > 0)) {
    message(sprintf(
      "Aligned %d samples (dropped %d from metadata, %d from abundance)",
      length(common_samples),
      length(missing_in_abundance),
      length(missing_in_metadata)
    ))
  }

  # Align data: subset and reorder to common samples
  abundance_aligned <- abundance[, common_samples, drop = FALSE]
  row_idx <- match(common_samples, metadata_samples)
  metadata_aligned <- as.data.frame(metadata[row_idx, , drop = FALSE])

  # Set rownames to sample names for downstream compatibility
  rownames(metadata_aligned) <- common_samples

  # Verify alignment
  if (!identical(colnames(abundance_aligned), as.character(metadata_aligned[[sample_col]]))) {
    stop("Internal error: sample alignment failed")
  }

  list(
    abundance = abundance_aligned,
    metadata = metadata_aligned,
    sample_col = sample_col,
    n_samples = length(common_samples)
  )
}

# =============================================================================
# Reference Data Loading
# =============================================================================

#' Load Reference Data from Package
#'
#' Unified function to load reference data from the package's data/ directory.
#' All reference data is lazy-loaded, so this function simply retrieves it
#' from the package namespace.
#'
#' @param ref_type Character string specifying the reference type:
#'   \itemize{
#'     \item "KO": KO annotation reference (id, description, pathway info)
#'     \item "EC": EC annotation reference (id, description)
#'     \item "MetaCyc": MetaCyc pathway reference (id, description)
#'     \item "KEGG": KEGG pathway reference (pathway, pathway_name)
#'     \item "ko_to_kegg": KO to KEGG pathway mapping
#'     \item "ko_to_go": KO to GO term mapping
#'     \item "metacyc_to_ec": MetaCyc to EC mapping
#'   }
#' @return Data frame with reference data
#' @noRd
load_reference_data <- function(ref_type) {
  # Mapping from user-friendly names to actual data object names

  ref_map <- list(
    "KO" = "ko_reference",
    "EC" = "ec_reference",
    "MetaCyc" = "metacyc_reference",
    "KEGG" = "kegg_pathway_reference",
    "ko_to_kegg" = "ko_to_kegg_reference",
    "ko_to_go" = "ko_to_go_reference",
    "metacyc_to_ec" = "metacyc_to_ec_reference"
  )

  ref_name <- ref_map[[ref_type]]
  if (is.null(ref_name)) {
    stop(sprintf(
      "Unknown reference type: '%s'. Valid options: %s",
      ref_type, paste(names(ref_map), collapse = ", ")
    ))
  }

  # Get from package namespace (lazy-loaded data)
  pkg_ns <- asNamespace("ggpicrust2")
  if (exists(ref_name, envir = pkg_ns)) {
    return(as.data.frame(get(ref_name, envir = pkg_ns)))
  }

  # Fallback: try using data() function

  env <- new.env()
 tryCatch({
    data(list = ref_name, package = "ggpicrust2", envir = env)
    if (exists(ref_name, envir = env)) {
      return(as.data.frame(get(ref_name, envir = env)))
    }
  }, error = function(e) NULL)

  stop(sprintf(
    "Reference data '%s' not found. Please reinstall the package.",
    ref_name
  ))
}

# =============================================================================
# Statistical Calculation Utilities
# =============================================================================

#' Calculate data-driven pseudocount for log transformation
#'
#' Calculates a pseudocount based on the data to avoid log(0) issues.
#' Uses half of the minimum non-zero value to ensure the pseudocount
#' is smaller than any real value in the data.
#'
#' @param values Numeric vector of abundance values
#' @return Pseudocount value (half of minimum non-zero value, or 1e-6 fallback)
#' @keywords internal
calculate_pseudocount <- function(values) {
  values <- as.numeric(values)
  non_zero <- values[values > 0 & !is.na(values) & !is.infinite(values)]
  if (length(non_zero) > 0) {
    min(non_zero) * 0.5
  } else {
    1e-6  # fallback for edge cases
 }
}

#' Calculate log2 fold change with consistent pseudocount handling
#'
#' Calculates log2 fold change between two groups with proper handling
#' of zero values using a data-driven pseudocount approach.
#'
#' @param mean1 Numeric. Mean abundance of group 1 (reference/control)
#' @param mean2 Numeric. Mean abundance of group 2 (comparison/treatment)
#' @param pseudocount Optional numeric. If NULL, calculated from reference_values
#' @param reference_values Optional numeric vector for calculating pseudocount
#' @return log2(mean2/mean1) with pseudocount protection
#' @keywords internal
#'
#' @details
#' The fold change direction is group2/group1, so:
#' - Positive values indicate higher abundance in group2
#' - Negative values indicate higher abundance in group1
calculate_log2_fold_change <- function(mean1, mean2, pseudocount = NULL,
                                        reference_values = NULL) {
  # Calculate pseudocount if not provided
  if (is.null(pseudocount)) {
    if (!is.null(reference_values)) {
      pseudocount <- calculate_pseudocount(reference_values)
    } else {
      # Use the means themselves to estimate pseudocount
      pseudocount <- calculate_pseudocount(c(mean1, mean2))
    }
  }

  # Calculate log2 fold change: group2 / group1
  log2((mean2 + pseudocount) / (mean1 + pseudocount))
}

# =============================================================================
# Input Validation Utilities
# =============================================================================

#' Validate Abundance Data
#'
#' Validates that abundance data is in the correct format for analysis.
#'
#' @param abundance Data frame or matrix to validate
#' @param min_samples Minimum number of samples required (default: 2)
#' @param require_numeric Whether all non-ID columns must be numeric (default: TRUE)
#' @param check_zero_columns If TRUE (default), reject sample columns whose
#'   total abundance is zero. Such samples would turn into NaN inside
#'   `x / sum(x)` and then be silently dropped by downstream
#'   `mean(..., na.rm = TRUE)`, producing stats computed from the wrong
#'   sample size without any warning.
#' @return TRUE if valid, otherwise stops with error
#' @noRd
validate_abundance <- function(abundance, min_samples = 2, require_numeric = TRUE,
                               check_zero_columns = TRUE) {
  if (!is.data.frame(abundance) && !is.matrix(abundance)) {
    stop("'abundance' must be a data frame or matrix")
  }

  if (ncol(abundance) < min_samples) {
    stop(sprintf("At least %d samples are required, found %d", min_samples, ncol(abundance)))
  }

  if (require_numeric && is.data.frame(abundance)) {
    # Check if numeric columns exist (excluding potential ID column)
    numeric_cols <- sapply(abundance, is.numeric)
    if (sum(numeric_cols) == 0) {
      stop("Abundance data must contain numeric values")
    }
  }

  if (check_zero_columns) {
    numeric_mat <- abundance_to_numeric_matrix(abundance)
    if (!is.null(numeric_mat)) {
      col_totals <- colSums(numeric_mat, na.rm = FALSE)
      bad <- which(!is.finite(col_totals) | col_totals == 0)
      if (length(bad) > 0) {
        offenders <- column_labels(numeric_mat, bad)
        stop(sprintf(
          "Invalid abundance data: %d sample column(s) have a total abundance of 0 or NA: %s.\n  Each sample must contain at least one non-zero feature. Drop or reprocess these samples before analysis.",
          length(bad),
          paste(head(offenders, 5), collapse = ", ")
        ), call. = FALSE)
      }
    }
  }

  TRUE
}

#' Extract the numeric sample-by-feature matrix from abundance input
#'
#' `validate_abundance()` accepts either a matrix or a data.frame that may
#' carry a non-numeric ID column. For the zero-column check we only care
#' about the numeric part; this helper normalizes both shapes.
#'
#' Returns NULL when no numeric columns exist (the caller has already
#' errored via `require_numeric` when appropriate).
#' @noRd
abundance_to_numeric_matrix <- function(abundance) {
  if (is.matrix(abundance)) {
    if (!is.numeric(abundance)) return(NULL)
    return(abundance)
  }
  # data.frame
  numeric_cols <- vapply(abundance, is.numeric, logical(1))
  if (!any(numeric_cols)) return(NULL)
  as.matrix(abundance[, numeric_cols, drop = FALSE])
}

#' Return human-readable column labels for a set of positions
#'
#' Prefer column names when they exist and are non-empty; otherwise fall
#' back to the position index so the error is still actionable.
#' @noRd
column_labels <- function(x, positions) {
  nm <- colnames(x)
  if (is.null(nm)) return(as.character(positions))
  labs <- nm[positions]
  labs[!nzchar(labs)] <- as.character(positions[!nzchar(labs)])
  labs
}

#' Compute Relative Abundance with Zero-Column Guard
#'
#' Divides each sample column by its column sum to produce relative
#' abundances. Replaces the duplicated `apply(., 2, function(x) x / sum(x))`
#' idiom and -- more importantly -- refuses zero-sum columns instead of
#' silently producing NaN that later gets dropped by `na.rm = TRUE` in
#' downstream aggregations (which would compute stats on fewer samples
#' than the user expects without any indication).
#'
#' @param abundance Numeric matrix or numeric-only data frame with features
#'   as rows and samples as columns.
#' @param context Short string describing the caller, used in the error
#'   message so the user knows where the zero-sum sample was detected.
#' @return A matrix of relative abundances with the same shape as the
#'   input.
#' @noRd
compute_relative_abundance <- function(abundance, context = "abundance") {
  mat <- if (is.matrix(abundance)) abundance else as.matrix(abundance)
  if (!is.numeric(mat)) {
    stop(sprintf("Cannot compute relative abundance: %s is not numeric.", context),
         call. = FALSE)
  }
  col_totals <- colSums(mat, na.rm = FALSE)
  bad <- which(!is.finite(col_totals) | col_totals == 0)
  if (length(bad) > 0) {
    offenders <- column_labels(mat, bad)
    stop(sprintf(
      "Cannot compute relative abundance for %s: %d sample column(s) have a total of 0 or NA: %s.\n  Drop or reprocess these samples before analysis.",
      context, length(bad),
      paste(head(offenders, 5), collapse = ", ")
    ), call. = FALSE)
  }
  sweep(mat, 2, col_totals, "/")
}

#' Summarize Per-Feature Mean and SD by Group
#'
#' Given a relative-abundance matrix (features x samples, already normalized)
#' and a parallel vector of per-sample group assignments, return a long-format
#' data frame with one row per (feature, group) combination containing the
#' mean and standard deviation across the samples in that group.
#'
#' This helper is the single source of truth for the per-feature-per-group
#' mean/sd that both `calculate_abundance_stats()` (used by `pathway_daa()`
#' and `pathway_errorbar_table()`) and `pathway_errorbar()` need. Previously
#' `pathway_errorbar()` recomputed these via a `pivot_longer + group_by +
#' summarise` chain without `na.rm`, which silently diverged from the table
#' path whenever any NA slipped through the pipeline.
#'
#' @param relative_abundance Numeric matrix (features x samples). Callers
#'   are expected to have already normalized with
#'   `compute_relative_abundance()` (or any equivalent) and validated
#'   samples with `validate_abundance()` so that zero-sum columns are
#'   already ruled out.
#' @param group_vector Character (or coercible) vector of length
#'   `ncol(relative_abundance)` assigning each sample column to a group.
#'   NA entries are dropped from the summary.
#' @return A data frame with columns `name`, `group`, `mean`, `sd` in long
#'   format. `name` preserves the row order of the input matrix; `group`
#'   appears in first-seen order from `group_vector`.
#' @noRd
summarize_abundance_by_group <- function(relative_abundance, group_vector) {
  if (!is.matrix(relative_abundance)) {
    relative_abundance <- as.matrix(relative_abundance)
  }
  n_samples <- ncol(relative_abundance)
  if (length(group_vector) != n_samples) {
    stop(sprintf(
      "summarize_abundance_by_group(): group_vector length (%d) must equal ncol(relative_abundance) (%d).",
      length(group_vector), n_samples
    ), call. = FALSE)
  }

  feature_names <- rownames(relative_abundance)
  if (is.null(feature_names)) {
    feature_names <- as.character(seq_len(nrow(relative_abundance)))
  }

  group_chr <- as.character(group_vector)
  # Preserve first-seen order; drop NA.
  present <- group_chr[!is.na(group_chr)]
  unique_groups <- unique(present)

  if (length(unique_groups) == 0) {
    return(data.frame(
      name = character(0), group = character(0),
      mean = numeric(0), sd = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # Compute per-group stats first (group-major intermediate), then
  # interleave to produce feature-major output: (feature1,gA), (feature1,gB),
  # (feature2,gA), (feature2,gB), .... This matches the row order that
  # `pathway_errorbar()` previously produced via `group_by(name, group) %>%
  # summarise()` and keeps downstream `match()`/`order()` code stable
  # regardless of the `order()` sort method.
  means <- matrix(NA_real_, nrow = nrow(relative_abundance),
                  ncol = length(unique_groups),
                  dimnames = list(feature_names, unique_groups))
  sds <- means
  for (g in unique_groups) {
    mask <- !is.na(group_chr) & group_chr == g
    sub <- relative_abundance[, mask, drop = FALSE]
    # `na.rm = TRUE` matches calculate_abundance_stats() so that the two
    # entry points agree on how to treat NA abundance cells. After the
    # zero-sum guard in compute_relative_abundance(), the only way NA
    # can reach here is if the upstream abundance already contained NA.
    means[, g] <- apply(sub, 1, mean, na.rm = TRUE)
    sds[, g] <- apply(sub, 1, stats::sd, na.rm = TRUE)
  }

  n_groups <- length(unique_groups)
  data.frame(
    name = rep(feature_names, each = n_groups),
    group = rep(unique_groups, times = nrow(relative_abundance)),
    mean = as.vector(t(means)),
    sd = as.vector(t(sds)),
    stringsAsFactors = FALSE
  )
}

#' Validate Metadata
#'
#' Validates that metadata is in the correct format.
#'
#' @param metadata Data frame to validate
#' @param required_cols Character vector of required column names (optional)
#' @return TRUE if valid, otherwise stops with error
#' @noRd
validate_metadata <- function(metadata, required_cols = NULL) {
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }

  if (nrow(metadata) == 0) {
    stop("'metadata' cannot be empty")
  }

  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, colnames(metadata))
    if (length(missing_cols) > 0) {
      stop(sprintf("Required column(s) not found in metadata: %s",
                   paste(missing_cols, collapse = ", ")))
    }
  }

  TRUE
}

#' Validate Group Column
#'
#' Validates that a group column exists and has appropriate values.
#'
#' @param metadata Data frame containing group column
#' @param group Column name for grouping variable
#' @param min_groups Minimum number of groups required (default: 2)
#' @param min_per_group Minimum samples per group (default: 1)
#' @return TRUE if valid, otherwise stops with error
#' @noRd
validate_group <- function(metadata, group, min_groups = 2, min_per_group = 1) {
  if (!group %in% colnames(metadata)) {
    stop(sprintf("Group column '%s' not found in metadata. Available columns: %s",
                 group, paste(colnames(metadata), collapse = ", ")))
  }

  group_values <- metadata[[group]]
  group_levels <- unique(group_values[!is.na(group_values)])

  if (length(group_levels) < min_groups) {
    stop(sprintf("At least %d groups are required, found %d: %s",
                 min_groups, length(group_levels),
                 paste(group_levels, collapse = ", ")))
  }

  if (min_per_group > 1) {
    group_counts <- table(group_values)
    small_groups <- names(group_counts)[group_counts < min_per_group]
    if (length(small_groups) > 0) {
      stop(sprintf("Groups with fewer than %d samples: %s",
                   min_per_group, paste(small_groups, collapse = ", ")))
    }
  }

  TRUE
}

# =============================================================================
# Extended Validation Utilities
# =============================================================================

#' Validate Numeric Matrix Quality
#'
#' Deep validation of numeric matrix including negative values, NA, and duplicates.
#'
#' @param mat Numeric matrix to validate
#' @param check_negative Check for negative values (default: TRUE)
#' @param check_na Check for NA values (default: TRUE, issues warning)
#' @param check_duplicates Check for duplicate column names (default: TRUE)
#' @return TRUE if valid, otherwise stops with error
#' @noRd
validate_numeric_matrix <- function(mat, check_negative = TRUE,
                                    check_na = TRUE, check_duplicates = TRUE) {
  if (!is.numeric(mat)) {
    stop("Abundance data contains non-numeric values. All sample columns must be numeric.")
  }

  if (check_negative) {
    neg_cols <- apply(mat, 2, function(x) any(x < 0, na.rm = TRUE))
    if (any(neg_cols)) {
      stop(sprintf("Negative values found in columns: %s",
                   paste(colnames(mat)[neg_cols], collapse = ", ")))
    }
  }

  if (check_na) {
    na_count <- sum(is.na(mat))
    if (na_count > 0) {
      warning(sprintf("Found %d NA values in data", na_count))
    }
  }

  if (check_duplicates && !is.null(colnames(mat))) {
    if (anyDuplicated(colnames(mat))) {
      stop("Duplicate column names found in data")
    }
  }

  TRUE
}

#' Validate Zero Abundance Issues
#'
#' Checks for empty data, all-zero data, and excessive zero rows.
#' Returns information about zero rows for downstream filtering.
#'
#' @param mat Numeric matrix (features x samples)
#' @param warn_threshold Proportion of zero rows to trigger warning (default: 0.9)
#' @return List with: valid (logical), zero_rows (logical vector), n_nonzero (integer)
#' @noRd
validate_zero_abundance <- function(mat, warn_threshold = 0.9) {
  # Check dimensions

  if (nrow(mat) == 0) {
    stop("No features in abundance data")
  }
  if (ncol(mat) == 0) {
    stop("No samples in abundance data")
  }

  # Check all-zero
  total <- sum(mat, na.rm = TRUE)
  if (total == 0) {
    stop("All abundance values are zero")
  }

  # Identify zero rows

  zero_rows <- rowSums(mat, na.rm = TRUE) == 0
  n_nonzero <- sum(!zero_rows)
  zero_prop <- mean(zero_rows)

  if (zero_prop > warn_threshold) {
    warning(sprintf("%.0f%% of features have zero abundance across all samples",
                    zero_prop * 100))
  }

  list(valid = TRUE, zero_rows = zero_rows, n_nonzero = n_nonzero)
}

#' Validate Feature ID Format
#'
#' Checks if feature IDs match expected format patterns.
#'
#' @param ids Character vector of feature IDs
#' @param type Expected ID type: "KO", "EC", "pathway", or "auto" (default: "auto")
#' @return TRUE (issues warning if format doesn't match)
#' @noRd
validate_feature_ids <- function(ids, type = "auto") {
  patterns <- list(
    KO = "^K\\d{5}$",
    EC = "^(EC:)?\\d+\\.\\d+\\.\\d+\\.(\\d+|-)$",
    pathway = "^(ko|map|ec)?\\d{5}$"
  )

  # Auto-detect type
  if (type == "auto") {
    for (t in names(patterns)) {
      if (mean(grepl(patterns[[t]], ids)) > 0.5) {
        type <- t
        break
      }
    }
    if (type == "auto") return(TRUE)  # Unknown format, skip validation
  }

  if (type %in% names(patterns)) {
    invalid <- ids[!grepl(patterns[[type]], ids)]
    if (length(invalid) == length(ids) && length(ids) > 0) {
      # Total mismatch: the input is almost certainly the wrong data type.
      # Fail loudly here so the caller gets an actionable message instead of
      # an empty downstream matrix.
      stop(sprintf(
        "None of the feature IDs match expected %s format (e.g., got '%s'). Check that your input is %s abundance data.",
        type, paste(head(invalid, 2), collapse = ", "), type))
    } else if (length(invalid) > 0) {
      warning(sprintf("%d %s IDs don't match expected format (e.g., %s)",
                      length(invalid), type, paste(head(invalid, 2), collapse = ", ")))
    }
  }

  TRUE
}

#' Validate Input for DAA Analysis
#'
#' Comprehensive validation for differential abundance analysis input.
#' Combines zero-abundance checks with method-specific requirements.
#'
#' @param mat Numeric matrix (features x samples)
#' @param method DAA method name (optional, for method-specific checks)
#' @param min_features Minimum required non-zero features (default: 2)
#' @param filter_zero Whether to return filtered matrix (default: FALSE)
#' @return If filter_zero=FALSE: TRUE. If filter_zero=TRUE: filtered matrix.
#' @noRd
validate_daa_input <- function(mat, method = NULL, min_features = 2, filter_zero = FALSE) {
  # Matrix-quality gate: numeric, no negatives, no duplicated sample names.
  # NA values only warn (some upstream pipelines emit them legitimately) so
  # downstream methods can still fail explicitly if they cannot tolerate NAs.
  validate_numeric_matrix(mat)

  # Basic zero checks
  result <- validate_zero_abundance(mat)

  # Method-specific requirements
  if (!is.null(method)) {
    if (method == "ALDEx2" && result$n_nonzero < 2) {
      stop(sprintf("ALDEx2 requires at least 2 non-zero features (found %d)", result$n_nonzero))
    }
    if (method == "LinDA" && result$n_nonzero < 10) {
      warning(sprintf("Only %d non-zero features for LinDA; results may be unreliable", result$n_nonzero))
    }
  }

  # General minimum check
  if (result$n_nonzero < min_features) {
    stop(sprintf("At least %d non-zero features required (found %d)", min_features, result$n_nonzero))
  }

  if (filter_zero) {
    return(mat[!result$zero_rows, , drop = FALSE])
  }

  TRUE
}

#' Require Package with Consistent Error Message
#'
#' Checks if a package is available and provides consistent installation instructions.
#'
#' @param pkg Package name
#' @param purpose Optional description of why the package is needed
#' @param bioc Whether it's a Bioconductor package (default: TRUE)
#' @return TRUE if available, otherwise stops with error
#' @noRd
require_package <- function(pkg, purpose = NULL, bioc = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install_cmd <- if (bioc) {
      sprintf("BiocManager::install('%s')", pkg)
    } else {
      sprintf("install.packages('%s')", pkg)
    }
    msg <- sprintf("Package '%s' is required", pkg)
    if (!is.null(purpose)) {
      msg <- sprintf("%s for %s", msg, purpose)
    }
    stop(sprintf("%s. Install with: %s", msg, install_cmd))
  }
  TRUE
}

#' Validate Choice Parameter
#'
#' Validates that a value is one of the allowed choices.
#'
#' @param value The value to validate
#' @param choices Vector of allowed choices
#' @param param_name Name of the parameter (for error message)
#' @return TRUE if valid, otherwise stops with error
#' @noRd
validate_choice <- function(value, choices, param_name = "value") {
  if (!value %in% choices) {
    stop(sprintf("%s must be one of: %s", param_name, paste(choices, collapse = ", ")))
  }
  TRUE
}

#' Normalize a scalar logical flag
#'
#' Public APIs historically accepted both logical TRUE/FALSE and the string
#' forms "TRUE"/"FALSE" for some flags. Normalize that contract once at the
#' boundary so downstream code can use plain logical tests without drifting
#' between `isTRUE(x)`, `x == TRUE`, and `if (x)` semantics.
#' @noRd
normalize_logical_flag <- function(value, param_name) {
  if (is.logical(value) && length(value) == 1 && !is.na(value)) {
    return(value)
  }

  if (is.character(value) && length(value) == 1 && !is.na(value)) {
    lowered <- tolower(trimws(value))
    if (lowered %in% c("true", "t", "1")) return(TRUE)
    if (lowered %in% c("false", "f", "0")) return(FALSE)
  }

  stop(sprintf("'%s' must be TRUE or FALSE.", param_name), call. = FALSE)
}

#' Validate Data Frame Input
#'
#' Validates that input is a data frame with required columns.
#'
#' @param df The data frame to validate
#' @param required_cols Character vector of required column names (optional)
#' @param param_name Name of the parameter (for error message)
#' @return TRUE if valid, otherwise stops with error
#' @noRd
validate_dataframe <- function(df, required_cols = NULL, param_name = "data") {
  if (!is.data.frame(df)) {
    stop(sprintf("'%s' must be a data frame", param_name))
  }
  if (!is.null(required_cols)) {
    missing <- setdiff(required_cols, colnames(df))
    if (length(missing) > 0) {
      stop(sprintf("Missing required columns in %s: %s", param_name, paste(missing, collapse = ", ")))
    }
  }
  TRUE
}

# =============================================================================
# Empty Result Frame Factories
# =============================================================================

#' Create Empty DAA Result Data Frame
#'
#' Creates an empty data frame with the correct structure for DAA results.
#'
#' @param method Character string specifying the DAA method used
#' @param include_log2fc Logical. If TRUE, includes log2_fold_change column
#'   instead of p_adjust and adj_method (for LinDA-style results)
#' @return Empty data frame with standard DAA result columns
#' @noRd
create_empty_daa_result <- function(method = "unknown", include_log2fc = FALSE) {
  if (include_log2fc) {
    data.frame(
      feature = character(0),
      method = character(0),
      group1 = character(0),
      group2 = character(0),
      p_values = numeric(0),
      log2_fold_change = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      feature = character(0),
      method = character(0),
      group1 = character(0),
      group2 = character(0),
      p_values = numeric(0),
      p_adjust = numeric(0),
      adj_method = character(0),
      stringsAsFactors = FALSE
    )
  }
}

#' Create Empty GSEA Result Data Frame
#'
#' Creates an empty data frame with the correct structure for GSEA results.
#'
#' @param method Character string specifying the GSEA method used
#' @param full Logical. If TRUE, includes ES, NES, leading_edge columns
#' @return Empty data frame with standard GSEA result columns
#' @noRd
create_empty_gsea_result <- function(method = "unknown", full = FALSE) {
  if (full) {
    data.frame(
      pathway_id = character(0),
      pathway_name = character(0),
      size = integer(0),
      ES = numeric(0),
      NES = numeric(0),
      pvalue = numeric(0),
      p.adjust = numeric(0),
      leading_edge = character(0),
      method = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      pathway_id = character(0),
      pathway_name = character(0),
      size = integer(0),
      direction = character(0),
      pvalue = numeric(0),
      p.adjust = numeric(0),
      method = character(0),
      NES = numeric(0),
      leading_edge = character(0),
      stringsAsFactors = FALSE
    )
  }
}

#' Validate DAA results data frame
#'
#' Validates that a DAA results data frame meets requirements for visualization
#' functions (single method, single group pair).
#'
#' @param daa_results_df Data frame containing DAA results
#' @param require_single_method Logical. If TRUE, requires exactly one method
#' @param require_single_group_pair Logical. If TRUE, requires exactly one group1/group2 pair
#' @return Invisible TRUE if validation passes
#' @keywords internal
validate_daa_results <- function(daa_results_df,
                                  require_single_method = TRUE,
                                  require_single_group_pair = TRUE) {
  if (require_single_method && length(unique(daa_results_df$method)) != 1) {
    stop("daa_results_df contains multiple methods. Filter to one method first.")
  }
  if (require_single_group_pair) {
    if (length(unique(daa_results_df$group1)) != 1 ||
        length(unique(daa_results_df$group2)) != 1) {
      stop("daa_results_df contains multiple group pairs. Filter to one pair first.")
    }
  }
  invisible(TRUE)
}

#' Require a column exists in a data frame
#'
#' Simple check that a single column exists. For checking multiple columns,
#' use validate_dataframe() with required_cols parameter.
#'
#' @param df Data frame to check
#' @param col Column name to require
#' @param param_name Name of the data frame parameter (for error message)
#' @return Invisible TRUE if column exists, otherwise stops with error
#' @keywords internal
require_column <- function(df, col, param_name = deparse(substitute(df))) {
  if (!col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in %s", col, param_name))
  }
  invisible(TRUE)
}
