#' Data Utilities for ggpicrust2
#'
#' Internal utility functions for data preprocessing, sample matching,
#' and format standardization.
#'
#' @name data_utils
#' @keywords internal
NULL

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
samples_match <- function(vec1, vec2, threshold = 0.5) {
  vec1 <- as.character(vec1)
  vec2 <- as.character(vec2)
  n_common <- length(intersect(vec1, vec2))
  min_length <- min(length(vec1), length(vec2))

  if (min_length == 0) return(FALSE)
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

  for (col in standard_names) {
    if (col %in% colnames(metadata)) {
      if (samples_match(metadata[[col]], abundance_samples)) {
        return(col)
      }
    }
  }

  # Priority 2: Any column that matches
  for (col in colnames(metadata)) {
    if (samples_match(metadata[[col]], abundance_samples)) {
      return(col)
    }
  }

  # Priority 3: Rownames
  if (samples_match(rownames(metadata), abundance_samples)) {
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
  metadata_aligned <- metadata[row_idx, , drop = FALSE]

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
