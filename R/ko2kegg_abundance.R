#' Convert KO abundance in picrust2 export files to KEGG pathway abundance
#'
#' This function takes a file containing KO (KEGG Orthology) abundance data in picrust2 export format and converts it to KEGG pathway abundance data.
#' The input file should be in .tsv, .txt, or .csv format.
#'
#' @param file A character string representing the file path of the input file containing KO abundance data in picrust2 export format. The input file should have KO identifiers in the first column and sample identifiers in the first row. The remaining cells should contain the abundance values for each KO-sample pair.
#' @param data An optional data.frame containing KO abundance data in the same format as the input file. If provided, the function will use this data instead of reading from the file. By default, this parameter is set to NULL.
#' @param method Method for calculating pathway abundance. One of:
#'   \itemize{
#'     \item \code{"abundance"}: (Default) PICRUSt2-style calculation using the mean of upper-half sorted KO abundances. This method is more robust and avoids inflating abundances for pathways with more KOs.
#'     \item \code{"sum"}: Simple summation of all KO abundances. This is the legacy method and may double-count KOs belonging to multiple pathways.
#'   }
#' @param filter_for_prokaryotes Logical. If TRUE (default), filters out KEGG pathways
#'   that are not relevant to prokaryotic (bacterial/archaeal) analysis. This removes
#'   pathways in categories such as:
#'   \itemize{
#'     \item Human diseases (cancer, neurodegenerative diseases, addiction, etc.)
#'     \item Organismal systems (immune system, nervous system, endocrine system, etc.)
#'   }
#'   Bacterial infection pathways and antimicrobial resistance pathways are retained.
#'   Set to FALSE to include all KEGG pathways (for eukaryotic analysis or custom filtering).
#'
#' @return
#' A data frame with KEGG pathway abundance values. Rows represent KEGG pathways, identified by their KEGG pathway IDs. Columns represent samples, identified by their sample IDs from the input file.
#'
#' @details
#' The default \code{"abundance"} method follows PICRUSt2's approach for calculating pathway abundance:
#' \enumerate{
#'   \item For each pathway, collect abundances of all associated KOs present in the data
#'   \item Sort the abundances in ascending order
#'   \item Take the upper half of the sorted values
#'   \item Calculate the mean as the pathway abundance
#' }
#'
#' This approach has several advantages over simple summation:
#' \itemize{
#'   \item Does not inflate abundances for pathways containing more KOs
#'   \item More robust to missing or low-abundance KOs
#'   \item Provides a more accurate representation of pathway activity
#' }
#'
#' The \code{"sum"} method is provided for backward compatibility and simply sums all KO abundances for each pathway.
#'
#' @section Pathway Filtering:
#' When \code{filter_for_prokaryotes = TRUE}, the function excludes KEGG pathways that are
#' biologically irrelevant to prokaryotic organisms. KEGG reference pathways include pathways
#' from all domains of life, and many human/animal-specific pathways would appear in bacterial
#' analysis simply because some KOs are shared across organisms.
#'
#' The following KEGG Level 2 categories are excluded:
#' \itemize{
#'   \item Cancer pathways (overview and specific types)
#'   \item Neurodegenerative diseases (Alzheimer's, Parkinson's, etc.)
#'   \item Substance dependence (addiction pathways)
#'   \item Cardiovascular diseases
#'   \item Endocrine and metabolic diseases
#'   \item Immune diseases
#'   \item Organismal systems (immune, nervous, endocrine, digestive, etc.)
#' }
#'
#' The following are RETAINED even with filtering:
#' \itemize{
#'   \item Infectious disease: bacterial (Salmonella, E. coli, Tuberculosis, etc.)
#'   \item Drug resistance: antimicrobial (antibiotic resistance)
#'   \item All Metabolism pathways
#'   \item Genetic/Environmental Information Processing
#'   \item Cellular Processes
#' }
#'
#' @examples
#' \dontrun{
#' library(ggpicrust2)
#' library(readr)
#'
#' # Example 1: Default - filtered for prokaryotic analysis
#' data(ko_abundance)
#' kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
#'
#' # Example 2: Include all pathways (for eukaryotic analysis)
#' kegg_abundance_all <- ko2kegg_abundance(data = ko_abundance, filter_for_prokaryotes = FALSE)
#'
#' # Example 3: Using legacy sum method with filtering
#' kegg_abundance_sum <- ko2kegg_abundance(data = ko_abundance, method = "sum")
#'
#' # Example 4: From file
#' input_file <- "path/to/your/picrust2/results/pred_metagenome_unstrat.tsv"
#' kegg_abundance <- ko2kegg_abundance(file = input_file)
#' }
#' @export
ko2kegg_abundance <- function (file = NULL, data = NULL, method = c("abundance", "sum"),
                               filter_for_prokaryotes = TRUE) {
  method <- match.arg(method)
  # Basic parameter validation
  if (is.null(file) & is.null(data)) {
    stop("Error: Please provide either a file or a data.frame.")
  }

  if (!is.null(file) && !is.null(data)) {
    warning("Both file and data provided. Using data and ignoring file.")
  }

  # File format validation function
  validate_file_format <- function(file_path) {
    # Check if file_path is actually a character string (file path)
    if (!is.character(file_path) || length(file_path) != 1) {
      stop("Error: 'file' parameter must be a character string representing a file path, not a data frame. If you have already loaded your data, please use the 'data' parameter instead.")
    }

    valid_extensions <- c(".txt", ".tsv", ".csv")
    ext <- tolower(tools::file_ext(file_path))

    # Handle case where file_ext returns empty string or multiple values
    if (length(ext) == 0 || any(ext == "")) {
      stop("Error: Unable to determine file extension. Please ensure the file has a valid extension (.txt, .tsv, or .csv).")
    }

    # Use any() to ensure we get a single logical value
    if (!any(paste0(".", ext) %in% valid_extensions)) {
      stop(sprintf("Error: Input file should be in %s format.",
                  paste(valid_extensions, collapse = ", ")))
    }
    return(paste0(".", ext[1]))  # Return first extension if multiple
  }
  
  # Data frame validation function
  validate_abundance_data <- function(df) {
    if (ncol(df) < 2) {
      stop("Error: Data must contain at least one KO column and one sample column.")
    }
    if (!is.character(df[[1]]) && !is.factor(df[[1]])) {
      stop("Error: First column must contain KO identifiers.")
    }
    if (!all(sapply(df[-1], is.numeric))) {
      stop("Error: All sample columns must contain numeric values.")
    }
  }
  
  # Input validation function
  validate_input <- function(data) {
    # Basic structure validation
    if (!is.data.frame(data)) {
      stop("The provided data must be a data.frame")
    }

    if (ncol(data) < 2) {
      stop("Data must contain at least one KO column and one sample column")
    }

    # Check first column (KO IDs)
    ko_col <- data[[1]]

    # Check KO ID format
    ko_pattern <- "^K\\d{5}$"
    invalid_kos <- ko_col[!grepl(ko_pattern, ko_col)]
    if (length(invalid_kos) > 0) {
      warning(sprintf(
        "Found %d KO IDs that don't match the expected format (K#####).\nFirst few invalid IDs: %s",
        length(invalid_kos),
        paste(head(invalid_kos, 3), collapse = ", ")
      ))
    }

    # Check numeric columns
    numeric_cols <- data[,-1, drop = FALSE]

    # Check if all sample columns are numeric
    non_numeric_cols <- names(numeric_cols)[!sapply(numeric_cols, is.numeric)]
    if (length(non_numeric_cols) > 0) {
      stop(sprintf(
        "The following columns contain non-numeric values: %s",
        paste(non_numeric_cols, collapse = ", ")
      ))
    }

    # Check for negative values
    neg_values <- sapply(numeric_cols, function(x) any(x < 0, na.rm = TRUE))
    if (any(neg_values)) {
      stop(sprintf(
        "Negative abundance values found in columns: %s",
        paste(names(numeric_cols)[neg_values], collapse = ", ")
      ))
    }

    # Check for missing values
    na_counts <- sapply(numeric_cols, function(x) sum(is.na(x)))
    cols_with_na <- names(na_counts[na_counts > 0])
    if (length(cols_with_na) > 0) {
      warning(sprintf(
        "Missing values found in columns: %s\nTotal NA count per column: %s",
        paste(cols_with_na, collapse = ", "),
        paste(paste(cols_with_na, na_counts[cols_with_na], sep = ": "), collapse = ", ")
      ))
    }

    # Check for all-zero rows
    zero_rows <- rowSums(numeric_cols == 0) == ncol(numeric_cols)
    if (any(zero_rows)) {
      warning(sprintf(
        "%d KOs have zero abundance across all samples",
        sum(zero_rows)
      ))
    }

    # Check for duplicate column names
    if (any(duplicated(colnames(data)))) {
      stop("Duplicate column names found in the input data")
    }

    return(TRUE)
  }

  # File format validation
  if (!is.null(file)) {
    file_format <- validate_file_format(file)
    message("Loading data from file...")
    abundance <- switch(
      file_format,
      ".txt" = readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
      ".tsv" = readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
      ".csv" = readr::read_delim(file, delim = ",", escape_double = FALSE, trim_ws = TRUE)
    )
  } else {
    # Validate that data is a data.frame before proceeding
    if (!is.data.frame(data)) {
      stop("The provided data must be a data.frame")
    }
    message("Processing provided data frame...")
    abundance <- data
  }

  # Validate data format
  if (ncol(abundance) < 2) {
    stop("Error: Data must have at least 2 columns (KO IDs and at least one sample)")
  }

  # PICRUSt2 compatibility: Accept various column name formats for KO IDs
  # Common formats: "function.", "function", "sequence", "#NAME"
  valid_ko_column_names <- c("function.", "function", "sequence", "#NAME")
  first_col_name <- colnames(abundance)[1]

  if (first_col_name %in% valid_ko_column_names) {
    # Standardize to "function." for internal processing
    if (first_col_name != "function.") {
      message(sprintf("Detected column name '%s', renaming to 'function.' for processing...", first_col_name))
      colnames(abundance)[1] <- "function."
    }
  } else {
    # Check if first column contains KO IDs (K##### format)
    first_col_values <- as.character(abundance[[1]])
    ko_pattern <- "(^K\\d{5}$)|(^ko:K\\d{5}$)"
    if (any(grepl(ko_pattern, first_col_values))) {
      message(sprintf("Column '%s' appears to contain KO IDs, renaming to 'function.' for processing...", first_col_name))
      colnames(abundance)[1] <- "function."
    } else {
      stop(sprintf("Error: First column '%s' does not appear to contain KO IDs. Expected formats: K##### or ko:K#####", first_col_name))
    }
  }

  if (nrow(abundance) == 0) {
    stop("Error: Data contains no rows")
  }

  # PICRUSt 2.6.2 compatibility: Clean KO ID format
  message("Checking KO ID format compatibility...")
  ko_ids <- abundance[[1]]

  # Detect and clean "ko:" prefix
  has_ko_prefix <- any(grepl("^ko:", ko_ids))
  if (has_ko_prefix) {
    message("Detected PICRUSt 2.6.2 format with 'ko:' prefix. Applying compatibility fix...")
    abundance[[1]] <- gsub("^ko:", "", abundance[[1]])
    message(sprintf("Cleaned %d KO IDs by removing 'ko:' prefix", sum(grepl("^ko:", ko_ids))))
  }

  # Run input validation
  tryCatch({
    validate_input(abundance)
  }, warning = function(w) {
    warning("Warning during data validation: ", w$message, call. = FALSE)
  }, error = function(e) {
    stop("Data validation failed: ", e$message, call. = FALSE)
  })
  
  message("Using KEGG reference data from package...")
  # Use internal data from sysdata.rda (automatically loaded in package namespace)
  # Objects in sysdata.rda are directly accessible within package functions
  # ko_to_kegg_reference is already defined in the package namespace from sysdata.rda
  # We just convert it to data.frame to ensure compatibility
  ko_to_kegg_reference <- as.data.frame(ko_to_kegg_reference)

  # Filter for prokaryotic pathways if requested
  if (filter_for_prokaryotes) {
    message("Filtering pathways for prokaryotic analysis...")

    # Define KEGG Level 2 categories to EXCLUDE (eukaryote-specific)
    eukaryote_specific_level2 <- c(
      # Human Diseases - exclude non-infectious diseases
      "9161 Cancer: overview",
      "9162 Cancer: specific types",
      "9163 Immune disease",
      "9164 Neurodegenerative disease",
      "9165 Substance dependence",
      "9166 Cardiovascular disease",
      "9167 Endocrine and metabolic disease",
      # Note: We KEEP "9171 Infectious disease: bacterial" and
      #       "9175 Drug resistance: antimicrobial" as they are prokaryote-relevant

      # Organismal Systems - exclude all (eukaryote-specific)
      "9149 Aging",
      "9151 Immune system",
      "9152 Endocrine system",
      "9153 Circulatory system",
      "9154 Digestive system",
      "9155 Excretory system",
      "9156 Nervous system",
      "9157 Sensory system",
      "9158 Development and regeneration",
      "9159 Environmental adaptation"
    )

    # Filter out eukaryote-specific pathways
    original_count <- length(unique(ko_to_kegg_reference$pathway_id))
    ko_to_kegg_reference <- ko_to_kegg_reference[
      !ko_to_kegg_reference$level2 %in% eukaryote_specific_level2,
    ]
    filtered_count <- length(unique(ko_to_kegg_reference$pathway_id))

    message(sprintf("  Removed %d eukaryote-specific pathways (%d -> %d pathways)",
                    original_count - filtered_count, original_count, filtered_count))
  }

  # Get all unique pathway IDs
  all_pathways <- unique(ko_to_kegg_reference$pathway_id)

  message(sprintf("Processing %d KEGG pathways using '%s' method...", length(all_pathways), method))

  # Create fast pathway → KOs lookup mapping
  message("Building pathway-KO index...")
  pathway_to_ko <- split(ko_to_kegg_reference$ko_id, ko_to_kegg_reference$pathway_id)

  # Create result data frame
  kegg_abundance <- data.frame(
    row.names = all_pathways  # Set pathways as row names
  )

  # Add columns for each sample
  sample_names <- colnames(abundance)[-1]
  for (sample in sample_names) {
    kegg_abundance[[sample]] <- 0
  }

  # Use tryCatch for error handling
  tryCatch({
    pb <- txtProgressBar(min = 0, max = nrow(kegg_abundance), style = 3)
    start_time <- Sys.time()

    total_matches <- 0

    for (i in seq_len(nrow(kegg_abundance))) {
      setTxtProgressBar(pb, i)
      current_kegg <- rownames(kegg_abundance)[i]

      # Get all KOs for current pathway from index
      relevant_kos <- pathway_to_ko[[current_kegg]]

      # Find matching KOs
      matching_rows <- abundance[[1]] %in% relevant_kos
      if (any(matching_rows)) {
        total_matches <- total_matches + sum(matching_rows)

        if (method == "sum") {
          # Legacy method: simple summation
          kegg_abundance[i, ] <- colSums(abundance[matching_rows, -1, drop = FALSE])
        } else {
          # PICRUSt2-style method: upper-half mean
          # For each sample, calculate the mean of the upper half of sorted KO abundances
          ko_abundances <- as.matrix(abundance[matching_rows, -1, drop = FALSE])

          for (j in seq_along(sample_names)) {
            sample_abun <- ko_abundances[, j]
            # Sort abundances in ascending order
            sorted_abun <- sort(sample_abun)
            # Calculate the starting index for upper half (ceiling handles odd numbers)
            half_i <- ceiling(length(sorted_abun) / 2)
            # Take upper half and calculate mean
            upper_half <- sorted_abun[half_i:length(sorted_abun)]
            kegg_abundance[i, j] <- mean(upper_half)
          }
        }
      }
    }
    
    end_time <- Sys.time()
    close(pb)

    # Remove zero-abundance pathways
    message("Removing KEGG pathways with zero abundance across all samples...")
    # Note: kegg_abundance only has sample columns, pathway IDs are in rownames
    zero_rows <- rowSums(kegg_abundance, na.rm = TRUE) == 0

    # Check if all pathways have zero abundance
    if (all(zero_rows)) {
      # All pathways are zero - this usually means no KO IDs matched any pathway
      # Return empty data frame
      message("All KEGG pathways have zero abundance. No KO IDs matched known pathways.")

      # Create empty but correctly structured data frame
      sample_names <- colnames(abundance)[-1]
      kegg_abundance <- data.frame(matrix(ncol = length(sample_names), nrow = 0))
      colnames(kegg_abundance) <- sample_names
    } else {
      # Normal filtering logic: remove zero rows
      kegg_abundance <- kegg_abundance[!zero_rows, , drop = FALSE]
    }

    message(sprintf("Final number of KEGG pathways: %d", nrow(kegg_abundance)))

    # Empty result is valid - it means no KO IDs matched any pathways
    # Simply return the empty data frame with correct structure
    
    return(kegg_abundance)
    
  }, error = function(e) {
    if (exists("pb")) close(pb)
    message("Error occurred during processing: ", e$message)
    stop(e)
  })
}
