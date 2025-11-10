#' Convert KO abundance in picrust2 export files to KEGG pathway abundance
#'
#' This function takes a file containing KO (KEGG Orthology) abundance data in picrust2 export format and converts it to KEGG pathway abundance data.
#' The input file should be in .tsv, .txt, or .csv format.
#'
#' @param file A character string representing the file path of the input file containing KO abundance data in picrust2 export format. The input file should have KO identifiers in the first column and sample identifiers in the first row. The remaining cells should contain the abundance values for each KO-sample pair.
#' @param data An optional data.frame containing KO abundance data in the same format as the input file. If provided, the function will use this data instead of reading from the file. By default, this parameter is set to NULL.
#'
#' @return
#' A data frame with KEGG pathway abundance values. Rows represent KEGG pathways, identified by their KEGG pathway IDs. Columns represent samples, identified by their sample IDs from the input file. Each cell contains the abundance of a specific KEGG pathway in a given sample, calculated by summing the abundances of the corresponding KOs in the input file.
#' @examples
#' \dontrun{
#' library(ggpicrust2)
#' library(readr)
#'
#' # Example 1: Demonstration with a hypothetical input file
#'
#' # Prepare an input file path
#' input_file <- "path/to/your/picrust2/results/pred_metagenome_unstrat.tsv"
#'
#' # Run ko2kegg_abundance function
#' kegg_abundance <- ko2kegg_abundance(file = input_file)
#'
#' # Alternatively, read the data from a file and use the data argument
#' file_path <- "path/to/your/picrust2/results/pred_metagenome_unstrat.tsv"
#' ko_abundance <- read_delim(file_path, delim = "\t")
#' kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
#'
#' # Example 2: Working with real data
#' # In this case, we're using an existing dataset from the ggpicrust2 package.
#'
#' # Load the data
#' data(ko_abundance)
#'
#' # Apply the ko2kegg_abundance function to our real dataset
#' kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
#' }
#' @export
ko2kegg_abundance <- function (file = NULL, data = NULL) {
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
    message("Processing provided data frame...")
    abundance <- data
  }

  # Validate data format
  if (ncol(abundance) < 2) {
    stop("Error: Data must have at least 2 columns (KO IDs and at least one sample)")
  }

  if (!("function." %in% colnames(abundance))) {
    if (colnames(abundance)[1] != "function.") {
      stop("Error: First column must be named 'function.' containing KO IDs")
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

  # Get all unique pathway IDs
  all_pathways <- unique(ko_to_kegg_reference$pathway_id)

  message(sprintf("Processing %d KEGG pathways...", length(all_pathways)))

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
        kegg_abundance[i, ] <- colSums(abundance[matching_rows, -1, drop = FALSE])
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
