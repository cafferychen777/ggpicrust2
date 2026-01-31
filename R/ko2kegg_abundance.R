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

  # Load data from file or use provided data
  if (!is.null(file)) {
    abundance <- read_abundance_file(file)
  } else {
    if (!is.data.frame(data)) {
      stop("'data' must be a data.frame")
    }
    abundance <- data
    if (ncol(abundance) < 2) {
      stop("Data must have at least 2 columns (KO IDs and samples)")
    }
  }

  # Standardize first column name (PICRUSt2 uses various formats)
  valid_ko_column_names <- c("function.", "function", "sequence", "#NAME")
  first_col_name <- colnames(abundance)[1]

  if (first_col_name %in% valid_ko_column_names) {
    colnames(abundance)[1] <- "function."
  } else {
    # Check if first column contains KO IDs
    first_col_values <- as.character(abundance[[1]])
    ko_pattern <- "(^K\\d{5}$)|(^ko:K\\d{5}$)"
    if (any(grepl(ko_pattern, first_col_values))) {
      colnames(abundance)[1] <- "function."
    } else {
      stop(sprintf("First column '%s' does not contain KO IDs (expected: K##### or ko:K#####)", first_col_name))
    }
  }

  if (nrow(abundance) == 0) {
    stop("Data contains no rows")
  }

  # Standardize KO ID format and validate
 abundance <- clean_ko_abundance(abundance, verbose = FALSE)

  # Validate numeric matrix (negative values, NA, duplicates)
  numeric_mat <- as.matrix(abundance[, -1, drop = FALSE])
  validate_numeric_matrix(numeric_mat)
  validate_feature_ids(abundance[[1]], type = "KO")
  
  # Load KEGG reference data using unified loader
  ko_to_kegg_reference <- load_reference_data("ko_to_kegg")

  # Filter for prokaryotic pathways if requested
  if (filter_for_prokaryotes) {
    # Exclude eukaryote-specific KEGG Level 2 categories
    eukaryote_specific_level2 <- c(
      "9161 Cancer: overview", "9162 Cancer: specific types",
      "9163 Immune disease", "9164 Neurodegenerative disease",
      "9165 Substance dependence", "9166 Cardiovascular disease",
      "9167 Endocrine and metabolic disease",
      "9149 Aging", "9151 Immune system", "9152 Endocrine system",
      "9153 Circulatory system", "9154 Digestive system",
      "9155 Excretory system", "9156 Nervous system",
      "9157 Sensory system", "9158 Development and regeneration",
      "9159 Environmental adaptation"
    )
    ko_to_kegg_reference <- ko_to_kegg_reference[
      !ko_to_kegg_reference$level2 %in% eukaryote_specific_level2,
    ]
  }

  # Get all unique pathway IDs
  all_pathways <- unique(ko_to_kegg_reference$pathway_id)

  # Create fast pathway → KOs lookup mapping
  pathway_to_ko <- split(ko_to_kegg_reference$ko_id, ko_to_kegg_reference$pathway_id)

  # Create result data frame
  kegg_abundance <- data.frame(
    row.names = all_pathways  # Set pathways as row names
  )

  sample_names <- colnames(abundance)[-1]
  kegg_abundance[sample_names] <- 0

  # Calculate pathway abundances with progress bar
  pb <- txtProgressBar(min = 0, max = nrow(kegg_abundance), style = 3)
  on.exit(close(pb), add = TRUE)

  for (i in seq_len(nrow(kegg_abundance))) {
    setTxtProgressBar(pb, i)
    current_kegg <- rownames(kegg_abundance)[i]
    relevant_kos <- pathway_to_ko[[current_kegg]]
    matching_rows <- abundance[[1]] %in% relevant_kos

    if (any(matching_rows)) {
      if (method == "sum") {
        kegg_abundance[i, ] <- colSums(abundance[matching_rows, -1, drop = FALSE])
      } else {
        # PICRUSt2-style: upper-half mean
        ko_abundances <- as.matrix(abundance[matching_rows, -1, drop = FALSE])
        kegg_abundance[i, ] <- apply(ko_abundances, 2, function(col) {
          sorted_col <- sort(col)
          mean(sorted_col[ceiling(length(sorted_col) / 2):length(sorted_col)])
        })
      }
    }
  }
  close(pb)

  # Remove zero-abundance pathways
  zero_rows <- rowSums(kegg_abundance, na.rm = TRUE) == 0
  if (all(zero_rows)) {
    warning("No KO IDs matched known pathways")
    kegg_abundance <- kegg_abundance[0, , drop = FALSE]
  } else {
    kegg_abundance <- kegg_abundance[!zero_rows, , drop = FALSE]
  }

  kegg_abundance
}
