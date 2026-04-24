#' Import Differential Abundance Analysis (DAA) results from MicrobiomeAnalyst
#'
#' This function imports DAA results from an external platform such as MicrobiomeAnalyst. It can be used to compare the results obtained from different platforms.
#' @name import_MicrobiomeAnalyst_daa_results
#' @param file_path a character string specifying the path to the CSV file containing the DAA results from MicrobiomeAnalyst. If this parameter is NULL and no data frame is provided, an error will be thrown. Default is NULL.
#' @param data a data frame containing the DAA results from MicrobiomeAnalyst. If this parameter is NULL and no file path is provided, an error will be thrown. Default is NULL.
#' @param method a character string specifying the method used for the DAA. This will be added as a new column in the returned data frame. Default is "MicrobiomeAnalyst".
#' @param group_levels a character vector specifying the group levels for the DAA. This will be added as new columns in the returned data frame. Default is c("control", "treatment").
#'
#' @return a data frame containing the DAA results from MicrobiomeAnalyst with additional columns for the method and group levels.
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
    data <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  }

  validate_dataframe(data, param_name = "data")

  if (!is.character(group_levels) || length(group_levels) < 1 || anyNA(group_levels)) {
    stop("group_levels must be a non-empty character vector without NA values.",
         call. = FALSE)
  }

  required_output_cols <- c("feature", "p_values", "p_adjust", "Statistics")
  if (ncol(data) < length(required_output_cols)) {
    stop("MicrobiomeAnalyst DAA results must contain at least four columns: ",
         paste(required_output_cols, collapse = ", "),
         call. = FALSE)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  names(data)[seq_along(required_output_cols)] <- required_output_cols

  # Create a new column for method
  data$method <- method

  # Create new columns for group levels
  for (i in seq_along(group_levels)) {
    data[paste0("group", i)] <- group_levels[i]
  }

  return(data)
}
