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
utils::globalVariables(c("read.csv"))
import_MicrobiomeAnalyst_daa_results <- function(file_path = NULL, data = NULL, method = "MicrobiomeAnalyst", group_levels = NULL) {
  # Check if a data frame is provided
  if (is.null(data)) {
    # Read the CSV file if data frame is not provided
    if (is.null(file_path)) {
      stop("Please provide either a file_path or a data frame.")
    }
    data <- read.csv(file_path, stringsAsFactors = FALSE)
  }

  if (is.null(group_levels)) {
    stop("Please provide group levels.")
  }

  # Rename the columns
  names(data) <- c("feature", "p_values", "p_adjust", "Statistics")

  # Create a new column for method
  data$method <- method

  # Create new columns for group levels
  for (i in seq_along(group_levels)) {
    data[paste0("group", i)] <- group_levels[i]
  }

  return(data)
}
