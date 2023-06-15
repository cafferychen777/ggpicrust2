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
  if (is.null(file) & is.null(data)) {
    stop("Error: Please provide either a file or a data.frame.")
  }

  if (!is.null(file)) {
    message("Loading data from file...")
    file_format <- substr(file, nchar(file) - 3, nchar(file))
    switch(
      file_format,
      ".txt" = abundance <- readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
      ".tsv" = abundance <- readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
      ".csv" = abundance <- readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
      stop("Error: Input file should be in .tsv, .txt or .csv format.\nThe optimal input file is the output file 'pred_metagenome_unstrat.tsv' from PICRUSt2.")
    )
  } else if (!is.null(data)) {
    if (is.data.frame(data)) {
      message("Processing provided data frame...")
      abundance <- data
    } else {
      stop("Error: The provided data must be a data.frame.")
    }
  }

  message("Loading KEGG reference data. This might take a while...")
  load(system.file("extdata", "kegg_reference.RData", package = "ggpicrust2"))

  message("Performing KO to KEGG conversion. Please be patient, this might take a while...")
  sample_names <- colnames(abundance)[-1]
  kegg_names <- ko_to_kegg_reference[, 1]

  kegg_abundance <- matrix(NA, nrow = nrow(kegg_names), ncol = length(sample_names))
  colnames(kegg_abundance) <- sample_names
  rownames(kegg_abundance) <- as.matrix(kegg_names)

  # 初始化一个文本进度条
  pb <- txtProgressBar(min = 0, max = nrow(kegg_abundance), style = 3)

  start_time <- Sys.time()  # 记录开始时间
  for (i in seq_len(nrow(kegg_abundance))) {
    for (j in seq_len(ncol(kegg_abundance))) {
      kegg_name <- rownames(kegg_abundance)[i]
      sample_name <- colnames(kegg_abundance)[j]
      ko_to_kegg <- ko_to_kegg_reference[ko_to_kegg_reference[,1] == kegg_name,-1]
      ko_to_kegg <- ko_to_kegg[!is.na(ko_to_kegg)]
      kegg_abundance[i, j] <- sum(abundance[as.matrix(abundance[,1]) %in% ko_to_kegg, sample_name])
    }
    # 更新进度条
    setTxtProgressBar(pb, i)
  }
  end_time <- Sys.time()  # 记录结束时间

  close(pb)  # 关闭进度条

  message(paste0("KO to KEGG conversion completed. Time elapsed: ", round(end_time - start_time, 2), " seconds."))

  message("Removing KEGG pathways with zero abundance across all samples...")
  kegg_abundance <- kegg_abundance[rowSums(kegg_abundance) != 0, ]
  kegg_abundance <- as.data.frame(kegg_abundance)

  message("KEGG abundance calculation completed successfully.")
  return(kegg_abundance)
}
