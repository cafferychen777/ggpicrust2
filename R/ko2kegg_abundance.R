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
  # 基本参数检查
  if (is.null(file) & is.null(data)) {
    stop("Error: Please provide either a file or a data.frame.")
  }
  
  if (!is.null(file) && !is.null(data)) {
    warning("Both file and data provided. Using data and ignoring file.")
  }
  
  # 文件格式验证函数
  validate_file_format <- function(file_path) {
    valid_extensions <- c(".txt", ".tsv", ".csv")
    ext <- tolower(tools::file_ext(file_path))
    if (!paste0(".", ext) %in% valid_extensions) {
      stop(sprintf("Error: Input file should be in %s format.",
                  paste(valid_extensions, collapse = ", ")))
    }
  }
  
  # 数据框验证函数
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
  
  # 添加新的输入验证函数
  validate_input <- function(data) {
    # 基本结构验证
    if (!is.data.frame(data)) {
      stop("The provided data must be a data.frame")
    }
    
    if (ncol(data) < 2) {
      stop("Data must contain at least one KO column and one sample column")
    }
    
    # 检查第一列（KO IDs）
    ko_col <- data[[1]]
    
    # 检查 KO ID 格式
    ko_pattern <- "^K\\d{5}$"
    invalid_kos <- ko_col[!grepl(ko_pattern, ko_col)]
    if (length(invalid_kos) > 0) {
      warning(sprintf(
        "Found %d KO IDs that don't match the expected format (K#####).\nFirst few invalid IDs: %s",
        length(invalid_kos),
        paste(head(invalid_kos, 3), collapse = ", ")
      ))
    }
    
    # 检查数值列
    numeric_cols <- data[,-1, drop = FALSE]
    
    # 检查是否所有样本列都是数值型
    non_numeric_cols <- names(numeric_cols)[!sapply(numeric_cols, is.numeric)]
    if (length(non_numeric_cols) > 0) {
      stop(sprintf(
        "The following columns contain non-numeric values: %s",
        paste(non_numeric_cols, collapse = ", ")
      ))
    }
    
    # 检查负值
    neg_values <- sapply(numeric_cols, function(x) any(x < 0, na.rm = TRUE))
    if (any(neg_values)) {
      stop(sprintf(
        "Negative abundance values found in columns: %s",
        paste(names(numeric_cols)[neg_values], collapse = ", ")
      ))
    }
    
    # 检查缺失值
    na_counts <- sapply(numeric_cols, function(x) sum(is.na(x)))
    cols_with_na <- names(na_counts[na_counts > 0])
    if (length(cols_with_na) > 0) {
      warning(sprintf(
        "Missing values found in columns: %s\nTotal NA count per column: %s",
        paste(cols_with_na, collapse = ", "),
        paste(paste(cols_with_na, na_counts[cols_with_na], sep = ": "), collapse = ", ")
      ))
    }
    
    # 检查全零行
    zero_rows <- rowSums(numeric_cols == 0) == ncol(numeric_cols)
    if (any(zero_rows)) {
      warning(sprintf(
        "%d KOs have zero abundance across all samples",
        sum(zero_rows)
      ))
    }
    
    # 检查列名是否唯一
    if (any(duplicated(colnames(data)))) {
      stop("Duplicate column names found in the input data")
    }
    
    return(TRUE)
  }
  
  # 文件格式验证
  if (!is.null(file)) {
    validate_file_format(file)
    message("Loading data from file...")
    file_format <- substr(file, nchar(file) - 3, nchar(file))
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
  
  # 运行输入验证
  tryCatch({
    validate_input(abundance)
  }, warning = function(w) {
    warning("Warning during data validation: ", w$message, call. = FALSE)
  }, error = function(e) {
    stop("Data validation failed: ", e$message, call. = FALSE)
  })
  
  message("Loading KEGG reference data. This might take a while...")
  load(system.file("extdata", "kegg_reference.RData", package = "ggpicrust2"))

  message("Performing KO to KEGG conversion. Please be patient, this might take a while...")
  sample_names <- colnames(abundance)[-1]
  kegg_names <- ko_to_kegg_reference[, 1]

  # 创建结果矩阵并确保正确设置列名
  kegg_abundance <- matrix(0, 
                          nrow = length(kegg_names), 
                          ncol = length(sample_names))
  colnames(kegg_abundance) <- sample_names  # 设置列名
  rownames(kegg_abundance) <- kegg_names    # 设置行名

  # 使用 tryCatch 进行错误处理
  tryCatch({
    pb <- txtProgressBar(min = 0, max = length(kegg_names), style = 3)
    start_time <- Sys.time()
    
    for (kegg_idx in seq_along(kegg_names)) {
      setTxtProgressBar(pb, kegg_idx)
      current_kegg <- kegg_names[kegg_idx]
      relevant_kos <- ko_to_kegg_reference[ko_to_kegg_reference[,1] == current_kegg, -1]
      relevant_kos <- relevant_kos[!is.na(relevant_kos)]
      
      matching_rows <- abundance[,1] %in% relevant_kos
      if (any(matching_rows)) {
        kegg_abundance[kegg_idx,] <- colSums(abundance[matching_rows, -1, drop = FALSE])
      }
    }
    
    end_time <- Sys.time()
    close(pb)
    
    message(sprintf("KO to KEGG conversion completed. Time elapsed: %.2f seconds.", 
                   as.numeric(end_time - start_time)))
    
    # 移除零丰度通路
    message("Removing KEGG pathways with zero abundance across all samples...")
    kegg_abundance <- kegg_abundance[rowSums(kegg_abundance) != 0, , drop = FALSE]
    
    # 转换为数据框时保持列名
    kegg_abundance <- as.data.frame(kegg_abundance, stringsAsFactors = FALSE)
    
    message("KEGG abundance calculation completed successfully.")
    return(kegg_abundance)
    
  }, error = function(e) {
    if (exists("pb")) close(pb)
    message("Error occurred during processing: ", e$message)
    stop(e)
  }, warning = function(w) {
    message("Warning during processing: ", w$message)
  }, finally = {
    if (exists("pb") && !inherits(pb, "try-error")) close(pb)
  })
}
