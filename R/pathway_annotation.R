#' Read and validate input file
#' @param file_path Path to the input file
#' @return Abundance data frame
#' @noRd
read_abundance_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  file_ext <- tolower(tools::file_ext(file_path))
  valid_formats <- c("txt", "tsv", "csv")
  
  if (!file_ext %in% valid_formats) {
    stop(
      "Invalid file format. Please input file in .tsv, .txt or .csv format. ",
      "The best input file format is the output file from PICRUSt2, ",
      "specifically 'pred_metagenome_unstrat.tsv'."
    )
  }
  
  delimiter <- if (file_ext == "csv") "," else "\t"
  
  tryCatch({
    abundance <- readr::read_delim(
      file_path,
      delim = delimiter,
      show_col_types = FALSE,
      progress = FALSE
    )
    
    if (ncol(abundance) < 2) {
      stop("Input file must contain at least two columns")
    }
    
    abundance %>% tibble::add_column(description = NA, .after = 1)
  }, 
  error = function(e) {
    stop("Error reading file: ", e$message)
  })
}

#' Load reference data for pathway annotation
#' @param pathway_type One of "KO", "EC", or "MetaCyc"
#' @return Reference data list
#' @noRd
load_reference_data <- function(pathway_type) {
  if (!pathway_type %in% c("KO", "EC", "MetaCyc")) {
    stop("Invalid pathway option. Please provide one of the following options: 'KO', 'EC', 'MetaCyc'.")
  }
  
  ref_file <- sprintf("%s_reference.RData", pathway_type)
  ref_path <- system.file("extdata", ref_file, package = "ggpicrust2")
  
  if (!file.exists(ref_path)) {
    stop("Reference data file not found: ", ref_file)
  }
  
  tryCatch({
    load(ref_path)
    ref_data <- get(paste0(pathway_type, "_reference"))
    
    if (pathway_type == "EC") {
      message("Note: EC description may appear to be duplicated due to shared EC numbers across different reactions.")
    }
    
    ref_data
  },
  error = function(e) {
    stop("Error loading reference data: ", e$message)
  })
}

#' Cache manager for KEGG annotations
#' @noRd
kegg_cache <- new.env(parent = emptyenv())

#' Get KEGG annotation with caching
#' @param ko_id KO identifier
#' @return KEGG annotation
#' @noRd
get_kegg_with_cache <- function(ko_id) {
  if (!exists(ko_id, envir = kegg_cache)) {
    # 添加请求间隔
    Sys.sleep(0.1)  # 100ms 延迟避免过快请求
    
    result <- with_retry({
      KEGGREST::keggGet(ko_id)
    })
    
    if (!is.null(result)) {
      assign(ko_id, result, envir = kegg_cache)
    }
  }
  get(ko_id, envir = kegg_cache, inherits = FALSE)
}

#' Retry mechanism for KEGG queries
#' @param expr Expression or function to retry
#' @param max_attempts Maximum number of retry attempts
#' @return Result or error
#' @noRd
with_retry <- function(expr, max_attempts = getOption("ggpicrust2.max_retries", 3)) {
  attempt <- 1
  last_error <- NULL
  
  while (attempt <= max_attempts) {
    result <- tryCatch({
      if (is.function(expr)) expr() else if (is.expression(expr)) eval(expr) else expr
    }, error = function(e) {
      # 特别处理 HTTP 404 错误
      if (grepl("HTTP 404", e$message)) {
        return(NULL)  # 直接返回 NULL，不再重试
      }
      
      last_error <- e
      if (attempt == max_attempts) {
        return(e)
      }
      log_message(sprintf("Attempt %d failed: %s, retrying...", attempt, e$message), "WARN")
      Sys.sleep(min(2^attempt, 30))  # 限制最大等待时间为30秒
      NULL
    })
    
    if (!is.null(result) && !inherits(result, "error")) {
      return(result)
    }
    
    attempt <- attempt + 1
  }
  
  NULL  # 如果所有尝试都失败，返回 NULL
}

#' Enhanced logging system
#' @noRd
log_message <- function(msg, level = "INFO") {
  if (getOption("ggpicrust2.verbose", default = TRUE)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] %s: %s", timestamp, level, msg))
  }
}

#' Create enhanced progress bar
#' @noRd
create_progress_bar <- function(total, format = NULL) {
  if (is.null(format)) {
    format <- paste0(
      "  [:bar] :percent | Elapsed: :elapsed | ETA: :eta\n",
      "  :current/:total (:rate/sec) | :what"
    )
  }
  
  progress::progress_bar$new(
    format = format,
    total = total,
    clear = FALSE,
    width = 80
  )
}

#' Validate input parameters
#' @noRd
validate_inputs <- function(file, pathway, daa_results_df, ko_to_kegg) {
  # File validation
  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1) {
      stop("'file' must be a single character string")
    }
  }
  
  # Pathway validation
  if (!is.null(pathway)) {
    valid_pathways <- c("KO", "EC", "MetaCyc")
    if (!pathway %in% valid_pathways) {
      stop(sprintf("'pathway' must be one of: %s", 
                  paste(valid_pathways, collapse = ", ")))
    }
  }
  
  # Data frame validation
  if (!is.null(daa_results_df)) {
    if (!is.data.frame(daa_results_df)) {
      stop("'daa_results_df' must be a data frame")
    }
    if (ko_to_kegg) {
      required_cols <- c("feature", "p_adjust")
      missing_cols <- setdiff(required_cols, colnames(daa_results_df))
      if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns: %s", 
                    paste(missing_cols, collapse = ", ")))
      }
    }
  }
}

#' Process KEGG annotations with improved error handling and caching
#' @param df Data frame with features to annotate
#' @return Annotated data frame
#' @noRd
process_kegg_annotations <- function(df) {
  if (nrow(df) == 0) {
    stop("Empty data frame provided for KEGG annotation")
  }
  
  filtered_df <- df[df$p_adjust < 0.05, ]
  if (nrow(filtered_df) == 0) {
    stop(
      "No statistically significant biomarkers found (p_adjust < 0.05).\n",
      "Consider using a less stringent threshold or reviewing your data."
    )
  }
  
  # 初始化新列
  new_cols <- c("pathway_name", "pathway_description", "pathway_class", "pathway_map")
  filtered_df[new_cols] <- NA_character_
  
  # 创建进度条
  total_features <- nrow(filtered_df)
  pb <- create_progress_bar(total_features)
  log_message("Starting KEGG annotation process")
  
  # 添加重试计数器
  retry_count <- 0
  max_retries <- 3
  
  # 处理每个特征
  for (i in seq_len(nrow(filtered_df))) {
    ko_id <- filtered_df$feature[i]
    
    # 更新进度条
    pb$tick(tokens = list(what = sprintf("Processing %s", ko_id)))
    
    # 获取 KEGG 注释
    entry <- tryCatch({
      result <- get_kegg_with_cache(ko_id)
      if (is.null(result)) {
        log_message(sprintf("No KEGG data found for %s", ko_id), "WARN")
        next
      }
      result
    }, error = function(e) {
      log_message(sprintf("Error processing %s: %s", ko_id, e$message), "ERROR")
      NULL
    })
    
    # 安全地提取数据
    if (!is.null(entry) && length(entry) > 0) {
      filtered_df$pathway_name[i] <- safe_extract(entry[[1]], "NAME", 1)
      filtered_df$pathway_description[i] <- safe_extract(entry[[1]], "DESCRIPTION", 1)
      filtered_df$pathway_class[i] <- safe_extract(entry[[1]], "CLASS", 1)
      filtered_df$pathway_map[i] <- safe_extract(entry[[1]], "PATHWAY_MAP", 1)
    }
  }
  
  log_message("KEGG annotation process completed")
  
  # 移除所有 NA 行
  filtered_df <- filtered_df[!is.na(filtered_df$pathway_name), ]
  
  if (nrow(filtered_df) == 0) {
    warning("No valid KEGG annotations found for any features")
    return(NULL)
  }
  
  filtered_df
}

#' Pathway information annotation
#'
#' @inheritParams pathway_annotation
#' @return Annotated data frame
#' @noRd
annotate_pathways <- function(abundance, pathway_type, ref_data) {
  if (nrow(abundance) == 0) {
    stop("Empty data frame provided for annotation")
  }
  
  # Vectorized operation instead of loop
  abundance$description <- ref_data[match(abundance[[1]], ref_data[[1]]), 
                                  if(pathway_type == "KO") 5 else 2]
  
  abundance
}

#' Pathway information annotation of "EC", "KO", "MetaCyc" pathway
#'
#' This function has two primary use cases:
#' 1. Annotating pathway information using the output file from PICRUSt2.
#' 2. Annotating pathway information from the output of `pathway_daa` function, and converting KO abundance to KEGG pathway abundance when `ko_to_kegg` is set to TRUE.
#'
#' @param file A character, address to store PICRUSt2 export files. Provide this parameter when using the function for the first use case.
#' @param pathway A character, consisting of "EC", "KO", "MetaCyc"
#' @param daa_results_df A data frame, output of pathway_daa. Provide this parameter when using the function for the second use case.
#' @param ko_to_kegg A logical, decide if convert KO abundance to KEGG pathway abundance. Default is FALSE. Set to TRUE when using the function for the second use case.
#'
#' @return A data frame containing pathway annotation information. The data frame has the following columns:
#' \itemize{
#'   \item \code{feature}: The feature ID of the pathway (e.g., KO, EC, or MetaCyc ID).
#'   \item \code{description}: The description or name of the pathway.
#'   \item Other columns depending on the input parameters and type of pathway.
#' }
#' If \code{ko_to_kegg} is set to TRUE, the output data frame will also include the following columns:
#' \itemize{
#'   \item \code{pathway_name}: The name of the KEGG pathway.
#'   \item \code{pathway_description}: The description of the KEGG pathway.
#'   \item \code{pathway_class}: The class of the KEGG pathway.
#'   \item \code{pathway_map}: The KEGG pathway map ID.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare the required input files and data frames
#' # Then, you can use the pathway_annotation function as follows:
#'
#' # Use case 1: Annotating pathway information using the output file from PICRUSt2
#' result1 <- pathway_annotation(file = "path/to/picrust2/export/file.txt",
#'                               pathway = "KO",
#'                               daa_results_df = NULL,
#'                               ko_to_kegg = FALSE)
#'
#' # Use case 2: Annotating pathway information from the output of pathway_daa function
#' # and converting KO abundance to KEGG pathway abundance
#' # This use case will be demonstrated using both a hypothetical example, and a real dataset.
#'
#' ## Hypothetical example
#' hypothetical_daa_results_df <- data.frame() # Replace this with your actual data frame
#' result2 <- pathway_annotation(file = NULL,
#'                               pathway = "KO",
#'                               daa_results_df = hypothetical_daa_results_df,
#'                               ko_to_kegg = TRUE)
#'
#' ## Real dataset example
#' # Load the real dataset
#' data(daa_results_df)
#' result3 <- pathway_annotation(file = NULL,
#'                               pathway = "KO",
#'                               daa_results_df = daa_results_df,
#'                               ko_to_kegg = TRUE)
#' }
pathway_annotation <- function(file = NULL,
                             pathway = NULL,
                             daa_results_df = NULL,
                             ko_to_kegg = FALSE) {
  
  # Input validation
  if (is.null(file) && is.null(daa_results_df)) {
    stop("Please input the picrust2 output or results of pathway_daa daa_results_df")
  }
  
  if (!is.null(daa_results_df) && nrow(daa_results_df) == 0) {
    stop("Input data frame is empty")
  }
  
  # Process file input
  if (!is.null(file)) {
    abundance <- read_abundance_file(file)
    ref_data <- load_reference_data(pathway)
    return(annotate_pathways(abundance, pathway, ref_data))
  }
  
  # Process DAA results
  if (!is.null(daa_results_df)) {
    if (!ko_to_kegg) {
      ref_data <- load_reference_data(pathway)
      return(annotate_pathways(daa_results_df, pathway, ref_data))
    } else {
      return(process_kegg_annotations(daa_results_df))
    }
  }
}

#' 添加安全提取函数
safe_extract <- function(list, field, index = 1) {
  tryCatch({
    if (is.null(list[[field]]) || length(list[[field]]) == 0) {
      NA_character_
    } else {
      as.character(list[[field]][index])
    }
  }, error = function(e) {
    NA_character_
  })
}
