# Helper function to suppress all messages and warnings during tests
quiet_run <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}

test_that("pathway_annotation basic functionality works", {
  # Setup test data
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    feature = c("K00001", "K00002"),
    sample1 = c(1, 2),
    sample2 = c(3, 4)
  )
  quiet_run(write.table(test_data, temp_file, sep = "\t", row.names = FALSE))

  result_ko <- quiet_run(
    pathway_annotation(
      file = temp_file,
      pathway = "KO",
      daa_results_df = NULL,
      ko_to_kegg = FALSE
    )
  )

  expect_s3_class(result_ko, "data.frame")
  expect_true("description" %in% colnames(result_ko))
  expect_equal(nrow(result_ko), 2)

  unlink(temp_file)
})

test_that("pathway_annotation handles invalid inputs correctly", {
  # Test missing both required inputs
  expect_error(
    pathway_annotation(),
    "Please input the picrust2 output or results of pathway_daa daa_results_df"
  )

  # Test invalid pathway type
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    feature = c("K00001", "K00002"),
    sample1 = c(1, 2),
    sample2 = c(3, 4)
  )
  quiet_run(write.table(test_data, temp_file, sep = "\t", row.names = FALSE))

  expect_error(
    pathway_annotation(file = temp_file, pathway = "INVALID"),
    "Invalid pathway option. Please provide one of the following options: 'KO', 'EC', 'MetaCyc'"
  )

  unlink(temp_file)
})

test_that("pathway_annotation works with daa_results_df input", {
  # Create mock daa_results_df
  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  result_daa <- quiet_run(
    pathway_annotation(
      file = NULL,
      pathway = "KO",
      daa_results_df = test_daa_df,
      ko_to_kegg = FALSE
    )
  )

  expect_s3_class(result_daa, "data.frame")
  expect_true("description" %in% colnames(result_daa))
  expect_equal(nrow(result_daa), 2)
})

test_that("pathway_annotation handles different file formats", {
  formats <- c(".tsv", ".txt", ".csv")

  for(fmt in formats) {
    temp_file <- tempfile(fileext = fmt)
    test_data <- data.frame(
      feature = c("K00001", "K00002"),
      sample1 = c(1, 2),
      sample2 = c(3, 4)
    )

    # 根据文件格式使用正确的分隔符
    delimiter <- if(fmt == ".csv") "," else "\t"
    quiet_run(
      write.table(
        test_data,
        temp_file,
        sep = delimiter,
        row.names = FALSE,
        quote = FALSE
      )
    )

    result <- quiet_run(
      pathway_annotation(
        file = temp_file,
        pathway = "KO",
        daa_results_df = NULL,
        ko_to_kegg = FALSE
      )
    )

    expect_s3_class(result, "data.frame")
    unlink(temp_file)
  }
})

test_that("ko_to_kegg functionality works", {
  skip_if_offline()

  # Create mock daa_results_df with significant results
  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  result_kegg <- quiet_run(
    pathway_annotation(
      file = NULL,
      pathway = "KO",
      daa_results_df = test_daa_df,
      ko_to_kegg = TRUE
    )
  )

  # Check KEGG specific columns exist
  expect_true(all(c(
    "pathway_name",
    "pathway_description",
    "pathway_class",
    "pathway_map"
  ) %in% colnames(result_kegg)))

  # Check that the result is a data frame
  expect_s3_class(result_kegg, "data.frame")

  # Check that we have the same number of rows as significant results
  expect_equal(nrow(result_kegg), nrow(test_daa_df))
})

test_that("pathway_annotation handles no significant results correctly", {
  # Create mock daa_results_df with no significant results
  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.06, 0.07),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_annotation(
      file = NULL,
      pathway = "KO",
      daa_results_df = test_daa_df,
      ko_to_kegg = TRUE
    ),
    "No statistically significant biomarkers found"
  )
})

test_that("pathway_annotation handles empty data frames correctly", {
  empty_df <- data.frame(
    feature = character(0),
    p_adjust = numeric(0),
    stringsAsFactors = FALSE
  )

  expect_error(
    quiet_run(
      pathway_annotation(
        file = NULL,
        pathway = "KO",
        daa_results_df = empty_df,
        ko_to_kegg = FALSE
      )
    ),
    "Input data frame is empty"
  )
})

test_that("pathway_annotation handles malformed input files", {
  # Test file with missing columns
  temp_file <- tempfile(fileext = ".tsv")
  bad_data <- data.frame(feature = c("K00001", "K00002"))
  quiet_run(write.table(bad_data, temp_file, sep = "\t", row.names = FALSE))

  expect_error(
    quiet_run(
      pathway_annotation(
        file = temp_file,
        pathway = "KO"
      )
    ),
    "Input file must contain at least two columns"
  )

  unlink(temp_file)

  # Test non-existent file
  expect_error(
    quiet_run(
      pathway_annotation(
        file = "nonexistent.tsv",
        pathway = "KO"
      )
    ),
    "File does not exist"
  )
})

test_that("pathway_annotation works with all pathway types", {
  test_daa_df <- data.frame(
    feature = c("K00001", "EC:1.1.1.1", "PWY-7219"),
    p_adjust = c(0.04, 0.03, 0.02),
    stringsAsFactors = FALSE
  )

  # Test KO pathway
  ko_result <- quiet_run(
    pathway_annotation(
      file = NULL,
      pathway = "KO",
      daa_results_df = test_daa_df,
      ko_to_kegg = FALSE
    )
  )
  expect_s3_class(ko_result, "data.frame")
  expect_true("description" %in% colnames(ko_result))

  # Test EC pathway
  ec_result <- quiet_run(
    pathway_annotation(
      file = NULL,
      pathway = "EC",
      daa_results_df = test_daa_df,
      ko_to_kegg = FALSE
    )
  )
  expect_s3_class(ec_result, "data.frame")
  expect_true("description" %in% colnames(ec_result))

  # Test MetaCyc pathway
  metacyc_result <- quiet_run(
    pathway_annotation(
      file = NULL,
      pathway = "MetaCyc",
      daa_results_df = test_daa_df,
      ko_to_kegg = FALSE
    )
  )
  expect_s3_class(metacyc_result, "data.frame")
  expect_true("description" %in% colnames(metacyc_result))
})

test_that("pathway_annotation preserves input data structure", {
  # Test with additional columns in input
  test_data <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    log2fc = c(1.5, -2.0),
    extra_col = c("A", "B"),
    stringsAsFactors = FALSE
  )

  result <- quiet_run(
    pathway_annotation(
      file = NULL,
      pathway = "KO",
      daa_results_df = test_data,
      ko_to_kegg = FALSE
    )
  )

  expect_true(all(c("log2fc", "extra_col") %in% colnames(result)))
  expect_equal(nrow(result), nrow(test_data))
})

test_that("pathway_annotation handles special characters in input", {
  # Test with special characters in feature IDs
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    feature = c("K00001-special", "K00002_test"),
    sample1 = c(1, 2),
    sample2 = c(3, 4)
  )
  quiet_run(write.table(test_data, temp_file, sep = "\t", row.names = FALSE))

  result <- quiet_run(
    pathway_annotation(
      file = temp_file,
      pathway = "KO",
      daa_results_df = NULL,
      ko_to_kegg = FALSE
    )
  )

  expect_s3_class(result, "data.frame")
  unlink(temp_file)
})

test_that("pathway_annotation handles KEGG query failures gracefully", {
  skip_if_offline()

  # 使用无效��� KO IDs
  test_daa_df <- data.frame(
    feature = c("K99999", "K88888"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  # 预期这个调用会返回 NULL 或抛出警告
  expect_warning(
    result <- pathway_annotation(
      file = NULL,
      pathway = "KO",
      daa_results_df = test_daa_df,
      ko_to_kegg = TRUE
    ),
    "No valid KEGG annotations found for any features"
  )

  # 如果返回结果为 NULL，测试通过
  expect_true(is.null(result))
})

test_that("pathway_annotation handles large datasets efficiently", {
  testthat::skip_if_not(identical(Sys.getenv("TEST_SLOW"), "true"),
                        "Skipping slow tests")

  # Create large test dataset
  large_data <- data.frame(
    feature = paste0("K", sprintf("%05d", 1:1000)),
    p_adjust = runif(1000, 0, 0.1),
    stringsAsFactors = FALSE
  )

  # Measure execution time
  start_time <- Sys.time()
  result <- quiet_run(
    pathway_annotation(
      file = NULL,
      pathway = "KO",
      daa_results_df = large_data,
      ko_to_kegg = FALSE
    )
  )
  end_time <- Sys.time()

  # Test should complete within reasonable time (e.g., 30 seconds)
  expect_true(as.numeric(end_time - start_time, units = "secs") < 30)
  expect_equal(nrow(result), nrow(large_data))
})

# 测试缓存系统
test_that("KEGG cache system works correctly", {
  skip_if_offline()

  # 清除现有缓存
  rm(list = ls(envir = kegg_cache), envir = kegg_cache)

  test_ko <- "K00001"

  # 第一次查询应该访问KEGG
  start_time <- Sys.time()
  result1 <- quiet_run(get_kegg_with_cache(test_ko))
  first_query_time <- as.numeric(Sys.time() - start_time, units = "secs")

  # 第二次查询应该使用缓存
  start_time <- Sys.time()
  result2 <- quiet_run(get_kegg_with_cache(test_ko))
  second_query_time <- as.numeric(Sys.time() - start_time, units = "secs")

  # 验证结果一致性
  expect_equal(result1, result2)
  # 验证缓存查询更快
  expect_true(second_query_time < first_query_time)
  # 验证数据已缓存
  expect_true(exists(test_ko, envir = kegg_cache))
})

# 测试重试机制
test_that("retry mechanism works correctly", {
  # 在测试环境中创建计数器
  test_env <- new.env()
  test_env$attempt_count <- 0

  # 创建测试函数
  test_fun <- function() {
    test_env$attempt_count <- test_env$attempt_count + 1
    if (test_env$attempt_count < 2) {
      stop("Simulated failure")
    }
    "success"
  }

  # 测试重试成功
  result <- quiet_run(with_retry(test_fun()))
  expect_equal(result, "success")
  expect_equal(test_env$attempt_count, 2)

  # 测试达到最大��试次数
  test_env$attempt_count <- 0
  always_fails <- function() {
    test_env$attempt_count <- test_env$attempt_count + 1
    stop("Always fails")
  }

  # 测试最大重试次数后返回 NULL
  result <- quiet_run(with_retry(always_fails(), max_attempts = 2))
  expect_null(result)
  expect_equal(test_env$attempt_count, 2)
})

# 测试日志系统
test_that("logging system works correctly", {
  # 测试不同日志级别
  expect_message(
    log_message("Test info message", "INFO"),
    "INFO: Test info message"
  )

  expect_message(
    log_message("Test warning message", "WARN"),
    "WARN: Test warning message"
  )

  # 测试日志禁用
  old_verbose <- options(ggpicrust2.verbose = FALSE)
  expect_no_message(log_message("This should not appear"))
  options(old_verbose)
})

# 测试进度条
test_that("progress bar works correctly", {
  pb <- create_progress_bar(100)
  expect_s3_class(pb, "progress_bar")
  expect_equal(pb$.__enclos_env__$private$total, 100)

  # 测试进度更新
  expect_no_error(pb$tick())
  expect_no_error(pb$tick(tokens = list(what = "Testing")))
})

# 测试输入验证
test_that("input validation works correctly", {
  # 测试文件路径验证
  expect_error(
    validate_inputs(file = 123, pathway = "KO", NULL, FALSE),
    "'file' must be a single character string"
  )

  # 测试pathway类型验证
  expect_error(
    validate_inputs(file = "test.txt", pathway = "invalid", NULL, FALSE),
    "must be one of: KO, EC, MetaCyc"
  )

  # 测试数据框验证
  invalid_df <- list(a = 1)
  expect_error(
    validate_inputs(NULL, "KO", invalid_df, TRUE),
    "'daa_results_df' must be a data frame"
  )
})

# 测试配置系统
test_that("configuration system works correctly", {
  # 保存原始设置
  old_options <- options()

  # 测试默认值
  expect_true(getOption("ggpicrust2.cache_enabled", TRUE))
  expect_equal(getOption("ggpicrust2.max_retries", 3), 3)

  # 测试自定义设置
  options(ggpicrust2.max_retries = 5)
  expect_equal(getOption("ggpicrust2.max_retries"), 5)

  # 恢复原始设置
  options(old_options)
})

# 测试KEGG注释处理
test_that("KEGG annotation processing works correctly", {
  skip_if_offline()

  # 测试有效数据
  valid_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.01, 0.02)
  )

  result <- quiet_run(process_kegg_annotations(valid_df))
  expect_true(all(c(
    "pathway_name",
    "pathway_description",
    "pathway_class",
    "pathway_map"
  ) %in% colnames(result)))

  # 测试无显著性结果
  nonsig_df <- data.frame(
    feature = "K00001",
    p_adjust = 0.06
  )
  expect_error(
    process_kegg_annotations(nonsig_df),
    "No statistically significant biomarkers found"
  )

  # 测试空数据框
  empty_df <- data.frame(
    feature = character(0),
    p_adjust = numeric(0)
  )
  expect_error(
    process_kegg_annotations(empty_df),
    "Empty data frame provided for KEGG annotation"
  )
})

# 测试大数据集处理
test_that("large dataset processing works efficiently", {
  skip_if_not(identical(Sys.getenv("TEST_SLOW"), "true"),
              "Skipping slow tests")

  # 创建大数据集
  n_rows <- 100
  large_df <- data.frame(
    feature = paste0("K", sprintf("%05d", 1:n_rows)),
    p_adjust = runif(n_rows, 0, 0.04)  # 确保有显著性结果
  )

  # 测试处理时间
  start_time <- Sys.time()
  result <- quiet_run(process_kegg_annotations(large_df))
  end_time <- Sys.time()

  # 验证处理时间在合理范围内（例如每个条目平均不超过1秒）
  processing_time <- as.numeric(end_time - start_time, units = "secs")
  expect_true(processing_time/n_rows < 1)

  # 验证结果完整性
  expect_equal(nrow(result), nrow(large_df))
  expect_true(all(c(
    "pathway_name",
    "pathway_description",
    "pathway_class",
    "pathway_map"
  ) %in% colnames(result)))
})
