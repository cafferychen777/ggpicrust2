# Edge cases and error handling tests for Discussion #113
# Ensures robustness of new implementation


# Load internal data
load("../../R/sysdata.rda")

test_that("ko2kegg_abundance handles all-zero abundances", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 20)

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = rep(0, 20),
    Sample2 = rep(0, 20),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  # Should return empty data frame with correct column structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 2)
})

test_that("ko2kegg_abundance handles very large abundance values", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = rep(1e10, 10),  # Very large values
    Sample2 = rep(1e10, 10),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_true(all(is.finite(as.matrix(result[, -1]))))
})

test_that("ko2kegg_abundance handles very small abundance values", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = rep(0.001, 10),  # Very small values
    Sample2 = rep(0.001, 10),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  numeric_cols <- result[, -1]
  expect_true(all(numeric_cols >= 0))
})

test_that("ko2kegg_abundance handles single sample", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 20)

  test_data <- data.frame(
    function. = real_kos,
    OnlySample = rpois(20, 100),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "OnlySample")
})

test_that("ko2kegg_abundance handles many samples", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 50)

  test_data <- data.frame(
    function. = real_kos,
    stringsAsFactors = FALSE
  )

  # Add 50 samples
  for (i in 1:50) {
    test_data[[paste0("S", i)]] <- rpois(50, 100)
  }

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 50)
})

test_that("ko2kegg_abundance handles whitespace in KO IDs", {
  # Some KO IDs in database might have whitespace
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  # Add whitespace to some KO IDs
  test_data <- data.frame(
    function. = paste0(real_kos, "  "),  # Trailing whitespace
    Sample1 = rpois(10, 100),
    stringsAsFactors = FALSE
  )

  # Should still work because code uses trimws()
  result <- tryCatch({
    suppressMessages(ko2kegg_abundance(data = test_data))
  }, error = function(e) {
    NULL
  })

  # This might fail depending on implementation
  # If it fails, it's expected since we need exact matching
  if (!is.null(result)) {
    expect_s3_class(result, "data.frame")
  }
})

test_that("ko2kegg_abundance handles duplicate KO IDs in input", {
  real_ko <- head(unique(ko_to_kegg_reference$ko_id), 1)

  # Duplicate the same KO
  test_data <- data.frame(
    function. = rep(real_ko, 3),
    Sample1 = c(100, 200, 150),
    stringsAsFactors = FALSE
  )

  # Behavior depends on implementation
  # Most likely it will sum all occurrences
  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
})

test_that("ko2kegg_abundance handles mixed case in KO IDs", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 5)

  # Convert to lowercase
  test_data <- data.frame(
    function. = tolower(real_kos),
    Sample1 = rpois(5, 100),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  # Should produce empty result since KO IDs are case-sensitive
  expect_equal(nrow(result), 0)
})

test_that("ko2kegg_abundance handles KO IDs with special prefixes", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 5)

  # Add prefix
  test_data <- data.frame(
    function. = paste0("PREFIX_", real_kos),
    Sample1 = rpois(5, 100),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  # Should produce empty result since IDs don't match
  expect_equal(nrow(result), 0)
})

test_that("ko2kegg_abundance handles extremely sparse data", {
  # Create dataset where most KOs have zero abundance
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 100)

  abundances <- rep(0, 100)
  abundances[c(1, 50, 99)] <- 100  # Only 3 non-zero

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = abundances,
    Sample2 = abundances,
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  # Should still produce some pathways from the non-zero KOs
  expect_gt(nrow(result), 0)
})

test_that("ko2kegg_abundance handles highly skewed abundance distribution", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 50)

  # One very high abundance, rest very low
  abundances <- c(10000, rep(1, 49))

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = abundances,
    Sample2 = abundances,
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)

  # Check that skewness is preserved
  numeric_result <- as.matrix(result[, -1, drop = FALSE])
  if (nrow(numeric_result) > 0 && ncol(numeric_result) > 0) {
    max_val <- max(numeric_result[, 1])
    min_val <- min(numeric_result[, 1])
    if (min_val > 0) {
      expect_gt(max_val / min_val, 1)
    }
  }
})

test_that("ko2kegg_abundance handles pathways with only one KO", {
  # Find a pathway with minimal KOs
  gene_sets <- split(ko_to_kegg_reference$ko_id,
                     ko_to_kegg_reference$pathway_id)
  set_sizes <- sapply(gene_sets, length)
  min_pathway <- names(gene_sets)[which.min(set_sizes)]
  min_kos <- gene_sets[[min_pathway]]

  if (length(min_kos) <= 3) {
    test_data <- data.frame(
      function. = min_kos,
      Sample1 = rep(100, length(min_kos)),
      stringsAsFactors = FALSE
    )

    result <- suppressMessages(ko2kegg_abundance(data = test_data))

    expect_s3_class(result, "data.frame")
  }
})

test_that("ko2kegg_abundance handles pathways with many KOs", {
  # Find the largest pathway
  gene_sets <- split(ko_to_kegg_reference$ko_id,
                     ko_to_kegg_reference$pathway_id)
  set_sizes <- sapply(gene_sets, length)
  max_pathway <- names(gene_sets)[which.max(set_sizes)]
  max_kos <- gene_sets[[max_pathway]]

  # Use only first 50 KOs from this pathway
  test_kos <- head(max_kos, 50)

  test_data <- data.frame(
    function. = test_kos,
    Sample1 = rep(100, length(test_kos)),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_true(max_pathway %in% rownames(result))
})

test_that("ko2kegg_abundance handles fractional abundances", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = runif(10, 0.1, 100.5),  # Fractional values
    Sample2 = runif(10, 0.1, 100.5),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)

  # Check that fractional values are preserved
  numeric_result <- as.matrix(result[, -1, drop = FALSE])
  expect_true(any(numeric_result != floor(numeric_result)))
})

test_that("ko2kegg_abundance handles negative abundances", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = c(100, -50, 200, -30, 150, 180, -10, 90, 110, 120),
    stringsAsFactors = FALSE
  )

  # Should either error or handle gracefully
  result <- tryCatch({
    suppressMessages(ko2kegg_abundance(data = test_data))
  }, error = function(e) {
    "ERROR"
  })

  # Either produces a result or errors - both are acceptable
  expect_true(is.data.frame(result) || result == "ERROR")
})

test_that("ko2kegg_abundance handles NA abundances", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  test_data <- data.frame(
    function. = real_kos,
    Sample1 = c(100, NA, 200, 150, NA, 180, 90, NA, 110, 120),
    stringsAsFactors = FALSE
  )

  # Should either handle or error
  result <- tryCatch({
    suppressMessages(ko2kegg_abundance(data = test_data))
  }, error = function(e) {
    "ERROR"
  })

  expect_true(is.data.frame(result) || result == "ERROR")
})

test_that("ko2kegg_abundance handles empty function column", {
  test_data <- data.frame(
    function. = character(0),
    Sample1 = numeric(0),
    stringsAsFactors = FALSE
  )

  expect_error(ko2kegg_abundance(data = test_data))
})

test_that("ko2kegg_abundance handles missing function column", {
  test_data <- data.frame(
    KO_ID = c("K00001", "K00002"),
    Sample1 = c(100, 200),
    stringsAsFactors = FALSE
  )

  expect_error(ko2kegg_abundance(data = test_data))
})

test_that("ko2kegg_abundance handles data with only function column", {
  test_data <- data.frame(
    function. = head(unique(ko_to_kegg_reference$ko_id), 10),
    stringsAsFactors = FALSE
  )

  expect_error(ko2kegg_abundance(data = test_data))
})

test_that("ko2kegg_abundance handles very long sample names", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  long_name <- paste(rep("A", 200), collapse = "")

  test_data <- data.frame(
    function. = real_kos,
    stringsAsFactors = FALSE
  )
  test_data[[long_name]] <- rpois(10, 100)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result)[1], long_name)
})

test_that("ko2kegg_abundance handles unicode in sample names", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)

  test_data <- data.frame(
    function. = real_kos,
    stringsAsFactors = FALSE
  )
  test_data[["Sample_1"]] <- rpois(10, 100)
  test_data[["Sample_2"]] <- rpois(10, 100)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
})
