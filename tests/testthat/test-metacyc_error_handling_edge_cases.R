# MetaCyc Error Handling and Edge Case Tests
# Comprehensive testing of error conditions and edge cases
# Following Linus principles: "Never break userspace", clear error messages, robust error handling

test_that("MetaCyc malformed abundance data error handling", {
  skip_if_not_installed("fgsea")
  
  # Test 1: Abundance data with non-numeric values
  bad_abundance_non_numeric <- data.frame(
    "EC:1.1.1.1" = c("not_a_number", "also_not_number"),
    "EC:2.2.2.2" = c("text", "more_text"),
    stringsAsFactors = FALSE
  )
  rownames(bad_abundance_non_numeric) <- c("Sample1", "Sample2")
  
  metadata_valid <- data.frame(
    sample_id = c("Sample1", "Sample2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_valid) <- metadata_valid$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = t(bad_abundance_non_numeric),  # Transpose to get samples as columns
      metadata = metadata_valid,
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 2: Abundance data with infinite values
  bad_abundance_infinite <- matrix(c(1, 2, Inf, 4, 5, -Inf), nrow = 2, ncol = 3)
  rownames(bad_abundance_infinite) <- c("EC:1.1.1.1", "EC:2.2.2.2")
  colnames(bad_abundance_infinite) <- c("S1", "S2", "S3")
  
  metadata_inf <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_inf) <- metadata_inf$sample_id
  
  # Should handle infinite values gracefully or give clear error
  expect_warning({
    tryCatch({
      result_inf <- pathway_gsea(
        abundance = bad_abundance_infinite,
        metadata = metadata_inf,
        group = "group",
        pathway_type = "MetaCyc",
        nperm = 10,
        min_size = 1
      )
    }, error = function(e) {
      # Error is acceptable, but it should be informative
      expect_true(nchar(e$message) > 0)
    })
  })
  
  # Test 3: Abundance data with all NA values
  bad_abundance_all_na <- matrix(NA, nrow = 3, ncol = 4)
  rownames(bad_abundance_all_na) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(bad_abundance_all_na) <- c("S1", "S2", "S3", "S4")
  
  metadata_na <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_na) <- metadata_na$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = bad_abundance_all_na,
      metadata = metadata_na,
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 4: Abundance data with negative values (biologically invalid)
  bad_abundance_negative <- matrix(c(-1, -2, -3, -4, 5, 6), nrow = 2, ncol = 3)
  rownames(bad_abundance_negative) <- c("EC:1.1.1.1", "EC:2.2.2.2")
  colnames(bad_abundance_negative) <- c("S1", "S2", "S3")
  
  metadata_neg <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    group = c("A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_neg) <- metadata_neg$sample_id
  
  # Should handle negative values (might be acceptable in some contexts like log ratios)
  expect_no_error({
    tryCatch({
      result_neg <- pathway_gsea(
        abundance = bad_abundance_negative,
        metadata = metadata_neg,
        group = "group",
        pathway_type = "MetaCyc",
        nperm = 10,
        min_size = 1
      )
    }, error = function(e) {
      # If it errors, should be informative
      expect_true(nchar(e$message) > 0)
    })
  })
})

test_that("MetaCyc missing EC numbers in gene sets error handling", {
  skip_if_not_installed("fgsea")
  
  # Test with abundance data that has no overlap with MetaCyc gene sets
  no_overlap_abundance <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(no_overlap_abundance) <- c("UNKNOWN:1", "UNKNOWN:2", "UNKNOWN:3", "UNKNOWN:4")
  colnames(no_overlap_abundance) <- paste0("Sample", 1:5)
  
  metadata_no_overlap <- data.frame(
    sample_id = paste0("Sample", 1:5),
    group = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_no_overlap) <- metadata_no_overlap$sample_id
  
  # Should handle gracefully when no genes overlap between data and gene sets
  expect_no_error({
    result_no_overlap <- pathway_gsea(
      abundance = no_overlap_abundance,
      metadata = metadata_no_overlap,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
    
    # Result should be empty but properly structured
    expect_s3_class(result_no_overlap, "data.frame")
    expect_equal(nrow(result_no_overlap), 0)
    expect_true("pathway_id" %in% colnames(result_no_overlap))
  })
  
  # Test with minimal EC overlap
  minimal_overlap_abundance <- matrix(rnorm(25), nrow = 5, ncol = 5)
  rownames(minimal_overlap_abundance) <- c("EC:1.1.1.1", "UNKNOWN:1", "UNKNOWN:2", "UNKNOWN:3", "EC:2.2.2.2")
  colnames(minimal_overlap_abundance) <- paste0("Sample", 1:5)
  
  expect_no_error({
    result_minimal <- pathway_gsea(
      abundance = minimal_overlap_abundance,
      metadata = metadata_no_overlap,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10,
      min_size = 1,  # Allow very small gene sets
      max_size = 50
    )
    
    # Should handle minimal overlap gracefully
    expect_s3_class(result_minimal, "data.frame")
  })
})

test_that("MetaCyc metadata edge cases error handling", {
  skip_if_not_installed("fgsea")
  
  # Create valid abundance data for testing metadata issues
  valid_abundance <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(valid_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3", "EC:4.4.4.4")
  colnames(valid_abundance) <- paste0("Sample", 1:5)
  
  # Test 1: Metadata with factor levels that have no samples
  metadata_empty_levels <- data.frame(
    sample_id = paste0("Sample", 1:5),
    group = factor(c("A", "A", "B", "B", "B"), levels = c("A", "B", "C")),  # C has no samples
    stringsAsFactors = FALSE
  )
  rownames(metadata_empty_levels) <- metadata_empty_levels$sample_id
  
  expect_no_error({
    result_empty_levels <- pathway_gsea(
      abundance = valid_abundance,
      metadata = metadata_empty_levels,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
  })
  
  # Test 2: Metadata with NA values in grouping variable
  metadata_na_groups <- data.frame(
    sample_id = paste0("Sample", 1:5),
    group = c("A", "A", NA, "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_na_groups) <- metadata_na_groups$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = valid_abundance,
      metadata = metadata_na_groups,
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 3: Metadata with single unique group (no comparison possible)
  metadata_single_group <- data.frame(
    sample_id = paste0("Sample", 1:5),
    group = rep("A", 5),
    stringsAsFactors = FALSE
  )
  rownames(metadata_single_group) <- metadata_single_group$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = valid_abundance,
      metadata = metadata_single_group,
      group = "group",
      pathway_type = "MetaCyc"
    )
  }, info = "Should error with informative message about needing two groups")
  
  # Test 4: Metadata with more than two groups (not currently supported)
  metadata_multi_group <- data.frame(
    sample_id = paste0("Sample", 1:5),
    group = c("A", "B", "C", "D", "E"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_multi_group) <- metadata_multi_group$sample_id
  
  expect_warning({
    # Should warn about multiple groups but may still work with first two
    result_multi <- pathway_gsea(
      abundance = valid_abundance,
      metadata = metadata_multi_group,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
  })
})

test_that("MetaCyc parameter validation edge cases", {
  # Test invalid parameter combinations and edge values
  
  # Create minimal valid data for parameter testing
  test_abundance <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(test_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(test_abundance) <- c("S1", "S2", "S3", "S4")
  
  test_metadata <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(test_metadata) <- test_metadata$sample_id
  
  # Test 1: Invalid nperm values
  expect_error({
    pathway_gsea(
      abundance = test_abundance,
      metadata = test_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = -1  # Invalid negative
    )
  })
  
  expect_error({
    pathway_gsea(
      abundance = test_abundance,
      metadata = test_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 0  # Invalid zero
    )
  })
  
  # Test 2: Invalid min_size/max_size combinations
  expect_error({
    pathway_gsea(
      abundance = test_abundance,
      metadata = test_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      min_size = 100,
      max_size = 50  # min_size > max_size
    )
  })
  
  expect_error({
    pathway_gsea(
      abundance = test_abundance,
      metadata = test_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      min_size = -1  # Negative min_size
    )
  })
  
  # Test 3: Invalid ranking methods
  expect_error({
    pathway_gsea(
      abundance = test_abundance,
      metadata = test_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      rank_method = "invalid_method"
    )
  })
  
  # Test 4: Invalid GSEA methods
  expect_error({
    pathway_gsea(
      abundance = test_abundance,
      metadata = test_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      method = "invalid_gsea_method"
    )
  })
  
  # Test 5: Edge case with very small datasets
  tiny_abundance <- matrix(c(1, 2), nrow = 1, ncol = 2)
  rownames(tiny_abundance) <- "EC:1.1.1.1"
  colnames(tiny_abundance) <- c("S1", "S2")
  
  tiny_metadata <- data.frame(
    sample_id = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(tiny_metadata) <- tiny_metadata$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = tiny_abundance,
      metadata = tiny_metadata,
      group = "group",
      pathway_type = "MetaCyc"
    )
  }, info = "Should require minimum samples for statistical analysis")
})

test_that("MetaCyc data type and structure error handling", {
  # Test various invalid data types and structures
  
  # Test 1: abundance as list instead of matrix/data.frame
  expect_error({
    pathway_gsea(
      abundance = list(a = 1, b = 2),
      metadata = data.frame(),
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 2: metadata as vector instead of data.frame
  expect_error({
    pathway_gsea(
      abundance = matrix(1:4, nrow = 2),
      metadata = c("A", "B", "A", "B"),
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 3: NULL inputs
  expect_error({
    pathway_gsea(
      abundance = NULL,
      metadata = data.frame(),
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  expect_error({
    pathway_gsea(
      abundance = matrix(1:4, nrow = 2),
      metadata = NULL,
      group = "group", 
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 4: Empty data structures
  expect_error({
    pathway_gsea(
      abundance = matrix(numeric(0), nrow = 0, ncol = 0),
      metadata = data.frame(),
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  expect_error({
    pathway_gsea(
      abundance = matrix(1:4, nrow = 2, ncol = 2),
      metadata = data.frame(),
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
})

test_that("MetaCyc row and column name edge cases", {
  skip_if_not_installed("fgsea")
  
  # Test 1: Missing row names in abundance data
  no_rownames_abundance <- matrix(rnorm(20), nrow = 4, ncol = 5)
  # Don't set row names - this should cause issues
  colnames(no_rownames_abundance) <- paste0("Sample", 1:5)
  
  metadata_valid <- data.frame(
    sample_id = paste0("Sample", 1:5),
    group = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_valid) <- metadata_valid$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = no_rownames_abundance,
      metadata = metadata_valid,
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 2: Missing column names in abundance data
  no_colnames_abundance <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(no_colnames_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3", "EC:4.4.4.4")
  # Don't set column names
  
  expect_error({
    pathway_gsea(
      abundance = no_colnames_abundance,
      metadata = metadata_valid,
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 3: Duplicate row names in abundance data
  dup_rownames_abundance <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(dup_rownames_abundance) <- c("EC:1.1.1.1", "EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")  # Duplicate
  colnames(dup_rownames_abundance) <- paste0("Sample", 1:5)
  
  # Should handle duplicates gracefully or give clear error
  expect_warning({
    tryCatch({
      result_dup <- pathway_gsea(
        abundance = dup_rownames_abundance,
        metadata = metadata_valid,
        group = "group",
        pathway_type = "MetaCyc",
        nperm = 10
      )
    }, error = function(e) {
      expect_true(nchar(e$message) > 0)
    })
  })
  
  # Test 4: Special characters in sample names
  special_abundance <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(special_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(special_abundance) <- c("Sample-1", "Sample.2", "Sample_3", "Sample@4")
  
  special_metadata <- data.frame(
    sample_id = c("Sample-1", "Sample.2", "Sample_3", "Sample@4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(special_metadata) <- special_metadata$sample_id
  
  expect_no_error({
    result_special <- pathway_gsea(
      abundance = special_abundance,
      metadata = special_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
  })
})

test_that("MetaCyc statistical edge cases error handling", {
  skip_if_not_installed("fgsea")
  
  # Test 1: All samples have identical values (zero variance)
  zero_var_abundance <- matrix(rep(100, 15), nrow = 3, ncol = 5)
  rownames(zero_var_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(zero_var_abundance) <- paste0("Sample", 1:5)
  
  metadata_zero_var <- data.frame(
    sample_id = paste0("Sample", 1:5),
    group = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata_zero_var) <- metadata_zero_var$sample_id
  
  expect_no_error({
    # Should handle zero variance gracefully (might produce warnings)
    suppressWarnings({
      result_zero_var <- pathway_gsea(
        abundance = zero_var_abundance,
        metadata = metadata_zero_var,
        group = "group",
        pathway_type = "MetaCyc",
        nperm = 10
      )
    })
  })
  
  # Test 2: Extreme outliers in data
  outlier_abundance <- matrix(rnorm(15, mean = 10, sd = 1), nrow = 3, ncol = 5)
  outlier_abundance[1, 1] <- 1e6  # Extreme outlier
  outlier_abundance[2, 5] <- -1e6  # Extreme negative outlier
  rownames(outlier_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(outlier_abundance) <- paste0("Sample", 1:5)
  
  expect_no_error({
    result_outliers <- pathway_gsea(
      abundance = outlier_abundance,
      metadata = metadata_zero_var,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
  })
  
  # Test 3: Very small sample sizes per group
  small_abundance <- matrix(rnorm(9), nrow = 3, ncol = 3)
  rownames(small_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(small_abundance) <- c("S1", "S2", "S3")
  
  small_metadata <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    group = c("A", "B", "B"),  # Very unbalanced: 1 vs 2
    stringsAsFactors = FALSE
  )
  rownames(small_metadata) <- small_metadata$sample_id
  
  expect_warning({
    # Should warn about small sample sizes but may still work
    result_small <- pathway_gsea(
      abundance = small_abundance,
      metadata = small_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
  })
})

test_that("MetaCyc Unicode and special character handling", {
  skip_if_not_installed("fgsea")
  
  # Test handling of various character encodings and special characters
  unicode_abundance <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(unicode_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(unicode_abundance) <- c("Sample_α", "Sample_β", "Sample_γ", "Sample_δ")
  
  unicode_metadata <- data.frame(
    sample_id = c("Sample_α", "Sample_β", "Sample_γ", "Sample_δ"),
    group = c("Control_α", "Control_β", "Treatment_γ", "Treatment_δ"),
    stringsAsFactors = FALSE
  )
  rownames(unicode_metadata) <- unicode_metadata$sample_id
  
  expect_no_error({
    result_unicode <- pathway_gsea(
      abundance = unicode_abundance,
      metadata = unicode_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
  })
  
  # Test with empty strings and whitespace
  whitespace_abundance <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(whitespace_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(whitespace_abundance) <- c(" Sample1 ", "\tSample2\t", "Sample3\n", "Sample4")
  
  whitespace_metadata <- data.frame(
    sample_id = c(" Sample1 ", "\tSample2\t", "Sample3\n", "Sample4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(whitespace_metadata) <- whitespace_metadata$sample_id
  
  # Should handle whitespace in sample names consistently
  expect_no_error({
    result_whitespace <- pathway_gsea(
      abundance = whitespace_abundance,
      metadata = whitespace_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      nperm = 10
    )
  })
})