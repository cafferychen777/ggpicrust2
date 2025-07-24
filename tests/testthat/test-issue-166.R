# Test for Issue #166: Add abundance statistics to pathway_daa and pathway_errorbar_table

test_that("pathway_daa include_abundance_stats parameter works correctly", {
  # Create test data
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(15, 25, 35),
    sample3 = c(30, 40, 50),
    sample4 = c(35, 45, 55),
    row.names = c("pathway1", "pathway2", "pathway3")
  )
  
  metadata <- data.frame(
    sample = paste0("sample", 1:4),
    group = factor(c("control", "control", "treatment", "treatment"))
  )
  
  # Test with include_abundance_stats = FALSE (default)
  result_basic <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2",
    include_abundance_stats = FALSE
  )
  
  # Should not have abundance statistics columns
  abundance_cols <- c("mean_rel_abundance_group1", "sd_rel_abundance_group1",
                     "mean_rel_abundance_group2", "sd_rel_abundance_group2", 
                     "log2_fold_change")
  
  expect_false(any(abundance_cols %in% colnames(result_basic)))
  
  # Test with include_abundance_stats = TRUE
  result_enhanced <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2",
    include_abundance_stats = TRUE
  )
  
  # Should have all abundance statistics columns
  expect_true(all(abundance_cols %in% colnames(result_enhanced)))
  
  # Check that values are numeric and not all NA
  for (col in abundance_cols) {
    expect_true(is.numeric(result_enhanced[[col]]))
    expect_false(all(is.na(result_enhanced[[col]])))
  }
  
  # Check that log2_fold_change values are reasonable
  expect_true(all(is.finite(result_enhanced$log2_fold_change)))
  
  # Check that standard deviations are non-negative
  expect_true(all(result_enhanced$sd_rel_abundance_group1 >= 0, na.rm = TRUE))
  expect_true(all(result_enhanced$sd_rel_abundance_group2 >= 0, na.rm = TRUE))
  
  # Check that relative abundances are between 0 and 1
  expect_true(all(result_enhanced$mean_rel_abundance_group1 >= 0 & 
                 result_enhanced$mean_rel_abundance_group1 <= 1, na.rm = TRUE))
  expect_true(all(result_enhanced$mean_rel_abundance_group2 >= 0 & 
                 result_enhanced$mean_rel_abundance_group2 <= 1, na.rm = TRUE))
})

test_that("calculate_abundance_stats helper function works correctly", {
  # Create test data
  abundance <- data.frame(
    sample1 = c(100, 200),
    sample2 = c(150, 250),
    sample3 = c(300, 400),
    sample4 = c(350, 450),
    row.names = c("feature1", "feature2")
  )
  
  metadata <- data.frame(
    sample = paste0("sample", 1:4),
    group_col = factor(c("A", "A", "B", "B"))
  )
  
  result <- calculate_abundance_stats(
    abundance = abundance,
    metadata = metadata,
    group = "group_col",
    features = rownames(abundance),
    group1 = "A",
    group2 = "B"
  )
  
  # Check structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  
  expected_cols <- c("feature", "mean_rel_abundance_group1", "sd_rel_abundance_group1",
                    "mean_rel_abundance_group2", "sd_rel_abundance_group2", 
                    "log2_fold_change")
  expect_true(all(expected_cols %in% colnames(result)))
  
  # Check that features are correct
  expect_equal(sort(result$feature), sort(rownames(abundance)))
  
  # Check that calculations are reasonable
  expect_true(all(is.finite(result$log2_fold_change)))
  expect_true(all(result$sd_rel_abundance_group1 >= 0))
  expect_true(all(result$sd_rel_abundance_group2 >= 0))
})

test_that("pathway_errorbar_table function works correctly", {
  # Create test data
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(15, 25, 35),
    sample3 = c(30, 40, 50),
    sample4 = c(35, 45, 55),
    row.names = c("pathway1", "pathway2", "pathway3")
  )
  
  metadata <- data.frame(
    sample = paste0("sample", 1:4),
    group = factor(c("control", "control", "treatment", "treatment"))
  )
  
  # Get DAA results
  daa_results <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2"
  )
  
  # Filter for single method
  daa_single_method <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]
  
  # Test pathway_errorbar_table
  result <- pathway_errorbar_table(
    abundance = abundance,
    daa_results_df = daa_single_method,
    Group = metadata$group,
    p_values_threshold = 1.0  # Include all features
  )
  
  # Check structure
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  
  expected_cols <- c("feature", "group1", "group2", 
                    "mean_rel_abundance_group1", "sd_rel_abundance_group1",
                    "mean_rel_abundance_group2", "sd_rel_abundance_group2", 
                    "log2_fold_change", "p_adjust")
  expect_true(all(expected_cols %in% colnames(result)))
  
  # Check that values are reasonable
  expect_true(all(is.finite(result$log2_fold_change)))
  expect_true(all(result$sd_rel_abundance_group1 >= 0))
  expect_true(all(result$sd_rel_abundance_group2 >= 0))
  expect_true(all(result$mean_rel_abundance_group1 >= 0 & result$mean_rel_abundance_group1 <= 1))
  expect_true(all(result$mean_rel_abundance_group2 >= 0 & result$mean_rel_abundance_group2 <= 1))
})

test_that("backward compatibility is maintained", {
  # Create test data
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(15, 25, 35),
    sample3 = c(30, 40, 50),
    sample4 = c(35, 45, 55),
    row.names = c("pathway1", "pathway2", "pathway3")
  )
  
  metadata <- data.frame(
    sample = paste0("sample", 1:4),
    group = factor(c("control", "control", "treatment", "treatment"))
  )
  
  # Test that default behavior is unchanged
  result_default <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2"
  )
  
  result_explicit_false <- pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2",
    include_abundance_stats = FALSE
  )
  
  # Results should be identical
  expect_equal(colnames(result_default), colnames(result_explicit_false))
  expect_equal(nrow(result_default), nrow(result_explicit_false))
  
  # Should have basic columns but not abundance stats
  basic_cols <- c("feature", "method", "group1", "group2", "p_values", "p_adjust", "adj_method")
  expect_true(all(basic_cols %in% colnames(result_default)))
  
  abundance_cols <- c("mean_rel_abundance_group1", "sd_rel_abundance_group1",
                     "mean_rel_abundance_group2", "sd_rel_abundance_group2", 
                     "log2_fold_change")
  expect_false(any(abundance_cols %in% colnames(result_default)))
})
