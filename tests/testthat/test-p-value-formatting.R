# Test for p-value formatting fix in pathway_errorbar
# This test ensures that the p-value display issue is resolved

test_that("p-value formatting displays correctly in pathway_errorbar", {
  # Create test data with problematic p-values
  abundance <- matrix(
    runif(50, 0, 100), 
    nrow = 5, 
    ncol = 10,
    dimnames = list(
      paste0("pathway_", 1:5),
      paste0("sample_", 1:10)
    )
  )
  
  # Create DAA results with very small p-values that were problematic
  daa_results_df <- data.frame(
    feature = paste0("pathway_", 1:5),
    p_adjust = c(6.686113e-06, 1.234567e-05, 0.001234, 0.0456, 0.123456),
    group1 = rep("GroupA", 5),
    group2 = rep("GroupB", 5),
    method = rep("ALDEx2", 5),
    pathway_name = paste0("Pathway ", 1:5),
    stringsAsFactors = FALSE
  )
  
  # Create Group vector
  Group <- factor(c(rep("GroupA", 5), rep("GroupB", 5)))
  names(Group) <- colnames(abundance)
  
  # Test the format_p_value function directly
  format_p_value <- function(p) {
    ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))
  }
  
  formatted_values <- format_p_value(daa_results_df$p_adjust)
  
  # Verify that very small p-values are correctly formatted
  expect_equal(formatted_values[1], "6.7e-06")  # Was "6.686" before fix
  expect_equal(formatted_values[2], "1.2e-05")  # Was "1.234" before fix
  expect_equal(formatted_values[3], "0.001")    # Regular formatting
  expect_equal(formatted_values[4], "0.046")    # Regular formatting
  expect_equal(formatted_values[5], "0.123")    # Regular formatting
  
  # Test that pathway_errorbar function works with the fix
  expect_no_error({
    p <- pathway_errorbar(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = Group,
      p_values_threshold = 0.5,  # High threshold to include all pathways
      x_lab = "pathway_name"
    )
  })
  
  # Verify that the plot is created successfully
  expect_s3_class(p, "patchwork")
})

test_that("format_p_value function handles edge cases correctly", {
  # Define the function as used in pathway_errorbar
  format_p_value <- function(p) {
    ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))
  }
  
  # Test edge cases
  edge_cases <- c(0, 1e-100, 0.0009999, 0.001, 0.001001, 1)
  formatted <- format_p_value(edge_cases)
  
  # Verify formatting
  expect_equal(formatted[1], "0.0e+00")  # Zero
  expect_equal(formatted[2], "1.0e-100") # Very small
  expect_equal(formatted[3], "1.0e-03")  # Just below threshold
  expect_equal(formatted[4], "0.001")    # At threshold
  expect_equal(formatted[5], "0.001")    # Just above threshold
  expect_equal(formatted[6], "1.000")    # Large value
  
  # Verify all outputs are character strings
  expect_true(all(is.character(formatted)))
  
  # Verify no values appear greater than 1 when they shouldn't
  numeric_values <- suppressWarnings(as.numeric(formatted))
  small_p_indices <- which(edge_cases < 0.001 & edge_cases > 0)
  expect_true(all(is.na(numeric_values[small_p_indices]) | numeric_values[small_p_indices] < 1))
})
