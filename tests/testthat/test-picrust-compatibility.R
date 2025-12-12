test_that("pathway_daa handles PICRUSt 2.6.2 compatibility issues", {
  # Create test data that simulates PICRUSt 2.6.2 output with many zero features
  abundance <- data.frame(
    sample1 = c(0, 0, 0, 10, 5),
    sample2 = c(0, 0, 0, 20, 8),
    sample3 = c(0, 0, 0, 15, 6),
    sample4 = c(0, 0, 0, 25, 10),
    row.names = c("pathway1", "pathway2", "pathway3", "pathway4", "pathway5")
  )

  metadata <- tibble::tibble(
    sample = paste0("sample", 1:4),
    group = c("control", "control", "treatment", "treatment")
  )

  # Test that the function handles mostly zero data gracefully
  expect_message(
    result <- pathway_daa(abundance, metadata, "group", daa_method = "ALDEx2"),
    "Filtering out.*features with zero abundance"
  )

  # Should still produce results for non-zero features
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
})

test_that("pathway_daa fails gracefully with all-zero data", {
  # Create test data with all zeros (extreme PICRUSt 2.6.2 issue)
  abundance <- data.frame(
    sample1 = c(0, 0, 0),
    sample2 = c(0, 0, 0),
    sample3 = c(0, 0, 0),
    sample4 = c(0, 0, 0),
    row.names = c("pathway1", "pathway2", "pathway3")
  )
  
  metadata <- tibble::tibble(
    sample = paste0("sample", 1:4),
    group = c("control", "control", "treatment", "treatment")
  )
  
  # Should fail with informative error message
  expect_error(
    pathway_daa(abundance, metadata, "group", daa_method = "ALDEx2"),
    "PICRUSt version compatibility issue"
  )
})

test_that("LinDA analysis handles empty data gracefully", {
  # Create test data that would cause LinDA to have limited features
  abundance <- data.frame(
    sample1 = c(0, 0, 10, 5, 8),
    sample2 = c(0, 0, 20, 8, 12),
    sample3 = c(0, 0, 15, 6, 10),
    sample4 = c(0, 0, 25, 10, 15),
    row.names = c("pathway1", "pathway2", "pathway3", "pathway4", "pathway5")
  )

  metadata <- tibble::tibble(
    sample = paste0("sample", 1:4),
    group = c("control", "control", "treatment", "treatment")
  )

  # Should handle the case where some features are filtered out
  expect_message(
    result <- pathway_daa(abundance, metadata, "group", daa_method = "LinDA"),
    "Filtering out.*features with zero abundance"
  )

  expect_true(is.data.frame(result))
})

test_that("ko2kegg_abundance handles PICRUSt 2.6.2 format", {
  # Create mock data that simulates PICRUSt 2.6.2 output issues
  # ko2kegg_abundance now requires 'function.' column format
  # Test with correct format (older PICRUSt 2.6.2 #NAME format is no longer supported)
  mock_ko_data <- data.frame(
    function. = c("K00001", "K00002", "K00003"),
    Sample1 = c(0, 0, 0),
    Sample2 = c(0, 0, 0),
    stringsAsFactors = FALSE
  )

  # Function should handle data with zero values
  expect_no_error({
    result <- ko2kegg_abundance(data = mock_ko_data)
  })
})
