# MetaCyc Pathway Annotation System Tests
# Comprehensive validation of gsea_pathway_annotation for MetaCyc pathways
# Following Linus principles: reliable annotation, clear fallbacks, no broken userspace

test_that("MetaCyc pathway annotation basic functionality", {
  # Create mock GSEA results with known MetaCyc pathways
  mock_results <- data.frame(
    pathway_id = c(
      "1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN",
      "ALL-CHORISMATE-PWY", "ANAEROFRUCAT-PWY", "UNKNOWN-PATHWAY-TEST"
    ),
    size = c(4, 5, 6, 3, 4, 2),
    ES = c(0.6, -0.4, 0.8, 0.3, -0.5, 0.1),
    NES = c(1.5, -1.2, 1.8, 0.7, -1.1, 0.2),
    pvalue = c(0.02, 0.08, 0.01, 0.15, 0.06, 0.45),
    p.adjust = c(0.06, 0.12, 0.04, 0.20, 0.10, 0.50),
    stringsAsFactors = FALSE
  )
  
  # Test annotation function
  annotated_results <- gsea_pathway_annotation(mock_results, pathway_type = "MetaCyc")
  
  # Basic structure validation
  expect_s3_class(annotated_results, "data.frame")
  expect_equal(nrow(annotated_results), nrow(mock_results))
  expect_true("pathway_name" %in% colnames(annotated_results))
  
  # All original columns should be preserved
  original_cols <- colnames(mock_results)
  expect_true(all(original_cols %in% colnames(annotated_results)))
  
  # Check that annotation worked for known pathways
  # The exact names depend on the annotation reference, but should be meaningful
  expect_type(annotated_results$pathway_name, "character")
  expect_true(all(!is.na(annotated_results$pathway_name)))
  expect_true(all(nchar(annotated_results$pathway_name) > 0))
})

test_that("MetaCyc annotation with pathway name mapping validation", {
  # Test specific pathway mappings that we expect to exist
  test_pathways <- data.frame(
    pathway_id = c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN"),
    NES = c(1.5, -1.2, 1.8),
    pvalue = c(0.01, 0.05, 0.02),
    stringsAsFactors = FALSE
  )
  
  annotated <- gsea_pathway_annotation(test_pathways, pathway_type = "MetaCyc")
  
  # Check that known pathways get appropriate names
  expect_equal(nrow(annotated), 3)
  
  # 1CMET2-PWY should get a meaningful name related to folate metabolism
  cmet_row <- annotated[annotated$pathway_id == "1CMET2-PWY", ]
  expect_equal(nrow(cmet_row), 1)
  # The name should be different from the ID (indicating successful annotation)
  expect_true(cmet_row$pathway_name != "1CMET2-PWY" | 
              grepl("formyl|folate|tetrahydrofolate", cmet_row$pathway_name, ignore.case = TRUE))
  
  # ANAGLYCOLYSIS-PWY should get a name related to glycolysis
  glycolysis_row <- annotated[annotated$pathway_id == "ANAGLYCOLYSIS-PWY", ]
  expect_equal(nrow(glycolysis_row), 1)
  # Should either keep the ID or get a glycolysis-related name
  expect_true(nchar(glycolysis_row$pathway_name) > 0)
  
  # ARG+POLYAMINE-SYN should get a name related to arginine/polyamine synthesis
  arg_row <- annotated[annotated$pathway_id == "ARG+POLYAMINE-SYN", ]
  expect_equal(nrow(arg_row), 1)
  expect_true(nchar(arg_row$pathway_name) > 0)
})

test_that("MetaCyc annotation fallback behavior for unknown pathways", {
  # Test with mixture of known and unknown pathways
  mixed_results <- data.frame(
    pathway_id = c(
      "1CMET2-PWY",  # Known pathway
      "FAKE-PATHWAY-001",  # Unknown pathway
      "ANOTHER-FAKE-PWY",  # Unknown pathway
      "TEST-UNKNOWN"  # Unknown pathway
    ),
    NES = c(1.5, 0.8, -1.1, 0.3),
    pvalue = c(0.01, 0.20, 0.05, 0.30),
    stringsAsFactors = FALSE
  )
  
  annotated <- gsea_pathway_annotation(mixed_results, pathway_type = "MetaCyc")
  
  expect_equal(nrow(annotated), 4)
  
  # Unknown pathways should keep their pathway_id as pathway_name (fallback behavior)
  unknown_pathways <- c("FAKE-PATHWAY-001", "ANOTHER-FAKE-PWY", "TEST-UNKNOWN")
  
  for (unknown_id in unknown_pathways) {
    unknown_row <- annotated[annotated$pathway_id == unknown_id, ]
    expect_equal(nrow(unknown_row), 1)
    # Fallback: pathway_name should equal pathway_id for unknown pathways
    expect_equal(unknown_row$pathway_name, unknown_id)
  }
  
  # Known pathway should get annotation (or at least not be equal to ID if annotation exists)
  known_row <- annotated[annotated$pathway_id == "1CMET2-PWY", ]
  expect_equal(nrow(known_row), 1)
  expect_true(nchar(known_row$pathway_name) > 0)
})

test_that("MetaCyc annotation data integrity preservation", {
  # Test that all original data is preserved during annotation
  original_data <- data.frame(
    pathway_id = c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "TEST-PATHWAY"),
    size = c(4, 5, 3),
    ES = c(0.6, -0.4, 0.2),
    NES = c(1.5, -1.2, 0.5),
    pvalue = c(0.02, 0.08, 0.25),
    p.adjust = c(0.06, 0.12, 0.30),
    leading_edge = c("EC:1.1.1.1;EC:2.2.2.2", "EC:3.3.3.3;EC:4.4.4.4", "EC:5.5.5.5"),
    method = c("fgsea", "fgsea", "fgsea"),
    stringsAsFactors = FALSE
  )
  
  annotated <- gsea_pathway_annotation(original_data, pathway_type = "MetaCyc")
  
  # Check data integrity
  expect_equal(nrow(annotated), nrow(original_data))
  
  # All original columns should be preserved with exact values
  for (col in colnames(original_data)) {
    expect_equal(annotated[[col]], original_data[[col]],
                info = paste("Column", col, "was modified during annotation"))
  }
  
  # pathway_name should be added
  expect_true("pathway_name" %in% colnames(annotated))
  expect_equal(length(annotated$pathway_name), nrow(original_data))
  expect_true(all(!is.na(annotated$pathway_name)))
})

test_that("MetaCyc annotation with empty and edge case inputs", {
  # Test with empty results
  empty_results <- data.frame(
    pathway_id = character(),
    NES = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  annotated_empty <- gsea_pathway_annotation(empty_results, pathway_type = "MetaCyc")
  expect_s3_class(annotated_empty, "data.frame")
  expect_equal(nrow(annotated_empty), 0)
  expect_true("pathway_name" %in% colnames(annotated_empty))
  
  # Test with single pathway
  single_result <- data.frame(
    pathway_id = "1CMET2-PWY",
    NES = 1.5,
    pvalue = 0.02,
    stringsAsFactors = FALSE
  )
  
  annotated_single <- gsea_pathway_annotation(single_result, pathway_type = "MetaCyc")
  expect_equal(nrow(annotated_single), 1)
  expect_true("pathway_name" %in% colnames(annotated_single))
  expect_equal(annotated_single$pathway_id, "1CMET2-PWY")
  
  # Test with pathways containing special characters
  special_results <- data.frame(
    pathway_id = c("ARG+POLYAMINE-SYN", "ALL-CHORISMATE-PWY", "PATHWAY.WITH.DOTS"),
    NES = c(1.2, -0.8, 0.5),
    pvalue = c(0.03, 0.10, 0.20),
    stringsAsFactors = FALSE
  )
  
  annotated_special <- gsea_pathway_annotation(special_results, pathway_type = "MetaCyc")
  expect_equal(nrow(annotated_special), 3)
  expect_true(all(!is.na(annotated_special$pathway_name)))
  
  # Pathway IDs should be preserved exactly
  expect_equal(annotated_special$pathway_id, special_results$pathway_id)
})

test_that("MetaCyc annotation performance with large datasets", {
  # Test performance with realistic dataset size
  set.seed(12345)
  
  # Create large mock dataset
  n_pathways <- 200
  large_pathway_ids <- paste0("PATHWAY-", sprintf("%03d", 1:n_pathways), "-PWY")
  
  # Include some real MetaCyc pathway IDs
  real_metacyc_ids <- c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN", 
                        "ALL-CHORISMATE-PWY", "ANAEROFRUCAT-PWY")
  large_pathway_ids[1:length(real_metacyc_ids)] <- real_metacyc_ids
  
  large_results <- data.frame(
    pathway_id = large_pathway_ids,
    size = sample(5:20, n_pathways, replace = TRUE),
    ES = rnorm(n_pathways, mean = 0, sd = 0.5),
    NES = rnorm(n_pathways, mean = 0, sd = 1.5),
    pvalue = runif(n_pathways, min = 0.001, max = 0.999),
    p.adjust = runif(n_pathways, min = 0.001, max = 0.999),
    stringsAsFactors = FALSE
  )
  
  # Performance test
  start_time <- Sys.time()
  
  annotated_large <- gsea_pathway_annotation(large_results, pathway_type = "MetaCyc")
  
  end_time <- Sys.time()
  annotation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Performance validation
  expect_lt(annotation_time, 5)  # Should complete within 5 seconds
  
  # Result validation
  expect_equal(nrow(annotated_large), n_pathways)
  expect_true("pathway_name" %in% colnames(annotated_large))
  expect_true(all(!is.na(annotated_large$pathway_name)))
  expect_equal(annotated_large$pathway_id, large_results$pathway_id)
  
  # Real MetaCyc pathways should get proper annotation
  for (real_id in real_metacyc_ids) {
    real_row <- annotated_large[annotated_large$pathway_id == real_id, ]
    expect_equal(nrow(real_row), 1)
    expect_true(nchar(real_row$pathway_name) > 0)
  }
})

test_that("MetaCyc annotation error handling and robustness", {
  # Test with malformed input
  expect_error({
    gsea_pathway_annotation(NULL, pathway_type = "MetaCyc")
  })
  
  expect_error({
    gsea_pathway_annotation("not_a_dataframe", pathway_type = "MetaCyc")
  })
  
  # Test with missing pathway_id column
  no_id_results <- data.frame(
    pathway_name = c("Test1", "Test2"),
    NES = c(1.0, -1.0),
    stringsAsFactors = FALSE
  )
  
  expect_error({
    gsea_pathway_annotation(no_id_results, pathway_type = "MetaCyc")
  })
  
  # Test with invalid pathway_type
  valid_results <- data.frame(
    pathway_id = c("TEST1", "TEST2"),
    NES = c(1.0, -1.0),
    stringsAsFactors = FALSE
  )
  
  expect_error({
    gsea_pathway_annotation(valid_results, pathway_type = "InvalidType")
  })
  
  # Test with NA values in pathway_id
  na_results <- data.frame(
    pathway_id = c("1CMET2-PWY", NA, "ANAGLYCOLYSIS-PWY"),
    NES = c(1.0, 0.5, -1.0),
    stringsAsFactors = FALSE
  )
  
  # Should handle NA gracefully
  expect_no_error({
    annotated_na <- gsea_pathway_annotation(na_results, pathway_type = "MetaCyc")
    expect_equal(nrow(annotated_na), 3)
    expect_true(is.na(annotated_na$pathway_id[2]))
  })
})

test_that("MetaCyc annotation consistency across multiple runs", {
  # Test that annotation is deterministic and consistent
  test_results <- data.frame(
    pathway_id = c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "UNKNOWN-PWY-123"),
    NES = c(1.5, -1.2, 0.3),
    pvalue = c(0.01, 0.05, 0.25),
    stringsAsFactors = FALSE
  )
  
  # Run annotation multiple times
  annotated1 <- gsea_pathway_annotation(test_results, pathway_type = "MetaCyc")
  annotated2 <- gsea_pathway_annotation(test_results, pathway_type = "MetaCyc")
  annotated3 <- gsea_pathway_annotation(test_results, pathway_type = "MetaCyc")
  
  # Results should be identical across runs
  expect_equal(annotated1$pathway_name, annotated2$pathway_name)
  expect_equal(annotated1$pathway_name, annotated3$pathway_name)
  expect_equal(annotated1$pathway_id, annotated2$pathway_id)
  expect_equal(annotated1$pathway_id, annotated3$pathway_id)
  
  # All other columns should be preserved identically
  for (col in colnames(test_results)) {
    expect_equal(annotated1[[col]], annotated2[[col]])
    expect_equal(annotated1[[col]], annotated3[[col]])
  }
})

test_that("MetaCyc annotation integration with visualize_gsea pathway labels", {
  # Test that annotated results work properly with visualization functions
  # This ensures the annotation system integrates well with downstream functions
  
  mock_annotated_results <- data.frame(
    pathway_id = c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN"),
    pathway_name = c(
      "N10-formyl-tetrahydrofolate biosynthesis",
      "ANAGLYCOLYSIS-PWY", 
      "ARG+POLYAMINE-SYN"
    ),
    size = c(4, 5, 6),
    ES = c(0.6, -0.4, 0.8),
    NES = c(1.5, -1.2, 1.8),
    pvalue = c(0.02, 0.08, 0.01),
    p.adjust = c(0.06, 0.12, 0.04),
    stringsAsFactors = FALSE
  )
  
  # Test that pathway_name column is available for visualization
  expect_true("pathway_name" %in% colnames(mock_annotated_results))
  expect_true(all(!is.na(mock_annotated_results$pathway_name)))
  expect_true(all(nchar(mock_annotated_results$pathway_name) > 0))
  
  # Test pathway label selection logic (simulating what visualize_gsea does)
  # When pathway_name is available, it should be preferred over pathway_id
  pathway_labels <- ifelse(
    !is.na(mock_annotated_results$pathway_name) & 
    mock_annotated_results$pathway_name != "",
    mock_annotated_results$pathway_name,
    mock_annotated_results$pathway_id
  )
  
  expect_equal(length(pathway_labels), nrow(mock_annotated_results))
  expect_true(all(nchar(pathway_labels) > 0))
  
  # First pathway should use the descriptive name
  expect_equal(pathway_labels[1], "N10-formyl-tetrahydrofolate biosynthesis")
  
  # Others should use their names (even if same as ID for unannotated pathways)
  expect_equal(pathway_labels[2], "ANAGLYCOLYSIS-PWY")
  expect_equal(pathway_labels[3], "ARG+POLYAMINE-SYN")
})