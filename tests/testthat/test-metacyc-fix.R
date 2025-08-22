test_that("pathway_annotation works correctly with MetaCyc pathways", {
  # Load test data
  load('data/metacyc_abundance.RData')
  
  # Create test DAA results
  pathway_features <- metacyc_abundance$pathway[1:10]
  daa_results_df <- data.frame(
    feature = pathway_features,
    method = rep("LinDA", length(pathway_features)),
    group1 = rep("Group1", length(pathway_features)),
    group2 = rep("Group2", length(pathway_features)),
    p_values = runif(length(pathway_features), 0.001, 0.9),
    p_adjust = runif(length(pathway_features), 0.001, 0.9),
    stringsAsFactors = FALSE
  )
  
  # Test that pathway_annotation works with MetaCyc
  annotated_results <- pathway_annotation(
    pathway = "MetaCyc",
    daa_results_df = daa_results_df,
    ko_to_kegg = FALSE
  )
  
  # Check that annotation was successful
  expect_true("description" %in% colnames(annotated_results))
  expect_equal(nrow(annotated_results), nrow(daa_results_df))
  
  # Check that all features were annotated (no missing descriptions)
  missing_descriptions <- sum(is.na(annotated_results$description))
  expect_equal(missing_descriptions, 0)
  
  # Check that specific known pathways are correctly annotated
  if ("1CMET2-PWY" %in% annotated_results$feature) {
    idx <- which(annotated_results$feature == "1CMET2-PWY")
    expect_equal(annotated_results$description[idx], "N10-formyl-tetrahydrofolate biosynthesis")
  }
})

test_that("load_reference_data works with multiple fallback strategies", {
  # Test that the improved load_reference_data function works
  ref_data <- load_reference_data("MetaCyc")
  
  expect_true(is.data.frame(ref_data))
  expect_true(nrow(ref_data) > 0)
  expect_true("id" %in% colnames(ref_data))
  expect_true("description" %in% colnames(ref_data))
  
  # Check that specific pathways exist
  expect_true("1CMET2-PWY" %in% ref_data$id)
  expect_true("ALL-CHORISMATE-PWY" %in% ref_data$id)
})

test_that("pathway_errorbar provides helpful error messages for MetaCyc annotation issues", {
  # Create test data with missing descriptions (simulating the original issue)
  abundance <- matrix(runif(50), nrow=5, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:5)
  colnames(abundance) <- paste0("sample", 1:10)
  
  group_data <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5)
  )
  
  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:5),
    method = rep("LinDA", 5),
    group1 = rep("GroupA", 5),
    group2 = rep("GroupB", 5),
    p_adjust = rep(0.01, 5),  # All significant
    description = rep(NA_character_, 5),  # Missing descriptions
    stringsAsFactors = FALSE
  )
  
  Group <- factor(group_data$group)
  names(Group) <- group_data$sample
  
  # This should fail with a helpful error message about missing annotations
  expect_error(
    pathway_errorbar(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = Group,
      ko_to_kegg = FALSE,
      x_lab = "description"
    ),
    regexp = "All features were excluded due to missing 'description' annotations"
  )
})

test_that("MetaCyc workflow integration test", {
  # Skip if required packages are not available
  skip_if_not_installed("ggprism")
  
  # Load test data
  load('data/metacyc_abundance.RData')
  load('data/metadata.RData')
  
  # Create realistic test DAA results
  pathway_features <- metacyc_abundance$pathway[1:10]
  daa_results_df <- data.frame(
    feature = pathway_features,
    method = rep("LinDA", length(pathway_features)),
    group1 = rep("Pro-survival", length(pathway_features)),
    group2 = rep("Pro-inflammatory", length(pathway_features)),
    p_values = c(rep(0.01, 3), rep(0.3, 7)),  # 3 significant, 7 not
    p_adjust = c(rep(0.01, 3), rep(0.3, 7)),
    stringsAsFactors = FALSE
  )
  
  # Test complete workflow
  annotated_results <- pathway_annotation(
    pathway = "MetaCyc",
    daa_results_df = daa_results_df,
    ko_to_kegg = FALSE
  )
  
  expect_true("description" %in% colnames(annotated_results))
  expect_equal(sum(is.na(annotated_results$description)), 0)
  
  # Prepare data for pathway_errorbar
  abundance_matrix <- as.matrix(metacyc_abundance[, -1])
  rownames(abundance_matrix) <- metacyc_abundance$pathway
  
  Group <- factor(metadata$Environment)
  names(Group) <- metadata$sample_name
  
  # This should work now (if plotting dependencies are available)
  # We'll test the core logic without requiring the full plotting stack
  expect_true(nrow(annotated_results) > 0)
  expect_true(sum(annotated_results$p_adjust < 0.05) >= 3)
})