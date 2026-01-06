# Test GO pathway support in GSEA functionality

test_that("ko_to_go_reference dataset is properly structured", {
  expect_no_error({
    data("ko_to_go_reference", package = "ggpicrust2")
  })

  expect_s3_class(ko_to_go_reference, "data.frame")
  expect_true(all(c("go_id", "go_name", "category", "ko_members") %in% colnames(ko_to_go_reference)))

  # Format validation
  expect_true(all(grepl("^GO:\\d{7}$", ko_to_go_reference$go_id)))
  expect_true(all(ko_to_go_reference$category %in% c("BP", "MF", "CC")))
  expect_true(length(unique(ko_to_go_reference$category)) == 3)
  expect_true(nrow(ko_to_go_reference) >= 100)
})

test_that("prepare_gene_sets works with GO pathway_type", {
  gene_sets_bp <- prepare_gene_sets("GO", go_category = "BP")
  gene_sets_mf <- prepare_gene_sets("GO", go_category = "MF")
  gene_sets_all <- prepare_gene_sets("GO", go_category = "all")

  expect_type(gene_sets_bp, "list")
  expect_true(length(gene_sets_bp) > 0)
  expect_true(all(grepl("^GO:\\d{7}$", names(gene_sets_bp))))

  # Categories should not overlap
  expect_true(length(intersect(names(gene_sets_bp), names(gene_sets_mf))) == 0)

  # All should include individual categories
  expect_true(all(names(gene_sets_bp) %in% names(gene_sets_all)))
})

test_that("pathway_gsea works with GO pathway_type", {
  data("ko_abundance", package = "ggpicrust2")
  data("metadata", package = "ggpicrust2")

  abundance_data <- ko_abundance %>%
    tibble::column_to_rownames("#NAME")
  abundance_data <- abundance_data[1:100, ]

  meta_df <- as.data.frame(metadata)
  rownames(meta_df) <- meta_df$sample_name

  expect_no_error({
    gsea_results <- pathway_gsea(
      abundance = abundance_data,
      metadata = meta_df,
      group = "Environment",
      pathway_type = "GO",
      go_category = "BP",
      method = "fgsea",
      nperm = 100,
      min_size = 2,
      max_size = 500
    )
  })
})
