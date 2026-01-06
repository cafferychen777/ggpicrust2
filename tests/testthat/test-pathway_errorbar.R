library(ggplot2)
library(patchwork)
library(dplyr)

# Note: ggprism is checked inside each test_that block

test_that("pathway_errorbar basic functionality works", {
  skip_if_not_installed("ggprism")
  library(ggprism)

  # Setup test data
  set.seed(123)
  abundance <- matrix(runif(100), nrow=10, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:10)
  colnames(abundance) <- paste0("sample", 1:10)

  # 创建一个数据框来存储分组信息
  metadata <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5)
  )

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:10),
    pathway_name = paste0("Pathway ", 1:10),
    description = paste0("Description ", 1:10),
    p_adjust = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1),
    method = rep("ALDEx2_Welch's t test", 10),
    group1 = rep("GroupA", 10),
    group2 = rep("GroupB", 10),
    stringsAsFactors = FALSE
  )

  Group <- metadata$group
  names(Group) <- metadata$sample

  # 添加数据检查
  testthat::expect_true(all(colnames(abundance) == names(Group)))
  testthat::expect_true(all(rownames(abundance) == daa_results_df$feature))

  p <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    select = paste0("pathway", 1:5),
    x_lab = "pathway_name"
  )

  expect_s3_class(p, "patchwork")
})

test_that("pathway_errorbar pathway_names_text_size parameter works", {
  skip_if_not_installed("ggprism")

  # Setup test data
  set.seed(123)
  abundance <- matrix(runif(100), nrow=10, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:10)
  colnames(abundance) <- paste0("sample", 1:10)

  metadata <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5)
  )

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:10),
    pathway_name = paste0("Pathway ", 1:10),
    description = paste0("Description ", 1:10),
    p_adjust = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1),
    method = rep("ALDEx2_Welch's t test", 10),
    group1 = rep("GroupA", 10),
    group2 = rep("GroupB", 10),
    stringsAsFactors = FALSE
  )

  Group <- metadata$group
  names(Group) <- metadata$sample

  # Test with auto text size
  p1 <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    select = paste0("pathway", 1:5),
    x_lab = "pathway_name",
    pathway_names_text_size = "auto"
  )

  expect_s3_class(p1, "patchwork")

  # Test with custom text size
  p2 <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    select = paste0("pathway", 1:5),
    x_lab = "pathway_name",
    pathway_names_text_size = 12
  )

  expect_s3_class(p2, "patchwork")
})

test_that("pathway_errorbar handles invalid inputs", {
  skip_if_not_installed("ggprism")

  # Test missing annotations
  abundance <- matrix(runif(20), nrow=2, ncol=10)
  rownames(abundance) <- c("pathway1", "pathway2")
  colnames(abundance) <- paste0("sample", 1:10)

  # 创建分组数据
  group_data <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5),
    stringsAsFactors = FALSE
  )

  daa_results_df_missing <- data.frame(
    feature = c("pathway1", "pathway2"),
    pathway_name = c(NA, "Pathway 2"),
    p_adjust = c(0.01, 0.02),
    method = rep("ALDEx2_Welch's t test", 2),
    group1 = rep("GroupA", 2),
    group2 = rep("GroupB", 2),
    stringsAsFactors = FALSE
  )

  # 创建 Group 向量，确保是因子类型
  Group <- factor(group_data$group, levels = c("GroupA", "GroupB"))
  names(Group) <- group_data$sample

  # 验证数据结构
  testthat::expect_equal(length(Group), ncol(abundance))
  testthat::expect_equal(names(Group), colnames(abundance))
  testthat::expect_equal(rownames(abundance), daa_results_df_missing$feature)

  # 使用 tryCatch 来捕获消息
  result <- tryCatch({
    pathway_errorbar(
      abundance = abundance,
      daa_results_df = daa_results_df_missing,
      Group = Group,
      x_lab = "pathway_name"
    )
  }, message = function(m) m)

  expect_true(grepl("The following pathways are missing annotations", result$message))
})

test_that("pathway_errorbar ordering works correctly", {
  skip_if_not_installed("ggprism")

  abundance <- matrix(runif(50), nrow=5, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:5)
  colnames(abundance) <- paste0("sample", 1:10)

  group_data <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5),
    stringsAsFactors = FALSE
  )

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:5),
    pathway_name = paste0("Pathway ", 1:5),
    pathway_class = c("Class1", "Class1", "Class2", "Class2", "Class3"),
    p_adjust = c(0.04, 0.01, 0.03, 0.02, 0.05),
    method = rep("ALDEx2_Welch's t test", 5),
    group1 = rep("GroupA", 5),
    group2 = rep("GroupB", 5),
    stringsAsFactors = FALSE
  )

  Group <- factor(group_data$group)
  names(Group) <- group_data$sample

  p1 <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    order = "p_values",
    x_lab = "pathway_name"
  )
  expect_s3_class(p1, "patchwork")
})

test_that("pathway_errorbar handles too many features", {
  skip_if_not_installed("ggprism")

  n_features <- 31  # 超过30个特征
  abundance <- matrix(runif(n_features * 10), nrow=n_features, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:n_features)
  colnames(abundance) <- paste0("sample", 1:10)

  # 创建分组数据

test_that("pathway_errorbar works with ko_to_kegg=TRUE and pathway_class (annotation_custom fix regression test)", {
  skip_if_not_installed("ggprism")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("GGally")
  
  # This test specifically verifies the fix for issue #184
  # The bug occurred when ko_to_kegg=TRUE with pathway_class annotations
  # Error: "no applicable method for 'rescale' applied to an object of class 'c('simpleUnit', 'unit', 'unit_v2')'"
  
  set.seed(123)
  abundance <- matrix(runif(60), nrow=6, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:6)
  colnames(abundance) <- paste0("sample", 1:10)
  
  metadata <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5)
  )
  
  # Create DAA results with pathway_class column (required for ko_to_kegg=TRUE scenario)
  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:6),
    pathway_name = paste0("Pathway ", 1:6),
    pathway_class = rep(c("Metabolism", "Genetic Information", "Environmental"), each=2),
    description = paste0("Description ", 1:6),
    p_adjust = c(0.001, 0.01, 0.02, 0.03, 0.04, 0.049),
    method = rep("ALDEx2_Welch's t test", 6),
    group1 = rep("GroupA", 6),
    group2 = rep("GroupB", 6),
    stringsAsFactors = FALSE
  )
  
  Group <- factor(metadata$group)
  names(Group) <- metadata$sample
  
  # Test with ko_to_kegg=TRUE which triggers pathway class background annotations
  # This would fail in ggplot2 4.0.0+ before the fix
  p <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    ko_to_kegg = TRUE,  # Triggers pathway class annotation code
    p_values_threshold = 0.05,
    order = "pathway_class",
    p_value_bar = TRUE,  # Creates 4-panel patchwork
    x_lab = "pathway_name"
  )
  
  # Should return a patchwork object with 4 panels
  expect_s3_class(p, "patchwork")
  
  # Test that plot can be rendered without the unit rescale error
  expect_no_error({
    # Render to null device to avoid display
    pdf(NULL)
    print(p)
    dev.off()
  })
  
  # Test with p_value_bar=FALSE (3-panel patchwork)
  p2 <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    p_values_threshold = 0.05,
    order = "pathway_class",
    p_value_bar = FALSE,  # Creates 3-panel patchwork
    x_lab = "pathway_name"
  )
  
  expect_s3_class(p2, "patchwork")
  
  # Verify rendering works
  expect_no_error({
    pdf(NULL)
    print(p2)
    dev.off()
  })
})

test_that("pathway_errorbar handles too many features", {
  skip_if_not_installed("ggprism")

  n_features <- 31  # 超过30个特征
  abundance <- matrix(runif(n_features * 10), nrow=n_features, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:n_features)
  colnames(abundance) <- paste0("sample", 1:10)

  # 创建分组数据
  group_data <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5),
    stringsAsFactors = FALSE
  )

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:n_features),
    pathway_name = paste0("Pathway ", 1:n_features),
    p_adjust = rep(0.01, n_features),  # 所有p值都显著
    method = rep("ALDEx2_Welch's t test", n_features),
    group1 = rep("GroupA", n_features),
    group2 = rep("GroupB", n_features),
    stringsAsFactors = FALSE
  )

  # 创建 Group 向量，确保是因子类型
  Group <- factor(group_data$group, levels = c("GroupA", "GroupB"))
  names(Group) <- group_data$sample

  # 验证数据结构
  testthat::expect_equal(length(Group), ncol(abundance))
  testthat::expect_equal(names(Group), colnames(abundance))
  testthat::expect_equal(rownames(abundance), daa_results_df$feature)

  # Function now warns instead of errors for too many features
  expect_warning(
    pathway_errorbar(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = Group,
      x_lab = "pathway_name"
    ),
    regexp = "The number of features with statistical significance exceeds 30"
  )
})

test_that("pathway_errorbar handles custom colors correctly", {
  skip_if_not_installed("ggprism")

  # Setup basic test data
  abundance <- matrix(runif(50), nrow=5, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:5)
  colnames(abundance) <- paste0("sample", 1:10)

  group_data <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5)
  )

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:5),
    pathway_name = paste0("Pathway ", 1:5),
    p_adjust = rep(0.01, 5),
    method = rep("ALDEx2_Welch's t test", 5),
    group1 = rep("GroupA", 5),
    group2 = rep("GroupB", 5)
  )

  Group <- factor(group_data$group)
  names(Group) <- group_data$sample

  # Test with custom colors
  custom_colors <- c("#FF0000", "#0000FF")
  p <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    colors = custom_colors,
    x_lab = "pathway_name"
  )

  expect_s3_class(p, "patchwork")
})

test_that("pathway_errorbar handles different ordering options", {
  skip_if_not_installed("ggprism")

  # Setup test data
  abundance <- matrix(runif(50), nrow=5, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:5)
  colnames(abundance) <- paste0("sample", 1:10)

  group_data <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5)
  )

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:5),
    pathway_name = paste0("Pathway ", 1:5),
    pathway_class = c("Class1", "Class1", "Class2", "Class2", "Class3"),
    p_adjust = c(0.04, 0.01, 0.03, 0.02, 0.05),
    method = rep("ALDEx2_Welch's t test", 5),
    group1 = rep("GroupA", 5),
    group2 = rep("GroupB", 5)
  )

  Group <- factor(group_data$group)
  names(Group) <- group_data$sample

  # Test different ordering options
  orders <- c("p_values", "name", "group", "pathway_class")
  for(order_type in orders) {
    p <- pathway_errorbar(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = Group,
      order = order_type,
      x_lab = "pathway_name"
    )
    expect_s3_class(p, "patchwork")
  }

  # Test invalid order type
  expect_error(
    pathway_errorbar(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = Group,
      order = "invalid_order",
      x_lab = "pathway_name"
    )
  )
})

test_that("pathway_errorbar handles p_value_bar parameter correctly", {
  skip_if_not_installed("ggprism")

  # Setup basic test data
  abundance <- matrix(runif(50), nrow=5, ncol=10)
  rownames(abundance) <- paste0("pathway", 1:5)
  colnames(abundance) <- paste0("sample", 1:10)

  group_data <- data.frame(
    sample = paste0("sample", 1:10),
    group = rep(c("GroupA", "GroupB"), each=5)
  )

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:5),
    pathway_name = paste0("Pathway ", 1:5),
    p_adjust = rep(0.01, 5),
    method = rep("ALDEx2_Welch's t test", 5),
    group1 = rep("GroupA", 5),
    group2 = rep("GroupB", 5)
  )

  Group <- factor(group_data$group)
  names(Group) <- group_data$sample

  # Test with p_value_bar = TRUE and FALSE
  p1 <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    p_value_bar = TRUE,
    x_lab = "pathway_name"
  )

  p2 <- pathway_errorbar(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = Group,
    p_value_bar = FALSE,
    x_lab = "pathway_name"
  )

  expect_s3_class(p1, "patchwork")
  expect_s3_class(p2, "patchwork")
  expect_false(identical(p1, p2))
})
