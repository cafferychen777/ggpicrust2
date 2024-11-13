# Create test data
test_that("pathway_daa works with basic inputs", {
  # Setup test data with more samples per group
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(20, 30, 40),
    sample3 = c(15, 25, 35),
    sample4 = c(30, 40, 50),
    sample5 = c(35, 45, 55),
    sample6 = c(25, 35, 45),
    row.names = c("pathway1", "pathway2", "pathway3")
  )

  metadata <- tibble::tibble(
    sample = paste0("sample", 1:6),
    group = c("control", "control", "control", "treatment", "treatment", "treatment")
  )

  # Test ALDEx2 method
  result <- pathway_daa(abundance, metadata, "group", daa_method = "ALDEx2")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("feature", "method", "p_values", "adj_method", "p_adjust") %in% colnames(result)))
  expect_equal(nrow(result), nrow(abundance))
  expect_equal(result$method[1], "ALDEx2")
  expect_true(all(!is.na(result$p_values)))
  expect_true(all(result$p_values >= 0 & result$p_values <= 1))
})

test_that("pathway_daa validates inputs correctly", {
  # 创建测试数据
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(20, 30, 40),
    sample3 = c(30, 40, 50),
    sample4 = c(40, 50, 60),  # 添加第四个样本
    row.names = c("pathway1", "pathway2", "pathway3")
  )

  # 创建不匹配的元数据
  metadata_mismatch <- data.frame(
    wrong_sample = c("wrong1", "wrong2", "wrong3", "wrong4"),
    group = c("control", "control", "treatment", "treatment")
  )

  # 测试样本名不匹配的情况
  expect_error(
    pathway_daa(abundance, metadata_mismatch, "group"),
    "No column in metadata matches the sample names in abundance data"
  )

  # 创建正确的元数据但组别数量不足
  metadata_single_group <- data.frame(
    sample = c("sample1", "sample2", "sample3", "sample4"),
    group = c("control", "control", "control", "control")
  )

  # 测试单一组别的情况
  expect_error(
    pathway_daa(abundance, metadata_single_group, "group"),
    "At least two groups are required"
  )

  # 测试样本数量不足的情况
  small_abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(20, 30, 40),
    sample3 = c(30, 40, 50),
    row.names = c("pathway1", "pathway2", "pathway3")
  )

  small_metadata <- data.frame(
    sample = c("sample1", "sample2", "sample3"),
    group = c("control", "control", "treatment")
  )

  expect_error(
    pathway_daa(small_abundance, small_metadata, "group"),
    "At least 4 samples are required for differential abundance analysis"
  )
})

test_that("pathway_daa methods produce expected results", {
  # 创建更大的测试数据集以避免 DESeq2 警告
  n_samples <- 10  # 增加样本数
  n_features <- 3

  # 创建模拟数据
  set.seed(123)
  abundance <- matrix(
    rpois(n_samples * n_features, lambda = 20),  # 使用泊松分布生成计数数据
    nrow = n_features,
    ncol = n_samples
  )
  rownames(abundance) <- paste0("pathway", 1:n_features)
  colnames(abundance) <- paste0("sample", 1:n_samples)
  abundance <- as.data.frame(abundance)

  # 创建平衡的元数据
  metadata <- data.frame(
    sample = paste0("sample", 1:n_samples),
    group = rep(c("control", "treatment"), each = n_samples/2)
  )

  # 测试所有方法
  for(method in c("ALDEx2", "DESeq2", "limma voom", "edgeR", "metagenomeSeq", "LinDA", "Maaslin2")) {
    result <- pathway_daa(abundance, metadata, "group", daa_method = method)

    # 确保结果是数据框
    expect_true(is.data.frame(result))

    # 验证必要的列存在
    expect_true(all(c("feature", "method", "p_values") %in% colnames(result)))

    # 验证行数正确
    expect_equal(nrow(result), n_features)

    # 验证 p 值在有效范围内
    expect_true(all(result$p_values >= 0 & result$p_values <= 1))
  }
})

test_that("pathway_daa handles sample selection correctly", {
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(20, 30, 40),
    sample3 = c(15, 25, 35),
    sample4 = c(30, 40, 50),
    row.names = c("pathway1", "pathway2", "pathway3")
  )

  metadata <- tibble::tibble(
    sample = paste0("sample", 1:4),
    group = c("control", "control", "treatment", "treatment")
  )

  # 测试有效的样本选择
  selected_samples <- c("sample1", "sample2", "sample3", "sample4")
  result <- pathway_daa(abundance, metadata, "group",
                       daa_method = "ALDEx2",
                       select = selected_samples)
  expect_equal(ncol(abundance[, selected_samples]), length(selected_samples))

  # 测试无效的样本选择
  expect_error(
    pathway_daa(abundance, metadata, "group",
                daa_method = "ALDEx2",
                select = c("sample1", "invalid_sample")),
    "Some selected samples are not present in the abundance data"
  )
})

test_that("pathway_daa handles multiple groups correctly", {
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(20, 30, 40),
    sample3 = c(30, 40, 50),
    sample4 = c(25, 35, 45),
    sample5 = c(15, 25, 35),
    sample6 = c(16, 26, 36),
    row.names = c("pathway1", "pathway2", "pathway3")
  )

  metadata <- tibble::tibble(
    sample = paste0("sample", 1:6),
    group = factor(c("control", "control", "treatment", "treatment", "other", "other"))
  )

  # 测试带参考组的多组比较
  suppressWarnings({
    result <- pathway_daa(abundance, metadata, "group",
                         daa_method = "limma voom",
                         reference = "control")
  })

  expect_s3_class(result, "data.frame")
  expect_true(all(!is.na(result$p_values)))
})

test_that("pathway_daa handles p-value adjustment correctly", {
  abundance <- data.frame(
    sample1 = c(10, 20, 30),
    sample2 = c(20, 30, 40),
    sample3 = c(30, 40, 50),
    sample4 = c(25, 35, 45),
    row.names = c("pathway1", "pathway2", "pathway3")
  )

  metadata <- tibble::tibble(
    sample = paste0("sample", 1:4),
    group = c("control", "control", "treatment", "treatment")
  )

  # 测试不同的p值调整方法
  methods <- c("BH", "holm", "bonferroni", "hochberg", "fdr")
  for(method in methods) {
    result <- pathway_daa(abundance, metadata, "group",
                         daa_method = "ALDEx2",
                         p.adjust = method)
    expect_true(all(result$adj_method == method))
    expect_true(all(!is.na(result$p_adjust)))
    expect_true(all(result$p_adjust >= 0 & result$p_adjust <= 1))
  }
})
