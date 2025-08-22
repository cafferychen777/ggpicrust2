#!/usr/bin/env Rscript

# 简化版图例和标注测试
cat("简化版图例和标注测试开始...\n")

# 加载必要的库
suppressPackageStartupMessages({
  library(ggpicrust2)
  library(dplyr)
  library(ggplot2)
})

# 创建测试结果文件夹
test_output_dir <- "legend_annotation_test_results"
if (!dir.exists(test_output_dir)) {
  dir.create(test_output_dir, recursive = TRUE)
}

# 加载所需文件
source("R/legend_annotation_utils.R")
source("R/color_themes.R")
source("R/pathway_errorbar.R")

# 创建简单测试数据
set.seed(42)
n_pathways <- 10
n_samples <- 12

# 创建丰度矩阵
abundance_test <- matrix(
  abs(rnorm(n_pathways * n_samples, mean = 50, sd = 25)),
  nrow = n_pathways,
  ncol = n_samples,
  dimnames = list(
    paste0("PWY", sprintf("%03d", 1:n_pathways)),
    paste0("S", sprintf("%02d", 1:n_samples))
  )
)

# 创建分组 (简化为2组)
Group_test <- rep(c("Control", "Treatment"), each = n_samples/2)

# 创建DAA结果数据
daa_results_test <- data.frame(
  feature = rownames(abundance_test),
  p_adjust = runif(n_pathways, 0.001, 0.1),
  method = "ALDEx2_Welch's t test",
  group1 = "Control",
  group2 = "Treatment",
  log_2_fold_change = rnorm(n_pathways, 0, 2),
  description = paste("Description for", rownames(abundance_test)),
  pathway_class = sample(c("Metabolism", "Signaling", "Transport"), 
                        n_pathways, replace = TRUE),
  stringsAsFactors = FALSE
)

cat("✓ 简化测试数据创建完成\n")

# 测试1: 基础功能测试
cat("测试1: 基础功能...\n")
tryCatch({
  p1 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    order = "p_values",
    color_theme = "nature"
  )
  
  filename1 <- file.path(test_output_dir, "01_basic_function.pdf")
  ggsave(filename1, p1, width = 12, height = 8)
  cat("✓ 基础功能测试成功 - 已保存:", basename(filename1), "\n")
  
}, error = function(e) {
  cat("✗ 基础功能测试失败:", e$message, "\n")
})

# 测试2: 图例位置测试
cat("测试2: 图例位置...\n")
tryCatch({
  p2 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    color_theme = "science",
    legend_position = "bottom",
    legend_direction = "horizontal",
    legend_title = "Sample Groups"
  )
  
  filename2 <- file.path(test_output_dir, "02_legend_position.pdf")
  ggsave(filename2, p2, width = 12, height = 8)
  cat("✓ 图例位置测试成功 - 已保存:", basename(filename2), "\n")
  
}, error = function(e) {
  cat("✗ 图例位置测试失败:", e$message, "\n")
})

# 测试3: P值格式化测试
cat("测试3: P值格式化...\n")
tryCatch({
  p3 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    color_theme = "cell",
    pvalue_format = "smart",
    pvalue_stars = TRUE,
    pvalue_size = 12
  )
  
  filename3 <- file.path(test_output_dir, "03_pvalue_format.pdf")
  ggsave(filename3, p3, width = 12, height = 8)
  cat("✓ P值格式化测试成功 - 已保存:", basename(filename3), "\n")
  
}, error = function(e) {
  cat("✗ P值格式化测试失败:", e$message, "\n")
})

# 测试4: P值颜色编码测试
cat("测试4: P值颜色编码...\n")
tryCatch({
  p4 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    color_theme = "viridis",
    pvalue_format = "smart",
    pvalue_stars = TRUE,
    pvalue_colors = TRUE,
    pvalue_angle = 0
  )
  
  filename4 <- file.path(test_output_dir, "04_pvalue_colors.pdf")
  ggsave(filename4, p4, width = 12, height = 8)
  cat("✓ P值颜色编码测试成功 - 已保存:", basename(filename4), "\n")
  
}, error = function(e) {
  cat("✗ P值颜色编码测试失败:", e$message, "\n")
})

# 测试5: Pathway Class标注测试
cat("测试5: Pathway Class标注...\n")
tryCatch({
  p5 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    color_theme = "nature",
    pathway_class_text_size = "auto",
    pathway_class_text_color = "auto", 
    pathway_class_text_face = "bold",
    pathway_class_text_angle = 0
  )
  
  filename5 <- file.path(test_output_dir, "05_pathway_class.pdf")
  ggsave(filename5, p5, width = 12, height = 8)
  cat("✓ Pathway Class标注测试成功 - 已保存:", basename(filename5), "\n")
  
}, error = function(e) {
  cat("✗ Pathway Class标注测试失败:", e$message, "\n")
})

# 测试6: 综合功能测试
cat("测试6: 综合功能...\n")
tryCatch({
  p6 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    color_theme = "nature",
    smart_colors = TRUE,
    legend_position = "top",
    legend_direction = "horizontal", 
    legend_title = "Treatment Groups",
    legend_title_size = 14,
    legend_text_size = 12,
    pvalue_format = "smart",
    pvalue_stars = TRUE,
    pvalue_colors = FALSE,
    pathway_class_text_size = "auto",
    pathway_class_text_color = "auto"
  )
  
  filename6 <- file.path(test_output_dir, "06_comprehensive.pdf")
  ggsave(filename6, p6, width = 12, height = 8)
  cat("✓ 综合功能测试成功 - 已保存:", basename(filename6), "\n")
  
}, error = function(e) {
  cat("✗ 综合功能测试失败:", e$message, "\n")
})

# 测试7: 不同主题测试
cat("测试7: 不同主题...\n")
themes_to_test <- c("nature", "science", "cell", "colorblind_friendly")

for (theme in themes_to_test) {
  tryCatch({
    p <- pathway_errorbar(
      abundance = abundance_test,
      daa_results_df = daa_results_test,
      Group = Group_test,
      p_values_threshold = 0.05,
      color_theme = theme,
      legend_title = paste(theme, "Theme"),
      pvalue_format = "smart",
      pvalue_stars = TRUE,
      pathway_class_text_color = "auto"
    )
    
    filename <- file.path(test_output_dir, paste0("07_theme_", theme, ".pdf"))
    ggsave(filename, p, width = 12, height = 8)
    cat("✓", theme, "主题测试成功 - 已保存:", basename(filename), "\n")
    
  }, error = function(e) {
    cat("✗", theme, "主题测试失败:", e$message, "\n")
  })
}

# 生成测试总结
cat("\n=======================================\n")
cat("简化版图例和标注测试完成!\n")
cat("=======================================\n")

files <- list.files(test_output_dir, pattern = "\\.pdf$", full.names = FALSE)
cat("生成的测试文件 (", length(files), "个):\n")
for (i in seq_along(files)) {
  cat(sprintf("%2d. %s\n", i, files[i]))
}

cat("\n所有测试图片已保存到:", test_output_dir, "文件夹\n")
cat("请检查生成的PDF文件以验证功能。\n")