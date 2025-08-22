#!/usr/bin/env Rscript

# 图例和标注美化系统完整测试脚本
# Comprehensive test script for legend and annotation beautification system

cat("========================================\n")
cat("图例和标注美化系统完整测试\n")
cat("Comprehensive Legend & Annotation Test\n")
cat("========================================\n\n")

# 加载必要的库
suppressPackageStartupMessages({
  library(ggpicrust2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# 创建测试结果文件夹
test_output_dir <- "legend_annotation_test_results"
if (!dir.exists(test_output_dir)) {
  dir.create(test_output_dir, recursive = TRUE)
}
cat("✓ 测试结果将保存到:", test_output_dir, "\n\n")

# 加载所需文件
source("R/legend_annotation_utils.R")
source("R/color_themes.R")
source("R/pathway_errorbar.R")

# 创建测试数据集
set.seed(42)
n_pathways <- 20
n_samples <- 24

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

# 创建分组
Group_test <- rep(c("Control", "Treatment_A", "Treatment_B"), each = n_samples/3)

# 创建DAA结果数据
pathway_classes <- c("Metabolism", "Signaling", "Transport", "Energy", "Biosynthesis")
daa_results_test <- data.frame(
  feature = rownames(abundance_test),
  p_adjust = c(
    runif(5, 0.0001, 0.0009),   # 极显著
    runif(5, 0.001, 0.009),     # 很显著
    runif(5, 0.01, 0.049),      # 显著
    runif(5, 0.05, 0.2)         # 不显著
  ),
  method = "ALDEx2_Welch's t test",
  group1 = "Control",
  group2 = "Treatment_A",
  log_2_fold_change = c(
    rnorm(10, 2, 1),    # 上调
    rnorm(10, -2, 1)    # 下调
  ),
  description = paste("Pathway description for", rownames(abundance_test)),
  pathway_class = sample(pathway_classes, n_pathways, replace = TRUE),
  stringsAsFactors = FALSE
)

cat("✓ 测试数据创建完成\n")
cat("  - 通路数量:", n_pathways, "\n")
cat("  - 样本数量:", n_samples, "\n")
cat("  - 分组:", paste(unique(Group_test), collapse = ", "), "\n")
cat("  - 显著通路数量:", sum(daa_results_test$p_adjust < 0.05), "\n\n")

# 测试函数
run_test <- function(test_name, plot_func, width = 12, height = 8) {
  cat("正在运行:", test_name, "\n")
  
  tryCatch({
    # 执行绘图函数
    p <- plot_func()
    
    # 保存图片
    filename <- file.path(test_output_dir, paste0(gsub("[^A-Za-z0-9_]", "_", test_name), ".pdf"))
    ggsave(filename, p, width = width, height = height, device = "pdf")
    
    cat("✓", test_name, "成功 - 已保存:", basename(filename), "\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗", test_name, "失败:", e$message, "\n")
    return(FALSE)
  })
}

# 开始测试
cat("开始图例和标注测试...\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

test_results <- list()

# 测试1: 基础图例位置测试
test_results$legend_positions <- run_test("1. 图例位置测试", function() {
  plots <- list()
  positions <- c("top", "bottom", "left", "right")
  
  for (pos in positions) {
    p <- pathway_errorbar(
      abundance = abundance_test,
      daa_results_df = daa_results_test,
      Group = Group_test,
      p_values_threshold = 0.05,
      order = "p_values",
      select = head(daa_results_test$feature, 12),
      color_theme = "nature",
      legend_position = pos,
      legend_title = paste("Legend", pos)
    )
    plots[[pos]] <- p + ggtitle(paste("Legend Position:", pos))
  }
  
  wrap_plots(plots, ncol = 2)
}, width = 16, height = 12)

# 测试2: 图例方向和样式
test_results$legend_styles <- run_test("2. 图例样式测试", function() {
  plots <- list()
  
  # 水平图例
  p1 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 10),
    color_theme = "science",
    legend_position = "top",
    legend_direction = "horizontal",
    legend_title = "Sample Groups",
    legend_title_size = 14,
    legend_text_size = 12,
    legend_key_size = 1.0
  ) + ggtitle("水平图例 (Horizontal Legend)")
  
  # 垂直图例
  p2 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 10),
    color_theme = "science",
    legend_position = "right",
    legend_direction = "vertical",
    legend_title = "Sample Groups",
    legend_title_size = 14,
    legend_text_size = 12,
    legend_key_size = 1.0
  ) + ggtitle("垂直图例 (Vertical Legend)")
  
  # 多列图例
  p3 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 10),
    color_theme = "cell",
    legend_position = "bottom",
    legend_direction = "horizontal",
    legend_title = "Treatment Groups",
    legend_ncol = 3,
    legend_text_size = 10
  ) + ggtitle("多列图例 (Multi-column Legend)")
  
  # 无图例
  p4 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 10),
    color_theme = "viridis",
    legend_position = "none"
  ) + ggtitle("无图例 (No Legend)")
  
  wrap_plots(p1, p2, p3, p4, ncol = 2)
}, width = 16, height = 12)

# 测试3: P值格式化测试
test_results$pvalue_formats <- run_test("3. P值格式化测试", function() {
  plots <- list()
  formats <- c("numeric", "scientific", "smart", "stars_only", "combined")
  
  for (i in 1:length(formats)) {
    fmt <- formats[i]
    p <- pathway_errorbar(
      abundance = abundance_test,
      daa_results_df = daa_results_test,
      Group = Group_test,
      p_values_threshold = 0.05,
      select = head(daa_results_test$feature, 8),
      color_theme = "nature",
      pvalue_format = fmt,
      pvalue_stars = TRUE,
      pvalue_size = 11
    ) + ggtitle(paste("P-value Format:", fmt))
    plots[[fmt]] <- p
  }
  
  wrap_plots(plots[1:4], ncol = 2) # 只显示前4个
}, width = 16, height = 12)

# 测试4: P值颜色编码测试
test_results$pvalue_colors <- run_test("4. P值颜色编码测试", function() {
  plots <- list()
  
  # 标准颜色编码
  p1 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 12),
    color_theme = "high_contrast",
    pvalue_format = "smart",
    pvalue_stars = TRUE,
    pvalue_colors = TRUE,
    pvalue_size = 12
  ) + ggtitle("P值颜色编码 (P-value Color Coding)")
  
  # 不同角度的P值显示
  p2 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 12),
    color_theme = "lancet",
    pvalue_format = "combined",
    pvalue_stars = TRUE,
    pvalue_colors = TRUE,
    pvalue_angle = 45,
    pvalue_size = 10
  ) + ggtitle("倾斜P值显示 (Angled P-values)")
  
  wrap_plots(p1, p2, ncol = 1)
}, width = 14, height = 16)

# 测试5: Pathway Class标注测试
test_results$pathway_class <- run_test("5. Pathway Class标注测试", function() {
  plots <- list()
  
  # 自动大小和颜色
  p1 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 15),
    color_theme = "nature",
    pathway_class_text_size = "auto",
    pathway_class_text_color = "auto",
    pathway_class_text_face = "bold",
    pathway_class_position = "left"
  ) + ggtitle("自动大小和颜色 (Auto Size & Color)")
  
  # 自定义样式
  p2 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 15),
    color_theme = "cell",
    pathway_class_text_size = 4.5,
    pathway_class_text_color = "#2166ac",
    pathway_class_text_face = "italic",
    pathway_class_text_angle = 15,
    pathway_class_position = "left"
  ) + ggtitle("自定义样式 (Custom Style)")
  
  # 右侧位置
  p3 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 15),
    color_theme = "science",
    pathway_class_text_size = 4,
    pathway_class_text_color = "#d62728",
    pathway_class_text_face = "bold",
    pathway_class_position = "right"
  ) + ggtitle("右侧位置 (Right Position)")
  
  wrap_plots(p1, p2, p3, ncol = 1)
}, width = 14, height = 18)

# 测试6: 主题集成测试
test_results$theme_integration <- run_test("6. 主题集成测试", function() {
  plots <- list()
  themes <- c("nature", "science", "cell", "viridis")
  
  for (theme in themes) {
    p <- pathway_errorbar(
      abundance = abundance_test,
      daa_results_df = daa_results_test,
      Group = Group_test,
      p_values_threshold = 0.05,
      select = head(daa_results_test$feature, 10),
      color_theme = theme,
      legend_position = "top",
      legend_title = paste(theme, "Theme"),
      pvalue_format = "smart",
      pvalue_stars = TRUE,
      pvalue_colors = TRUE,
      pathway_class_text_color = "auto"
    ) + ggtitle(paste("Theme:", theme))
    plots[[theme]] <- p
  }
  
  wrap_plots(plots, ncol = 2)
}, width = 16, height = 12)

# 测试7: 综合功能测试
test_results$comprehensive <- run_test("7. 综合功能测试", function() {
  plots <- list()
  
  # 学术期刊风格1 - Nature
  p1 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 12),
    color_theme = "nature",
    smart_colors = TRUE,
    legend_position = "top",
    legend_direction = "horizontal",
    legend_title = "Treatment Groups",
    legend_title_size = 13,
    legend_text_size = 11,
    pvalue_format = "smart",
    pvalue_stars = TRUE,
    pvalue_colors = FALSE,
    pathway_class_text_size = "auto",
    pathway_class_text_color = "auto"
  ) + ggtitle("Nature Journal Style")
  
  # 学术期刊风格2 - Science
  p2 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 12),
    color_theme = "science",
    legend_position = "bottom",
    legend_direction = "horizontal",
    legend_title = "Sample Types",
    legend_ncol = 3,
    pvalue_format = "combined",
    pvalue_stars = TRUE,
    pvalue_colors = TRUE,
    pathway_class_text_size = 4,
    pathway_class_text_face = "bold"
  ) + ggtitle("Science Journal Style")
  
  wrap_plots(p1, p2, ncol = 1)
}, width = 14, height = 16)

# 测试8: 可访问性测试
test_results$accessibility <- run_test("8. 可访问性测试", function() {
  plots <- list()
  
  # 色盲友好
  p1 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 12),
    color_theme = "colorblind_friendly",
    accessibility_mode = TRUE,
    legend_position = "top",
    pvalue_format = "smart",
    pvalue_stars = TRUE
  ) + ggtitle("色盲友好设计 (Colorblind Friendly)")
  
  # 高对比度
  p2 <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    select = head(daa_results_test$feature, 12),
    color_theme = "high_contrast",
    legend_position = "top",
    pvalue_format = "stars_only",
    pvalue_stars = TRUE,
    pvalue_colors = TRUE
  ) + ggtitle("高对比度设计 (High Contrast)")
  
  wrap_plots(p1, p2, ncol = 1)
}, width = 14, height = 16)

# 生成测试总结
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("测试完成总结\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

success_count <- sum(unlist(test_results))
total_tests <- length(test_results)

cat("总测试数量:", total_tests, "\n")
cat("成功测试数量:", success_count, "\n")
cat("测试成功率:", round(success_count/total_tests * 100, 1), "%\n\n")

cat("生成的测试文件:\n")
files <- list.files(test_output_dir, pattern = "\\.pdf$", full.names = FALSE)
for (i in seq_along(files)) {
  cat(sprintf("%2d. %s\n", i, files[i]))
}

cat("\n图例和标注美化系统测试完成！\n")
cat("所有测试图片已保存到:", test_output_dir, "文件夹\n")
cat("请检查生成的PDF文件以验证功能正确性。\n")

# 创建测试报告
report_file <- file.path(test_output_dir, "test_report.md")
cat("# 图例和标注美化系统测试报告\n\n", file = report_file)
cat("## 测试概述\n", file = report_file, append = TRUE)
cat("- 测试日期:", Sys.Date(), "\n", file = report_file, append = TRUE)
cat("- 总测试数:", total_tests, "\n", file = report_file, append = TRUE)
cat("- 成功测试:", success_count, "\n", file = report_file, append = TRUE)
cat("- 成功率:", round(success_count/total_tests * 100, 1), "%\n\n", file = report_file, append = TRUE)

cat("## 测试内容\n", file = report_file, append = TRUE)
cat("1. 图例位置测试 - 验证top/bottom/left/right位置\n", file = report_file, append = TRUE)
cat("2. 图例样式测试 - 验证水平/垂直方向，多列布局\n", file = report_file, append = TRUE)
cat("3. P值格式化测试 - 验证numeric/scientific/smart/stars_only格式\n", file = report_file, append = TRUE)
cat("4. P值颜色编码测试 - 验证显著性颜色和角度显示\n", file = report_file, append = TRUE) 
cat("5. Pathway Class标注测试 - 验证自动/自定义样式和位置\n", file = report_file, append = TRUE)
cat("6. 主题集成测试 - 验证与不同颜色主题的集成\n", file = report_file, append = TRUE)
cat("7. 综合功能测试 - 验证期刊风格综合应用\n", file = report_file, append = TRUE)
cat("8. 可访问性测试 - 验证色盲友好和高对比度设计\n", file = report_file, append = TRUE)

cat("\n测试报告已保存:", report_file, "\n")