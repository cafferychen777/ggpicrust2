#\!/usr/bin/env Rscript

# 调试测试脚本
cat("开始调试pathway_errorbar函数...\n")

# 加载必要的库
suppressPackageStartupMessages({
  library(ggpicrust2)
  library(dplyr)
  library(ggplot2)
})

# 加载所需文件
source("R/legend_annotation_utils.R")
source("R/color_themes.R")
source("R/pathway_errorbar.R")

# 创建最简单的测试数据
set.seed(42)
n_pathways <- 5
n_samples <- 6

# 创建丰度矩阵
abundance_test <- matrix(
  abs(rnorm(n_pathways * n_samples, mean = 50, sd = 25)),
  nrow = n_pathways,
  ncol = n_samples,
  dimnames = list(
    paste0("PWY", 1:n_pathways),
    paste0("S", 1:n_samples)
  )
)

# 创建分组
Group_test <- rep(c("Control", "Treatment"), each = n_samples/2)

# 创建DAA结果数据
daa_results_test <- data.frame(
  feature = rownames(abundance_test),
  p_adjust = c(0.001, 0.01, 0.03, 0.07, 0.1),
  method = "ALDEx2_Welch's t test",
  group1 = "Control",
  group2 = "Treatment",
  log_2_fold_change = c(2, -1.5, 1, -0.5, 0.3),
  description = paste("Description for", rownames(abundance_test)),
  pathway_class = c("Metabolism", "Signaling", "Transport", "Energy", "Biosynthesis"),
  stringsAsFactors = FALSE
)

cat("测试数据创建完成\n")

# 先测试原始函数是否工作
cat("\n测试1: 原始pathway_errorbar函数...\n")
tryCatch({
  p_original <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05,
    order = "p_values"
  )
  cat("✓ 原始函数执行成功\n")
  
  # 尝试保存图片
  ggsave("debug_original.pdf", p_original, width = 10, height = 6)
  cat("✓ 原始图片保存成功\n")
  
}, error = function(e) {
  cat("✗ 原始函数执行失败:\n")
  cat("错误信息:", e$message, "\n")
})

cat("\n调试测试完成。\n")
EOF < /dev/null