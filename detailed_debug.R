#\!/usr/bin/env Rscript

cat("详细调试测试...\n")

suppressPackageStartupMessages({
  library(ggpicrust2)
  library(dplyr)
  library(ggplot2)
})

source("R/legend_annotation_utils.R")
source("R/color_themes.R")
source("R/pathway_errorbar.R")

# 创建测试数据
set.seed(42)
abundance_test <- matrix(
  abs(rnorm(15, mean = 50, sd = 25)),
  nrow = 3,
  ncol = 5,
  dimnames = list(
    c("PWY001", "PWY002", "PWY003"),
    c("S1", "S2", "S3", "S4", "S5")
  )
)

Group_test <- c("Control", "Control", "Treatment", "Treatment", "Treatment")

daa_results_test <- data.frame(
  feature = c("PWY001", "PWY002", "PWY003"),
  p_adjust = c(0.001, 0.03, 0.08),
  method = "ALDEx2_Welch's t test",
  group1 = "Control",
  group2 = "Treatment",
  log_2_fold_change = c(2.1, -1.2, 0.5),
  description = c("Desc1", "Desc2", "Desc3"),
  pathway_class = c("Metabolism", "Signaling", "Transport"),
  stringsAsFactors = FALSE
)

cat("数据准备完成\n")
cat("abundance维度:", dim(abundance_test), "\n")
cat("Group长度:", length(Group_test), "\n")
cat("daa_results行数:", nrow(daa_results_test), "\n")

# 测试最简版本
cat("\n尝试最简单的调用...\n")
tryCatch({
  p <- pathway_errorbar(
    abundance = abundance_test,
    daa_results_df = daa_results_test,
    Group = Group_test,
    p_values_threshold = 0.05
  )
  cat("✓ 最简单调用成功\n")
  print(summary(p))
  
  # 检查图片是否为空
  if (inherits(p, "ggplot")) {
    cat("✓ 返回了ggplot对象\n")
    # 尝试保存
    ggsave("test_simple.pdf", p, width = 8, height = 6)
    cat("✓ 图片保存成功\n")
  } else {
    cat("✗ 返回的不是ggplot对象\n")
  }
  
}, error = function(e) {
  cat("✗ 最简单调用失败:\n")
  cat("错误:", e$message, "\n")
  if(exists("traceback")) {
    cat("调用栈:\n")
    traceback()
  }
})

cat("\n调试完成\n")
END < /dev/null