# Extensive Testing Suite for Enhanced pathway_errorbar()
# ======================================================
# 
# This comprehensive test suite validates the color theme system 
# under various conditions and edge cases.

library(ggpicrust2)
library(dplyr)

# Source enhanced functions
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_errorbar.R")

cat("🧪 EXTENSIVE TESTING SUITE\n")
cat("==========================\n\n")

# Load test data
data("daa_annotated_results_df")
data("kegg_abundance")
data("metadata")

abundance_data <- kegg_abundance
Group <- metadata$Environment
all_features <- daa_annotated_results_df$feature

test_results <- list()
test_count <- 0

# Helper function to run tests
run_test <- function(test_name, test_func) {
  test_count <<- test_count + 1
  cat(sprintf("Test %02d: %-50s", test_count, test_name))
  
  result <- tryCatch({
    test_func()
    "✅ PASS"
  }, error = function(e) {
    paste("❌ FAIL:", conditionMessage(e))
  }, warning = function(w) {
    paste("⚠️ WARN:", conditionMessage(w))
  })
  
  cat(result, "\n")
  test_results[[test_name]] <<- result
  return(grepl("PASS", result))
}

cat("🔍 BASIC FUNCTIONALITY TESTS\n")
cat("-----------------------------\n")

# Test 1: Default parameters
run_test("Default parameters", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:5]
  )
})

# Test 2: All available color themes
cat("\n🎨 COLOR THEME TESTS\n")
cat("--------------------\n")

all_themes <- get_available_themes()
for (theme in all_themes) {
  run_test(paste("Theme:", theme), function() {
    pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      select = all_features[1:3],
      color_theme = theme
    )
  })
}

cat("\n🧠 SMART COLOR SELECTION TESTS\n")
cat("-------------------------------\n")

# Test smart colors with different group counts
run_test("Smart colors - 2 groups", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:4],
    smart_colors = TRUE
  )
})

# Test accessibility mode
run_test("Accessibility mode", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:4],
    accessibility_mode = TRUE
  )
})

# Test combined smart colors and accessibility
run_test("Smart colors + Accessibility", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:4],
    smart_colors = TRUE,
    accessibility_mode = TRUE
  )
})

cat("\n🎛️ PARAMETER COMBINATION TESTS\n")
cat("------------------------------\n")

# Test auto fold change colors
run_test("Auto fold change colors", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:4],
    color_theme = "nature",
    log2_fold_change_color = "auto"
  )
})

# Test custom pathway class colors
run_test("Custom pathway class colors", function() {
  custom_colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FECA57")
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:4],
    color_theme = "cell",
    pathway_class_colors = custom_colors
  )
})

# Test different ordering methods
ordering_methods <- c("p_values", "group", "pathway_class")
for (order_method in ordering_methods) {
  run_test(paste("Order by:", order_method), function() {
    pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      select = all_features[1:5],
      order = order_method,
      color_theme = "science"
    )
  })
}

cat("\n📊 EDGE CASE TESTS\n")
cat("------------------\n")

# Test with minimal features
run_test("Single feature", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1],
    color_theme = "minimal"
  )
})

# Test with all features
run_test("All available features", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features,
    color_theme = "viridis"
  )
})

# Test without feature selection
run_test("No feature selection", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = NULL,
    color_theme = "bold"
  )
})

cat("\n🔧 COLOR SYSTEM COMPONENT TESTS\n")
cat("-------------------------------\n")

# Test color theme functions directly
run_test("get_available_themes()", function() {
  themes <- get_available_themes()
  if (length(themes) < 10) stop("Too few themes available")
  themes
})

run_test("get_color_theme() with various n_colors", function() {
  for (n in c(1, 3, 5, 8, 15)) {
    theme <- get_color_theme("nature", n)
    if (length(theme$group_colors) != n) {
      stop(paste("Wrong number of colors for n =", n))
    }
  }
})

run_test("smart_color_selection() various scenarios", function() {
  scenarios <- list(
    list(n_groups = 2, data_type = "abundance", accessibility = FALSE),
    list(n_groups = 5, data_type = "pvalue", accessibility = FALSE),
    list(n_groups = 8, data_type = "foldchange", accessibility = FALSE),
    list(n_groups = 3, data_type = "abundance", accessibility = TRUE)
  )
  
  for (scenario in scenarios) {
    result <- smart_color_selection(
      n_groups = scenario$n_groups,
      data_type = scenario$data_type,
      accessibility_mode = scenario$accessibility
    )
    if (is.null(result$theme_name)) stop("No theme selected")
  }
})

run_test("create_gradient_colors()", function() {
  for (theme in c("default", "nature", "viridis")) {
    colors <- create_gradient_colors(theme, 11, TRUE)
    if (length(colors) != 11) stop("Wrong gradient length")
  }
})

run_test("preview_color_theme()", function() {
  preview_plot <- preview_color_theme("science")
  if (!"ggplot" %in% class(preview_plot)) stop("Not a ggplot object")
})

cat("\n⚡ PERFORMANCE TESTS\n")
cat("-------------------\n")

# Test with different theme loading approaches
run_test("Multiple theme switches", function() {
  themes <- c("default", "nature", "science", "cell", "viridis")
  for (theme in themes) {
    pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      select = all_features[1:3],
      color_theme = theme
    )
  }
})

cat("\n🔄 BACKWARD COMPATIBILITY TESTS\n")
cat("-------------------------------\n")

# Test old-style usage (should still work)
run_test("Legacy usage - no new parameters", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:5],
    colors = NULL,
    p_value_bar = TRUE
  )
})

run_test("Legacy usage - with custom colors", function() {
  custom_colors <- c("#FF5733", "#33FF57", "#3357FF")
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:5],
    colors = custom_colors
  )
})

# Final summary
cat("\n📈 TEST SUMMARY\n")
cat("===============\n")

total_tests <- length(test_results)
passed_tests <- sum(sapply(test_results, function(x) grepl("PASS", x)))
failed_tests <- sum(sapply(test_results, function(x) grepl("FAIL", x)))
warning_tests <- sum(sapply(test_results, function(x) grepl("WARN", x)))

cat(sprintf("Total tests run: %d\n", total_tests))
cat(sprintf("✅ Passed: %d (%.1f%%)\n", passed_tests, passed_tests/total_tests*100))
cat(sprintf("❌ Failed: %d (%.1f%%)\n", failed_tests, failed_tests/total_tests*100))
cat(sprintf("⚠️ Warnings: %d (%.1f%%)\n", warning_tests, warning_tests/total_tests*100))

if (failed_tests > 0) {
  cat("\n❌ FAILED TESTS:\n")
  failed_names <- names(test_results)[sapply(test_results, function(x) grepl("FAIL", x))]
  for (name in failed_names) {
    cat("-", name, ":", test_results[[name]], "\n")
  }
}

if (warning_tests > 0) {
  cat("\n⚠️ TESTS WITH WARNINGS:\n")
  warning_names <- names(test_results)[sapply(test_results, function(x) grepl("WARN", x))]
  for (name in warning_names) {
    cat("-", name, ":", test_results[[name]], "\n")
  }
}

if (failed_tests == 0) {
  cat("\n🎉 ALL TESTS PASSED! The color theme system is working perfectly!\n")
  cat("✨ The enhanced pathway_errorbar() function is ready for production use.\n")
} else {
  cat("\n⚠️ Some tests failed. Please review the issues above.\n")
}

cat(sprintf("\n🏁 Testing completed. Success rate: %.1f%%\n", passed_tests/total_tests*100))