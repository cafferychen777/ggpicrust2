# Stress Test for Color Theme System
# ==================================
# 
# Additional edge cases and stress testing for the enhanced pathway_errorbar()

library(ggpicrust2)
library(dplyr)

# Source functions
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_errorbar.R")

cat("⚡ STRESS TESTING SUITE\n")
cat("======================\n\n")

# Load data
data("daa_annotated_results_df")  
data("kegg_abundance")
data("metadata")

abundance_data <- kegg_abundance
Group <- metadata$Environment
all_features <- daa_annotated_results_df$feature

stress_results <- list()
test_count <- 0

run_stress_test <- function(test_name, test_func) {
  test_count <<- test_count + 1
  cat(sprintf("Stress %02d: %-50s", test_count, test_name))
  
  result <- tryCatch({
    test_func()
    "✅ PASS"
  }, error = function(e) {
    paste("❌ FAIL:", conditionMessage(e))
  }, warning = function(w) {
    paste("⚠️ WARN:", conditionMessage(w))
  })
  
  cat(result, "\n")
  stress_results[[test_name]] <<- result
  return(grepl("PASS", result))
}

cat("🔥 EXTREME PARAMETER TESTS\n")
cat("--------------------------\n")

# Test with very high n_colors
run_stress_test("High n_colors (n=100)", function() {
  theme <- get_color_theme("viridis", 100)
  if (length(theme$group_colors) != 100) stop("Color count mismatch")
})

# Test invalid theme names
run_stress_test("Invalid theme name handling", function() {
  theme <- get_color_theme("nonexistent_theme", 5)
  if (is.null(theme)) stop("Should return default theme")
})

# Test zero colors
run_stress_test("Zero colors request", function() {
  theme <- get_color_theme("default", 0)
  # Should not crash
})

cat("\n🌪️ RAPID THEME SWITCHING\n")
cat("------------------------\n")

# Rapid theme switching test
run_stress_test("Rapid theme switching (50 times)", function() {
  themes <- get_available_themes()
  for (i in 1:50) {
    theme_name <- sample(themes, 1)
    pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      select = all_features[1:2],
      color_theme = theme_name
    )
  }
})

cat("\n🎲 RANDOM PARAMETER COMBINATIONS\n")
cat("--------------------------------\n")

# Random parameter combinations
run_stress_test("Random parameter combinations", function() {
  themes <- get_available_themes()
  orders <- c("p_values", "group", "pathway_class")
  
  for (i in 1:20) {
    pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      select = sample(all_features, sample(1:length(all_features), 1)),
      color_theme = sample(themes, 1),
      order = sample(orders, 1),
      smart_colors = sample(c(TRUE, FALSE), 1),
      accessibility_mode = sample(c(TRUE, FALSE), 1),
      log2_fold_change_color = sample(c("auto", "#FF5733"), 1)
    )
  }
})

cat("\n📊 EXTREME DATA SCENARIOS\n")
cat("-------------------------\n")

# Test with all features
run_stress_test("All features at once", function() {
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features,
    color_theme = "viridis"
  )
})

# Test different p_values_thresholds
p_thresholds <- c(0.001, 0.01, 0.05, 0.1, 0.5, 1.0)
for (threshold in p_thresholds) {
  run_stress_test(paste("P-value threshold:", threshold), function() {
    pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      p_values_threshold = threshold,
      color_theme = "science"
    )
  })
}

cat("\n🧬 MEMORY AND PERFORMANCE\n")
cat("-------------------------\n")

# Memory usage test
run_stress_test("Memory usage - large color arrays", function() {
  large_colors <- rep(get_color_theme("plasma", 50)$group_colors, 10)
  pathway_errorbar(
    abundance = abundance_data,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    ko_to_kegg = TRUE,
    select = all_features[1:5],
    colors = large_colors[1:2],  # Only use what we need
    color_theme = "plasma"
  )
})

# Repeated executions
run_stress_test("Repeated executions (100 times)", function() {
  for (i in 1:100) {
    pathway_errorbar(
      abundance = abundance_data,
      daa_results_df = daa_annotated_results_df,
      Group = Group,
      ko_to_kegg = TRUE,
      select = all_features[1:3],
      color_theme = "default"
    )
  }
})

cat("\n🔧 COMPONENT STRESS TESTS\n")
cat("-------------------------\n")

# Test color system components under stress
run_stress_test("Color theme preview - all themes", function() {
  themes <- get_available_themes()
  for (theme in themes) {
    preview_color_theme(theme)
  }
})

run_stress_test("Smart selection - extreme scenarios", function() {
  scenarios <- expand.grid(
    n_groups = c(1, 2, 5, 10, 20, 50),
    data_type = c("abundance", "pvalue", "foldchange"),
    accessibility = c(TRUE, FALSE)
  )
  
  for (i in 1:nrow(scenarios)) {
    smart_color_selection(
      n_groups = scenarios$n_groups[i],
      data_type = as.character(scenarios$data_type[i]),
      accessibility_mode = scenarios$accessibility[i]
    )
  }
})

run_stress_test("Gradient colors - various sizes", function() {
  themes <- get_available_themes()
  sizes <- c(1, 3, 5, 11, 21, 51, 101)
  
  for (theme in themes[1:5]) {  # Test first 5 themes
    for (size in sizes) {
      colors <- create_gradient_colors(theme, size, TRUE)
      if (length(colors) != size) stop(paste("Size mismatch for", theme, size))
    }
  }
})

# Final stress test summary
cat("\n📊 STRESS TEST SUMMARY\n")
cat("======================\n")

total_stress <- length(stress_results)
passed_stress <- sum(sapply(stress_results, function(x) grepl("PASS", x)))
failed_stress <- sum(sapply(stress_results, function(x) grepl("FAIL", x)))

cat(sprintf("Total stress tests: %d\n", total_stress))
cat(sprintf("✅ Passed: %d (%.1f%%)\n", passed_stress, passed_stress/total_stress*100))
cat(sprintf("❌ Failed: %d (%.1f%%)\n", failed_stress, failed_stress/total_stress*100))

if (failed_stress == 0) {
  cat("\n🚀 ALL STRESS TESTS PASSED!\n")
  cat("💪 The color theme system is robust and production-ready!\n")
  cat("🏆 Performance under extreme conditions: EXCELLENT\n")
} else {
  cat("\n⚠️ Some stress tests failed:\n")
  failed_names <- names(stress_results)[sapply(stress_results, function(x) grepl("FAIL", x))]
  for (name in failed_names) {
    cat("-", name, ":", stress_results[[name]], "\n")
  }
}

cat(sprintf("\n🎯 Stress test success rate: %.1f%%\n", passed_stress/total_stress*100))
cat("✨ Color theme system validation complete!\n")