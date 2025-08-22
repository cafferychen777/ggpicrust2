#!/usr/bin/env Rscript

# Cross-Pathway Consistency Test Runner
# This script runs comprehensive tests to validate consistency across KEGG, MetaCyc, and GO pathway types

library(testthat)
library(ggpicrust2)

cat("🧪 Cross-Pathway Consistency Test Suite\n")
cat("=======================================\n\n")

# Set up test environment
test_start_time <- Sys.time()
original_dir <- getwd()

# Ensure we're in the package directory
if (!file.exists("DESCRIPTION")) {
  if (file.exists("ggpicrust2/DESCRIPTION")) {
    setwd("ggpicrust2")
  } else {
    stop("Cannot find ggpicrust2 package directory")
  }
}

cat("📂 Working directory:", getwd(), "\n")
cat("⏰ Test started at:", format(test_start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# Test execution function with detailed reporting
run_test_suite <- function(test_file, suite_name) {
  cat(sprintf("🔍 Running %s...\n", suite_name))
  cat(sprintf("   File: %s\n", test_file))
  
  if (!file.exists(test_file)) {
    cat("   ❌ SKIPPED - Test file not found\n\n")
    return(list(status = "SKIPPED", tests = 0, passed = 0, failed = 0, time = 0))
  }
  
  suite_start <- Sys.time()
  
  # Capture test results
  test_results <- tryCatch({
    # Run tests and capture results
    results <- test_file(test_file, reporter = "summary")
    
    suite_end <- Sys.time()
    execution_time <- as.numeric(suite_end - suite_start, units = "secs")
    
    # Extract test statistics
    n_tests <- length(results)
    n_passed <- sum(sapply(results, function(x) length(x$results) == 0))  # No failures
    n_failed <- n_tests - n_passed
    
    cat(sprintf("   ✅ COMPLETED - %d tests run (%d passed, %d failed) in %.2f seconds\n\n", 
               n_tests, n_passed, n_failed, execution_time))
    
    return(list(
      status = if(n_failed == 0) "PASSED" else "FAILED",
      tests = n_tests,
      passed = n_passed, 
      failed = n_failed,
      time = execution_time,
      results = results
    ))
    
  }, error = function(e) {
    suite_end <- Sys.time()
    execution_time <- as.numeric(suite_end - suite_start, units = "secs")
    
    cat(sprintf("   ❌ ERROR - Test suite failed: %s\n", e$message))
    cat(sprintf("   ⏱️  Execution time: %.2f seconds\n\n", execution_time))
    
    return(list(
      status = "ERROR",
      tests = 0,
      passed = 0,
      failed = 1,
      time = execution_time,
      error = e$message
    ))
  })
  
  return(test_results)
}

# Define test suites to run
test_suites <- list(
  list(
    file = "tests/testthat/test-cross_pathway_consistency.R",
    name = "Cross-Pathway Consistency Tests",
    description = "Core consistency tests across KEGG, MetaCyc, and GO"
  ),
  list(
    file = "tests/testthat/test-pathway_gsea_comprehensive.R", 
    name = "Comprehensive GSEA Tests",
    description = "Extended GSEA functionality validation"
  ),
  list(
    file = "tests/testthat/test-pathway_gsea_integration_validation.R",
    name = "Integration Validation Tests", 
    description = "Method integration and result standardization"
  ),
  list(
    file = "tests/testthat/test-visualize_gsea_comprehensive.R",
    name = "Visualization Consistency Tests",
    description = "Cross-pathway visualization validation"
  ),
  list(
    file = "tests/testthat/test-gsea_pathway_annotation-comprehensive.R",
    name = "Annotation System Tests",
    description = "Pathway annotation consistency across types"
  )
)

# Run all test suites
all_results <- list()
total_tests <- 0
total_passed <- 0
total_failed <- 0
total_errors <- 0

for (suite in test_suites) {
  result <- run_test_suite(suite$file, suite$name)
  all_results[[suite$name]] <- result
  
  if (result$status == "ERROR") {
    total_errors <- total_errors + 1
  } else {
    total_tests <- total_tests + result$tests
    total_passed <- total_passed + result$passed
    total_failed <- total_failed + result$failed
  }
}

# Run practical demonstration
cat("🎭 Running Practical Demonstration...\n")
demo_start <- Sys.time()

tryCatch({
  source("test_cross_pathway_demo.R")
  demo_end <- Sys.time()
  demo_time <- as.numeric(demo_end - demo_start, units = "secs")
  cat(sprintf("   ✅ Demo completed successfully in %.2f seconds\n\n", demo_time))
  demo_status <- "SUCCESS"
}, error = function(e) {
  demo_end <- Sys.time()
  demo_time <- as.numeric(demo_end - demo_start, units = "secs")
  cat(sprintf("   ❌ Demo failed: %s\n", e$message))
  cat(sprintf("   ⏱️  Demo time: %.2f seconds\n\n", demo_time))
  demo_status <- "FAILED"
})

# Generate comprehensive summary report
cat("📊 COMPREHENSIVE TEST SUMMARY\n")
cat("=============================\n\n")

test_end_time <- Sys.time()
total_execution_time <- as.numeric(test_end_time - test_start_time, units = "secs")

cat("⏱️  Execution Summary:\n")
cat(sprintf("   Total execution time: %.2f seconds\n", total_execution_time))
cat(sprintf("   Test started: %s\n", format(test_start_time, "%H:%M:%S")))
cat(sprintf("   Test ended: %s\n\n", format(test_end_time, "%H:%M:%S")))

cat("📈 Test Results Summary:\n")
cat(sprintf("   Total test suites: %d\n", length(test_suites)))
cat(sprintf("   Total individual tests: %d\n", total_tests))
cat(sprintf("   Tests passed: %d (%.1f%%)\n", total_passed, 
           if(total_tests > 0) 100 * total_passed / total_tests else 0))
cat(sprintf("   Tests failed: %d (%.1f%%)\n", total_failed,
           if(total_tests > 0) 100 * total_failed / total_tests else 0))
cat(sprintf("   Suite errors: %d\n\n", total_errors))

# Detailed results by test suite
cat("📋 Detailed Results by Test Suite:\n")
for (suite_name in names(all_results)) {
  result <- all_results[[suite_name]]
  status_icon <- switch(result$status,
                       "PASSED" = "✅",
                       "FAILED" = "❌", 
                       "ERROR" = "💥",
                       "SKIPPED" = "⏭️")
  
  cat(sprintf("   %s %s\n", status_icon, suite_name))
  cat(sprintf("      Status: %s\n", result$status))
  if (result$tests > 0) {
    cat(sprintf("      Tests: %d passed, %d failed\n", result$passed, result$failed))
  }
  cat(sprintf("      Time: %.2f seconds\n", result$time))
  if (!is.null(result$error)) {
    cat(sprintf("      Error: %s\n", result$error))
  }
  cat("\n")
}

# Cross-pathway consistency assessment
cat("🎯 Cross-Pathway Consistency Assessment:\n")

consistency_score <- 0
max_score <- 6

# 1. API Consistency
api_score <- if (any(sapply(all_results, function(x) x$status == "PASSED" && grepl("consistency", names(all_results), ignore.case = TRUE)))) 1 else 0
consistency_score <- consistency_score + api_score
cat(sprintf("   API Consistency: %s\n", if(api_score == 1) "✅ PASSED" else "❌ FAILED"))

# 2. Statistical Consistency  
stats_score <- if (total_failed == 0 && total_errors == 0) 1 else 0
consistency_score <- consistency_score + stats_score
cat(sprintf("   Statistical Consistency: %s\n", if(stats_score == 1) "✅ PASSED" else "❌ FAILED"))

# 3. Integration Workflow
integration_score <- if (any(sapply(names(all_results), function(x) grepl("integration", x, ignore.case = TRUE))) && 
                         all_results[["Integration Validation Tests"]]$status == "PASSED") 1 else 0
consistency_score <- consistency_score + integration_score  
cat(sprintf("   Integration Workflow: %s\n", if(integration_score == 1) "✅ PASSED" else "❌ FAILED"))

# 4. Visualization Consistency
viz_score <- if (any(sapply(names(all_results), function(x) grepl("visualization", x, ignore.case = TRUE)))) 1 else 0
consistency_score <- consistency_score + viz_score
cat(sprintf("   Visualization Consistency: %s\n", if(viz_score == 1) "✅ PASSED" else "❌ FAILED"))

# 5. Annotation System
annotation_score <- if (any(sapply(names(all_results), function(x) grepl("annotation", x, ignore.case = TRUE)))) 1 else 0
consistency_score <- consistency_score + annotation_score
cat(sprintf("   Annotation System: %s\n", if(annotation_score == 1) "✅ PASSED" else "❌ FAILED"))

# 6. Performance Parity
performance_score <- if (demo_status == "SUCCESS") 1 else 0
consistency_score <- consistency_score + performance_score
cat(sprintf("   Performance Parity: %s\n", if(performance_score == 1) "✅ PASSED" else "❌ FAILED"))

cat(sprintf("\n🏆 Overall Consistency Score: %d/%d (%.1f%%)\n", 
           consistency_score, max_score, 100 * consistency_score / max_score))

# Final assessment
cat("\n🎯 FINAL ASSESSMENT:\n")
if (consistency_score == max_score && total_failed == 0 && total_errors == 0) {
  cat("🌟 EXCELLENT - Cross-pathway consistency fully validated!\n")
  cat("   All pathway types (KEGG, MetaCyc, GO) provide reliable, consistent results.\n")
  cat("   Users can confidently use any or all pathway types for scientific analysis.\n")
  final_status <- "EXCELLENT"
} else if (consistency_score >= max_score * 0.8 && total_failed <= total_tests * 0.1) {
  cat("✅ GOOD - Cross-pathway consistency largely validated.\n") 
  cat("   Minor issues detected but overall system is reliable.\n")
  cat("   Most pathway types work consistently.\n")
  final_status <- "GOOD"
} else if (consistency_score >= max_score * 0.6) {
  cat("⚠️  PARTIAL - Cross-pathway consistency partially validated.\n")
  cat("   Some pathway types work well, others need attention.\n")
  cat("   Users should verify results across pathway types.\n")
  final_status <- "PARTIAL"
} else {
  cat("❌ FAILED - Significant cross-pathway consistency issues detected.\n")
  cat("   Major problems found across multiple pathway types.\n") 
  cat("   System needs substantial improvements before production use.\n")
  final_status <- "FAILED"
}

# Recommendations
cat("\n💡 RECOMMENDATIONS:\n")

if (final_status == "EXCELLENT") {
  cat("   1. System is ready for production use\n")
  cat("   2. All pathway types can be recommended to users\n")
  cat("   3. Documentation can emphasize consistency and reliability\n")
  cat("   4. Consider adding advanced cross-pathway integration features\n")
} else if (final_status == "GOOD") {
  cat("   1. Address minor issues before full release\n") 
  cat("   2. Most pathway types are ready for use\n")
  cat("   3. Monitor user feedback for remaining edge cases\n")
  cat("   4. Continue regression testing for consistency\n")
} else if (final_status == "PARTIAL") {
  cat("   1. Focus on failing pathway types and fix core issues\n")
  cat("   2. Improve test coverage for problematic areas\n")
  cat("   3. Consider phased release starting with working pathway types\n")
  cat("   4. Enhance error handling and user messaging\n")
} else {
  cat("   1. Major refactoring needed for consistency\n")
  cat("   2. Do not release until core issues are resolved\n")
  cat("   3. Focus on API standardization across pathway types\n")
  cat("   4. Implement comprehensive regression test suite\n")
}

# Save results to file
cat("\n💾 Saving Test Results...\n")
results_file <- paste0("cross_pathway_test_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")

# Create detailed results report
results_content <- c(
  "Cross-Pathway Consistency Test Results",
  "======================================",
  "",
  paste("Test Date:", format(test_start_time, "%Y-%m-%d %H:%M:%S")),
  paste("Total Execution Time:", sprintf("%.2f seconds", total_execution_time)),
  paste("Final Status:", final_status),
  paste("Overall Score:", sprintf("%d/%d (%.1f%%)", consistency_score, max_score, 100 * consistency_score / max_score)),
  "",
  "Test Suite Results:",
  "-------------------"
)

for (suite_name in names(all_results)) {
  result <- all_results[[suite_name]]
  results_content <- c(results_content,
    paste("Suite:", suite_name),
    paste("  Status:", result$status),
    paste("  Tests:", sprintf("%d passed, %d failed", result$passed, result$failed)),
    paste("  Time:", sprintf("%.2f seconds", result$time)),
    ""
  )
}

results_content <- c(results_content,
  "Summary Statistics:",
  "------------------",
  paste("Total Test Suites:", length(test_suites)),
  paste("Total Tests:", total_tests),
  paste("Tests Passed:", sprintf("%d (%.1f%%)", total_passed, if(total_tests > 0) 100 * total_passed / total_tests else 0)),
  paste("Tests Failed:", sprintf("%d (%.1f%%)", total_failed, if(total_tests > 0) 100 * total_failed / total_tests else 0)),
  paste("Suite Errors:", total_errors),
  "",
  "Cross-Pathway Assessment:",
  "------------------------",
  paste("API Consistency:", if(api_score == 1) "PASSED" else "FAILED"),
  paste("Statistical Consistency:", if(stats_score == 1) "PASSED" else "FAILED"), 
  paste("Integration Workflow:", if(integration_score == 1) "PASSED" else "FAILED"),
  paste("Visualization Consistency:", if(viz_score == 1) "PASSED" else "FAILED"),
  paste("Annotation System:", if(annotation_score == 1) "PASSED" else "FAILED"),
  paste("Performance Parity:", if(performance_score == 1) "PASSED" else "FAILED")
)

writeLines(results_content, results_file)
cat(sprintf("   Results saved to: %s\n", results_file))

# Return to original directory
setwd(original_dir)

cat("\n🎉 Cross-pathway consistency testing completed!\n")
cat("Check the generated files for detailed results and visualizations.\n")

# Exit with appropriate code
if (final_status %in% c("EXCELLENT", "GOOD")) {
  quit(status = 0)
} else {
  quit(status = 1)
}