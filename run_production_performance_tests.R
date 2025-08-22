#!/usr/bin/env Rscript

# Production Performance Test Runner
# 
# Quick execution script for comprehensive GSEA performance testing
# Run this script to validate production deployment readiness

# Suppress startup messages for cleaner output
options(warn = 1)
suppressMessages(library(ggpicrust2))

cat("🚀 PRODUCTION PERFORMANCE TESTING\n")
cat("=================================\n\n")

cat("Loading performance test suite...\n")

# Source the master test suite
source("master_performance_test_suite.R")

cat("✓ Test suite loaded\n\n")

cat("🎯 QUICK PERFORMANCE VALIDATION\n")
cat("Running essential tests for production deployment...\n\n")

# Run basic performance tests (fastest, most essential)
start_time <- Sys.time()

results <- run_master_performance_suite(
  test_levels = c("basic"), 
  output_report = TRUE,
  save_results = TRUE
)

end_time <- Sys.time()
duration <- end_time - start_time

cat("\n" * 2)
cat("=" * 60, "\n")
cat("QUICK PERFORMANCE VALIDATION COMPLETE\n")
cat("=" * 60, "\n")
cat(sprintf("Duration: %.1f minutes\n", as.numeric(duration, units = "mins")))

# Quick summary
readiness <- results$production_readiness
cat(sprintf("Overall Score: %d%%\n", readiness$overall_score))
cat(sprintf("Status: %s\n", readiness$deployment_status))

if (readiness$deployment_status == "APPROVED") {
  cat("\n✅ PRODUCTION READY\n")
  cat("System meets all performance targets for deployment.\n")
  
} else if (readiness$deployment_status == "CONDITIONAL") {
  cat("\n⚠️  CONDITIONAL APPROVAL\n")
  cat("System mostly ready, some optimizations recommended.\n")
  if (length(readiness$recommendations) > 0) {
    cat("Key recommendations:\n")
    for (rec in head(readiness$recommendations, 3)) {
      cat(sprintf("• %s\n", rec))
    }
  }
  
} else {
  cat("\n❌ OPTIMIZATION NEEDED\n")
  cat("Performance issues detected. Address before production deployment.\n")
  if (length(readiness$recommendations) > 0) {
    cat("Critical issues:\n")
    for (rec in readiness$recommendations) {
      cat(sprintf("• %s\n", rec))
    }
  }
}

cat("\n" * 2)
cat("🔬 For comprehensive testing including stress tests:\n")
cat("source('master_performance_test_suite.R')\n")
cat("full_results <- run_master_performance_suite(test_levels = c('basic', 'comprehensive', 'stress'))\n")

cat("\n📊 Results saved to file and available in 'results' variable.\n")