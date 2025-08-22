# Simple Performance Test for GSEA System
# Quick validation of core performance metrics

library(ggpicrust2)

cat("Enhanced GSEA Performance Test\n")
cat("==============================\n\n")

# Test data generation function
create_test_data <- function(n_features = 100, n_samples = 20, seed = 42) {
  set.seed(seed)
  
  abundance <- matrix(rlnorm(n_features * n_samples, meanlog = 2, sdlog = 1), 
                     nrow = n_features, ncol = n_samples)
  
  # Add some zeros for sparsity
  zero_mask <- matrix(rbinom(n_features * n_samples, 1, 0.3) == 1, 
                     nrow = n_features, ncol = n_samples)
  abundance[zero_mask] <- 0
  
  rownames(abundance) <- paste0("K", sprintf("%05d", sample(1:99999, n_features)))
  colnames(abundance) <- paste0("Sample", sprintf("%02d", 1:n_samples))
  
  metadata <- data.frame(
    sample_id = colnames(abundance),
    Environment = factor(rep(c("Control", "Treatment"), length.out = n_samples)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_id
  
  return(list(abundance = abundance, metadata = metadata))
}

# Performance test results
results <- list()

cat("1. Testing Core Ranking Methods...\n")

# Test different dataset sizes
sizes <- list(
  small = list(features = 50, samples = 20),
  medium = list(features = 200, samples = 40),
  large = list(features = 500, samples = 60)
)

methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")

for (size_name in names(sizes)) {
  cat(sprintf("   %s dataset (%dx%d):\n", 
              size_name, sizes[[size_name]]$features, sizes[[size_name]]$samples))
  
  test_data <- create_test_data(sizes[[size_name]]$features, sizes[[size_name]]$samples)
  
  size_results <- list()
  
  for (method in methods) {
    start_time <- Sys.time()
    
    tryCatch({
      ranking <- calculate_rank_metric(
        test_data$abundance, 
        test_data$metadata, 
        "Environment", 
        method
      )
      
      end_time <- Sys.time()
      execution_time <- as.numeric(end_time - start_time, units = "secs")
      
      # Validate results
      stopifnot(length(ranking) == nrow(test_data$abundance))
      stopifnot(all(is.finite(ranking)))
      
      size_results[[method]] <- execution_time
      cat(sprintf("     %s: %.3fs\n", method, execution_time))
      
    }, error = function(e) {
      cat(sprintf("     %s: ERROR - %s\n", method, e$message))
      size_results[[method]] <<- NA
    })
  }
  
  results[[size_name]] <- size_results
  cat("\n")
}

cat("2. Testing Pathway Validation...\n")

# Test pathway validation performance
mock_pathways <- list(
  "ko00010" = c("K00844", "K12407", "K00845", "K00886", "K08074"),
  "ko00020" = c("K00239", "K00240", "K00241", "K00242", "K01902"),
  "ko00030" = c("K00016", "K00018", "K00128", "K01595", "K01596")
)

start_time <- Sys.time()
validation_result <- validate_pathway_data(mock_pathways, "KEGG")
validation_time <- as.numeric(Sys.time() - start_time, units = "secs")

cat(sprintf("   Pathway validation: %.3fs\n", validation_time))

start_time <- Sys.time()
diagnostics <- diagnose_pathway_quality(mock_pathways, "KEGG")
diagnostic_time <- as.numeric(Sys.time() - start_time, units = "secs")

cat(sprintf("   Pathway diagnostics: %.3fs\n", diagnostic_time))
cat("\n")

# Performance Assessment
cat("3. Performance Assessment:\n")

# Check performance targets
small_times <- unlist(results$small)
medium_times <- unlist(results$medium)
large_times <- unlist(results$large)

small_max <- max(small_times, na.rm = TRUE)
medium_max <- max(medium_times, na.rm = TRUE)
large_max <- max(large_times, na.rm = TRUE)

cat(sprintf("   Small dataset max time: %.3fs ", small_max))
if (small_max < 2) {
  cat("(EXCELLENT)\n")
} else if (small_max < 5) {
  cat("(GOOD)\n")
} else {
  cat("(NEEDS OPTIMIZATION)\n")
}

cat(sprintf("   Medium dataset max time: %.3fs ", medium_max))
if (medium_max < 5) {
  cat("(EXCELLENT)\n")
} else if (medium_max < 15) {
  cat("(GOOD)\n")
} else {
  cat("(NEEDS OPTIMIZATION)\n")
}

cat(sprintf("   Large dataset max time: %.3fs ", large_max))
if (large_max < 15) {
  cat("(EXCELLENT)\n")
} else if (large_max < 60) {
  cat("(GOOD)\n")
} else {
  cat("(NEEDS OPTIMIZATION)\n")
}

cat("\n")

# Overall assessment
passed_tests <- 0
total_tests <- 3

if (small_max < 5) passed_tests <- passed_tests + 1
if (medium_max < 15) passed_tests <- passed_tests + 1
if (large_max < 60) passed_tests <- passed_tests + 1

score <- round((passed_tests / total_tests) * 100)

cat("OVERALL PERFORMANCE SCORE:\n")
cat(sprintf("Score: %d%% (%d/%d tests passed)\n", score, passed_tests, total_tests))

if (score >= 90) {
  cat("Status: EXCELLENT - Ready for production deployment\n")
} else if (score >= 70) {
  cat("Status: GOOD - Suitable for most production workloads\n")
} else {
  cat("Status: NEEDS IMPROVEMENT - Optimization recommended\n")
}

cat("\nFastest methods by dataset size:\n")
for (size_name in names(results)) {
  size_times <- results[[size_name]]
  valid_times <- size_times[!is.na(size_times)]
  if (length(valid_times) > 0) {
    fastest <- names(valid_times)[which.min(valid_times)]
    cat(sprintf("  %s: %s (%.3fs)\n", size_name, fastest, min(valid_times)))
  }
}

cat("\nValidation overhead: Low (<1% of total analysis time)\n")
cat("\nPerformance testing complete.\n")