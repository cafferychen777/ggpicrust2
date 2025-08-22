# Master Performance Test Suite for Enhanced GSEA System
#
# Comprehensive performance testing for production deployment and GitHub release
# Integrates all performance benchmarks with detailed reporting
#
# Following Linus principles:
# - "Performance that matters is what users actually experience"
# - "Measure twice, optimize once" 
# - "Real problems require real tests"

library(ggpicrust2)

# =============================================================================
# MASTER TEST ORCHESTRATION
# =============================================================================

#' Execute complete performance test suite
#' 
#' This is the main entry point for comprehensive GSEA system performance testing.
#' It orchestrates all individual benchmark suites and generates a unified report.
run_master_performance_suite <- function(test_levels = c("basic", "comprehensive", "stress"), 
                                        output_report = TRUE,
                                        save_results = TRUE) {
  
  cat("\n")
  cat("◆" * 80, "\n")
  cat("MASTER PERFORMANCE TEST SUITE FOR ENHANCED GSEA SYSTEM\n")
  cat("Production Deployment & GitHub Release Validation\n")
  cat("◆" * 80, "\n")
  cat(sprintf("Started: %s\n", Sys.time()))
  cat(sprintf("Test levels: %s\n", paste(test_levels, collapse = ", ")))
  cat("◆" * 80, "\n\n")
  
  # Initialize master results container
  master_results <- list(
    metadata = list(
      start_time = Sys.time(),
      test_levels = test_levels,
      system_info = get_system_info(),
      package_version = utils::packageVersion("ggpicrust2")
    ),
    benchmarks = list(),
    summary = list()
  )
  
  # ==========================================================================
  # BASIC PERFORMANCE TESTS (Always run)
  # ==========================================================================
  
  if ("basic" %in% test_levels || TRUE) {  # Basic tests always run
    cat("█ BASIC PERFORMANCE TESTS\n")
    cat("─" * 50, "\n\n")
    
    # Core algorithm performance
    cat("1. Core Algorithm Performance...\n")
    source("comprehensive_performance_benchmarks.R")
    master_results$benchmarks$core <- run_comprehensive_benchmarks()
    
    cat("\n2. Essential Workflow Performance...\n")
    source("end_to_end_workflow_benchmarks.R")
    master_results$benchmarks$workflow <- benchmark_kegg_workflow("medium")
    
    cat("✓ Basic performance tests completed\n\n")
  }
  
  # ==========================================================================
  # COMPREHENSIVE PERFORMANCE TESTS
  # ==========================================================================
  
  if ("comprehensive" %in% test_levels) {
    cat("█ COMPREHENSIVE PERFORMANCE TESTS\n")
    cat("─" * 50, "\n\n")
    
    # Full workflow benchmarks
    cat("1. Full Workflow Benchmarks...\n")
    if (!exists("run_workflow_benchmarks")) {
      source("end_to_end_workflow_benchmarks.R")
    }
    master_results$benchmarks$full_workflow <- run_workflow_benchmarks()
    
    # Visualization performance
    cat("2. Visualization Performance...\n")
    source("visualization_performance_benchmarks.R")
    master_results$benchmarks$visualization <- run_visualization_benchmarks()
    
    cat("✓ Comprehensive performance tests completed\n\n")
  }
  
  # ==========================================================================
  # STRESS TESTS (Optional, intensive)
  # ==========================================================================
  
  if ("stress" %in% test_levels) {
    cat("█ STRESS TESTING\n")
    cat("─" * 50, "\n\n")
    
    # Maximum load testing
    cat("1. Maximum Dataset Stress Test...\n")
    stress_results <- run_stress_tests()
    master_results$benchmarks$stress <- stress_results
    
    # Memory leak testing
    cat("2. Memory Leak Detection...\n")
    memory_leak_results <- test_memory_leaks()
    master_results$benchmarks$memory_leaks <- memory_leak_results
    
    # Concurrent load testing
    cat("3. Concurrent Load Testing...\n")
    concurrent_results <- test_concurrent_load()
    master_results$benchmarks$concurrent_load <- concurrent_results
    
    cat("✓ Stress testing completed\n\n")
  }
  
  # ==========================================================================
  # RESULTS ANALYSIS AND REPORTING
  # ==========================================================================
  
  cat("█ RESULTS ANALYSIS\n")
  cat("─" * 50, "\n\n")
  
  master_results$metadata$end_time <- Sys.time()
  master_results$metadata$total_duration <- master_results$metadata$end_time - master_results$metadata$start_time
  
  # Generate comprehensive analysis
  master_results$summary <- analyze_performance_results(master_results$benchmarks)
  
  # Production readiness assessment
  master_results$production_readiness <- assess_production_readiness(master_results$summary)
  
  # Generate detailed report
  if (output_report) {
    generate_master_performance_report(master_results)
  }
  
  # Save results if requested
  if (save_results) {
    save_performance_results(master_results)
  }
  
  cat("\n")
  cat("◆" * 80, "\n")
  cat(sprintf("MASTER PERFORMANCE SUITE COMPLETED: %s\n", Sys.time()))
  cat(sprintf("Total Duration: %.1f minutes\n", 
              as.numeric(master_results$metadata$total_duration, units = "mins")))
  cat("◆" * 80, "\n")
  
  return(master_results)
}

# =============================================================================
# SYSTEM INFORMATION COLLECTION
# =============================================================================

#' Collect system information for performance context
get_system_info <- function() {
  list(
    r_version = R.version.string,
    platform = R.version$platform,
    os = Sys.info()["sysname"],
    cpu_cores = parallel::detectCores(),
    memory_gb = round(as.numeric(gsub("[^0-9]", "", system("sysctl hw.memsize", intern = TRUE))) / 1024^3, 1),
    timestamp = Sys.time(),
    locale = Sys.getlocale(),
    packages = list(
      ggplot2 = utils::packageVersion("ggplot2"),
      dplyr = if (requireNamespace("dplyr", quietly = TRUE)) utils::packageVersion("dplyr") else "Not installed",
      fgsea = if (requireNamespace("fgsea", quietly = TRUE)) utils::packageVersion("fgsea") else "Not installed"
    )
  )
}

# =============================================================================
# STRESS TESTING FUNCTIONS
# =============================================================================

#' Run maximum load stress tests
run_stress_tests <- function() {
  cat("  Maximum dataset size testing...\n")
  
  stress_configs <- list(
    extreme = list(features = 20000, samples = 500, sparsity = 0.95),
    maximum = list(features = 50000, samples = 1000, sparsity = 0.98)
  )
  
  results <- list()
  
  for (config_name in names(stress_configs)) {
    config <- stress_configs[[config_name]]
    cat(sprintf("    Testing %s configuration (%d×%d)...\n", 
                config_name, config$features, config$samples))
    
    tryCatch({
      start_time <- Sys.time()
      
      # Create massive dataset
      stress_data <- create_benchmark_data(config$features, config$samples, config$sparsity)
      
      # Test ranking calculation
      ranking <- calculate_rank_metric(stress_data$abundance, stress_data$metadata,
                                     "Environment", "signal2noise")
      
      end_time <- Sys.time()
      execution_time <- as.numeric(end_time - start_time, units = "secs")
      
      results[[config_name]] <- list(
        success = TRUE,
        time = execution_time,
        features = config$features,
        samples = config$samples,
        data_size_mb = object.size(stress_data$abundance) / 1024^2
      )
      
      cat(sprintf("    ✓ Completed in %.1f seconds\n", execution_time))
      
      # Cleanup
      rm(stress_data, ranking)
      gc()
      
    }, error = function(e) {
      results[[config_name]] <<- list(
        success = FALSE,
        error = e$message,
        features = config$features,
        samples = config$samples
      )
      cat(sprintf("    ✗ Failed: %s\n", e$message))
    })
  }
  
  return(results)
}

#' Test for memory leaks in repeated operations
test_memory_leaks <- function() {
  cat("  Memory leak detection over repeated operations...\n")
  
  # Create baseline dataset
  test_data <- create_benchmark_data(500, 50)
  
  # Record initial memory
  gc()
  initial_memory <- sum(gc()[, 2])
  
  memory_timeline <- numeric(20)
  
  # Run multiple iterations
  for (i in 1:20) {
    # Perform analysis
    ranking <- calculate_rank_metric(test_data$abundance, test_data$metadata,
                                   "Environment", "signal2noise")
    
    # Force cleanup
    rm(ranking)
    gc()
    
    # Record memory usage
    memory_timeline[i] <- sum(gc()[, 2])
  }
  
  # Analyze memory trend
  final_memory <- tail(memory_timeline, 1)
  memory_increase <- final_memory - initial_memory
  
  # Linear regression to detect trend
  memory_trend <- lm(memory_timeline ~ seq_along(memory_timeline))
  trend_slope <- coef(memory_trend)[2]
  
  cat(sprintf("    Memory change: %.1f MB\n", memory_increase))
  cat(sprintf("    Trend slope: %.3f MB/iteration\n", trend_slope))
  
  leak_detected <- trend_slope > 0.1  # More than 0.1 MB increase per iteration
  
  if (leak_detected) {
    cat("    ⚠ Potential memory leak detected\n")
  } else {
    cat("    ✓ No memory leaks detected\n")
  }
  
  return(list(
    leak_detected = leak_detected,
    memory_increase = memory_increase,
    trend_slope = trend_slope,
    memory_timeline = memory_timeline
  ))
}

#' Test concurrent load handling
test_concurrent_load <- function() {
  cat("  Concurrent analysis load testing...\n")
  
  # Create multiple datasets for concurrent processing
  datasets <- list()
  for (i in 1:5) {
    datasets[[i]] <- create_benchmark_data(300, 40, seed = i * 111)
  }
  
  # Sequential baseline
  start_time <- Sys.time()
  sequential_results <- list()
  for (i in 1:5) {
    sequential_results[[i]] <- calculate_rank_metric(
      datasets[[i]]$abundance, datasets[[i]]$metadata, "Environment", "signal2noise"
    )
  }
  sequential_time <- as.numeric(Sys.time() - start_time, units = "secs")
  
  cat(sprintf("    Sequential processing: %.1f seconds\n", sequential_time))
  
  # Test data isolation (simulated concurrent access)
  # Modify one dataset and ensure others are unaffected
  original_value <- datasets[[1]]$abundance[1, 1]
  datasets[[1]]$abundance[1, 1] <- original_value * 10
  
  # Recalculate
  modified_result <- calculate_rank_metric(
    datasets[[1]]$abundance, datasets[[1]]$metadata, "Environment", "signal2noise"
  )
  unchanged_result <- calculate_rank_metric(
    datasets[[2]]$abundance, datasets[[2]]$metadata, "Environment", "signal2noise"
  )
  
  # Check isolation
  isolation_maintained <- identical(sequential_results[[2]], unchanged_result)
  modification_detected <- !identical(sequential_results[[1]], modified_result)
  
  cat(sprintf("    Data isolation: %s\n", if (isolation_maintained) "✓ PASSED" else "✗ FAILED"))
  cat(sprintf("    Modification detection: %s\n", if (modification_detected) "✓ PASSED" else "✗ FAILED"))
  
  return(list(
    sequential_time = sequential_time,
    isolation_maintained = isolation_maintained,
    modification_detected = modification_detected,
    datasets_tested = length(datasets)
  ))
}

# =============================================================================
# RESULTS ANALYSIS
# =============================================================================

#' Analyze performance results across all benchmarks
analyze_performance_results <- function(benchmark_results) {
  
  analysis <- list()
  
  # Core performance analysis
  if (!is.null(benchmark_results$core)) {
    analysis$core_performance <- analyze_core_performance(benchmark_results$core)
  }
  
  # Workflow performance analysis
  if (!is.null(benchmark_results$workflow) || !is.null(benchmark_results$full_workflow)) {
    workflow_data <- benchmark_results$full_workflow %||% benchmark_results$workflow
    analysis$workflow_performance <- analyze_workflow_performance(workflow_data)
  }
  
  # Visualization performance analysis
  if (!is.null(benchmark_results$visualization)) {
    analysis$visualization_performance <- analyze_visualization_performance(benchmark_results$visualization)
  }
  
  # Stress testing analysis
  if (!is.null(benchmark_results$stress)) {
    analysis$stress_performance <- analyze_stress_performance(benchmark_results$stress)
  }
  
  return(analysis)
}

#' Analyze core algorithm performance
analyze_core_performance <- function(core_results) {
  
  # Extract ranking method performance across dataset sizes
  ranking_results <- core_results$ranking
  
  methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  sizes <- names(ranking_results)
  
  # Performance matrix: methods × sizes
  perf_matrix <- matrix(NA, nrow = length(methods), ncol = length(sizes),
                       dimnames = list(methods, sizes))
  
  for (size in sizes) {
    for (method in methods) {
      if (!is.null(ranking_results[[size]][[method]])) {
        perf_matrix[method, size] <- ranking_results[[size]][[method]]$median_time
      }
    }
  }
  
  # Find fastest/slowest methods
  fastest_method <- methods[which.min(rowMeans(perf_matrix, na.rm = TRUE))]
  slowest_method <- methods[which.max(rowMeans(perf_matrix, na.rm = TRUE))]
  
  # Scaling analysis
  scaling_coefficients <- numeric(length(methods))
  names(scaling_coefficients) <- methods
  
  for (method in methods) {
    method_times <- perf_matrix[method, ]
    if (sum(!is.na(method_times)) >= 3) {
      # Extract data sizes
      data_sizes <- sapply(sizes, function(s) {
        config <- BENCHMARK_CONFIGS[[s]]
        config$features * config$samples
      })
      
      # Fit scaling model
      valid_idx <- !is.na(method_times)
      if (sum(valid_idx) >= 3) {
        scaling_model <- lm(log(method_times[valid_idx]) ~ log(data_sizes[valid_idx]))
        scaling_coefficients[method] <- coef(scaling_model)[2]
      }
    }
  }
  
  return(list(
    performance_matrix = perf_matrix,
    fastest_method = fastest_method,
    slowest_method = slowest_method,
    scaling_coefficients = scaling_coefficients,
    average_scaling = mean(scaling_coefficients, na.rm = TRUE)
  ))
}

#' Analyze workflow performance
analyze_workflow_performance <- function(workflow_results) {
  
  # Extract workflow timings if available
  if (is.list(workflow_results) && "times" %in% names(workflow_results)) {
    # Single workflow result
    times <- workflow_results$times
    
    total_time <- times$total
    steps <- names(times)[names(times) != "total"]
    
    # Calculate step proportions
    step_proportions <- sapply(steps, function(step) times[[step]] / total_time)
    bottleneck_step <- names(step_proportions)[which.max(step_proportions)]
    
    return(list(
      total_time = total_time,
      step_times = times[steps],
      step_proportions = step_proportions,
      bottleneck_step = bottleneck_step,
      efficiency_score = calculate_efficiency_score(total_time)
    ))
    
  } else {
    # Multiple workflow results
    return(list(
      message = "Multiple workflow results - detailed analysis requires individual workflow data"
    ))
  }
}

#' Analyze visualization performance  
analyze_visualization_performance <- function(viz_results) {
  
  max_times <- list()
  
  # Enrichment plots
  if (!is.null(viz_results$enrichment)) {
    enrich_times <- unlist(lapply(viz_results$enrichment, function(x) lapply(x, function(y) y$median_time)))
    max_times$enrichment <- max(enrich_times, na.rm = TRUE)
  }
  
  # Heatmaps
  if (!is.null(viz_results$heatmap)) {
    heatmap_times <- unlist(lapply(viz_results$heatmap, function(x) lapply(x, function(y) y$time)))
    max_times$heatmap <- max(heatmap_times, na.rm = TRUE)
  }
  
  # Network plots
  if (!is.null(viz_results$network)) {
    network_times <- sapply(viz_results$network, function(x) x$time)
    max_times$network <- max(network_times, na.rm = TRUE)
  }
  
  # Interactive plots
  if (!is.null(viz_results$interactive)) {
    interactive_times <- sapply(viz_results$interactive, function(x) x$time)
    max_times$interactive <- max(interactive_times, na.rm = TRUE)
  }
  
  return(list(
    max_times = max_times,
    overall_max = max(unlist(max_times), na.rm = TRUE),
    responsive_threshold = 5.0  # seconds
  ))
}

#' Analyze stress test performance
analyze_stress_performance <- function(stress_results) {
  
  analysis <- list()
  
  # Maximum load analysis
  if (!is.null(stress_results$extreme) || !is.null(stress_results$maximum)) {
    max_success <- any(sapply(stress_results, function(x) x$success %||% FALSE))
    analysis$max_load_handled = max_success
    
    if (max_success) {
      successful_tests <- stress_results[sapply(stress_results, function(x) x$success %||% FALSE)]
      max_features <- max(sapply(successful_tests, function(x) x$features))
      max_samples <- max(sapply(successful_tests, function(x) x$samples))
      
      analysis$max_dataset_size <- paste(max_features, "×", max_samples)
    }
  }
  
  # Memory leak analysis
  if (!is.null(stress_results$memory_leaks)) {
    analysis$memory_leak_status <- if (stress_results$memory_leaks$leak_detected) "DETECTED" else "NONE"
  }
  
  # Concurrent load analysis
  if (!is.null(stress_results$concurrent_load)) {
    analysis$concurrent_capability <- stress_results$concurrent_load$isolation_maintained && 
                                    stress_results$concurrent_load$modification_detected
  }
  
  return(analysis)
}

# =============================================================================
# PRODUCTION READINESS ASSESSMENT
# =============================================================================

#' Assess production deployment readiness
assess_production_readiness <- function(performance_summary) {
  
  readiness <- list(
    overall_score = 0,
    criteria = list(),
    recommendations = list(),
    deployment_status = "UNKNOWN"
  )
  
  total_criteria <- 0
  passed_criteria <- 0
  
  # Core Performance Criteria
  if (!is.null(performance_summary$core_performance)) {
    total_criteria <- total_criteria + 3
    
    core <- performance_summary$core_performance
    
    # Criterion 1: Small dataset speed (<5s)
    small_time <- core$performance_matrix[, "small"]
    if (any(!is.na(small_time)) && max(small_time, na.rm = TRUE) < 5) {
      passed_criteria <- passed_criteria + 1
      readiness$criteria$small_dataset_speed <- "PASS"
    } else {
      readiness$criteria$small_dataset_speed <- "FAIL"
      readiness$recommendations <- append(readiness$recommendations, "Optimize small dataset performance")
    }
    
    # Criterion 2: Large dataset speed (<120s)
    if ("large" %in% colnames(core$performance_matrix)) {
      large_time <- core$performance_matrix[, "large"]
      if (any(!is.na(large_time)) && max(large_time, na.rm = TRUE) < 120) {
        passed_criteria <- passed_criteria + 1
        readiness$criteria$large_dataset_speed <- "PASS"
      } else {
        readiness$criteria$large_dataset_speed <- "FAIL"
        readiness$recommendations <- append(readiness$recommendations, "Optimize large dataset performance")
      }
    }
    
    # Criterion 3: Reasonable scaling (exponent < 2.0)
    if (!is.na(core$average_scaling) && core$average_scaling < 2.0) {
      passed_criteria <- passed_criteria + 1
      readiness$criteria$scaling_efficiency <- "PASS"
    } else {
      readiness$criteria$scaling_efficiency <- "FAIL"
      readiness$recommendations <- append(readiness$recommendations, "Improve algorithmic scaling")
    }
  }
  
  # Workflow Performance Criteria
  if (!is.null(performance_summary$workflow_performance)) {
    total_criteria <- total_criteria + 1
    
    workflow <- performance_summary$workflow_performance
    
    # Criterion: Complete workflow under 2 minutes
    if (!is.null(workflow$total_time) && workflow$total_time < 120) {
      passed_criteria <- passed_criteria + 1
      readiness$criteria$workflow_speed <- "PASS"
    } else {
      readiness$criteria$workflow_speed <- "FAIL"
      readiness$recommendations <- append(readiness$recommendations, "Optimize workflow bottlenecks")
    }
  }
  
  # Visualization Performance Criteria
  if (!is.null(performance_summary$visualization_performance)) {
    total_criteria <- total_criteria + 1
    
    viz <- performance_summary$visualization_performance
    
    # Criterion: All visualizations under 5s
    if (!is.null(viz$overall_max) && viz$overall_max < 5) {
      passed_criteria <- passed_criteria + 1
      readiness$criteria$visualization_responsiveness <- "PASS"
    } else {
      readiness$criteria$visualization_responsiveness <- "FAIL"
      readiness$recommendations <- append(readiness$recommendations, "Optimize slow visualizations")
    }
  }
  
  # Stress Testing Criteria
  if (!is.null(performance_summary$stress_performance)) {
    total_criteria <- total_criteria + 2
    
    stress <- performance_summary$stress_performance
    
    # Criterion: Handle large datasets
    if (!is.null(stress$max_load_handled) && stress$max_load_handled) {
      passed_criteria <- passed_criteria + 1
      readiness$criteria$stress_handling <- "PASS"
    } else {
      readiness$criteria$stress_handling <- "FAIL"
      readiness$recommendations <- append(readiness$recommendations, "Improve large dataset handling")
    }
    
    # Criterion: No memory leaks
    if (!is.null(stress$memory_leak_status) && stress$memory_leak_status == "NONE") {
      passed_criteria <- passed_criteria + 1
      readiness$criteria$memory_stability <- "PASS"
    } else {
      readiness$criteria$memory_stability <- "FAIL"
      readiness$recommendations <- append(readiness$recommendations, "Fix memory leaks")
    }
  }
  
  # Calculate overall score
  if (total_criteria > 0) {
    readiness$overall_score <- round((passed_criteria / total_criteria) * 100)
  }
  
  # Determine deployment status
  if (readiness$overall_score >= 90) {
    readiness$deployment_status <- "APPROVED"
  } else if (readiness$overall_score >= 75) {
    readiness$deployment_status <- "CONDITIONAL"
  } else {
    readiness$deployment_status <- "NOT_READY"
  }
  
  return(readiness)
}

# =============================================================================
# REPORTING FUNCTIONS
# =============================================================================

#' Generate comprehensive master performance report
generate_master_performance_report <- function(results) {
  cat("\n")
  cat("◆" * 80, "\n")
  cat("MASTER PERFORMANCE REPORT\n")
  cat("Enhanced GSEA System - Production Deployment Assessment\n")
  cat("◆" * 80, "\n\n")
  
  # System Information
  cat("SYSTEM CONFIGURATION:\n")
  cat("─" * 30, "\n")
  sys_info <- results$metadata$system_info
  cat(sprintf("• R Version: %s\n", sys_info$r_version))
  cat(sprintf("• Platform: %s\n", sys_info$platform))
  cat(sprintf("• CPU Cores: %d\n", sys_info$cpu_cores))
  cat(sprintf("• Memory: %.1f GB\n", sys_info$memory_gb))
  cat(sprintf("• Package Version: %s\n", results$metadata$package_version))
  cat(sprintf("• Test Duration: %.1f minutes\n", as.numeric(results$metadata$total_duration, units = "mins")))
  cat("\n")
  
  # Production Readiness Assessment
  readiness <- results$production_readiness
  cat("PRODUCTION READINESS ASSESSMENT:\n")
  cat("─" * 40, "\n")
  cat(sprintf("Overall Score: %d%%\n", readiness$overall_score))
  cat(sprintf("Deployment Status: %s\n", readiness$deployment_status))
  cat("\n")
  
  # Detailed criteria
  cat("Performance Criteria:\n")
  for (criterion in names(readiness$criteria)) {
    status <- readiness$criteria[[criterion]]
    symbol <- if (status == "PASS") "✓" else "✗"
    cat(sprintf("  %s %s: %s\n", symbol, gsub("_", " ", toupper(criterion)), status))
  }
  cat("\n")
  
  # Recommendations
  if (length(readiness$recommendations) > 0) {
    cat("RECOMMENDATIONS:\n")
    cat("─" * 20, "\n")
    for (i in seq_along(readiness$recommendations)) {
      cat(sprintf("%d. %s\n", i, readiness$recommendations[[i]]))
    }
    cat("\n")
  }
  
  # Performance Summary by Component
  summary <- results$summary
  
  if (!is.null(summary$core_performance)) {
    cat("CORE ALGORITHM PERFORMANCE:\n")
    cat("─" * 35, "\n")
    core <- summary$core_performance
    cat(sprintf("• Fastest Method: %s\n", core$fastest_method))
    cat(sprintf("• Slowest Method: %s\n", core$slowest_method))
    cat(sprintf("• Average Scaling: O(n^%.2f)\n", core$average_scaling))
    cat("\n")
  }
  
  if (!is.null(summary$workflow_performance)) {
    cat("WORKFLOW PERFORMANCE:\n")
    cat("─" * 25, "\n")
    workflow <- summary$workflow_performance
    if (!is.null(workflow$total_time)) {
      cat(sprintf("• Total Workflow Time: %.1fs\n", workflow$total_time))
      cat(sprintf("• Primary Bottleneck: %s\n", workflow$bottleneck_step))
      cat(sprintf("• Efficiency Score: %.1f/10\n", workflow$efficiency_score))
    }
    cat("\n")
  }
  
  if (!is.null(summary$visualization_performance)) {
    cat("VISUALIZATION PERFORMANCE:\n")
    cat("─" * 30, "\n")
    viz <- summary$visualization_performance
    cat(sprintf("• Maximum Render Time: %.2fs\n", viz$overall_max))
    cat(sprintf("• Responsiveness Target: %.1fs\n", viz$responsive_threshold))
    status <- if (viz$overall_max < viz$responsive_threshold) "MEETS TARGET" else "EXCEEDS TARGET"
    cat(sprintf("• Status: %s\n", status))
    cat("\n")
  }
  
  # Final Assessment
  cat("FINAL ASSESSMENT:\n")
  cat("─" * 20, "\n")
  
  if (readiness$deployment_status == "APPROVED") {
    cat("✅ APPROVED FOR PRODUCTION DEPLOYMENT\n")
    cat("   • All critical performance criteria met\n")
    cat("   • System ready for GitHub release\n")
    cat("   • Can handle expected production workloads\n")
    
  } else if (readiness$deployment_status == "CONDITIONAL") {
    cat("⚠️  CONDITIONAL APPROVAL\n")
    cat("   • Most performance criteria met\n")
    cat("   • Some optimizations recommended before release\n")
    cat("   • Monitor performance in production environment\n")
    
  } else {
    cat("❌ NOT READY FOR PRODUCTION\n")
    cat("   • Critical performance issues identified\n")
    cat("   • Optimization required before deployment\n")
    cat("   • Re-test after improvements\n")
  }
  
  cat("\n")
  cat("◆" * 80, "\n")
}

#' Save performance results to file
save_performance_results <- function(results, filename = NULL) {
  if (is.null(filename)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- sprintf("gsea_performance_results_%s.RData", timestamp)
  }
  
  save(results, file = filename)
  cat(sprintf("Performance results saved to: %s\n", filename))
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Calculate efficiency score for workflow performance
calculate_efficiency_score <- function(total_time) {
  # Score based on total time (10 = excellent, 1 = poor)
  if (total_time < 30) return(10)
  if (total_time < 60) return(8)
  if (total_time < 120) return(6)
  if (total_time < 300) return(4)
  if (total_time < 600) return(2)
  return(1)
}

#' Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if (interactive() || !exists("master_performance_results")) {
  cat("MASTER PERFORMANCE TEST SUITE\n")
  cat("============================\n\n")
  
  cat("This comprehensive test suite will evaluate:\n")
  cat("• Core algorithm performance across dataset sizes\n")
  cat("• End-to-end workflow efficiency\n")
  cat("• Visualization rendering performance\n")
  cat("• System stress handling and memory management\n")
  cat("• Production deployment readiness\n\n")
  
  cat("Test levels available:\n")
  cat("• 'basic': Essential performance tests (~5-10 minutes)\n")
  cat("• 'comprehensive': Full benchmark suite (~15-20 minutes)\n") 
  cat("• 'stress': Maximum load and stress tests (~10-15 minutes)\n\n")
  
  # Run with default comprehensive testing
  master_performance_results <- run_master_performance_suite(
    test_levels = c("basic", "comprehensive"),
    output_report = TRUE,
    save_results = TRUE
  )
  
  cat("\nMaster performance results available in 'master_performance_results' variable.\n")
  cat("Detailed report generated above.\n")
  cat("Results saved to timestamped file.\n")
}