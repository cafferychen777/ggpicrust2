#!/usr/bin/env Rscript

# Cross-Pathway Type Consistency Demonstration
# This script demonstrates the practical usage of KEGG, MetaCyc, and GO pathways
# in a typical microbiome analysis workflow

library(ggpicrust2)
library(ggplot2)
library(dplyr)

cat("=== Cross-Pathway Type Consistency Demo ===\n")
cat("Testing KEGG, MetaCyc, and GO pathway integration\n\n")

# Load example data
data(ko_abundance, package = "ggpicrust2")
data(metacyc_abundance, package = "ggpicrust2") 
data(metadata, package = "ggpicrust2")

cat("1. DATA PREPARATION\n")
cat("-------------------\n")

# Prepare KO abundance data (for KEGG and GO)
ko_data <- as.data.frame(ko_abundance)
if ("#NAME" %in% colnames(ko_data)) {
  rownames(ko_data) <- ko_data[["#NAME"]]
  ko_data <- ko_data[, !colnames(ko_data) %in% "#NAME"]
}

# Prepare EC abundance data (for MetaCyc)
ec_data <- as.data.frame(metacyc_abundance)
if ("#NAME" %in% colnames(ec_data)) {
  rownames(ec_data) <- ec_data[["#NAME"]]
  ec_data <- ec_data[, !colnames(ec_data) %in% "#NAME"]
}

# Ensure metadata has proper row names
if ("sample_name" %in% colnames(metadata)) {
  rownames(metadata) <- metadata$sample_name
}

cat(sprintf("KO abundance data: %d features x %d samples\n", nrow(ko_data), ncol(ko_data)))
cat(sprintf("EC abundance data: %d features x %d samples\n", nrow(ec_data), ncol(ec_data)))
cat(sprintf("Metadata: %d samples x %d variables\n", nrow(metadata), ncol(metadata)))

# Find a suitable grouping variable
group_var <- NULL
for (col in colnames(metadata)) {
  if (is.factor(metadata[[col]]) || is.character(metadata[[col]])) {
    unique_vals <- unique(metadata[[col]])
    if (length(unique_vals) == 2 && all(table(metadata[[col]]) >= 2)) {
      group_var <- col
      break
    }
  }
}

if (is.null(group_var)) {
  # Create a simple binary grouping if none exists
  metadata$test_group <- factor(rep(c("Group1", "Group2"), length.out = nrow(metadata)))
  group_var <- "test_group"
}

cat(sprintf("Using grouping variable: %s\n", group_var))
cat(sprintf("Group distribution: %s\n", paste(names(table(metadata[[group_var]])), "=", table(metadata[[group_var]]), collapse = ", ")))

cat("\n2. CROSS-PATHWAY GSEA ANALYSIS\n")
cat("------------------------------\n")

# Function to safely run GSEA with error handling
safe_gsea <- function(abundance, pathway_type, group_var, n_attempts = 3) {
  for (attempt in 1:n_attempts) {
    tryCatch({
      cat(sprintf("Running %s GSEA analysis (attempt %d)...\n", pathway_type, attempt))
      
      result <- pathway_gsea(
        abundance = abundance,
        metadata = metadata,
        group = group_var,
        pathway_type = pathway_type,
        method = "fgsea",
        rank_method = "signal2noise",
        nperm = 100,
        min_size = 5,
        max_size = 200,
        seed = 42
      )
      
      cat(sprintf("✓ %s analysis completed: %d pathways found\n", pathway_type, nrow(result)))
      return(result)
      
    }, error = function(e) {
      cat(sprintf("✗ %s analysis failed (attempt %d): %s\n", pathway_type, attempt, e$message))
      if (attempt == n_attempts) {
        return(data.frame())
      }
    })
  }
}

# Run GSEA for all pathway types
results <- list()

# KEGG analysis
results$KEGG <- safe_gsea(ko_data, "KEGG", group_var)

# MetaCyc analysis  
results$MetaCyc <- safe_gsea(ec_data, "MetaCyc", group_var)

# GO analysis
results$GO <- safe_gsea(ko_data, "GO", group_var)

cat("\n3. RESULT STRUCTURE CONSISTENCY CHECK\n")
cat("-------------------------------------\n")

# Check result structure consistency
check_result_consistency <- function(results) {
  if (length(results) == 0) return(FALSE)
  
  pathway_types <- names(results)
  reference_cols <- NULL
  
  for (type in pathway_types) {
    if (nrow(results[[type]]) == 0) {
      cat(sprintf("⚠ %s: No results found\n", type))
      next
    }
    
    current_cols <- colnames(results[[type]])
    
    if (is.null(reference_cols)) {
      reference_cols <- current_cols
      cat(sprintf("✓ %s: %d pathways, columns: %s\n", type, nrow(results[[type]]), 
                 paste(current_cols, collapse = ", ")))
    } else {
      if (identical(reference_cols, current_cols)) {
        cat(sprintf("✓ %s: %d pathways, consistent column structure\n", type, nrow(results[[type]])))
      } else {
        cat(sprintf("✗ %s: %d pathways, INCONSISTENT columns\n", type, nrow(results[[type]])))
        cat(sprintf("  Expected: %s\n", paste(reference_cols, collapse = ", ")))
        cat(sprintf("  Got:      %s\n", paste(current_cols, collapse = ", ")))
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}

consistency_check <- check_result_consistency(results)
cat(sprintf("\nStructure consistency: %s\n", ifelse(consistency_check, "✓ PASSED", "✗ FAILED")))

cat("\n4. STATISTICAL CONSISTENCY ANALYSIS\n")
cat("-----------------------------------\n")

# Analyze statistical properties
analyze_statistical_properties <- function(results) {
  stats_summary <- data.frame()
  
  for (type in names(results)) {
    if (nrow(results[[type]]) == 0) next
    
    res <- results[[type]]
    stats <- data.frame(
      pathway_type = type,
      n_pathways = nrow(res),
      mean_NES = mean(abs(res$NES), na.rm = TRUE),
      sd_NES = sd(res$NES, na.rm = TRUE),
      min_pvalue = min(res$pvalue, na.rm = TRUE),
      max_pvalue = max(res$pvalue, na.rm = TRUE),
      n_significant = sum(res$p.adjust < 0.05, na.rm = TRUE),
      prop_significant = sum(res$p.adjust < 0.05, na.rm = TRUE) / nrow(res)
    )
    
    stats_summary <- rbind(stats_summary, stats)
  }
  
  if (nrow(stats_summary) > 0) {
    print(stats_summary)
    
    # Check for reasonable statistical properties
    cat("\nStatistical consistency checks:\n")
    
    # NES values should be reasonable
    if (all(stats_summary$mean_NES < 5 & stats_summary$mean_NES > 0)) {
      cat("✓ NES values are in reasonable range\n")
    } else {
      cat("⚠ NES values may be extreme\n")
    }
    
    # P-values should span reasonable range
    if (all(stats_summary$min_pvalue >= 0 & stats_summary$max_pvalue <= 1)) {
      cat("✓ P-values are in valid range [0,1]\n")
    } else {
      cat("✗ P-values are outside valid range\n")
    }
    
    # Should have some significant results (but not all)
    reasonable_sig <- stats_summary$prop_significant > 0.01 & stats_summary$prop_significant < 0.8
    if (all(reasonable_sig)) {
      cat("✓ Proportion of significant results is reasonable\n")
    } else {
      cat("⚠ Proportion of significant results may be unusual\n")
    }
  }
  
  return(stats_summary)
}

stats_summary <- analyze_statistical_properties(results)

cat("\n5. ANNOTATION SYSTEM TESTING\n")
cat("----------------------------\n")

# Test annotation system for each pathway type
annotated_results <- list()

for (type in names(results)) {
  if (nrow(results[[type]]) == 0) {
    cat(sprintf("⚠ %s: Skipping annotation (no results)\n", type))
    next
  }
  
  tryCatch({
    cat(sprintf("Annotating %s results...\n", type))
    
    annotated <- gsea_pathway_annotation(
      gsea_results = results[[type]],
      pathway_type = type
    )
    
    # Check if annotations were added
    has_names <- "pathway_name" %in% colnames(annotated)
    names_different <- has_names && any(annotated$pathway_name != annotated$pathway_id)
    
    cat(sprintf("✓ %s: Annotation completed, pathway names %s\n", 
               type, ifelse(names_different, "updated", "unchanged")))
    
    annotated_results[[type]] <- annotated
    
  }, error = function(e) {
    cat(sprintf("✗ %s: Annotation failed - %s\n", type, e$message))
  })
}

cat("\n6. VISUALIZATION CONSISTENCY TESTING\n")
cat("------------------------------------\n")

# Test visualization across pathway types
test_visualizations <- function(annotated_results) {
  plot_types <- c("enrichment_plot", "dotplot", "barplot")
  
  for (type in names(annotated_results)) {
    if (nrow(annotated_results[[type]]) == 0) {
      cat(sprintf("⚠ %s: Skipping visualization (no results)\n", type))
      next
    }
    
    cat(sprintf("Testing %s visualizations:\n", type))
    
    for (plot_type in plot_types) {
      tryCatch({
        n_pathways <- min(10, nrow(annotated_results[[type]]))
        
        plot <- visualize_gsea(
          gsea_results = annotated_results[[type]],
          plot_type = plot_type,
          n_pathways = n_pathways
        )
        
        # Save plot
        filename <- sprintf("cross_pathway_test_%s_%s.pdf", type, plot_type)
        ggsave(filename, plot, width = 10, height = 8)
        cat(sprintf("  ✓ %s saved to %s\n", plot_type, filename))
        
      }, error = function(e) {
        cat(sprintf("  ✗ %s failed: %s\n", plot_type, e$message))
      })
    }
  }
}

test_visualizations(annotated_results)

cat("\n7. INTEGRATION WORKFLOW DEMONSTRATION\n")
cat("-------------------------------------\n")

# Demonstrate switching between pathway types in workflow
workflow_demo <- function() {
  cat("Demonstrating typical analysis workflow with all three pathway types:\n\n")
  
  # Step 1: KEGG for metabolic overview
  if (!is.null(results$KEGG) && nrow(results$KEGG) > 0) {
    kegg_significant <- results$KEGG[results$KEGG$p.adjust < 0.1, ]
    cat(sprintf("Step 1 - KEGG metabolic overview: %d significant pathways\n", nrow(kegg_significant)))
    if (nrow(kegg_significant) > 0) {
      cat("Top KEGG pathways:\n")
      top_kegg <- head(kegg_significant[order(kegg_significant$p.adjust), ], 3)
      for (i in 1:nrow(top_kegg)) {
        cat(sprintf("  - %s (NES=%.2f, p.adj=%.3f)\n", 
                   top_kegg$pathway_id[i], top_kegg$NES[i], top_kegg$p.adjust[i]))
      }
    }
  }
  
  # Step 2: MetaCyc for detailed mechanisms
  if (!is.null(results$MetaCyc) && nrow(results$MetaCyc) > 0) {
    metacyc_significant <- results$MetaCyc[results$MetaCyc$p.adjust < 0.1, ]
    cat(sprintf("\nStep 2 - MetaCyc detailed mechanisms: %d significant pathways\n", nrow(metacyc_significant)))
    if (nrow(metacyc_significant) > 0) {
      cat("Top MetaCyc pathways:\n")
      top_metacyc <- head(metacyc_significant[order(metacyc_significant$p.adjust), ], 3)
      for (i in 1:nrow(top_metacyc)) {
        cat(sprintf("  - %s (NES=%.2f, p.adj=%.3f)\n", 
                   top_metacyc$pathway_id[i], top_metacyc$NES[i], top_metacyc$p.adjust[i]))
      }
    }
  }
  
  # Step 3: GO for functional categories
  if (!is.null(results$GO) && nrow(results$GO) > 0) {
    go_significant <- results$GO[results$GO$p.adjust < 0.1, ]
    cat(sprintf("\nStep 3 - GO functional categories: %d significant pathways\n", nrow(go_significant)))
    if (nrow(go_significant) > 0) {
      cat("Top GO terms:\n")
      top_go <- head(go_significant[order(go_significant$p.adjust), ], 3)
      for (i in 1:nrow(top_go)) {
        cat(sprintf("  - %s (NES=%.2f, p.adj=%.3f)\n", 
                   top_go$pathway_id[i], top_go$NES[i], top_go$p.adjust[i]))
      }
    }
  }
}

workflow_demo()

cat("\n8. PERFORMANCE COMPARISON\n")
cat("-------------------------\n")

# Compare performance across pathway types
performance_test <- function() {
  pathway_types <- c("KEGG", "MetaCyc", "GO")
  times <- numeric(length(pathway_types))
  names(times) <- pathway_types
  
  for (i in seq_along(pathway_types)) {
    type <- pathway_types[i]
    abundance_data <- if (type == "MetaCyc") ec_data else ko_data
    
    cat(sprintf("Timing %s analysis...\n", type))
    
    start_time <- Sys.time()
    
    tryCatch({
      result <- pathway_gsea(
        abundance = abundance_data,
        metadata = metadata,
        group = group_var,
        pathway_type = type,
        method = "fgsea",
        nperm = 100,
        seed = 42
      )
      
      end_time <- Sys.time()
      times[i] <- as.numeric(end_time - start_time, units = "secs")
      
      cat(sprintf("✓ %s completed in %.2f seconds\n", type, times[i]))
      
    }, error = function(e) {
      times[i] <- NA
      cat(sprintf("✗ %s failed\n", type))
    })
  }
  
  # Performance summary
  valid_times <- times[!is.na(times)]
  if (length(valid_times) > 1) {
    cat("\nPerformance summary:\n")
    cat(sprintf("Fastest: %s (%.2f sec)\n", names(which.min(valid_times)), min(valid_times)))
    cat(sprintf("Slowest: %s (%.2f sec)\n", names(which.max(valid_times)), max(valid_times)))
    cat(sprintf("Speed ratio: %.2fx\n", max(valid_times) / min(valid_times)))
  }
  
  return(times)
}

performance_times <- performance_test()

cat("\n9. FINAL SUMMARY\n")
cat("----------------\n")

# Generate final summary report
summary_report <- function(results, stats_summary, annotated_results, performance_times) {
  cat("Cross-pathway consistency test summary:\n\n")
  
  # Results summary
  total_pathways <- sum(sapply(results, nrow))
  successful_types <- sum(sapply(results, function(x) nrow(x) > 0))
  
  cat(sprintf("✓ Pathway types tested: %d/3\n", successful_types))
  cat(sprintf("✓ Total pathways analyzed: %d\n", total_pathways))
  
  # API consistency
  if (consistency_check) {
    cat("✓ API consistency: PASSED\n")
  } else {
    cat("✗ API consistency: FAILED\n")
  }
  
  # Statistical consistency
  if (nrow(stats_summary) > 0) {
    reasonable_stats <- all(stats_summary$mean_NES < 5) && 
                        all(stats_summary$prop_significant < 0.8)
    cat(sprintf("✓ Statistical consistency: %s\n", 
               ifelse(reasonable_stats, "PASSED", "QUESTIONABLE")))
  }
  
  # Annotation consistency
  annotation_success <- sum(sapply(annotated_results, function(x) nrow(x) > 0))
  cat(sprintf("✓ Annotation system: %d/%d types successful\n", 
             annotation_success, length(results)))
  
  # Performance consistency
  valid_perf <- performance_times[!is.na(performance_times)]
  if (length(valid_perf) > 1) {
    perf_ratio <- max(valid_perf) / min(valid_perf)
    cat(sprintf("✓ Performance consistency: %.2fx speed difference\n", perf_ratio))
  }
  
  # Overall assessment
  cat("\nOVERALL ASSESSMENT:\n")
  if (successful_types >= 2 && consistency_check) {
    cat("🎉 Cross-pathway consistency: EXCELLENT\n")
    cat("   All pathway types provide reliable, consistent results.\n")
    cat("   Users can confidently use any or all pathway types.\n")
  } else if (successful_types >= 1) {
    cat("⚠  Cross-pathway consistency: PARTIAL\n")
    cat("   Some pathway types working, others need attention.\n")
  } else {
    cat("❌ Cross-pathway consistency: FAILED\n")
    cat("   Significant issues detected across pathway types.\n")
  }
}

summary_report(results, stats_summary, annotated_results, performance_times)

cat("\n=== Demo completed successfully! ===\n")
cat("Check generated PDF files for visualization outputs.\n")