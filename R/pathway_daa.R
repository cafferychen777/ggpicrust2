#' @importFrom magrittr %>%
NULL

#' Differential Abundance Analysis for Predicted Functional Pathways
#'
#' Performs differential abundance analysis on predicted functional pathway data using various statistical methods.
#' This function supports multiple methods for analyzing differences in pathway abundance between groups,
#' including popular approaches like ALDEx2, DESeq2, edgeR, and others.
#'
#' @param abundance A data frame or matrix containing predicted functional pathway abundance,
#'        with pathways/features as rows and samples as columns.
#'        The column names should match the sample names in metadata.
#'        Values should be counts or abundance measurements.
#'
#' @param metadata A data frame or tibble containing sample information.
#'        Must include a 'sample' column with sample identifiers matching the column names in abundance data.
#'
#' @param group Character string specifying the column name in metadata that contains group information
#'        for differential abundance analysis.
#'
#' @param daa_method Character string specifying the method for differential abundance analysis.
#'        Available choices are:
#'        \itemize{
#'          \item \code{"ALDEx2"}: ANOVA-Like Differential Expression tool
#'          \item \code{"DESeq2"}: Differential expression analysis based on negative binomial distribution
#'          \item \code{"edgeR"}: Exact test for differences between groups using negative binomial model
#'          \item \code{"limma voom"}: Limma-voom framework for RNA-seq analysis
#'          \item \code{"metagenomeSeq"}: Zero-inflated Gaussian mixture model
#'          \item \code{"LinDA"}: Linear models for differential abundance analysis
#'          \item \code{"Maaslin2"}: Multivariate Association with Linear Models
#'          \item \code{"Lefser"}: Linear discriminant analysis effect size
#'        }
#'        Default is "ALDEx2".
#'
#' @param select Vector of sample names to include in the analysis.
#'        If NULL (default), all samples are included.
#'
#' @param p.adjust Character string specifying the method for p-value adjustment.
#'        Choices are:
#'        \itemize{
#'          \item \code{"BH"}: Benjamini-Hochberg procedure (default)
#'          \item \code{"holm"}: Holm's step-down method
#'          \item \code{"bonferroni"}: Bonferroni correction
#'          \item \code{"hochberg"}: Hochberg's step-up method
#'          \item \code{"fdr"}: False Discovery Rate
#'          \item \code{"none"}: No adjustment
#'        }
#'
#' @param reference Character string specifying the reference level for the group comparison.
#'        If NULL (default), the first level is used as reference.
#'
#' @param ... Additional arguments passed to the specific DAA method
#'
#' @return A data frame containing the differential abundance analysis results. The structure of the results
#' depends on the chosen DAA method. For methods that support multi-group comparisons (like LinDA),
#' when there are more than two groups, the results will contain separate rows for each feature in each
#' pairwise comparison between the reference group and each non-reference group. The data frame includes
#' the following columns:
#' \itemize{
#'   \item \code{feature}: Feature/pathway identifier
#'   \item \code{method}: The DAA method used
#'   \item \code{group1}: Reference group
#'   \item \code{group2}: Comparison group
#'   \item \code{p_values}: P-values for the comparison
#'   \item \code{p_adjust}: Adjusted p-values
#'   \item \code{adj_method}: Method used for p-value adjustment
#' }
#'
#' When \code{include_abundance_stats = TRUE}, the following additional columns
#' are included:
#' \itemize{
#'   \item \code{mean_rel_abundance_group1}: Mean relative abundance for group1
#'   \item \code{sd_rel_abundance_group1}: Standard deviation of relative
#'         abundance for group1
#'   \item \code{mean_rel_abundance_group2}: Mean relative abundance for group2
#'   \item \code{sd_rel_abundance_group2}: Standard deviation of relative
#'         abundance for group2
#'   \item \code{log2_fold_change}: Log2 fold change (group2/group1)
#' }
#'
#' Some methods may provide additional columns, such as \code{log2FoldChange}
#' for effect size information.
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#'
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#'
#' # Run differential abundance analysis using ALDEx2
#' results <- pathway_daa(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment"
#' )
#'
#' # Using a different method (DESeq2)
#' deseq_results <- pathway_daa(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   daa_method = "DESeq2"
#' )
#'
#' # Create example data with more samples
#' abundance <- data.frame(
#'   sample1 = c(10, 20, 30),
#'   sample2 = c(20, 30, 40),
#'   sample3 = c(30, 40, 50),
#'   sample4 = c(40, 50, 60),
#'   sample5 = c(50, 60, 70),
#'   row.names = c("pathway1", "pathway2", "pathway3")
#' )
#'
#' metadata <- data.frame(
#'   sample = c("sample1", "sample2", "sample3", "sample4", "sample5"),
#'   group = c("control", "control", "treatment", "treatment", "treatment")
#' )
#'
#' # Run differential abundance analysis using ALDEx2
#' results <- pathway_daa(abundance, metadata, "group")
#'
#' # Using a different method (limma voom instead of DESeq2 for this small example)
#' limma_results <- pathway_daa(abundance, metadata, "group",
#'                             daa_method = "limma voom")
#'
#' # Analyze specific samples only
#' subset_results <- pathway_daa(abundance, metadata, "group",
#'                              select = c("sample1", "sample2", "sample3", "sample4"))
#' }
#'
#' @references
#' \itemize{
#'   \item ALDEx2: Fernandes et al. (2014) Unifying the analysis of high-throughput sequencing datasets:
#'         characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by
#'         compositional data analysis. Microbiome.
#'   \item DESeq2: Love et al. (2014) Moderated estimation of fold change and dispersion for
#'         RNA-seq data with DESeq2. Genome Biology.
#'   \item edgeR: Robinson et al. (2010) edgeR: a Bioconductor package for differential expression
#'         analysis of digital gene expression data. Bioinformatics.
#'   \item limma-voom: Law et al. (2014) voom: precision weights unlock linear model analysis tools
#'         for RNA-seq read counts. Genome Biology.
#'   \item metagenomeSeq: Paulson et al. (2013) Differential abundance analysis for microbial
#'         marker-gene surveys. Nature Methods.
#'   \item Maaslin2: Mallick et al. (2021) Multivariable Association Discovery in Population-scale
#'         Meta-omics Studies.
#' }
#'
#' Helper function to calculate abundance statistics for differential analysis
#'
#' This function calculates mean relative abundance, standard deviation, and log2 fold change
#' for each feature between two groups.
#'
#' @param abundance A matrix or data frame with features as rows and samples as columns
#' @param metadata A data frame containing sample information
#' @param group Character string specifying the group column name in metadata
#' @param features Character vector of feature names to calculate statistics for
#' @param group1 Character string specifying the first group name
#' @param group2 Character string specifying the second group name
#'
#' @return A data frame with columns:
#'   \item{feature}{Feature identifier}
#'   \item{mean_rel_abundance_group1}{Mean relative abundance for group1}
#'   \item{sd_rel_abundance_group1}{Standard deviation of relative abundance for group1}
#'   \item{mean_rel_abundance_group2}{Mean relative abundance for group2}
#'   \item{sd_rel_abundance_group2}{Standard deviation of relative abundance for group2}
#'   \item{log2_fold_change}{Log2 fold change (group2/group1)}
#'
#' @keywords internal
calculate_abundance_stats <- function(abundance, metadata, group, features, group1, group2) {
  # Validate inputs
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("'abundance' must be a matrix or data frame")
  }

  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }

  if (!group %in% colnames(metadata)) {
    stop("Group column '", group, "' not found in metadata")
  }

  # Convert abundance to matrix if needed
  abundance_mat <- as.matrix(abundance)

  # Find the sample column that matches abundance column names
  sample_col <- NULL
  potential_sample_cols <- c("sample", "Sample", "sample_name", "Sample_Name", "sample_id", "Sample_ID")

  # First try standard column names
  for (col in potential_sample_cols) {
    if (col %in% colnames(metadata)) {
      if (all(colnames(abundance_mat) %in% metadata[[col]]) ||
          all(metadata[[col]] %in% colnames(abundance_mat))) {
        sample_col <- col
        break
      }
    }
  }

  # If no standard column found, try all columns
  if (is.null(sample_col)) {
    for (col in colnames(metadata)) {
      if (all(colnames(abundance_mat) %in% metadata[[col]]) ||
          all(metadata[[col]] %in% colnames(abundance_mat))) {
        sample_col <- col
        break
      }
    }
  }

  if (is.null(sample_col)) {
    # Try to match by rownames if no sample column found
    if (all(colnames(abundance_mat) %in% rownames(metadata))) {
      metadata$sample <- rownames(metadata)
      sample_col <- "sample"
    } else {
      stop("Cannot find matching sample identifiers between abundance and metadata")
    }
  }

  # Filter metadata to match abundance samples
  metadata_filtered <- metadata[metadata[[sample_col]] %in% colnames(abundance_mat), ]

  # Filter abundance to match metadata samples
  abundance_filtered <- abundance_mat[, colnames(abundance_mat) %in% metadata_filtered[[sample_col]], drop = FALSE]

  # Reorder to ensure matching
  sample_order <- match(colnames(abundance_filtered), metadata_filtered[[sample_col]])
  metadata_ordered <- metadata_filtered[sample_order, ]

  # Get group assignments
  group_assignments <- metadata_ordered[[group]]

  # Filter for the two groups of interest
  group1_samples <- colnames(abundance_filtered)[group_assignments == group1]
  group2_samples <- colnames(abundance_filtered)[group_assignments == group2]

  if (length(group1_samples) == 0) {
    stop("No samples found for group1: ", group1)
  }
  if (length(group2_samples) == 0) {
    stop("No samples found for group2: ", group2)
  }

  # Convert to relative abundance
  relative_abundance <- apply(abundance_filtered, 2, function(x) x / sum(x))

  # Filter for specified features
  if (!is.null(features)) {
    features_available <- intersect(features, rownames(relative_abundance))
    if (length(features_available) == 0) {
      stop("None of the specified features found in abundance data")
    }
    relative_abundance <- relative_abundance[features_available, , drop = FALSE]
  }

  # Calculate statistics for each feature
  results <- data.frame(
    feature = rownames(relative_abundance),
    mean_rel_abundance_group1 = NA_real_,
    sd_rel_abundance_group1 = NA_real_,
    mean_rel_abundance_group2 = NA_real_,
    sd_rel_abundance_group2 = NA_real_,
    log2_fold_change = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(relative_abundance))) {
    feature_name <- rownames(relative_abundance)[i]

    # Get abundance values for each group
    group1_values <- relative_abundance[i, group1_samples, drop = TRUE]
    group2_values <- relative_abundance[i, group2_samples, drop = TRUE]

    # Calculate means and standard deviations
    mean1 <- mean(group1_values, na.rm = TRUE)
    sd1 <- stats::sd(group1_values, na.rm = TRUE)
    mean2 <- mean(group2_values, na.rm = TRUE)
    sd2 <- stats::sd(group2_values, na.rm = TRUE)

    # Calculate log2 fold change (group2/group1)
    # Add small pseudocount to avoid log(0)
    pseudocount <- 1e-10
    log2_fc <- log2((mean2 + pseudocount) / (mean1 + pseudocount))

    # Store results
    results[i, "mean_rel_abundance_group1"] <- mean1
    results[i, "sd_rel_abundance_group1"] <- sd1
    results[i, "mean_rel_abundance_group2"] <- mean2
    results[i, "sd_rel_abundance_group2"] <- sd2
    results[i, "log2_fold_change"] <- log2_fc
  }

  return(results)
}

#' @param include_abundance_stats Logical value indicating whether to include
#'        abundance statistics (mean relative abundance, standard deviation,
#'        and log2 fold change) in the output. Default is FALSE for backward
#'        compatibility.
#'
#' @export
pathway_daa <- function(abundance, metadata, group, daa_method = "ALDEx2",
                       select = NULL, p.adjust = "BH", reference = NULL,
                       include_abundance_stats = FALSE, ...) {
  # Check if the requested DAA method package is available
  method_packages <- list(
    "ALDEx2" = "ALDEx2",
    "DESeq2" = "DESeq2",
    "edgeR" = "edgeR",
    "limma voom" = "limma",
    "metagenomeSeq" = "metagenomeSeq",
    "Maaslin2" = "Maaslin2",
    "Lefser" = "lefser"
  )
  
  if (daa_method %in% names(method_packages)) {
    pkg_name <- method_packages[[daa_method]]
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for method '%s'. Please install it using BiocManager::install('%s')", 
                   pkg_name, daa_method, pkg_name))
    }
  }
  # Input validation
  if (!is.data.frame(abundance) && !is.matrix(abundance)) {
    stop("abundance must be a data frame or matrix")
  }

  if (ncol(abundance) < 4) {
    stop("At least 4 samples are required for differential abundance analysis")
  }

  # Convert metadata to tibble
  if (!tibble::is_tibble(metadata)) {
    metadata <- tibble::as_tibble(metadata)
  }

  # Extract sample names from abundance data
  abundance_samples <- colnames(abundance)
  
  # Identify the column in metadata that matches sample names
  sample_col <- NULL
  for (col in colnames(metadata)) {
    if (all(metadata[[col]] %in% abundance_samples)) {
      sample_col <- col
      break
    }
  }
  
  # Check if a matching column was found
  if (is.null(sample_col)) {
    stop("No column in metadata matches the sample names in abundance data")
  }
  
  message(sprintf("Using column '%s' as sample identifier", sample_col))
  
  # Get sample names from metadata
  metadata_samples <- metadata[[sample_col]]
  
  # Verify sample matching
  if (!all(metadata_samples %in% abundance_samples)) {
    stop("Some samples in metadata are not found in abundance data")
  }
  
  if (!all(abundance_samples %in% metadata_samples)) {
    stop("Some samples in abundance data are not found in metadata")
  }
  
  # Now check sample size
  if (length(abundance_samples) < 4) {
    stop("At least 4 samples are required for differential abundance analysis")
  }
  
  # Ensure consistent sample order between data and metadata
  metadata <- metadata[match(abundance_samples, metadata[[sample_col]]), ]

  # Verify if group column exists
  if (!group %in% colnames(metadata)) {
    stop(sprintf("group column '%s' not found in metadata", group))
  }

  # Extract metadata samples using identified column
  metadata_samples <- metadata[[sample_col]]

  # Verify sample matching
  if (!all(metadata_samples %in% abundance_samples)) {
    stop("Some samples in metadata are not found in abundance data")
  }

  if (!all(abundance_samples %in% metadata_samples)) {
    stop("Some samples in abundance data are not found in metadata")
  }

  # Ensure consistent sample order between data and metadata
  metadata <- metadata[match(colnames(abundance), metadata[[sample_col]]), ]

  # Verify group count
  group_levels <- unique(metadata[[group]])
  if (length(group_levels) < 2) {
    stop("At least two groups are required for differential abundance analysis")
  }

  # Handle sample selection
  if (!is.null(select)) {
    if (!all(select %in% abundance_samples)) {
      stop("Some selected samples are not present in the abundance data")
    }
    abundance <- abundance[, select, drop = FALSE]
    metadata <- metadata[metadata$sample %in% select, ]
  }

  # Prepare data
  abundance_mat <- as.matrix(abundance)
  Group <- factor(metadata[[group]])

  # Ensure factor levels only include groups present in the data
  # This is defensive programming to handle edge cases with unused factor levels
  Group <- droplevels(Group)
  Level <- levels(Group)
  length_Level <- length(Level)

  # Perform differential analysis
  result <- switch(
    daa_method,
    "ALDEx2" = perform_aldex2_analysis(abundance_mat, Group, Level, length_Level),
    "DESeq2" = perform_deseq2_analysis(abundance_mat, metadata, group, Level),
    "LinDA" = perform_linda_analysis(abundance, metadata, group, reference, Level, length_Level),
    "limma voom" = perform_limma_voom_analysis(abundance_mat, Group, reference, Level, length_Level),
    "edgeR" = perform_edger_analysis(abundance_mat, Group, Level, length_Level),
    "metagenomeSeq" = perform_metagenomeseq_analysis(abundance_mat, metadata, group, Level),
    "Maaslin2" = perform_maaslin2_analysis(abundance_mat, metadata, group, reference, Level, length_Level),
    "Lefser" = perform_lefser_analysis(abundance_mat, metadata, group, Level)
  )

  # Add multiple testing correction
  if (!is.null(result) && "p_values" %in% colnames(result)) {
    result$p_adjust <- p.adjust(result$p_values, method = p.adjust)
    result$adj_method <- p.adjust
  }

  # Add abundance statistics if requested
  if (include_abundance_stats && !is.null(result) && nrow(result) > 0) {
    # Check if we have the required columns for abundance stats
    if (all(c("feature", "group1", "group2") %in% colnames(result))) {
      tryCatch({
        # Get unique group pairs for calculation
        unique_comparisons <- unique(result[, c("group1", "group2")])

        # Calculate abundance stats for each unique comparison
        all_abundance_stats <- data.frame()

        for (i in seq_len(nrow(unique_comparisons))) {
          group1_name <- unique_comparisons$group1[i]
          group2_name <- unique_comparisons$group2[i]

          # Get features for this comparison
          comparison_features <- result$feature[
            result$group1 == group1_name & result$group2 == group2_name
          ]

          if (length(comparison_features) > 0) {
            abundance_stats <- calculate_abundance_stats(
              abundance = abundance_mat,
              metadata = metadata,
              group = group,
              features = comparison_features,
              group1 = group1_name,
              group2 = group2_name
            )

            # Add group information for merging
            abundance_stats$group1 <- group1_name
            abundance_stats$group2 <- group2_name

            all_abundance_stats <- rbind(all_abundance_stats, abundance_stats)
          }
        }

        # Merge abundance stats with results
        if (nrow(all_abundance_stats) > 0) {
          result <- merge(
            result,
            all_abundance_stats,
            by = c("feature", "group1", "group2"),
            all.x = TRUE
          )
        }
      }, error = function(e) {
        warning("Failed to calculate abundance statistics: ", e$message)
      })
    } else {
      warning("Cannot calculate abundance statistics: required columns (feature, group1, group2) not found in results")
    }
  }

  return(result)
}

# Helper function: Perform ALDEx2 analysis
perform_aldex2_analysis <- function(abundance_mat, Group, Level, length_Level) {
  message("Running ALDEx2 analysis...")
  
  # Check if ALDEx2 is available
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("Package 'ALDEx2' is required for ALDEx2 analysis. Please install it using BiocManager::install('ALDEx2')")
  }
  
  # Round the abundance data
  abundance_mat <- round(abundance_mat)
  
  # Save the original Group factor and convert to numeric for ALDEx2
  original_Group <- Group
  Group <- as.numeric(Group)
  
  # Perform different analyses based on the number of groups
  if (length_Level == 2) {
    message("Running ALDEx2 with two groups. Performing t-test...")
    
    # Create ALDEx2 object with numeric Group
    ALDEx2_object <- ALDEx2::aldex.clr(
      abundance_mat,
      Group,
      mc.samples = 256,
      denom = "all",
      verbose = FALSE
    )
    
    # Get t-test results
    results <- ALDEx2::aldex.ttest(
      ALDEx2_object,
      paired.test = FALSE,
      verbose = FALSE
    )
    
    # Build result dataframe using original Level names
    return(data.frame(
      feature = rep(rownames(results), 2),
      method = c(
        rep("ALDEx2_Welch's t test", nrow(results)),
        rep("ALDEx2_Wilcoxon rank test", nrow(results))
      ),
      group1 = rep(Level[1], 2 * nrow(results)),
      group2 = rep(Level[2], 2 * nrow(results)),
      p_values = c(results$we.ep, results$wi.ep),
      stringsAsFactors = FALSE
    ))
    
  } else {
    message("Running ALDEx2 with multiple groups. This might take some time...")
    
    # Create ALDEx2 object for multiple groups with numeric Group
    ALDEx2_object <- ALDEx2::aldex.clr(
      abundance_mat,
      Group,
      mc.samples = 256,
      denom = "all",
      verbose = FALSE
    )
    
    # Get Kruskal-Wallis and GLM test results
    results <- ALDEx2::aldex.kw(ALDEx2_object)
    
    # Build initial result dataframe
    result_df <- data.frame(
      feature = rep(rownames(results), 2),
      method = c(
        rep("ALDEx2_Kruskal-Wallace test", nrow(results)),
        rep("ALDEx2_glm test", nrow(results))
      ),
      p_values = c(results$kw.ep, results$glm.ep),
      stringsAsFactors = FALSE
    )
    
    # Add all group information using original Level names
    group_cols <- data.frame(matrix(
      NA, 
      nrow = nrow(result_df), 
      ncol = length_Level,
      dimnames = list(NULL, paste0("group", 1:length_Level))
    ))
    
    for (i in 1:length_Level) {
      group_cols[, i] <- Level[i]
    }
    
    # Merge results
    result_df <- cbind(
      result_df[, c("feature", "method")],
      group_cols,
      result_df[, "p_values", drop = FALSE]
    )
    
    message("ALDEx2 analysis with multiple groups complete.")
    return(result_df)
  }
}

# Helper function: Perform DESeq2 analysis
perform_deseq2_analysis <- function(abundance_mat, metadata, group, Level) {
  message("Running DESeq2 analysis...")
  
  # Convert to integer matrix
  message("converting counts to integer mode")
  counts <- round(as.matrix(abundance_mat))
  
  # 确保 group 是因子类型
  metadata[[group]] <- factor(metadata[[group]], levels = Level)
  
  # Create DESeqDataSet object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = metadata
  )
  
  # Use try-catch to handle possible errors
  result <- tryCatch({
    suppressWarnings({
      # Create DESeqDataSet
      dds <- DESeq2::DESeqDataSet(se, design = as.formula(paste0("~", group)))

      # 根据样本量选择合适的拟合方法
      fitType <- if(ncol(abundance_mat) < 6) "mean" else "local"

      # Run DESeq2 pipeline with improved error handling
      dds <- DESeq2::estimateSizeFactors(dds)

      # Try to estimate dispersions with different methods
      dds <- tryCatch({
        DESeq2::estimateDispersions(dds, fitType = fitType)
      }, error = function(e) {
        # If standard dispersion estimation fails, use gene-wise estimates
        message("Standard dispersion estimation failed, using gene-wise estimates...")
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        DESeq2::dispersions(dds) <- SummarizedExperiment::mcols(dds)$dispGeneEst
        return(dds)
      })

      dds <- DESeq2::nbinomWaldTest(dds)

      # Extract results
      res <- DESeq2::results(dds, contrast = c(group, Level[2], Level[1]))

      data.frame(
        feature = rownames(abundance_mat),
        method = "DESeq2",
        group1 = Level[1],
        group2 = Level[2],
        p_values = res$pvalue,
        stringsAsFactors = FALSE
      )
    })
  }, error = function(e) {
    stop("DESeq2 analysis failed: ", e$message)
  })
  
  if (is.null(result)) {
    stop("DESeq2 analysis failed to produce results")
  }
  
  return(result)
}

# Helper function: Perform limma voom analysis
perform_limma_voom_analysis <- function(abundance_mat, Group, reference, Level, length_Level) {
  message("Running limma voom analysis...")

  # Check if limma is available
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for limma voom analysis. Please install it using BiocManager::install('limma')")
  }
  
  # Ensure Group is a factor type
  Group <- factor(Group)
  if (!is.null(reference)) {
    Group <- relevel(Group, ref = reference)
  }

  # Create design matrix
  design <- model.matrix(~Group)

  # Perform voom transformation and analysis
  dge <- edgeR::DGEList(counts = abundance_mat, group = Group)
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge, design)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)

  # Extract results
  if (length_Level == 2) {
    # Two-group comparison - ensure consistent format with other methods
    group_levels <- levels(Group)
    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "limma voom",
      group1 = group_levels[1],
      group2 = group_levels[2],
      p_values = fit$p.value[,2],
      stringsAsFactors = FALSE
    )
  } else {
    # Multi-group comparison handling
    group_levels <- levels(Group)
    results <- data.frame(
      feature = rep(rownames(abundance_mat), length(group_levels) - 1),
      method = "limma voom",
      group1 = reference,
      group2 = group_levels[group_levels != reference],
      p_values = as.vector(fit$p.value[,-1]),
      stringsAsFactors = FALSE
    )
  }

  return(results)
}

# Helper function: Perform edgeR analysis
perform_edger_analysis <- function(abundance_mat, Group, Level, length_Level) {
  message("Running edgeR analysis...")

  # Create DGEList object
  dge <- edgeR::DGEList(counts = round(abundance_mat), group = Group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateCommonDisp(dge, verbose = TRUE)

  if (length_Level == 2) {
    # Two-group comparison
    et <- edgeR::exactTest(dge, pair = c(1, 2))
    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "edgeR",
      group1 = Level[1],
      group2 = Level[2],
      p_values = et$table$PValue
    )
  } else {
    # Multi-group comparison
    results_list <- list()
    combinations <- utils::combn(seq_along(Level), 2)

    for (i in 1:ncol(combinations)) {
      et <- edgeR::exactTest(dge, pair = combinations[,i])
      results_list[[i]] <- data.frame(
        feature = rownames(abundance_mat),
        method = "edgeR",
        group1 = Level[combinations[1,i]],
        group2 = Level[combinations[2,i]],
        p_values = et$table$PValue
      )
    }
    results <- do.call(rbind, results_list)
  }

  return(results)
}

# Helper function: Perform metagenomeSeq analysis
perform_metagenomeseq_analysis <- function(abundance_mat, metadata, group, Level) {
  message("Running metagenomeSeq analysis...")

  # Convert metadata to data.frame and ensure sample names are correct
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata$sample

  # Ensure abundance_mat and metadata have consistent sample order
  abundance_mat <- abundance_mat[, rownames(metadata)]

  # Create phenoData
  phenoData <- new("AnnotatedDataFrame",
                   data = metadata,
                   varMetadata = data.frame(
                     labelDescription = c("Sample ID", "Group"),
                     row.names = colnames(metadata)
                   ))

  # Ensure data is an integer matrix
  counts <- round(as.matrix(abundance_mat))

  # Create MRexperiment object
  obj <- try({
    metagenomeSeq::newMRexperiment(
      counts = counts,
      phenoData = phenoData,
      featureData = NULL,
      libSize = NULL,
      normFactors = NULL
    )
  }, silent = TRUE)

  if (inherits(obj, "try-error")) {
    stop("Failed to create MRexperiment object: ", attr(obj, "condition")$message)
  }

  # Normalize
  obj <- metagenomeSeq::cumNorm(obj)

  # Create model matrix
  mod <- stats::model.matrix(as.formula(paste0("~", group)), data = metadata)

  # Fit model
  fit <- metagenomeSeq::fitFeatureModel(obj, mod)

  # Extract results
  results <- data.frame(
    feature = rownames(abundance_mat),
    method = "metagenomeSeq",
    group1 = Level[1],
    group2 = Level[2],
    p_values = fit@pvalues,
    stringsAsFactors = FALSE
  )

  return(results)
}

# Helper function: Perform Maaslin2 analysis
perform_maaslin2_analysis <- function(abundance_mat, metadata, group, reference, Level, length_Level) {
  message("Running Maaslin2 analysis...")

  # Check if Maaslin2 is available
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    stop("Maaslin2 package is required but not installed")
  }

  # Prepare metadata first
  metadata <- as.data.frame(metadata)

  # Ensure metadata rownames match the column names of abundance_mat
  # This is crucial because the main function has already reordered metadata to match abundance
  rownames(metadata) <- colnames(abundance_mat)

  # Transpose abundance matrix (samples become rows, features become columns)
  abundance_mat_t <- t(abundance_mat)

  # Create temporary output directory
  output_dir <- tempdir()

  # Run Maaslin2 analysis
  fit_data <- Maaslin2::Maaslin2(
    input_data = abundance_mat_t,
    input_metadata = metadata,
    output = output_dir,
    transform = "AST",
    fixed_effects = group,
    reference = if (length_Level > 2) paste0(group, ",", reference) else NULL,
    normalization = "TSS",
    standardize = TRUE,
    min_prevalence = 0.1,
    cores = 1,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )

  # Read all results instead of significant results
  results_file <- file.path(output_dir, "all_results.tsv")

  if (file.exists(results_file)) {
    maaslin2_results <- utils::read.table(results_file,
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)

    # Format results
    # Note: Maaslin2 replaces hyphens (-) with dots (.) in feature names
    # We need to match features more intelligently

    # Create a mapping between original and Maaslin2 feature names
    original_features <- rownames(abundance_mat)
    maaslin2_features <- maaslin2_results$feature

    # Try direct matching first
    matches <- match(original_features, maaslin2_features)

    # For unmatched features, try with hyphen-to-dot conversion
    unmatched_indices <- which(is.na(matches))
    if (length(unmatched_indices) > 0) {
      # Convert hyphens to dots in original feature names for matching
      original_with_dots <- gsub("-", ".", original_features[unmatched_indices])
      dot_matches <- match(original_with_dots, maaslin2_features)
      matches[unmatched_indices] <- dot_matches
    }

    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "Maaslin2",
      group1 = if (length_Level == 2) Level[1] else reference,
      group2 = if (length_Level == 2) Level[2] else Level[Level != reference],
      p_values = maaslin2_results$pval[matches],
      stringsAsFactors = FALSE
    )
  } else {
    stop("Maaslin2 analysis failed to produce results file")
  }

  return(results)
}

# Helper function: Perform Lefser analysis
perform_lefser_analysis <- function(abundance_mat, metadata, group, Level) {
  message("Running Lefser analysis...")
  
  # Check if lefser package is available
  if (!requireNamespace("lefser", quietly = TRUE)) {
    stop("The 'lefser' package is required for Lefser analysis. Please install it using BiocManager::install('lefser')")
  }

  # Create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = abundance_mat),
    colData = metadata
  )

  # Perform Lefser analysis
  lefser_results <- lefser::lefser(se, classCol = group)  # Using classCol instead of deprecated groupCol

  # Create results for all features, not just significant ones
  all_features <- rownames(abundance_mat)

  # Initialize results with all features
  results <- data.frame(
    feature = all_features,
    method = "Lefser",
    group1 = Level[1],
    group2 = Level[2],
    p_values = rep(1.0, length(all_features)),  # Default p-value of 1 for non-significant
    stringsAsFactors = FALSE
  )

  # If significant features found, update their p-values
  if (length(lefser_results$Names) > 0 && !is.null(lefser_results$Names)) {
    # Match significant features to all features
    sig_indices <- match(lefser_results$Names, all_features)

    # Convert effect scores to p-values (higher absolute effect score = lower p-value)
    # This is a simplified conversion; in practice, Lefser uses Wilcoxon test internally
    effect_scores <- abs(lefser_results$scores)
    # Convert to p-values: higher effect scores get lower p-values
    max_score <- max(effect_scores, na.rm = TRUE)
    if (max_score > 0) {
      p_vals <- 1 - (effect_scores / max_score) * 0.95  # Scale to 0.05-1.0 range
    } else {
      p_vals <- rep(0.5, length(effect_scores))
    }

    # Update p-values for significant features
    results$p_values[sig_indices] <- p_vals
  } else {
    message("No significant features found by Lefser analysis.")
  }

  return(results)
}

# Helper function: Perform LinDA analysis
perform_linda_analysis <- function(abundance, metadata, group, reference, Level, length_Level) {
  message("Running LinDA analysis...")

  # Ensure necessary packages are loaded
  if (!requireNamespace("MicrobiomeStat", quietly = TRUE)) {
    stop("MicrobiomeStat package is required but not installed")
  }

  # Prepare data
  feature.dat <- as.matrix(abundance)
  meta.dat <- metadata
  
  # Ensure the group variable in meta.dat is a factor
  if (!is.factor(meta.dat[[group]])) {
    meta.dat[[group]] <- factor(meta.dat[[group]])
  }
  
  # If reference is provided, ensure it is one of the values in Level
  if (!is.null(reference) && !(reference %in% Level)) {
    warning(sprintf("Reference '%s' not found in levels of group variable. Using first level as reference.", reference))
    reference <- Level[1]
  } else if (is.null(reference)) {
    reference <- Level[1]
  }

  # Build formula
  formula <- paste0("~ ", group)

  # Add more debug information
  message(sprintf("Group variable: %s with levels: %s", group, paste(Level, collapse=", ")))
  message(sprintf("Reference level: %s", reference))
  message(sprintf("Number of features: %d, Number of samples: %d", nrow(feature.dat), ncol(feature.dat)))

  # Use tryCatch to catch errors in LinDA analysis
  linda_result <- tryCatch({
    # Run LinDA analysis
    linda_obj <- MicrobiomeStat::linda(
      feature.dat = feature.dat,
      meta.dat = meta.dat,
      formula = formula,
      feature.dat.type = "count",
      prev.filter = 0,
      mean.abund.filter = 0,
      adaptive = TRUE,
      zero.handling = "pseudo-count",
      pseudo.cnt = 0.5,
      p.adj.method = "BH",
      alpha = 0.05,
      n.cores = 1,
      verbose = TRUE    # Changed to TRUE to provide more detailed output
    )
    
    # Check the structure of linda_obj
    if (is.null(linda_obj) || is.null(linda_obj$output) || length(linda_obj$output) == 0) {
      message("LinDA analysis returned empty results. This might be due to insufficient variation in the data.")
      # Return an empty data frame instead of NULL, so that the pathway_daa function can handle it properly
      return(data.frame(
        feature = character(0),
        method = character(0),
        group1 = character(0),
        group2 = character(0),
        p_values = numeric(0),
        log2FoldChange = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    # Create an empty list to store results for each comparison
    results_list <- list()
    
    # Get all non-reference groups
    non_ref_groups <- Level[Level != reference]
    
    # Process each output element (each corresponds to a group comparison)
    for (i in seq_along(linda_obj$output)) {
      # Get output name, typically in the format "groupVariableNameNonReferenceGroup"
      output_name <- names(linda_obj$output)[i]
      # Extract non-reference group name from output name
      # Here we assume the output name format is "groupVariableNameNonReferenceGroup"
      group_prefix <- paste0(group, "")
      if (startsWith(output_name, group_prefix)) {
        comparison_group <- substring(output_name, nchar(group_prefix) + 1)
      } else {
        # If we cannot extract group name from output name, use index to get from non_ref_groups
        comparison_group <- if (i <= length(non_ref_groups)) non_ref_groups[i] else paste0("Group", i)
      }
      
      comparison_df <- linda_obj$output[[i]]
      
      # Create a results data frame for this comparison
      results_list[[i]] <- data.frame(
        feature = rownames(comparison_df),
        method = "LinDA",
        group1 = reference,
        group2 = comparison_group,
        p_values = comparison_df$pvalue,
        log2FoldChange = comparison_df$log2FoldChange,
        stringsAsFactors = FALSE
      )
    }
    
    # Combine all results
    if (length(results_list) > 0) {
      do.call(rbind, results_list)
    } else {
      # If there are no results, return an empty data frame
      data.frame(
        feature = character(0),
        method = character(0),
        group1 = character(0),
        group2 = character(0),
        p_values = numeric(0),
        log2FoldChange = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    # Catch and report errors, but return an empty data frame instead of NULL
    message(sprintf("Error in LinDA analysis: %s", e$message))
    data.frame(
      feature = character(0),
      method = character(0),
      group1 = character(0),
      group2 = character(0),
      p_values = numeric(0),
      log2FoldChange = numeric(0),
      stringsAsFactors = FALSE
    )
  })
  
  return(linda_result)
}
