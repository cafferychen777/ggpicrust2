#' @importFrom magrittr %>%
NULL

#' Differential Abundance Analysis for Predicted Functional Pathways
#'
#' @description
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
#' @return A data frame containing the differential abundance analysis results
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#'
#' # Run differential abundance analysis using ALDEx2
#' results <- pathway_daa(
#'   abundance = ko_abundance %>% column_to_rownames("#NAME"),
#'   metadata = metadata,
#'   group = "Environment"
#' )
#'
#' # Using a different method (DESeq2)
#' deseq_results <- pathway_daa(
#'   abundance = ko_abundance %>% column_to_rownames("#NAME"),
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
#' # Using a different method (DESeq2)
#' deseq_results <- pathway_daa(abundance, metadata, "group",
#'                             daa_method = "DESeq2")
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
#' @export
pathway_daa <- function(abundance, metadata, group, daa_method = "ALDEx2",
                       select = NULL, p.adjust = "BH", reference = NULL) {
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

  return(result)
}

# Helper function: Perform ALDEx2 analysis
perform_aldex2_analysis <- function(abundance_mat, Group, Level, length_Level) {
  message("Running ALDEx2 analysis...")
  
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
      
      # Run DESeq2 pipeline
      dds <- DESeq2::estimateSizeFactors(dds)
      dds <- DESeq2::estimateDispersions(dds, fitType = fitType)
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
    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "limma voom",
      p_values = fit$p.value[,2]
    )
  } else {
    # Multi-group comparison handling
    group_levels <- levels(Group)
    results <- data.frame(
      feature = rep(rownames(abundance_mat), length(group_levels) - 1),
      method = "limma voom",
      group1 = reference,
      group2 = group_levels[group_levels != reference],
      p_values = as.vector(fit$p.value[,-1])
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

  # Ensure necessary packages are loaded
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    stop("Maaslin2 package is required but not installed")
  }

  # Transpose abundance matrix
  abundance_mat_t <- t(abundance_mat)

  # Prepare metadata
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata$sample

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
    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "Maaslin2",
      group1 = if (length_Level == 2) Level[1] else reference,
      group2 = if (length_Level == 2) Level[2] else Level[Level != reference],
      p_values = maaslin2_results$pval[match(rownames(abundance_mat), maaslin2_results$feature)],
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

  # Create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = abundance_mat),
    colData = metadata
  )

  # Perform Lefser analysis
  lefser_results <- lefser::lefser(se, groupCol = group)

  # Extract results
  results <- data.frame(
    feature = lefser_results$Names,
    method = "Lefser",
    group1 = Level[1],
    group2 = Level[2],
    effect_scores = lefser_results$scores
  )

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

  # Build formula
  formula <- paste0("~ ", group)

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
    verbose = FALSE
  )

  # Extract results
  if (length(linda_obj$output) > 0) {
    # Get results of the first variable (usually group comparison)
    first_comparison <- linda_obj$output[[1]]

    results <- data.frame(
      feature = rownames(first_comparison),
      method = "LinDA",
      group1 = if (length_Level == 2) Level[1] else reference,
      group2 = if (length_Level == 2) Level[2] else Level[Level != reference],
      p_values = first_comparison$pvalue,
      stringsAsFactors = FALSE
    )
  } else {
    stop("LinDA analysis failed to produce results")
  }

  return(results)
}
