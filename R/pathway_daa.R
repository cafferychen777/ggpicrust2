#' Predictional functional patwhay differential abundance (DA)
#'
#' @param abundance a data frame containing predicted functional pathway abundance, with pathways/features as rows and samples as columns. The column names of abundance should match the sample names in metadata. Pathway abundance values should be counts
#' @param metadata a tibble containing samples information
#' @param group a character specifying the group name for differential abundance analysis
#' @param daa_method a character specifying the method for differential abundance analysis, choices are:
#' - "ALDEx2": ANOVA-Like Differential Expression tool for high throughput sequencing data
#' - "DESeq2": Differential expression analysis based on the negative binomial distribution using DESeq2
#' - "edgeR": Exact test for differences between two groups of negative-binomially distributed counts using edgeR
#' - "limma voom": Limma-voom framework for the analysis of RNA-seq data
#' - "metagenomeSeq": Fit logistic regression models to test for differential abundance between groups using metagenomeSeq
#' - "LinDA": Linear models for differential abundance analysis of microbiome compositional data
#' - "Maaslin2": Multivariate Association with Linear Models (MaAsLin2) for differential abundance analysis
#' - "Lefser": Linear discriminant analysis (LDA) effect size algorithm for high-dimensional microbiome data
#' @default "ALDEx2"
#' @param select a vector containing sample names for analysis, if NULL all samples are included. This parameter can be used to specify which samples are included in the differential abundance analysis. Default is NULL.
#' @param p.adjust a character specifying the method for p-value adjustment, choices are:
#'- "BH": Benjamini-Hochberg correction
#'- "holm": Holm's correction
#'- "bonferroni": Bonferroni correction
#'- "hochberg": Hochberg's correction
#'- "fdr": False discovery rate correction
#'- "none": No p-value adjustment.
#' @default "BH"
#' @param reference a character specifying the reference group level, required for several differential abundance analysis methods such as LinDA, limme voom and Maaslin2. This parameter is used to specify the reference group when there are more than two groups. Default is NULL.
#'
#' @return  a data frame containing the differential abundance analysis results.
#' @value
#' A data frame containing the differential abundance analysis results. The data frame has the following columns:
#' \itemize{
#'   \item \code{feature}: The feature ID of the pathway.
#'   \item \code{statistic}: The test statistic for the differential abundance analysis.
#'   \item \code{p_value}: The raw p-value for the differential abundance analysis.
#'   \item \code{p_adjust}: The adjusted p-value for the differential abundance analysis.
#' }
#' @export
#'
#' @examples
#' \donttest{
#' library(ggpicrust2)
#' library(MicrobiomeStat)
#' library(tibble)
#' library(magrittr)
#' abundance <- data.frame(sample1 = c(10, 20, 30),
#' sample2 = c(20, 30, 40),
#' sample3 = c(30, 40, 50),
#' row.names = c("pathway1", "pathway2", "pathway3"))
#'
#' metadata <- tibble::tibble(sample = paste0("sample", 1:3),
#' group = c("control", "control", "treatment"))
#'
#' #Run pathway_daa function
#' result <- pathway_daa(abundance = abundance, metadata = metadata, group = "group",
#' daa_method = "LinDA")
#'
#' data(metacyc_abundance)
#' data(metadata)
#' daa_results_df <- pathway_daa(metacyc_abundance %>%
#' column_to_rownames("pathway"), metadata, "Environment", daa_method = "Lefser")
#' }
pathway_daa <-
  function(abundance,
           metadata,
           group,
           daa_method = "ALDEx2",
           select = NULL,
           p.adjust = "BH",
           reference = NULL) {

    # Define valid differential abundance analysis methods
    valid_methods <- c("ALDEx2", "DESeq2", "edgeR", "limma voom", "metagenomeSeq", "LinDA", "Maaslin2", "Lefser")

    # Check if the provided differential abundance analysis method is valid
    if (!daa_method %in% valid_methods) {
      stop(paste("Invalid differential abundance analysis method. Please choose from:",
                 paste(valid_methods, collapse = ", ")))
    }

    # Convert metadata to a tibble if it is not already
    if (!tibble::is_tibble(metadata)) {
      message("Converting metadata to tibble...")
      metadata <- tibble::as_tibble(metadata)
    }

    # Extract the sample names
    sample_names <- colnames(abundance)
    message("Sample names extracted.")

    # Identify the columns in metadata that match the sample names
    message("Identifying matching columns in metadata...")
    matches <- base::lapply(metadata, function(x) {
      intersect(sample_names, x)
    })
    matching_columns <- names(metadata)[sapply(matches, function(x) {
      length(x) == length(sample_names)
    })][1]

    # Highlight the importance of matching_columns
    if (!is.null(matching_columns)) {
      message(paste("Matching columns identified:", matching_columns,
                    ". This is important for ensuring data consistency."))
    } else {
      stop("No matching columns found. Please check your input data.")
    }

    # If select is null, use all columns, otherwise filter the abundance and metadata
    switch(is.null(select),
           "TRUE" = {
             message("Using all columns in abundance.")
           },
           "FALSE" = {
             message("Filtering the abundance and metadata...")
             abundance <- abundance[, colnames(abundance) %in% select]
             metadata <- metadata[as.matrix(metadata[, matching_columns]) %in% select,]
           })

    # Reextract sample names after filtering
    sample_names <- colnames(abundance)

    # Convert abundance to a matrix
    message("Converting abundance to a matrix...")
    abundance_mat <- as.matrix(abundance)

    # Reorder metadata to match the sample names in abundance
    message("Reordering metadata...")
    metadata_order <- match(sample_names, as.matrix(metadata[, matching_columns]))
    metadata <- metadata[metadata_order,]

    # Convert metadata to a matrix and data frame
    message("Converting metadata to a matrix and data frame...")
    metadata_mat <- as.matrix(metadata)
    metadata_df <- as.data.frame(metadata)

    # Extract group information and calculate the number of levels
    message("Extracting group information...")
    Group <- factor(metadata_mat[, group])
    Level <- levels(Group)
    length_Level <- length(Level)

    switch(
      daa_method,
      "ALDEx2" = {
        # Round the abundance data
        ALDEx2_abundance <- round(abundance)

        # Check if there are two levels and execute the relevant analysis
        if (length_Level == 2){
          message("Running ALDEx2 with two groups. Performing t-test...")

          # Creating ALDEx2 object for two groups
          ALDEx2_object <- ALDEx2::aldex.clr(
            ALDEx2_abundance,
            Group,
            mc.samples = 256,
            denom = "all",
            verbose = FALSE
          )

          # Retrieving t-test results
          ALDEx2_results <- ALDEx2::aldex.ttest(ALDEx2_object,
                                                paired.test = FALSE,
                                                verbose = FALSE)
          # Building the results dataframe
          p_values_df <- data.frame(
            feature = rep(rownames(ALDEx2_results), 2),
            method = c(
              rep("ALDEx2_Welch's t test", nrow(ALDEx2_results)),
              rep("ALDEx2_Wilcoxon rank test", nrow(ALDEx2_results))
            ),
            group1 = rep(Level[1], 2 * nrow(ALDEx2_results)),
            group2 = rep(Level[2], 2 * nrow(ALDEx2_results)),
            p_values = c(ALDEx2_results$we.ep, ALDEx2_results$wi.ep)
          )

          message("ALDEx2 analysis with two groups complete.")

        } else {
          message("Running ALDEx2 with multiple groups. This might take some time, please wait patiently...")

          # Creating ALDEx2 object for multiple groups
          ALDEx2_object <- ALDEx2::aldex.clr(
            ALDEx2_abundance,
            Group,
            mc.samples = 256,
            denom = "all",
            verbose = FALSE
          )

          # Retrieving Kruskal-Wallis test results
          ALDEx2_results <- ALDEx2::aldex.kw(ALDEx2_object)

          # Building the initial results dataframe
          p_values_df <- data.frame(
            feature = rep(rownames(ALDEx2_results), 2),
            method = c(
              rep("ALDEx2_Kruskal-Wallace test", nrow(ALDEx2_results)),
              rep("ALDEx2_glm test", nrow(ALDEx2_results))
            ),
            p_values = c(ALDEx2_results$kw.ep, ALDEx2_results$glm.ep)
          )

          # Loop through levels to complete the results dataframe
          for (k in 1:length_Level) {
            p_values_df <- cbind(p_values_df[, 1:(1 + k)],
                                 rep(Level[k], 2 * nrow(ALDEx2_results)),
                                 p_values_df$p_values)
          }
          colnames(p_values_df)[3:(3 + length_Level)] <- c(paste0("group", 1:length_Level), "p_values")

          message("ALDEx2 analysis with multiple groups complete.")
        }
      },
      "DESeq2" = {
        # Inform the user about the requirements of DESeq2
        message("Running DESeq2. Note: DESeq2 is only suitable for comparison between two groups.")

        # Rename columns in metadata for DESeq2
        DESeq2_metadata <- metadata_df
        DESeq2_abundance_mat <- abundance_mat
        DESeq2_colnames <- colnames(DESeq2_metadata)
        DESeq2_colnames[DESeq2_colnames == group] <- "Group_group_nonsense"
        colnames(DESeq2_metadata) <- DESeq2_colnames
        DESeq2_metadata[, "Group_group_nonsense"] <- factor(DESeq2_metadata[, "Group_group_nonsense"])

        # Generate combinations of groups for comparison
        DESeq2_combinations <- utils::combn(unique(DESeq2_metadata[, "Group_group_nonsense"]), 2)
        DESeq2_results <- list()

        # Loop through combinations and perform DESeq2 analysis
        message("Performing pairwise comparisons with DESeq2...")
        for (i in seq_len(ncol(DESeq2_combinations))) {
          j <- DESeq2_combinations[, i]

          # Subsetting the data for the current combination of groups
          DESeq2_sub_group <- DESeq2_metadata$Group_group_nonsense %in% j
          DESeq2_metadata_sub <- DESeq2_metadata[DESeq2_sub_group,]
          DESeq2_abundance_mat_sub <- DESeq2_abundance_mat[, DESeq2_sub_group]
          DESeq2_abundance_mat_sub <- round(DESeq2_abundance_mat_sub)

          # Creating DESeq2 object and performing analysis
          DESeq2_object <- DESeq2::DESeqDataSetFromMatrix(
            countData = DESeq2_abundance_mat_sub,
            colData = DESeq2_metadata_sub,
            design = ~ Group_group_nonsense
          )
          DESeq2_object <- BiocGenerics::estimateSizeFactors(DESeq2_object, type = "poscounts")
          DESeq2_object_finish <- DESeq2::DESeq(DESeq2_object)
          DESeq2_results[[i]] <- DESeq2::results(DESeq2_object_finish)
        }

        # Compile the results into a data frame
        message("Compiling DESeq2 results...")
        DESeq2_results_nrow <- unlist(lapply(DESeq2_results, function(x) nrow(x)))
        DESeq2_results_matrix <- as.matrix(do.call(rbind, DESeq2_results))
        DESeq2_combinations_matrix_t <- t(DESeq2_combinations)
        DESeq2_group_matrix <- matrix(ncol = 2)

        # Loop through to build the final results matrix
        for (i in 1:length(DESeq2_results_nrow)) {
          DESeq2_group_matrix <- rbind(matrix(
            rep(DESeq2_combinations_matrix_t[i,], times = DESeq2_results_nrow[i]),
            ncol = 2,
            byrow = TRUE
          ), DESeq2_group_matrix)
        }

        # Cleanup and format the results matrix
        DESeq2_group_matrix <- DESeq2_group_matrix[stats::complete.cases(DESeq2_group_matrix),]
        colnames(DESeq2_group_matrix) <- c("group1", "group2")
        p_values_matrix <- cbind(
          feature = rownames(DESeq2_results_matrix),
          method = "DESeq2",
          DESeq2_group_matrix,
          p_values = as.vector(DESeq2_results_matrix[, "pvalue"])
        )
        p_values_df <- as.data.frame(p_values_matrix)

        message("DESeq2 analysis complete.")
      },
      "Maaslin2" = {
        # Inform the user that Maaslin2 is starting
        message("Running Maaslin2 analysis...")

        # Preparing the data for Maaslin2
        Maaslin2_abundance_mat <- abundance_mat
        Maaslin2_abundance_mat <- t(Maaslin2_abundance_mat)
        Maaslin2_metadata_df <- metadata_df
        rownames(Maaslin2_metadata_df) <- Maaslin2_metadata_df[, matching_columns]
        Maaslin2_metadata_df <- dplyr::select(Maaslin2_metadata_df, -matching_columns)

        # Perform Maaslin2 analysis based on the number of levels
        message("Performing Maaslin2 analysis...")
        switch(length_Level == 2,
               "TRUE" = {
                 Maaslin2_results <- Maaslin2::Maaslin2(
                   Maaslin2_abundance_mat,
                   Maaslin2_metadata_df,
                   output = paste0("Maaslin2_results_", group),
                   transform = "AST",
                   fixed_effects = group,
                   reference = Level,
                   normalization = "TSS",
                   standardize = TRUE,
                   min_prevalence = 0.1
                 )
                 p_values_matrix <- cbind(
                   feature = Maaslin2_results$results$feature,
                   method = "Maaslin2",
                   group1 = Level[1],
                   group2 = Level[2],
                   p_values = Maaslin2_results$results$pval
                 )
                 p_values_df <- as.data.frame(p_values_matrix)
               },
               {
                 Maaslin2_results <- Maaslin2::Maaslin2(
                   Maaslin2_abundance_mat,
                   Maaslin2_metadata_df,
                   output = paste0("Maaslin2_results_", group),
                   transform = "AST",
                   fixed_effects = group,
                   reference = paste0(group, ",", reference),
                   normalization = "TSS",
                   standardize = TRUE,
                   min_prevalence = 0.1
                 )
                 p_values_matrix <- cbind(
                   feature = Maaslin2_results$results$feature,
                   method = "Maaslin2",
                   group1 = Maaslin2_results$results$value,
                   group2 = reference,
                   p_values = Maaslin2_results$results$pval
                 )
                 p_values_df <- as.data.frame(p_values_matrix)
               })

        # Inform the user of where to find the results
        message(paste0(
          "Maaslin2 analysis complete. You can view the full analysis results and logs in the current default file location: ",
          getwd(),
          "/Maaslin2_results_",
          group
        ))
      },
      "LinDA" = {
        # Inform the user that LinDA is starting
        message("Running LinDA analysis...")

        # Preparing the data for LinDA
        LinDA_metadata_df <- metadata_df
        LinDA_colnames <- colnames(LinDA_metadata_df)
        LinDA_colnames[LinDA_colnames == group] <- "Group_group_nonsense_"
        colnames(LinDA_metadata_df) <- LinDA_colnames
        rownames(LinDA_metadata_df) <- LinDA_metadata_df[, matching_columns]
        LinDA_metadata_df <- dplyr::select(LinDA_metadata_df, -matching_columns)
        LinDA_metadata_df$Group_group_nonsense_ <- factor(LinDA_metadata_df$Group_group_nonsense_)

        # Check conditions for LinDA analysis
        if (length_Level != 2) {
          if (is.null(reference)) {
            stop("Error: A reference group is required when using LinDA or limma voom for comparisons among more than two groups. Please specify a reference group.")
          }
          LinDA_metadata_df$Group_group_nonsense_ <- stats::relevel(LinDA_metadata_df$Group_group_nonsense_, ref = reference)
        }

        # Perform LinDA analysis
        message("Performing LinDA analysis...")
        LinDA_results <- MicrobiomeStat::linda(abundance,
                                               LinDA_metadata_df,
                                               formula = "~Group_group_nonsense_",
                                               alpha = 0.05,
        )$output

        # Processing LinDA results
        message("Processing LinDA results...")
        if (length_Level != 2){
          for (i in 1:length(LinDA_results)) {
            LinDA_results[[i]] <- cbind(
              feature = rownames(LinDA_results[[i]]),
              method = "LinDA",
              group1 = substr(names(LinDA_results)[i], 22, stop = nchar(names(LinDA_results)[i])),
              group2 = reference,
              p_values = LinDA_results[[i]]$pvalue
            )
          }
        }else{
          for (i in 1:length(LinDA_results)) {
            LinDA_results[[i]] <- cbind(
              feature = rownames(LinDA_results[[i]]),
              method = "LinDA",
              group1 = Level[1],
              group2 = Level[2],
              p_values = LinDA_results[[i]]$pvalue
            )
          }
        }

        # Inform the user that LinDA analysis is complete
        message("LinDA analysis is complete.")

        # Combining the results
        p_values_matrix <- as.matrix(do.call(rbind, LinDA_results))
        p_values_df <- as.data.frame(p_values_matrix)
      },
      "edgeR" = {
        message("Processing data with edgeR method...")

        # Rounding the abundance matrix values
        edgeR_abundance_mat <- round(abundance_mat)

        # Initializing edgeR object
        message("Initializing edgeR object...")
        edgeR_object <- edgeR::DGEList(counts = edgeR_abundance_mat, group = Group)

        # Normalization
        message("Calculating normalization factors...")
        edgeR_object <- edgeR::calcNormFactors(edgeR_object)

        # Estimating common dispersion
        message("Estimating common dispersions...")
        edgeR_object <- edgeR::estimateCommonDisp(edgeR_object, verbose = TRUE)

        # Check if there are only two levels
        switch(length_Level == 2,
               "TRUE" = {
                 # For two groups, perform exact test
                 message("Performing exact test for two groups...")
                 edgeR_results <- edgeR::exactTest(edgeR_object, pair = c(1, 2))

                 # Extracting and structuring results
                 p_values_matrix <-
                   cbind(
                     feature = rownames(edgeR_abundance_mat),
                     method = "edgeR",
                     group1 = Level[1],
                     group2 = Level[2],
                     p_values = edgeR_results$table$PValue
                   )
                 p_values_df <- as.data.frame(p_values_matrix)
               },
               {
                 # For more than two groups, perform multiple tests
                 message("Performing multiple exact tests for multiple groups...")

                 edgeR_results <- list()
                 edgeR_numbers <- seq_along(Level)
                 edgeR_combinations <- utils::combn(edgeR_numbers, 2)

                 for (j in 1:sum(1:(length_Level - 1))) {
                   message(paste("Processing combination", j, "of", sum(1:(length_Level - 1)), "..."))
                   DGEExact <- edgeR::exactTest(edgeR_object, pair = edgeR_combinations[, j])
                   edgeR_results[[j]] <- cbind(
                     feature = rownames(DGEExact$table),
                     method = "edgeR",
                     group1 = DGEExact$comparison[1],
                     group2 = DGEExact$comparison[2],
                     p_values = DGEExact$table$PValue
                   )
                 }

                 # Combining all results
                 message("Combining results...")
                 p_values_matrix <- as.matrix(do.call(rbind, edgeR_results))
                 p_values_df <- as.data.frame(p_values_matrix)
               })

        message("edgeR processing completed.")
      },
      "limma voom" = {
        # Inform the user that limma voom analysis is starting
        message("Starting limma voom analysis...")

        # Preparing the data for limma voom analysis
        limma_voom_abundance_mat <- round(abundance_mat)

        # Check the number of levels and set the reference group accordingly
        if (length_Level != 2) {
          if (is.null(reference)) {
            stop("Error: A reference group is required when using limma voom for comparisons among more than two groups. Please specify a reference group.")
          }
          Group <- stats::relevel(Group, ref = reference)
          message(paste0("Using '", reference, "' as the reference group for limma voom analysis."))
        } else {
          reference <- Level[1]
          Group <- stats::relevel(Group, ref = reference)
          message(paste0("Using '", reference, "' as the default reference group for limma voom analysis with two levels."))
        }

        # Perform limma voom analysis
        message("Performing limma voom analysis...")
        limma_voom_object <- edgeR::DGEList(counts = limma_voom_abundance_mat, group = Group)
        limma_voom_object <- edgeR::calcNormFactors(limma_voom_object)
        limma_voom_Fit <- limma::lmFit(limma::voom(limma_voom_object))
        limma_voom_Fit <- limma::eBayes(limma_voom_Fit)

        # Collect p-values
        p_values_matrix <- cbind(
          feature = rep(rownames(abundance), length_Level - 1),
          method = "limma voom",
          group1 = reference,
          group2 = rep(substr(
            colnames(limma_voom_Fit$p.value)[-1], 6, nchar(colnames(limma_voom_Fit$p.value)[-1])
          ), each = nrow(abundance)),
          p_values = c(limma_voom_Fit$p.value[,-1])
        )

        # Convert to a data frame
        p_values_df <- as.data.frame(p_values_matrix)

        # Inform the user that the analysis is complete
        message("Limma voom analysis complete.")
      },
      "metagenomeSeq" = {
        # Inform the user that metagenomeSeq analysis is starting
        message("Starting metagenomeSeq analysis...")

        # Preparing the data for metagenomeSeq analysis
        metagenomeSeq_combinations <- utils::combn(Level, 2)
        metagenomeSeq_results_list <- list()

        # Loop through the combinations
        for (i in seq_len(ncol(metagenomeSeq_combinations))) {
          message(paste0("Processing combination ", i, " of ", ncol(metagenomeSeq_combinations), "..."))

          # Subset the data for the current combination
          sub_metadata_df <- metadata_df[metadata_df[, group] %in% metagenomeSeq_combinations[, i], ]
          sub_abundance <- abundance[, metadata_df[, group] %in% metagenomeSeq_combinations[, i]]

          # Perform metagenomeSeq analysis
          metagenomeSeq_list <- list(sub_abundance, as.data.frame(rownames(sub_abundance)))
          metagenomeSeq_metadata_df <- sub_metadata_df
          rownames(metagenomeSeq_metadata_df) <- metagenomeSeq_metadata_df[, matching_columns]
          metagenomeSeq_metadata_df <- dplyr::select(metagenomeSeq_metadata_df, -matching_columns)
          metagenomeSeq_colnames <- colnames(metagenomeSeq_metadata_df)
          metagenomeSeq_colnames[metagenomeSeq_colnames == group] <- "Group_group_nonsense_"
          colnames(metagenomeSeq_metadata_df) <- metagenomeSeq_colnames
          metagenomeSeq_phenotypeData <- Biobase::AnnotatedDataFrame(metagenomeSeq_metadata_df)
          metagenomeSeq_object <- metagenomeSeq::newMRexperiment(metagenomeSeq_list[[1]], phenoData = metagenomeSeq_phenotypeData)
          metagenomeSeq_object <- metagenomeSeq::cumNorm(metagenomeSeq_object, p = metagenomeSeq::cumNormStatFast(metagenomeSeq_object))
          metagenomeSeq_mod <- stats::model.matrix(~ Group_group_nonsense_, data = Biobase::pData(metagenomeSeq_object))
          metagenomeSeq_results <- metagenomeSeq::fitFeatureModel(metagenomeSeq_object, metagenomeSeq_mod)

          # Store the results
          metagenomeSeq_results_list[[i]] <- cbind(
            feature = names(metagenomeSeq_results@pvalues),
            method = "metagenomeSeq",
            group1 = metagenomeSeq_combinations[, i][1],
            group2 = metagenomeSeq_combinations[, i][2],
            p_values = as.vector(metagenomeSeq_results@pvalues)
          )
        }

        # Combine the results into a matrix
        message("Combining results...")
        p_values_matrix <- as.matrix(do.call(rbind, metagenomeSeq_results_list))

        # Convert to a data frame
        p_values_df <- as.data.frame(p_values_matrix)

        # Inform the user that the analysis is complete
        message("MetagenomeSeq analysis complete.")
      }
      ,
      "Lefser" = {
        # Inform the user that Lefser analysis is starting
        message("Starting Lefser analysis...")

        # Preparing the data for Lefser analysis
        Lefser_combinations <- utils::combn(Level, 2)
        Lefser_results <- list()
        Lefser_metadata_df <- metadata_df

        # Loop through the combinations
        for (i in seq_len(ncol(Lefser_combinations))) {
          message(paste0("Processing combination ", i, " of ", ncol(Lefser_combinations), "..."))

          # Subset the data for the current combination
          Lefser_sub_metadata_df <- Lefser_metadata_df[Lefser_metadata_df[, group] %in% Lefser_combinations[, i], ]
          Lefser_sub_abundance <- abundance[, Lefser_metadata_df[, group] %in% Lefser_combinations[, i]]
          colnames(Lefser_sub_metadata_df)[colnames(Lefser_sub_metadata_df) == group] <- "Group_group_nonsense_"
          Lefser_sub_metadata_df$Group_group_nonsense_ <- factor(Lefser_sub_metadata_df$Group_group_nonsense_)

          # Perform Lefser analysis
          Lefser_object <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(Lefser_sub_abundance)),
                                                                      colData = Lefser_sub_metadata_df)
          Lefser_kw_filter <- apply(Lefser_object@assays@data$counts, 1L, function(x) {
            stats::kruskal.test(x ~ as.numeric(Lefser_sub_metadata_df[, "Group_group_nonsense_"]) -
                                  1)[["p.value"]]
          })
          if (!sum(stats::na.omit(as.numeric(Lefser_kw_filter < 0.05)))) {
            next
          }
          Lefser_results[[i]] <- cbind(
            feature = lefser::lefser(Lefser_object, groupCol = "Group_group_nonsense_")$Names,
            method = "Lefser",
            group1 = Lefser_combinations[, i][1],
            group2 = Lefser_combinations[, i][2],
            effect_scores = lefser::lefser(Lefser_object, groupCol = "Group_group_nonsense_")$scores
          )
        }

        # Combine the results into a matrix
        message("Combining results...")
        p_values_matrix <- as.matrix(do.call(rbind, Lefser_results))

        # Convert to a data frame
        p_values_df <- as.data.frame(p_values_matrix)

        # Inform the user that Lefser analysis is complete
        message("Lefser analysis complete.")

        # Warn the user that Lefser results are not suitable for pathway_errorbar visualization
        message("Warning: Lefser results are not suitable for visualization with pathway_errorbar, as Lefser does not support outputting p-values.")

        # Return the results
        return(p_values_df)
      }
    )

    valid_p_adjust <- c("BH", "holm", "bonferroni",  "hochberg", "fdr", "none")

    if (!p.adjust %in% valid_p_adjust) {
      stop(paste("Invalid p.adjust method. Please choose from:",
                 paste(valid_p_adjust, collapse = ", ")))
    }

      switch(
        p.adjust,
        "BH" = {
          adjusted_p_values <- p.adjust(p_values_df$p_values, method = "BH")
        },
        "holm" = {
          adjusted_p_values <- p.adjust(p_values_df$p_values, method = "holm")
        },
        "bonferroni" = {
          adjusted_p_values <-
            p.adjust(p_values_df$p_values, method = "bonferroni")
        },
        "hochberg" = {
          adjusted_p_values <-
            p.adjust(p_values_df$p_values, method = "hochberg")
        },
        "fdr" = {
          adjusted_p_values <-
            p.adjust(p_values_df$p_values, method = "fdr")
        },
        "none" = {
          adjusted_p_values <-
            p.adjust(p_values_df$p_values, method = "none")
        }
      )
      daa_results_df <-
        cbind(p_values_df, adj_method = p.adjust, p_adjust = adjusted_p_values)
      return(daa_results_df)
  }
