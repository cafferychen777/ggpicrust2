#' Predictional functional patwhay differential abundance (DA)
#'
#' @param abundance a data frame containing predicted functional pathway abundance
#' @param metadata a tibble containing samples information
#' @param group a character specifying the group name for differential abundance analysis
#' @param daa_method a character specifying the method for differential abundance analysis, default is "ALDEx2"
#' @param select a vector containing pathway names for analysis, if NULL all pathways are included, default is NULL
#' @param p.adjust a character specifying the method for p-value adjustment, default is "BH"
#' @param reference a character specifying the reference group level, required for several differential abundance analysis methods such as LinDA, limme voom and Maaslin2, default is NULL
#'
#' @return  a data frame containing the differential abundance analysis results.
#' @export
#'
#' @examples
#' library(ggpicrust2)
#' library(MicrobiomeStat)
#' abundance <- data.frame(sample1 = c(10, 20, 30),
#' sample2 = c(20, 30, 40),
#' sample3 = c(30, 40, 50),
#' row.names = c("pathway1", "pathway2", "pathway3"))
#'
#' metadata <- tibble::tibble(sample = paste0("sample", 1:3),
#' group = c("control", "control", "treatment"))
#'
#' # Run pathway_daa function
#' # result <- pathway_daa(abundance = abundance, metadata = metadata, group = "group",
#' # daa_method = "LinDA")
pathway_daa <-
  function(abundance,
           metadata,
           group,
           daa_method = "ALDEx2",
           select = NULL,
           p.adjust = "BH",
           reference = NULL) {
    if (!tibble::is_tibble(metadata)) {
      metadata <- tibble::as_tibble(metadata)
    }
    sample_names <- colnames(abundance)
    matches <-
      base::lapply(metadata, function(x) {
        intersect(sample_names, x)
      })
    matching_columns <-
      names(metadata)[sapply(matches, function(x) {
        length(x) == length(sample_names)
      })]
    switch(is.null(select),
           "TRUE" = {
             abundance <- abundance
           },
           "FALSE" = {
             abundance <- abundance[, colnames(abundance) %in% select]
             metadata <-
               metadata[as.matrix(metadata[, matching_columns]) %in% select,]
           })
    sample_names <- colnames(abundance)
    abundance_mat <- as.matrix(abundance)
    metadata_order <-
      match(sample_names, as.matrix(metadata[, matching_columns]))
    metadata <- metadata[metadata_order,]
    metadata_mat <- as.matrix(metadata)
    metadata_df <- as.data.frame(metadata)
    Group <- factor(metadata_mat[, group])
    Level <- levels(Group)
    length_Level <- length(Level)
    switch(
      daa_method,
      "ALDEx2" = {
        ALDEx2_abundance <- round(abundance)
        if (length_Level == 2){
          ALDEx2_object <-
            ALDEx2::aldex.clr(
              ALDEx2_abundance,
              Group,
              mc.samples = 256,
              denom = "all",
              verbose = FALSE
            )
          ALDEx2_results <-
            ALDEx2::aldex.ttest(ALDEx2_object,
                        paired.test = FALSE,
                        verbose = FALSE)
          p_values_df <<-
            data.frame(
              feature = rep(rownames(ALDEx2_results), 2),
              method = c(
                rep("ALDEx2_Welch's t test", nrow(ALDEx2_results)),
                rep("ALDEx2_Wilcoxon rank test", nrow(ALDEx2_results))
              ),
              group1 = rep(Level[1], 2 * nrow(ALDEx2_results)),
              group2 = rep(Level[2], 2 * nrow(ALDEx2_results)),
              p_values = c(ALDEx2_results$we.ep, ALDEx2_results$wi.ep)
            )
        } else {
          #messgae("ALDEx2 takes a long time to complete the calculation, please wait patiently.")
          ALDEx2_object <-
            ALDEx2::aldex.clr(
              ALDEx2_abundance,
              Group,
              mc.samples = 256,
              denom = "all",
              verbose = FALSE
            )
          ALDEx2_results <- ALDEx2::aldex.kw(ALDEx2_object)
          p_values_df <<-
            data.frame(
              feature = rep(rownames(ALDEx2_results), 2),
              method = c(
                rep("ALDEx2_Kruskal-Wallace test", nrow(ALDEx2_results)),
                rep("ALDEx2_glm test", nrow(ALDEx2_results))
              ),
              p_values = c(ALDEx2_results$kw.ep, ALDEx2_results$glm.ep)
            )
          for (k in 1:length_Level) {
            p_values_df <<-
              cbind(p_values_df[, 1:(1 + k)],
                    rep(Level[k], 2 * nrow(ALDEx2_results)),
                    p_values_df$p_values)
          }
          colnames(p_values_df)[3:(3 + length_Level)] <-
            c(paste0("group", 1:length_Level), "p_values")
        }
      },
      "DESeq2" = {
        message("DESeq2 is only suitable for comparison between two groups.")
        DESeq2_metadata <- metadata_df
        DESeq2_abundance_mat <- abundance_mat
        DESeq2_colnames <- colnames(DESeq2_metadata)
        DESeq2_colnames[DESeq2_colnames == group] <-
          "Group_group_nonsense"
        colnames(DESeq2_metadata) <- DESeq2_colnames
        DESeq2_metadata[, "Group_group_nonsense"] <-
          factor(DESeq2_metadata[, "Group_group_nonsense"])
        DESeq2_combinations <-
          utils::combn(unique(DESeq2_metadata[, "Group_group_nonsense"]), 2)
        DESeq2_results <- list()
        for (i in seq_len(ncol(DESeq2_combinations))) {
          j <- DESeq2_combinations[, i]
          DESeq2_sub_group <-
            DESeq2_metadata$Group_group_nonsense %in% j
          DESeq2_metadata_sub <- DESeq2_metadata[DESeq2_sub_group,]
          DESeq2_abundance_mat_sub <-
            DESeq2_abundance_mat[, DESeq2_sub_group]
          DESeq2_abundance_mat_sub <- round(DESeq2_abundance_mat_sub)
          DESeq2_object <-
            DESeq2::DESeqDataSetFromMatrix(
              countData = DESeq2_abundance_mat_sub,
              colData = DESeq2_metadata_sub,
              design = ~ Group_group_nonsense
            )
          DESeq2_object <-
            BiocGenerics::estimateSizeFactors(DESeq2_object, type = "poscounts")
          DESeq2_object_finish <- DESeq2::DESeq(DESeq2_object)
          DESeq2_results[[i]] <- DESeq2::results(DESeq2_object_finish)
        }
        DESeq2_results_nrow <-
          unlist(lapply(DESeq2_results, function(x)
            nrow(x)))
        DESeq2_results_matrix <-
          as.matrix(do.call(rbind, DESeq2_results))
        DESeq2_combinations_matrix_t <- t(DESeq2_combinations)
        DESeq2_group_matrix <- matrix(ncol = 2)
        for (i in 1:length(DESeq2_results_nrow)) {
          DESeq2_group_matrix <-
            rbind(matrix(
              rep(DESeq2_combinations_matrix_t[i,], times = DESeq2_results_nrow[i]),
              ncol = 2,
              byrow = TRUE
            ), DESeq2_group_matrix)
        }
        DESeq2_group_matrix <-
          DESeq2_group_matrix[stats::complete.cases(DESeq2_group_matrix),]
        colnames(DESeq2_group_matrix) <- c("group1", "group2")
        p_values_matrix <-
          cbind(
            feature = rownames(DESeq2_results_matrix),
            method = "DESeq2",
            DESeq2_group_matrix,
            p_values = as.vector(DESeq2_results_matrix[, "pvalue"])
          )
        p_values_df <<- as.data.frame(p_values_matrix)
      },
      "Maaslin2" = {
        Maaslin2_abundance_mat <- abundance_mat
        Maaslin2_abundance_mat <- t(Maaslin2_abundance_mat)
        Maaslin2_metadata_df <- metadata_df
        rownames(Maaslin2_metadata_df) <-
          Maaslin2_metadata_df[, matching_columns]
        Maaslin2_metadata_df <-
          dplyr::select(Maaslin2_metadata_df, -matching_columns)
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
                 p_values_matrix <-
                   cbind(
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
                 p_values_matrix <-
                   cbind(
                     feature = Maaslin2_results$results$feature,
                     method = "Maaslin2",
                     group1 = Maaslin2_results$results$value,
                     group2 = reference,
                     p_values = Maaslin2_results$results$pval
                   )
                 p_values_df <<- as.data.frame(p_values_matrix)
               })
        message(
          paste0(
            "You can view the full analysis results and logs in the current default file location: ",
            getwd(),
            "/Maaslin2_results_",
            group
          )
        )
      },
      "LinDA" = {
        LinDA_metadata_df <- metadata_df
        LinDA_colnames <- colnames(LinDA_metadata_df)
        LinDA_colnames[LinDA_colnames == group] <-
          "Group_group_nonsense_"
        colnames(LinDA_metadata_df) <- LinDA_colnames
        rownames(LinDA_metadata_df) <-
          LinDA_metadata_df[, matching_columns]
        LinDA_metadata_df <-
          dplyr::select(LinDA_metadata_df, -matching_columns)
        LinDA_metadata_df$Group_group_nonsense_ <-
          factor(LinDA_metadata_df$Group_group_nonsense_)
        if (length_Level != 2) {
          if (is.null(reference)) {
            stop("If you use the LinDA or limma voom, you should give a reference.")
          }
          LinDA_metadata_df$Group_group_nonsense_ <-
            stats::relevel(LinDA_metadata_df$Group_group_nonsense_, ref = reference)
        }
        LinDA_results <- MicrobiomeStat::linda(abundance,
                                               LinDA_metadata_df,
                                               formula = "~Group_group_nonsense_",
                                               alpha = 0.05,
        )$output
        if (length_Level != 2){
          for (i in 1:length(LinDA_results)) {
            LinDA_results[[i]] <-
              cbind(
                feature = rownames(LinDA_results[[i]]),
                method = "LinDA",
                group1 = substr(names(LinDA_results)[i], 22, stop = nchar(names(
                  LinDA_results
                )[i])),
                group2 = reference,
                p_values = LinDA_results[[i]]$pvalue
              )
          }
        }else{
          for (i in 1:length(LinDA_results)) {
            LinDA_results[[i]] <-
              cbind(
                feature = rownames(LinDA_results[[i]]),
                method = "LinDA",
                group1 = Level[1],
                group2 = Level[2],
                p_values = LinDA_results[[i]]$pvalue
              )
          }
        }
        p_values_matrix <- as.matrix(do.call(rbind, LinDA_results))
        p_values_df <<- as.data.frame(p_values_matrix)
      },
      "edgeR" = {
        edgeR_abundance_mat <- round(abundance_mat)
        edgeR_object <-
          edgeR::DGEList(counts = edgeR_abundance_mat, group = Group)
        edgeR_object <- edgeR::calcNormFactors(edgeR_object)
        edgeR_object <-
          edgeR::estimateCommonDisp(edgeR_object, verbose = TRUE)
        switch(length_Level == 2,
               "TRUE" = {
                 edgeR_results <- edgeR::exactTest(edgeR_object, pair = c(1, 2))
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
                 edgeR_results <- list()
                 edgeR_numbers <- seq_along(Level)
                 edgeR_combinations <- utils::combn(edgeR_numbers, 2)
                 for (j in 1:sum(1:(length_Level - 1))) {
                   DGEExact <- edgeR::exactTest(edgeR_object, pair = edgeR_combinations[, j])
                   edgeR_results[[j]] <- cbind(
                     feature = rownames(DGEExact$table),
                     method = "edgeR",
                     group1 = DGEExact$comparison[1],
                     group2 = DGEExact$comparison[2],
                     p_values = DGEExact$table$PValue
                   )
                 }
                 p_values_matrix <-
                   as.matrix(do.call(rbind, edgeR_results))
                 p_values_df <<- as.data.frame(p_values_matrix)
               })
      },
      "limma voom" = {
        limma_voom_abundance_mat <- round(abundance_mat)
        if (length_Level != 2){
          if (is.null(reference)) {
            stop("If you use the LinDA or limma voom, you should give a reference.")
          }
          Group <- stats::relevel(Group, ref = reference)
        }else{
          reference <- Level[1]
          Group <- stats::relevel(Group, ref = reference)
        }
        limma_voom_object <-
          edgeR::DGEList(counts = limma_voom_abundance_mat, group = Group)
        limma_voom_object <-
          edgeR::calcNormFactors(limma_voom_object)
        limma_voom_Fit <- limma::lmFit(limma::voom(limma_voom_object))
        limma_voom_Fit <- limma::eBayes(limma_voom_Fit)
        p_values_matrix <-
          cbind(
            feature = rep(rownames(abundance), length_Level - 1),
            method = "limma voom",
            group1 = reference,
            group2 = rep(substr(
              colnames(limma_voom_Fit$p.value)[-1], 6, nchar(colnames(limma_voom_Fit$p.value)[-1])
            ), each = nrow(abundance)),
            p_values = c(limma_voom_Fit$p.value[,-1])
          )
        p_values_df <<- as.data.frame(p_values_matrix)
      },
      "metagenomeSeq" = {
        metagenomeSeq_combinations <- utils::combn(Level, 2)
        metagenomeSeq_results_list <- list()
        for (i in seq_len(ncol(metagenomeSeq_combinations))) {
          sub_metadata_df <-
            metadata_df[metadata_df[, group] %in% metagenomeSeq_combinations[, i],]
          sub_abundance <-
            abundance[, metadata_df[, group] %in% metagenomeSeq_combinations[, i]]
          metagenomeSeq_list <- list()
          metagenomeSeq_list[[1]] <- sub_abundance
          metagenomeSeq_list[[2]] <-
            as.data.frame(rownames(sub_abundance))
          metagenomeSeq_metadata_df <- sub_metadata_df
          rownames(metagenomeSeq_metadata_df) <-
            metagenomeSeq_metadata_df[, matching_columns]
          metagenomeSeq_metadata_df <-
            select(metagenomeSeq_metadata_df,-matching_columns)
          metagenomeSeq_colnames <-
            colnames(metagenomeSeq_metadata_df)
          metagenomeSeq_colnames[metagenomeSeq_colnames == group] <-
            "Group_group_nonsense_"
          colnames(metagenomeSeq_metadata_df) <-
            metagenomeSeq_colnames
          metagenomeSeq_phenotypeData <-
            Biobase::AnnotatedDataFrame(metagenomeSeq_metadata_df)
          metagenomeSeq_object <-
            metagenomeSeq::newMRexperiment(metagenomeSeq_list[[1]], phenoData = metagenomeSeq_phenotypeData)
          metagenomeSeq_object <-
            metagenomeSeq::cumNorm(metagenomeSeq_object, p = metagenomeSeq::cumNormStatFast(metagenomeSeq_object))
          metagenomeSeq_mod <-
            stats::model.matrix(~ Group_group_nonsense_, data = Biobase::pData(metagenomeSeq_object))
          metagenomeSeq_results <-
            metagenomeSeq::fitFeatureModel(metagenomeSeq_object, metagenomeSeq_mod)
          metagenomeSeq_results_list[[i]] <-
            cbind(
              feature = names(metagenomeSeq_results@pvalues),
              method = "metagenomeSeq",
              group1 = metagenomeSeq_combinations[, i][1],
              group2 = metagenomeSeq_combinations[, i][2],
              p_values = as.vector(metagenomeSeq_results@pvalues)
            )
        }
        p_values_matrix <-
          as.matrix(do.call(rbind, metagenomeSeq_results_list))
        p_values_df <<- as.data.frame(p_values_matrix)
      }
      ,
      "Lefser" = {
        Lefser_combinations <- utils::combn(Level, 2)
        Lefser_results <- list()
        Lefser_metadata_df <- metadata_df
        for (i in seq_len(ncol(Lefser_combinations))) {
          Lefser_sub_metadata_df <-
            Lefser_metadata_df[Lefser_metadata_df[, group] %in% Lefser_combinations[, i],]
          Lefser_sub_abundance <-
            abundance[, Lefser_metadata_df[, group] %in% Lefser_combinations[, i]]
          colnames(Lefser_sub_metadata_df)[colnames(Lefser_sub_metadata_df) == group] <-
            "Group_group_nonsense_"
          Lefser_sub_metadata_df$Group_group_nonsense_ <-
            factor(Lefser_sub_metadata_df$Group_group_nonsense_)
          Lefser_object <-
            SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(Lefser_sub_abundance)),
                                 colData = Lefser_sub_metadata_df)
          Lefser_kw_filter <-
            apply(Lefser_object@assays@data$counts, 1L, function(x) {
              stats::kruskal.test(x ~ as.numeric(Lefser_sub_metadata_df[, "Group_group_nonsense_"]) -
                             1)[["p.value"]]
            })
          if (!sum(stats::na.omit(as.numeric(Lefser_kw_filter < 0.05)))) {
            next
          }
          Lefser_results[[i]] <-
            cbind(
              feature = lefser::lefser(Lefser_object, groupCol = "Group_group_nonsense_")$Names,
              method = "Lefser",
              group1 = Lefser_combinations[, i][1],
              group2 = Lefser_combinations[, i][2],
              effect_scores = lefser::lefser(Lefser_object, groupCol = "Group_group_nonsense_")$scores
            )
        }
        p_values_matrix <-
          as.matrix(do.call(rbind, Lefser_results))
        p_values_df <<- as.data.frame(p_values_matrix)
        return(p_values_df)
      }
    )
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
