#' Title
#'
#' @param abundance
#' @param metadata
#' @param group
#' @param daa_method
#' @param select
#' @param p.adjust
#' @param reference
#'
#' @return
#' @export
#'
#' @examples
pathway_daa <-
  function(abundance,
           metadata,
           group,
           daa_method,
           select = NULL,
           p.adjust = "BH",
           reference = NULL) {
    if (!is_tibble(metadata)) {
      metadata <- as_tibble(metadata)
    }
    sample_names <- colnames(abundance)
    matches <-
      lapply(metadata, function(x) {
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
    Level <- names(table(Group))
    switch(
      daa_method,
      "ALDEx2" = {
        ALDEx2_abundance <- round(abundance)
        switch(length(Level) == 2,
               "TRUE" = {
                 ALDEx2_object <-
                   aldex.clr(
                     ALDEx2_abundance,
                     Group,
                     mc.samples = 256,
                     denom = "all",
                     verbose = F
                   )
                 ALDEx2_results <-
                   aldex.ttest(ALDEx2_object,
                               paired.test = FALSE,
                               verbose = FALSE)
                 p_values_df <-
                   data.frame(
                     feature = rep(rownames(ALDEx2_results), 2),
                     method = c(
                       rep("Welch’s t test", nrow(ALDEx2_results)),
                       rep("Wilcoxon rank test", nrow(ALDEx2_results))
                     ),
                     group1 = rep(names(table(Group))[1], 2 * nrow(ALDEx2_results)),
                     group2 = rep(names(table(Group))[2], 2 * nrow(ALDEx2_results)),
                     p_values = c(ALDEx2_results$we.ep, ALDEx2_results$wi.ep)
                   )
               },
               {
                 messgae("ALDEx2 takes a long time to complete the calculation, please wait patiently.")
                 ALDEx2_object <-
                   aldex.clr(
                     ALDEx2_abundance,
                     Group,
                     mc.samples = 256,
                     denom = "all",
                     verbose = FALSE
                   )
                 ALDEx2_results <- aldex.kw(ALDEx2_object)
                 p_values_df <-
                   data.frame(
                     feature = rep(rownames(ALDEx2_results), 2),
                     method = c(
                       rep("Kruskal-Wallace test", nrow(ALDEx2_results)),
                       rep("glm test", nrow(ALDEx2_results))
                     ),
                     p_values = c(ALDEx2_results$kw.ep, ALDEx2_results$glm.ep)
                   )
                 for (k in 1:length(Level)) {
                   p_values_df <-
                     cbind(p_values_df[, 1:(1 + k)],
                           rep(names(table(Group))[k], 2 * nrow(ALDEx2_results)),
                           p_values_df$p_values)
                 }
                 colnames(p_values_df)[3:(3 + length(Level))] <-
                   c(paste0("group", 1:length(Level)), "p_values")
               })
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
          combn(unique(DESeq2_metadata[, "Group_group_nonsense"]), 2)
        DESeq2_results <- list()
        for (i in 1:ncol(DESeq2_combinations)) {
          j <- DESeq2_combinations[, i]
          DESeq2_sub_group <-
            DESeq2_metadata$Group_group_nonsense %in% j
          DESeq2_metadata_sub <- DESeq2_metadata[DESeq2_sub_group, ]
          DESeq2_abundance_mat_sub <-
            DESeq2_abundance_mat[, DESeq2_sub_group]
          DESeq2_object <-
            DESeqDataSetFromMatrix(
              countData = DESeq2_abundance_mat_sub,
              colData = DESeq2_metadata_sub,
              design = ~ Group_group_nonsense
            )
          DESeq2_object <-
            estimateSizeFactors(DESeq2_object, type = "poscounts")
          DESeq2_object_finish <- DESeq(DESeq2_object)
          DESeq2_results[[i]] <- results(DESeq2_object_finish)
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
              rep(DESeq2_combinations_matrix_t[i, ], times = DESeq2_results_nrow[i]),
              ncol = 2,
              byrow = TRUE
            ), DESeq2_group_matrix)
        }
        DESeq2_group_matrix <-
          DESeq2_group_matrix[complete.cases(DESeq2_group_matrix), ]
        colnames(DESeq2_group_matrix) <- c("group1", "group2")
        p_values_matrix <-
          cbind(
            feature = rownames(DESeq2_results_matrix),
            DESeq2_group_matrix,
            p_values = as.vector(DESeq2_results_matrix[, "pvalue"])
          )
        p_values_df <- as.data.frame(p_values_matrix)
      },
      "Maaslin2" = {
        Maaslin2_abundance_mat <- abundance_mat
        Maaslin2_abundance_mat <- t(Maaslin2_abundance_mat)
        Maaslin2_metadata_df <- metadata_df
        rownames(Maaslin2_metadata_df) <-
          Maaslin2_metadata_df[, matching_columns]
        Maaslin2_metadata_df <-
          select(Maaslin2_metadata_df,-matching_columns)
        switch(length(Level) == 2,
               "TRUE" = {
                 Maaslin2_results <- Maaslin2(
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
                     group1 = names(table(Group))[1],
                     group2 = names(table(Group))[2],
                     p_values = Maaslin2$results$pval
                   )
                 p_values_df <- as.data.frame(p_values_matrix)
               },
               {
                 Maaslin2_results <- Maaslin2(
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
                     group1 = Maaslin2_results$results$value,
                     group2 = reference,
                     p_values = Maaslin2$results$pval
                   )
                 p_values_df <- as.data.frame(p_values_matrix)
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
          select(LinDA_metadata_df,-matching_columns)
        LinDA_metadata_df$Group_group_nonsense_ <-
          factor(LinDA_metadata_df$Group_group_nonsense_)
        LinDA_metadata_df$Group_group_nonsense_ <-
          relevel(LinDA_metadata_df$Group_group_nonsense_, ref = reference)
        LinDA_results <- MicrobiomeStat::linda(abundance,
                                               LinDA_metadata_df,
                                               formula = "~Group_group_nonsense_",
                                               alpha = 0.05,)$output
        for (i in 1:length(LinDA_results)){
          LinDA_results[[i]] <- cbind(feature = rownames(LinDA_results[[i]]),method= "LinDA",group1 = substr(names(LinDA_results)[i],22,stop=nchar(names(LinDA_results)[i])),group2 = reference,p_values = LinDA_results[[i]]$pvalue)
        }

        p_values_matrix <- as.matrix(do.call(rbind, LinDA_results))
        p_values_df <- as.data.frame(p_values_matrix)
      },
      "edgeR" = {
        edgeR_abundance_mat <- t(round(abundance_mat))
        edgeR_object <-
          DGEList(counts = edgeR_abundance_mat, group = Group)
        edgeR_object <- calcNormFactors(edgeR_object)
        edgeR_object <-
          estimateCommonDisp(edgeR_object, verbose = T)
        length_Level <- length(Level)
        switch(length_Level == 2,
               "TRUE" = {
                 edgeR_results <- exactTest(edgeR_object, pair = c(1, 2))
                 p_values <-
                   as.matrix(edgeR_results$table)[,-c(1, 2)]
               },
               {
                 edgeR_results <- list()
                 numbers <- seq_along(Level)
                 combinations <- combn(numbers, 2)
                 for (j in 1:sum(1:(length(Level) - 1))) {
                   DGEExact <- exactTest(edgeR_object, pair = combinations[, j])
                   edgeR_results[[j]] <- tibble(
                     feature = rownames(DGEExact$table),
                     Group1 = DGEExact$comparison[1],
                     Group2 = DGEExact$comparison[2],
                     PValue = DGEExact$table$PValue
                   )
                 }
                 edgeR_results <- bind_rows(edgeR_results)
               })
      },
      "limma voom" = {
        edgeR_abundance_mat <- round(abundance_mat)
        edgeR_object <-
          DGEList(counts = edgeR_abundance_mat, group = Group)
        edgeR_object <- calcNormFactors(edgeR_object)
        limma_voom_Fit <- lmFit(voom(edgeR_object))
        limma_voom_Fit <- eBayes(limma_voom_Fit)
        p_values <- limma_voom_Fit$p.value[, 2]
      },
      "metagenomeSeq" = {
        if (length(Level) != 2) {
          stop("MetagenomeSeq only support two groups comparison.")
        }
        metagenomeSeq_list <- list()
        metagenomeSeq_list[[1]] <- abundance
        metagenomeSeq_list[[2]] <-
          as.data.frame(rownames(abundance))
        metagenomeSeq_metadata_df <- metadata_df
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
          AnnotatedDataFrame(metagenomeSeq_metadata_df)
        metagenomeSeq_taxa <- data.frame(OTU = rownames(abundance))
        rownames(metagenomeSeq_taxa) <- rownames(abundance)
        metagenomeSeq_OTUdata <-
          AnnotatedDataFrame(metagenomeSeq_taxa)
        metagenomeSeq_object <-
          newMRexperiment(metagenomeSeq_list[[1]], phenoData = metagenomeSeq_phenotypeData)
        metagenomeSeq_object <-
          cumNorm(metagenomeSeq_object, p = cumNormStatFast(metagenomeSeq_object))
        metagenomeSeq_mod <-
          model.matrix( ~ Group_group_nonsense_, data = pData(metagenomeSeq_object))
        metagenomeSeq_results <-
          fitFeatureModel(metagenomeSeq_object, metagenomeSeq_mod)
        p_values <- metagenomeSeq_results@pvalues
      },
      "Lefser" = {
        # Lefser only for two groups.
        if (length(Level) != 2) {
          stop("Lefser only support two groups comparison.")
        }
        BiocManager::install("lefser")
        Lefser_object <-
          SummarizedExperiment(assays = list(counts = abundance), colData = metadata_df)
        Lefser_results <- lefser(se1, groupCol = "Enviroment")
        lefserPlot(
          Lefser_results,
          colors = c("#7fb1d3", "#fdb462"),
          trim.names = FALSE
        )
      }
    )




    switch(
      p.adjust,
      "BH" = {
        adjusted_p_values <- p.adjust(p_values, method = "BH")
      },
      "holm" = {
        adjusted_p_values <- p.adjust(p_values, method = "holm")
      },
      "bonferroni" = {
        adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
      },
      "Bayes" = {
        library(rstan)
        model_string <- "
  data {
    int<lower=0> N;
    real p[N];
  }
  parameters {
    real<lower=0,upper=1> mu;
  }
  model {
    for (i in 1:N) {
      p[i] ~ beta(mu, 1);
    }
  }
"
        model <- stan_model(model_code = model_string)
        posterior <-
          sampling(model, data = list(N = length(p_values), p = p_values))
      },
      "Random Effect" = {
        # 加载必要的包
        library(glmnet)
        # 将 p 值转换为二进制类别
        y <- ifelse(p_values < 0.05, 1, 0)
        # 建立模型
        model <- cv.glmnet(x, y, family = "binomial")
        # 计算系数
        coefs <- coef(model, s = "lambda.min")
      }
    )
  }
