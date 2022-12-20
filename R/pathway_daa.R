

metadata <- read_delim("~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)



pathway_daa <-
  function(abundance,
           metadata,
           group,
           daa_method,
           select = c(),
           p.adjust = "BH",
           maaslin2_reference = NULL) {
    sample_names <- colnames(abundance)
    matches <-
      lapply(metadata, function(x)
        intersect(sample_names, x))
    matching_columns <-
      names(metadata)[sapply(matches, function(x)
        length(x) == length(sample_names))]
    switch(length(select),
           0 = {
             abundance <- abundance
           },
           {
             abundance <- abundance[, colnames(abundance) %in% select]
             metadata <-
               metadata[as.matrix(metadata[, matching_columns]) %in% select,]
           })
    sample_names <- colnames(abundance)
    abundance_mat <- as.matrix(abundance)
    metadata_order <-
      match(sample_names, as.matrix(metadata[, matching_columns]))
    metadata <- metadata[metadata_order, ]
    metadata_mat <- as.matrix(metadata)
    metadata_df <- as.data.frame(metadata)
    Group <- factor(as.matrix(metadata[, group]))
    Level <- names(table(Group))
    switch(
      daa_method,
      "ALDEx2" = {
        library(ALDEx2)
        abundance <- round(abundance)
        switch(length(Level),
               2 = {
                 clr_abundance <-
                   aldex.clr(
                     abundance,
                     Group,
                     mc.samples = 256,
                     denom = "all",
                     verbose = F
                   )
                 p_values <-
                   aldex.ttest(clr_abundance,
                               paired.test = FALSE,
                               verbose = FALSE)[, c(1, 3)]
               }, {
                 clr_abundance <-
                   aldex.clr(
                     abundance,
                     Group,
                     mc.samples = 32,
                     denom = "all",
                     verbose = F
                   )
                 p_values <- aldex.kw(clr_abundance)[, c(1, 3)]
               })
      },
      "ANCOMBC" = {
        library(SummarizedExperiment)
        assays = SimpleList(counts = b)
        metadata_df <- as.data.frame(metadata)
        tax_tab = matrix(rep("unknown", nrow(b) * 7),
                         nrow = nrow(b),
                         ncol = 7)
        rownames(tax_tab) = rownames(b)
        colnames(tax_tab) = c("Kingdom",
                              "Phylum",
                              "Class",
                              "Order",
                              "Family",
                              "Genus",
                              "Species")
        tax_tab = DataFrame(tax_tab)
        tse = TreeSummarizedExperiment(assays = assays,
                                       colData = metadata_df,
                                       rowData = tax_table)
        pseq = makePhyloseqFromTreeSummarizedExperiment(tse)
      },
      "DESeq2" = {
        if (!require("DESeq2")) {
          install.packages("DESeq2")
          require("DESeq2")
        }
        DESeq2_metadata <- metadata_df
        DESeq2_abundance_mat <- abundance_mat
        DESeq2_colnames <- colnames(DESeq2_metadata)
        DESeq2_colnames[DESeq2_colnames == group] = "Group_group_nonsense"
        colnames(DESeq2_metadata) <-  DESeq2_colnames
        DESeq2_object <-
          DESeqDataSetFromMatrix(
            countData = DESeq2_abundance_mat,
            colData = DESeq2_metadata,
            design = ~ Group_group_nonsense
          )
        DESeq2_object <-
          estimateSizeFactors(DESeq2_object, type = "poscounts")
        DESeq2_object_finish <- DESeq(DESeq2_object)
        DESeq2_results <- results(DESeq2_object_finish)
        DESeq2_results_tibble <- as_tibble(DESeq2_results)
        DESeq2_results_tibble <-
          DESeq2_results_tibble %>% add_column(.,
                                               pathway = rownames(DESeq2_results),
                                               .before = 1)
        DESeq2_results_tibble <- DESeq2_results_tibble[, c(1, 6)]
        DESeq2_results_mat <- as.matrix(DESeq2_results_tibble)
        rownames(DESeq2_results_mat) <- DESeq2_results_mat[, 1]
        DESeq2_results_mat <- DESeq2_results_mat[, -1]
        p_values <- DESeq2_results_mat
      },
      "Maaslin2" = {
        Maaslin2_abundance_mat <- abundance_mat
        Maaslin2_abundance_mat <- t(Maaslin2_abundance_mat)
        Maaslin2_metadata_df <- metadata_df
        rownames(Maaslin2_metadata_df) <-
          Maaslin2_metadata_df[, matching_columns]
        Maaslin2_metadata_df <- Maaslin2_metadata_df[, -1]
        switch(length(Level),
               2 = {
                 Maaslin2 <- Maaslin2(
                   Maaslin2_abundance_mat,
                   metadata_df,
                   output = "DAA example",
                   transform = "AST",
                   fixed_effects = group,
                   reference = names(table(metadata_df[, group])),
                   normalization = "TSS",
                   standardize = TRUE,
                   min_prevalence = 0.1
                 )
               }, {
                 Maaslin2 <- Maaslin2(
                   abundance_mat,
                   metadata_df,
                   output = "DAA example",
                   transform = "AST",
                   fixed_effects = group,
                   reference = paste0(group, maaslin2_reference),
                   #In case of multiple groups, be sure to specify the baseline reference
                   normalization = "TSS",
                   standardize = TRUE,
                   min_prevalence = 0.1
                 )
               })
        p_values <- Maaslin2$results$pval
      },
      "LinDA" = {
        LinDA_metadata_df <- metadata_df
        LinDA_colnames <- colnames(LinDA_metadata_df)
        LinDA_colnames[LinDA_colnames == group] = "Group_group_nonsense_"
        colnames(LinDA_metadata_df) <-  LinDA_colnames
        LinDA_results <- linda(
          abundance,
          LinDA_metadata_df,
          formula = '~Group_group_nonsense_',
          alpha = 0.05,
          prev.filter = 0,
          mean.abund.filter = 0
        )
        LinDA_results$output
      },
      "edgeR" = {
        library(edgeR)
        edgeR_abundance_mat <- t(round(abundance_mat))
        edgeR_object <-
          DGEList(counts = edgeR_abundance_mat, group = Group)
        edgeR_object <- calcNormFactors(edgeR_object)
        edgeR_object <-
          estimateCommonDisp(edgeR_object, verbose = T)
        length_Level <- length(Level)
        switch(length_Level,
               2 = {
                 edgeR_results <- exactTest(edgeR_object, pair = c(1, 2))
                 p_values <-
                   as.matrix(edgeR_results$table)[, -c(1, 2)]
               },
               {
                 edgeR_results <- list()
                 numbers <- 1:length(Level)
                 combinations <- combn(numbers, 2)
                 for (j in 1:sum(1:(length(Level)-1))){
                   DGEExact <- exactTest(edgeR_object, pair = combinations[,j])
                   edgeR_results[[j]] <- tibble(
                     pathway = rownames(DGEExact$table),
                     Group1 = DGEExact$comparison[1],
                     Group2 = DGEExact$comparison[2],
                     PValue = DGEExact$table$PValue
                   )
                 }
                 edgeR_results <- bind_rows(edgeR_results)
               })

      },
      "limma voom" = {
        library(edgeR)
        library(limma)
        edgeR_abundance_mat <- t(round(abundance_mat))
        edgeR_object <-
          DGEList(counts = edgeR_abundance_mat, group = Group)
        edgeR_object <- calcNormFactors(edgeR_object)
        limma_voom_Fit <- lmFit(voom(edgeR_object))
        limma_voom_Fit <- eBayes(limma_voom_Fit)
        p_values <- limma_voom_Fit$p.value
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
