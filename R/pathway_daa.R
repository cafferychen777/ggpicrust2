

# metadata <- read_delim("~/Microbiome/C9orf72/Code And Data/new_metadata.txt",delim = "\t", escape_double = FALSE,trim_ws = TRUE)


#abundance: df, rownames:pathway, colnames:samples
#metadata: tibble, sample_names single a column
pathway_daa <-
  function(abundance,
           metadata,
           group,
           daa_method,
           select = NULL,
           p.adjust = "BH",
           maaslin2_reference = NULL) {
    if (!is_tibble(metadata)){
      metadata <- as_tibble(metadata)
    }
    sample_names <- colnames(abundance)
    matches <-
      lapply(metadata, function(x)
        intersect(sample_names, x))
    matching_columns <-
      names(metadata)[sapply(matches, function(x)
        length(x) == length(sample_names))]
    length_select <- length(select)
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
    metadata <- metadata[metadata_order, ]
    metadata_mat <- as.matrix(metadata)
    metadata_df <- as.data.frame(metadata)
    Group <- factor(metadata_mat[,group])
    Level <- names(table(Group))
    switch(
      daa_method,
      "ALDEx2" = {
        if (!require("ALDEx2")) {
          install.packages("ALDEx2")
          require("ALDEx2")
        }
        ALDEx2_abundance <- round(abundance)
        switch(length(Level)==2,
               "TRUE" = {
                 ALDEx2_object <-
                   aldex.clr(
                     ALDEx2_abundance,
                     Group,
                     mc.samples = 256,
                     denom = "all",
                     verbose = F
                   )
                 p_values <-
                   aldex.ttest(ALDEx2_object,
                               paired.test = FALSE,
                               verbose = FALSE)[, c(1, 3)]
               }, {
                 ALDEx2_object <-
                   aldex.clr(
                     ALDEx2_abundance,
                     Group,
                     mc.samples = 256,
                     denom = "all",
                     verbose = F
                   )
                 p_values <- aldex.kw(ALDEx2_object)[, c(1, 3)]
               })
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
        if (!require("ALDEx2")) {
          if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install("Maaslin2")
          require("Maaslin2")
        }
        Maaslin2_abundance_mat <- abundance_mat
        Maaslin2_abundance_mat <- t(Maaslin2_abundance_mat)
        Maaslin2_metadata_df <- metadata_df
        rownames(Maaslin2_metadata_df) <-
          Maaslin2_metadata_df[, matching_columns]
        Maaslin2_metadata_df <- select(Maaslin2_metadata_df,-matching_columns)

        switch(length(Level)==2,
               "TRUE" = {
                 Maaslin2 <- Maaslin2(
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
               }, {
                 Maaslin2 <- Maaslin2(
                   Maaslin2_abundance_mat,
                   Maaslin2_metadata_df,
                   output = paste0("Maaslin2_results_", group),
                   transform = "AST",
                   fixed_effects = group,
                   reference = paste0(group,",",maaslin2_reference),
                   #In case of multiple groups, be sure to specify the baseline reference
                   normalization = "TSS",
                   standardize = TRUE,
                   min_prevalence = 0.1
                 )
               })
        message(
          paste0(
            "You can view the full analysis results and logs in the current default file location: ",
            getwd(),
            "/Maaslin2_results_",
            group
          )
        )
        p_values <- Maaslin2$results$pval
      },
      "LinDA" = {
        if (!require("LinDA")) {
          if (!require("devtools")) {
            install.packages("devtools")
            require("devtools")
          }
          devtools::install_github("zhouhj1994/LinDA")
          require("LinDA")
        }
        LinDA_metadata_df <- metadata_df
        LinDA_colnames <- colnames(LinDA_metadata_df)
        LinDA_colnames[LinDA_colnames == group] = "Group_group_nonsense_"
        colnames(LinDA_metadata_df) <-  LinDA_colnames
        rownames(LinDA_metadata_df) <-
          LinDA_metadata_df[, matching_columns]
        LinDA_metadata_df <- select(LinDA_metadata_df,-matching_columns)
        LinDA_results <- linda(
          abundance,
          LinDA_metadata_df,
          formula = '~Group_group_nonsense_',
          alpha = 0.05,
        )
        LinDA_results$output
        #在多组时需要整理输出结论
      },
      "edgeR" = {
        if (!require("edgeR")) {
          if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install("edgeR")
          require("edgeR")
        }
        edgeR_abundance_mat <- t(round(abundance_mat))
        edgeR_object <-
          DGEList(counts = edgeR_abundance_mat, group = Group)
        edgeR_object <- calcNormFactors(edgeR_object)
        edgeR_object <-
          estimateCommonDisp(edgeR_object, verbose = T)
        length_Level <- length(Level)
        switch(length_Level==2,
               "TRUE" = {
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
        if (!require("edgeR")) {
          if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install("edgeR")
          require("edgeR")
        }
        if (!require("limma")) {
          if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install("limma")
          require("limma")
        }
        edgeR_abundance_mat <- round(abundance_mat)
        edgeR_object <-
          DGEList(counts = edgeR_abundance_mat, group = Group)
        edgeR_object <- calcNormFactors(edgeR_object)
        limma_voom_Fit <- lmFit(voom(edgeR_object))
        limma_voom_Fit <- eBayes(limma_voom_Fit)
        p_values <- limma_voom_Fit$p.value[,2]
      },
      "metagenomeSeq" = {

      },
      "Lefser" = { #Lefser only for two groups.
        if (length(Level)!=2){
          stop("Lefser only support two groups comparison.")
        }
        if (!require("lefser")) {
          if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install("lefser")
          require("lefser")
        }
        Lefser_object <- SummarizedExperiment(assays = list(counts = abundance),colData = metadata_df)
        Lefser_results <- lefser(se1, groupCol = "Enviroment")
        lefserPlot(Lefser_results,colors = c("#7fb1d3", "#fdb462"),trim.names=FALSE)
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
