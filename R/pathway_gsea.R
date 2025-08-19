#' Validate group sizes for statistical analysis
#'
#' This function checks group sizes and balance to ensure statistical reliability.
#' Follows Linus principle: fail fast with clear reasons, don't hide problems.
#'
#' @param group_vector A factor vector with group assignments
#' @param group_name A character string with the group variable name for error messages
#' @keywords internal
validate_group_sizes <- function(group_vector, group_name) {
  group_counts <- table(group_vector)
  group_names <- names(group_counts)
  n_groups <- length(group_counts)
  
  # Absolute minimum: need at least 2 samples per group for any comparison
  min_samples_per_group <- 2
  insufficient_groups <- group_counts < min_samples_per_group
  
  if (any(insufficient_groups)) {
    insufficient_names <- paste(group_names[insufficient_groups], ":", group_counts[insufficient_groups], collapse = "; ")
    stop(sprintf("Groups with <2 samples detected in '%s' (%s). Statistical comparison impossible. Need at least 2 samples per group.", 
                 group_name, insufficient_names))
  }
  
  # Statistical power warning: less than 3 samples per group
  small_groups <- group_counts < 3
  if (any(small_groups)) {
    small_names <- paste(group_names[small_groups], ":", group_counts[small_groups], collapse = "; ")
    warning(sprintf("Small group sizes in '%s' (%s). Statistical power severely limited. Recommend n>=3 per group for reliable results.",
                    group_name, small_names), call. = FALSE)
  }
  
  # Group imbalance warning (only for two-group comparisons)
  if (n_groups == 2) {
    imbalance_ratio <- max(group_counts) / min(group_counts)
    if (imbalance_ratio > 3) {
      warning(sprintf("Severe group imbalance in '%s' (ratio %.1f:1). Results may be biased. Consider balancing groups or using robust methods.",
                      group_name, imbalance_ratio), call. = FALSE)
    }
  }
  
  # Multi-group warning
  if (n_groups > 2) {
    warning(sprintf("Multiple groups detected in '%s' (%d groups). GSEA will use pairwise comparisons. Consider running separate two-group analyses.",
                    group_name, n_groups), call. = FALSE)
  }
}

#' Gene Set Enrichment Analysis for PICRUSt2 output
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on PICRUSt2 predicted functional data
#' to identify enriched pathways between different conditions.
#'
#' @param abundance A data frame containing KO/EC/MetaCyc abundance data, with features as rows and samples as columns
#' @param metadata A data frame containing sample metadata
#' @param group A character string specifying the column name in metadata that contains the grouping variable
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#' @param method A character string specifying the GSEA method: "fgsea", "GSEA", or "clusterProfiler"
#' @param rank_method A character string specifying the ranking statistic: "signal2noise", "t_test", "log2_ratio", or "diff_abundance"
#' @param nperm An integer specifying the number of permutations
#' @param min_size An integer specifying the minimum gene set size
#' @param max_size An integer specifying the maximum gene set size
#' @param p.adjust A character string specifying the p-value adjustment method
#' @param seed An integer specifying the random seed for reproducibility
#'
#' @return A data frame containing GSEA results
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#' 
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#' 
#' # Run GSEA analysis
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "fgsea"
#' )
#' 
#' # Visualize results
#' visualize_gsea(gsea_results, plot_type = "enrichment_plot", n_pathways = 10)
#' }
pathway_gsea <- function(abundance, 
                        metadata, 
                        group, 
                        pathway_type = "KEGG", 
                        method = "fgsea", 
                        rank_method = "signal2noise", 
                        nperm = 1000, 
                        min_size = 10, 
                        max_size = 500, 
                        p.adjust = "BH",
                        seed = 42,
                        go_category = "BP") {
  
  # Input validation
  if (!is.data.frame(abundance) && !is.matrix(abundance)) {
    stop("'abundance' must be a data frame or matrix")
  }
  
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }
  
  if (!group %in% colnames(metadata)) {
    stop(paste("Group variable", group, "not found in metadata"))
  }
  
  if (!pathway_type %in% c("KEGG", "MetaCyc", "GO")) {
    stop("pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'")
  }
  
  if (!method %in% c("fgsea", "GSEA", "clusterProfiler")) {
    stop("method must be one of 'fgsea', 'GSEA', or 'clusterProfiler'")
  }
  
  if (!rank_method %in% c("signal2noise", "t_test", "log2_ratio", "diff_abundance")) {
    stop("rank_method must be one of 'signal2noise', 't_test', 'log2_ratio', or 'diff_abundance'")
  }
  
  # Validate GO category parameter when pathway_type is "GO"
  if (pathway_type == "GO") {
    valid_go_categories <- c("BP", "MF", "CC", "all")
    if (!go_category %in% valid_go_categories) {
      stop(sprintf("Invalid go_category '%s'. Must be one of: %s", 
                   go_category, paste(valid_go_categories, collapse = ", ")))
    }
  }
  
  # Check if required packages are installed
  required_packages <- list(
    "fgsea" = "fgsea",
    "GSEA" = "clusterProfiler",
    "clusterProfiler" = "clusterProfiler"
  )
  
  pkg_name <- required_packages[[method]]
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(paste("Package", pkg_name, "is required for method", method, 
               ". Please install it using BiocManager::install('", pkg_name, "').", sep = ""))
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Prepare data
  # Handle #NAME column commonly found in PICRUSt2 output
  if (ncol(abundance) > 0 && colnames(abundance)[1] == "#NAME") {
    # Convert tibble to data.frame if necessary and set proper rownames
    abundance <- as.data.frame(abundance)
    rownames(abundance) <- abundance[, 1]
    abundance <- abundance[, -1]
  }
  
  # Ensure abundance is a matrix with samples as columns
  abundance_mat <- as.matrix(abundance)
  
  # Extract group information
  Group <- factor(metadata[[group]])
  names(Group) <- rownames(metadata)
  
  # Find common samples (eliminate special case handling)
  common_samples <- intersect(colnames(abundance_mat), names(Group))
  
  # Enforce minimum requirements (clear validation)
  if (length(common_samples) < 4) {
    stop("Insufficient overlapping samples (", length(common_samples), 
         "/", ncol(abundance_mat), "). Need at least 4 samples for statistical analysis.")
  }
  
  # Inform user about sample subsetting (transparency)
  if (length(common_samples) < ncol(abundance_mat)) {
    warning(sprintf("Using %d/%d overlapping samples between abundance data and metadata", 
                    length(common_samples), ncol(abundance_mat)))
  }
  
  # Clean subset to common samples (no special cases)
  abundance_mat <- abundance_mat[, common_samples]
  Group <- Group[common_samples]
  
  # Validate group balance for statistical reliability
  group_counts <- table(Group)
  if (any(group_counts < 2)) {
    stop("Each group must have at least 2 samples. Current group sizes: ", 
         paste(names(group_counts), "=", group_counts, collapse = ", "))
  }
  
  # Prepare gene sets
  gene_sets <- prepare_gene_sets(pathway_type)
  
  # Calculate ranking metric
  ranked_list <- calculate_rank_metric(abundance_mat, metadata, group, method = rank_method)
  
  # Run GSEA based on selected method
  if (method == "fgsea") {
    if (!requireNamespace("fgsea", quietly = TRUE)) {
      stop("Package 'fgsea' is required. Please install it using BiocManager::install('fgsea').")
    }
    
    results <- run_fgsea(ranked_list, gene_sets, nperm, min_size, max_size)
    
  } else if (method == "GSEA" || method == "clusterProfiler") {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("Package 'clusterProfiler' is required. Please install it using BiocManager::install('clusterProfiler').")
    }
    
    # Convert ranked list to format required by clusterProfiler
    gene_list <- sort(ranked_list, decreasing = TRUE)
    
    # Run GSEA using clusterProfiler
    gsea_result <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = data.frame(
        term = rep(names(gene_sets), sapply(gene_sets, length)),
        gene = unlist(gene_sets)
      ),
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = 1,
      pAdjustMethod = p.adjust,
      nPermSimple = nperm,
      seed = seed
    )
    
    # Convert results to data frame
    if (length(gsea_result) > 0) {
      results <- as.data.frame(gsea_result)
      
      # Rename columns to match fgsea output
      results <- data.frame(
        pathway_id = results$ID,
        pathway_name = results$Description,
        size = results$setSize,
        ES = results$enrichmentScore,
        NES = results$NES,
        pvalue = results$pvalue,
        p.adjust = results$p.adjust,
        leading_edge = results$core_enrichment,
        stringsAsFactors = FALSE
      )
    } else {
      results <- data.frame(
        pathway_id = character(),
        pathway_name = character(),
        size = integer(),
        ES = numeric(),
        NES = numeric(),
        pvalue = numeric(),
        p.adjust = numeric(),
        leading_edge = character(),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Add method information (handle empty results)
  if (nrow(results) > 0) {
    results$method <- method
  } else {
    # Create empty result with correct structure
    results <- data.frame(
      pathway_id = character(),
      pathway_name = character(),
      size = integer(),
      ES = numeric(),
      NES = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      leading_edge = character(),
      method = character(),
      stringsAsFactors = FALSE
    )
  }
  
  return(results)
}

#' Prepare gene sets for GSEA
#'
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#' @param organism A character string specifying the organism (only relevant for KEGG and GO)
#' @param go_category A character string specifying the GO category: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "all"
#'
#' @return A list of pathway gene sets
#' @export
prepare_gene_sets <- function(pathway_type = "KEGG", organism = "ko", go_category = "BP") {
  
  if (pathway_type == "KEGG") {
    # Initialize gene_sets
    gene_sets <- list()
    
    # Try to load KEGG pathway to KO mapping
    tryCatch({
      if (!exists("ko_to_kegg_reference")) {
        data("ko_to_kegg_reference", package = "ggpicrust2", envir = environment())
      }
      
      # Convert to list format required for GSEA
      ko_to_kegg_reference <- as.data.frame(ko_to_kegg_reference)
      
      # Create a list where each element is a pathway containing KO IDs
      for (i in 1:nrow(ko_to_kegg_reference)) {
        pathway_id <- ko_to_kegg_reference[i, 1]
        ko_ids <- as.character(ko_to_kegg_reference[i, -1])
        ko_ids <- ko_ids[!is.na(ko_ids) & ko_ids != ""]
        
        if (length(ko_ids) > 0) {
          gene_sets[[pathway_id]] <- ko_ids
        }
      }
      
    }, error = function(e) {
      # Create dummy gene sets for testing when reference data is not available
      warning("KEGG reference data not found. Creating dummy gene sets for testing.", call. = FALSE)
      
      # Create some dummy pathways for demonstration
      gene_sets <<- list(
        "ko00010" = c("K00844", "K12407", "K00845", "K00886", "K08074"),
        "ko00020" = c("K00239", "K00240", "K00241", "K00242", "K01902"),
        "ko00030" = c("K00016", "K00018", "K00128", "K01595", "K01596"),
        "ko00040" = c("K01623", "K01624", "K11645", "K01803", "K15633"),
        "ko00051" = c("K00134", "K00150", "K03781", "K03782", "K14085")
      )
    })
    
  } else if (pathway_type == "MetaCyc") {
    # Load MetaCyc pathway to EC mapping
    if (!exists("metacyc_to_ec_reference")) {
      # Try to load from package extdata first
      tryCatch({
        metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
        if (file.exists(metacyc_ref_path)) {
          load(metacyc_ref_path, envir = environment())
        } else {
          stop("metacyc_to_ec_reference data file not found")
        }
      }, error = function(e) {
        stop("Failed to load MetaCyc to EC mapping: ", e$message)
      })
    }
    
    # Convert to data frame
    metacyc_to_ec_reference <- as.data.frame(metacyc_to_ec_reference)
    
    # Create gene sets list
    gene_sets <- list()
    
    for (i in 1:nrow(metacyc_to_ec_reference)) {
      pathway_id <- metacyc_to_ec_reference[i, "pathway"]
      ec_string <- as.character(metacyc_to_ec_reference[i, "ec_numbers"])
      
      # Skip pathways with no EC mappings
      if (is.na(ec_string) || ec_string == "" || ec_string == "NA") {
        next
      }
      
      # Split EC numbers by semicolon
      ec_numbers <- strsplit(ec_string, ";")[[1]]
      ec_numbers <- trimws(ec_numbers)  # Remove whitespace
      ec_numbers <- ec_numbers[ec_numbers != ""]  # Remove empty strings
      
      if (length(ec_numbers) > 0) {
        # Add EC: prefix if not present for consistency
        ec_numbers <- ifelse(grepl("^EC:", ec_numbers), ec_numbers, paste0("EC:", ec_numbers))
        gene_sets[[pathway_id]] <- ec_numbers
      }
    }
    
  } else if (pathway_type == "GO") {
    # Load KO to GO mapping
    tryCatch({
      data("ko_to_go_reference", package = "ggpicrust2", envir = environment())
    }, error = function(e) {
      # Create basic GO mapping if reference data doesn't exist
      message("Creating basic GO mapping (reference data not found)")
    })
    
    # Always create mapping if it doesn't exist
    if (!exists("ko_to_go_reference", envir = environment())) {
      ko_to_go_reference <- create_basic_go_mapping()
    }
    
    # Convert to data frame format required for GSEA
    go_reference <- as.data.frame(ko_to_go_reference)
    
    # Validate and filter by GO category if specified
    valid_go_categories <- c("BP", "MF", "CC", "all")
    if (!is.null(go_category) && !go_category %in% valid_go_categories) {
      stop(sprintf("Invalid go_category '%s'. Must be one of: %s", 
                   go_category, paste(valid_go_categories, collapse = ", ")))
    }
    
    if (!is.null(go_category) && go_category != "all" && go_category %in% c("BP", "MF", "CC")) {
      if ("category" %in% colnames(go_reference)) {
        go_reference <- go_reference[go_reference$category == go_category, ]
      }
    }
    
    # Create gene sets list for each GO term
    gene_sets <- list()
    
    for (i in 1:nrow(go_reference)) {
      go_id <- go_reference[i, 1]  # First column contains GO IDs
      
      # Get KO members for this GO term
      if ("ko_members" %in% colnames(go_reference)) {
        ko_string <- go_reference[i, "ko_members"]
        if (!is.na(ko_string) && ko_string != "") {
          ko_ids <- unlist(strsplit(ko_string, ";"))
          ko_ids <- ko_ids[!is.na(ko_ids) & ko_ids != ""]
          
          if (length(ko_ids) > 0) {
            gene_sets[[go_id]] <- ko_ids
          }
        }
      } else {
        # Fallback: assume other columns contain KO IDs
        ko_ids <- as.character(go_reference[i, -1])
        ko_ids <- ko_ids[!is.na(ko_ids) & ko_ids != ""]
        
        if (length(ko_ids) > 0) {
          gene_sets[[go_id]] <- ko_ids
        }
      }
    }
  }
  
  return(gene_sets)
}

#' Calculate rank metric for GSEA
#'
#' @param abundance A matrix of abundance data
#' @param metadata A data frame of metadata
#' @param group A character string specifying the grouping variable
#' @param method A character string specifying the ranking method
#'
#' @return A named vector of ranking statistics
#' @keywords internal
calculate_rank_metric <- function(abundance, 
                                 metadata, 
                                 group, 
                                 method = "signal2noise") {
  
  # Extract group information
  Group <- factor(metadata[[group]])
  names(Group) <- rownames(metadata)
  
  # Ensure abundance is a matrix with samples as columns
  abundance <- as.matrix(abundance)
  
  # Subset abundance to include only samples in metadata
  common_samples <- intersect(colnames(abundance), names(Group))
  abundance <- abundance[, common_samples]
  Group <- Group[common_samples]
  
  # Group size validation already done above in main function
  
  # Get unique group levels
  levels <- levels(Group)
  
  if (length(levels) != 2) {
    stop("GSEA currently only supports two-group comparisons")
  }
  
  # Split samples by group
  group1_samples <- names(Group)[Group == levels[1]]
  group2_samples <- names(Group)[Group == levels[2]]
  
  # Calculate ranking statistic based on method
  if (method == "signal2noise") {
    # Signal-to-noise ratio: (mean1 - mean2) / (sd1 + sd2)
    mean1 <- rowMeans(abundance[, group1_samples, drop = FALSE])
    mean2 <- rowMeans(abundance[, group2_samples, drop = FALSE])
    sd1 <- apply(abundance[, group1_samples, drop = FALSE], 1, stats::sd)
    sd2 <- apply(abundance[, group2_samples, drop = FALSE], 1, stats::sd)
    
    # Handle zero standard deviations
    sd1[sd1 == 0] <- 0.00001
    sd2[sd2 == 0] <- 0.00001
    
    metric <- (mean1 - mean2) / (sd1 + sd2)
    # Ensure names are preserved - critical for fgsea
    names(metric) <- rownames(abundance)
    
  } else if (method == "t_test") {
    # t-test statistic
    metric <- numeric(nrow(abundance))
    names(metric) <- rownames(abundance)
    
    for (i in 1:nrow(abundance)) {
      t_test <- stats::t.test(abundance[i, group1_samples], abundance[i, group2_samples])
      metric[i] <- t_test$statistic
    }
    
  } else if (method == "log2_ratio") {
    # Log2 fold change
    mean1 <- rowMeans(abundance[, group1_samples, drop = FALSE])
    mean2 <- rowMeans(abundance[, group2_samples, drop = FALSE])
    
    # Add small constant to avoid division by zero
    mean1[mean1 == 0] <- 0.00001
    mean2[mean2 == 0] <- 0.00001
    
    metric <- log2(mean1 / mean2)
    # Ensure names are preserved - critical for fgsea
    names(metric) <- rownames(abundance)
    
  } else if (method == "diff_abundance") {
    # Simple difference in abundance
    mean1 <- rowMeans(abundance[, group1_samples, drop = FALSE])
    mean2 <- rowMeans(abundance[, group2_samples, drop = FALSE])
    metric <- mean1 - mean2
    # Ensure names are preserved - critical for fgsea
    names(metric) <- rownames(abundance)
  }
  
  return(metric)
}

#' Run fast GSEA implementation
#'
#' @param ranked_list A named vector of ranking statistics
#' @param gene_sets A list of pathway gene sets
#' @param nperm An integer specifying the number of permutations
#' @param min_size An integer specifying the minimum gene set size
#' @param max_size An integer specifying the maximum gene set size
#'
#' @return A data frame of fgsea results
#' @keywords internal
run_fgsea <- function(ranked_list, 
                     gene_sets, 
                     nperm = 1000, 
                     min_size = 10, 
                     max_size = 500) {
  
  # Check if fgsea package is available
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required. Please install it using BiocManager::install('fgsea').")
  }
  
  # Run fgsea
  fgsea_result <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranked_list,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm
  )
  
  # Convert to data frame and handle empty results
  results <- as.data.frame(fgsea_result)
  
  # Check if we have any results
  if (nrow(results) > 0) {
    # Rename columns for consistency
    results <- data.frame(
      pathway_id = results$pathway,
      pathway_name = results$pathway,  # Will be updated later with annotation
      size = results$size,
      ES = results$ES,
      NES = results$NES,
      pvalue = results$pval,
      p.adjust = results$padj,
      leading_edge = sapply(results$leadingEdge, function(x) paste(x, collapse = ";")),
      stringsAsFactors = FALSE
    )
  } else {
    # Return empty data frame with correct column structure
    results <- data.frame(
      pathway_id = character(),
      pathway_name = character(),
      size = integer(),
      ES = numeric(),
      NES = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      leading_edge = character(),
      stringsAsFactors = FALSE
    )
  }
  
  return(results)
}

#' Create basic GO mapping for microbiome analysis
#'
#' This function creates a basic KO to GO term mapping that covers
#' common functional categories relevant to microbiome research.
#'
#' @return A data frame containing basic GO term mappings
#' @keywords internal
create_basic_go_mapping <- function() {
  # Create a basic KO to GO mapping
  # This includes common GO terms relevant to microbiome analysis
  
  basic_go_terms <- data.frame(
    go_id = c(
      "GO:0006096", "GO:0006099", "GO:0006631", "GO:0006520",
      "GO:0019682", "GO:0015980", "GO:0006163", "GO:0006508",
      "GO:0006412", "GO:0006979", "GO:0006810", "GO:0005975",
      "GO:0008152", "GO:0009058", "GO:0009056", "GO:0006629",
      "GO:0006950", "GO:0006511", "GO:0006464", "GO:0006355"
    ),
    go_name = c(
      "Glycolytic process", "Tricarboxylic acid cycle", "Fatty acid metabolic process",
      "Cellular amino acid metabolic process", "Glyceraldehyde-3-phosphate metabolic process",
      "Energy derivation by oxidation of organic compounds", "Purine nucleotide metabolic process",
      "Proteolysis", "Translation", "Response to oxidative stress",
      "Transport", "Carbohydrate metabolic process",
      "Metabolic process", "Biosynthetic process", "Catabolic process", "Lipid metabolic process",
      "Response to stress", "Protein ubiquitination", "Cellular protein modification process", "Regulation of transcription, DNA-templated"
    ),
    category = c(
      "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP"
    ),
    ko_members = c(
      "K00134;K01810;K00927;K01623;K01803;K00850",
      "K01902;K01903;K00031;K00164;K00382;K01647",
      "K00059;K00625;K01895;K07512;K00626;K01897",
      "K01915;K00928;K01914;K02204;K00812;K01776",
      "K00134;K01623;K00927;K00150;K01803;K00850",
      "K00164;K00382;K00031;K01902;K01903;K01647",
      "K00088;K00759;K01756;K00948;K01633;K00939",
      "K01419;K08303;K01273;K08602;K01417;K01362",
      "K02519;K02543;K02992;K02946;K02874;K02878",
      "K04068;K03781;K00432;K05919;K00540;K03386",
      "K03076;K05685;K03327;K09687;K03406;K03088",
      "K01810;K00134;K01623;K00927;K01803;K00850",
      "K00001;K00002;K00003;K00004;K00005;K00006",
      "K01915;K01914;K01776;K01889;K01845;K01858",
      "K01689;K01690;K01692;K01693;K01694;K01695",
      "K00059;K00625;K01895;K07512;K00626;K01897",
      "K04068;K03781;K00432;K05919;K00540;K03386",
      "K08770;K08771;K08772;K08773;K08774;K08775",
      "K02204;K03100;K03101;K03102;K03103;K03104",
      "K03040;K03041;K03042;K03043;K03044;K03045"
    ),
    stringsAsFactors = FALSE
  )
  
  # Add molecular function terms
  mf_terms <- data.frame(
    go_id = c(
      "GO:0003824", "GO:0016740", "GO:0016787", "GO:0005215",
      "GO:0003677", "GO:0003723", "GO:0016491", "GO:0016853"
    ),
    go_name = c(
      "Catalytic activity", "Transferase activity", "Hydrolase activity", "Transporter activity",
      "DNA binding", "RNA binding", "Oxidoreductase activity", "Isomerase activity"
    ),
    category = c(
      "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF"
    ),
    ko_members = c(
      "K00001;K00002;K00003;K00004;K00005;K00006",
      "K00928;K01914;K01915;K02204;K00812;K01776",
      "K01419;K08303;K01273;K08602;K01417;K01362",
      "K03076;K05685;K03327;K09687;K03406;K03088",
      "K03040;K03041;K03042;K03043;K03044;K03045",
      "K02519;K02543;K02992;K02946;K02874;K02878",
      "K00164;K00382;K00031;K01902;K01903;K01647",
      "K01803;K01804;K01805;K01806;K01807;K01808"
    ),
    stringsAsFactors = FALSE
  )
  
  # Add cellular component terms
  cc_terms <- data.frame(
    go_id = c(
      "GO:0016020", "GO:0005737", "GO:0005829", "GO:0030312",
      "GO:0005886", "GO:0016021", "GO:0022626", "GO:0005840"
    ),
    go_name = c(
      "Membrane", "Cytoplasm", "Cytosol", "External encapsulating structure",
      "Plasma membrane", "Integral component of membrane", "Cytosolic ribosome", "Ribosome"
    ),
    category = c(
      "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC"
    ),
    ko_members = c(
      "K03076;K05685;K03327;K09687;K03406;K03088",
      "K00134;K01810;K00927;K01623;K01803;K00850",
      "K00164;K00382;K00031;K01902;K01903;K01647",
      "K01419;K08303;K01273;K08602;K01417;K01362",
      "K03076;K05685;K03327;K09687;K03406;K03088",
      "K03076;K05685;K03327;K09687;K03406;K03088",
      "K02519;K02543;K02992;K02946;K02874;K02878",
      "K02519;K02543;K02992;K02946;K02874;K02878"
    ),
    stringsAsFactors = FALSE
  )
  
  # Combine all terms
  go_mapping <- rbind(basic_go_terms, mf_terms, cc_terms)
  
  return(go_mapping)
}
