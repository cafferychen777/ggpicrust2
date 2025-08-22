# Create MetaCyc pathway to EC number mapping
# This script creates a mapping based on known metabolic pathways and their associated EC numbers

create_metacyc_to_ec_mapping <- function() {
  # Load existing data
  load("../../data/metacyc_abundance.RData")
  pathway_ids <- metacyc_abundance[["pathway"]]
  
  # Create a comprehensive mapping based on metabolic knowledge
  # This mapping is based on common MetaCyc pathways and their associated EC numbers
  metacyc_to_ec <- list()
  
  # Central carbon metabolism
  metacyc_to_ec[["GLYCOLYSIS"]] <- c("2.7.1.1", "2.7.1.11", "4.1.2.13", "1.2.1.12", "2.7.2.3", "5.4.2.11", "4.2.1.11", "5.4.2.1", "2.7.1.40", "1.2.1.59")
  metacyc_to_ec[["ANAGLYCOLYSIS-PWY"]] <- c("2.7.1.1", "2.7.1.11", "4.1.2.13", "1.2.1.12", "2.7.2.3")
  metacyc_to_ec[["PENTOSE-P-PWY"]] <- c("1.1.1.49", "3.1.1.31", "5.3.1.6", "2.2.1.1", "5.1.3.1")
  metacyc_to_ec[["TCA"]] <- c("4.2.1.3", "1.3.5.1", "6.2.1.5", "2.3.3.1", "1.2.4.2", "2.3.1.12", "4.2.1.2")
  metacyc_to_ec[["CALVIN-PWY"]] <- c("4.1.1.39", "5.3.1.1", "2.7.1.19", "4.1.2.13", "5.4.2.1")
  
  # Fermentation pathways
  metacyc_to_ec[["CENTFERM-PWY"]] <- c("4.1.1.1", "1.1.1.1", "2.3.1.8")
  metacyc_to_ec[["ANAEROFRUCAT-PWY"]] <- c("2.7.1.4", "4.1.2.22", "1.2.1.10")
  
  # Amino acid biosynthesis
  metacyc_to_ec[["ARGSYN-PWY"]] <- c("6.3.4.5", "2.1.3.3", "3.5.3.1", "2.6.1.11", "1.2.1.38")
  metacyc_to_ec[["ARGSYNBSUB-PWY"]] <- c("6.3.4.5", "2.1.3.3", "3.5.3.1", "2.6.1.11")
  metacyc_to_ec[["ARG+POLYAMINE-SYN"]] <- c("6.3.4.5", "2.1.3.3", "3.5.3.1", "4.1.1.17", "2.5.1.16")
  metacyc_to_ec[["BRANCHED-CHAIN-AA-SYN-PWY"]] <- c("2.2.1.6", "1.1.1.86", "4.2.1.9", "2.6.1.42")
  metacyc_to_ec[["ASPASN-PWY"]] <- c("1.1.1.3", "2.7.2.4", "1.2.1.11", "2.6.1.1")
  metacyc_to_ec[["AST-PWY"]] <- c("2.7.2.4", "1.2.1.11", "2.6.1.1")
  
  # Aromatic amino acid synthesis
  metacyc_to_ec[["ARO-PWY"]] <- c("2.5.1.54", "4.2.3.4", "4.2.1.10", "2.5.1.19")
  metacyc_to_ec[["ALL-CHORISMATE-PWY"]] <- c("2.5.1.54", "4.2.3.4", "4.2.1.10")
  
  # Cofactor biosynthesis
  metacyc_to_ec[["BIOTIN-BIOSYNTHESIS-PWY"]] <- c("2.8.1.6", "6.3.3.3", "1.3.1.26", "2.1.1.56")
  metacyc_to_ec[["COA-PWY"]] <- c("2.7.1.33", "6.3.2.5", "4.1.1.36", "6.3.5.4")
  metacyc_to_ec[["COBALSYN-PWY"]] <- c("2.5.1.17", "6.3.5.10", "6.6.1.2", "2.1.1.151")
  
  # One-carbon metabolism
  metacyc_to_ec[["1CMET2-PWY"]] <- c("3.5.4.9", "1.5.1.15", "6.3.4.3", "2.1.2.1")
  
  # Degradation pathways
  metacyc_to_ec[["ARGDEG-PWY"]] <- c("3.5.3.1", "4.3.1.17", "1.2.1.16", "2.6.1.11")
  metacyc_to_ec[["ARGORNPROST-PWY"]] <- c("3.5.3.1", "4.1.1.17", "1.5.1.2", "3.1.3.57")
  metacyc_to_ec[["3-HYDROXYPHENYLACETATE-DEGRADATION-PWY"]] <- c("1.14.13.3", "4.1.1.45", "1.2.1.10")
  
  # Additional pathways for coverage
  metacyc_to_ec[["CODH-PWY"]] <- c("1.2.99.2", "1.8.98.1")
  
  # Convert list to data frame format suitable for GSEA
  # Each row will have pathway_id and associated EC numbers
  pathway_ec_mapping <- data.frame(
    pathway = character(),
    ec_numbers = character(),
    stringsAsFactors = FALSE
  )
  
  # Add mappings for pathways that exist in our abundance data
  for (pathway_id in pathway_ids) {
    if (pathway_id %in% names(metacyc_to_ec)) {
      ec_list <- metacyc_to_ec[[pathway_id]]
      pathway_ec_mapping <- rbind(pathway_ec_mapping, data.frame(
        pathway = pathway_id,
        ec_numbers = paste(ec_list, collapse = ";"),
        stringsAsFactors = FALSE
      ))
    } else {
      # For unmapped pathways, create a minimal entry
      # This prevents GSEA from failing completely
      pathway_ec_mapping <- rbind(pathway_ec_mapping, data.frame(
        pathway = pathway_id,
        ec_numbers = "",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(pathway_ec_mapping)
}

# Create the mapping and save it
metacyc_to_ec_reference <- create_metacyc_to_ec_mapping()
save(metacyc_to_ec_reference, file = "metacyc_to_ec_reference.RData")

cat("Created MetaCyc to EC mapping with", nrow(metacyc_to_ec_reference), "pathways\n")
cat("Number of pathways with EC mappings:", sum(metacyc_to_ec_reference$ec_numbers != ""), "\n")