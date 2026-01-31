#!/usr/bin/env Rscript
# ============================================================================
# ggpicrust2 Reference Data Orchestrator
# ============================================================================
#
# Rebuilds all reference datasets used by the package.
# Large datasets delegate to standalone build scripts;
# small datasets (MetaCyc) are built inline.
#
# Usage:
#   source("data-raw/build_all_reference_data.R")
#   build_all()                         # rebuild everything
#   build_all(go = FALSE)               # skip GO (slow, ~35 min)
#   build_all(kegg = TRUE, go = FALSE)  # only rebuild KEGG
#
# Datasets built:
#   ko_to_kegg_reference  <- build_ko_to_kegg_reference.R (BRITE CSV)
#   ko_to_go_reference    <- build_ko_to_go_reference.R   (KEGG + QuickGO)
#   metacyc_reference     <- inline (curated pathway list)
#
# Note: data/ko_reference.rda and data/ec_reference.rda exist but are not
# consumed by any package code. They are legacy files and not rebuilt here.

source("data-raw/build_utils.R")
source("data-raw/build_ko_to_kegg_reference.R")
source("data-raw/build_ko_to_go_reference.R")

# ============================================================================
# Inline builder: MetaCyc pathway reference
# ============================================================================

#' Build metacyc_reference from curated pathway list.
#'
#' MetaCyc requires a registered account for bulk download, so this uses
#' a curated list of major metabolic pathways (based on MetaCyc 29.5).
#'
#' @return data.frame: id, description, category (invisible).
build_metacyc_reference <- function() {
  message("=== Building metacyc_reference ===")
  message("Source: curated MetaCyc 29.5 pathway list\n")

  metacyc_reference <- data.frame(
    id = c(
      # Central carbon metabolism
      "GLYCOLYSIS", "ANAGLYCOLYSIS-PWY", "TCA", "TCA-GLYOX-BYPASS",
      "PENTOSE-P-PWY", "NONOXIPENT-PWY",
      # Fermentation
      "FERMENTATION-PWY", "PWY-5100", "CENTFERM-PWY",
      # Amino acid biosynthesis
      "ARGSYN-PWY", "ASPASN-PWY", "BRANCHED-CHAIN-AA-SYN-PWY",
      "TRPSYN-PWY", "TYRSYN", "PHESYN", "HISTSYN-PWY",
      "LEUSYN-PWY", "ILEUSYN-PWY", "VALSYN-PWY", "METSYN-PWY",
      "GLNSYN-PWY", "GLUTSYN-PWY", "PROSYN-PWY", "SERSYN-PWY",
      "GLYSYN-PWY", "CYSTSYN-PWY", "THRESYN-PWY", "LYSINE-AMINOAD-PWY",
      # Nucleotide biosynthesis
      "PWY-7219", "PWY-7220", "PWY-6123", "PWY-6122",
      "DENOVOPURINE2-PWY", "PWY0-162",
      # Cofactor biosynthesis
      "BIOTIN-BIOSYNTHESIS-PWY", "COA-PWY", "FOLSYN-PWY",
      "NAD-BIOSYNTHESIS-II", "RIBOSYN2-PWY", "THISYNARA-PWY",
      "PYRIDOXSYN-PWY", "PWY-6268", "HEMESYN2-PWY",
      # Fatty acid metabolism
      "FAO-PWY", "FASYN-ELONG-PWY", "PWY-5971",
      # Cell wall biosynthesis
      "PEPTIDOGLYCANSYN-PWY", "PWY-6387",
      # Electron transport
      "PWY-3781", "PWY0-1329", "PWY0-1334"
    ),
    description = c(
      # Central carbon metabolism
      "glycolysis I (from glucose 6-phosphate)",
      "glycolysis III (from glucose)",
      "TCA cycle I (prokaryotic)",
      "superpathway of glyoxylate bypass and TCA",
      "pentose phosphate pathway",
      "pentose phosphate pathway (non-oxidative branch)",
      # Fermentation
      "mixed acid fermentation",
      "pyruvate fermentation to acetate and lactate II",
      "pyruvate fermentation to butanoate",
      # Amino acid biosynthesis
      "L-arginine biosynthesis I (via L-ornithine)",
      "superpathway of L-aspartate and L-asparagine biosynthesis",
      "superpathway of branched amino acid biosynthesis",
      "L-tryptophan biosynthesis",
      "L-tyrosine biosynthesis I",
      "L-phenylalanine biosynthesis I",
      "L-histidine biosynthesis",
      "L-leucine biosynthesis",
      "L-isoleucine biosynthesis I",
      "L-valine biosynthesis",
      "L-methionine biosynthesis I",
      "L-glutamine biosynthesis I",
      "L-glutamate biosynthesis I",
      "L-proline biosynthesis I",
      "L-serine biosynthesis",
      "glycine biosynthesis I",
      "L-cysteine biosynthesis I",
      "L-threonine biosynthesis",
      "L-lysine biosynthesis IV",
      # Nucleotide biosynthesis
      "adenosine ribonucleotides de novo biosynthesis",
      "adenosine deoxyribonucleotides de novo biosynthesis II",
      "inosine-5'-phosphate biosynthesis I",
      "5-aminoimidazole ribonucleotide biosynthesis II",
      "superpathway of purine nucleotides de novo biosynthesis II",
      "superpathway of pyrimidine ribonucleotides de novo biosynthesis",
      # Cofactor biosynthesis
      "biotin biosynthesis I",
      "coenzyme A biosynthesis I",
      "superpathway of tetrahydrofolate biosynthesis",
      "NAD biosynthesis II (from tryptophan)",
      "flavin biosynthesis I (bacteria and plants)",
      "thiamine diphosphate biosynthesis I",
      "pyridoxal 5'-phosphate biosynthesis I",
      "adenosylcobalamin salvage from cobalamin",
      "heme b biosynthesis I (aerobic)",
      # Fatty acid metabolism
      "fatty acid beta-oxidation I",
      "fatty acid elongation -- saturated",
      "palmitate biosynthesis II (bacteria and plants)",
      # Cell wall biosynthesis
      "peptidoglycan biosynthesis I (meso-diaminopimelate containing)",
      "UDP-N-acetylmuramoyl-pentapeptide biosynthesis I",
      # Electron transport
      "aerobic respiration I (cytochrome c)",
      "succinate to cytochrome bd oxidase electron transfer",
      "NADH to cytochrome bd oxidase electron transfer I"
    ),
    category = c(
      rep("Central Carbon Metabolism", 6),
      rep("Fermentation", 3),
      rep("Amino Acid Biosynthesis", 19),
      rep("Nucleotide Biosynthesis", 6),
      rep("Cofactor Biosynthesis", 9),
      rep("Fatty Acid Metabolism", 3),
      rep("Cell Wall Biosynthesis", 2),
      rep("Electron Transport", 3)
    ),
    stringsAsFactors = FALSE
  )

  validate_reference(metacyc_reference,
                     c("id", "description", "category"),
                     "metacyc_reference")
  save_reference(metacyc_reference, "metacyc_reference")

  invisible(metacyc_reference)
}

# ============================================================================
# Orchestrator
# ============================================================================

#' Rebuild all reference datasets.
#'
#' @param kegg Rebuild ko_to_kegg_reference? (default TRUE)
#' @param metacyc Rebuild metacyc_reference? (default TRUE)
#' @param go Rebuild ko_to_go_reference? (default TRUE; slow, ~35 min)
build_all <- function(kegg = TRUE, metacyc = TRUE, go = TRUE) {
  start <- Sys.time()
  message("=== ggpicrust2 Reference Data Build ===")
  message(sprintf("  Started: %s\n", format(start, "%Y-%m-%d %H:%M:%S")))

  results <- list()

  if (kegg) {
    results$ko_to_kegg <- tryCatch(
      build_ko_to_kegg_reference(),
      error = function(e) { message("  KEGG build failed: ", e$message); NULL }
    )
  }

  if (metacyc) {
    results$metacyc <- tryCatch(
      build_metacyc_reference(),
      error = function(e) { message("  MetaCyc build failed: ", e$message); NULL }
    )
  }

  if (go) {
    results$ko_to_go <- tryCatch(
      build_ko_to_go_reference(),
      error = function(e) { message("  GO build failed: ", e$message); NULL }
    )
  }

  # Summary
  elapsed <- difftime(Sys.time(), start, units = "mins")
  message(sprintf("\n=== Build complete (%.1f min) ===", as.numeric(elapsed)))
  for (name in names(results)) {
    status <- if (!is.null(results[[name]])) "OK" else "FAILED"
    message(sprintf("  %s: %s", name, status))
  }

  invisible(results)
}
