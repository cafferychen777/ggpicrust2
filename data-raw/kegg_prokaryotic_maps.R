# Build `kegg_prokaryotic_maps` (internal data in R/sysdata.rda): the KEGG
# pathway maps kept by the default prokaryote filter
# (`filter_for_prokaryotes = TRUE`).
#
# Rule. A map is kept if
#   1. it is not an Organismal Systems (09150) or Human Diseases (09160) map,
#      except bacterial infection (09171) and antimicrobial resistance (09175),
#      and
#   2. KEGG draws an organism-specific version of it (e.g. eco02026) for at
#      least one sampled bacterial or archaeal genome.
# KEGG generates an organism-specific map only when the map lies within the
# organism's taxonomic scope and the genome carries its genes. Rule 2 therefore
# removes eukaryotic signalling, cell-cycle, cell-death, cell-junction, glycan
# and plant or animal secondary-metabolism maps that rule 1 leaves in place.
# Some eukaryote-named maps (e.g. ko04142 Lysosome, ko04146 Peroxisome) pass,
# because KEGG draws them for prokaryotes whose genomes carry their enzymes.
#
# Inputs (in data-raw/):
#   kegg_prokaryotic_sample.tsv  The 575 sampled KEGG prokaryotes: the first
#     listed and one random genome of each of the 332 prokaryotic orders in the
#     KEGG Organisms hierarchy (br:br08601), the first genome of 62 gut genera,
#     and 15 reference organisms named by lineage-specific maps (biofilm,
#     Caulobacter cell cycle, bacterial infections, clavulanic acid), which
#     KEGG draws only within the named lineage.
#   kegg_prokaryotic_maps.tsv  One row per map of the bundled KO-to-KEGG
#     reference (KEGG BRITE Release 117.0) and of the KEGG REST reference,
#     with the number of sampled prokaryotes for which KEGG draws the map
#     (`list/pathway/<org>`, retrieved 2026-10-23), counts by sample stratum,
#     six eukaryote controls (hsa, mmu, dme, cel, sce, ath), and the
#     resulting `default_filter_keep` flag.
#
# Usage, from the package root:
#   Rscript data-raw/kegg_prokaryotic_maps.R          # rebuild R/sysdata.rda
#   KEGG_REFRESH=true Rscript data-raw/kegg_prokaryotic_maps.R
# The second form re-queries KEGG REST for the same sampled organisms and
# reports maps whose prokaryote status has changed; it does not overwrite the
# TSV. KEGG asks for modest request rates, so requests are throttled.

maps <- utils::read.delim("data-raw/kegg_prokaryotic_maps.tsv",
                          stringsAsFactors = FALSE)

stopifnot(
  !anyDuplicated(maps$pathway_id),
  identical(maps$default_filter_keep,
            maps$category_filter_keep & maps$n_sampled_prokaryotes >= 1L)
)

if (identical(tolower(Sys.getenv("KEGG_REFRESH")), "true")) {
  sample <- utils::read.delim("data-raw/kegg_prokaryotic_sample.tsv",
                              stringsAsFactors = FALSE)
  drawn <- lapply(sample$org, function(org) {
    Sys.sleep(0.4)
    ids <- sub("\t.*", "", readLines(
      url(paste0("https://rest.kegg.jp/list/pathway/", org)), warn = FALSE))
    unique(paste0("ko", substring(ids, nchar(org) + 1L)))
  })
  counts <- table(unlist(drawn))
  now <- maps$pathway_id %in% names(counts)
  changed <- maps$pathway_id[now != (maps$n_sampled_prokaryotes >= 1L)]
  message(length(changed), " maps changed prokaryote status: ",
          paste(changed, collapse = ", "))
}

kegg_prokaryotic_maps <- maps[maps$default_filter_keep,
                              c("pathway_id", "pathway_name",
                                "n_sampled_prokaryotes")]
rownames(kegg_prokaryotic_maps) <- NULL
attr(kegg_prokaryotic_maps, "sampled_prokaryotes") <- 575L
attr(kegg_prokaryotic_maps, "retrieved") <- as.Date("2026-10-23")

usethis::use_data(kegg_prokaryotic_maps, internal = TRUE, overwrite = TRUE)
