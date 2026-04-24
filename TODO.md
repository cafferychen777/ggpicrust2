# TODO

## First-principles audit fixes

- [x] `ko2kegg_abundance(file=..., data=...)` warning and behavior must agree.
- [x] `ggpicrust2(ko_to_kegg=...)` must normalize logical-like inputs before validation.
- [x] `pathway_volcano()` must handle all-zero p-values without producing infinite plot data.
- [x] `calculate_rank_metric(..., method="t_test")` must tolerate constant features.
- [x] `calculate_rank_metric()` must align samples through the same metadata contract as public entry points.
- [x] `prepare_gene_sets("KEGG")` must exclude BRITE / non-pathway KEGG buckets.
- [x] `compare_gsea_daa()` significant pathway counts must use set semantics.
- [x] `pathway_gsea()` must normalize PICRUSt2 `#NAME` inputs before abundance validation.
- [x] Remove redundant wrapper/core alignment by making alignment ownership explicit.
- [x] `ko2kegg_abundance()` progress output must be opt-in or interactive-only.

## Follow-up core audit fixes

- [x] `import_MicrobiomeAnalyst_daa_results()` documentation and `group_levels` default must agree.
- [x] `import_MicrobiomeAnalyst_daa_results()` must validate/normalize input columns instead of creating `NA` column names.
- [x] `import_MicrobiomeAnalyst_daa_results(file_path=..., data=...)` must follow the package-wide data-wins warning convention.
- [x] `pathway_errorbar(order="group")` must handle tied maximum groups deterministically without replacement warnings.
- [x] `gsea_pathway_annotation()` must preserve GSEA result order after annotation.
- [x] `visualize_gsea(plot_type="network")` must fail fast when required network columns are missing.
- [x] `visualize_gsea()` network parameters must be validated before building delayed ggplot objects.
- [x] `visualize_gsea()` must honor `edge_width_by` or remove the dead API promise.
- [x] `visualize_gsea()` heatmap annotation colors must handle explicit `NULL` consistently.
- [x] `pathway_errorbar()` logical parameters must use the shared logical normalization contract.
- [x] `pathway_errorbar_table()` must remove or honor its unused `ko_to_kegg` parameter.
- [x] `compare_metagenome_results()` must check optional heatmap namespaces before using them.
