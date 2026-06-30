# ggpicrust2 2.5.17

## Bug Fixes

* Count-like controls such as `n_pathways`, `top_n`, and
  `correlation_permutations` now reject values larger than R's integer range
  without emitting coercion warnings. `visualize_gsea()` also applies the same
  no-coercion integer check to GSEA `size` values before plotting.
* `taxa_contribution_bar(show_percentage = TRUE)` now treats absent
  sample/function rows as zero-total combinations when validating percentage
  denominators, instead of silently dropping those sample/function bars from
  the plot. Explicit `function_ids` are also checked after sample alignment so
  partially missing requests fail with the missing IDs instead of being
  silently reduced to the matching subset.
* `pathway_annotation()` now normalizes logical-like `ko_to_kegg` strings,
  preventing `ko_to_kegg = "TRUE"` from silently taking the local-reference
  branch instead of the requested KEGG annotation branch.
* `pathway_gsea(method = "camera"|"fry")` now rejects non-finite or
  rank-deficient design matrices before limma is called. Constant covariates
  or covariates perfectly confounded with the group variable now fail with an
  actionable error instead of entering an unestimable camera/fry contrast.
* `pathway_gsea(method = "camera"|"fry")` now builds limma design formulas
  from literal metadata column names, so non-syntactic group or covariate
  names such as `"treatment group"` and `"age years"` work without requiring
  users to rename their metadata.
* Fixed R scoping bug in `pathway_annotation()`: error counting inside
  `tryCatch()` error handler used `<-` (local assignment) instead of `<<-`,
  so `error_count` and `error_ids` were silently never updated when KEGG
  API calls failed.
* Fixed row-specific KEGG annotation merge-back in `pathway_annotation()`
  when the same feature appears in multiple DAA rows. Non-significant rows
  with the same feature ID as a significant row now keep `NA` annotation
  fields instead of inheriting annotations by feature-name matching.
* `ko2kegg_abundance()` now rejects duplicated KO identifiers after cleaning
  optional `ko:` prefixes. Duplicate KO rows previously entered pathway
  aggregation as repeated evidence and could distort both upper-half mean and
  legacy sum abundances.
* `ko2kegg_abundance()` now rejects missing or non-finite KO abundance values
  before pathway aggregation. The previous upper-half mean path could silently
  drop missing KO values during sorting, changing the number of KOs
  contributing to a pathway abundance.
* `visualize_gsea()` now rejects missing, empty, or duplicate `pathway_id`
  values in the selected rows for every plot type, preventing low-level
  duplicate row-name errors in heatmaps and ambiguous multi-method/contrast
  visualizations. Missing or empty pathway labels now fall back to
  `pathway_id` row-wise.
* `pathway_ridgeplot()` now calculates fold changes using explicit
  group-level semantics. Two-group inputs use factor-level order, while
  multi-group inputs must supply `comparison = c(group1, group2)`, preventing
  sample-column order from silently changing the log2 fold-change direction
  or plotting a comparison that does not match the GSEA contrast.
* Fixed dead auto-adjust logic in `create_legend_theme()`: `missing(direction)`
  was checked after `match.arg(direction, ...)` which always evaluates the

  formal, so the auto-adjust based on legend position never executed.
* Aligned `run_fgsea()` internal `min_size` default (was 10) with the
  public `pathway_gsea()` API default of 5.
* Aligned `perform_aldex2_analysis()` internal `include_effect_size` default
  (was `FALSE`) with the public `pathway_daa()` API default of `TRUE`.
* Two-group ALDEx2 analyses now stop if explicitly requested
  `aldex.effect()` output fails or is malformed, instead of warning and
  returning p-value-only rows without the documented effect-size columns.
  Effect rows are matched to test rows by feature identifier and the required
  finite `effect`, `diff.btw`, `rab.all`, and probability-valued `overlap`
  columns are validated before use.
* Paired ALDEx2 results from `compare_metagenome_results()` now use the same
  feature-ID alignment and effect-output validation. Previously, a reordered
  `aldex.effect()` table could attach one feature's effect size and fold change
  to another feature's p-values by row position.
* Fixed `run_fgsea()` empty-result path: missing method argument to
  `create_empty_gsea_result()` produced `method = "unknown"` instead of
  `"fgsea"`.
* Collapsed dead contrast-selection branches in `run_limma_gsea()` where
  both if/else paths assigned the same value; the unreachable final `else`
  now raises an informative error for unexpected `contrast` types.
* Fixed `pathway_gsea(method = "camera"|"fry")` contrast resolution:
  multi-group designs now require an explicit contrast, character contrasts
  must exactly match a design column or non-reference group level, and
  numeric contrasts are validated against the design matrix. This prevents
  substring matching from silently testing the wrong coefficient.
  Named numeric contrast vectors are now matched to design column names before
  being passed to limma, preventing out-of-order named vectors from silently
  testing the wrong coefficient combination.
* Fixed `pathway_daa(daa_method = "ALDEx2", reference = ...)`: ALDEx2
  condition levels are now releveled before converting them to numeric
  conditions, so `diff_btw`/`log2_fold_change` direction and `group1`/`group2`
  labels honor the requested reference level.
* Fixed `pathway_daa(daa_method = "Lefser", reference = ...)`: Lefser now
  relevels the grouping factor before analysis, converts counts to lefser's
  expected relative-abundance assay with `lefser::relativeAb()`, and reports
  all-feature Kruskal-Wallis p-values from that same relative-abundance
  scale. Previously raw-count library-size differences could drive the
  reported p-values and the user-supplied reference was ignored.
  Actual Lefser backend errors now fail the call instead of returning a table
  with missing LDA scores.
* The top-level `ggpicrust2()` wrapper no longer rejects
  `daa_method = "Lefser"` with the stale claim that Lefser does not output
  p-values. `pathway_daa()` now owns Lefser validation and supplies
  p-values/adjusted p-values, while `pathway_errorbar()` supplies the
  display-only log2 fold-change fallback when no method-native log2 fold
  change is returned.
* `pathway_daa(daa_method = "Lefser")` now fails clearly if a per-feature
  Kruskal-Wallis p-value cannot be computed on the Lefser relative-abundance
  scale, instead of silently converting the failed test to `p = 1`.
* `pathway_daa(daa_method = "LinDA")` now surfaces backend errors as
  errors instead of returning an empty DAA table, so method failures cannot
  be mistaken for "no differential features".
* `pathway_daa()` now validates `reference` once after sample alignment and
  `select` filtering. Invalid, empty, or filtered-out reference levels now
  fail with an actionable error instead of silently falling back to the first
  observed group level.
* `pathway_daa(select = ...)` now requires unique, non-empty sample names and
  rechecks that at least four samples remain after filtering, preventing
  duplicated or undersized selected sample sets from reaching backend fitting.
* `pathway_daa()` now rejects missing/empty group labels and requires at
  least two samples per retained group after sample alignment and `select`
  filtering. This avoids fitting DAA models or tests on singleton groups
  where within-group variation cannot be assessed.
* `pathway_daa()` now rejects missing, non-finite, or negative abundance
  values through shared DAA input validation before backend fitting. Invalid
  abundance cells no longer reach different DAA backends with backend-specific
  low-level failures or partial statistics.
* `pathway_daa()` now preserves method-native adjusted p-values for DESeq2
  (`padj` from `results()`), LinDA (`padj`), and Maaslin2 (`qval`) instead
  of discarding them and recomputing a generic wrapper-level adjustment.
  The requested `p_adjust_method` is forwarded to those backends where
  supported.
* `pathway_daa(daa_method = "DESeq2")` no longer suppresses all backend
  warnings. DESeq2 diagnostics such as identical values across every sample
  now reach the caller instead of being hidden while a result table is
  returned.
* `pathway_daa(daa_method = "DESeq2"|"limma voom"|"edgeR")` now validates
  backend feature identifiers and aligns p-values, adjusted p-values, and
  log fold changes by feature ID before constructing the result table. A
  reordered, incomplete, or malformed backend result now fails clearly
  instead of attaching statistics to features by row position.
* `pathway_daa(daa_method = "metagenomeSeq")` no longer silently replaces
  all CSS normalization-quantile failures with `p = 0.5`. Samples with only
  one positive feature now fail with an actionable error, matching
  metagenomeSeq's own `cumNormStatFast()` requirement, while the existing
  degenerate-search fallback to `p = 0.5` now emits a warning.
* `pathway_daa(daa_method = "Maaslin2")` now mirrors MaAsLin2's
  `make.names()` feature-name sanitization when mapping results back to
  original feature IDs, rejects ambiguous sanitized feature IDs before model
  fitting, and validates MaAsLin2 result columns before returning them.
* `pathway_daa()` now validates raw and adjusted p-value columns at the DAA
  wrapper boundary, including method-specific adjusted p-values supplied by
  individual backends, so malformed probability columns fail instead of
  propagating into downstream summaries.
* `pathway_daa(daa_method = "DESeq2"|"metagenomeSeq"|"LinDA"|"Maaslin2")`
  now models against an internal syntactic metadata column, so valid
  user-supplied group columns with names such as `"treatment group"` work
  without forcing users to rename their metadata. Output `group1` and
  `group2` labels continue to use the original group levels.
* `pathway_daa(daa_method = "LinDA")` now validates the documented
  MicrobiomeStat::linda output columns consumed by ggpicrust2
  (`pvalue`, `padj`, and `log2FoldChange`) and fails on missing, malformed,
  non-finite, or out-of-range values instead of filling them with `NA`.
* Wrapper-computed DAA adjusted p-values are now calculated within each
  method/contrast (`method`, `group1`, `group2`) instead of across all
  pairwise contrasts at once, matching the per-contrast result semantics of
  DESeq2, limma, and edgeR-style outputs.
* Count-based DAA backends that require or assume integer counts (ALDEx2,
  DESeq2, edgeR, and metagenomeSeq) now warn when non-integer abundance
  values are rounded before fitting, replacing previously silent data
  mutation.
* `pathway_daa()` now requires explicit, non-empty, unique feature identifiers
  in abundance row names (or a leading feature-ID column that is normalized to
  row names). Missing IDs fail before backend fitting, and duplicated IDs no
  longer propagate into ambiguous DAA result rows or downstream annotation
  merges.
* `pathway_daa()` now rejects additional arguments supplied through `...`
  instead of silently ignoring them. The previous documentation claimed these
  arguments were passed to backend DAA methods, which could mislead users into
  believing formula, fixed-effect, or covariate-adjustment arguments had been
  applied when the fitted model was actually unchanged.
* `pathway_gsea()` now validates `min_size`, `max_size`, `nperm`, `seed`,
  `p_adjust_method`, and `inter.gene.cor` at the API boundary, including
  rejecting impossible gene-set size ranges before downstream limma/fgsea/
  clusterProfiler calls.
* `pathway_gsea()` now revalidates group structure and complete design
  variables after sample alignment. Missing group labels or camera/fry
  covariates now fail before ranking/model fitting, preventing
  `model.matrix()` from dropping samples and preventing preranked methods
  from silently using a different sample set. This validation is also applied
  when preranked methods use an explicit `comparison = c(group1, group2)`, so
  missing aligned group labels are not silently dropped from the ranked list.
* `pathway_gsea(method = "camera"|"fry")` now passes raw non-negative counts
  directly to `limma::voom()`. The previous pre-voom `+0.5` pseudocount
  changed library sizes and could distort voom's mean-variance trend.
  A failed voom transformation now stops the analysis instead of silently
  switching to a log2 transform with unit weights, which discarded voom's
  observation-level mean-variance weights while still reporting camera/fry
  results.
* Camera/fry GSEA output now explicitly labels its legacy `NES` compatibility
  column as `score_type = "signed_log10_pvalue"` with
  `score_label = "Signed -log10(p-value)"`. limma camera/fry do not estimate
  a true normalized enrichment score, and `visualize_gsea()` now uses the
  explicit score label for axes and legends to avoid misinterpretation.
* `pathway_gsea(method = "fgsea"|"GSEA"|"clusterProfiler")` now accepts
  `comparison = c(group1, group2)` to define preranked GSEA ranking
  direction explicitly. Positive ranking statistics and positive ES/NES now
  map to features higher in `group1`; multi-group preranked analyses must
  specify `comparison` instead of relying on implicit factor-level order.
* `pathway_gsea()` now rejects design-only arguments that cannot be honored
  by preranked GSEA methods. Supplying `covariates` to `fgsea`, `GSEA`, or
  `clusterProfiler`, or supplying `contrast` to a preranked method, now fails
  with instructions to use limma-based `camera`/`fry` or the preranked
  `comparison` argument instead of silently running an unadjusted analysis.
* `pathway_gsea(method = "fgsea")` now honors the requested
  `p_adjust_method` by recomputing adjusted p-values from fgsea raw
  p-values. Previously the fgsea branch always returned fgsea's built-in
  BH-adjusted `padj` values even when callers requested another correction.
* Preranked GSEA ranking vectors are now validated before fgsea or
  clusterProfiler execution. Ranking statistics must be a named numeric vector
  with unique feature names, finite non-missing values, and at least two
  distinct statistics. This prevents all-tied rankings (for example all-zero
  group differences) from producing enrichment results driven by arbitrary
  input order rather than biological signal.
* `pathway_gsea(method = "GSEA"|"clusterProfiler")` now returns a standard
  empty GSEA result schema when no gene sets overlap the ranked feature list
  after size filtering, instead of returning a partial data frame missing
  columns such as `leading_edge` and `method`.
* `pathway_gsea()` now requires a numeric abundance matrix with explicit,
  non-duplicated feature row names (or a leading non-numeric feature-ID
  column) before gene-set matching, so malformed input fails at the API
  boundary instead of producing empty or low-level GSEA errors.
* `pathway_gsea()` and its internal preranked/limma GSEA helpers now reject
  negative, missing, or non-finite count-like abundance values instead of
  replacing them with zero, preventing silent changes to ranking statistics
  and voom's mean-variance model.
* `pathway_gsea(pathway_type = "MetaCyc")` now rejects pathway-level MetaCyc
  identifiers in the abundance row names. MetaCyc GSEA requires EC-level input
  for EC-to-pathway gene sets; pathway-level MetaCyc abundance tables should be
  analyzed with `pathway_daa()` instead of being treated as empty enrichment
  evidence.
* `pathway_gsea()`, `prepare_gene_sets()`, `run_fgsea()`, and
  `run_limma_gsea()` now validate scalar method/pathway/ranking choices and
  named gene-set lists before enrichment testing. Missing, empty, or duplicated
  gene-set names are rejected because those names become `pathway_id` values;
  this prevents downstream R/limma name repair from silently changing pathway
  identifiers such as `set1` into `set1.1`.
* Shared choice-parameter validation now requires a single non-missing
  supported value. This gives clear errors for invalid `plot_type`, `sort_by`,
  `order`, method, and pathway choices instead of low-level R condition-length
  errors when callers pass vectors or `NA`.
* Added shared validation for `p_adjust_method` across DAA/GSEA-facing
  entry points using R's `stats::p.adjust.methods`, so unsupported
  adjustment names fail even when a backend returns method-specific adjusted
  p-values.
* Shared abundance validation now requires numeric matrix input and numeric
  sample columns in data frames, while treating a leading non-numeric feature
  ID column as metadata rather than a sample. This prevents malformed
  abundance tables from reaching low-level `colSums()`, `rowMeans()`,
  `scale()`, or backend model-fitting calls with misleading errors.
* Abundance-consuming entry points now normalize a leading non-numeric feature
  ID column into row names before sample alignment and matrix conversion. This
  prevents feature/pathway IDs from being dropped or replaced by default row
  numbers in `pathway_daa()`, `pathway_gsea()`, heatmaps, PCA, ridge plots, and
  error-bar summaries.
* `ggpicrust2(data = ..., ko_to_kegg = FALSE)` now preserves abundance inputs
  that already store feature IDs in row names. The wrapper no longer
  unconditionally treats the first sample column as an ID column, preventing
  the first sample from being dropped and abundance values from becoming
  feature labels.
* `pathway_heatmap()` now handles zero-variance rows or sample profiles when
  row/column clustering uses correlation or Spearman distance, instead of
  sending `NA` distances to `hclust()`.
* `pathway_heatmap()` now accepts finite transformed inputs with negative
  values or zero-sum sample columns, while explicitly rejecting missing or
  non-finite values before z-score scaling.
* `pathway_heatmap()` now revalidates primary and secondary grouping variables
  after sample alignment and rejects missing aligned group labels, preventing
  extra metadata rows from making an otherwise invalid heatmap grouping appear
  valid.
* `ko2kegg_abundance(method = "abundance")` now matches PICRUSt2's
  unstructured pathway upper-half indexing for even numbers of matched KOs
  (`floor(n / 2) + 1` in R). The previous `ceiling(n / 2)` start included
  one lower-half KO for even pathway sizes and could underestimate pathway
  abundance.
* Clarified `ko2kegg_abundance()` documentation: the default method is an
  offline KO-to-KEGG aggregation approximation based on PICRUSt2's
  unstructured pathway rule, not the full PICRUSt2 pathway pipeline with
  MinPath and structured MetaCyc pathway inference.
* Fixed `pathway_errorbar()` documentation defaults for `pvalue_format`
  (was "smart", code uses "numeric") and `pathway_class_position` (was
  "left", code uses "right").
* Fixed significance star/color assignment in `get_significance_stars()` and
  `get_significance_colors()`: iteration order now goes from least to most
  significant threshold so the tightest match wins.
* `pathway_pca()` now correctly requires at least 2 groups (`min_groups = 2`)
  instead of allowing single-group input that produces a degenerate PCA.
* `pathway_pca()` now rejects missing or empty group labels after sample
  alignment. Previously a metadata table could retain ungrouped samples in
  the PCA coordinates as long as the non-missing labels still contained two
  groups, making color and marginal-density interpretation incomplete.
* `pathway_pca()` now draws confidence ellipses only for groups with at least
  four samples, matching ggplot2's `stat_ellipse()` implementation. Smaller
  groups remain in the PCA scatter plot and trigger an
  explicit warning instead of relying on ggplot2's low-level
  "Too few points to calculate an ellipse" message.
* `pathway_pca()` no longer rejects finite zero-sum sample columns and no
  longer drops samples whose values are constant across pathways. In
  `prcomp(t(abundance))`, samples are observations and pathways are variables;
  only zero-variance pathways are removed before scaling.
* `pathway_ridgeplot()` now guards against `NA` values in direction and NES
  columns instead of silently propagating them into ggplot aesthetics.
* `pathway_ridgeplot()` now aligns metadata to abundance columns by sample ID
  before calculating displayed gene/KO log2 fold changes, validates GSEA
  p-value/FDR/NES columns and display parameters, and ranks `sort_by = "NES"`
  by absolute effect size.
* `pathway_ridgeplot()` now rejects missing, empty, or duplicated
  `pathway_id` values and invalid `pathway_type` values before reference
  mapping. Missing or empty `pathway_name` labels now fall back to
  `pathway_id` instead of causing low-level string-length errors.
* `visualize_gsea()` now validates required GSEA statistic columns before
  plotting, including finite NES values, probability-valued p-values/FDRs,
  and positive integer gene set sizes.
* `gsea_pathway_annotation()` now rejects missing or empty `pathway_id`
  values and invalid `pathway_type` values at the API boundary, preventing
  annotated GSEA tables with blank or `NA` pathway labels.
* `visualize_gsea(plot_type = "network"|"heatmap")` now parses missing or
  empty `leading_edge` values as empty sets, preventing missing leading-edge
  annotations from creating false pathway similarity edges.
* `visualize_gsea(plot_type = "heatmap")` now recognizes abundance input with
  a leading non-numeric feature-ID column, validates the aligned abundance
  matrix, rejects missing group annotations, and fails clearly when non-empty
  leading-edge genes do not match abundance row names instead of drawing an
  all-zero heatmap.
* `import_MicrobiomeAnalyst_daa_results()` now parses MicrobiomeAnalyst result
  columns by semantic names such as `Pvalues`, `FDR`, `Statistics`, and `log2FC`
  instead of assuming a fixed four-column order. Feature IDs may come from a
  feature/name column or from row names, preventing method-specific DE outputs
  from being silently misread as p-values/FDR values. Imported feature IDs,
  p-values, FDR values, optional statistics/fold changes, method labels, and
  group labels are validated before returning a ggpicrust2-style DAA table.
* P-value annotation helpers now validate probability values, significance
  thresholds, star symbols, and colors. `pathway_errorbar()` validates its
  p-value display parameters at the API boundary.
* `pathway_volcano()` now validates that `fc_threshold` is a non-negative
  number.
* `pathway_volcano()` and `compare_gsea_daa()` now validate probability-valued
  p-value/FDR columns and p-value thresholds, and keep zero adjusted p-values
  finite in scatter/volcano log-scale displays.
* `compare_gsea_daa(plot_type = "scatter")` now rejects duplicated pathway
  effect-size rows before merging GSEA and DAA results, preventing many-to-many
  joins from inflating plotted associations. Pathway identifiers must also be
  non-empty.
* `compare_gsea_daa(plot_type = "scatter")` now requires explicit GSEA and
  DAA group-direction columns and aligns DAA `log2_fold_change` values to the
  GSEA-positive NES direction before plotting. This prevents default
  two-group analyses from showing an artificial sign reversal between
  preranked GSEA (`group1` vs `group2`) and DAA (`group2/group1`) effects.
* `compare_gsea_daa(plot_type = "venn"|"upset")` now rejects direction-aware
  inputs that contain multiple or incompatible `group1`/`group2` pairs,
  preventing pathways from different biological contrasts from being counted
  as method agreement.
* `compare_daa_results()` now compares multi-group DAA discoveries as
  feature/group-pair units rather than feature IDs alone, preventing agreement
  from being overstated when different methods flag the same feature in
  different pairwise contrasts. Group pairs are canonicalized as unordered for
  this set-level comparison, so `A vs B` and `B vs A` are treated as the same
  biological comparison while `A vs B` and `A vs C` remain distinct. Method
  labels and feature/group identifiers are validated for unambiguous output.
* `compare_metagenome_results()` now drops undefined per-feature Spearman
  correlations before summarizing and fails with a clear error when all shared
  features are constant across aligned samples, instead of surfacing a
  low-level `wilcox.test()` error.
* `compare_metagenome_results()` now preserves the paired-sample design during
  differential abundance analysis. The function uses paired ALDEx2 tests by
  default and also supports paired Wilcoxon signed-rank tests on relative
  abundance. Independent-group backends are rejected because treating repeated
  quantifications of the same samples as independent replicates produces
  pseudoreplication.
* Correlation p-values in `compare_metagenome_results()` now use joint
  sample-label permutations of the median per-feature Spearman correlation,
  preserving dependence among features while breaking cross-metagenome sample
  correspondence. Diagonal p-values are `NA`, Monte Carlo p-values use the
  plus-one correction, and a multiplicity-adjusted `p_adjust_matrix` plus the
  contributing-feature counts are returned.
* `compare_metagenome_results()` now validates metagenome labels, feature row
  names, sample column names, and finite non-negative matrix values before
  intersecting matrices. Duplicate feature/sample IDs now fail instead of
  creating ambiguous by-name alignments.
* `compare_daa_results()`, `pathway_errorbar()`, `pathway_errorbar_table()`,
  `pathway_annotation()`, and the top-level `ggpicrust2()` wrapper now use
  shared validation for adjusted p-value columns and significance thresholds
  before filtering significant features.
* `pathway_errorbar_table()` now aligns named `Group` vectors to abundance
  sample columns when metadata is not provided, matching `pathway_errorbar()`
  behavior and preventing order-dependent mean/log2 fold-change errors.
  Both functions now reject duplicated sample names in `Group`.
* `pathway_errorbar()` and `pathway_errorbar_table()` now reject duplicated
  feature IDs within a selected DAA method/group pair before merging abundance
  statistics back to annotations, preventing duplicated plot/table rows from
  many-to-one DAA result inputs.
* `pathway_errorbar()` and `pathway_errorbar_table()` now reject missing,
  empty, or incompatible `Group` labels before calculating group means and
  standard deviations. Previously samples with `NA` group labels could be
  silently dropped from plotted/table abundance summaries, and mismatched
  group labels could produce summaries unrelated to the DAA contrast.
* `read_contrib_file()`, `read_strat_file()`, and
  `aggregate_taxa_contributions()` now validate contribution identifier keys
  (`sample`, `function_id`, and `taxon`) before aggregation. Missing or empty
  IDs, duplicated wide-format stratified sample columns, and invalid requested
  pathway/function IDs now fail fast instead of being silently dropped or
  mislabeled by downstream aggregation.
* `read_contrib_file()` and `aggregate_taxa_contributions()` now reject
  contribution tables that mix gene-family-level identifiers (such as KOs/ECs)
  with pathway-level identifiers in the same `function_id` column. This avoids
  silently treating direct pathway contributions and KEGG pathway-to-KO
  expansions as one comparable biological unit.
* `aggregate_taxa_contributions()` now validates DAA adjusted p-values,
  significance thresholds, `top_n`, and contribution metric values. Missing,
  infinite, or negative contribution values now fail fast instead of being
  silently dropped or propagated into taxa contribution summaries.
* `aggregate_taxa_contributions()` now ranks `top_n` taxa by total
  contribution mass instead of mean contribution over observed rows, avoiding
  over-selection of sparsely observed taxa with a single high contribution.
* `taxa_contribution_bar()` and `taxa_contribution_heatmap()` now validate
  contribution values, identifier keys, requested `function_ids`, and
  `n_functions` before plotting.
* `taxa_contribution_bar(show_percentage = TRUE)` now rejects zero-total
  sample/function combinations instead of displaying undefined relative
  contributions as all-zero percentage bars.
* `taxa_contribution_heatmap()` now treats absent sparse contribution
  combinations as zero when computing sample means, matching PICRUSt2's
  sparse contribution outputs and avoiding inflated mean heatmap intensities
  for taxa observed in only a subset of samples.
* `taxa_contribution_heatmap(annotation_data = ...)` now rejects conflicting
  labels for the same plotted function ID instead of silently choosing the
  first duplicate annotation label.
* `taxa_contribution_bar()` now ranks default `n_functions` by between-sample
  variance in total function contribution, after summing taxa and filling
  absent sparse sample/function combinations with zero. This prevents taxon
  composition heterogeneity within a function from being mistaken for
  between-sample functional variation.
* Corrected the multi-group ALDEx2 method label from
  `ALDEx2_Kruskal-Wallace test` to `ALDEx2_Kruskal-Wallis test`; package
  internals still recognize the legacy spelling as an alias for
  backward compatibility.

## Internal

* Renamed `order` variable to `sort_idx` in `pathway_errorbar()` to avoid
  shadowing the function parameter; removed unreachable default `switch()`
  branch.
* Removed redundant pre-sort in `visualize_gsea()` enrichment plot
  (`reorder()` in the aesthetic handles display ordering).
* Trimmed verbose per-sample/per-row log messages in `pathway_heatmap()`.

# ggpicrust2 2.5.16

## New Features

* Added pathway-level taxa contribution support for PICRUSt2
  `path_abun_contrib.tsv` output via `read_pathway_contrib_file()` and
  expanded `read_contrib_file(type = "auto")` parsing.
* `aggregate_taxa_contributions()` now supports contribution tables that do
  not contain `norm_taxon_function_contrib`, using available PICRUSt2
  contribution metrics without requiring a differential abundance step.
* `pathway_annotation()` now accepts data-frame input, including
  rowname-based `ko2kegg_abundance()` output, and can annotate local KEGG
  pathway IDs with `pathway = "KEGG"`.

# ggpicrust2 2.5.14

## Behavior Changes

* `pathway_daa()` now defaults `include_effect_size = TRUE`. ALDEx2 results
  include `effect_size`, `diff_btw`, `log2_fold_change`, `rab_all`, and
  `overlap` columns by default, aligning ALDEx2 output with DESeq2, edgeR,
  limma voom, LinDA, Maaslin2, and metagenomeSeq, which all return log2 fold
  changes without an opt-in flag. The extra `ALDEx2::aldex.effect()` call
  reuses the already-computed CLR object, so the cost is modest. Pass
  `include_effect_size = FALSE` to restore the p-value-only output
  (requested in #181).
* Multi-group ALDEx2 runs no longer warn that "effect size only available
  for two-group comparisons"; the flag is silently ignored when
  `aldex.effect()` cannot apply.
* `pathway_daa(daa_method = "DESeq2", reference = ...)` now honors the
  user-supplied reference level and supports multi-group designs (one
  row-block per non-reference contrast), matching the shape returned by
  `edgeR` and `limma voom`. Previously the `reference` argument was
  silently ignored and multi-group input yielded only a single
  `Level[2]` vs `Level[1]` contrast.

## Bug Fixes

* `ko2kegg_abundance(filter_for_prokaryotes = TRUE)` now actually removes
  eukaryotic / human-system pathway classes from the bundled KEGG hierarchy.
  The previous filter matched obsolete Level 2 labels without the `091xx`
  BRITE prefixes used in `ko_to_kegg_reference`, so it removed no rows.
  `ko2kegg_abundance()` also now excludes KEGG BRITE hierarchies and
  "Not Included in Pathway or Brite" buckets such as `ko99980` before
  abundance calculation, because they are not pathway maps and cannot be
  consistently annotated as pathways. Bacterial infection and antimicrobial
  resistance pathways remain available in the default prokaryotic mode.
* `pathway_daa()` with `include_abundance_stats = TRUE` no longer produces
  `log2_fold_change.x` / `log2_fold_change.y` columns when the DAA method
  already returns its own `log2_fold_change` (ALDEx2 with effect size,
  DESeq2, edgeR, limma voom, LinDA, Maaslin2, metagenomeSeq). The
  method-native log2 fold change is kept and the relative-abundance ratio
  is not recomputed, so model-based and ratio-based effect sizes are never
  conflated under the same column name.
* `pathway_daa(include_abundance_stats = TRUE)` now fails if the requested
  abundance summaries cannot be calculated for the returned feature/group
  pairs, instead of warning and returning a result table that omits the
  user-requested summary columns.
* `pathway_daa(select = ...)` now reorders metadata rows to match the
  reordered abundance columns. Previously the `select` branch reordered
  abundance but only filtered metadata with `%in%`, leaving group labels
  desynchronized from samples whenever `select` was not in natural order,
  which silently produced wrong p-values and log fold changes.
* `pathway_daa(daa_method = "limma voom")` with three or more groups now
  labels `group2` correctly. Previously the length-(k-1) contrast vector
  was recycled into a result of `n_features * (k-1)` rows, producing
  interleaved `B,C,B,C,...` labels that no longer corresponded to the
  `as.vector()`-flattened p-values and coefficients.
* `pathway_daa(daa_method = "Maaslin2")` with three or more groups now
  emits one row per (feature, non-reference level) contrast instead of
  flattening the output via `match()` to a single row per feature. The
  `group2` column is taken from Maaslin2's own `value` column rather than
  being recycled from a length-(k-1) vector.
* `pathway_daa(daa_method = "metagenomeSeq")` no longer hardcodes
  `metadata$sample`; sample-column autodetection from `align_samples()`
  is honored so metadata with a non-default sample identifier (e.g.
  `SampleID`) works end-to-end.
* `pathway_daa(daa_method = "metagenomeSeq")` now extracts p-values and
  log-fold changes from the complete `fitFeatureModel()` result using feature
  identifiers. The previous code called `MRcoefs()`, a top-table display
  helper that can return sorted or partial feature tables and currently fails
  on `fitFeatureModelResults` in supported metagenomeSeq versions; ggpicrust2
  then returned all-`NA` `log2_fold_change` values with only a warning.
* `pathway_daa(daa_method = "metagenomeSeq")` with three or more groups
  now emits one row-block per (reference, non-reference) contrast --
  `(k - 1) * n_features` rows total -- matching the shape returned by
  DESeq2 / edgeR / limma voom / LinDA / Maaslin2. Previously the
  function built a full k-column model matrix, called
  `fitFeatureModel()` once, read `coef = 2`, and hard-coded the labels
  as `group1 = Level[1] / group2 = Level[2]`, so any contrast beyond
  the first non-reference level was silently dropped while the output
  shape looked like a two-group result. `fitFeatureModel()` is also
  metagenomeSeq's documented two-group entry point (it tests a single
  coefficient and returns one p-value per feature), so each pairwise
  contrast is now refit on the subset of samples in the two levels of
  interest. The `reference` argument governs which level is held fixed
  as `group1` across all contrasts.
* Calling `pathway_daa(daa_method = "Maaslin2")` twice in the same R
  session no longer fails with `cannot open the connection`. Stale
  handlers left by Maaslin2's `logging` package after its first-call
  tempdir is cleaned up are now cleared before each invocation.
* `pathway_daa()` now rejects abundance matrices containing negative
  values or duplicate sample identifiers at the validation layer,
  instead of letting them propagate into method-specific failures with
  cryptic messages.
* `pathway_daa()` and `pathway_errorbar()` now refuse sample columns
  with a total abundance of zero (or NA) instead of silently producing
  NaN inside the `x / sum(x)` relative-abundance step. The NaN used to
  be absorbed by downstream `mean(..., na.rm = TRUE)` aggregations, so
  group statistics and error-bar plots were computed from fewer samples
  than supplied with no warning. The error now names the offending
  sample(s). The shared `compute_relative_abundance()` helper replaces
  the duplicated `apply(., 2, function(x) x / sum(x))` idiom, and
  `validate_abundance()` gained a `check_zero_columns` gate so every
  entry point shares the same contract.
* `pathway_daa()` now validates `daa_method` against the supported set
  up front and suggests the canonical spelling for common typos
  (e.g. `"linDA"` -> `"LinDA"`, `"Lefse"` -> `"Lefser"`, `"aldex"` ->
  `"ALDEx2"`). Previously the method dispatch fell through `switch()`
  with no default branch and returned `NULL` silently, breaking
  downstream annotation and plotting with opaque errors.
* `compare_metagenome_results()` now aligns every input metagenome on
  the shared feature set by row name before the `cbind` and per-feature
  Spearman correlation steps. Previously both steps indexed rows by
  position, so two metagenomes with identical row names but different
  row orders were compared feature-by-position and produced meaningless
  (often negative) correlations. A by-name alignment now correctly
  returns correlation = 1 for identical inputs regardless of row order.
  The function also errors cleanly when a metagenome is missing row
  names or when the metagenomes share no features.
* `compare_metagenome_results()` now aligns every input metagenome on
  the shared **sample** set by column name in addition to the existing
  feature alignment. Previously the per-feature Spearman correlation
  computed `stats::cor(m1[k, ], m2[k, ])` by indexing columns by
  position, so two metagenomes with identical content but columns in
  different orders were compared "sample i of metagenome A" against
  "sample i of metagenome B" as if they were the same biological
  sample, producing median correlations that could be strongly negative
  (e.g. `-0.6`) on data that was really identical under by-name
  alignment. Mismatched column counts also used to fall through to
  `stats::cor()` and abort mid-loop with "incompatible dimensions"; the
  function now stops at the boundary with an actionable message when
  column names are missing, sample sets are disjoint, and warns when
  the intersection drops samples. Per-sample cross-metagenome
  correlation is only defined for parallel samples; this makes that
  contract explicit instead of enforcing it by happy-path coincidence.
* `find_sample_column()` now requires the Priority 1 standard-named
  column (e.g. `sample`, `Sample`, `sample_id`, `sample_name`) to
  contain unique values that match the abundance sample IDs. Previously
  a standard-named column with duplicate values (e.g. `sample =
  c("S1","S1","S2","S2")`) was picked purely on its name, silently
  misaligning every downstream function that relies on
  `align_samples()`.
* `ggpicrust2()`'s internal call to `pathway_daa()` now uses the modern
  `p_adjust_method` argument. Previously it still forwarded the
  deprecated `p.adjust` argument, so every normal `ggpicrust2()` call
  emitted the deprecation warning that is meant to fire only when a
  user explicitly supplies the legacy name. The legacy `p.adjust`
  parameter remains accepted with a deprecation warning for backward
  compatibility.
* `pathway_errorbar_table()` now delegates sample-column detection and
  Group reordering to the shared `align_samples()` helper, removing a
  parallel implementation of the same logic as well as a duplicated
  `length(Group) != ncol(abundance)` check. Accepted metadata shapes
  are now identical to `pathway_daa()` and `ggpicrust2()`.
* `DESCRIPTION` now lists `Maaslin2` and `metagenomeSeq` under
  `Suggests`. Both packages are documented as supported
  \code{daa_method} values in \code{pathway_daa()} and are dispatched
  via \code{require_package()}, but neither appeared in `Suggests`, so
  the declared dependency graph, the public method list, and the
  runtime dispatch had drifted apart. Declaring them brings the three
  sources back into alignment.
* `pathway_gsea(organism = ...)` and
  `prepare_gene_sets(organism = ...)` now warn when a non-default value
  is supplied. Both functions advertised an `organism` argument but
  the KEGG and GO branches read KO-based reference tables
  (`ko_to_kegg_reference`, `ko_to_go_reference`) that are
  organism-independent by construction, so a caller passing e.g.
  `organism = "hsa"` silently got the same gene sets as `"ko"`. The
  argument is retained for signature compatibility with a deprecation
  warning, its documentation now records the no-op, and the parameter
  will be removed in a future release. Callers that rely on the
  default value are unaffected.
* `ggpicrust2()` no longer forwards its `select` argument into
  `pathway_daa()`. The wrapper's `@param select` documents a vector of
  **pathway names** for plot-time feature selection, but the same value
  was also being passed into `pathway_daa()` where `select` means
  **sample names**, so any call with `select = <pathway names>` aborted
  immediately inside `pathway_daa()` with
  "Some selected samples not in abundance data". The user-supplied
  `select` is now only forwarded to `pathway_errorbar()` -- where
  feature-level filtering for the figure actually happens -- and
  `pathway_daa()` runs on the full sample set as documented.
* `pathway_errorbar()` no longer silently overwrites a method-native
  `log2_fold_change` column with a relative-abundance mean ratio. The
  previous code added the column as NA only when missing, then
  unconditionally overwrote every row inside a `for` loop, so the
  effect size displayed in the side panel disagreed with the model
  output that produced the p_adjust shown next to it. The bar is now
  taken as-is when the DAA method supplies `log2_fold_change` (DESeq2,
  edgeR, limma voom, LinDA, Maaslin2, metagenomeSeq, and ALDEx2 with
  `include_effect_size = TRUE`), so both panels of the same figure
  report the same model-based estimate. The mean-ratio fallback still
  runs when no `log2_fold_change` column is supplied (ALDEx2 with
  `include_effect_size = FALSE`, Lefser, or custom DAA frames).
* `pathway_errorbar()` and `pathway_errorbar_table()` now require every
  displayed significant DAA feature to be present in the abundance row names.
  The previous intersection-based subsetting could silently drop significant
  DAA rows from abundance summaries or leave plot panels describing different
  feature sets. `pathway_errorbar()` also validates displayed method-native
  `log2_fold_change` values, so non-finite effect sizes fail before plotting.
* `pathway_errorbar_table()` no longer derives its two group names via
  `unique(daa_results_filtered_sub_df$group1)[1]` /
  `unique(daa_results_filtered_sub_df$group2)[1]`. The surrounding
  `validate_daa_results()` call already hard-rejects multi-contrast
  input, so the `unique(...)[1]` idiom was dead defensive code --
  but it was also shaped exactly like "silently pick the first of
  many", which would have masked any future validator bypass by
  collapsing a `(k-1) * n_features` multi-contrast DAA result (as
  produced by `pathway_daa()` for >=3 groups with DESeq2 / edgeR /
  limma voom / LinDA / Maaslin2 / metagenomeSeq) down to a single
  contrast without warning. Both names are now read via direct `[1]`
  indexing, keeping the fast-fail path routed through the validator
  where contract violations belong.
* `pathway_errorbar()` and `pathway_errorbar_table()` (via
  `calculate_abundance_stats()`) now derive per-feature, per-group mean
  and standard deviation from a single shared helper,
  `summarize_abundance_by_group()`. Previously `pathway_errorbar()`
  rolled its own `pivot_longer() %>% group_by(name, group) %>%
  summarise(mean(value), sd(value))` path without `na.rm = TRUE`, so the
  same abundance matrix could produce different bar heights in the plot
  versus the companion table if any NA slipped through the pipeline.
  Unifying the aggregation removes that latent divergence and guarantees
  both entry points evolve together.

## Internal

* Documented the intentional duplicate `align_samples()` call in
  `ggpicrust2()` and `pathway_daa()`. The wrapper must pre-align
  abundance/metadata before Step 4 builds `Group_vec` by positional
  zipping of `metadata[[group]]` with `colnames(abundance)`;
  `pathway_daa()` must also align independently to honor its
  standalone-caller contract. The two call sites are deliberately
  invoked with identical arguments, and `align_samples()` is
  deterministic and idempotent, so running it twice on the same
  inputs is a cheap no-op and cannot drift. A regression test in
  `test-data_utils.R` now locks the idempotency invariant so any
  future change that breaks it fails loudly at test time.
* Removed the `"nonsense"` placeholder columns and values from
  `pathway_errorbar()`'s internal data frames. The log2-fold-change bar
  now sets its fill directly on `geom_bar()` instead of routing a single
  color through `aes(fill = group_nonsense)` + `scale_fill_manual()`,
  and the pathway-class / p-value side panels anchor all labels at a
  single x via `aes(x = "")` instead of padding each data frame with a
  constant dummy column. A dead `$group2 <- "nonsense"` column that was
  written but never read downstream is also gone. No user-visible
  change.
* `pathway_annotation(file = ..., ko_to_kegg = FALSE)` now actually
  populates the `description` column. The file-mode branch previously
  extracted features from sample column names (skipping columns 1 and 2
  and taking the rest as IDs), so the description column was always
  filled with `NA`. Feature IDs are now read from the first column, in
  one unified code path shared with the DAA-results branch.
* `pathway_daa(daa_method = "metagenomeSeq")` no longer aborts with
  `missing value where TRUE/FALSE needed` on small or near-uniform
  inputs (e.g. the minimum 4-sample / 2-group case). The normalization
  quantile returned by `cumNormStatFast()` is now checked for NA/NaN
  and falls back to metagenomeSeq's documented default (`p = 0.5`).
* `pathway_errorbar_table(sample_col = ...)` now defaults to
  auto-detection via the same logic `align_samples()` uses, so
  metadata using `sample`, `Sample`, `sample_id`, etc. works without
  passing `sample_col` explicitly. Previously the default was
  hardcoded to `"sample_name"`, which caused the common `sample`
  convention to fail with `Column 'sample_name' not found in metadata`.
* `pathway_daa(daa_method = "LinDA", reference = ...)` now actually
  contrasts against the user-specified reference level. The formula
  passed to `MicrobiomeStat::linda()` did not relevel the grouping
  factor, so LinDA kept the factor's natural first level as reference
  while the `group1` label in the result used the user's `reference`
  argument. That produced rows with `group1 == group2` and an
  `log2_fold_change` whose sign did not reflect the requested
  direction.
* `pathway_daa(daa_method = "Maaslin2", reference = ...)` now honors
  the requested reference in the two-group case. Previously the
  `reference` argument was only forwarded to Maaslin2 for k > 2
  groups and the two-group branch passed `NULL`, so Maaslin2 fell
  back to its alphabetical default while the result labeled
  `group1 = <user reference>` -- producing the same `group1 == group2`
  /no-sign-flip symptom as the LinDA bug above.
* `pathway_daa()` now re-validates the group count after sample
  alignment and `select` filtering. A narrow `select =` that removes
  every sample of a level, or `align_samples()` dropping the only
  samples for a group, would previously let a single-group dataset
  reach the backends with a less actionable downstream error.
* `pathway_annotation(ko_to_kegg = TRUE)` now returns the full input
  `daa_results_df` with annotation columns, populating `pathway_name`,
  `pathway_description`, `pathway_class`, and `pathway_map` only for
  rows where `p_adjust < p_adjust_threshold` and leaving the rest as
  `NA`. Previously it returned just the significant subset, silently
  dropping non-significant rows that downstream code (e.g.
  `ggpicrust2()`'s `plot_result_list$daa_results_df`) expected to
  still be present.
* `pathway_annotation(ko_to_kegg = TRUE, organism = ...)` no longer
  picks the first organism-specific gene linked to a KO as a
  "representative" and fetches its record in place of the KO entry.
  Gene order from `KEGGREST::keggLink()` is not semantically
  meaningful, so isozymes or paralogs participating in different
  pathways could yield different annotations across KEGG builds. The
  function now always fetches the generic KO entry (the authoritative
  KO-level record) and rewrites pathway IDs from the `ko` prefix to
  the organism prefix (e.g. `ko00010` → `hsa00010`), which is
  KEGG's own convention for organism-specific pathway projection.
* `find_sample_column()` (internal) no longer picks up categorical
  columns or columns with only partial overlap when auto-detecting
  the sample identifier column in metadata. The scan-every-column
  fallback now requires unique values and >= 90% overlap with the
  abundance sample names; standard-named columns (`sample`,
  `sample_id`, etc.) retain the previous lenient threshold because
  the column name is itself strong evidence. This prevents
  `align_samples()` from mistakenly treating a `subject_id` or
  `batch` column as the sample ID when it happens to share a few
  strings with the sample names.
* `pathway_daa(daa_method = "edgeR", reference = ...)` now honors the
  user-supplied reference. edgeR's `exactTest()` took `pair = c(1, 2)`
  against the raw factor order, so `reference` was silently ignored
  and the result always labeled `group1 = Level[1]` / `group2 = Level[2]`.
  The grouping factor is now releveled so the tested contrast is
  `log(non-ref / ref)` and the labels reflect the requested direction.
  Multi-group edgeR runs now emit one block per (reference,
  non-reference) contrast rather than every pairwise combination,
  matching the shape returned by DESeq2 / limma voom / LinDA / Maaslin2.
* `pathway_daa(daa_method = "metagenomeSeq", reference = ...)` now
  honors the user-supplied reference. The model matrix used the raw
  factor order and the result labels were hardcoded to
  `Level[1]` / `Level[2]`, so flipping `reference` left both labels
  and coefficients unchanged. The grouping factor is releveled before
  `fitFeatureModel()` so the contrast and labels track the requested
  direction.
* `taxa_contribution_heatmap(annotation_data = ...)` now accepts the
  column shape actually produced by `pathway_annotation()`. The
  heatmap used to look up `annotation_data$pathway` / `$description`,
  but `pathway_annotation()` emits `feature` (ID) and `description`
  (non-ko_to_kegg) or `pathway_name` (ko_to_kegg). The lookup silently
  returned all-NA, leaving raw IDs on the axis. The heatmap now
  detects `feature`/`description` and `pathway`/`pathway_name`
  column pairs and relabels correctly.

# ggpicrust2 2.5.13

## Bug Fixes

* `ggpicrust2()` now rejects the incompatible combination of
  `ko_to_kegg = TRUE` with `pathway = "EC"` or `"MetaCyc"` up front, with an
  actionable error message. Previously this misuse produced a cryptic
  `No features in abundance data` error thrown from deep inside
  `pathway_daa()` (reported in #198).
* `ko2kegg_abundance()` now errors when none of the input feature IDs match
  the expected KO format (e.g. when EC numbers are passed in), instead of
  silently producing an empty result. Total-mismatch detection in the
  internal `validate_feature_ids()` helper was previously a blind spot.
* `ko2kegg_abundance()` also errors when the input contains only KO IDs that
  are absent from the KEGG reference, rather than returning an empty data
  frame that would break downstream DAA.

# ggpicrust2 2.5.12

## CRAN Resubmission

* Removed transient social-media URLs from the README in favor of stable
  project links.
* Made the `metagenomeSeq` workflow a truly optional runtime dependency
  instead of a declared package-level suggestion.

# ggpicrust2 2.5.11

## Documentation

* Unified the package website URL around the canonical
  `https://cafferyang.com/ggpicrust2/`.
* Simplified README citation guidance to use the paper DOI and
  `citation("ggpicrust2")` instead of redundant publisher-specific links.
* Replaced a moved LinDA reference URL with its DOI-based link and cleaned a
  KO-to-GO manual encoding warning.

# ggpicrust2 2.5.10

## Bug Fixes

* Removed `Maaslin2` from `Suggests` to avoid CRAN/BioC availability failures
  on special check flavors where the package is not in mainstream repos.
* Switched the internal MaAsLin2 call in `pathway_daa()` to dynamic lookup so
  the method remains optional without forcing repository availability checks.

# ggpicrust2 2.5.9

## Bug Fixes

* Refactored `compare_daa_results()` examples to use minimal in-memory DAA-like
  result tables instead of running external method pipelines.
* Removed hard dependency on optional method packages in examples (notably
  `Maaslin2`) so CRAN special `donttest` checks can run without failing.

# ggpicrust2 2.5.8

## Bug Fixes

* Fixed CRAN pretest failure in `tests/testthat/test-pathway_daa.R` by splitting
  default/core method coverage from optional extended methods that require
  non-mainstream dependencies (e.g., `Maaslin2`).
* Kept extended DAA method coverage available behind
  `GGPICRUST2_RUN_EXTENDED_DAA_TESTS=true`.
* Corrected `metacyc_reference` documentation to match actual data columns
  (`id`, `description`), resolving `codoc` mismatch warnings.

# ggpicrust2 2.5.7

## Bug Fixes

* Added a regression test for `pathway_errorbar()` to explicitly cover the
  `ko_to_kegg = TRUE` + `order = "pathway_class"` path.
* This guards against reintroducing the historical
  `tibble::column_to_rownames()` / `Can't find column '.'` failure mode.

# ggpicrust2 2.5.6

## Breaking Changes

* **Removed `ggpicrust2_extended()` function**:
  - This wrapper function provided minimal value over calling `ggpicrust2()` and `pathway_gsea()` separately
  - Users can achieve the same functionality by calling the individual functions directly
  - This change reduces maintenance burden and improves code clarity

## Major Features

### Covariate Adjustment & Improved Statistical Methods for pathway_gsea() (#193)

* **Added limma camera and fry methods as new GSEA options**:
  - `method = "camera"` (now default): Competitive gene set test using limma's camera function
  - `method = "fry"`: Fast rotation gene set test (self-contained)
  - Both methods account for inter-gene correlations, providing more reliable p-values than preranked GSEA

* **Added covariate adjustment support**:
  - New `covariates` parameter to adjust for confounding factors (age, sex, BMI, etc.)
  - Covariates are incorporated into the design matrix for proper statistical adjustment
  - Essential for microbiome studies where host factors can confound results

* **New parameters**:
  - `covariates`: Character vector of covariate column names from metadata
  - `contrast`: For multi-group comparisons, specify the contrast to test
  - `inter.gene.cor`: Inter-gene correlation for camera method (default: 0.01)

* **Scientific background**:
  - Wu et al. (2012) demonstrated that preranked GSEA methods can produce "spectacularly wrong p-values" due to not accounting for inter-gene correlations
  - The camera and fry methods from limma address this limitation
  - Reference: Wu, D., & Smyth, G. K. (2012). Nucleic Acids Research, 40(17), e133.

* **Backward compatibility**:
  - Existing `fgsea` and `clusterProfiler` methods remain available
  - Added informational message when using preranked methods about p-value reliability

* **New internal functions**:
  - `run_limma_gsea()`: Core implementation for camera/fry methods
  - `build_design_matrix()`: Constructs design matrix with covariates

* **Updated documentation and vignettes**:
  - Method selection guide with comparison table
  - Covariate adjustment examples
  - Updated gsea_analysis.Rmd vignette

### New Visualization Functions

* **Added `pathway_volcano()` function**:
  - Creates publication-quality volcano plots for differential abundance analysis
  - Visualizes both statistical significance (-log10 p-value) and effect size (log2 fold change)
  - Smart label placement using ggrepel to avoid overlapping labels
  - Color-coded significance categories (Up/Down/Not Significant)
  - Automatic handling of NA pathway names (won't display "NA" labels)
  - Handles infinite p-values (when p = 0) gracefully
  - Customizable thresholds, colors, and appearance

* **Added `pathway_ridgeplot()` function**:
  - Creates ridge plots (joy plots) for GSEA results interpretation
  - Shows distribution of gene abundances/fold changes within enriched pathways
  - Color-coded by enrichment direction (Up/Down)
  - Automatic pathway-KO mapping using built-in ko_to_kegg_reference data
  - Supports KEGG and GO pathway types
  - Helps identify whether pathways are predominantly up- or down-regulated
  - Requires ggridges package (added to Suggests)

* **New dependencies added to Suggests**:
  - `ggridges`: Required for ridge plot visualization
  - `ggrepel`: Used for smart label placement in volcano plots

---

# ggpicrust2 2.5.5

## Major Features

### Prokaryote-Specific Pathway Filtering (#191)

* **New `filter_for_prokaryotes` parameter in `ko2kegg_abundance()`**:
  - Defaults to TRUE, automatically filtering out eukaryote-specific pathways
  - Removes biologically irrelevant pathways from bacterial/archaeal analysis:
    - Cancer pathways (overview and specific types)
    - Neurodegenerative diseases (Alzheimer's, Parkinson's, etc.)
    - Substance dependence (addiction pathways)
    - Cardiovascular diseases
    - Endocrine and metabolic diseases (human-specific)
    - Immune diseases (human-specific)
    - Organismal systems (immune, nervous, endocrine, digestive, etc.)
  - Retains prokaryote-relevant pathways:
    - All Metabolism pathways
    - Infectious disease: bacterial (Salmonella, E. coli, Tuberculosis, etc.)
    - Drug resistance: antimicrobial (antibiotic resistance)
    - Genetic/Environmental Information Processing
    - Cellular Processes
  - Set `filter_for_prokaryotes = FALSE` for eukaryotic analysis or to include all pathways
  - Reduces pathway count from ~370 to ~290 for typical bacterial analyses

### Local KEGG Database Implementation (#113)

* **Replaced KEGG API dependency with local database**:
  - Implemented comprehensive local KO-to-KEGG pathway mapping (61,655 mappings)
  - Covers 557 pathways and 27,127 KO IDs
  - Eliminates dependency on external KEGG API
  - Provides 100% data coverage with 0% missing values (vs 84% NA in previous format)
  - Significantly improved performance and reliability

* **Enhanced data structure**:
  - Migrated from wide format (306×326 matrix) to long format (61,655×9 table)
  - Added rich pathway metadata including hierarchical classification (Level1-3)
  - Includes pathway names, KO descriptions, and EC numbers
  - Optimized with fast lookup index for O(N) performance

* **Updated functions**:
  - `ko2kegg_abundance()`: Now uses internal database instead of KEGG API
  - `pathway_gsea()`: Updated to use long-format data for gene set enrichment
  - Both functions maintain backward compatibility
  - Added comprehensive input validation and error handling

* **PICRUSt 2.6.2 compatibility**:
  - Automatic detection and cleaning of "ko:" prefix in KO IDs
  - Handles both old (K##### format) and new (ko:K##### format) PICRUSt2 outputs
  - Seamless migration path for existing users

* **Data quality improvements**:
  - Validates KO ID format (K##### pattern)
  - Detects negative abundance values
  - Reports missing values with detailed statistics
  - Identifies all-zero KOs across samples
  - Checks for duplicate column names

* **Performance**:
  - Processing speed: 2,145-6,452 KOs/second
  - Handles large datasets (1000+ KOs, 50+ samples) efficiently
  - Progress bar for long-running operations

* **Testing**:
  - Added comprehensive test suite with 64 tests
  - Covers data structure, functionality, performance, and edge cases
  - All tests passing with 100% coverage of new features

This resolves Discussion #113 and provides a robust, API-independent solution for KEGG pathway analysis.

# ggpicrust2 2.5.4

## Improvements

### Enhanced Error Messages for pathway_annotation() (#142)

* **Significantly improved diagnostic messages when no significant pathways are found**:
  - Provides clear explanation when p_adjust < 0.05 filter removes all pathways
  - Shows detailed statistics (total features, significant count, minimum p-value)
  - Lists possible reasons (sample size, effect size, variability, method choice)
  - Offers concrete recommendations (data quality checks, alternative methods)
  - Suggests using ko_to_kegg = FALSE for local annotation without KEGG API

* **Enhanced error messages for KEGG API failures**:
  - Detailed diagnostic information when all queries fail
  - Distinguishes between "not found (HTTP 404)" and "network errors"
  - Lists affected KO IDs for troubleshooting
  - Provides specific recommendations based on error type
  - Suggests local annotation as reliable alternative

* **Updated documentation**:
  - Clarified p_adjust < 0.05 filtering behavior in function description
  - Added note about NA columns when no significant pathways exist
  - Improved parameter descriptions for ko_to_kegg parameter
  - Enhanced return value documentation with filtering behavior

* **User experience improvements**:
  - All messages follow R package standards (warning(), stop() instead of cat())
  - No emoji characters (professional text output)
  - Clear formatting for readability
  - Actionable information to help users resolve issues quickly

This addresses Discussion #142 where users received NA annotations without understanding why.

## Bug Fixes

### pathway_gsea() MetaCyc Support Fix (#174)

* **Fixed undefined organism parameter bug**:
  - Added `organism = "ko"` parameter with default value
  - Prevents "object 'organism' not found" error
  - Properly passes organism to prepare_gene_sets() function

* **Added input validation for MetaCyc data**:
  - Detects when MetaCyc pathway IDs are provided instead of EC numbers
  - Issues clear warning directing users to use pathway_daa() for pathway-level data
  - Helps users understand the difference between gene-level and pathway-level analysis
  
* **Enhanced documentation**:
  - Clarified that GSEA requires gene-level data (EC numbers for MetaCyc)
  - Added cross-reference to pathway_daa() for pathway abundance analysis
  - Improved parameter descriptions to prevent data type confusion

* **Scientific integrity maintained**:
  - Ensures correct analysis method for each data type
  - Prevents misleading results from incorrect data usage
  - Guides users to appropriate functions based on their data

### Critical annotation_custom() Fix (#184)

* **Fixed annotation_custom() parameter type issue in pathway_errorbar()**:
  - Removed unit object wrappers from `annotation_custom()` position parameters
  - Now uses numeric values directly for xmin, xmax, ymin, ymax parameters
  - Resolves "no applicable method for 'rescale' applied to an object of class 'c('simpleUnit', 'unit', 'unit_v2')'" error
  - Fixes pathway class background color rendering failures

* **Root cause identified and resolved**:
  - `annotation_custom()` expects numeric values, not unit objects for position parameters
  - Previous fix attempt (changing ggplot2::unit to grid::unit) was incorrect
  - Both ggplot2::unit() and grid::unit() return identical objects - the issue was using unit objects at all
  - Theme-related unit usage (legend.key.size, plot.margin) remains unchanged and correct

This fix resolves the critical rendering issue where users encountered errors when generating
pathway error bar plots with pathway class backgrounds, particularly with R 4.4+ and ggplot2 4.0.0.

# ggpicrust2 2.5.3

## Previous Release
* Initial attempt to fix issue #184 (superseded by v2.5.4)

# ggpicrust2 2.5.2

## Bug Fixes

### Critical Edge Case Resolution

* **Fixed pathway_errorbar empty data handling**:
  - `pathway_errorbar()` now returns NULL instead of crashing when no annotation data is available
  - Main `ggpicrust2()` function gracefully handles NULL plot objects
  - Users still receive complete results data even when plots cannot be generated
  - Improved warning messages for better user experience

* **Enhanced robustness for datasets with no significant pathways**:
  - Prevents crashes when all pathways have p_adjust > 0.05
  - Returns meaningful results with empty annotation columns for visualization
  - Maintains data integrity throughout the analysis pipeline
  - Works seamlessly with all PICRUSt2 versions including 2.6.2

* **Function signature consistency**:
  - Fixed `pathway_annotation()` function calls in main function
  - Ensures proper parameter passing throughout the workflow
  - Resolves compatibility issues between internal functions

These fixes ensure the package works reliably with all types of microbiome data, 
including edge cases where no statistically significant pathways are found.

# ggpicrust2 2.5.1

## Bug Fixes

### PICRUSt 2.6.2 Compatibility (#174) - Complete Resolution

* **Added automatic KO ID format detection and conversion**:
  - Automatically detects PICRUSt 2.6.2 format with "ko:" prefixes
  - Transparently removes "ko:" prefixes during data loading
  - Maintains full backward compatibility with PICRUSt 2.5.2 format
  - Eliminates "subscript out of bounds" errors caused by format mismatches

* **Enhanced core functions for seamless compatibility**:
  - Updated `ko2kegg_abundance()` with automatic format detection
  - Updated `ggpicrust2()` main function with compatibility layer
  - Added clear informational messages about format conversion
  - Zero manual preprocessing required for users

* **Comprehensive testing and validation**:
  - Tested with real PICRUSt 2.6.2 output files (5,000+ KO features)
  - Verified compatibility with all major DAA methods
  - Confirmed backward compatibility with existing workflows
  - Performance optimized for large datasets

# ggpicrust2 2.5.0

## Bug Fixes

### PICRUSt 2.6.2 Compatibility (#174)

* **Fixed compatibility issues with PICRUSt 2.6.2 output**:
  - Added comprehensive data validation for PICRUSt compatibility
  - Improved zero-abundance data filtering with fallback strategies
  - Enhanced error handling with specific PICRUSt version guidance
  - Added graceful handling of sparse data scenarios

* **Enhanced ALDEx2 error detection**:
  - Better error messages for insufficient or invalid data
  - Improved handling of edge cases (single features, all-zero data)
  - Added PICRUSt version-specific troubleshooting guidance

* **Improved LinDA analysis robustness**:
  - Better filtering of zero-abundance features before analysis
  - Enhanced data validation and error reporting
  - Added warnings for very sparse data scenarios

* **Enhanced ko2kegg_abundance function**:
  - Added fallback strategies for zero-abundance KO data
  - Improved compatibility warnings and error messages
  - Better handling of PICRUSt format variations

* **Added comprehensive compatibility guide**:
  - Created detailed troubleshooting documentation
  - Provided alternative analysis strategies
  - Added version-specific recommendations

# ggpicrust2 2.4.0

## New Features

### Enhanced Multi-Grouping Support for pathway_heatmap()

* **Added secondary_groups parameter (#171)**:
  - Enables multi-level grouping in pathway heatmaps
  - Supports nested faceting with multiple grouping variables
  - Maintains full backward compatibility with existing code
  - Automatic color generation for complex grouping structures

* **Deprecated facet_by parameter**:
  - Replaced with more flexible secondary_groups parameter
  - Shows deprecation warning with migration guidance
  - Will be removed in future major version

* **Enhanced Examples and Documentation**:
  - Added comprehensive examples for multi-grouping usage
  - Migration guide from facet_by to secondary_groups
  - Updated function documentation with new parameter descriptions

* **Comprehensive Testing**:
  - Added unit tests for all new functionality
  - Integration tests with real data
  - Backward compatibility validation
  - Error handling and edge case testing

## Bug Fixes

* **Improved color allocation for multi-level grouping**:
  - Smart color generation based on grouping complexity
  - Automatic color repetition when insufficient colors provided
  - Better handling of color assignment in nested structures

---

# ggpicrust2 2.3.3

## Major Bug Fixes and Improvements

### DAA Methods Comprehensive Fixes

* **Fixed limma voom compatibility with compare_daa_results (#163)**:
  - Added missing group1 and group2 columns for two-group comparisons
  - Ensures consistent output format with other DAA methods
  - Resolves "Unknown comparison type" error when using limma voom with compare_daa_results

* **Fixed Maaslin2 multiple critical issues (#164)**:
  - Fixed sample name matching between abundance matrix and metadata
  - Implemented intelligent feature name matching to handle Maaslin2's hyphen-to-dot conversion
  - Resolves "Unable to find samples in data and metadata files" error
  - Eliminates NA p-values caused by feature name mismatches
  - Robust handling of various feature naming conventions (hyphens, dots, underscores)

* **Enhanced DESeq2 robustness**:
  - Added fallback mechanism for dispersion estimation failures
  - Uses gene-wise estimates when standard dispersion estimation fails
  - Improved compatibility with small sample sizes and low-variance data

* **Standardized Lefser implementation**:
  - Added lefser package to method detection list
  - Standardized output format to match other DAA methods
  - Returns results for all features, not just significant ones
  - Converts effect scores to p-values for consistency

### Comprehensive Testing and Validation

* **100% success rate achieved for all 8 DAA methods**:
  - ALDEx2, DESeq2, edgeR, limma voom, metagenomeSeq, Maaslin2, LinDA, Lefser
  - All methods now fully compatible with compare_daa_results function
  - Consistent output format across all methods
  - Extensive testing with various feature naming patterns and edge cases

### Technical Improvements

* Enhanced error handling and robustness across all DAA methods
* Improved feature name matching algorithms
* Better sample-metadata alignment validation
* Defensive programming against edge cases with factor levels

# ggpicrust2 2.3.2

## Bug Fixes

* Fixed MetaCyc pathway annotation NA description issue (#154):
  - Standardized MetaCyc reference data column names from 'X1'/'X2' to 'id'/'description'
  - Applied fix in all three loading paths of load_reference_data function
  - Resolves column name mismatch that caused all MetaCyc annotations to return NA
  - Tested with 100% success rate on sample data
  - Maintains backward compatibility with KO and EC pathway types

# ggpicrust2 2.3.1

## Bug Fixes

* Fixed MetaCyc reference data loading issue:
  - Enhanced the file search mechanism in the `load_reference_data` function
  - Added multiple search paths for reference data files
  - Improved error messages with more diagnostic information
  - Fixed "Reference data file not found" error that some users encountered

# ggpicrust2 2.3.0

## Bug Fixes

* Fixed NA handling in pathway_errorbar function:
  - Improved handling of NA values in the feature column
  - Added robust error checking for group ordering option
  - Prevents errors when processing MetaCyc pathway data with missing annotations

* Enhanced pathway_pca function to handle zero variance data:
  - Automatically detects and filters out rows (pathways) with zero variance
  - Keeps sample profiles as PCA observations, including samples with zero
    variance across pathways
  - Provides informative warnings about removed pathways
  - Improves error handling with clear diagnostic messages

* Fixed file extension handling in ko2kegg_abundance function:
  - Resolved "the condition has length > 1" error when processing files
  - Improved extension detection for different file types
  - Enhanced robustness for handling various input formats

# ggpicrust2 2.2.2

## Bug Fixes

* Fixed column name compatibility issues in reference data files:
  - Standardized column names in KO_reference.RData and EC_reference.RData
  - Resolved warnings when using pathway_annotation() with ko_to_kegg=FALSE
  - Improved backward compatibility with existing code

# ggpicrust2 2.2.1

## Bug Fixes

* Added aggregate_by_group parameter to pathway_heatmap function:
  - Allows displaying representative samples (e.g., mean) for each group instead of all individual samples
  - Improved handling of NA values in heatmap visualization
  - Added aggregate_fun parameter for customizing the aggregation function

# ggpicrust2 2.1.4

## Reference Data Updates

* Updated reference databases for improved pathway annotation:
  - EC reference data updated from 3,180 to 8,371 entries (163% increase)
  - KO reference data updated from 23,917 to 27,531 unique KO IDs (15.4% increase)
  - These updates provide more comprehensive and accurate pathway annotations

# ggpicrust2 2.1.1

## Bug Fixes

* Improved error handling for KEGG database connections (Issue #138):
  - Added robust error handling for HTTP 404 errors when connecting to the KEGG database
  - Function now continues processing other KO IDs even when some IDs return HTTP 404 errors
  - Added detailed logging about which KO IDs were not found
  - Only throws a fatal error when all KO IDs fail to be processed
  - Includes summary statistics about successful, not found, and error counts

* Fixed the LinDA analysis for multi-group comparisons (Issue #144):
  - Modified the `perform_linda_analysis` function to handle multi-group comparisons correctly
  - The function now creates separate result entries for each pairwise comparison
  - Added the `log2FoldChange` column to the results for effect size information

# ggpicrust2 2.1.0

## Major Changes

* Added Gene Set Enrichment Analysis (GSEA) functionality with the following new functions:
  - `pathway_gsea()`: Performs GSEA analysis, supporting KEGG, MetaCyc, and GO pathways
  - `visualize_gsea()`: Creates visualizations of GSEA results, including enrichment plots, dotplots, barplots, network plots, and heatmaps
  - `compare_gsea_daa()`: Compares GSEA and Differential Abundance Analysis (DAA) results
  - `gsea_pathway_annotation()`: Adds pathway annotations to GSEA results

* Improved network and heatmap visualization capabilities with richer parameter options and better error handling

* Added preliminary support for MetaCyc and GO pathways

* Fixed various bugs and optimized code structure

# ggpicrust2 2.0.1

## Major Changes

* Fixed a bug in the pathway_annotation function, resolving issues caused by changes in the KEGG API response structure

# ggpicrust2 2.0.0

## 主要变更

* Refactored the package dependencies, moving most Bioconductor packages from Imports to Suggests, reducing mandatory dependencies.

* Added conditional checks to ensure that packages are only used when they are available.

* Fixed the function export issue, ensuring that all public API functions are correctly exported.

* Updated the example code, improving its stability and compatibility.

* Updated the documentation format to comply with the latest CRAN standards.

* Fixed invalid URL links in the README.

* Optimized code quality, removing unused variables.

* Added missing import declarations to ensure package integrity.

# ggpicrust2 1.7.5

# ggpicrust2 1.7.2

# ggpicrust2 1.7.1

# ggpicrust2 1.7.0

# ggpicrust2 1.6.6

# ggpicrust2 1.6.5

# ggpicrust2 1.6.4

# ggpicrust2 1.6.3

# ggpicrust2 1.6.2

# ggpicrust2 1.6.1

# ggpicrust2 1.6.0

# ggpicrust2 1.5.1

# ggpicrust2 1.5.0

# ggpicrust2 1.4.12

# ggpicrust2 1.4.11

# ggpicrust2 1.4.10

# ggpicrust2 1.4.9

# ggpicrust2 1.4.8

# ggpicrust2 1.4.7

# ggpicrust2 1.4.6

# ggpicrust2 1.4.5

# ggpicrust2 1.4.4

# ggpicrust2 1.4.3

# ggpicrust2 1.4.2

# ggpicrust2 1.4.1

# ggpicrust2 1.4.0

# ggpicrust2 1.3.1

# ggpicrust2 1.3.0

# ggpicrust2 1.2.3

# ggpicrust2 1.2.2

# ggpicrust2 1.2.1

# ggpicrust2 1.2.0

# ggpicrust2 1.1.1

# ggpicrust2 1.1.0

# ggpicrust2 1.0.1

# ggpicrust2 1.0.0

# ggpicrust2 0.2.0

* Added a `NEWS.md` file to track changes to the package.
