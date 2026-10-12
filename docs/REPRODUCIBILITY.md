# Reproducibility and provenance

## Required inputs

- `Processed/expr/crispr.exp.rds`: genes × ModelIDs. Original dependency release, scoring algorithm and Sanger Score inclusion remain unconfirmed. ScreenGeneEffect is not an interchangeable input.
- `Processed/ph/cell_line.ph.rds`: model metadata; depmap_id, lineage_1/2 and display names.
- `Processed/mut/mutations_NAsdropped.csv`: historical variants used by bootstrap stages.
- Default 26Q1 `OmicsSomaticMutations.csv` and `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`.
- `Processed/protein/Proteomics_subsetted.csv`, `protein_quant_current_normalized.csv.gz`, `CCLE_RPPA_20181003.csv`, `Protein_Array_subsetted.csv`.
- `Processed/RNA_SEQ/r_data/Cell_lines_annotations_20181226.txt`, `CCLE_RNAseq_genes_counts_20180929.gct`, `CCLE_RNAseq_genes_rpkm_20180929.gct.gz`.
- Regional methylation CSV in `Processed/Methylation/` and the CN/LoH/profile/segment/global-signature files referenced by the audit scripts. The optional new-genomics audit also uses molecular-subtype annotations.

These processed filenames do not establish the original release provenance. Download dates, RNA log base and some original dataset identifiers still require confirmation. Do not replace model-level CRISPR values with screen-level values.

## Historical sequence

1. `output/pdf/build_fig1e.py` creates the legacy discovery master; its old MUT comparison is not the final LOF analysis.
2. `output/pdf/build_fig1h.py` creates the historical functional-class source table.
3. `output/audit_multiomics.py`, then `output/audit_multiomics_overlap.py`, generate mapped molecular targets and the earlier master.
4. `output/audit_new_26q1.py` applies the current default APC/RNA definitions and produces `output/multiomics_audit_26Q1/master_audit_26Q1.csv`.
5. `output/fig2_resources.py` downloads resources. Compare downloads against the historical URLs/hashes in `resource_sources.json`; the live PROGENy endpoint may change. Access and resource license requirements apply separately.
6. `output/fig2_analysis.py` runs broad multi-omics analyses. It additionally requires `output/multiomics_audit_26Q1/APC_CN_mutation_sequencing_match.csv`. Construction of this same-sequencing mutation/CN intermediate is not fully captured by these scripts: this is a known reproduction gap. Do not synthesize allelic classifications to bypass it.
7. `output/fig2_validate.py` independently checks selected FWL fits, computes leave-one-LOF-out models and beta-ranked GSEA. `fig2_plot.py` and `pdf/build_refig2.py` reproduce earlier exploratory figure layouts.
8. `output/cofig1_analysis.py` uses `output/Fig2_updated/tables/master_analysis.csv` plus the CRISPR matrix for final Fig1 screening and regression.
9. `output/pdf/build_cofig1.py` also reads dependency_associations, APC_multivariable_models, methylation_associations, group_comparisons and codependency_ranking tables under `output/Fig2_updated/tables/`, and the newly generated `output/coFig1_analysis/` tables.

The matched master and computed tables are not included in this code-only disclosure. Use a fresh output directory for changed inputs; otherwise caches can retain stale results. The historical multi-omics analyses require more inputs than the final subset of manuscript panels.

## Method and verification notes

Final classes: WT1005/LOF76/missense66. Genome-wide tests require ≥3 finite values/group; 17930 genes are tested. Lineage models use OLS/HC1 with t inference and lineage_1. Family bootstrap: 5000 resamples, seed20260923. Wnt/co-dependency/GSEA uses seed20260922 where applicable. PROGENy is a custom weighted footprint score. The full complete-case multi-omic model did not establish an abundance-independent APC effect.

Checks performed during packaging: AST syntax, literal local-drive paths, recognized GitHub-token/private-key signatures and hashes. No new numerical analysis was run. The working environment used for packaging no longer has all original analysis libraries installed; clean installation and numerical reproducibility remain to be tested once all inputs are supplied.
