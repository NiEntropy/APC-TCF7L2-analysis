# APC–TCF7L2 dependency analysis

Research code for Figure 1/S1–S3 and associated exploratory multi-omics analyses, plus experimental source data and final statistical reanalysis for Figures 2–5/S4–S12.

## Scope

APC functional classification, genome-wide LOF-versus-WT screening, lineage-adjusted regression, TCF/Wnt-family specificity, Hallmark/PROGENy-derived/seven-gene Wnt scores, RNA/protein associations, co-dependencies and editable scientific figures.

Experimental numerical source data and final main-figure statistical reanalysis for Figures 2–5/S4–S12 are included in [experimental/](experimental/README.md), together with RNA-seq result tables and historical figure builders. Original instrument/image files and a complete FASTQ-to-DESeq2 pipeline are not bundled; this repository is not complete code for every experimental step.

## Installation

Use an isolated Python environment and `python -m pip install -r requirements.txt`. The pinned versions come from the historical `environment-recorded.json`; the original Python, reportlab and openpyxl versions were not captured. A fresh dependency installation and end-to-end numerical run have not been verified for this public package. Arial was used for the figures; font substitution requires layout review. The archived supplementary compositor uses font discovery instead of a machine-specific Windows path.

## Data and reproduction

Original DepMap matrices, raw experimental instrument/image files, result figures and unpublished manuscript drafts are not bundled. Processed experimental source data are included under `experimental/data/`. Obtain datasets under their original terms. See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for inputs and known gaps.

After preparing the required matched inputs and upstream tables:

```sh
python output/cofig1_analysis.py
python output/pdf/build_cofig1.py
```

Keep the original `Processed/` and `output/` directory layout. This is an analysis-code archive, not a self-contained reproduction from a new DepMap download. Some input provenance and intermediate construction remain unresolved. Scripts reuse existing caches; a cache hit must not be interpreted as reanalysis of changed data.

## Interpretation

APC classes are variant-based annotations, not proof of biallelic functional loss. The genome-wide screen is unadjusted; lineage adjustment is a separate model. Genome-wide and eight-gene FDR families differ. RNA/protein group differences do not establish lineage-independent abundance differences, causality or synthetic lethality.

## Validation, citation and licensing

Packaging checks cover Python syntax, obvious credential/local-drive-path patterns and hashes (`code_manifest.json`), not numerical reproduction or a security audit. Cite the actual repository URL and commit/release after publication. No DOI or publication metadata is asserted. The authors have not selected a reuse license; public visibility alone is not an open-source license. Third-party data and gene sets retain their own terms.

## Experimental analyses

See [experimental/README.md](experimental/README.md) for data units and the final statistical-method map, and [experimental/RNAseq_methods.md](experimental/RNAseq_methods.md) for sequencing/GO/GSEA provenance. Run `python experimental/reproduce_statistics.py` after installing `experimental/requirements.txt`. Agreement with archived final main-figure tests and DEG/overlap counts was verified on 2026-10-12. Historical figure builders and upstream RNA-seq steps were not rerun.
