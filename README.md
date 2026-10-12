# APC–TCF7L2 dependency analysis

Research code for the merged Figure 1 and supplementary Figures S1–S3, with associated exploratory multi-omics analyses.

## Scope

APC functional classification, genome-wide LOF-versus-WT screening, lineage-adjusted regression, TCF/Wnt-family specificity, Hallmark/PROGENy-derived/seven-gene Wnt scores, RNA/protein associations, co-dependencies and editable scientific figures.

Wet-laboratory RNA-seq, animal, flow-cytometry and other experimental-analysis code maintained outside this workspace is not included. This repository must not be cited as the complete code for every manuscript experiment.

## Installation

Use an isolated Python environment and `python -m pip install -r requirements.txt`. The pinned versions come from the historical `environment-recorded.json`; the original Python, reportlab and openpyxl versions were not captured. A fresh dependency installation and end-to-end numerical run have not been verified for this public package. Arial was used for the figures; font substitution requires layout review. The archived supplementary compositor uses font discovery instead of a machine-specific Windows path.

## Data and reproduction

Original matrices, experimental data, result figures and unpublished manuscript drafts are not bundled. Obtain datasets under their original terms. See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for inputs and known gaps.

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
