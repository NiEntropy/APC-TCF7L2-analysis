# Experimental analyses: Figures 2–5 and S4–S12

This additive package complements the Fig1/S1–S3 DepMap code at the repository root. It shares numerical figure inputs, final statistical outputs, transcriptomic result tables, and historical figure builders. It is not a complete archive of original instrument files or a FASTQ-to-manuscript pipeline.

## Run the verified numerical reanalysis

From the repository root, in an isolated Python environment:

```sh
python -m pip install -r experimental/requirements.txt
python experimental/reproduce_statistics.py
```

The script reads only the packaged tables, recomputes final Fig2F, Fig3D/F, Fig4C/D/F, Fig5G statistics and DEG/overlap counts, checks agreement with archived results, and writes `experimental/recomputed/`. It does not read caches or require the original Chinese Windows directory. Versions in `environment-tested.json` describe the environment in which this check passed on 2026-10-12; they are not asserted to be the historical sequencing-provider environment. A fresh dependency installation was not separately tested.

## Source-data map

| Panels | Files under `data/` | Experimental unit and interpretation |
|---|---|---|
| Fig2F, S4–S6 SRB | `Fig2/Fig2_growth_data.csv`, `Fig2_source_data.csv` | Raw recorded OD575 per well at 12, 24, 48, 72 and 96 h; relative 96-h signal is KD / matched NC mean ×100. Six technical wells per group; do not claim six independent biological experiments. |
| S4–S6 WB, qPCR, colonies | `Fig2/FigS4_S6_source_data.csv` | Source-normalized observations. The observation index alone does not establish biological independence. |
| Source-supplied stars | `Fig2/Significance_source_map.csv` | Preserved annotations; this table does not document or verify the original statistical test. |
| Fig3D/F and S7/S8 numerical observations | `Fig3/plotted_data.csv` | Flow quadrant percentages, protein ratios and normalized qPCR. Fig3D total apoptosis = early + late; Fig3F ratios normalized to the NC mean. Three biological replicate records, paired by experimental batch. |
| Fig4C | `Fig4/Fig4C_individual_volumes.csv` | Eight animals/group, ten time points. Day0 zeros come from the Prism plotting baseline, not measured tumour volumes. |
| Fig4D | `Fig4/Fig4D_Prism_weights.csv` | Final Prism endpoint weights in grams, eight animals/group. This is the source selected for the final main figure. |
| Fig4F | `Fig4/TUNEL_fields.csv` | Three fields per tissue; average fields first, then compare eight tissue samples/group. Fields are not independent animals. |
| S9 and source comparison | `Fig4/animal_measurements_Excel.csv` | Excel volume/body-weight record and Excel endpoint weights. Excel endpoint weights differ from final Fig4D Prism inputs; this file must not replace the final Fig4D source. |
| Fig5G | `Fig5/Fig5G_independent_samples.csv` | Six independent biological samples/group, source-normalized mRNA expression; no new Ct normalization is performed. |
| Fig5, S10–S12 | `RNAseq/*DEG_full.csv.gz`, `*GO_original.csv.gz`, `*GSEA_original.csv.gz`, PCA/correlation tables, shared genes/GO, heatmaps | Provider-derived DEG/GO/full-GO GSEA tables and manuscript-derived intersections/heatmaps. PCA coordinates/correlations are supplied results, not recalculated from raw counts by the verification script. |

The older qPCR averaging scheme, log-transformed Fig5G tests, mixed/repeated-measures Fig4 analyses and alternative endpoint-source tests are intentionally excluded from the current statistics package. Historical builders may contain earlier intermediate analyses; final numerical tables and the method table below take precedence.

## Final main-figure statistical methods

| Panel | Final method | Multiple comparisons |
|---|---|---|
| Fig2F | Two-sided unpaired Welch t tests on raw 96-h OD575, seven comparisons | Holm across the seven tests; exploratory technical-well inference |
| Fig3D/F | Two-sided paired t tests, three biological pairs matched by batch | Unadjusted P values |
| Fig4C | Ordinary two-way ANOVA, group × time interaction | Šídák adjustment across ten between-group time-point comparisons using the ANOVA pooled residual variance |
| Fig4D | Two-sided unpaired equal-variance Student t test | Unadjusted |
| Fig4F | Two-sided unpaired Welch t test on tissue-level field means | Unadjusted |
| Fig5G | Two-sided unpaired equal-variance Student t tests, six biological samples/group | Unadjusted, separately for seven genes |

Ordinary Fig4C ANOVA follows the final main figure; it does not model within-animal longitudinal correlation. Source-supplied stars in Fig2B–E/S4–S6/S7C–D/S8A–D have not been independently assigned a statistical method. A missing test must not be inferred from a star. See the companion RNA-seq methods record for FDR families and GSEA limitations.

## Historical figure builders

`archive/` preserves the relevant original figure-composition and analysis scripts, with source hashes in `archive/provenance.json`. These are provenance records, not standalone public entry points: they depend on original Excel/Prism files, microscopy/blot images, zipped sequencing reports and the original project folder layout; some reuse caches. A historical cache hit is not new analysis. Plotting requires additional packages (matplotlib, openpyxl, xlrd, Pillow, pymupdf) and suitable fonts. Figure assembly has not been rerun from a fresh public checkout. `reproduce_statistics.py` is the portable, numerically checked entry point.

## Access and licensing

RNA sequencing is deposited under [GEO GSE348858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348858). An accession is not proof that all files are already publicly accessible; confirm release status when submitting. FASTQ/BAM, original FCS/FlowJo files, uncropped blot images and full microscopy archives are not included here. An appropriate data archive can host those raw files and provide a stable accession/DOI linked from this repository.

No unpublished manuscript drafts, contact sheets, private GEO reviewer links, credentials or third-party article PDFs are included. No reuse license has been selected; existing repository and third-party resource terms remain applicable. SHA-256 checksums in `manifest.json` identify this package's files; they do not certify original biological provenance.

## 中文说明

本目录补充 Fig2–Fig5 及 S4–S12 的数值源数据、最终统计结果、转录组分析结果与原始作图脚本。主图统计方法以本页表格为准。已验证复算脚本能独立读取本目录数据，复现主图已明确的统计结果和差异基因/交集数量；不宣称已从 FASTQ 重新执行上游测序分析。Fig2F 的六个孔是技术重复；Fig3D/F 为按批次配对的三组生物学重复；Fig5G 为每组六个独立生物学样本。旧版 Fig4/Fig5 统计不作为当前结果。原始 FCS、完整 WB 和显微图像仍需另行归档共享。
