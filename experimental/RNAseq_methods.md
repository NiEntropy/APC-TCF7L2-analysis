# RNA-seq analysis provenance

## Samples and upstream processing

HT29 and SW620 were analysed separately, with three independent biological samples per condition: shNC versus HT29 shTCF7L2 #1 or SW620 shTCF7L2 #6 (12 samples total). Sequencing data are deposited as GSE348858. Poly(A)-enriched, non-stranded libraries were sequenced using paired-end 150-bp reads on NovaSeq X Plus.

Provider records identify fastp 0.19.7 (`-g -q 5 -u 50 -n 15 -l 150`), HISAT2 2.2.1 (`--phred33 --dta`), SAMtools 1.20 and featureCounts 2.0.6. Reference annotation was GRCh38.p12 / GCA_000001405.27, gene annotation updated July 2018. An Ensembl release number is not inferred from the assembly name. Counting excluded low-mapping-quality (<10), multimapping and unpaired reads according to the provider method record; a complete counting command is not retained here.

DESeq2 1.42.0 used raw counts separately within each cell line, retaining genes with total counts >1 across the six samples. The exact design formula is not retained in this package and is not reconstructed by assumption. Differential expression: BH-adjusted P ≤0.05 and |log2FC| ≥1, with positive log2FC indicating higher expression after TCF7L2 knockdown. Packaged tables are delivered analysis results, not the output of a new DESeq2 run.

The reproducibility script verifies HT29 375 up/195 down, SW620 1257 up/664 down, with 99 shared upregulated and 17 shared downregulated Ensembl gene IDs. Heatmaps use log2(FPKM+1), gene-wise z scores separately within each cell line, clipped to ±2 (zero variance mapped to zero). Selected Fig5F genes need not all meet differential-expression thresholds in both cell lines.

## GO enrichment

Provider GO analysis used clusterProfiler 4.8.1 and hypergeometric over-representation. BP background sizes were 14285 (HT29) and 14688 (SW620). BH correction was applied within each ontology (BP/CC/MF), rather than treating all GO categories as a single undocumented family.

Shared-gene BP enrichment used the intersection of genes with nonmissing adjusted P values in both DEG tables and BP annotations: 13209 background genes. Of the shared genes, 89/99 upregulated and 14/17 downregulated genes mapped to BP. Retain terms containing 10–500 background genes, including zero-hit terms in the testing family. Compute `hypergeom.sf(k-1, N, M, n)` and BH-adjust all retained terms separately for shared-up and shared-down lists. Gene ratio = k/n. The historical builder contains the implementation; re-running this enrichment also requires the original provider `go.xls` annotation mapping, which is not bundled. Packaged shared-GO tables preserve the original completed analysis.

## Distinguish the two GSEA analyses

| Record | Fig3 custom apoptosis analysis | S11 provider full-GO analysis |
|---|---|---|
| Input | Archived `express.gct` | Archived `readcount.gct` |
| Gene sets | Six HT29 and two SW620 selected cell-death sets | Full GO collection; 6667 HT29 / 6763 SW620 sets in archived GMT files |
| Minimum set size | HT29 1; SW620 15 | 15 in both cell lines |
| Maximum set size | 5000 | 5000 |
| Software version | Custom rerun version not independently retained | GSEA 4.3.2 in provider report |

Both use 1000 phenotype permutations, Signal2Noise ranking in descending order, weighted enrichment, mean-division normalization and NoCollapse; seeds were timestamp-derived. With only three samples per phenotype, the number of distinct label partitions is limited; 1000 requested permutations do not represent 1000 independent biological assignments. Nominal P=0 should not be read as an exact probability of zero.

Fig3 running enrichment curves were reconstructed from archived ranks and hit indices, with ES/NES/nominal P/FDR taken from the saved results; this is not a new GSEA permutation run. Exploratory GSEA uses q<0.25, while DEG/GO use the stated 0.05 thresholds. Nominal P is not interchangeable with FDR q. FDR controls an expected false-discovery proportion within a specified family; it does not guarantee that no more than that percentage of a particular selected list is false, nor is q a posterior probability that an individual finding is false.

Archived gene-set SHA-256 identifiers (source identification, not a redistribution license):

| GMT | SHA-256 |
|---|---|
| HT29 custom | `502237b9aad6ade2bf88de674a3c36c7d9f896d2ded9aa3e9c9971f6605432ef` |
| SW620 custom | `1c2ca96245571aba8e479d27b15b901105fd69897f5e31524a5153930e633289` |
| HT29 full GO | `9f34caf0d3a23aeb7e8e2a4ba9d4ed91720e2a68abda9361215a9f76bc606230` |
| SW620 full GO | `3c00a999231ed83130e7694a5e4a6ceb24926b77a1350e7feae8b9fdf0aa4c0c` |

## Remaining reproduction limits

FASTQ alignment/counting and DESeq2 were not rerun when packaging. Full command logs, complete annotation/GMT files and a fully specified DESeq2 design would be needed for an exact upstream rerun. Original sequencing-provider report identifiers: HT29 X101SC25102166-Z01-J001 (2025-10-25), SW620 X101SC25112550-Z01-J001 (2025-11-28). No private submission token or reviewer access URL is disclosed.
