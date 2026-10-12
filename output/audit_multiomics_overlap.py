from pathlib import Path
import re,json
import numpy as np
import pandas as pd
import pyreadr

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/multiomics_audit'
ann=pd.read_csv(ROOT/'Processed/RNA_SEQ/r_data/Cell_lines_annotations_20181226.txt',sep='\t')
assert not ann.depMapID.dropna().duplicated().any()
mapping=ann.set_index('CCLE_ID').depMapID
base=pd.read_csv(ROOT/'output/pdf/Fig1H_source_data.csv').set_index('depmap_id')
master=base[['cell_line_display_name','lineage_1','lineage_2','APC_class','TCF7L2_gene_effect']].copy()
ph=next(iter(pyreadr.read_r(str(ROOT/'Processed/ph/cell_line.ph.rds')).values())).set_index('depmap_id')
for col in ['GDSC_MSI','CCLE_MSI','TP53_status']:master[col]=ph[col]
ex=next(iter(pyreadr.read_r(str(ROOT/'Processed/expr/crispr.exp.rds')).values()))
master['CTNNB1_gene_effect']=ex.loc['CTNNB1']
rna=pd.read_csv(OUT/'CCLE_RNAseq_genes_rpkm_20180929.gct.gz.targets.csv').set_index('Description').drop(columns='Name').T
rna.index=rna.index.map(mapping)
assert not rna.index.dropna().duplicated().any()
rna=rna.loc[rna.index.notna()].apply(pd.to_numeric,errors='coerce')
for g in rna.columns:master[f'RNA_RPKM_{g}']=rna[g]
prot=pd.read_csv(OUT/'Proteomics_subsetted.csv.targets.csv',index_col=0)
pc=[c for c in prot if c.startswith('TCF7L2 ')][0];bc=[c for c in prot if c.startswith('CTNNB1 ')][0]
raw=pd.read_csv(OUT/'raw_protein_targets.csv')
samplecols=[c for c in raw if re.search(r'_TenPx\d+$',c)]
samplemap=pd.DataFrame({'protein_sample':samplecols,'CCLE_ID':[re.sub(r'_TenPx\d+$','',c) for c in samplecols]})
samplemap['depmap_id']=samplemap.CCLE_ID.map(mapping)
assert samplemap.depmap_id.notna().all()
rv=raw.loc[raw.Gene_Symbol.eq('TCF7L2'),samplecols].iloc[0].to_numpy(float)
assert np.allclose(prot[pc].to_numpy(float),rv,equal_nan=True)
samplemap['processed_id']=prot.index
samplemap.to_csv(OUT/'protein_sample_map.csv',index=False)
prot.index=samplemap.depmap_id
# Aggregate duplicate normalized measurements within the same actual cell line.
prot=prot.groupby(level=0).mean()
master['MS_TCF7L2']=prot[pc];master['MS_CTNNB1']=prot[bc]
rp=pd.read_csv(OUT/'Protein_Array_subsetted.csv.targets.csv',index_col=0)
master['RPPA_beta_catenin']=rp['beta-Catenin'];master['RPPA_beta_catenin_pT41_S45']=rp['beta-Catenin_pT41_S45']
meth=pd.read_csv(OUT/'methylation_TCF7L2.csv',index_col=0)
for col in meth:master['METH_'+col]=meth[col]
mc=['METH_'+c for c in meth]
master['RNA_available']=master.RNA_RPKM_TCF7L2.notna()
master['Protein_available']=master.MS_TCF7L2.notna()
master['Meth_any']=master[mc].notna().any(axis=1)
master['Meth_all3']=master[mc].notna().all(axis=1)
master['RNA_Protein']=master.RNA_available&master.Protein_available
master['RNA_Protein_Meth_any']=master.RNA_Protein&master.Meth_any
master['RNA_Protein_Meth_all3']=master.RNA_Protein&master.Meth_all3
for col in mc:master['Full_'+col]=master.RNA_Protein&master[col].notna()
features=['TCF7L2_gene_effect','CTNNB1_gene_effect','RNA_RPKM_TCF7L2','MS_TCF7L2','MS_CTNNB1','RPPA_beta_catenin']+mc+['Meth_any','Meth_all3','RNA_Protein','RNA_Protein_Meth_any','RNA_Protein_Meth_all3']+['Full_'+c for c in mc]
rows=[]
for scope,df in [('Pan-cancer',master),('CRC',master.loc[master.lineage_2.eq('Colorectal Adenocarcinoma')])]:
    for col in features:
        mask=df[col] if df[col].dtype==bool else df[col].notna()
        v=df.loc[mask];ct=v.APC_class.value_counts()
        rows.append(dict(scope=scope,feature=col,n_total=len(v),WT=int(ct.get('WT',0)),Missense_only=int(ct.get('Missense only',0)),LOF=int(ct.get('Truncating/LOF',0)),Other=int(ct.get('Other mutation',0)),n_WT_LOF=int(ct.get('WT',0)+ct.get('Truncating/LOF',0)),lineages=v.lineage_1.nunique()))
coverage=pd.DataFrame(rows);coverage.to_csv(OUT/'coverage_by_APC_group.csv',index=False)
master.to_csv(OUT/'master_availability_and_targets.csv')
# Exact published CCLE_ID -> DepMap ID mapping, no fuzzy cell-name matching.
rid=pd.read_csv(OUT/'CCLE_RNAseq_genes_rpkm_20180929.gct.gz.ids.csv').iloc[:,0]
mapped=rid.map(mapping)
stats={'RNA_samples':len(rid),'RNA_mapped':int(mapped.notna().sum()),'RNA_overlap_1150':int(mapped.isin(master.index).sum()),'MS_samples':len(samplecols),'MS_unique_lines':len(prot),'MS_TCF7L2_unique_nonmissing':int(prot[pc].notna().sum()),'MS_targets_equal_raw_processed':True,'protein_id_rows_repaired':int(samplemap.depmap_id.ne(samplemap.processed_id).sum()),'RPPA_features':214,'MS_features':12755,'METH_regions':len(meth.columns),'RNA_Wnt_targets':[g for g in ['AXIN2','MYC','CCND1','LGR5','ASCL2','RNF43','SP5'] if g in rna],'RNA_missing_mapping':rid[mapped.isna()].tolist(),'MS_TCF7L2_raw_range':[float(np.nanmin(rv)),float(np.nanmax(rv))]}
# Design estimability only: no outcome association tests or coefficient selection.
designs=[]
for target in mc+['all3']:
    required=['RNA_RPKM_TCF7L2','MS_TCF7L2']+(mc if target=='all3' else [target])
    df=master.loc[master.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=required+['lineage_1','TCF7L2_gene_effect'])
    X=pd.DataFrame({'intercept':1.,'APC_LOF':df.APC_class.eq('Truncating/LOF').astype(float),'RNA':np.log2(df.RNA_RPKM_TCF7L2+1),'Protein':df.MS_TCF7L2},index=df.index)
    for c in required[2:]:X[c]=df[c]
    X=pd.concat([X,pd.get_dummies(df.lineage_1,drop_first=True,dtype=float)],axis=1)
    rank=np.linalg.matrix_rank(X)
    designs.append({'methylation':target,'n':len(df),'WT':int(df.APC_class.eq('WT').sum()),'LOF':int(df.APC_class.eq('Truncating/LOF').sum()),'lineages':df.lineage_1.nunique(),'columns':X.shape[1],'rank':int(rank),'residual_df':len(df)-int(rank),'CRC_LOF':int((df.lineage_2.eq('Colorectal Adenocarcinoma')&df.APC_class.eq('Truncating/LOF')).sum())})
pd.DataFrame(designs).to_csv(OUT/'full_model_design_audit.csv',index=False)
(OUT/'mapping_summary.json').write_text(json.dumps(stats,indent=2))
print(json.dumps(stats,indent=2));print(coverage.to_string(index=False));print(pd.DataFrame(designs).to_string(index=False))
