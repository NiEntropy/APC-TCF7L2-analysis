from pathlib import Path
import csv,gzip,json
import numpy as np
import pandas as pd
import openpyxl

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/multiomics_audit';OUT.mkdir(exist_ok=True)
def emit(s):print(s,flush=True)
ann=pd.read_csv(ROOT/'Processed/RNA_SEQ/r_data/Cell_lines_annotations_20181226.txt',sep='\t')
emit({'annotation_shape':ann.shape,'duplicate_CCLE_ID':int(ann.CCLE_ID.duplicated().sum()),'duplicate_depmap':int(ann.depMapID.duplicated().sum()),'annotation_examples':ann[['CCLE_ID','depMapID']].head().to_dict('records')})
for fname in ['CCLE_RPPA_20181003.csv','Protein_Array_subsetted.csv','Proteomics_subsetted.csv']:
    d=pd.read_csv(ROOT/'Processed/protein'/fname,index_col=0)
    targets=[x for x in d.columns if any(g in x.upper() for g in ['TCF7L2','CTNNB1','CATENIN','TCF4'])]
    emit({'file':fname,'shape':d.shape,'sample_ids':d.index[:5].tolist(),'targets':targets,'nonnull':{x:int(d[x].notna().sum()) for x in targets}})
    d[targets].to_csv(OUT/(fname+'.targets.csv'))
    pd.Series(d.index).to_csv(OUT/(fname+'.ids.csv'),index=False)
raw=ROOT/'Processed/protein/protein_quant_current_normalized.csv.gz'
d=pd.read_csv(raw)
emit({'file':raw.name,'shape':d.shape,'columns_first':d.columns[:10].tolist(),'columns_middle':d.columns[40:65].tolist(),'columns_last':d.columns[-10:].tolist()})
mask=d.Gene_Symbol.astype(str).str.contains(r'(^|[;, ])(TCF7L2|CTNNB1|TCF4)([;, ]|$)',regex=True)
emit({'raw_protein_targets':d.loc[mask].iloc[:,:10].to_dict('records')})
d.loc[mask].to_csv(OUT/'raw_protein_targets.csv',index=False)
pd.Series(d.columns).to_csv(OUT/'raw_protein_columns.csv',index=False)
mfile=next((ROOT/'Processed/Methylation').glob('*.csv'))
cols=pd.read_csv(mfile,nrows=0).columns
targets=[x for x in cols if x.startswith('TCF7L2_')]
m=pd.read_csv(mfile,usecols=[cols[0]]+targets,index_col=0)
emit({'methylation_shape_target':m.shape,'n_features':len(cols)-1,'sample_ids':m.index[:5].tolist(),'target_regions':targets,'nonmissing':m.notna().sum().to_dict(),'ranges':{x:[float(m[x].min()),float(m[x].max())] for x in targets}})
m.to_csv(OUT/'methylation_TCF7L2.csv')
genes={'TCF7L2','CTNNB1','AXIN2','MYC','CCND1','LGR5','ASCL2','RNF43','SP5'}
for fname in ['CCLE_RNAseq_genes_counts_20180929.gct','CCLE_RNAseq_genes_rpkm_20180929.gct.gz']:
    p=ROOT/'Processed/RNA_SEQ/r_data'/fname
    f=gzip.open(p,'rt') if p.suffix=='.gz' else open(p)
    version=f.readline().strip();dims=f.readline().strip();rd=csv.reader(f,delimiter='\t');header=next(rd);found=[];symbols=[];n=0
    for row in rd:
        n+=1;symbols.append(row[:2])
        if any(v in genes for v in row[:2]):found.append(row)
    emit({'RNA':fname,'dimensions':dims,'observed_rows':n,'first_columns':header[:7],'target_rows':[x[:2] for x in found],'samples':len(header)-2})
    pd.DataFrame(found,columns=header).to_csv(OUT/(fname+'.targets.csv'),index=False)
    pd.DataFrame(symbols,columns=header[:2]).to_csv(OUT/(fname+'.genes.csv'),index=False)
    pd.Series(header[2:]).to_csv(OUT/(fname+'.ids.csv'),index=False)
    f.close()
for p in (ROOT/'Processed/protein').glob('*.xlsx'):
    w=openpyxl.load_workbook(p,read_only=True,data_only=True)
    emit({'workbook':p.name,'guide':[r for r in w['Guide'].iter_rows(values_only=True)]})
    w.close()
