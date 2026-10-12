from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];P=ROOT/'Processed';O=ROOT/'output/multiomics_audit_26Q1';O.mkdir(exist_ok=True)
old=pd.read_csv(ROOT/'output/multiomics_audit/master_availability_and_targets.csv',index_col=0)
master=old.copy();inventory=[];extract={}
targets=['APC','TCF7L2','CTNNB1','AXIN2','MYC','CCND1','LGR5','ASCL2','RNF43','SP5','KRAS','BRAF','TP53']
for name in ['OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv','OmicsAbsoluteCNGene.csv','OmicsLoH.csv','OmicsCNGeneWGS.csv','OmicsSomaticMutationsMatrixDamaging.csv','ScreenGeneEffect.csv']:
    cols=pd.read_csv(P/name,nrows=0).columns
    meta=[c for c in cols[:6] if '(' not in c]
    gc=[c for c in cols if any(c.startswith(g+' (') for g in targets)]
    df=pd.read_csv(P/name,usecols=list(dict.fromkeys(meta+gc)))
    if 'IsDefaultEntryForModel' in df:
        dd=df.loc[df.IsDefaultEntryForModel.eq('Yes')].copy();idcol='ModelID'
    else:dd=df.copy();idcol=cols[0]
    dup=int(dd[idcol].duplicated().sum())
    inventory.append({'file':name,'rows':len(df),'feature_columns':len(cols)-len(meta),'default_rows':len(dd),'unique_ids':dd[idcol].nunique(),'duplicate_default_ids':dup,'overlap_old':int(dd[idcol].isin(old.index).sum())})
    dd.to_csv(O/(name+'.targets.csv'),index=False);extract[name]=dd
    print(inventory[-1],flush=True)
    if dup:continue
    dd=dd.set_index(idcol)
    if name.startswith('OmicsExpression'):
        for c in gc:master['RNA26_'+c.split(' (')[0]]=dd[c]
    elif name=='OmicsAbsoluteCNGene.csv':master['APC_absolute_CN']=dd['APC (324)']
    elif name=='OmicsLoH.csv':master['APC_LoH_code']=dd['APC (324)']
    elif name=='OmicsCNGeneWGS.csv':master['APC_relative_CN_WGS']=dd['APC (324)']
    elif name=='OmicsSomaticMutationsMatrixDamaging.csv':master['APC_damaging_matrix_code']=dd['APC (324)']

ms=[];covered=set()
uc=['SequencingID','ModelID','IsDefaultEntryForModel','Chrom','Pos','Ref','Alt','AF','DP','GT','PS','VariantInfo','ProteinChange','HugoSymbol','LikelyLoF','VepImpact','VepBiotype','EnsemblFeatureID']
for ch in pd.read_csv(P/'OmicsSomaticMutations.csv',usecols=uc,chunksize=250000):
    ch=ch.loc[ch.IsDefaultEntryForModel.eq('Yes')];covered.update(ch.ModelID)
    ms.append(ch.loc[ch.HugoSymbol.eq('APC')])
apc=pd.concat(ms).drop_duplicates();apc['is_LOF']=apc.LikelyLoF.eq(True)|apc.VariantInfo.str.contains('stop_gained|frameshift_variant',na=False)
# Match prior coding-only grouping convention; do not turn intronic variants into coding MUT.
apc['coding_consequence']=apc.VariantInfo.str.contains('missense_variant|stop_gained|frameshift_variant|inframe_deletion|inframe_insertion|start_lost|stop_lost|splice_acceptor_variant|splice_donor_variant',na=False)
apc.to_csv(O/'APC_default_variants_26Q1.csv',index=False)
lof=set(apc.loc[apc.is_LOF,'ModelID']);coding=set(apc.loc[apc.coding_consequence,'ModelID'])
mis=set(apc.loc[apc.VariantInfo.str.contains('missense_variant',na=False),'ModelID'])-lof
other=set(apc.loc[apc.coding_consequence&~apc.VariantInfo.str.contains('missense_variant',na=False),'ModelID'])-lof
mis-=other
master['APC_class_26Q1_provisional']=[('Unavailable' if i not in covered else 'Truncating/LOF' if i in lof else 'Missense only' if i in mis else 'Other coding' if i in coding else 'WT') for i in master.index]
changes=master.loc[master.APC_class.ne(master.APC_class_26Q1_provisional),['cell_line_display_name','APC_class','APC_class_26Q1_provisional']]
changes.to_csv(O/'APC_class_changes_review.csv')
sig=pd.read_csv(P/'OmicsGlobalSignatures.csv');sig=sig.loc[sig.IsDefaultEntryForModel.eq('Yes')]
print('signature duplicates',sig.ModelID.duplicated().sum(),flush=True)
if not sig.ModelID.duplicated().any():
    for c in ['MSIScore','Ploidy','LoHFraction','WGD']:master['26Q1_'+c]=sig.set_index('ModelID')[c]

# Profile-level allele-specific copy number: check join coverage, not automatic biallelic calls.
pr=pd.read_csv(P/'OmicsProfiles.csv');seg=pd.read_csv(P/'OmicsAbsoluteCNSegmentsProfile.csv')
segids=pd.Index(seg.ProfileID.unique());matches=pr.loc[pr.ProfileID.isin(segids)]
print('absolute segment profiles',len(segids),'mapped',matches.ProfileID.nunique(),'default profiles',len(matches.loc[matches.IsDefaultEntryForModel.eq('Yes')]),flush=True)
matches.to_csv(O/'absolute_CN_profile_mapping.csv',index=False)
# Match chr5 segments to observed APC variant positions, which are hg38 per README.
chr5=seg.loc[seg.Chromosome.astype(str).isin(['5','chr5'])]
mapped=chr5.merge(pr[['ProfileID','ModelID','IsDefaultEntryForModel','DataType']],on='ProfileID',how='left')
default=mapped.loc[mapped.IsDefaultEntryForModel.eq('Yes')]
apcsegments=[]
for _,v in apc.loc[apc.ModelID.isin(old.index)&apc.is_LOF].iterrows():
    hits=default.loc[default.ModelID.eq(v.ModelID)&default.Start.le(v.Pos)&default.End.ge(v.Pos)]
    for _,r in hits.iterrows():apcsegments.append({'ModelID':v.ModelID,'mutation_sequencing_id':v.SequencingID,'Pos':v.Pos,'ProteinChange':v.ProteinChange,'AF':v.AF,**r.to_dict()})
hits=pd.DataFrame(apcsegments);hits.to_csv(O/'APC_LOF_variant_position_CN_evidence.csv',index=False)

mc=[c for c in master if c.startswith('METH_')]
master['RNA26_protein']=master.RNA26_TCF7L2.notna()&master.MS_TCF7L2.notna()
master['Full26_all3']=master.RNA26_protein&master[mc].notna().all(axis=1)
for c in mc:master['Full26_'+c]=master.RNA26_protein&master[c].notna()
features=['RNA_RPKM_TCF7L2','RNA26_TCF7L2','RNA_Protein','RNA26_protein','RNA_Protein_Meth_all3','Full26_all3']+['Full26_'+c for c in mc]+['APC_absolute_CN','APC_LoH_code','APC_relative_CN_WGS','APC_damaging_matrix_code','26Q1_MSIScore']
rows=[]
for groupcol in ['APC_class','APC_class_26Q1_provisional']:
    for scope,dd in [('Pan-cancer',master),('CRC',master.loc[master.lineage_2.eq('Colorectal Adenocarcinoma')])]:
        for c in features:
            mask=dd[c] if dd[c].dtype==bool else dd[c].notna();v=dd.loc[mask];ct=v[groupcol].value_counts()
            rows.append({'group_definition':groupcol,'scope':scope,'feature':c,'n_total':len(v),'WT':int(ct.get('WT',0)),'LOF':int(ct.get('Truncating/LOF',0)),'Missense':int(ct.get('Missense only',0))})
coverage=pd.DataFrame(rows);coverage.to_csv(O/'coverage_comparison.csv',index=False)
master.to_csv(O/'master_audit_26Q1.csv');pd.DataFrame(inventory).to_csv(O/'new_file_inventory.csv',index=False)
model=pd.read_csv(P/'Model.csv');ls=model.loc[model.StrippedCellLineName.str.upper().eq('LS174T')]
summary={'old_class_counts':master.APC_class.value_counts().to_dict(),'new_provisional_class_counts':master.APC_class_26Q1_provisional.value_counts().to_dict(),'changed_labels':len(changes),'CN_distribution':master.APC_absolute_CN.value_counts().to_dict(),'LoH_code_distribution':master.APC_LoH_code.value_counts().to_dict(),'segments_fields':seg.columns.tolist(),'absolute_segment_profiles':len(segids),'profile_mapping_matches':matches.ProfileID.nunique(),'apc_LOF_variant_CN_models':hits.ModelID.nunique() if len(hits) else 0,'LS174T':ls[['ModelID','CellLineName']].to_dict('records')}
(O/'summary.json').write_text(json.dumps(summary,indent=2))
print(coverage.loc[coverage.group_definition.eq('APC_class')&coverage.scope.eq('Pan-cancer')].to_string(index=False),flush=True);print(json.dumps(summary,indent=2),flush=True)
