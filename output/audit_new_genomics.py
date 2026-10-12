from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];P=ROOT/'Processed'
OUT=ROOT/'output/multiomics_audit/update_20260922';OUT.mkdir(exist_ok=True)
m=pd.read_csv(ROOT/'output/multiomics_audit/master_availability_and_targets.csv',index_col='depmap_id')
inventory=[]
for name,prefix in [('OmicsAbsoluteCNGene.csv','AbsoluteCN'),('OmicsLoH.csv','LoH'),('OmicsCNGeneWGS.csv','RelativeCN_WGS')]:
    h=pd.read_csv(P/name,nrows=0).columns
    targets=[c for c in h if c.startswith(('APC (','TCF7L2 (','CTNNB1 ('))]
    meta=[c for c in h if c in ['Unnamed: 0','ModelID','SequencingID','ModelConditionID','IsDefaultEntryForModel','IsDefaultEntryForMC']]
    d=pd.read_csv(P/name,usecols=meta+targets)
    original_n=len(d)
    if 'ModelID' not in d:d=d.rename(columns={h[0]:'ModelID'})
    if 'IsDefaultEntryForModel' in d:d=d.loc[d.IsDefaultEntryForModel.eq('Yes')]
    assert d.ModelID.is_unique,(name,'duplicate models')
    d=d.set_index('ModelID');d.to_csv(OUT/(prefix+'_target_genes.csv'))
    for c in targets:m[prefix+'_'+c.split(' (')[0]]=d[c]
    inventory.append(dict(file=name,rows=original_n,default_models=len(d),cohort_matched=int(d.index.isin(m.index).sum()),targets=targets))
profiles=pd.read_csv(P/'OmicsProfiles.csv')
profiles.to_csv(OUT/'profiles_metadata_snapshot.csv',index=False)
seg=pd.read_csv(P/'OmicsAbsoluteCNSegmentsProfile.csv',dtype={'Chromosome':str})
profile_ids=set(seg.ProfileID)
hit=profiles.loc[profiles.ProfileID.isin(profile_ids)].copy()
print('ABS_SEG_PROFILE_MAPPING',{'segments':len(seg),'profiles':len(profile_ids),'matched_profiles':hit.ProfileID.nunique(),'matched_models':hit.ModelID.nunique(),'duplicate_ProfileID':int(hit.ProfileID.duplicated().sum()),'default_values':hit.IsDefaultEntryForModel.value_counts(dropna=False).to_dict(),'LoHStatus':seg.LoHStatus.value_counts(dropna=False).to_dict()},flush=True)
pd.Series(sorted(profile_ids-set(hit.ProfileID))).to_csv(OUT/'unmapped_absolute_segment_profiles.csv',index=False)
hit.to_csv(OUT/'absolute_segment_profile_mapping.csv',index=False)
# Limit chromosome-5 segments for downstream locus-specific review; no unverified gene coordinates are used.
seg.loc[seg.Chromosome.isin(['5','chr5'])].to_csv(OUT/'absolute_segments_chr5.csv',index=False)
inventory.append(dict(file='OmicsAbsoluteCNSegmentsProfile.csv',rows=len(seg),profiles=len(profile_ids),mapped_profiles=hit.ProfileID.nunique(),mapped_models=hit.ModelID.nunique(),cohort_matched=hit.loc[hit.ModelID.isin(m.index)].ModelID.nunique()))
global_=pd.read_csv(P/'OmicsGlobalSignatures.csv')
gd=global_.loc[global_.IsDefaultEntryForModel.eq('Yes')].set_index('ModelID')
assert gd.index.is_unique
for c in ['MSIScore','LoHFraction','WGD','CIN','Ploidy','Aneuploidy']:m['New_'+c]=gd[c]
inventory.append(dict(file='OmicsGlobalSignatures.csv',rows=len(global_),default_models=len(gd),cohort_matched=int(gd.index.isin(m.index).sum())))
sub=pd.read_csv(P/'OmicsInferredMolecularSubtypes.csv').set_index('ModelID')
assert sub.index.is_unique
for c in ['APC_LoF','MSI','KRAS p.G12','KRAS p.G13','KRAS p.Q61','BRAF p.V600E']:
    if c in sub:m['New_subtype_'+c]=sub[c]
inventory.append(dict(file='OmicsInferredMolecularSubtypes.csv',rows=len(sub),default_models=len(sub),cohort_matched=int(sub.index.isin(m.index).sum())))
fields=[c for c in m if c.startswith(('AbsoluteCN_','LoH_','RelativeCN_','New_'))]
rows=[]
for scope,df in [('Pan-cancer',m),('CRC',m.loc[m.lineage_2.eq('Colorectal Adenocarcinoma')])]:
    for c in fields:
        v=df.loc[df[c].notna()];ct=v.APC_class.value_counts()
        rows.append(dict(scope=scope,field=c,n=len(v),WT=int(ct.get('WT',0)),Missense=int(ct.get('Missense only',0)),LOF=int(ct.get('Truncating/LOF',0)),Other=int(ct.get('Other mutation',0))))
coverage=pd.DataFrame(rows);coverage.to_csv(OUT/'new_data_coverage.csv',index=False)
pd.DataFrame(inventory).to_csv(OUT/'new_file_inventory.csv',index=False)
joint=m[['AbsoluteCN_APC','LoH_APC']].notna().all(axis=1)
m['APC_CN_LOH_covered']=joint
# Observed-data descriptions, not inferred mono/biallelic loss labels.
m['APC_LOF_with_LoH_code_1']=m.APC_class.eq('Truncating/LOF')&m.LoH_APC.eq(1)
m['APC_LOF_with_LoH_code_0']=m.APC_class.eq('Truncating/LOF')&m.LoH_APC.eq(0)
m.to_csv(OUT/'master_with_new_genomics.csv')
cross=pd.crosstab(m.APC_class,m.LoH_APC,dropna=False);cross.to_csv(OUT/'APC_class_by_LoH_code.csv')
compare=pd.crosstab(m.APC_class,m.New_subtype_APC_LoF,dropna=False);compare.to_csv(OUT/'APC_class_by_new_inferred_LoF.csv')
summary={'joint_APC_CN_LOH':int(joint.sum()),'joint_CRC':int((joint&m.lineage_2.eq('Colorectal Adenocarcinoma')).sum()),'APC_absolute_values':m.AbsoluteCN_APC.value_counts(dropna=False).sort_index().to_dict(),'APC_LOH_codes':m.LoH_APC.value_counts(dropna=False).to_dict(),'segment_profiles_total':len(profile_ids),'segment_profiles_mapped':hit.ProfileID.nunique()}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2,default=str))
print(pd.DataFrame(inventory).to_string(index=False),flush=True)
print(coverage.loc[coverage.field.isin(['AbsoluteCN_APC','LoH_APC','RelativeCN_WGS_APC','New_MSIScore','New_Ploidy','New_subtype_APC_LoF','New_subtype_MSI'])].to_string(index=False),flush=True)
print('APC_class_by_LoH',cross.to_string(),flush=True)
print('APC_class_by_new_subtype',compare.to_string(),flush=True)
print(json.dumps(summary,indent=2,default=str),flush=True)
