from pathlib import Path
import json,re,os,warnings
import numpy as np
import pandas as pd
import pyreadr,gseapy as gp
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1];O=ROOT/'output/Fig2_updated';T=O/'tables';T.mkdir(exist_ok=True,parents=True);R=O/'resources'
SEED=20260922
def say(x):print(x,flush=True)
def bh(v):
    a=np.asarray(v,float);out=np.full(len(a),np.nan);m=np.isfinite(a);out[m]=multipletests(a[m],method='fdr_bh')[1];return out
def save(df,name):df.to_csv(T/(name+'.csv'),index=False)
def readgmt(p):
    return {r[0]:r[2:] for l in p.read_text().splitlines() if len(r:=l.split('\t'))>=3}
def z(v):return (v-v.mean())/v.std(ddof=1)
m=pd.read_csv(ROOT/'output/multiomics_audit_26Q1/master_audit_26Q1.csv',index_col=0)
m['APC_class']=m.APC_class_26Q1_provisional;m['LOF']=m.APC_class.eq('Truncating/LOF').astype(float)
classes=['WT','Missense only','Truncating/LOF'];meth=[c for c in m if c.startswith('METH_')]
assert m.index.is_unique

# Full new RNA matrix; select one default biological profile per ModelID.
cache=O/'RNA_default.pkl'
if cache.exists():rna=pd.read_pickle(cache)
else:
    rd=pd.read_csv(ROOT/'Processed/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv')
    rd=rd.loc[rd.IsDefaultEntryForModel.eq('Yes')].set_index('ModelID')
    assert rd.index.is_unique
    cols=[c for c in rd if re.search(r' \(\d+\)$',c)]
    rna=rd[cols].rename(columns=lambda x:x.rsplit(' (',1)[0]).T.groupby(level=0).mean().T
    rna=rna.loc[rna.index.intersection(m.index)].astype(float);rna.to_pickle(cache)
say(f'RNA {rna.shape}')
hall=readgmt(R/'h.all.v2025.1.Hs.symbols.gmt');key='HALLMARK_WNT_BETA_CATENIN_SIGNALING'
scorefile=T/'hallmark_ssgsea_scores.csv'
if scorefile.exists():ss=pd.read_csv(scorefile,index_col=0)
else:
    ssres=gp.ssgsea(data=rna.T,gene_sets=hall,outdir=None,min_size=10,max_size=500,sample_norm_method='rank',permutation_num=0,no_plot=True,threads=4,seed=SEED)
    ss=ssres.res2d.pivot(index='Name',columns='Term',values='ES').astype(float);ss.to_csv(scorefile)
m['Wnt_Hallmark']=z(ss[key])
custom=['AXIN2','MYC','CCND1','LGR5','ASCL2','RNF43','SP5'];m['Wnt_custom']=rna[custom].apply(z).mean(axis=1)
scoremeta={'Hallmark_WNT_genes_total':len(hall[key]),'Hallmark_WNT_genes_measured':len(set(hall[key])&set(rna.columns)),'custom_genes':custom}

# PROGENy footprint model: top 100 targets by resource p-value, weighted gene-z score.
pro=R/'PROGENy_omnipath.tsv'
if pro.exists():
    q=pd.read_csv(pro,sep='\t');say(('PROGENy columns',q.columns.tolist()))
    if {'record_id','label','value','genesymbol'}.issubset(q.columns):
        q=q.pivot_table(index=['record_id','genesymbol'],columns='label',values='value',aggfunc='first').reset_index()
        pathway='pathway' if 'pathway' in q else 'source'
        pp=q.loc[q[pathway].astype(str).str.upper().eq('WNT')].copy()
        pp['p_value']=pd.to_numeric(pp['p_value']);pp['weight']=pd.to_numeric(pp['weight'])
        pp=pp.sort_values('p_value').drop_duplicates('genesymbol').head(100)
        pp=pp.loc[pp.genesymbol.isin(rna.columns)];pp.to_csv(R/'PROGENy_WNT_top100_used.csv',index=False)
        w=pp.set_index('genesymbol').weight;m['Wnt_PROGENy']=z(rna[w.index].apply(z).mul(w).sum(axis=1)/np.sqrt((w*w).sum()))
        scoremeta['PROGENy_WNT_targets_used']=len(w)
(T/'score_definitions.json').write_text(json.dumps(scoremeta,indent=2))

def model(df,y,continuous=(),apc=False,lineage=True):
    need=[y]+list(continuous)+(['lineage_1'] if lineage else [])
    a=df.dropna(subset=need).copy()
    if apc:a=a.loc[a.APC_class.isin(['WT','Truncating/LOF'])]
    X=pd.DataFrame({'const':np.ones(len(a))},index=a.index)
    if apc:X['LOF']=a.LOF
    for c in continuous:X[c]=z(a[c])
    if lineage:X=pd.concat([X,pd.get_dummies(a.lineage_1,drop_first=True,dtype=float)],axis=1)
    assert np.linalg.matrix_rank(X)==X.shape[1],(y,continuous,len(a))
    return sm.OLS(a[y].astype(float),X.astype(float)).fit(cov_type='HC1',use_t=True),a
def termrow(f,a,term,**kw):
    lo,hi=f.conf_int().loc[term]
    return dict(**kw,n=len(a),WT=int(a.APC_class.eq('WT').sum()),LOF=int(a.LOF.sum()),beta=f.params[term],SE=f.bse[term],CI_low=lo,CI_high=hi,P=f.pvalues[term],df=f.df_resid)

# Pairwise descriptive differences: independent bootstrap, MW and within-variable Holm.
rng=np.random.default_rng(SEED);pw=[];summ=[]
features=['Wnt_Hallmark','Wnt_custom']+(['Wnt_PROGENy'] if 'Wnt_PROGENy' in m else [])+['RNA26_TCF7L2','MS_TCF7L2']+meth+['MS_CTNNB1','RPPA_beta_catenin']
for feat in features:
    vals={g:m.loc[m.APC_class.eq(g),feat].dropna().to_numpy() for g in classes}
    for g,v in vals.items():summ.append(dict(feature=feat,group=g,n=len(v),median=np.median(v),mean=np.mean(v),sd=np.std(v,ddof=1)))
    start=len(pw)
    for ref,test in [('WT','Missense only'),('WT','Truncating/LOF'),('Missense only','Truncating/LOF')]:
        w,v=vals[ref],vals[test];p=stats.mannwhitneyu(v,w,alternative='two-sided',method='asymptotic').pvalue
        boots=[]
        for _ in range(10):boots.extend(np.median(rng.choice(v,(500,len(v))),axis=1)-np.median(rng.choice(w,(500,len(w))),axis=1))
        ci=np.quantile(boots,[.025,.975]);pw.append(dict(feature=feat,comparison=f'{test} - {ref}',n_reference=len(w),n_test=len(v),median_difference=np.median(v)-np.median(w),CI_low=ci[0],CI_high=ci[1],P=p))
    qs=multipletests([r['P'] for r in pw[start:]],method='holm')[1]
    for rr,qv in zip(pw[start:],qs):rr['P_Holm']=qv
pw=pd.DataFrame(pw);pw['FDR_all_feature_comparisons']=bh(pw.P);save(pw,'group_comparisons');save(pd.DataFrame(summ),'group_summary')

# Primary scatter associations in WT/LOF, all-class and CRC sensitivity.
assoc=[];groupmodels=[]
for feat in features:
    for scope,df in [('WT_LOF',m.loc[m.APC_class.isin(['WT','Truncating/LOF'])]),('All_classes',m.loc[m.APC_class.isin(classes)]),('CRC',m.loc[m.lineage_2.eq('Colorectal Adenocarcinoma')&m.APC_class.isin(['WT','Truncating/LOF'])])]:
        f,a=model(df,'TCF7L2_gene_effect',[feat],lineage=scope!='CRC')
        rho,p=stats.spearmanr(a[feat],a.TCF7L2_gene_effect)
        assoc.append(termrow(f,a,feat,feature=feat,scope=scope,rho=rho,Spearman_P=p))
    f,a=model(m,feat,apc=True);groupmodels.append(termrow(f,a,'LOF',feature=feat))
assoc=pd.DataFrame(assoc);assoc['FDR']=bh(assoc.P);assoc['Spearman_FDR']=bh(assoc.Spearman_P);save(assoc,'dependency_associations')
gm=pd.DataFrame(groupmodels);gm['FDR']=bh(gm.P);save(gm,'lineage_adjusted_group_effects')

# Six models on identical complete cases; all available and single-region sensitivities.
modelspec=[('Unadjusted',[],False),('+ Lineage',[],True),('+ RNA',['RNA26_TCF7L2'],True),('+ Protein',['MS_TCF7L2'],True),('+ Methylation',meth,True),('Full model',['RNA26_TCF7L2','MS_TCF7L2']+meth,True)]
fullvars=['RNA26_TCF7L2','MS_TCF7L2']+meth
cc=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=fullvars+['TCF7L2_gene_effect','lineage_1'])
fr=[];fits={}
for mode,df in [('Complete_case',cc),('All_available',m)]:
    for name,covs,lin in modelspec:
        f,a=model(df,'TCF7L2_gene_effect',covs,apc=True,lineage=lin);fr.append(termrow(f,a,'LOF',mode=mode,model=name))
        if mode=='Complete_case':fits[name]=f
for region in meth:
    a=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=['RNA26_TCF7L2','MS_TCF7L2',region])
    for name,covs,lin in modelspec:
        cv=[region if x==meth[0] else x for x in covs if x not in meth[1:]]
        f,dd=model(a,'TCF7L2_gene_effect',cv,apc=True,lineage=lin);fr.append(termrow(f,dd,'LOF',mode=region,model=name))
save(pd.DataFrame(fr),'APC_multivariable_models')
influence=fits['Full model'].get_influence();diag=pd.DataFrame({'ModelID':cc.index,'APC_class':cc.APC_class,'cooks_distance':influence.cooks_distance[0],'hat':influence.hat_matrix_diag});save(diag,'full_model_influence')
say(('Complete cases',len(cc),'LOF',int(cc.LOF.sum()),'full APC',fr[5]))

# Methylation correlations to RNA, protein and dependency, within each APC class too.
mr=[]
for reg in meth:
    for target in ['RNA26_TCF7L2','MS_TCF7L2','TCF7L2_gene_effect']:
        for group in ['WT_LOF','WT','Truncating/LOF']:
            df=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])] if group=='WT_LOF' else m.loc[m.APC_class.eq(group)]
            f,a=model(df,target,[reg]);rho,p=stats.spearmanr(a[reg],a[target]);mr.append(termrow(f,a,reg,region=reg,target=target,group=group,rho=rho,Spearman_P=p))
mr=pd.DataFrame(mr);mr['FDR']=bh(mr.P);mr['Spearman_FDR']=bh(mr.Spearman_P);save(mr,'methylation_associations')

# Exploratory mediation, lineage-stratified bootstrap; no causal interpretation.
meddf=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=['Wnt_Hallmark','TCF7L2_gene_effect']).copy()
def medcalc(df):
    cols=['LOF','Wnt_Hallmark','TCF7L2_gene_effect'];v=df[cols]-df.groupby('lineage_1')[cols].transform('mean');a0=v.LOF.to_numpy();mv=v.Wnt_Hallmark.to_numpy();yv=v.TCF7L2_gene_effect.to_numpy()
    aa=(a0@mv)/(a0@a0);bb=np.linalg.lstsq(np.column_stack([a0,mv]),yv,rcond=None)[0][1];total=a0@yv/(a0@a0);return aa,bb,aa*bb,total
ma,mb,ind,total=medcalc(meddf);parts=[g for _,g in meddf.groupby('lineage_1')];boots=[]
for _ in range(2000):boots.append(medcalc(pd.concat([g.iloc[rng.integers(0,len(g),len(g))] for g in parts]))[2])
ci=np.quantile(boots,[.025,.975]);save(pd.DataFrame([dict(n=len(meddf),a=ma,b=mb,indirect=ind,CI_low=ci[0],CI_high=ci[1],total=total,bootstrap=2000)]),'exploratory_mediation')

# Evidence-based LOF+LOH stratification; not mono/biallelic assignment.
ev=pd.read_csv(ROOT/'output/multiomics_audit_26Q1/APC_CN_mutation_sequencing_match.csv')
ev=ev.loc[ev.mutation_sequencing_id.eq(ev.SequencingID)]
minor=ev.groupby('ModelID').MinorAlleleAbsoluteCN.agg(['min','max'])
m['LOH_evidence']=pd.NA
m.loc[m.APC_class.eq('WT')&m.APC_absolute_CN.notna(),'LOH_evidence']='WT with CN'
m.loc[m.index.isin(minor.index)&m.APC_class.eq('Truncating/LOF'),'LOH_evidence']='LOF; minor CN > 0'
m.loc[m.index.isin(minor.index[minor['min'].eq(0)])&m.APC_class.eq('Truncating/LOF'),'LOH_evidence']='LOF; minor CN = 0'
evrows=[]
egs=['WT with CN','LOF; minor CN > 0','LOF; minor CN = 0']
for ref,test in [(egs[0],egs[1]),(egs[0],egs[2]),(egs[1],egs[2])]:
    dd=m.loc[m.LOH_evidence.isin([ref,test])].copy();dd['APC_class']=np.where(dd.LOH_evidence.eq(test),'Truncating/LOF','WT');dd['LOF']=dd.APC_class.eq('Truncating/LOF').astype(float)
    f,a=model(dd,'TCF7L2_gene_effect',apc=True);evrows.append(termrow(f,a,'LOF',comparison=f'{test} vs {ref}'))
evrows=pd.DataFrame(evrows);evrows['P_Holm']=multipletests(evrows.P,method='holm')[1];save(evrows,'LOH_evidence_comparisons')
m.to_csv(T/'master_analysis.csv');say('Targeted analyses complete')

# Genome-wide association, FWL OLS with HC1 robust SE and per-feature sample selection.
def wide_assoc(mat,kind):
    out=[]
    for j,g in enumerate(mat.columns):
        v=pd.concat([mat[g].rename('x'),m.TCF7L2_gene_effect.rename('y'),m.lineage_1.rename('lin')],axis=1).dropna()
        if len(v)<100 or v.x.std()<1e-8:continue
        v['x']=z(v.x);xr=v.x-v.groupby('lin').x.transform('mean');yr=v.y-v.groupby('lin').y.transform('mean');ssx=xr@xr
        if ssx<1e-9:continue
        beta=xr@yr/ssx;res=yr-beta*xr;df=len(v)-v.lin.nunique()-1
        var=len(v)/df*float(((xr*res)**2).sum())/ssx**2;se=np.sqrt(var);tv=beta/se;p=2*stats.t.sf(abs(tv),df)
        out.append(dict(gene=g,n=len(v),beta_per_SD=beta,SE=se,t=tv,P=p,rank_metric=-tv))
        if j%3000==0:say(f'{kind}: {j}/{len(mat.columns)}')
    out=pd.DataFrame(out);out['FDR']=bh(out.P);save(out,f'{kind}_association');return out
def cod():
    p=T/'codependency_ranking.csv'
    if p.exists():return pd.read_csv(p)
    ex=next(iter(pyreadr.read_r(str(ROOT/'Processed/expr/crispr.exp.rds')).values())).T
    out=[]
    for j,g in enumerate(ex.columns):
        if g=='TCF7L2':continue
        v=pd.concat([ex[g].rename('x'),m.TCF7L2_gene_effect.rename('y'),m.lineage_1.rename('lin')],axis=1).dropna()
        if len(v)<100 or v.x.std()<1e-8:continue
        rx=v.x-v.groupby('lin').x.transform('mean');ry=v.y-v.groupby('lin').y.transform('mean')
        rho=stats.spearmanr(rx,ry).statistic;df=len(v)-v.lin.nunique()-1
        pval=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-12)),df)
        out.append(dict(gene=g,n=len(v),rho=rho,P=pval,rank_metric=rho))
        if j%4000==0:say(f'Co-dependency: {j}/{len(ex.columns)}')
    out=pd.DataFrame(out);out['FDR']=bh(out.P);out=out.sort_values('rho',ascending=False);save(out,'codependency_ranking');return out
co=cod()
rafile=T/'RNA_association.csv'
if rafile.exists():ra=pd.read_csv(rafile)
else:
    # Detectable in >=10% of samples and variable, then standardized continuous regression.
    keep=(rna.gt(1).mean()>=.1)&(rna.std()>1e-6);ra=wide_assoc(rna.loc[:,keep],'RNA')
pafile=T/'Protein_association.csv'
if pafile.exists():pa=pd.read_csv(pafile)
else:
    pr=pd.read_csv(ROOT/'Processed/protein/Proteomics_subsetted.csv',index_col=0)
    pm=pd.read_csv(ROOT/'output/multiomics_audit/protein_sample_map.csv');pr.index=pm.depmap_id;pr=pr.groupby(level=0).mean()
    pr=pr.T.groupby(pr.columns.str.split(' ').str[0]).mean().T
    pa=wide_assoc(pr,'Protein')
# Hallmark + Reactome + KEGG MEDICUS; FDR across sets within each ranking.
sets={}
for path in sorted(R.glob('*.gmt')):sets.update(readgmt(path))
for label,rank in [('codependency',co),('RNA',ra),('Protein',pa)]:
    dest=T/f'{label}_GSEA.csv'
    if dest.exists():continue
    rr=rank[['gene','rank_metric']].dropna().sort_values(['rank_metric','gene'],ascending=[False,True])
    result=gp.prerank(rnk=rr,gene_sets=sets,min_size=15,max_size=500,permutation_num=1000,weight=1,threads=4,seed=SEED,outdir=None,no_plot=True)
    result.res2d.to_csv(dest,index=False);say(f'{label} GSEA complete')
say('All analyses complete')
