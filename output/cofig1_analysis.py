from pathlib import Path
import numpy as np,pandas as pd,pyreadr,json
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'output/coFig1_analysis';O.mkdir(exist_ok=True)
m=pd.read_csv(ROOT/'output/Fig2_updated/tables/master_analysis.csv',index_col=0)
e=next(iter(pyreadr.read_r(str(ROOT/'Processed/expr/crispr.exp.rds')).values()))
wtids=m.index[m.APC_class.eq('WT')];loids=m.index[m.APC_class.eq('Truncating/LOF')]
w=e[wtids].to_numpy(float);l=e[loids].to_numpy(float);nw=np.isfinite(w).sum(axis=1);nl=np.isfinite(l).sum(axis=1)
pv=np.full(len(e),np.nan)
for i in range(0,len(e),256):pv[i:i+256]=stats.mannwhitneyu(l[i:i+256],w[i:i+256],axis=1,nan_policy='omit',method='asymptotic',alternative='two-sided').pvalue
ok=(nw>=3)&(nl>=3)&np.isfinite(pv);pv[~ok]=np.nan;q=np.full(len(e),np.nan);q[ok]=multipletests(pv[ok],method='fdr_bh')[1]
s=pd.DataFrame({'gene':e.index,'n_WT':nw,'n_LOF':nl,'mean_WT':np.nanmean(w,axis=1),'mean_LOF':np.nanmean(l,axis=1),'mean_difference':np.nanmean(l,axis=1)-np.nanmean(w,axis=1),'P':pv,'FDR':q,'tested':ok})
s=s.sort_values(['mean_difference','FDR','gene']).reset_index(drop=True);s.to_csv(O/'genomewide_LOF_vs_WT.csv',index=False)
top=s.loc[s.FDR.lt(.05)&s.mean_difference.lt(0)].head(15).copy();top.to_csv(O/'top15_LOF_dependencies.csv',index=False)
genes=['TCF7L2','TCF7','TCF7L1','LEF1','BCL9','BCL9L','PYGO1','PYGO2'];family=s.set_index('gene').loc[genes].reset_index();family['FDR_family8']=multipletests(family.P,method='fdr_bh')[1]
rng=np.random.default_rng(20260923);cis=[]
for g in genes:
    a=e.loc[g,loids].dropna().to_numpy();b=e.loc[g,wtids].dropna().to_numpy();bs=[]
    for _ in range(10):bs.extend(rng.choice(a,(500,len(a))).mean(axis=1)-rng.choice(b,(500,len(b))).mean(axis=1))
    cis.append(np.quantile(bs,[.025,.975]))
family['CI_low']=[v[0] for v in cis];family['CI_high']=[v[1] for v in cis];family.to_csv(O/'family_LOF_vs_WT.csv',index=False)
d=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=['TCF7L2_gene_effect','lineage_1']).copy()
X=pd.concat([pd.Series(1.,index=d.index,name='const'),d.APC_class.eq('Truncating/LOF').astype(float).rename('LOF'),pd.get_dummies(d.lineage_1,drop_first=True,dtype=float)],axis=1)
f=sm.OLS(d.TCF7L2_gene_effect,X).fit(cov_type='HC1',use_t=True);ci=f.conf_int().loc['LOF']
adj={'n':len(d),'WT':int(d.APC_class.eq('WT').sum()),'LOF':int(d.APC_class.eq('Truncating/LOF').sum()),'beta':f.params.LOF,'SE':f.bse.LOF,'CI_low':ci[0],'CI_high':ci[1],'P':f.pvalues.LOF,'df':f.df_resid}
pd.DataFrame([adj]).to_csv(O/'lineage_adjusted_LOF_effect.csv',index=False)
c=d.loc[d.lineage_2.eq('Colorectal Adenocarcinoma')].sort_values(['TCF7L2_gene_effect','cell_line_display_name']).copy();c['rank']=np.arange(1,len(c)+1);c.to_csv(O/'CRC_dependency_landscape.csv')
v=[m.loc[m.APC_class.eq(g),'TCF7L2_gene_effect'].dropna() for g in ['WT','Missense only','Truncating/LOF']]
prs=[]
for i,j in [(0,1),(0,2),(1,2)]:prs.append({'comparison':f'{i} vs {j}','P':stats.mannwhitneyu(v[i],v[j],alternative='two-sided',method='asymptotic').pvalue})
pp=pd.DataFrame(prs);pp['P_Holm']=multipletests(pp.P,method='holm')[1];pp.to_csv(O/'dependency_APC_pairwise.csv',index=False)
m.to_csv(O/'master_classification.csv');summary={'WT':len(wtids),'LOF':len(loids),'genes_tested':int(ok.sum()),'significant':int((q<.05).sum()),'negative_hits':int((s.FDR.lt(.05)&s.mean_difference.lt(0)).sum()),'adjusted':adj,'CRC_n':len(c),'Kruskal_Wallis_P':stats.kruskal(*v).pvalue}
(O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
print(family[['gene','mean_difference','FDR_family8']].to_string(index=False),flush=True)
