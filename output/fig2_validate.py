from pathlib import Path
import json, importlib.metadata as im
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
import gseapy as gp

O=Path(__file__).resolve().parent/'Fig2_updated'; T=O/'tables'; R=O/'resources'
m=pd.read_csv(T/'master_analysis.csv',index_col=0)
def fit(a,cols):
    X=(a[cols]-a[cols].mean())/a[cols].std()
    X=pd.concat([pd.Series(1.,index=a.index,name='const'),a.LOF,X,pd.get_dummies(a.lineage_1,drop_first=True,dtype=float)],axis=1)
    return sm.OLS(a.TCF7L2_gene_effect,X).fit(cov_type='HC1',use_t=True)
cols=['RNA26_TCF7L2','MS_TCF7L2']+[c for c in m if c.startswith('METH_')]
cc=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=cols+['lineage_1','TCF7L2_gene_effect'])
f=fit(cc,cols);old=pd.read_csv(T/'APC_multivariable_models.csv');ref=old.loc[old['mode'].eq('Complete_case')&old.model.eq('Full model')].iloc[0]
assert np.isclose(f.params.LOF,ref.beta,atol=1e-12)
out=[]
for ix in cc.loc[cc.LOF.eq(1)].index:
    ff=fit(cc.drop(index=ix),cols); ci=ff.conf_int().loc['LOF']
    out.append(dict(excluded_ModelID=ix,beta=ff.params.LOF,P=ff.pvalues.LOF,CI_low=ci[0],CI_high=ci[1]))
pd.DataFrame(out).to_csv(T/'full_model_leave_one_LOF_out.csv',index=False)
# Verify FWL genome-wide regression against an independently assembled dummy-variable design.
rna=pd.read_pickle(O/'RNA_default.pkl'); rank=pd.read_csv(T/'RNA_association.csv');checks=[]
for gene in ['TCF7L2','CTNNB1','AXIN2','MYC']:
    a=pd.concat([rna[gene].rename('x'),m.TCF7L2_gene_effect.rename('y'),m.lineage_1.rename('lin')],axis=1).dropna()
    X=pd.concat([pd.Series(1.,index=a.index,name='const'),((a.x-a.x.mean())/a.x.std()).rename('x'),pd.get_dummies(a.lin,drop_first=True,dtype=float)],axis=1)
    ff=sm.OLS(a.y,X).fit(cov_type='HC1',use_t=True);rr=rank.loc[rank.gene.eq(gene)].iloc[0]
    assert np.isclose(ff.params.x,rr.beta_per_SD,atol=1e-10) and np.isclose(ff.bse.x,rr.SE,atol=1e-10)
    checks.append(dict(gene=gene,beta=ff.params.x,SE=ff.bse.x,verified=True))
pd.DataFrame(checks).to_csv(T/'independent_regression_checks.csv',index=False)
versions={p:im.version(p) for p in ['numpy','pandas','scipy','statsmodels','gseapy','matplotlib','pymupdf','pyreadr']}
(O/'software_versions.json').write_text(json.dumps(versions,indent=2))
print('Independent regression verification passed; LOF leave-one-out:',pd.DataFrame(out)[['beta','P']].agg(['min','max']).to_dict(),flush=True)
# The requested effect-size ranking complements the precision-weighted t-statistic ranking.
sets={}
for p in sorted(R.glob('*.gmt')):
    for l in p.read_text().splitlines():
        r=l.split('\t');sets[r[0]]=r[2:]
for label in ['RNA','Protein']:
    dest=T/f'{label}_GSEA_beta.csv'
    if dest.exists(): continue
    rr=pd.read_csv(T/f'{label}_association.csv');rr['rank_beta']=-rr.beta_per_SD
    rr=rr[['gene','rank_beta']].dropna().sort_values(['rank_beta','gene'],ascending=[False,True])
    r=gp.prerank(rnk=rr,gene_sets=sets,min_size=15,max_size=500,permutation_num=1000,weight=1,threads=4,seed=20260922,outdir=None,no_plot=True)
    r.res2d.to_csv(dest,index=False);print(label,'beta-ranked GSEA complete',flush=True)
