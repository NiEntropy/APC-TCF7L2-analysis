"""Recompute final main-figure statistics from packaged numerical inputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'recomputed'
OUT.mkdir(exist_ok=True)
def read(folder, name):
    return pd.read_csv(ROOT / 'data' / folder / name)
def check(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-12)

growth = read('Fig2', 'Fig2_growth_data.csv')
expected = read('Fig2', 'Fig2F_statistics.csv')
ps = []
for row in expected.itertuples():
    d = growth[(growth.cell_line == row.cell_line) & (growth.time_h == 96)]
    a = d[d.shRNA == row.shRNA].absorbance_575nm
    b = d[d.shRNA == 'NC'].absorbance_575nm
    ps.append(stats.ttest_ind(a, b, equal_var=False).pvalue)
check(ps, expected.p_raw)
order = np.argsort(ps)
adj = np.empty(len(ps))
adj[order] = np.minimum(1, np.maximum.accumulate(np.array(ps)[order] * np.arange(len(ps), 0, -1)))
check(adj, expected.p_Holm)
pd.DataFrame({'cell_line':expected.cell_line,'shRNA':expected.shRNA,'p_raw':ps,'p_Holm':adj}).to_csv(OUT/'Fig2F.csv',index=False)

d = read('Fig3', 'plotted_data.csv')
expected = read('Fig3', 'Fig3D_F_paired_statistics.csv')
ps = []
for row in expected.itertuples():
    v = d[(d.cell == row.cell) & (d.metric == row.metric)].pivot(index='record', columns='group', values='value')
    ps.append(stats.ttest_rel(v.shTCF7L2, v.shNC).pvalue)
check(ps, expected.p_unadjusted)
expected.assign(p_recomputed=ps).to_csv(OUT/'Fig3D_F.csv',index=False)

d = read('Fig4', 'Fig4C_individual_volumes.csv')
days = sorted(d.day.unique())
groups = ['shNC','shTCF7L2']
a = np.stack([d[d.group==g].pivot(index='animal_ID_as_recorded',columns='day',values='volume_mm3')[days].to_numpy() for g in groups])
n = a.shape[1]; nt = len(days)
m = a.mean(axis=1); grand = a.mean(); gm = m.mean(axis=1); tm = m.mean(axis=0)
ss = [n*np.sum((m-gm[:,None]-tm+grand)**2),n*2*np.sum((tm-grand)**2),n*nt*np.sum((gm-grand)**2)]
dfs = [nt-1,nt-1,1]; errdf = 2*nt*(n-1)
sse = np.sum((a-m[:,None,:])**2); mse = sse/errdf
anova = pd.DataFrame({'effect':['Interaction','Time','Group'],'SS':ss,'df':dfs})
anova['F'] = anova.SS/anova.df/mse
anova['pvalue'] = stats.f.sf(anova.F,anova.df,errdf)
expected = read('Fig4','Fig4C_ordinary_two_way_ANOVA.csv').iloc[:3]
check(anova.F,expected.F); check(anova.pvalue,expected.pvalue)
anova.to_csv(OUT/'Fig4C_ANOVA.csv',index=False)
raw = 2*stats.t.sf(abs((m[0]-m[1])/np.sqrt(mse*2/n)),errdf)
with np.errstate(divide='ignore'):
    sidak = -np.expm1(nt*np.log1p(-raw))
check(sidak,read('Fig4','Fig4C_Sidak_comparisons.csv').p_Sidak_10_days)
pd.DataFrame({'day':days,'p_raw':raw,'p_Sidak':sidak}).to_csv(OUT/'Fig4C_Sidak.csv',index=False)

w = read('Fig4','Fig4D_Prism_weights.csv')
pD = stats.ttest_ind(w.shNC_g,w.shTCF7L2_g,equal_var=True).pvalue
fields = read('Fig4','TUNEL_fields.csv')
means = fields.groupby(['group','tissue_id_as_recorded']).positive_area_percent.mean()
pF = stats.ttest_ind(means['shNC'],means['shTCF7L2'],equal_var=False).pvalue
check([pD,pF],read('Fig4','Fig4_endpoint_tests.csv').pvalue)
pd.DataFrame({'panel':['Fig4D','Fig4F'],'p_unadjusted':[pD,pF]}).to_csv(OUT/'Fig4_endpoints.csv',index=False)

d = read('Fig5','Fig5G_independent_samples.csv')
expected = read('Fig5','Fig5G_Student_t.csv')
ps = []
for gene in expected.gene:
    v = d[d.gene==gene]
    ps.append(stats.ttest_ind(v[v.group=='shNC'].expression,v[v.group=='shTCF7L2'].expression,equal_var=True).pvalue)
check(ps,expected.P)
expected.assign(p_recomputed=ps).to_csv(OUT/'Fig5G.csv',index=False)

counts = {}
sets = {}
for cell in ['HT29','SW620']:
    d = pd.read_csv(ROOT/'data'/'RNAseq'/f'{cell}_DEG_full.csv.gz')
    sig = d[d.padj.le(.05) & d.log2FoldChange.abs().ge(1)]
    sets[cell] = {direction:set(sig.loc[(sig.log2FoldChange>0) if direction=='up' else (sig.log2FoldChange<0),'gene_id']) for direction in ['up','down']}
    counts[cell] = {k:len(v) for k,v in sets[cell].items()}
assert counts == {'HT29':{'up':375,'down':195},'SW620':{'up':1257,'down':664}}, counts
for direction, n in [('up',99),('down',17)]:
    overlap = sets['HT29'][direction] & sets['SW620'][direction]
    assert len(overlap)==n
    counts['shared_'+direction] = n
(OUT/'DEG_counts.json').write_text(json.dumps(counts,indent=2))
print('PASS: final Fig2F, Fig3D/F, Fig4C/D/F, Fig5G statistics and DEG/overlap counts match archived outputs.')
