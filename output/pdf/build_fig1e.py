from pathlib import Path
import json
import numpy as np
import pandas as pd
import pyreadr
import statsmodels.api as sm
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
import pymupdf

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
ph = next(iter(pyreadr.read_r(str(ROOT/'Processed/ph/cell_line.ph.rds')).values()))
expr = next(iter(pyreadr.read_r(str(ROOT/'Processed/expr/crispr.exp.rds')).values()))
covered, mutant = set(), set()
for chunk in pd.read_csv(ROOT/'Processed/mut/mutations_NAsdropped.csv', usecols=['gene','depmap_id'], chunksize=200000):
    covered.update(chunk.depmap_id.dropna())
    mutant.update(chunk.loc[chunk.gene.eq('APC'),'depmap_id'])
d = ph[['depmap_id','cell_line_display_name','lineage_1','lineage_2']].copy().reset_index(drop=True)
d['TCF7L2_gene_effect'] = d.depmap_id.map(expr.loc['TCF7L2'])
d['mutation_table_covered'] = d.depmap_id.isin(covered)
d['APC_MUT'] = d.depmap_id.isin(mutant).astype(int)
d['included'] = d.mutation_table_covered & d.lineage_1.notna() & d.TCF7L2_gene_effect.notna()
d.to_csv(OUT/'Fig1E_source_data.csv',index=False)
a=d.loc[d.included].copy()
assert a.depmap_id.is_unique
lineage = pd.get_dummies(a.lineage_1.astype(str),prefix='lineage',drop_first=True,dtype=float)
X0 = sm.add_constant(a[['APC_MUT']].astype(float))
X1 = pd.concat([X0,lineage],axis=1)
y = a.TCF7L2_gene_effect.astype(float)
results=[]
for label,X in [('Unadjusted',X0),('Lineage-adjusted',X1)]:
    assert np.linalg.matrix_rank(X)==X.shape[1]
    fit=sm.OLS(y,X).fit(cov_type='HC1',use_t=True)
    b=float(fit.params['APC_MUT']);lo,hi=fit.conf_int().loc['APC_MUT']
    results.append(dict(model=label,beta=b,ci_low=float(lo),ci_high=float(hi),p=float(fit.pvalues['APC_MUT']),se=float(fit.bse['APC_MUT']),n=int(fit.nobs),rank=X.shape[1],residual_df=int(fit.df_resid),r_squared=float(fit.rsquared)))
    pd.DataFrame({'term':fit.params.index,'beta':fit.params.values,'SE_HC1':fit.bse.values,'CI_low':fit.conf_int()[0].values,'CI_high':fit.conf_int()[1].values,'P':fit.pvalues.values}).to_csv(OUT/f'Fig1E_{label.lower()}_coefficients.csv',index=False)
    (OUT/f'Fig1E_{label.lower()}_model.txt').write_text(str(fit.summary()),encoding='utf-8')
# Independently verify adjusted APC coefficient by within-lineage demeaning.
xd=a.APC_MUT-a.groupby('lineage_1',observed=True).APC_MUT.transform('mean')
yd=y-a.groupby('lineage_1',observed=True).TCF7L2_gene_effect.transform('mean')
assert np.isclose(results[1]['beta'],float(xd@yd/(xd@xd)),atol=1e-10)
counts=a.groupby('lineage_1',observed=True).APC_MUT.agg(n='size',MUT='sum')
counts['WT']=counts.n-counts.MUT
counts.to_csv(OUT/'Fig1E_lineage_counts.csv')
pd.DataFrame(results).to_csv(OUT/'Fig1E_statistics.csv',index=False)
n=len(a);nm=int(a.APC_MUT.sum());nl=len(counts)

pdf=OUT/'Fig1E_Lineage_adjusted_multivariable_analysis.pdf'
c=canvas.Canvas(str(pdf),pagesize=(500,268))
c.setTitle('Figure 1E | Lineage-adjusted multivariable analysis')
ink=HexColor('#24282E');grey=HexColor('#66717D');red=HexColor('#C55358');blue=HexColor('#377DA8')
def text(x,y,s,size=8,font='Helvetica',color=ink,align='left'):
    c.setFont(font,size);c.setFillColor(color)
    {'left':c.drawString,'right':c.drawRightString,'center':c.drawCentredString}[align](x,y,s)
def beta_text(x,y,size=9):
    text(x,y,'\u03b2',size,'Symbol')
text(14,242,'E',16,'Helvetica-Bold')
text(38,244,'Lineage-adjusted multivariable analysis',12,'Helvetica-Bold')
text(38,228,f'Pan-cancer | n = {n:,} | {nl} tissue lineages | APC WT {n-nm:,} / MUT {nm}',8.5,color=grey)
text(18,196,'APC effect (MUT vs WT)',8.5,'Helvetica-Bold')
beta_text(352,196,8.5);text(360,196,' [95% CI]',8.5,'Helvetica-Bold')
text(477,196,'P value',8.5,'Helvetica-Bold',align='right')
left,right=169,331
low=min(r['ci_low'] for r in results);high=max(r['ci_high'] for r in results)
xmin=min(-.1,np.floor((low-.015)/.1)*.1);xmax=max(.1,np.ceil((high+.015)/.1)*.1)
def x(v):return left+(v-xmin)/(xmax-xmin)*(right-left)
c.setFillColor(HexColor('#F6F0F1'));c.roundRect(12,117,476,40,3,fill=1,stroke=0)
c.setStrokeColor(HexColor('#AAB3BC'));c.setLineWidth(.7);c.setDash(2,3);c.line(x(0),114,x(0),184);c.setDash()
for r,yy,col in zip(results,[171,136],[blue,red]):
    text(18,yy-3,r['model'],9,'Helvetica-Bold' if yy==136 else 'Helvetica')
    c.setStrokeColor(col);c.setFillColor(col);c.setLineWidth(1.4)
    c.line(x(r['ci_low']),yy,x(r['ci_high']),yy)
    c.line(x(r['ci_low']),yy-3,x(r['ci_low']),yy+3)
    c.line(x(r['ci_high']),yy-3,x(r['ci_high']),yy+3)
    c.circle(x(r['beta']),yy,3.4,stroke=0,fill=1)
    text(346,yy-3,f"{r['beta']:.3f} [{r['ci_low']:.3f}, {r['ci_high']:.3f}]",7.7)
    p=f"{r['p']:.3g}" if r['p']>=.001 else f"{r['p']:.1e}"
    text(477,yy-3,p,8,align='right')
c.setStrokeColor(ink);c.setLineWidth(.65);c.line(left,108,right,108)
for v in np.arange(xmin,xmax+.001,.1):
    c.line(x(v),108,x(v),105);text(x(v),94,f'{v:.1f}'.replace('-0.0','0.0'),7.5,align='center')
text((left+right)/2,79,'APC regression coefficient',8.5,align='center')
beta_text(305,79,8.5)
text((left+right)/2,65,'Negative: stronger dependency in APC MUT',7,color=grey,align='center')
text(18,42,'Adjusted model: TCF7L2 gene effect ~ APC status + tissue lineage',8)
text(18,28,'Points: OLS coefficients; bars: 95% CI (HC1 robust SE); two-sided P values.',7.5,color=grey)
text(18,15,'APC MUT: recorded coding mutation; WT: no recorded APC mutation in covered samples.',7,color=grey)
c.showPage();c.save()
doc=pymupdf.open(pdf)
doc[0].get_pixmap(matrix=pymupdf.Matrix(3,3)).save(OUT/'Fig1E_preview.png')
assert len(doc)==1 and not doc[0].get_images()
assert 'Lineage-adjusted' in doc[0].get_text()
assert all(b[0]>=0 and b[1]>=0 and b[2]<=500 and b[3]<=268 for b in doc[0].get_text('blocks'))
r=results[1]
notes=f'''# Figure 1E: lineage-adjusted multivariable analysis

Sources: Processed/ph/cell_line.ph.rds (lineage_1), Processed/expr/crispr.exp.rds (TCF7L2 row), Processed/mut/mutations_NAsdropped.csv (gene and depmap_id). Join by unique depmap_id. No score transformation; release and score algorithm are not documented in these source files and are not inferred.

Same APC definition as panel D: any APC record in the supplied processed coding mutation table gives MUT; absence gives WT only if the sample appears in the mutation table. WT therefore means no recorded APC coding mutation, not independently confirmed wild type. Copy-number alterations are not included. {len(d)} source cell lines; {n} complete cases; {len(d)-n} excluded for missing mutation-table coverage, lineage or score. {nm} MUT and {n-nm} WT, {nl} tissue lineages. Non-cancerous annotated models present in the provided 1,150-line discovery cohort are retained; this reproduces the supplied discovery cohort rather than silently narrowing it.

Primary model: TCF7L2_gene_effect ~ intercept + APC_MUT + categorical lineage_1. APC coefficient is MUT minus WT adjusted for tissue lineage. Lineage uses {nl-1} dummy variables with alphabetically first observed lineage as reference. No interactions. All available lineages retained, including singleton categories. The plot shows only the APC coefficient; lineage coefficients are nuisance adjustment terms provided in the full coefficient table. The unadjusted model uses exactly the same rows and serves as context.

OLS with HC1 heteroskedasticity-robust sandwich covariance and Student t inference using residual degrees of freedom. HC1 rather than HC3 avoids division by zero in singleton-lineage observations (unit leverage). Robust intervals remain approximate with sparse groups. Two-sided tests; no FDR adjustment for this prespecified APC contrast. Full-rank design verified; the adjusted APC coefficient independently reproduced using within-lineage demeaned outcome and exposure. These are mean regression effects, whereas D reports a median difference and a rank test.

The full model summary may warn that the omnibus robust covariance of constraints has reduced rank because singleton categories have exactly fitted observations. The omnibus model F test is not used or interpreted here; the reported APC contrast has finite positive robust variance.

Primary adjusted APC beta = {r['beta']:.8f}, 95% CI [{r['ci_low']:.8f}, {r['ci_high']:.8f}], P = {r['p']:.8g}. The interval {'excludes' if r['ci_high']<0 or r['ci_low']>0 else 'includes'} zero. Negative values indicate stronger dependency in APC MUT. Adjustment addresses measured tissue-lineage differences only; it does not establish causality or eliminate unmeasured confounding. MSI and other mutations are not included in this requested basic lineage-adjustment model.

Consistency caveat from D remains: the supplied data yield 142 APC MUT among all 1,150 annotated lines using this definition, not the 35 described in Fig1.md. No original 35-sample list or filtering script is available. Align upstream mutation criteria before assembling a manuscript.

Output PDF uses editable Helvetica/Symbol text and native vector graphical elements with zero raster image objects. Figure source data, lineage counts, both full coefficient tables, model summaries, statistics and this reproducible script are saved alongside it.
'''
(OUT/'Fig1E_methods.md').write_text(notes,encoding='utf-8')
print(json.dumps({'models':results,'n_lineages':nl,'singleton_lineages':counts[counts.n.eq(1)].index.tolist(),'excluded':len(d)-n},indent=2))
