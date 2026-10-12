from pathlib import Path
from itertools import combinations
import json
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import mannwhitneyu,kruskal
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor,Color,white
import pymupdf

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
d=pd.read_csv(OUT/'Fig1E_source_data.csv')
cols=['gene','depmap_id','variant_info','likely_lof','protein_change']
a=pd.concat([x.loc[x.gene.eq('APC')] for x in pd.read_csv(ROOT/'Processed/mut/mutations_NAsdropped.csv',usecols=cols,chunksize=200000)],ignore_index=True).drop_duplicates()
a=a.loc[a.depmap_id.isin(d.depmap_id)].copy()
a['truncating']=a.variant_info.str.contains('stop_gained|frameshift_variant',regex=True)
a['LOF_record']=a.truncating | a.likely_lof.eq(True)
mut=set(a.depmap_id);lof=set(a.loc[a.LOF_record,'depmap_id'])
miss=set(a.loc[a.variant_info.eq('missense_variant'),'depmap_id'])-lof
miss={i for i in miss if a.loc[a.depmap_id.eq(i),'variant_info'].eq('missense_variant').all()}
d['APC_class']=d.depmap_id.map(lambda i:'Truncating/LOF' if i in lof else 'Missense only' if i in miss else 'Other mutation' if i in mut else 'WT')
d['APC_variants']=d.depmap_id.map(a.groupby('depmap_id').apply(lambda x:'; '.join(x.variant_info+':'+x.protein_change),include_groups=False)).fillna('')
d['H_included']=d.included & d.APC_class.ne('Other mutation') & d.TCF7L2_gene_effect.notna()
d.to_csv(OUT/'Fig1H_source_data.csv',index=False)
a.to_csv(OUT/'Fig1H_APC_variant_annotations.csv',index=False)
data=d.loc[d.H_included].copy()
groups=['WT','Missense only','Truncating/LOF']
vals={g:data.loc[data.APC_class.eq(g),'TCF7L2_gene_effect'].to_numpy() for g in groups}
assert sum(map(len,vals.values()))==len(data)
overall=kruskal(*vals.values())
pd.DataFrame([dict(group=g,n=len(v),median=np.median(v),q1=np.quantile(v,.25),q3=np.quantile(v,.75)) for g,v in vals.items()]).to_csv(OUT/'Fig1H_group_summary.csv',index=False)
rng=np.random.default_rng(20260921)
rows=[]
for ref,test in combinations(groups,2):
    w,m=vals[ref],vals[test]
    mw=mannwhitneyu(m,w,alternative='two-sided',method='asymptotic')
    boots=[]
    for _ in range(40):
        boots.extend(np.median(rng.choice(m,(500,len(m))),axis=1)-np.median(rng.choice(w,(500,len(w))),axis=1))
    lo,hi=np.quantile(boots,[.025,.975])
    rows.append(dict(comparison=f'{test} - {ref}',median_difference=np.median(m)-np.median(w),CI_low=lo,CI_high=hi,P=mw.pvalue,Cliffs_delta=2*mw.statistic/(len(m)*len(w))-1))
stats=pd.DataFrame(rows);stats['P_Holm']=multipletests(stats.P,method='holm')[1]
stats.to_csv(OUT/'Fig1H_pairwise_statistics.csv',index=False)
# Sensitivity analysis using the same lineage adjustment as E.
X=sm.add_constant(pd.concat([pd.DataFrame({'Missense':data.APC_class.eq('Missense only').astype(float),'LOF':data.APC_class.eq('Truncating/LOF').astype(float)},index=data.index),pd.get_dummies(data.lineage_1,prefix='lineage',drop_first=True,dtype=float)],axis=1))
assert np.linalg.matrix_rank(X)==X.shape[1]
fit=sm.OLS(data.TCF7L2_gene_effect,X).fit(cov_type='HC1',use_t=True)
adj=[]
for label,mcoef,lcoef in [('Missense - WT',1,0),('Truncating/LOF - WT',0,1),('Truncating/LOF - Missense',-1,1)]:
    contrast=np.zeros(X.shape[1]);contrast[X.columns.get_loc('Missense')]=mcoef;contrast[X.columns.get_loc('LOF')]=lcoef
    r=fit.t_test(contrast);ci=r.conf_int()[0]
    adj.append(dict(comparison=label,beta=float(np.asarray(r.effect).item()),CI_low=ci[0],CI_high=ci[1],P=float(np.asarray(r.pvalue).item())))
adj=pd.DataFrame(adj);adj['P_Holm']=multipletests(adj.P,method='holm')[1];adj.to_csv(OUT/'Fig1H_lineage_adjusted_statistics.csv',index=False)

pdf=OUT/'Fig1H_APC_LOF_TCF7L2_dependency.pdf'
c=canvas.Canvas(str(pdf),pagesize=(460,490));c.setTitle('Figure 1H | APC loss-of-function status and TCF7L2 dependency')
ink=HexColor('#24282E');grey=HexColor('#66717D');colors=[HexColor('#377DA8'),HexColor('#BE9348'),HexColor('#C55358')]
def text(x,y,s,size=8,font='Helvetica',color=ink,align='left'):
    c.setFont(font,size);c.setFillColor(color)
    {'left':c.drawString,'center':c.drawCentredString,'right':c.drawRightString}[align](x,y,s)
def pformat(p):return f'{p:.3f}' if p>=.001 else f'{p:.1e}'
text(14,463,'H',16,'Helvetica-Bold')
text(38,465,'APC loss-of-function status',12,'Helvetica-Bold')
text(38,449,f'TCF7L2 dependency | Pan-cancer (n = {len(data):,})',9,color=grey)
text(38,433,f'Kruskal-Wallis P = {pformat(overall.pvalue)}',8.5)
left,right,bottom,top=67,438,208,407
ymin=np.floor(data.TCF7L2_gene_effect.min()*2)/2;ymax=np.ceil(data.TCF7L2_gene_effect.max()*2)/2
def y(v):return bottom+(v-ymin)/(ymax-ymin)*(top-bottom)
c.setStrokeColor(ink);c.setLineWidth(.65);c.line(left,bottom,left,top);c.line(left,bottom,right,bottom)
for v in np.arange(ymin,ymax+.001,.5):
    c.line(left-3,y(v),left,y(v));text(left-7,y(v)-3,f'{v:.1f}',8,align='right')
c.setStrokeColor(HexColor('#CCD3D9'));c.setDash(2,3);c.line(left,y(0),right,y(0));c.setDash()
c.saveState();c.translate(21,(top+bottom)/2);c.rotate(90);text(0,0,'TCF7L2 dependency score',9,align='center');c.restoreState()
jitter=np.random.default_rng(42)
for xx,g,col in zip([129,250,371],groups,colors):
    v=vals[g];q1,med,q3=np.quantile(v,[.25,.5,.75]);iqr=q3-q1
    lo=v[v>=q1-1.5*iqr].min();hi=v[v<=q3+1.5*iqr].max()
    for value,j in zip(v,jitter.uniform(-30,30,len(v))):
        c.setStrokeColor(Color(col.red,col.green,col.blue,alpha=.40 if g=='WT' else .7));c.setFillColor(Color(col.red,col.green,col.blue,alpha=.18 if g=='WT' else .55));c.setLineWidth(.35);c.circle(xx+j,y(value),1.6 if g=='WT' else 2.2,stroke=1,fill=1)
    c.setStrokeColor(col);c.setLineWidth(1);c.line(xx,y(lo),xx,y(hi));c.line(xx-8,y(lo),xx+8,y(lo));c.line(xx-8,y(hi),xx+8,y(hi))
    c.setFillColor(Color(col.red,col.green,col.blue,alpha=.09));c.rect(xx-24,y(q1),48,y(q3)-y(q1),fill=1,stroke=1)
    c.setStrokeColor(ink);c.setLineWidth(1.6);c.line(xx-24,y(med),xx+24,y(med))
    text(xx,192,'APC WT' if g=='WT' else g,9,'Helvetica-Bold',col,'center')
    text(xx,178,f'n = {len(v):,}',8,color=grey,align='center')
text(252,161,'Lower score indicates stronger dependency',7.5,color=grey,align='center')
c.setFillColor(HexColor('#F4F6F8'));c.roundRect(14,62,432,86,4,stroke=0,fill=1)
text(24,132,'Comparison',8,'Helvetica-Bold');text(235,132,'Median diff. [95% CI]',8,'Helvetica-Bold');text(434,132,'P (Holm)',8,'Helvetica-Bold',align='right')
for i,r in stats.iterrows():
    yy=113-i*18
    text(24,yy,r.comparison,7.8)
    text(235,yy,f'{r.median_difference:.3f} [{r.CI_low:.3f}, {r.CI_high:.3f}]',7.8)
    text(434,yy,pformat(r.P_Holm),7.8,align='right')
text(18,45,'LOF takes priority over missense; 1 in-frame-deletion-only cell line excluded.',7.2,color=grey)
text(18,31,'Boxes: median and IQR; whiskers: 1.5 IQR; points: individual cell lines.',7.2,color=grey)
text(18,17,'Pairwise tests: two-sided Mann-Whitney; CI: 20,000 bootstrap resamples.',7.2,color=grey)
c.showPage();c.save()
doc=pymupdf.open(pdf);doc[0].get_pixmap(matrix=pymupdf.Matrix(3,3)).save(OUT/'Fig1H_preview.png')
assert not doc[0].get_images()
assert all(b[0]>=0 and b[1]>=0 and b[2]<=460 and b[3]<=490 for b in doc[0].get_text('blocks'))
notes='''# Figure 1H: APC loss-of-function status and TCF7L2 dependency

Data and cohort: Fig1E_source_data.csv contains the original 1,150 annotated cell lines, their TCF7L2 scores and coverage status from the local RDS/processed mutation files. All covered lines with finite score enter classification. APC records are read from Processed/mut/mutations_NAsdropped.csv; identical selected-field records are deduplicated. The same provided discovery cohort and mutation definition as D-F are retained (including non-cancerous annotated models); no release or score algorithm is inferred.

Mutually exclusive classification, in precedence order:
1. Truncating/LOF: at least one stop_gained or frameshift_variant consequence (including compound consequences), or likely_lof=True. Mixed missense/truncating lines go here. All LOF-classified lines in these data also have a truncating variant; the two duplicate likely_lof missense records belong to C84, which also has truncations.
2. Missense only: APC mutation records are exclusively missense_variant, without any LOF record.
3. WT: no APC mutation record, among lines covered by the mutation table. This is not independently confirmed wild type.
4. Other mutation: EN (ACH-000978), APC p.E2617del, inframe_deletion, likely_lof=False. Recorded in source data and excluded from the three-group comparison, never relabeled WT or missense.

LOF is a predicted sequence-annotation category, not a demonstration of biallelic loss, absence of protein, or functional knockout. Copy-number loss, allele-specific status and second hits are not evaluated.

Main panel: pan-cancer comparison of the three groups. Points are individual lines with deterministic jitter (seed 42); WT uses smaller translucent points for visibility at n=1,008. Boxes show quartiles and median; whiskers are observed extremes within 1.5 IQR. All observations, including outliers, are shown. Overall Kruskal-Wallis and three pairwise two-sided asymptotic Mann-Whitney U tests (tie and continuity corrected). Pairwise P values receive Holm correction across all three comparisons. Median differences are first group minus second group as labeled; percentile 95% bootstrap intervals use 20,000 independent within-group resamples (seed 20260921). Intervals are pointwise, not familywise-adjusted; rank tests do not specifically test median equality. Cliff's delta is saved in the pairwise statistics.

Sensitivity table: OLS TCF7L2 score ~ APC class + categorical lineage_1, WT as class reference, HC1 robust covariance with residual-df t inference. Three APC-class contrasts and separate Holm adjustment. Full-rank design checked. HC1 avoids HC3 division by zero in singleton-lineage observations, though robust inference is approximate for sparse groups. Adjusted regression effects are mean differences, not the main panel's median differences.

This is observational association, not proof of a dose-response or causality. Tissue composition may explain unadjusted differences. No ordered trend is assumed or imposed. As in D-F, the original 35/1,115 split in Fig1.md is not reproduced: the source data give 142 APC-mutant lines (73 LOF, 68 missense-only, 1 other) and 1,008 without a recorded APC mutation.

PDF uses editable text and vector marks with zero raster images. Reproduction code, cell-line and variant data, group summaries, pairwise and lineage-adjusted statistics accompany the figure.

'''
notes+=f'Overall Kruskal-Wallis H={overall.statistic:.8g}, P={overall.pvalue:.8g}.\n\nGroup summaries:\n'+pd.read_csv(OUT/'Fig1H_group_summary.csv').to_string(index=False)+'\n\nPairwise results:\n'+stats.to_string(index=False)+'\n\nLineage-adjusted sensitivity:\n'+adj.to_string(index=False)
(OUT/'Fig1H_methods.md').write_text(notes,encoding='utf-8')
print(pd.read_csv(OUT/'Fig1H_group_summary.csv').to_string(index=False));print(stats.to_string(index=False));print('Adjusted:');print(adj.to_string(index=False));print('Excluded:',d.loc[~d.H_included,['cell_line_display_name','APC_class']].to_string(index=False))
