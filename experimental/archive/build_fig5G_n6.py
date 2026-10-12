from pathlib import Path
import sys
ROOT=Path.cwd();sys.path.insert(0,str(ROOT/'tmp/pdfs/python_libs'))
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
import pymupdf as fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json,shutil,zipfile
SRC=ROOT/'output/pdf/Fig4_Fig5_S9_S12'
OUT=ROOT/'output/pdf/Fig5G_n6_20261008';OUT.mkdir(parents=True,exist_ok=True)
DATA=OUT/'data';DATA.mkdir(exist_ok=True)
import xml.etree.ElementTree as ET
source=next((ROOT/'2026 崔珣 合成致死课题原始数据'/'Figure 4').rglob('HT29 shTCF7L2 细胞 qpcr.pzfx'))
root=ET.parse(source).getroot();ns={'p':'htt<LOCAL_PATH_REMOVED>'}
genes=['CCL5','CXCL10','IFNB1','IL18','HLA-B','HLA-C','ULBP-2'];arrays={};rows=[];inputs=[]
for table in root.findall('p:Table',ns):
    gene=table.findtext('p:Title',namespaces=ns)
    if gene not in genes:continue
    arrays[gene]=[np.array([float(d.text) for d in col.findall('.//p:d',ns) if d.text]) for col in table.findall('p:YColumn',ns)]
for gene in genes:
    a,b=arrays[gene];assert len(a)==len(b)==6
    t=ttest_ind(a,b,equal_var=True)
    symbol='****' if t.pvalue<.0001 else '***' if t.pvalue<.001 else '**' if t.pvalue<.01 else '*' if t.pvalue<.05 else 'ns'
    rows.append(dict(gene=gene,n_NC=6,n_KD=6,mean_NC=a.mean(),mean_KD=b.mean(),t=t.statistic,df=t.df,P=t.pvalue,symbol=symbol))
    for group,values in zip(['shNC','shTCF7L2'],[a,b]):
        for i,value in enumerate(values):inputs.append(dict(gene=gene,group=group,sample=i+1,expression=value))
result=pd.DataFrame(rows)
assert abs(result.set_index('gene').loc['HLA-B','P']-.005120653089)<1e-10
result.to_csv(DATA/'Fig5G_Student_t.csv',index=False)
pd.DataFrame(inputs).to_csv(DATA/'Fig5G_independent_samples.csv',index=False)
plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,'xtick.major.size':2,'ytick.major.size':2,'axes.labelpad':3})
f=plt.figure(figsize=(210/25.4,270/25.4));ax=f.add_axes([.65,.17,.31,.17]);colors=['#777D83','#8A3048']
for j,gene in enumerate(genes):
    for gi,v in enumerate(arrays[gene]):
        x=j+(gi-.5)*.30;ax.scatter(x+np.linspace(-.065,.065,6),v,s=12,c=colors[gi],edgecolor='white',lw=.25)
        ax.errorbar(x,v.mean(),yerr=v.std(ddof=1)/np.sqrt(6),fmt='_',ms=8,c=colors[gi],capsize=2,lw=.7)
    y=max(max(v.max(),v.mean()+v.std(ddof=1)/np.sqrt(6)) for v in arrays[gene])+.27
    ax.plot([j-.18,j-.18,j+.18,j+.18],[y,y+.10,y+.10,y],c='black',lw=.55)
    ax.text(j,y+.14,result.symbol.iloc[j],ha='center',va='bottom',fontsize=7)
ax.set(xticks=range(7),xticklabels=genes,ylabel='Source-normalized mRNA expression',ylim=(0,7));ax.tick_params(axis='x',rotation=45,labelsize=6);ax.spines[['top','right']].set_visible(False)
f.text(.65,.105,'6 independent biological samples / group; mean +/- SEM.\nTwo-sided unpaired Student t tests; unadjusted P values.',fontsize=5.8,linespacing=1.35)
tmp=ROOT/'tmp/pdfs/Fig5G_overlay.pdf';f.savefig(tmp);plt.close(f)
old=fitz.open(SRC/'Fig5.pdf');bounds=old[0].rect
r=fitz.Rect(.59*bounds.width,(1-.35)*bounds.height,.99*bounds.width,(1-.08)*bounds.height)
doc=fitz.open();p=doc.new_page(width=bounds.width,height=bounds.height)
for clip in [fitz.Rect(0,0,bounds.width,r.y0),fitz.Rect(0,r.y1,bounds.width,bounds.height),fitz.Rect(0,r.y0,r.x0,r.y1),fitz.Rect(r.x1,r.y0,bounds.width,r.y1)]:
    p.show_pdf_page(clip,old,0,clip=clip)
overlay=fitz.open(tmp);p.show_pdf_page(r,overlay,0,clip=r)
doc.save(OUT/'Fig5.pdf',garbage=4,deflate=True)
p.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(str(OUT/'Fig5_preview.png'))
p.get_pixmap(matrix=fitz.Matrix(3,3),clip=r).save(str(OUT/'Fig5G_preview.png'))
(OUT/'Fig5.svg').write_text(p.get_svg_image(text_as_path=False),encoding='utf-8')
notes='Fig5G statistical revision\n\nUser confirmed 6 independent biological samples per group. Source: '+str(source)+'\nTwo-sided unpaired equal-variance Student t test on source-normalized expression, separately for 7 genes. P values are unadjusted, matching the supplied Prism analysis. No log transformation. Points show all 6 samples; mean +/- SEM.\nEarlier interpretation of paired wells is superseded by user confirmation. Source cached Ct discrepancies have not been resolved; this reproduces Prism normalized inputs, without recalculating Ct.\n\n'+result.to_string(index=False)+'\n\nLegend: * P<0.05; ** P<0.01; *** P<0.001; **** P<0.0001; ns P>=0.05.\n'
(OUT/'Analysis_notes.txt').write_text(notes,encoding='utf-8')
shutil.copy2(__file__,OUT/'reproduce.py')
print(result.to_string(index=False),flush=True)
