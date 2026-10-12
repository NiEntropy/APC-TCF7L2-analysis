from pathlib import Path
import json, textwrap, zipfile, shutil, hashlib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import pymupdf as fitz

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/pdf/coFig1_final'; OUT.mkdir(parents=True,exist_ok=True)
PAN=OUT/'individual_panels'; PAN.mkdir(exist_ok=True)
TMP=ROOT/'tmp/pdfs/coFig1';TMP.mkdir(parents=True,exist_ok=True)
T=ROOT/'output/Fig2_updated/tables'
m=pd.read_csv(T/'master_analysis.csv',index_col=0)
ass=pd.read_csv(T/'dependency_associations.csv'); fr=pd.read_csv(T/'APC_multivariable_models.csv')
mr=pd.read_csv(T/'methylation_associations.csv');pw=pd.read_csv(T/'group_comparisons.csv')
co=pd.read_csv(T/'codependency_ranking.csv');meth=[c for c in m if c.startswith('METH_')]
COL={'WT':'#377DA8','Missense only':'#BE9348','Truncating/LOF':'#C55358'}
INK='#24282E'; GREY='#66717D'; LIGHT='#CDD4DA';MM=72/25.4
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.titlesize':7,'axes.labelsize':6.5,'xtick.labelsize':6,'ytick.labelsize':6,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':2,'ytick.major.size':2,'xtick.major.pad':1.7,'ytick.major.pad':1.7,'axes.labelpad':1.5,'savefig.facecolor':'white'})
manifest=[];meta={}; bounds=[]
def p(v):return f'{v:.3f}' if v>=.001 else f'{v:.1e}'
def label(fig,x,y,s,size=6.5,weight='normal',ha='left',**kw):
    w,h=fig.get_size_inches()*25.4
    return fig.text(x/w,1-y/h,s,fontsize=size,fontweight=weight,ha=ha,va='top',color=INK,**kw)
def new(letter,title,w,h):
    fig=plt.figure(figsize=(w/25.4,h/25.4))
    label(fig,.6,.65,letter,9.5,'bold');label(fig,5,.95,title,7,'bold')
    return fig
def axmm(fig,x,y,w,h):
    fw,fh=fig.get_size_inches()*25.4
    return fig.add_axes([x/fw,1-(y+h)/fh,w/fw,h/fh])
def chart(fig,w,h,top=10,left=10,bottom=9,right=2):return axmm(fig,left,top,w-left-right,h-top-bottom)
def save(fig,key,w,h,info=''):
    # Check artists before PDF export as PDF text extraction cannot detect already-clipped text.
    for axis in fig.axes:
        for direction in ['x','y']:
            lim=getattr(axis,'get_'+direction+'lim')();lo,hi=sorted(lim)
            ticks=getattr(axis,'get_'+direction+'ticks')()
            if len(ticks) and any(t<lo or t>hi for t in ticks):
                getattr(axis,'set_'+direction+'ticks')([t for t in ticks if lo<=t<=hi])
                getattr(axis,'set_'+direction+'lim')(lim)
    fig.canvas.draw(); renderer=fig.canvas.get_renderer();bb=fig.bbox
    for tx in fig.findobj(matplotlib.text.Text):
        if not tx.get_visible() or not tx.get_text():continue
        r=tx.get_window_extent(renderer)
        if r.x0 < -1 or r.y0 < -1 or r.x1>bb.width+1 or r.y1>bb.height+1:
            bounds.append({'file':key,'text':tx.get_text(),'bbox':[r.x0,r.y0,r.x1,r.y1]})
    dest=PAN/(key+'.pdf');fig.savefig(dest);plt.close(fig)
    meta[key]={'size_mm':[w,h],'description':info}
    verify(dest,w,h);return key
def verify(dest,w,h):
    d=fitz.open(dest);pg=d[0]
    assert len(d)==1 and not pg.get_images(),dest
    assert abs(pg.rect.width/MM-w)<.01 and abs(pg.rect.height/MM-h)<.01
    assert pg.get_text().strip() and pg.get_drawings()
    pg.get_pixmap(matrix=fitz.Matrix(3,3)).save(TMP/(dest.stem+'.png'))
    manifest.append({'file':str(dest.relative_to(OUT)),'width_mm':w,'height_mm':h,'pages':1,'raster_images':len(pg.get_images()),'fonts':';'.join(sorted(set(x[2] for x in pg.get_fonts(full=True))))})
def box(ax,feature,ylabel,df=m,counts=True):
    vs=[];rng=np.random.default_rng(71)
    for i,g in enumerate(COL):
        v=df.loc[df.APC_class.eq(g),feature].dropna().to_numpy();vs.append(v)
        ax.scatter(i+rng.uniform(-.19,.19,len(v)),v,s=5.5,c=COL[g],alpha=.20 if len(v)>200 else .55,linewidths=0)
        ax.boxplot([v],positions=[i],widths=.49,showfliers=False,patch_artist=True,boxprops={'facecolor':'none','edgecolor':COL[g],'linewidth':.75},medianprops={'color':INK,'linewidth':.8},whiskerprops={'color':COL[g],'linewidth':.65},capprops={'color':COL[g],'linewidth':.65})
    ax.set_xlim(-.55,2.55);ax.set_xticks(range(3),[f'{g}\nn={len(v)}' if counts else g for g,v in zip(['WT','Missense','LOF'],vs)])
    ax.set_ylabel(ylabel);ax.yaxis.set_major_locator(MaxNLocator(4));return stats.kruskal(*vs).pvalue
def scatter(ax,feature,xlabel,scope='WT_LOF',df=None):
    if df is None:df=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])]
    d=df.dropna(subset=[feature,'TCF7L2_gene_effect'])
    for g in ['WT','Truncating/LOF']:
        v=d.loc[d.APC_class.eq(g)];ax.scatter(v[feature],v.TCF7L2_gene_effect,s=6,c=COL[g],alpha=.30 if len(d)>300 else .65,linewidths=0)
    ax.axhline(0,lw=.5,color=LIGHT,ls='--');ax.set_xlabel(xlabel);ax.set_ylabel('TCF7L2 gene effect')
    ax.xaxis.set_major_locator(MaxNLocator(4));ax.yaxis.set_major_locator(MaxNLocator(4))
    return ass.loc[ass.feature.eq(feature)&ass.scope.eq(scope)].iloc[0]
def scatterstats(fig,row,x,y,scope='WT_LOF'):
    pp='P_OLS' if scope=='CRC' else 'P_lineage'
    label(fig,x,y,f'rho = {row.rho:.3f}; {pp} = {p(row.P)}',6)
def group_panel(key,letter,title,feat,ylabel,w,h,df=m):
    f=new(letter,title,w,h);a=chart(f,w,h,top=11,bottom=9)
    pv=box(a,feat,ylabel,df);label(f,10,6,f'Kruskal-Wallis P = {p(pv)}',6)
    return save(f,key,w,h,f'{feat}; group comparison; n printed below classes')
def scatter_panel(key,letter,title,feat,xlabel,w,h,scope='WT_LOF',df=None):
    f=new(letter,title,w,h);a=chart(f,w,h,top=11,bottom=9)
    r=scatter(a,feat,xlabel,scope,df);scatterstats(f,r,10,6,scope)
    return save(f,key,w,h,f'n={int(r.n)}; Spearman P={r.Spearman_P:.6g}; beta/SD={r.beta:.6g}; HC1 P={r.P:.6g}; scope={scope}')
def schematic(fig,x,y,w,h,text,edge=GREY,fs=6.5):
    fw,fh=fig.get_size_inches()*25.4
    patch=FancyBboxPatch((x/fw,1-(y+h)/fh),w/fw,h/fh,boxstyle='round,pad=0.002',transform=fig.transFigure,facecolor='#F6F7F9',edgecolor=edge,lw=.7)
    fig.add_artist(patch)
    label(fig,x+w/2,y+h/2,text,fs,ha='center',multialignment='center').set_verticalalignment('center')
def link(fig,x1,y1,x2,y2,arrow=False,dash=False):
    fw,fh=fig.get_size_inches()*25.4
    if arrow:
        aa=axmm(fig,0,0,fw,fh);aa.axis('off');aa.annotate('',xy=(x2/fw,1-y2/fh),xytext=(x1/fw,1-y1/fh),arrowprops={'arrowstyle':'->','lw':.65,'color':GREY,'linestyle':'--' if dash else '-'})
    else:fig.add_artist(Line2D([x1/fw,x2/fw],[1-y1/fh,1-y2/fh],transform=fig.transFigure,lw=.65,color=GREY,ls='--' if dash else '-'))

AN=ROOT/'output/coFig1_analysis'
s=pd.read_csv(AN/'genomewide_LOF_vs_WT.csv');fam=pd.read_csv(AN/'family_LOF_vs_WT.csv');adj=pd.read_csv(AN/'lineage_adjusted_LOF_effect.csv').iloc[0]
summary=json.loads((AN/'summary.json').read_text())
f=new('A','Workflow',32,46)
for y,h,txt in [(6,5,'DepMap CRISPR'),(14,8,'APC WT / missense\n/ truncating-LOF'),(25,5,'Genome-wide screen'),(33,5,'TCF7L2 candidate'),(41,4,'RNA / protein / Wnt')]:
    schematic(f,1,y,30,h,txt,fs=6)
    if y!=41:link(f,16,y+h,16,y+h+3,True)
save(f,'Fig1A',32,46,'Integrated screen and molecular profiling; RNA/protein/Wnt are characterized in F-H.')

f=new('B','Genome-wide differential dependency',94,46);a=chart(f,94,46,top=10,bottom=9)
t=s.loc[s.tested].copy();t['logq']=-np.log10(t.FDR.clip(lower=1e-300));sig=t.FDR.lt(.05)
for mask,col,sz in [(~sig,LIGHT,3),(sig&t.mean_difference.lt(0),COL['Truncating/LOF'],4),(sig&t.mean_difference.ge(0),COL['WT'],4)]:a.scatter(t.loc[mask,'mean_difference'],t.loc[mask,'logq'],s=sz,c=col,alpha=.55,linewidths=0)
a.axhline(-np.log10(.05),lw=.6,ls='--',color=GREY);a.axvline(0,lw=.5,ls='--',color=LIGHT)
a.set_xlabel('Mean gene-effect difference (APC LOF - WT)');a.set_ylabel('-log10(FDR)');a.set_ylim(-.4,t.logq.max()*1.16)
for gene,dx,dy in [('TCF7L2',-.15,2.8),('CTNNB1',.09,1.0)]:
    r=t.loc[t.gene.eq(gene)].iloc[0];a.annotate(gene,(r.mean_difference,r.logq),xytext=(r.mean_difference+dx,r.logq+dy),fontsize=6.5,fontweight='bold',color=COL['Truncating/LOF'],arrowprops={'arrowstyle':'-','lw':.6,'color':GREY})
label(f,10,5.5,f'WT n=1005 / LOF n=76; {int(t.tested.sum()):,} genes',6)
save(f,'Fig1B',94,46,'26Q1 APC classes; two-sided MW and genome-wide BH; unadjusted screen; missense excluded')

f=new('C','APC functional class',48,46);a=chart(f,48,46,top=12,bottom=9)
pv=box(a,'TCF7L2_gene_effect','TCF7L2 gene effect');label(f,10,6,f'KW P = {p(pv)}',6)
save(f,'Fig1C',48,46,'All three APC classes; sample counts beneath class labels')

f=new('D','Lineage-adjusted\nAPC LOF effect',34,44);a=axmm(f,8,16,24,15)
a.errorbar(adj.beta,0,xerr=[[adj.beta-adj.CI_low],[adj.CI_high-adj.beta]],fmt='o',ms=3,lw=.8,capsize=2,color=COL['Truncating/LOF']);a.axvline(0,color=GREY,ls='--',lw=.6)
a.set_xlim(-.42,.06);a.set_ylim(-1,1);a.set_yticks([]);a.set_xticks([-.4,-.2,0]);a.set_xlabel('APC LOF beta (95% CI)',fontsize=6)
label(f,17,11,f'P = {p(adj.P)}',6,ha='center');label(f,17,38.5,f'n={int(adj.n)}; WT vs LOF',6,ha='center')
save(f,'Fig1D',34,44,'GE ~ APC_LOF + lineage; HC1 t inference; no missense models')

f=new('E','TCF/Wnt family specificity',54,44);a=axmm(f,13,8,28,25)
for i,r in fam.iterrows():
    col=COL['Truncating/LOF'] if r.gene=='TCF7L2' else GREY
    a.errorbar(r.mean_difference,i,xerr=[[r.mean_difference-r.CI_low],[r.CI_high-r.mean_difference]],fmt='o',ms=2.4,lw=.7,capsize=1.4,color=col)
    a.text(1.04,i,'<0.001' if r.FDR_family8<.001 else f'{r.FDR_family8:.3f}',transform=a.get_yaxis_transform(),fontsize=6,va='center')
a.set_yticks(range(8),fam.gene);a.invert_yaxis();a.axvline(0,lw=.6,ls='--',color=GREY);a.set_xlabel('Mean difference\n(LOF - WT)',fontsize=6.5);a.xaxis.set_major_locator(MaxNLocator(3));a.text(1.04,1.05,'FDR',transform=a.transAxes,fontsize=6)
save(f,'Fig1E',54,44,'Mean differences and 5000-bootstrap pointwise CIs; MW tests with BH across eight genes')

group_panel('Fig1F','F','Wnt activity by APC class','Wnt_Hallmark','Hallmark Wnt ssGSEA score (z)',86,44)
f=new('G','TCF7L2 abundance across APC states',88,50)
for j,(feat,yl) in enumerate([('RNA26_TCF7L2','TCF7L2 RNA (log-TPM)'),('MS_TCF7L2','TCF7L2 protein')]):
    x=j*46;a=axmm(f,x+10,12,30,28);pv=box(a,feat,yl);label(f,x+10,7,f'KW P = {p(pv)}',6)
save(f,'Fig1G',88,50,'Two 42-mm slots and 4-mm gap; unadjusted RNA/protein comparisons; lineage-adjusted differences not significant')
f=new('H','TCF7L2 abundance vs dependency',88,50)
for j,(feat,xl) in enumerate([('RNA26_TCF7L2','TCF7L2 RNA (log-TPM)'),('MS_TCF7L2','TCF7L2 protein')]):
    x=j*46;a=axmm(f,x+10,14,30,27);r=scatter(a,feat,xl);label(f,x+10,6,f'rho = {r.rho:.3f}',6);label(f,x+10,9,f'P_lineage = {p(r.P)}',6)
save(f,'Fig1H',88,50,'RNA n1016; protein n270; Spearman rho and lineage-adjusted HC1 P; two 42-mm slots, 4-mm gap')

f=new('A','Top APC LOF-associated dependencies',88,50);a=axmm(f,15,10,70,32)
top=pd.read_csv(AN/'top15_LOF_dependencies.csv').iloc[::-1]
a.hlines(np.arange(len(top)),top.mean_difference,0,color=LIGHT,lw=.7);a.scatter(top.mean_difference,np.arange(len(top)),s=9,c=[COL['Truncating/LOF'] if g in ['TCF7L2','CTNNB1'] else GREY for g in top.gene])
a.set_yticks(range(len(top)),top.gene);a.tick_params(axis='y',length=0);a.set_xlabel('Mean gene-effect difference (LOF - WT)');a.axvline(0,lw=.6,color=GREY);a.xaxis.set_major_locator(MaxNLocator(4))
label(f,15,5.8,'Top 15 by mean difference among FDR <0.05',6)
save(f,'FigS1A',88,50,'Top 15 negative hits sorted by effect among genome-wide FDR<.05; candidates not forced into ranking')

f=new('B','CRC dependency landscape',88,50);a=axmm(f,10,10,75,31)
c=pd.read_csv(AN/'CRC_dependency_landscape.csv');names=['SW620','HT29','LOVO','RKO']
for group,col in [('WT',COL['WT']),('Truncating/LOF',COL['Truncating/LOF'])]:
    v=c.loc[c.APC_class.eq(group)];a.vlines(v['rank'],0,v.TCF7L2_gene_effect,lw=.45,color=col,alpha=.3);a.scatter(v['rank'],v.TCF7L2_gene_effect,s=7,c=col,linewidths=0)
a.axhline(0,lw=.5,ls='--',color=LIGHT);a.set_ylim(-2,.3);a.set_xlim(0,len(c)+1);a.set_xlabel('CRC rank (stronger to weaker dependency)');a.set_ylabel('TCF7L2 gene effect');a.xaxis.set_major_locator(MaxNLocator(5))
slots={'SW620':(8,-1.85),'HT29':(23,-1.6),'LOVO':(40,-1.1),'RKO':(51,-.75)}
for n,(x,y) in slots.items():
    r=c.loc[c.cell_line_display_name.eq(n)].iloc[0];a.annotate(n,(r['rank'],r.TCF7L2_gene_effect),xytext=(x,y),fontsize=6.5,fontweight='bold',color=COL[r.APC_class],arrowprops={'arrowstyle':'-','lw':.6,'color':COL[r.APC_class]})
label(f,10,5.6,f'WT n={c.APC_class.eq("WT").sum()} / LOF n={c.APC_class.eq("Truncating/LOF").sum()}',6)
save(f,'FigS1B',88,50,'CRC WT/LOF only; SW620, HT29, LOVO and RKO labeled with actual values')

f=new('C','Top 15 TCF7L2 co-dependencies',180,58);a=axmm(f,17,8,159,42)
top=co.sort_values(['rho','gene'],ascending=[False,True]).head(15).iloc[::-1];a.hlines(np.arange(15),0,top.rho,color=LIGHT,lw=.7);a.scatter(top.rho,np.arange(15),s=10,c=[COL['Truncating/LOF'] if g in ['CTNNB1','FZD5','BCL9L'] else GREY for g in top.gene])
a.set_yticks(range(15),top.gene);a.tick_params(axis='y',length=0);a.set_xlim(0,top.rho.max()*1.05);a.set_xlabel('Lineage-residualized Spearman rho');a.xaxis.set_major_locator(MaxNLocator(6))
save(f,'FigS1C',180,58,'Co-dependency residualized within lineage; TCF7L2 excluded; full ranking retained')

group_panel('FigS2A','A','PROGENy WNT by APC class','Wnt_PROGENy','PROGENy WNT weighted score (z)',88,52)
scatter_panel('FigS2B','B','PROGENy WNT vs dependency','Wnt_PROGENy','PROGENy WNT weighted score (z)',88,52)
group_panel('FigS2C','C','Seven-gene score by APC class','Wnt_custom','Seven-gene mean-z score',88,52)
scatter_panel('FigS2D','D','Seven-gene score vs dependency','Wnt_custom','Seven-gene mean-z score',88,52)
crc=m.loc[m.lineage_2.eq('Colorectal Adenocarcinoma')];crcwl=crc.loc[crc.APC_class.isin(['WT','Truncating/LOF'])]
group_panel('FigS3A','A','CRC: Hallmark Wnt by APC class','Wnt_Hallmark','Hallmark Wnt score (z)',88,52,crc)
for letter,feat,xl in [('B','Wnt_Hallmark','Hallmark Wnt score (z)'),('C','RNA26_TCF7L2','TCF7L2 RNA (log-TPM)'),('D','MS_TCF7L2','TCF7L2 protein')]:scatter_panel('FigS3'+letter,letter,'CRC: '+('Wnt' if letter=='B' else 'RNA' if letter=='C' else 'protein')+' vs dependency',feat,xl,88,52,'CRC',crcwl)

def composite(name,w,h,placements):
    d=fitz.open();pg=d.new_page(width=w*MM,height=h*MM)
    for key,x,y in placements:
        pw,ph=meta[key]['size_mm'];src=fitz.open(PAN/(key+'.pdf'));pg.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+pw)*MM,(y+ph)*MM),src,0,keep_proportion=False)
    dest=OUT/(name+'.pdf');d.save(dest,garbage=4,deflate=True);d.close();verify(dest,w,h)
composite('Fig1_composite',180,146,[('Fig1A',0,0),('Fig1B',35,0),('Fig1C',132,0),('Fig1D',0,49),('Fig1E',37,49),('Fig1F',94,49),('Fig1G',.5,96),('Fig1H',91.5,96)])
composite('FigS1_composite',180,112,[('FigS1A',0,0),('FigS1B',92,0),('FigS1C',0,54)])
for nm in ['FigS2','FigS3']:composite(nm+'_composite',180,108,[(nm+chr(65+i),(i%2)*92,(i//2)*56) for i in range(4)])
pd.DataFrame(manifest).to_csv(OUT/'PDF_manifest.csv',index=False);(OUT/'panel_metadata.json').write_text(json.dumps(meta,indent=2))
(TMP/'text_bounds.json').write_text(json.dumps(bounds,indent=2))
print('PDFs:',len(manifest),'Bounds:',len(bounds),flush=True)
for b in bounds:print(b['file'],ascii(b['text']))
assert len(manifest)==23 and not bounds

