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
OUT=ROOT/'output/pdf/reFig2_final'; OUT.mkdir(parents=True,exist_ok=True)
PAN=OUT/'individual_panels'; PAN.mkdir(exist_ok=True)
TMP=ROOT/'tmp/pdfs/reFig2';TMP.mkdir(parents=True,exist_ok=True)
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

# Main figure: exact prescribed 42/58/74, 78/99, 103/74 mm geometry.
f=new('A','Study design',42,44)
schematic(f,1.3,7,39.4,9,'DepMap cohort\nAPC WT / Missense / LOF',fs=6.5)
for x,w,s in [(1.3,10,'RNA'),(13.1,11,'Protein'),(25.9,14.8,'DNA\nmethylation')]:
    schematic(f,x,21,w,8,s,fs=6);link(f,21,16,x+w/2,21,True);link(f,x+w/2,29,21,34,True)
schematic(f,2.3,34,37.4,8,'TCF7L2 molecular state\n& dependency',edge=COL['WT'],fs=6.5)
save(f,'Fig2A',42,44,'Navigation schematic; 26Q1 RNA/mutation integrated with existing protein/methylation/CRISPR')
f=new('B','Wnt activity by APC class',58,44);a=chart(f,58,44,top=13,bottom=9)
pv=box(a,'Wnt_Hallmark','Hallmark Wnt score (z)');r=pw.loc[pw.feature.eq('Wnt_Hallmark')&pw.comparison.eq('Truncating/LOF - WT')].iloc[0]
label(f,10,5.6,f'Kruskal-Wallis P = {p(pv)}',6);label(f,10,8.6,f'LOF vs WT: Holm P = {p(r.P_Holm)}',6)
save(f,'Fig2B',58,44,'Hallmark ssGSEA; WT944/missense65/LOF72')
scatter_panel('Fig2C','C','Wnt activity vs dependency','Wnt_PROGENy','PROGENy WNT weighted score (z)',74,44)
f=new('D','TCF7L2 abundance across APC states',78,50)
for j,(feat,yl) in enumerate([('RNA26_TCF7L2','TCF7L2 RNA (log-TPM)'),('MS_TCF7L2','TCF7L2 protein')]):
    x=j*41;a=axmm(f,x+9,12,27,28);pv=box(a,feat,yl);label(f,x+9,7,f'KW P = {p(pv)}',6)
save(f,'Fig2D',78,50,'Two 37-mm slots separated by 4 mm; normalized MS protein; no methylation')
f=new('E','TCF7L2 abundance vs dependency',99,50)
for j,(feat,xl) in enumerate([('RNA26_TCF7L2','TCF7L2 RNA (log-TPM)'),('MS_TCF7L2','TCF7L2 protein')]):
    x=j*52;a=axmm(f,x+10,14,35,27);r=scatter(a,feat,xl)
    label(f,x+10,6,f'rho = {r.rho:.3f}',6);label(f,x+10,9,f'P_lineage = {p(r.P)}',6)
save(f,'Fig2E',99,50,'Two 47-mm slots separated by 5 mm; RNA n1016/protein n270; OLS HC1 with lineage')
f=new('F','Top 15 lineage-adjusted co-dependencies',103,43);a=axmm(f,16,6,84,30)
top=co.sort_values(['rho','gene'],ascending=[False,True]).head(15).iloc[::-1]
a.hlines(np.arange(15),0,top.rho,color=LIGHT,lw=.7)
a.scatter(top.rho,np.arange(15),s=8,c=[COL['Truncating/LOF'] if x in ['CTNNB1','FZD5','BCL9L'] else GREY for x in top.gene])
a.set_yticks(np.arange(15),top.gene);a.tick_params(axis='y',length=0);a.set_ylim(-.7,14.7);a.set_xlim(0,top.rho.max()*1.07);a.set_xlabel('Lineage-residualized Spearman rho');a.xaxis.set_major_locator(MaxNLocator(5))
save(f,'Fig2F',103,43,'Top 15 by positive residual Spearman; TCF7L2 excluded; red highlights CTNNB1/FZD5/BCL9L, not significance')
f=new('G','Evidence-based working model',74,43)
label(f,37,5.2,'Observed associations',6,ha='center')
schematic(f,25,9,24,5,'APC truncating/LOF',fs=6)
schematic(f,1.5,17,29,6,'Wnt-active state',fs=6.5);schematic(f,37,17,35,6,'TCF7L2 RNA/protein',fs=6.5)
link(f,30,14,16,17);link(f,44,14,54,17)
schematic(f,21,26,34,5,'TCF7L2 dependency',edge=COL['Truncating/LOF'],fs=6.5)
link(f,16,23,28,26);link(f,54,23,48,26)
label(f,37,33,'Causal hypothesis: APC loss -?-> dependency',6,ha='center')
label(f,37,37.4,'To test: isogenic APC-WT / APC-KO',6,ha='center')
save(f,'Fig2G',74,43,'Lines are observational associations; abundance group differences are unadjusted. No claim that isogenic validation is completed.')

# Supplementary S2: 86 x 50 panels, 4 mm internal gaps, 2 mm external margins.
group_panel('FigS2A','A','PROGENy WNT by APC class','Wnt_PROGENy','PROGENy WNT weighted score (z)',86,50)
scatter_panel('FigS2B','B','PROGENy WNT vs dependency','Wnt_PROGENy','PROGENy WNT weighted score (z)',86,50)
group_panel('FigS2C','C','Seven-gene score by APC class','Wnt_custom','Seven-gene mean-z score',86,50)
scatter_panel('FigS2D','D','Seven-gene score vs dependency','Wnt_custom','Seven-gene mean-z score',86,50)

# S3 has shared row/column headings in the composite; the independent PDFs identify their variables by axes.
for i,reg in enumerate(meth):
    for j,(tar,yl) in enumerate([('RNA26_TCF7L2','RNA (log-TPM)'),('MS_TCF7L2','Protein abundance'),('TCF7L2_gene_effect','TCF7L2 gene effect')]):
        letter=chr(65+i*3+j);f=new(letter,'',55,46);a=chart(f,55,46,top=10,bottom=9)
        dd=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=[reg,tar])
        for gg in ['WT','Truncating/LOF']:
            v=dd.loc[dd.APC_class.eq(gg)];a.scatter(v[reg],v[tar],s=5.5,c=COL[gg],alpha=.4,linewidths=0)
        r=mr.loc[mr.region.eq(reg)&mr.target.eq(tar)&mr.group.eq('WT_LOF')].iloc[0]
        label(f,10,2,f'rho = {r.rho:.2f}',6);label(f,10,5.3,f'FDR = {r.FDR:.3f}',6)
        a.set_xlabel(f'Region {i+1} methylation');a.set_ylabel(yl);a.xaxis.set_major_locator(MaxNLocator(3));a.yaxis.set_major_locator(MaxNLocator(4))
        save(f,f'FigS3{letter}',55,46,f'{reg} versus {tar}; n={int(r.n)}; FDR of lineage OLS across 27 models')

def forest_panel(key,letter,title,mode):
    f=new(letter,title,86,52);a=axmm(f,23,11,50,32);dd=fr.loc[fr['mode'].eq(mode)].reset_index(drop=True)
    for i,r in dd.iterrows():a.errorbar(r.beta,i,xerr=[[r.beta-r.CI_low],[r.CI_high-r.beta]],fmt='o',ms=2.5,lw=.7,capsize=1.5,color=COL['Truncating/LOF'] if i==5 else COL['WT'])
    a.axvline(0,color=GREY,lw=.6,ls='--');a.set_yticks(range(6),dd.model);a.invert_yaxis();a.set_xlabel('APC LOF coefficient (95% CI)');a.xaxis.set_major_locator(MaxNLocator(4))
    for i,r in dd.iterrows():a.text(1.04,i,p(r.P),transform=a.get_yaxis_transform(),fontsize=6,va='center')
    a.text(1.04,1.07,'P',transform=a.transAxes,fontsize=6)
    label(f,23,6,f'Full model: n={int(dd.iloc[-1].n)}, LOF={int(dd.iloc[-1].LOF)}',6)
    save(f,key,86,52,f'{mode}; all six coefficients/CI/sample sizes in source table')
for letter,mode,title in zip('ABCD',['All_available']+meth,['All available samples','Region 1 complete cases','Region 2 complete cases','Region 3 complete cases']):forest_panel('FigS4'+letter,letter,title,mode)

# S5: 86 x 48 panels in upper 2 x 2; 86 x 40 enrichment panels below.
for group_letter,scat_letter,feat,title in [('A','B','MS_CTNNB1','CTNNB1 mass spectrometry'),('C','D','RPPA_beta_catenin','Total beta-catenin RPPA')]:
    group_panel('FigS5'+group_letter,group_letter,title+' by APC',feat,'Normalized protein abundance',86,48)
    scatter_panel('FigS5'+scat_letter,scat_letter,title+' vs dependency',feat,title,86,48)
for letter,omics in [('E','RNA'),('F','Protein')]:
    f=new(letter,omics+' Hallmark pathways',86,40);a=axmm(f,43,8,32,25)
    gs=pd.read_csv(T/f'{omics}_GSEA_beta.csv');v=gs.loc[gs.Term.str.startswith('HALLMARK_')].sort_values(['FDR q-val','Term']).head(5).sort_values('NES')
    a.barh(np.arange(len(v)),v.NES,height=.65,color=[COL['Truncating/LOF'] if n>0 else COL['WT'] for n in v.NES],alpha=.8)
    a.set_yticks(np.arange(len(v)),[textwrap.fill(x.replace('HALLMARK_','').replace('_',' '),26) for x in v.Term]);a.tick_params(axis='y',length=0,labelsize=6)
    a.axvline(0,color=GREY,lw=.6);a.set_xlabel('NES');a.xaxis.set_major_locator(MaxNLocator(4))
    for i,q in enumerate(v['FDR q-val']):a.text(1.03,i,'<0.001' if q<.001 else f'{q:.3f}',transform=a.get_yaxis_transform(),fontsize=6,va='center')
    a.text(1.03,1.05,'FDR',transform=a.transAxes,fontsize=6)
    save(f,'FigS5'+letter,86,40,'Five lowest-FDR Hallmark sets to preserve >=6pt text in 40mm height; -standardized beta ranking; full tables retained')

f=new('A','APC LOF and minor-copy-number status',88,70);a=chart(f,88,70,top=16,bottom=13)
egs=['WT with CN','LOF; minor CN > 0','LOF; minor CN = 0'];rng=np.random.default_rng(53)
for i,g in enumerate(egs):
    v=m.loc[m.LOH_evidence.eq(g),'TCF7L2_gene_effect'].dropna();c=[COL['WT'],COL['Missense only'],COL['Truncating/LOF']][i]
    a.scatter(i+rng.uniform(-.18,.18,len(v)),v,s=6,c=c,alpha=.3,linewidths=0);a.boxplot([v],positions=[i],widths=.5,showfliers=False,medianprops={'color':INK,'linewidth':.8},boxprops={'linewidth':.7},whiskerprops={'linewidth':.7},capprops={'linewidth':.7})
a.set_xticks([0,1,2],['WT\nn=803','LOF; minor CN >0\nn=27','LOF; minor CN =0\nn=26']);a.tick_params(axis='x',labelsize=6);a.set_ylabel('TCF7L2 gene effect')
label(f,10,6,'Minor CN=0 vs WT: Holm P=0.0004',6);label(f,10,9.2,'Within LOF, minor CN=0 vs >0: P=0.331',6)
label(f,44,66.4,'Allelic-loss evidence; not confirmed biallelic LOF',6,ha='center')
save(f,'FigS6A',88,70,'LOH evidence strata; OLS HC1 lineage adjusted, Holm three comparisons')
med=pd.read_csv(T/'exploratory_mediation.csv').iloc[0]
f=new('B','Exploratory statistical mediation',88,70)
for x,s in [(1.5,'APC LOF'),(31.5,'Hallmark\nWnt score'),(61.5,'TCF7L2\ngene effect')]:schematic(f,x,19,25,12,s,fs=6.5)
link(f,26.5,25,31.5,25,True);link(f,56.5,25,61.5,25,True)
label(f,29,34,f'a={med.a:.3f}',6,ha='center');label(f,59,34,f'b={med.b:.3f}',6,ha='center')
label(f,44,43,f'Indirect association = {med.indirect:.3f}\n95% CI [{med.CI_low:.3f}, {med.CI_high:.3f}]\nn=1016; 2,000 lineage-stratified resamples',6.5,ha='center',linespacing=1.35)
label(f,44,61.5,'Cross-sectional observational data\ndo not establish causal mediation.',6,ha='center')
save(f,'FigS6B',88,70,'Exploratory non-causal mediation; same lineage residualization as source analysis')
crc=m.loc[m.lineage_2.eq('Colorectal Adenocarcinoma')];crcwl=crc.loc[crc.APC_class.isin(['WT','Truncating/LOF'])]
group_panel('FigS7A','A','CRC: Hallmark Wnt by APC class','Wnt_Hallmark','Hallmark Wnt score (z)',86,50,crc)
for letter,feat,xl in [('B','Wnt_Hallmark','Hallmark Wnt score (z)'),('C','RNA26_TCF7L2','TCF7L2 RNA (log-TPM)'),('D','MS_TCF7L2','TCF7L2 protein')]:scatter_panel('FigS7'+letter,letter,'CRC: '+('Wnt' if letter=='B' else 'RNA' if letter=='C' else 'protein')+' vs dependency',feat,xl,86,50,'CRC',crcwl)

def composite(name,w,h,placements,headings=None):
    d=fitz.open();pg=d.new_page(width=w*MM,height=h*MM)
    for key,x,y in placements:
        pw,ph=meta[key]['size_mm'];src=fitz.open(PAN/(key+'.pdf'))
        pg.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+pw)*MM,(y+ph)*MM),src,0,keep_proportion=False)
    if headings:
        pg.insert_font(fontname='Arial',fontfile=matplotlib.font_manager.findfont('Arial'))
        for x,y,txt,rot in headings:pg.insert_text((x*MM,y*MM),txt,fontname='Arial',fontsize=6.5,rotate=rot,color=(.14,.16,.18))
    dest=OUT/(name+'.pdf');d.save(dest,garbage=4,deflate=True);d.close();verify(dest,w,h)
composite('Fig2_composite',180,143,[('Fig2A',0,0),('Fig2B',45,0),('Fig2C',106,0),('Fig2D',0,47),('Fig2E',81,47),('Fig2F',0,100),('Fig2G',106,100)])
grid=lambda pre,h:[(pre+chr(65+i),2+(i%2)*90,2+(i//2)*(h+4)) for i in range(4)]
composite('FigS2_composite',180,108,grid('FigS2',50))
heads=[(20,4.5,'RNA expression',0),(77,4.5,'Protein abundance',0),(132,4.5,'TCF7L2 dependency',0)]
heads +=[(3,41+i*50,'TSS region '+str(i+1),90) for i in range(3)]
composite('FigS3_composite',180,158,[(f'FigS3{chr(65+i*3+j)}',7+j*59,8+i*50) for i in range(3) for j in range(3)],heads)
composite('FigS4_composite',180,112,grid('FigS4',52))
composite('FigS5_composite',180,150,[(f'FigS5{chr(65+i)}',2+(i%2)*90,3+(i//2)*52) for i in range(4)]+[('FigS5E',2,107),('FigS5F',92,107)])
composite('FigS6_composite',180,75,[('FigS6A',0,2.5),('FigS6B',92,2.5)])
composite('FigS7_composite',180,108,grid('FigS7',50))
pd.DataFrame(manifest).to_csv(OUT/'PDF_manifest.csv',index=False)
(OUT/'panel_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
(TMP/'text_bounds.json').write_text(json.dumps(bounds,indent=2),encoding='utf-8')
print('PDF count:',len(manifest),'out-of-page text artists:',len(bounds),flush=True)
for b in bounds:print('BOUNDS',b['file'],ascii(b['text']),flush=True)
assert len(manifest)==43
assert not bounds,'Out-of-page text requires correction'
