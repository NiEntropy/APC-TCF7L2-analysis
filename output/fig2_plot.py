from pathlib import Path
import json,textwrap,zipfile
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import pymupdf

ROOT=Path(__file__).resolve().parents[1];O=ROOT/'output/Fig2_updated';T=O/'tables';P=ROOT/'output/pdf/Fig2_updated';P.mkdir(exist_ok=True,parents=True);V=O/'previews';V.mkdir(exist_ok=True)
m=pd.read_csv(T/'master_analysis.csv',index_col=0)
pw=pd.read_csv(T/'group_comparisons.csv');assoc=pd.read_csv(T/'dependency_associations.csv');fr=pd.read_csv(T/'APC_multivariable_models.csv');gm=pd.read_csv(T/'lineage_adjusted_group_effects.csv')
meth=[c for c in m if c.startswith('METH_')];colors={'WT':'#377DA8','Missense only':'#BE9348','Truncating/LOF':'#C55358'};INK='#24282E';GREY='#66717D'
groups=list(colors)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'savefig.facecolor':'white'})
manifest=[]
def pstr(v):return f'{v:.3f}' if v>=.001 else f'{v:.1e}'
def base(letter,title,w=100,h=85,nrows=1,ncols=1):
    fig,axes=plt.subplots(nrows,ncols,figsize=(w/25.4,h/25.4),squeeze=False)
    fig.subplots_adjust(left=.16 if ncols==1 else .08,right=.97,top=.77,bottom=.20,wspace=.42,hspace=.73)
    fig.text(.025,.955,letter,fontsize=13,fontweight='bold',va='top',color=INK)
    fig.text(.09 if letter[0]!='S' else .11,.949,title,fontsize=10,fontweight='bold',va='top',color=INK)
    return fig,axes
def save(fig,name,note=''):
    if note:fig.text(.04,.025,note,fontsize=6.6,color=GREY,va='bottom')
    path=P/(name+'.pdf');fig.savefig(path);plt.close(fig)
    d=pymupdf.open(path);page=d[0];assert len(d)==1 and not page.get_images()
    spans=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
    over=[s['text'] for s in spans if s['bbox'][0]<-.5 or s['bbox'][1]<-.5 or s['bbox'][2]>page.rect.width+.5 or s['bbox'][3]>page.rect.height+.5]
    if over:print('OVERFLOW',name,over,flush=True)
    page.get_pixmap(matrix=pymupdf.Matrix(3,3)).save(V/(name+'.png'))
    manifest.append(dict(file=path.name,width_mm=page.rect.width*25.4/72,height_mm=page.rect.height*25.4/72,images=len(page.get_images()),overflow=';'.join(over)))
def groupplot(ax,feature,label,df=m,showtests=True):
    vals=[df.loc[df.APC_class.eq(g),feature].dropna().to_numpy() for g in groups]
    rng=np.random.default_rng(101)
    for i,(g,v) in enumerate(zip(groups,vals)):
        ax.scatter(i+rng.uniform(-.20,.20,len(v)),v,s=6 if len(v)>200 else 10,c=colors[g],alpha=.18 if len(v)>200 else .55,linewidths=0,rasterized=False,zorder=1)
        ax.boxplot([v],positions=[i],widths=.50,showfliers=False,patch_artist=True,boxprops={'facecolor':'none','edgecolor':colors[g],'linewidth':1},medianprops={'color':INK,'linewidth':1.3},whiskerprops={'color':colors[g]},capprops={'color':colors[g]})
    ax.set_xticks(range(3),[f'{g}\nn = {len(v)}' for g,v in zip(['WT','Missense','LOF'],vals)])
    ax.set_ylabel(label);ax.set_xlim(-.55,2.55)
    if showtests:
        q=stats.kruskal(*vals).pvalue;ax.set_title(f'Kruskal-Wallis P = {pstr(q)}',fontsize=8)
def scatter(ax,feature,label,scope='WT_LOF',df=None,report=True):
    if df is None:df=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])]
    a=df.dropna(subset=[feature,'TCF7L2_gene_effect'])
    for g in ['WT','Truncating/LOF']:
        v=a.loc[a.APC_class.eq(g)];ax.scatter(v[feature],v.TCF7L2_gene_effect,s=9,c=colors[g],alpha=.38 if len(a)>300 else .65,linewidths=0,rasterized=False)
    ax.set_xlabel(label);ax.set_ylabel('TCF7L2 gene effect')
    ax.axhline(0,color='#CCD3D9',lw=.6,ls='--')
    if report:
        r=assoc.loc[assoc.feature.eq(feature)&assoc.scope.eq(scope)].iloc[0]
        ax.set_title(f'n = {int(r.n)}; Spearman rho = {r.rho:.3f}\n'
                     f'P = {pstr(r.Spearman_P)}; '+('OLS' if scope=='CRC' else 'adjusted')+f' beta = {r.beta:.3f}\n'
                     +('OLS' if scope=='CRC' else 'adjusted')+f' P = {pstr(r.P)}',fontsize=7.5,loc='left')
def legend(fig,y=.845):
    fig.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor=colors[g],markeredgecolor='none',markersize=4,label='APC '+('LOF' if g!='WT' else 'WT')) for g in ['WT','Truncating/LOF']],loc='upper center',bbox_to_anchor=(.53,y),ncol=2,frameon=False,fontsize=7)
def forest(ax,dd,title):
    dd=dd.reset_index(drop=True);y=np.arange(len(dd))
    for i,r in dd.iterrows():ax.errorbar(r.beta,i,xerr=[[r.beta-r.CI_low],[r.CI_high-r.beta]],fmt='o',color='#C55358' if i==len(dd)-1 else '#377DA8',markersize=4,capsize=2,lw=1)
    ax.axvline(0,color=GREY,lw=.7,ls='--');ax.set_yticks(y,dd['model']);ax.invert_yaxis();ax.set_xlabel('APC LOF coefficient (95% CI)');ax.set_title(title,fontsize=8)
    for i,r in dd.iterrows():ax.text(1.03,i,f'{pstr(r.P)}',transform=ax.get_yaxis_transform(),fontsize=7,va='center')
    ax.text(1.03,1.02,'P',transform=ax.transAxes,fontsize=7)
def schematic_box(ax,xy,w,h,lines,col='#377DA8',size=9):
    x,y=xy;ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012',edgecolor=col,facecolor='#F6F7F9',linewidth=1,clip_on=False))
    ax.text(x+w/2,y+h/2,lines,ha='center',va='center',fontsize=size,color=INK)
def arrow(ax,x1,y1,x2,y2):ax.annotate('',(x2,y2),(x1,y1),arrowprops={'arrowstyle':'->','color':GREY,'lw':1})

fig,aa=base('A','Multi-omics study design',110,85);ax=aa[0,0];ax.axis('off');ax.set_position([.04,.09,.92,.78]);ax.set_xlim(0,1);ax.set_ylim(0,1)
schematic_box(ax,(.23,.79),.54,.16,'Model-matched DepMap cohort\nAPC WT / missense / LOF',size=8)
for x,w,label,col in [(.01,.29,'RNA-seq\nn = 1,082','#377DA8'),(.355,.29,'TCF7L2 protein\nn = 294','#BE9348'),(.70,.29,'Methylation\n3 TSS regions','#C55358')]:
    schematic_box(ax,(x,.42),w,.20,label,col,8);arrow(ax,.5,.78,x+w/2,.64);arrow(ax,x+w/2,.40,.5,.25)
schematic_box(ax,(.12,.03),.76,.21,'TCF7L2 molecular associations\nFull WT/LOF complete cases: 190 (16 LOF)',size=8)
save(fig,'Fig2A_Study_design','26Q1 RNA/mutations + existing model-level CRISPR, proteomics and methylation.')

fig,aa=base('B','Wnt activity by APC class',100,90);groupplot(aa[0,0],'Wnt_Hallmark','Hallmark Wnt ssGSEA score (z)')
r=pw.loc[pw.feature.eq('Wnt_Hallmark')&pw.comparison.eq('Truncating/LOF - WT')].iloc[0]
fig.text(.16,.86,f'LOF vs WT: Holm P = {pstr(r.P_Holm)}',fontsize=8)
save(fig,'Fig2B_Wnt_by_APC','MSigDB Hallmark WNT/BETA-CATENIN; score is standardized across RNA samples.')

fig,aa=base('C','Wnt activity vs dependency',110,95);fig.subplots_adjust(top=.69,bottom=.20);scatter(aa[0,0],'Wnt_Hallmark','Hallmark Wnt ssGSEA score (z)');legend(fig,.885)
save(fig,'Fig2C_Wnt_dependency','Adjusted beta: gene-effect change per 1 SD score, controlling for lineage (HC1).')

fig,aa=base('D','TCF7L2 abundance across APC states',180,90,ncols=3);fig.subplots_adjust(top=.74,bottom=.24,wspace=.47)
groupplot(aa[0,0],'RNA26_TCF7L2','RNA (provided log-TPM)');groupplot(aa[0,1],'MS_TCF7L2','Protein (normalized abundance)')
ax=aa[0,2];rng=np.random.default_rng(43)
for j,reg in enumerate(meth):
    for k,g in enumerate(groups):
        v=m.loc[m.APC_class.eq(g),reg].dropna();pos=j*4+k
        ax.boxplot([v],positions=[pos],widths=.65,showfliers=False,patch_artist=True,boxprops={'facecolor':colors[g]+'22','edgecolor':colors[g]},medianprops={'color':INK},whiskerprops={'color':colors[g]},capprops={'color':colors[g]})
ax.set_xticks([1,5,9],['Region 1','Region 2','Region 3']);ax.set_ylabel('TSS-region methylation');ax.set_title('Three TSS regions\nNo outcome-based selection',fontsize=7.5)
fig.legend(handles=[Line2D([],[],color=colors[g],lw=2,label=g) for g in groups],loc='lower right',bbox_to_anchor=(.98,.08),ncol=3,fontsize=6.5,frameon=False)
save(fig,'Fig2D_TCF7L2_abundance','Methylation: existing 1-kb upstream-TSS regions; region coordinates and all tests in tables.')

fig,aa=base('E','TCF7L2 abundance vs dependency',180,95,ncols=2);fig.subplots_adjust(top=.69,bottom=.22,wspace=.35)
scatter(aa[0,0],'RNA26_TCF7L2','TCF7L2 RNA (provided log-TPM)');scatter(aa[0,1],'MS_TCF7L2','TCF7L2 protein (normalized)');legend(fig,.89)
save(fig,'Fig2E_Abundance_dependency','Spearman: unadjusted; adjusted beta/P: lineage-adjusted OLS with HC1 robust inference.')

fig,aa=base('F','Multi-omic adjusted APC effect',125,95);fig.subplots_adjust(left=.29,right=.84,top=.77,bottom=.19)
dd=fr.loc[fr['mode'].eq('Complete_case')];forest(aa[0,0],dd,'Same 190 cell lines in every model (174 WT / 16 LOF)')
save(fig,'Fig2F_Complete_case_APC_effect','+RNA / +protein / +methylation each add to lineage; full model contains all three regions.')

co=pd.read_csv(T/'codependency_ranking.csv');gs=pd.read_csv(T/'codependency_GSEA.csv')
def gsplot(ax,df,n=8,hallmark=True):
    v=df.loc[df.Term.str.startswith('HALLMARK_')] if hallmark else df
    v=v.assign(q=pd.to_numeric(v['FDR q-val']),nes=pd.to_numeric(v.NES)).sort_values(['q','Term']).head(n).sort_values('nes')
    labels=[textwrap.fill(x.replace('HALLMARK_','').replace('_',' '),29 if hallmark else 36) for x in v.Term]
    ax.barh(np.arange(len(v)),v.nes,color=['#C55358' if x>0 else '#377DA8' for x in v.nes],alpha=.8)
    ax.set_yticks(np.arange(len(v)),labels);ax.tick_params(axis='y',labelsize=7);ax.axvline(0,color=GREY,lw=.6);ax.set_xlabel('Normalized enrichment score')
    for i,r in enumerate(v.itertuples()):ax.text(1.03,i,('<0.001' if r.q<.001 else f'{r.q:.3f}'),transform=ax.get_yaxis_transform(),fontsize=6.5,va='center')
    ax.text(1.03,.995,'FDR',transform=ax.transAxes,fontsize=6.5)
fig,aa=base('G','TCF7L2 co-dependency and pathways',180,105,ncols=2);fig.subplots_adjust(left=.12,right=.90,top=.77,bottom=.15,wspace=1.18)
top=co.sort_values(['rho','gene'],ascending=[False,True]).head(10).iloc[::-1]
ax=aa[0,0];ax.barh(top.gene,top.rho,color=['#C55358' if g=='CTNNB1' else '#377DA8' for g in top.gene]);ax.set_xlabel('Residual Spearman rho');ax.set_title('Top 10 co-dependent genes',fontsize=8)
gsplot(aa[0,1],gs,n=7);aa[0,1].set_title('Hallmark pathways (top FDR)',fontsize=8)
save(fig,'Fig2G_Codependency_GSEA','Lineage-residualized gene effects; TCF7L2 excluded from ranking. Full pathway results supplied.')

fig,aa=base('H','Evidence-based working model',150,80);ax=aa[0,0];ax.axis('off');ax.set_position([.03,.15,.94,.70]);ax.set_xlim(0,1);ax.set_ylim(0,1)
schematic_box(ax,(.01,.53),.25,.29,'APC truncating/LOF',size=9)
schematic_box(ax,(.365,.53),.27,.29,'Wnt transcriptional\nstate',size=9)
schematic_box(ax,(.745,.53),.245,.29,'TCF7L2\ndependency',size=9)
ax.plot([.27,.355],[.675,.675],color=GREY,lw=1.3);ax.plot([.645,.735],[.675,.675],color=GREY,lw=1.3)
r1=gm.loc[gm.feature.eq('Wnt_Hallmark')].iloc[0];r2=assoc.loc[assoc.feature.eq('Wnt_Hallmark')&assoc.scope.eq('WT_LOF')].iloc[0]
ax.text(.31,.38,f'Lineage-adjusted\nP = {pstr(r1.P)}',ha='center',fontsize=8);ax.text(.69,.38,f'Lineage-adjusted\nP = {pstr(r2.P)}',ha='center',fontsize=8)
rf=dd.iloc[-1];ax.text(.5,.17,f'Full complete-case APC effect: beta = {rf.beta:.3f}, P = {pstr(rf.P)}\nAbundance-independent APC effect is not established.',ha='center',fontsize=9,color='#C55358')
save(fig,'Fig2H_Working_model','Connecting lines indicate associations, not causal arrows. Isogenic experiments remain necessary.')

# Supplementary figures for the remaining supported analyses.
fig,aa=base('S1','Alternative Wnt scores',180,145,nrows=2,ncols=2);fig.subplots_adjust(top=.82,bottom=.12,hspace=.85,wspace=.40)
for j,(feat,label) in enumerate([('Wnt_PROGENy','PROGENy WNT weighted score (z)'),('Wnt_custom','Seven-target mean-z score')]):
    groupplot(aa[j,0],feat,label);scatter(aa[j,1],feat,label)
save(fig,'Fig2S1_Wnt_score_sensitivity','PROGENy: top 100 footprint targets, weighted standardized expression; custom: 7 prespecified genes.')

mr=pd.read_csv(T/'methylation_associations.csv');fig,aa=base('S2','TCF7L2 methylation associations',180,165,nrows=3,ncols=3);fig.subplots_adjust(top=.86,bottom=.10,hspace=.77,wspace=.50)
for i,reg in enumerate(meth):
    for j,(tar,label) in enumerate([('RNA26_TCF7L2','RNA (log-TPM)'),('MS_TCF7L2','Protein'),('TCF7L2_gene_effect','Gene effect')]):
        ax=aa[i,j];v=m.loc[m.APC_class.isin(['WT','Truncating/LOF'])].dropna(subset=[reg,tar])
        for g in ['WT','Truncating/LOF']:
            a=v.loc[v.APC_class.eq(g)];ax.scatter(a[reg],a[tar],s=5,c=colors[g],alpha=.45,linewidths=0)
        rr=mr.loc[mr.region.eq(reg)&mr.target.eq(tar)&mr.group.eq('WT_LOF')].iloc[0]
        ax.set_title(f'n={int(rr.n)}; rho={rr.rho:.2f}\nAdjusted FDR={rr.FDR:.3f}',fontsize=7)
        ax.set_xlabel(f'Region {i+1} methylation');ax.set_ylabel(label)
save(fig,'Fig2S2_Methylation_associations','Adjusted FDR: lineage-adjusted regional tests, BH across 27 region/target/stratum models. Red LOF; blue WT.')

fig,aa=base('S3','APC-effect sensitivity models',180,155,nrows=2,ncols=2);fig.subplots_adjust(left=.15,right=.91,top=.84,bottom=.12,wspace=.70,hspace=.65)
for ax,mode,title in zip(aa.flat,['All_available']+meth,['All available samples','Region 1 complete cases','Region 2 complete cases','Region 3 complete cases']):
    ff=fr.loc[fr['mode'].eq(mode)];forest(ax,ff,title+f'\nFull n={int(ff.iloc[-1].n)}, LOF={int(ff.iloc[-1].LOF)}')
save(fig,'Fig2S3_APC_model_sensitivity','All-available models use different samples; coefficient changes cannot be attributed to adjustment alone.')

fig,aa=base('S4','CTNNB1 protein as a molecular bridge',180,145,nrows=2,ncols=2);fig.subplots_adjust(top=.82,bottom=.12,hspace=.82,wspace=.4)
for j,(feat,label) in enumerate([('MS_CTNNB1','CTNNB1 mass spectrometry'),('RPPA_beta_catenin','Total beta-catenin RPPA')]):groupplot(aa[j,0],feat,label);scatter(aa[j,1],feat,label)
save(fig,'Fig2S4_CTNNB1_protein','Platforms analyzed separately. Total protein does not measure nuclear active beta-catenin.')

rg=pd.read_csv(T/'RNA_GSEA_beta.csv');pg=pd.read_csv(T/'Protein_GSEA_beta.csv')
fig,aa=base('S5','Dependency-associated RNA and protein pathways',180,145,nrows=1,ncols=2);fig.subplots_adjust(left=.28,right=.90,top=.81,bottom=.15,wspace=2.6)
gsplot(aa[0,0],rg,10);gsplot(aa[0,1],pg,10);aa[0,0].set_title('RNA Hallmark pathways',fontsize=8);aa[0,1].set_title('Protein Hallmark pathways',fontsize=8)
save(fig,'Fig2S5_RNA_protein_pathways','Rank = negative standardized beta; positive NES denotes expression linked to stronger dependency.')
ccg=rg.loc[rg.Term.str.startswith('HALLMARK_'),['Term','NES','FDR q-val']].merge(pg[['Term','NES','FDR q-val']],on='Term',suffixes=('_RNA','_Protein'));ccg.to_csv(T/'RNA_protein_Hallmark_concordance.csv',index=False)

fig,aa=base('S6','APC LOF and allelic-loss evidence',130,95);ax=aa[0,0];egs=['WT with CN','LOF; minor CN > 0','LOF; minor CN = 0'];rng=np.random.default_rng(23)
for i,g in enumerate(egs):
    v=m.loc[m.LOH_evidence.eq(g),'TCF7L2_gene_effect'].dropna();col=['#377DA8','#BE9348','#C55358'][i]
    ax.scatter(i+rng.uniform(-.2,.2,len(v)),v,s=7,c=col,alpha=.3,linewidths=0);ax.boxplot([v],positions=[i],widths=.5,showfliers=False,medianprops={'color':INK})
ax.set_xticks(range(3),[f'{g.replace("; ",chr(10))}\nn={m.LOH_evidence.eq(g).sum()}' for g in egs]);ax.set_ylabel('TCF7L2 gene effect');ax.set_title('Lineage-adjusted LOF vs WT: minor CN = 0, Holm P = 0.0004\nWithin LOF: minor CN = 0 vs > 0, Holm P = 0.331',fontsize=7.5)
save(fig,'Fig2S6_APC_LOH_evidence','Same-sequencing mutation/CN matches; phasing and absolute-CN quality require further validation.')

md=pd.read_csv(T/'exploratory_mediation.csv').iloc[0]
fig,aa=base('S7','Exploratory statistical mediation',135,85);ax=aa[0,0];ax.axis('off');ax.set_position([.04,.12,.92,.70])
schematic_box(ax,(.02,.55),.23,.26,'APC LOF',size=9);schematic_box(ax,(.385,.55),.23,.26,'Hallmark Wnt',size=9);schematic_box(ax,(.75,.55),.23,.26,'Gene effect',size=9)
arrow(ax,.26,.68,.375,.68);arrow(ax,.625,.68,.74,.68)
ax.text(.315,.45,f'a = {md.a:.3f}',ha='center',fontsize=8);ax.text(.685,.45,f'b = {md.b:.3f}',ha='center',fontsize=8)
ax.text(.5,.20,f'Indirect association = {md.indirect:.3f}\n95% bootstrap CI [{md.CI_low:.3f}, {md.CI_high:.3f}]\nn = {int(md.n)}; 2,000 lineage-stratified resamples',ha='center',fontsize=9)
save(fig,'Fig2S7_Exploratory_mediation','Path diagram is a statistical model; cross-sectional observational data do not identify causal mediation.')

crc=m.loc[m.lineage_2.eq('Colorectal Adenocarcinoma')];fig,aa=base('S8','CRC-restricted sensitivity analyses',180,145,nrows=2,ncols=2);fig.subplots_adjust(top=.82,bottom=.12,hspace=.85,wspace=.4)
groupplot(aa[0,0],'Wnt_Hallmark','Hallmark Wnt score (z)',df=crc)
for ax,feat,label in [(aa[0,1],'Wnt_Hallmark','Hallmark Wnt score (z)'),(aa[1,0],'RNA26_TCF7L2','TCF7L2 RNA'),(aa[1,1],'MS_TCF7L2','TCF7L2 protein')]:scatter(ax,feat,label,scope='CRC',df=crc.loc[crc.APC_class.isin(['WT','Truncating/LOF'])])
save(fig,'Fig2S8_CRC_sensitivity','CRC models contain no lineage covariate. Small WT/protein groups limit precision.')
pd.DataFrame(manifest).to_csv(O/'PDF_manifest_QA.csv',index=False)
print(pd.DataFrame(manifest).to_string(index=False),flush=True)
