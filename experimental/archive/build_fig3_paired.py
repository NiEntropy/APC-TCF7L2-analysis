from pathlib import Path
import xml.etree.ElementTree as ET,zipfile,json,csv,sys,shutil
sys.path.insert(0,str(Path.cwd()/'tmp/pdfs/python_libs'))
import numpy as np,openpyxl,pymupdf as fitz
from scipy.stats import ttest_rel
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0,str(Path.cwd()/'tmp/pdfs/python_libs'))
import xlrd
ROOT=Path.cwd(); SRC=ROOT/'2026 崔珣 合成致死课题原始数据/Figure 4'; OUT=ROOT/'output/pdf/Fig3_paired_statistics';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,'xtick.major.size':2,'ytick.major.size':2,'axes.labelpad':2})
colors=['#777D83','#8A3048']; cells=['HT29','SW620']; sources={}; rows=[]
flow={'HT29':{'early':np.array([[2.90,3.72,3.35],[8.89,7.91,8.38]]),'late':np.array([[4.71,5.16,4.47],[11.3,12.3,12.5]])},'SW620':{'early':np.array([[5.04,4.93,4.59],[11.1,13.3,12.4]]),'late':np.array([[1.45,1.26,1.40],[4.25,5.18,3.55]])}}
for c in cells:flow[c]['total']=flow[c]['early']+flow[c]['late']
flowpdf={'HT29':next((SRC/'Figure 4-2/HT29流式数据/细胞凋亡').glob('*Layout.pdf')),'SW620':next((SRC/'Figure 4-2/SW620流式数据/细胞凋亡').rglob('SW620 凋亡流式.pdf'))}
# Native export rectangles include the unmodified density plot, axes and quadrant labels.
clips={'HT29':[(70,102,230,257),(234,103,395,258),(395,103,558,258),(77,304,237,460),(242,304,403,460),(408,304,570,460)],'SW620':[(75,96,235,243),(239,96,400,243),(403,96,565,243),(75,292,235,447),(239,292,400,447),(403,292,565,447)]}
gsea={}
for c in cells:
 zpath=next((SRC/'Figure 4-1'/f'{c} RNA-sec GSEA分析').glob('*.zip'))
 with zipfile.ZipFile(zpath) as z:
  rank=np.array([float(l.split('\t')[1]) for l in z.read(next(n for n in z.namelist() if n.endswith('express.rnk'))).decode().splitlines()])
  edb=ET.fromstring(z.read(next(n for n in z.namelist() if n.endswith('results.edb'))))
  dtg=next(x for x in edb.iter('DTG') if 'GO:0097190' in x.get('GENESET',''))
  hits=np.array(dtg.get('HIT_INDICES').split(),int); a=np.full(len(rank),-1/(len(rank)-len(hits)));a[hits]=abs(rank[hits])/abs(rank[hits]).sum();running=np.cumsum(a)
  assert abs(running.max()-float(dtg.get('ES')))<.002,(c,running.max(),dtg.get('ES'))
  assert len(dtg.get('RND_ES','').split()) == 1000 and float(dtg.get('NP')) == 0
  gsea[c]=(rank,hits,running,{k:float(dtg.get(k)) for k in ['ES','NES','NP','FDR','FWER']})
  sources[c]={'gsea':str(zpath.relative_to(ROOT)),'flow':str(flowpdf[c].relative_to(ROOT))}
wb={};qpcr={};raw={};crops={}
b= openpyxl.load_workbook(next((SRC/'Figure 4-3/3 凋亡相关WB检测').rglob('*.xlsx')),data_only=True).active
for c,start in [('HT29',3),('SW620',12)]:
 vals={k:np.array([b.cell(i,col).value for i in range(start,start+6)],float) for k,col in [('GAPDH',3),('BAX',6),('BCL-2',11),('Caspase-3',19),('Cleaved caspase-3',24)]}
 wb[c]={}
 for k in ['BCL-2','BAX','Caspase-3','Cleaved caspase-3']:
  v=vals[k]/vals['GAPDH'];wb[c][k]=(v/v[:3].mean()).reshape(2,3)
 for k,n,d in [('BAX / BCL-2','BAX','BCL-2'),('Cleaved / total caspase-3','Cleaved caspase-3','Caspase-3')]:
  v=vals[n]/vals[d];wb[c][k]=(v/v[:3].mean()).reshape(2,3)
 qpath=next((SRC/'Figure 4-3/2 热图相关凋亡因子qPCR检测'/c).glob('*.xls'));s=xlrd.open_workbook(qpath).sheet_by_name('Sheet1');qpcr[c]={}
 for startq in range(1,78,13):
  gene=s.cell_value(startq,5).replace('H-','');qpcr[c][gene]=np.array([s.cell_value(startq+j,12) for j in range(6)]).reshape(2,3)
 sources[c]['qPCR']=str(qpath.relative_to(ROOT));raw[c]={}
 for sub in (SRC/'Figure 4-3/3 凋亡相关WB检测'/c).iterdir():
  if sub.is_dir() and list(sub.glob('*signal*.jpg')):raw[c][sub.name]=list(sub.glob('*signal*.jpg'))[0]
crops={'HT29':{'Bax':[530,482,834,525],'Bcl-2':[545,479,845,520],'Caspase3':[457,535,759,581],'Cleaved-Caspase3':[378,590,681,630],'GAPDH':[562,516,868,566]},'SW620':{'Bax':[423,473,741,513],'Bcl2':[558,510,860,555],'Caspase3':[457,533,768,590],'Cleaved-Caspase3':[609,623,1013,668],'GAPDH':[432,525,742,575]}}
for typ,data in [('flow',flow),('WB',wb),('qPCR',qpcr)]:
 for c,metrics in data.items():
  for metric,arr in metrics.items():
   for gi,g in enumerate(['shNC','shTCF7L2']):
    for j,v in enumerate(arr[gi]):rows.append([typ,c,metric,g,j+1,float(v)])
with (OUT/'plotted_data.csv').open('w',newline='',encoding='utf-8-sig') as f:
 w=csv.writer(f);w.writerow(['assay','cell','metric','group','record','value']);w.writerows(rows)
(OUT/'source_map.json').write_text(json.dumps({'sources':sources,'WB_crops_pixels':crops,'flow_clips_points':clips,'GSEA':{c:gsea[c][3] for c in cells}},ensure_ascii=False,indent=2),encoding='utf-8')

statrows=[];stars={}
for panelname,metric,data in [('D','total',flow),('F','BAX / BCL-2',wb),('F','Cleaved / total caspase-3',wb)]:
 for c in cells:
  a=data[c][metric];test=ttest_rel(a[1],a[0]);pval=float(test.pvalue)
  star='****' if pval<.0001 else '***' if pval<.001 else '**' if pval<.01 else '*' if pval<.05 else 'ns'
  stars[(c,metric)]=star
  diff=a[1]-a[0]
  statrows.append({'panel':panelname,'cell':c,'metric':metric,'n_pairs':3,'test':'Two-sided paired t-test','t_KD_minus_NC':float(test.statistic),'df':2,'p_unadjusted':pval,'significance':star,'mean_NC':float(a[0].mean()),'mean_KD':float(a[1].mean()),'mean_paired_difference':float(diff.mean())})
with (OUT/'Fig3D_F_paired_statistics.csv').open('w',newline='',encoding='utf-8-sig') as sf:
 writer=csv.DictWriter(sf,fieldnames=statrows[0].keys());writer.writeheader();writer.writerows(statrows)
print(json.dumps(statrows,ensure_ascii=False,indent=2))

def txt(fig,x,y,text,**kw):fig.text(x,y,text,**kw)
def panel(fig,x,y,letter,title=''):
 txt(fig,x,y,letter,fontsize=11,fontweight='bold');txt(fig,x+.033,y,title,fontsize=8,fontweight='bold')
def chart(ax,arrs,labels,ylabel,stars=None):
 maxv=max(float(np.max(a)) for a in arrs);n=len(arrs)
 for i,a in enumerate(arrs):
  for g in range(2):
   x=i+(g-.5)*.28;v=a[g];ax.bar(x,v.mean(),width=.25,color=colors[g],alpha=.32,edgecolor=colors[g],linewidth=.7)
   ax.errorbar(x,v.mean(),yerr=v.std(ddof=1)/np.sqrt(len(v)),fmt='none',color=colors[g],capsize=2,lw=.8)
   ax.scatter(x+np.linspace(-.06,.06,len(v)),v,s=13,color=colors[g],edgecolor='white',lw=.25,zorder=4)
  if stars and stars[i]:
   y=a.max()+maxv*.10;ax.plot([i-.15,i-.15,i+.15,i+.15],[y,y+maxv*.025,y+maxv*.025,y],color='black',lw=.6);ax.text(i,y+maxv*.035,stars[i],ha='center',fontsize=8)
 ax.set(xticks=range(n),xticklabels=labels,ylabel=ylabel,ylim=(0,maxv*1.35),xlim=(-.6,n-.4));ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=6.5)
def legend(fig,x=.74,y=.97):
 for i,g in enumerate(['shNC','shTCF7L2']):txt(fig,x+i*.115,y,'● '+g,color=colors[i],fontsize=7)
def save(fig,name,overlays=[],append=None):
 temp=OUT/(name+'_base.pdf');fig.savefig(temp);plt.close(fig);doc=fitz.open(temp);p=doc[0]
 for source,clip,box in overlays:
  sd=fitz.open(source); x,y,w,h=box; rect=fitz.Rect(x*p.rect.width,(1-y-h)*p.rect.height,(x+w)*p.rect.width,(1-y)*p.rect.height);p.show_pdf_page(rect,sd,0,clip=fitz.Rect(clip))
 if append:
  for q in append:doc.insert_pdf(fitz.open(q))
 doc.save(OUT/(name+'.pdf'),garbage=4,deflate=True);doc.close();temp.unlink()
 final=fitz.open(OUT/(name+'.pdf'))
 for i,p in enumerate(final):
  suffix='' if len(final)==1 else f'_p{i+1}'
  (OUT/(name+suffix+'.svg')).write_text(p.get_svg_image(text_as_path=False),encoding='utf-8')
  p.get_pixmap(matrix=fitz.Matrix(1.6,1.6)).save(str(ROOT/'tmp/pdfs'/f'{name}{suffix}.png'))
 return final
# FIGURE 3: 180 x 205 mm, three evidence rows.
f=plt.figure(figsize=(180/25.4,205/25.4));legend(f,.72,.982)
for ci,c in enumerate(cells):
 x=.09+ci*.49;panel(f,x-.06,.947,chr(65+ci),c+' | Apoptotic signaling pathway')
 txt(f,x,.93,'GO:0097190',fontsize=6,color='#555555');rank,hits,es,stats=gsea[c]; ax=f.add_axes([x,.80,.385,.115]);ax.plot(np.arange(len(es)),es,color=colors[1],lw=1);ax.axhline(0,color='#999999',lw=.5);ax.set(xlim=(0,len(es)),ylabel='Enrichment score');ax.set_xticks([]);ax.spines[['top','right']].set_visible(False)
 ax.text(.98,.90,f"NES = {stats['NES']:.4f}\nNominal P < 0.001\nFDR q = {stats['FDR']:.4f}",transform=ax.transAxes,ha='right',va='top',fontsize=7)
 a=f.add_axes([x,.775,.385,.017]);a.vlines(hits,0,1,color='#444444',lw=.3);a.set(xlim=(0,len(es)));a.axis('off')
 a=f.add_axes([x,.735,.385,.035]);a.fill_between(np.arange(len(rank)),rank,color='#A8A8A8',lw=0);a.set(xlim=(0,len(rank)),xlabel='Gene rank',yticks=[]);a.tick_params(labelsize=6);a.spines[['top','right','left']].set_visible(False)
 txt(f,x,.708,'shTCF7L2-enriched',fontsize=6,color=colors[1]);txt(f,x+.27,.708,'shNC-enriched',fontsize=6,color=colors[0])
panel(f,.03,.675,'C','Annexin V-APC / DAPI');panel(f,.66,.675,'D','Annexin V-positive cells')
ov=[]
for ri,c in enumerate(cells):
 y=.49-ri*.17;txt(f,.045,y+.11,c,rotation=90,fontsize=8,fontweight='bold')
 for gi,idx in enumerate(([2,4] if c=='HT29' else [1,5])):
  x=.075+gi*.275;ov.append((flowpdf[c],clips[c][idx],(x,y,.245,.157)))
  if ri==0:txt(f,x+.122,y+.158,['shNC','shTCF7L2'][gi],ha='center',fontsize=7,color=colors[gi])
chart(f.add_axes([.72,.385,.235,.23]),[flow[c]['total'] for c in cells],cells,'Annexin V+ cells (%)',[stars[(c,'total')] for c in cells])
# Main WB rows retain all six lanes, no lane splicing.
panel(f,.03,.285,'E','Apoptosis-associated proteins');panel(f,.64,.285,'F','Protein ratios')
proteins=[('BCL-2','Bcl-2','Bcl2'),('BAX','Bax','Bax'),('Cleaved\ncaspase-3','Cleaved-Caspase3','Cleaved-Caspase3'),('Caspase-3','Caspase3','Caspase3'),('GAPDH','GAPDH','GAPDH')]
for ci,c in enumerate(cells):
 x=.155+ci*.23;txt(f,x+.09,.257,c,ha='center',fontweight='bold',fontsize=7)
 txt(f,x+.043,.242,'shNC',ha='center',fontsize=6);txt(f,x+.14,.242,'shTCF7L2',ha='center',fontsize=6)
 for pi,(label,k1,k2) in enumerate(proteins):
  key=k1 if ci==0 else k2;y=.206-pi*.036
  a=f.add_axes([x,y,.19,.026]);a.imshow(Image.open(raw[c][key]).crop(crops[c][key]),cmap='gray',aspect='auto');a.axis('off')
  if ci==0:txt(f,.035,y+.008,label,fontsize=6,va='center')
for i,key in enumerate(['BAX / BCL-2','Cleaved / total caspase-3']):
 ax=f.add_axes([.69+i*.17,.085,.135,.145]);chart(ax,[wb[c][key] for c in cells],cells,'NC-normalized ratio' if i==0 else '',[stars[(c,key)] for c in cells]);ax.set_title('BAX / BCL-2' if i==0 else 'Cleaved / total\ncaspase-3',fontsize=6.5)
txt(f,.03,.018,'n = 3 biological replicates; mean ± SEM. D, F: two-sided paired t-tests (unadjusted P).',fontsize=6)
save(f,'Fig3',ov)

(OUT/'统计方法说明.md').write_text("""# Fig3D、F 配对显著性分析

用户已确认：每组 3 个记录为独立生物学重复，NC 与 shTCF7L2 按批次一一配对。按原始记录顺序 1 对 1、2 对 2、3 对 3，进行双侧配对 t 检验，n=3 对，df=2；图示均值±SEM。

D：逐样本总 Annexin V+ = Q2 + Q3。F：同一记录的 BAX/BCL-2 和 cleaved/total caspase-3 比值，除以各自 NC 组均值归一化后检验；共同缩放不改变配对 t 检验的 P 值。

六项比较显示未校正 P 值；未做多重比较校正，星号仅表示各项比较的名义显著性。* P<0.05，** P<0.01，*** P<0.001，**** P<0.0001；否则 ns。配对 t 检验假定不同批次相互独立、批次内差值近似正态；n=3 无法充分评估正态性。

所有星号均重新计算，未沿用旧图星号。仅更新 Fig3D、F 及其图注；S7、S8 未重新开展统计。详细数值见 Fig3D_F_paired_statistics.csv。原版 Fig3 保留于 Fig3_S7_S8 目录，本目录为更新版本。
""",encoding='utf-8')
shutil.copy2(__file__,OUT/'build_fig3_paired.py')

(OUT/'GSEA_nominal_P说明.md').write_text('Fig3A（HT29）与 Fig3B（SW620）原始 results.edb 的 NP 均为 0.0000，各含 1000 个置换 ES。按照 GSEA 官方说明，图中显示 Nominal P < 0.001，而不是 P = 0。NES、FDR 及 D/F 配对统计保持不变。\n参考：http<LOCAL_PATH_REMOVED>',encoding='utf-8')
