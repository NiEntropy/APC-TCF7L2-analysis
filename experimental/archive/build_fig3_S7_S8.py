from pathlib import Path
import xml.etree.ElementTree as ET,zipfile,json,csv,sys,shutil
import numpy as np,openpyxl,pymupdf as fitz
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0,str(Path.cwd()/'tmp/pdfs/python_libs'))
import xlrd
ROOT=Path.cwd(); SRC=ROOT/'2026 崔珣 合成致死课题原始数据/Figure 4'; OUT=ROOT/'output/pdf/Fig3_S7_S8';OUT.mkdir(parents=True,exist_ok=True)
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
  gsea[c]=(rank,hits,running,{k:float(dtg.get(k)) for k in ['ES','NES','FDR','FWER']})
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
 ax.text(.98,.90,f"NES = {stats['NES']:.4f}\nFDR q = {stats['FDR']:.4f}",transform=ax.transAxes,ha='right',va='top',fontsize=7)
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
chart(f.add_axes([.72,.385,.235,.23]),[flow[c]['total'] for c in cells],cells,'Annexin V+ cells (%)')
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
 ax=f.add_axes([.69+i*.17,.085,.135,.145]);chart(ax,[wb[c][key] for c in cells],cells,'NC-normalized ratio' if i==0 else '',None if i==0 else ['*','*']);ax.set_title('BAX / BCL-2' if i==0 else 'Cleaved / total\ncaspase-3',fontsize=6.5)
txt(f,.03,.018,'Points: recorded samples / lanes; mean ± SEM. Original significance retained for the cleaved/total ratio.',fontsize=6)
save(f,'Fig3',ov)
# S7: all final quadrant exports and corrected early / late outcomes.
f=plt.figure(figsize=(210/25.4,270/25.4));panel(f,.035,.965,'A','HT29 | All exported samples');legend(f,.73,.97);ov=[]
for ci,c in enumerate(cells):
 top=.89-ci*.40
 if ci:panel(f,.035,top+.075,'B','SW620 | All exported samples')
 for j in range(6):
  col=j%3;rr=j//3;x=.09+col*.295;y=top-rr*.17-.13
  ov.append((flowpdf[c],clips[c][j],(x,y,.255,.158)));txt(f,x+.125,y+.158,('shNC' if j<3 else 'shTCF7L2')+f' | sample {col+1}',ha='center',fontsize=7,color=colors[rr])
for ci,c in enumerate(cells):
 panel(f,.045+ci*.49,.175,chr(67+ci),c+' | Apoptosis fractions')
 chart(f.add_axes([.10+ci*.49,.055,.35,.10]),[flow[c]['early'],flow[c]['late']],['Early (APC+ / DAPI-)','Late (APC+ / DAPI+)'],'Cells (%)',['****','****'] if ci==0 else ['****','**'])
save(f,'FigS7_quadrants',ov)
f=plt.figure(figsize=(210/25.4,180/25.4));panel(f,.05,.94,'E','Gating hierarchy and source consistency')
for ci,c in enumerate(cells):
 y=.82-ci*.31;txt(f,.06,y,c,fontweight='bold',fontsize=10)
 labels=['All events','Cells / lymphocytes','Single cells']+(['GFP / SSC subset'] if c=='SW620' else [])+['Annexin V-APC / DAPI']
 for j,label in enumerate(labels):
  x=.055+j*.186;txt(f,x,y-.085,label,fontsize=7,bbox=dict(boxstyle='round,pad=.5',fc='#F2F2F2',ec='#999999')); 
  if j<len(labels)-1:txt(f,x+.145,y-.085,'→',fontsize=11)
 txt(f,.06,y-.155,'Hierarchy transcribed from the saved FlowJo workspace; boxes are schematic, not event plots.',fontsize=7)
 msg='HT29: saved workspace counts differ from the final exported plots.\nThe export is used for Fig3 and S7; matching upstream gate plots cannot be verified.' if c=='HT29' else 'SW620: saved workspace population counts agree with the exported quadrant percentages.\nQ3 = early apoptosis; Q2 = late apoptosis. Original Prism column titles were reversed.'
 txt(f,.06,y-.23,msg,fontsize=8,linespacing=1.5)
txt(f,.06,.12,'Quadrant percentages and total Annexin V+ values are tabulated in plotted_data.csv.\nTotal Annexin V+ = Q2 + Q3; no original star is transferred to this newly summed outcome.\nThe three observations per group are recorded samples, not verified independent biological repeats.',fontsize=8,linespacing=1.6)
save(f,'FigS7_gating')
d=fitz.open(OUT/'FigS7_quadrants.pdf');d.insert_pdf(fitz.open(OUT/'FigS7_gating.pdf'));d.save(OUT/'FigS7.pdf');d.close()
# S8: individual proteins and qPCR from cached numeric analysis tables.
f=plt.figure(figsize=(210/25.4,230/25.4));legend(f,.73,.97)
keys=['BCL-2','BAX','Caspase-3','Cleaved caspase-3'];genes=['CDKN1A','TNFSF10','TNF','DIXDC1','SOX12','MYB']
for ci,c in enumerate(cells):
 panel(f,.04,.91-ci*.42,chr(65+ci*2),c+' | Individual protein abundance')
 chart(f.add_axes([.10,.71-ci*.42,.36,.16]),[wb[c][k] for k in keys],['BCL-2','BAX','Total C3','Cleaved C3'],'Protein / GAPDH (NC = 1)',['ns','ns','*','*'] if ci==0 else ['ns','*','ns','****'])
 panel(f,.53,.91-ci*.42,chr(66+ci*2),c+' | Apoptosis-associated transcripts')
 chart(f.add_axes([.59,.71-ci*.42,.36,.16]),[qpcr[c][g] for g in genes],genes,'Relative mRNA expression',['****','***','****','**','*','ns'] if ci==0 else ['****','ns','****','****','****','****'])
 for ax in f.axes[-1:]:ax.tick_params(axis='x',labelrotation=45,labelsize=6)
txt(f,.06,.16,'Points show three recorded values per group; bars show mean ± SEM.\nStars / ns are transcribed from the matching original analysis figures; tests were not rerun.\nqPCR uses the cached, sample-averaged normalized values in the original workbook.\nC3: caspase-3. Full source WB signal images follow on pages 2 and 3.',fontsize=8,linespacing=1.7)
save(f,'FigS8_quantification')
full=[]
for ci,c in enumerate(cells):
 f=plt.figure(figsize=(210/25.4,297/25.4));panel(f,.04,.965,chr(69+ci),c+' | Full source WB signal images')
 for j,(label,k1,k2) in enumerate(proteins):
  key=k1 if ci==0 else k2;x=.045+(j%2)*.49;y=.67-(j//2)*.30
  a=f.add_axes([x,y,.43,.255]);a.imshow(Image.open(raw[c][key]));a.axis('off');a.set_title(label.replace('\n',' '),fontsize=9)
  txt(f,x,y-.012,raw[c][key].name,fontsize=5.5)
 txt(f,.55,.10,'Full signal files, uncropped.\nMain-figure strips use contiguous\nsix-lane crops from these images.\nLane order: shNC × 3, shTCF7L2 × 3.\nNo brightness / contrast adjustment.',fontsize=8,linespacing=1.6)
 name=f'FigS8_full_{c}';save(f,name);full.append(OUT/(name+'.pdf'))
d=fitz.open(OUT/'FigS8_quantification.pdf')
for p in full:d.insert_pdf(fitz.open(p))
d.save(OUT/'FigS8.pdf');d.close()
# Export unified SVG pages and previews, and keep only three final PDFs.
for name in ['FigS7','FigS8']:
 d=fitz.open(OUT/(name+'.pdf'))
 for i,p in enumerate(d):
  (OUT/f'{name}_p{i+1}.svg').write_text(p.get_svg_image(text_as_path=False),encoding='utf-8');p.get_pixmap(matrix=fitz.Matrix(1.4,1.4)).save(str(ROOT/'tmp/pdfs'/f'{name}_p{i+1}.png'))
for stem in ['FigS7_quadrants','FigS7_gating','FigS8_quantification','FigS8_full_HT29','FigS8_full_SW620']:
 for ext in ['.pdf','.svg']:(OUT/(stem+ext)).unlink(missing_ok=True)
notes='''# Fig3、S7、S8 交付说明\n\n- PDF 和 SVG 中的图表、文字、坐标与统计符号保留矢量；WB 信号图和流式密度点阵保留原始位图，不伪造为矢量实验数据。\n- Fig3：两个细胞系相同 GO:0097190 的 GSEA，流式代表图与总 Annexin V+，五种 WB 蛋白与两个比值。代表流式图按各组总 Annexin V+ 中位数选取。\n- GSEA 的 FDR 为 HT29 0.0831、SW620 0.0722，均不小于 0.05；不标显著性星号。曲线按原 ranked list 和 hit indices 重建，并核对 ES。\n- 总凋亡率采用原导出图 Q2+Q3；BAX/BCL-2 按同一记录的原始灰度先作比值，再除以 NC 组均值。两者不直接继承原星号。\n- S7 共 2 页：全部 12 个导出四象限图、早晚期定量、门控层级及版本说明。HT29 WSP 与已导出图计数不一致，无法把该工作区上游门控当作导出图的匹配门控。S7 第二页是门控示意，不是完整实测门控散点图；投稿前仍需补齐匹配版本的门控导出。\n- SW620 原 Prism 早晚期列名颠倒，按 Q3（APC+/DAPI-）=早期、Q2（APC+/DAPI+）=晚期纠正。实际染料为 DAPI，不是 PI。\n- S8 共 3 页：单蛋白、六种转录本和两细胞系完整 WB 信号文件。沿用源工作簿缓存的 qPCR 归一化结果，不重新处理 Ct。\n- 星号与 ns 来自对应原 PDF，未重新开展检验，原检验方法/独立生物学重复数仍需作者核对；3 个记录不宣称为 3 次独立实验。* P<0.05, ** P<0.01, *** P<0.001, **** P<0.0001 为通常定义，精确阈值需以原统计设置为准。\n- HT29 流式文件名含 shTCF7L2-6，而 WB/qPCR 分析图含 shTCF7L2-1；SW620 部分表格样本名也有不一致。图中统一写 shTCF7L2，不将其擅自解释为相同构建。\n- 同批凋亡 WB 未提供 TCF7L2 条带，未从其他实验拼接。无 LOVO 凋亡数据。\n- 原始输入未修改；裁剪坐标见 source_map.json，所有绘制数值见 plotted_data.csv。\n'''
(OUT/'说明与待核对事项.md').write_text(notes,encoding='utf-8')
shutil.copy2(__file__,OUT/'build_fig3.py')
for name in ['Fig3','FigS7','FigS8']:
 d=fitz.open(OUT/(name+'.pdf'));print(name,'pages',len(d),'text',sum(len(p.get_text()) for p in d),'vector_paths',sum(len(p.get_drawings()) for p in d))
print(OUT)


