from pathlib import Path
import sys,ast,pickle,shutil,json
ROOT=Path.cwd();sys.path.insert(0,str(ROOT/'tmp/pdfs/python_libs'))
s=(ROOT/'tmp/pdfs/build_fig45.py').read_text(encoding='utf-8')
exec(s[:s.index('ROOT=Path.cwd();')])
OUT=ROOT/'output/pdf/Fig5_BC_rasterized_20261008';OUT.mkdir(exist_ok=True)
DATA=OUT/'data';DATA.mkdir(exist_ok=True)
CELLS=['HT29','SW620'];NC='#777D83';KD='#8A3048';BLUE='#3978A4';COLORS=[NC,KD];FIGURES={}
D=pickle.loads((ROOT/'tmp/pdfs/fig45_data.pkl').read_bytes())
plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,'xtick.major.size':2,'ytick.major.size':2,'axes.labelpad':3,'savefig.facecolor':'white'})
for node in ast.parse(s).body:
    if isinstance(node,ast.FunctionDef):exec(compile(ast.Module(body=[node],type_ignores=[]),'<helpers>','exec'))
block=s[s.index('# Figure 5: seven panels.'):s.index('# S10: source PCA')]
latest=(ROOT/'tmp/pdfs/revise_fig5g.py').read_text(encoding='utf-8')
df=pd.read_csv(ROOT/'output/pdf/Fig5G_n6_20261008/data/Fig5G_independent_samples.csv')
result=pd.read_csv(ROOT/'output/pdf/Fig5G_n6_20261008/data/Fig5G_Student_t.csv')
genes=result.gene.tolist();arrays={g:[df[(df.gene==g)&(df.group==group)].expression.to_numpy() for group in ['shNC','shTCF7L2']] for g in genes}
start=block.index("ax=f.add_axes([.65,.17,.31,.17])")
end=block.index("text(f,.045,.066",start)
g=latest[latest.index('for j,gene in enumerate(genes):'):latest.index("tmp=ROOT/'tmp/pdfs/Fig5G_overlay.pdf'")]
block=block[:start]+"ax=f.add_axes([.65,.17,.31,.17]);colors=COLORS\n"+g+block[end:]
SAMPLECOLS=['shNC_1','shNC_2','shNC_3','shTCF7L2_1','shTCF7L2_2','shTCF7L2_3']
exec(block)
# Rasterize only volcano point collections, at 600 dpi. Other artists remain vector.
for ax in f.axes[:2]:
    assert len(ax.collections)==3
    for collection in ax.collections:collection.set_rasterized(True)
f.savefig(OUT/'Fig5.pdf',dpi=600)
f.savefig(OUT/'Fig5.svg',dpi=600)
plt.close(f)
doc=fitz.open(OUT/'Fig5.pdf');p=doc[0]
p.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(str(OUT/'Fig5_preview.png'))
assert 'HLA-B' in p.get_text() and '6 independent biological samples' in p.get_text()
images=p.get_images(full=True);assert len(images)>=2
(OUT/'rendering_notes.txt').write_text('Panels B and C: scatter collections rasterized at 600 dpi. Axes, labels, annotations, threshold lines and other panels remain vector. Panel G retains n=6 independent biological samples/group, two-sided unpaired Student t tests and unadjusted P values.\n'+str([(i[2],i[3]) for i in images]),encoding='utf-8')
for file in (ROOT/'output/pdf/Fig5G_n6_20261008/data').glob('*.csv'):shutil.copy2(file,DATA/file.name)
shutil.copy2(__file__,OUT/'reproduce.py')
print({'images':[(i[2],i[3]) for i in images],'PDF_bytes':(OUT/'Fig5.pdf').stat().st_size})

