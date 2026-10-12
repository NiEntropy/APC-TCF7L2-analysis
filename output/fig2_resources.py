from pathlib import Path
import urllib.request,json,hashlib
O=Path(__file__).resolve().parent/'Fig2_updated';R=O/'resources';R.mkdir(parents=True,exist_ok=True)
urls={f'{prefix}.v2025.1.Hs.symbols.gmt':f'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/{prefix}.v2025.1.Hs.symbols.gmt' for prefix in ['h.all','c2.cp.reactome','c2.cp.kegg_medicus']}
urls['PROGENy_omnipath.tsv']='https://omnipathdb.org/annotations?databases=PROGENy&format=tsv&license=academic'
log=[]
for name,url in urls.items():
    try:
        p=R/name
        if not p.exists():
            with urllib.request.urlopen(url,timeout=50) as r:p.write_bytes(r.read())
        log.append(dict(file=name,url=url,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        print(name,p.stat().st_size,flush=True)
    except Exception as e:log.append(dict(file=name,url=url,error=str(e)));print(name,str(e),flush=True)
(R/'sources.json').write_text(json.dumps(log,indent=2))
