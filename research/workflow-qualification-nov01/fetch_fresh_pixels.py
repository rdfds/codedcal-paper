"""Finite one-field original-TIFF qualification reader; no comparative fits."""
from pathlib import Path
import os,sys,json,time,hashlib,shutil,urllib.request,urllib.error,csv,io,re,fcntl
import numpy as np
import tifffile
H=Path(__file__).resolve().parent;CACHE=H.parents[1].parent/'workflow-qualification-fresh.nosync';CACHE.mkdir(exist_ok=True)
lock=(CACHE/'collector.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
state=CACHE/'state.json';assert not state.exists(),'Diagnose a prior attempt explicitly before resuming'
entries=json.loads((H/'public-metadata/metadata-1.txt').read_text())['data']['latestVersion']['files']
originals=json.loads((H/'fresh-originals/receipts.json').read_text())['records'];assert len(originals)==3 and all(r['publisher_md5_matches'] for r in originals)
organization=list(csv.DictReader(io.StringIO((H/'fresh-originals/12555511.bin').read_text())))
organization=[r for r in organization if r.get('readoutName','').startswith('RS')]
assert len(organization)==16 and sorted(int(r['bitNumber']) for r in organization)==list(range(1,17))
plan=[]
for roundno in range(1,10):
    name=f'merFISH_merged_{roundno:02}_000.tiff';e=next(x for x in entries if x['dataFile']['filename']==name);f=e['dataFile'];assert not e['restricted'] and f['filesize']==232717320
    plan.append(dict(round=roundno,name=name,id=f['id'],md5=f['checksum']['value'],bytes=f['filesize'],url=f"https://dataverse.harvard.edu/api/access/datafile/{f['id']}"))
(H/'fresh-pixel-pins.json').write_text(json.dumps(plan,indent=2)+'\n')
records=[];total=0;start=time.time()
def save(stage,error=None):
    d=dict(pid=os.getpid(),command=str(Path(__file__).resolve()),stage=stage,records=records,transport_bytes=total,started_epoch=start,error=error,comparative_outcomes_examined=False)
    tmp=state.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(state)
save('running')
try:
    for i,p in enumerate(plan):
        if i:time.sleep(10)
        assert shutil.disk_usage(CACHE).free>=1073741824+p['bytes'],'Disk reserve'
        row=dict(**p,stage='requested');records.append(row);save('running');part=CACHE/'current.tiff.partial';assert not part.exists()
        digest=hashlib.sha256();md=hashlib.md5();n=0
        with urllib.request.urlopen(urllib.request.Request(p['url'],headers={'User-Agent':'codedcal-public-qualification/0.1'}),timeout=180) as r:
            row.update(status=r.status,headers={k:v for k,v in r.headers.items() if k.lower() in ['content-length','content-type','etag','last-modified','retry-after']})
            assert r.status==200 and int(r.headers['Content-Length'])==p['bytes'] and not r.headers.get('Retry-After')
            with part.open('xb') as out:
                while block:=r.read(min(1048576,p['bytes']-n+1)):
                    n+=len(block);total+=len(block);out.write(block);digest.update(block);md.update(block)
                    assert n<=p['bytes'] and total<=2094455880
                    assert time.time()-start<=3600,'Finite worker time exceeded'
        row.update(bytes_received=n,sha256=digest.hexdigest(),actual_md5=md.hexdigest());save('running')
        assert n==p['bytes'] and md.hexdigest()==p['md5'],'Original TIFF length/MD5 mismatch'
        with tifffile.TiffFile(part) as t:
            assert len(t.pages)==45 and t.pages[0].shape==(1608,1608)
            arrays={};mapping=[]
            for rna in organization:
                if int(rna['imagingRound'])!=p['round']:continue
                frames=[int(v) for v in re.findall(r'\d+',rna['frame'])];assert len(frames)==15
                frame=frames[7];a=t.pages[frame].asarray();a=a.T[::-1,::-1];a=a[676:932,676:932].copy()
                key='bit'+str(int(rna['bitNumber']));arrays[key]=a;mapping.append(dict(bit=int(rna['bitNumber']),readout=rna['readoutName'],color=int(rna['color']),z_index=7,frame_index=frame,shape=list(a.shape)))
            fid=t.pages[6].asarray().T[::-1,::-1][676:932,676:932].copy();arrays['fiducial_frame6']=fid
            dest=CACHE/f'round-{p["round"]:02}.npz';np.savez_compressed(dest,**arrays)
            row.update(tiff_pages=45,mapping=mapping,extraction=str(dest),extraction_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),publisher_md5_verified=True,stage='verified_and_extracted')
        save('running');part.unlink();print(p['round'],row['stage'],len(mapping),total,flush=True)
    save('raw_input_retained')
except Exception as e:
    if isinstance(e,urllib.error.HTTPError):records[-1].update(status=e.code,retry_after=e.headers.get('Retry-After'))
    save('stopped_with_error',type(e).__name__+': '+str(e));raise
