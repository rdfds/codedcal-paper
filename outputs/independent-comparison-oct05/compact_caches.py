"""Retain exact calls at every scored threshold, remove verified ignored tensors.

Raw author HDFs, all scientific parameters and results remain. Recreate density
caches with the frozen runners to repeat tensor-level checks in analyze.py.
This changes storage, not the frozen scientific implementation or outcomes.
"""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.ndimage import maximum_filter

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,obj):
    p=Path(p);temp=p.with_suffix(p.suffix+'.tmp')
    temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temp.replace(p)
def portable(path,xy,words):
    out=io.StringIO();w=csv.writer(out,lineterminator='\n');w.writerow(['x','y','word']);w.writerows(zip(xy[:,0],xy[:,1],words))
    path.write_bytes(gzip.compress(out.getvalue().encode(),mtime=0))

def main():
    ledgerpath=HERE/'cache-projection-ledger.json'
    ledger=read(ledgerpath) if ledgerpath.exists() else []
    prior={a['density_path'] for a in ledger}
    for item in read(HERE/'density-and-call-manifest.json'):
        if item['path'] in prior:continue
        path=HERE/item['path'];primary=HERE/item['calls_path']
        assert subprocess.run(['git','check-ignore','--quiet',str(path.relative_to(ROOT))],cwd=ROOT).returncode==0
        assert sha(path)==item['sha256'] and sha(primary)==item['calls_sha256']
        condition,policy,decoder=path.stem.split('__')
        result=read(path.parent/(condition+'.json'))
        arm=result['arms'][policy+'__'+decoder]
        evidence=dict(density_path=item['path'],density_file_sha256=item['sha256'],
            primary_calls_path=item['calls_path'],primary_calls_sha256=item['calls_sha256'],thresholds={},bytes_removed=path.stat().st_size)
        with np.load(path) as data:
            density=np.asarray(data['density'],np.float32)
            maxima=maximum_filter(density,size=(1,3,3,1),mode='constant',cval=np.finfo(np.float32).min)
            for factor in [.5,1.,2.]:
                cut=np.float32(arm['threshold']*factor*(float(density.max()) if decoder=='regression' else 1.))
                z,y,x,j=np.where((density==maxima)&(maxima>cut));order=np.lexsort((x,y,z,j))
                xy=np.c_[x[order],y[order]].astype(float);words=j[order]
                inside=np.all((xy>=24)&(xy<176),axis=1)
                expected=arm['metrics']['5']['calls'] if factor==1 else arm['threshold_sensitivity'][str(factor)]['calls']
                assert int(inside.sum())==expected
                target=primary if factor==1 else path.with_suffix(f'.calls-factor-{factor}.csv.gz')
                if factor==1:
                    assert np.array_equal(data['xy'],xy) and np.array_equal(data['words'],words)
                    with gzip.open(primary,'rt') as f:rows=list(csv.DictReader(f))
                    existing=np.array([[float(r['x']),float(r['y']),int(r['word'])] for r in rows]).reshape(-1,3)
                    assert np.array_equal(existing,np.c_[xy,words])
                else:portable(target,xy,words)
                evidence['thresholds'][str(factor)]=dict(path=str(target.relative_to(HERE)),sha256=sha(target),all_image_calls=len(xy),interior_calls=int(inside.sum()),absolute_float32_cut=float(cut))
            evidence['density_max']=float(density.max())
            evidence['density_array_sha256']=hashlib.sha256(density.tobytes()).hexdigest()
        ledger.append(evidence);write(ledgerpath,ledger)
        # Every output is written and verified before removing this derivative.
        path.unlink()
    write(HERE/'cache-projection-summary.json',dict(verified_density_caches=len(ledger),
        reclaimed_bytes=sum(a['bytes_removed'] for a in ledger),scored_thresholds_per_cache=3,
        exact_primary_coordinates_and_words=True,raw_sources_retained=True,
        script_sha256=sha(__file__),reproduction='Frozen runners regenerate tensors; portable call tables reproduce scored endpoints without them'))
    print('Projected and removed',len(ledger),'verified density derivatives',flush=True)

if __name__=='__main__':main()
