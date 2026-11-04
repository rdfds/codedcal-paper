"""Recount completed primary endpoints using retained calls, no dense TF cache."""
import ast
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import h5py
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.evaluation import spatial_match
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def inside(xy):return np.all((xy>=24)&(xy<176),axis=1)
namespace=dict(np=np,inside=inside,spatial_match=spatial_match)
node=next(n for n in ast.parse((HERE/'run.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='score')
exec(compile(ast.Module(body=[node],type_ignores=[]),'frozen score','exec'),namespace)
score=namespace['score']

def main():
    references={};verified=0
    for split in ['development','validation']:
        source=read(HERE/split/'source.json')
        path=ROOT.parent/'bardensr-upstream.nosync/bardensr/benchmarks'/(source['name']+'.hdf5')
        assert sha(path)==source['sha256']
        with h5py.File(path) as f:
            xy=np.c_[f['rolonies/m2'][:],f['rolonies/m1'][:]];raw_words=f['rolonies/j'][:]
            assert np.array_equal(raw_words,raw_words.astype(int))
            words=raw_words.astype(int)
            mask=(f['rolonies/status'][:].astype('U')=='good')&inside(xy)
            references[split]=(xy[mask],words[mask],f['codebook'].shape[-1])
    for item in read(HERE/'density-and-call-manifest.json'):
        path=HERE/item['calls_path'];assert sha(path)==item['calls_sha256']
        with gzip.open(path,'rt') as f:rows=list(csv.DictReader(f))
        calls=np.array([[float(r['x']),float(r['y'])] for r in rows]).reshape(-1,2)
        ids=np.array([int(r['word']) for r in rows],dtype=int)
        condition,policy,decoder=path.name.split('__');decoder=decoder.split('.')[0]
        split=path.parent.name.split('-strong')[0]
        arm=read(path.parent/(condition+'.json'))['arms'][policy+'__'+decoder]
        for radius in [1.5,3,5]:
            actual=score(*references[split][:2],calls,ids,references[split][2],radius)
            expected=arm['metrics'][str(radius)]
            for key in ['calls','correct_labels','spatial_matches','identity_f1','spatial_f1','label_constrained_f1','wrong_label_fraction','composition_tv']:
                a,b=actual[key],expected[key]
                assert (a is None and b is None) or (a is not None and b is not None and abs(a-b)<1e-12),(item['calls_path'],radius,key,a,b)
        verified+=1
    assert verified==312
    if (HERE/'cache-projection-ledger.json').exists():
        for entry in read(HERE/'cache-projection-ledger.json'):
            for value in entry['thresholds'].values():assert sha(HERE/value['path'])==value['sha256']
    result=dict(completed_primary_arms=verified,radii=[1.5,3,5],exact_metric_recounts=True,
        author_raw_sources_verified=True,script_sha256=sha(__file__),frozen_score_source_sha256=sha(HERE/'run.py'))
    (HERE/'portable-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)

if __name__=='__main__':main()
