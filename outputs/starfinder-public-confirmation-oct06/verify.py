"""Independent retained-output recount; no TensorFlow or original metric function."""
import collections,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent.parent
HERE=Path(__file__).resolve().parent

def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while b:=f.read(1024**2):h.update(b)
    return h.hexdigest()
def verify():
    records=[];refs=read(HERE/'probe/references.json')
    for folder,expected in [('probe',54),('source-comparators',27)]:
        paths=list((HERE/folder).glob('*__*.json'));assert len(paths)==expected
        for p in paths:
            r=read(p);assert r['status']=='completed'
            pred=r['prediction'];words=pred['words'];n=sum(w<998 for w in words);b=sum(w>=998 for w in words)
            assert n==pred['coding'] and b==pred['negative'] and len(words)==len(pred['xyz'])
            assert len(words)+pred['abstained']==pred['candidate_objects']
            native=refs[r['decoder']];ca=collections.Counter(w for w in native['words'] if w<998);cb=collections.Counter(w for w in words if w<998)
            tv=sum(abs(ca[w]/sum(ca.values())-cb[w]/sum(cb.values())) for w in ca.keys()|cb.keys())/2 if cb else None
            for m in r['metrics'].values():
                assert m['reference_calls']==len(native['words']) and m['predicted_calls']==len(words)
                if tv is not None:assert abs(m['composition_tv']-tv)<1e-12
                assert abs(m['identity_f1']-2*m.get('correct_labels',0)/(len(native['words'])+len(words)))<1e-12
                assert 0<=m['identity_f1']<=1
            if r['metrics']['3']['identity_f1']==1:assert pred['xyz']==native['xyz'] and words==native['words']
            records.append(dict(file=str(p.relative_to(ROOT)),sha256=sha(p),coding=n,negative=b,identity_f1=r['metrics']['3']['identity_f1']))
    # All original frozen libraries/core and source manifests remain byte-identical.
    libraries=[]
    for frozen in ['freeze.json','local-trace/freeze.json']:
        for rel,digest in read(HERE/frozen)['files'].items():assert sha(ROOT/rel)==digest;libraries.append(rel)
    core=ROOT/'outputs/barseq3-native-repair-oct05/global_repair.py';assert sha(core)=='82407386fccc55da8b745fd5a31c83ed27d41e4c4fcb2eb1b3267d5c14df4fdc'
    init=list((HERE/'initialization').glob('starfinder_*.json'))+list((HERE/'initialization').glob('barseq3_*.json'));assert len(init)==30 and read(HERE/'initialization/state.json')['completed']==30
    for p in init:assert read(p)['status']=='completed';records.append(dict(file=str(p.relative_to(ROOT)),sha256=sha(p)))
    result=dict(passed=True,primary_records=81,initialization_records=30,verified_records=len(records),frozen_files_verified=len(set(libraries)),unchanged_core_sha256=sha(core),counts_composition_and_identity_scalars_verified=True,perfect_records_exact_catalog_equality=True,independent_RNA_truth=False,independent_preparation_workflow_confirmation=False,records=records)
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',len(records),'retained scientific rows')
if __name__=='__main__':verify()
