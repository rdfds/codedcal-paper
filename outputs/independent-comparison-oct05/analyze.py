"""Read every arm, preserve failures, freeze the simple preflight, verify calls."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment

HERE=Path(__file__).resolve().parent
POLICIES=['none','common_round','common_iterative','population_template','barcode_model','oracle']
DECODERS=['correlation','regression']

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def write(path,value): Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def records(split):
    source=read(HERE/split/'source.json')
    found=[]
    for path in sorted((HERE/split).glob('*-d*-s*.json')):
        obj=read(path)
        assert obj['source_sha256']==source['sha256']
        assert obj['settings_sha256']==sha(HERE/'settings-freeze.json')
        assert obj['amendment_sha256']==sha(HERE/'PREVALIDATION_AMENDMENT.md')
        extra_path=HERE/(split+'-strong-common')/path.name
        if not extra_path.exists():raise ValueError('strong common-round comparison incomplete')
        extra=read(extra_path)
        assert extra['original_condition_sha256']==sha(path)
        assert extra['amendment_sha256']==sha(HERE/'STRONG_COMMON_AMENDMENT.md')
        assert extra['strong_script_sha256']==sha(HERE/'strong_common.py')
        obj['arms'].update(extra['arms'])
        obj['failures'].update({'common_iterative__'+k:v for k,v in extra['failures'].items()})
        obj['corrections_yx']['common_iterative']=extra['correction_yx']
        obj['_path']=str(path)
        found.append(obj)
    return found

def feature(obj,comparator='common_iterative'):
    a=np.array(obj['corrections_yx']['population_template'])
    b=np.array(obj['corrections_yx'][comparator])
    return float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))

def confusion(truth,pred):
    a,b=np.asarray(truth,bool),np.asarray(pred,bool)
    tp=int(np.sum(a&b)); tn=int(np.sum(~a&~b))
    fp=int(np.sum(~a&b)); fn=int(np.sum(a&~b))
    bal=(tp/(tp+fn)+tn/(tn+fp))/2 if (tp+fn) and (tn+fp) else None
    return dict(tp=tp,tn=tn,fp=fp,fn=fn,balanced_accuracy=bal)

def freeze():
    data=records('development')
    assert len(data)==13, 'development incomplete'
    settings=read(HERE/'settings-freeze.json')
    checks={}
    for decoder in DECODERS:
        baseline=settings['thresholds'][decoder]
        zero=next(d for d in data if d['condition'].startswith('zero'))
        zero_f1=zero['arms']['none__'+decoder]['metrics']['5']['identity_f1']
        oracle_drops=[zero_f1-d['arms']['oracle__'+decoder]['metrics']['5']['identity_f1'] for d in data]
        checks[decoder]=dict(unperturbed_f1=zero_f1,max_oracle_drop=max(oracle_drops),
                             oracle_qualification=max(oracle_drops)<=.05)
    nonzero=[d for d in data if not d['condition'].startswith('zero')]
    x=np.array([feature(d) for d in nonzero])
    frozen=dict(extension_sha256=sha(HERE/'DECISION_EXTENSION.md'),
                analysis_sha256=sha(__file__),settings_sha256=sha(HERE/'settings-freeze.json'),
                development_sha256={Path(d['_path']).name:sha(d['_path']) for d in data},
                strong_amendment_sha256=sha(HERE/'STRONG_COMMON_AMENDMENT.md'),
                oracle_checks=checks,rules={},phase_only_rules={})
    for comparator,output in [('common_iterative','rules'),('common_round','phase_only_rules')]:
        x=np.array([feature(d,comparator) for d in nonzero])
        for decoder in DECODERS:
            gains=np.array([d['arms']['population_template__'+decoder]['metrics']['5']['identity_f1']-
                            d['arms'][comparator+'__'+decoder]['metrics']['5']['identity_f1'] for d in nonzero])
            truth=gains>=.05
            curve=[]
            for threshold in [0,.25,.5,.75,1,1.5,2,3,4,8]:
                curve.append(dict(threshold=threshold,**confusion(truth,x>=threshold)))
            selected=max(curve,key=lambda v:(v['balanced_accuracy'] if v['balanced_accuracy'] is not None else -1,
                                            -v['fp'],v['threshold']))
            qualified=(int(truth.sum())>=3 and int((~truth).sum())>=3 and
                       selected['balanced_accuracy'] is not None and selected['balanced_accuracy']>=.8)
            frozen[output][decoder]=dict(qualified=qualified,selected=selected,curve=curve,
                cases=[dict(condition=d['condition'],feature=float(xx),gain=float(g),positive=bool(t))
                       for d,xx,g,t in zip(nonzero,x,gains,truth)])
    frozen['validation_comparison_qualified']=all(d['oracle_qualification'] for d in checks.values())
    # A bad predictor does not block the separately frozen scientific comparison.
    write(HERE/'preflight-freeze.json',frozen)
    print(json.dumps({k:frozen[k] for k in ['oracle_checks','validation_comparison_qualified']},indent=2))
    print({k:v['selected']|{'qualified':v['qualified']} for k,v in frozen['rules'].items()})

def simple_global_match(a,b,radius):
    """Independent all-points assignment; checks cardinality, not tie label choice."""
    if not len(a) or not len(b):return 0
    nr,nc=len(a),len(b)
    if nr+nc>5000:raise ValueError('global verification memory ceiling exceeded')
    dist=np.linalg.norm(a[:,None]-b[None],axis=2)
    reward=(min(nr,nc)+1)*radius
    cost=np.zeros((nr+nc,nr+nc))
    cost[:nr,:nc]=np.where(dist<=radius,dist-reward,reward)
    rr,cc=linear_sum_assignment(cost)
    keep=(rr<nr)&(cc<nc)
    return int(np.sum(dist[rr[keep],cc[keep]]<=radius))

def report():
    # Avoid importing the TensorFlow scientific runner in this read-only analysis.
    import h5py
    frozen=read(HERE/'preflight-freeze.json')
    summaries={};rows=[];manifest=[];all_failures=[];verified=0
    for split in ['development','validation']:
        data=records(split)
        assert len(data)==13, split+' incomplete'
        assert all(d['run_sha256']==sha(HERE/'run.py') for d in data), 'executed source changed'
        source=read(HERE/split/'source.json')
        with h5py.File(HERE.parent.parent.parent/'bardensr-upstream.nosync/bardensr/benchmarks'/(source['name']+'.hdf5')) as f:
            statuses=f['rolonies/status'][:].astype('U')
            xy=np.c_[f['rolonies/m2'][:],f['rolonies/m1'][:]].astype(float)
            keep=(statuses=='good')&np.all((xy>=24)&(xy<176),axis=1)
            reference=xy[keep]
        selected_rows={decoder:[] for decoder in DECODERS}
        for obj in data:
            condition=obj['condition'];kind=condition.split('-')[0]
            all_failures.extend([dict(split=split,condition=condition,arm=k,traceback=v)
                                 for k,v in obj['failures'].items()])
            for policy in POLICIES:
                for decoder in DECODERS:
                    key=policy+'__'+decoder
                    if key not in obj['arms']:
                        rows.append(dict(split=split,condition=condition,kind=kind,policy=policy,decoder=decoder,status='failed'))
                        continue
                    arm=obj['arms'][key]
                    metrics=arm['metrics']['5']
                    arm_folder=HERE/(split+'-strong-common') if policy=='common_iterative' else HERE/split
                    path=arm_folder/(condition+'__'+key+'.npz')
                    with np.load(path) as raw:
                        xy_pred=raw['xy'];ids=raw['words']
                        mask=np.all((xy_pred>=24)&(xy_pred<176),axis=1)
                        pred=xy_pred[mask]
                        assert len(pred)==metrics['calls']
                        card=simple_global_match(reference,pred,5)
                        assert card==metrics['spatial_matches'], 'independent assignment cardinality mismatch'
                        verified+=1
                        stream=io.StringIO();out=csv.writer(stream,lineterminator='\n');out.writerow(['x','y','word'])
                        out.writerows(zip(xy_pred[:,0],xy_pred[:,1],ids))
                        csvpath=arm_folder/(condition+'__'+key+'.calls.csv.gz')
                        csvpath.write_bytes(gzip.compress(stream.getvalue().encode(),mtime=0))
                    manifest.append(dict(path=str(path.relative_to(HERE)),sha256=sha(path),
                                         calls_path=str(csvpath.relative_to(HERE)),calls_sha256=sha(csvpath)))
                    runtime=arm['runtime']
                    rows.append(dict(split=split,condition=condition,kind=kind,
                        imposed_rms=obj['imposed_rms'],feature=feature(obj),policy=policy,
                        decoder=decoder,status='completed',**{k:v for k,v in metrics.items() if k!='matching'},
                        residual_imposed_rms=arm['residual_imposed_rms'],
                        total_core_s=sum(runtime[k] for k in ['preprocessing_s','decoder_s','calibration_s','warp_s']),
                        ambiguity=metrics['matching']['ambiguous_reference']))
            if kind!='zero':
                for decoder in DECODERS:
                    rule=frozen['rules'][decoder]
                    if not rule['qualified']: continue
                    recommend='population_template' if feature(obj)>=rule['selected']['threshold'] else 'common_iterative'
                    a=obj['arms'];f=lambda p:a[p+'__'+decoder]['metrics']['5']['identity_f1']
                    gain=f('population_template')-f('common_iterative')
                    selected_rows[decoder].append(dict(condition=condition,feature=feature(obj),
                        recommend=recommend,truth=gain>=.05,predict=recommend=='population_template',
                        selected_f1=f(recommend),common_f1=f('common_iterative'),phase_only_f1=f('common_round'),template_f1=f('population_template'),
                        model_f1=f('barcode_model'),harm_vs_common=f(recommend)-f('common_iterative')<=-.05))
        decisions={}
        for decoder, cases in selected_rows.items():
            decisions[decoder]=dict(qualified=frozen['rules'][decoder]['qualified'],cases=cases,
                confusion=confusion([c['truth'] for c in cases],[c['predict'] for c in cases]) if cases else None,
                mean_f1={key:float(np.mean([c[key] for c in cases])) for key in
                         ['selected_f1','common_f1','template_f1','model_f1']} if cases else None,
                harms_vs_common=sum(c['harm_vs_common'] for c in cases))
        summaries[split]=dict(source=source,conditions=len(data),decision=decisions)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (HERE/'all-arms.csv').open('w') as f:
        out=csv.DictWriter(f,fieldnames=fields);out.writeheader();out.writerows(rows)
    write(HERE/'density-and-call-manifest.json',manifest)
    write(HERE/'comparison-summary.json',dict(preflight_freeze_sha256=sha(HERE/'preflight-freeze.json'),
        studies=summaries,expected_arms=312,successful_arms=verified,failures=all_failures,
        verification='Primary call counts and spatial cardinality independently checked by global Hungarian assignment; label tie ambiguity remains disclosed.'))
    print('Verified arms',verified,'failures',len(all_failures))
    for split,value in summaries.items():
        print(split,{k:v['confusion'] for k,v in value['decision'].items()})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','report'])
    args=parser.parse_args()
    (freeze if args.stage=='freeze' else report)()
