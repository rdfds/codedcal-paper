"""One finite streaming worker; never modifies Source Data 69/70."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import time
import platform
import scipy
from collections import Counter
import numpy as np
from methods import arithmetic_target, geometric_target, contains, estimate

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'temporal-validation-oct09'
TAGS = ('npc','edf10rep1','edf10rep2','fig4untreated','rtx_rep1','rtx_rep2','rtx_rep3','sdc_corrected','sdc_raw')
MODELS = ('gaussian','residual')
PRIMARY = ('zero','homogeneous','group_common','heterogeneous')
WORLDS = PRIMARY+('noise_shift','source_noise_shift','code_offset')
Q = (.1,.15)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024**2), b''):h.update(block)
    return h.hexdigest()


def write(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n'); temporary.replace(path)


def freeze():
    assert not (HERE/'freeze.json').exists(), 'Do not replace an existing freeze'
    files = [HERE/'PROTOCOL.md', HERE/'methods.py', HERE/'run.py', HERE/'controls.py', HERE/'controls.json',
             HERE/'verify.py', SOURCE/'artifact-manifest.json', SOURCE/'summary.json', SOURCE/'freeze.json',
             HERE.parents[1]/'tool/src/codedcal/uncertainty.py']
    write(HERE/'freeze.json', dict(protocol='retained-observation comparison',
          files={str(p.relative_to(HERE.parents[1])):sha(p) for p in files},
          environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
          original_outcomes_already_exposed=True, no_new_simulations=True))


def assert_frozen():
    frozen = json.loads((HERE/'freeze.json').read_text())
    for name, digest in frozen['files'].items():
        assert sha(HERE.parents[1]/name)==digest, name


def case(groups, variant, target, observable, original, arithmetic_physical):
    result = estimate(groups, variant)
    available = result['ratio'] is not None and np.isfinite(result['ratio'])
    result.update(physical_target=target, observable_target=observable,
                  arithmetic_physical_target=arithmetic_physical,
                  available=bool(available), covers=contains(result['intervals'], target),
                  observable_covers=contains(result['intervals'], observable),
                  lower=(result['intervals'][0]['lower'] if result['usable'] else None),
                  upper=(result['intervals'][0]['upper'] if result['usable'] else None),
                  D=original['D'])
    if variant=='arithmetic':
        assert result['kind']==original['kind'], (result,original)
        assert result['usable']==original['usable']
        assert result['covers']==original['covers']
        if available:
            assert np.isclose(result['ratio'],original['ratio'],rtol=1e-10,atol=1e-12)
    return result


def execute():
    assert_frozen()
    lock = HERE/'worker.lock'
    descriptor = os.open(lock, os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.write(descriptor,str(os.getpid()).encode());os.close(descriptor)
    start = time.monotonic(); total_fields=0; total_cases=0; sources={}
    manifest=json.loads((SOURCE/'artifact-manifest.json').read_text())['artifacts']
    try:
        for tag in TAGS:
            for model in MODELS:
                for world in WORLDS:
                    label=f'{tag}-{model}-{world}'; folder=SOURCE/'worlds'/label
                    state_path=folder/'state.json'; state=json.loads(state_path.read_text())
                    assert state['complete'] and state['completed']==state['expected']
                    out=HERE/'worlds'/f'{label}.jsonl.gz'; marker=out.with_suffix('.json')
                    if marker.exists():
                        old=json.loads(marker.read_text());assert sha(out)==old['output_sha256']
                        assert old['freeze_sha256']==sha(HERE/'freeze.json')
                        total_fields+=old['fields'];total_cases+=old['cases'];sources.update(old['sources']);continue
                    assert not out.exists(), f'Unmarked outcome {out}; preserve/diagnose before resuming'
                    out.parent.mkdir(exist_ok=True)
                    truth_path=folder/'truth.npz'
                    assert sha(state_path)==manifest[str(state_path.relative_to(SOURCE))]
                    assert sha(truth_path)==manifest[str(truth_path.relative_to(SOURCE))]
                    with np.load(truth_path) as data:pop=data['population']
                    local_sources={str(p.relative_to(SOURCE)):sha(p) for p in (state_path,truth_path)}
                    count=0
                    with gzip.open(str(out)+'.tmp','wt') as stream:
                        for entry in state['replicates']:
                            path=folder/entry['filename']; rp=folder/f"field-{entry['rep']:03d}.json"
                            assert sha(path)==entry['sha256'] and sha(rp)==entry['record_sha256']
                            local_sources[str(path.relative_to(SOURCE))]=entry['sha256']
                            local_sources[str(rp.relative_to(SOURCE))]=entry['record_sha256']
                            original=json.loads(rp.read_text())
                            with np.load(path) as arrays:
                                for old in original['cases']:
                                    data=(arrays['coded'][old['assignment']] if old['method']=='coded' else arrays['known_source'])
                                    j=Q.index(old['q']); groups=data[:old['group_count'],j]
                                    for variant, functional in [('arithmetic',arithmetic_target),('geometric',geometric_target)]:
                                        result=case(groups,variant,functional(pop[0,j]),functional(pop[1,j]),old,arithmetic_target(pop[0,j]))
                                        result.update(tag=tag, model=model, world=world, variant=variant,
                                                      rep=old['rep'], scope=old['scope'], group_count=old['group_count'],
                                                      method=old['method'],assignment=old['assignment'],q=old['q'])
                                        stream.write(json.dumps(result,allow_nan=False)+'\n');count+=1
                    Path(str(out)+'.tmp').replace(out)
                    write(marker,dict(freeze_sha256=sha(HERE/'freeze.json'),output_sha256=sha(out),
                                      fields=state['completed'],cases=count,sources=local_sources))
                    total_fields+=state['completed'];total_cases+=count;sources.update(local_sources)
                    write(HERE/'execution-state.json',dict(stage='running',last_world=label,fields=total_fields,
                          cases=total_cases,pid=os.getpid(),elapsed_s=time.monotonic()-start))
                    print(label,count,flush=True)
        write(HERE/'source-hashes.json',sources)
        write(HERE/'execution-state.json',dict(stage='complete',fields=total_fields,cases=total_cases,
              worlds=len(TAGS)*len(MODELS)*len(WORLDS),elapsed_s=time.monotonic()-start,
              freeze_sha256=sha(HERE/'freeze.json')))
    finally:
        lock.unlink(missing_ok=True)


def summarize():
    assert_frozen()
    assert json.loads((HERE/'execution-state.json').read_text())['stage']=='complete'
    summaries=[]
    for file in sorted((HERE/'worlds').glob('*.jsonl.gz')):
        buckets={}
        with gzip.open(file,'rt') as stream:
            for line in stream:
                r=json.loads(line);key=(r['scope'],r['method'],r['q'],r['variant'])
                buckets.setdefault(key,[]).append(r)
        for key, rows in buckets.items():
            rows.sort(key=lambda r:(r['rep'],r['assignment']))
            available=np.array([r['available'] for r in rows]);usable=np.array([r['usable'] for r in rows])
            covers=np.array([r['covers'] for r in rows]);joint=usable&covers
            ratio=np.array([r['ratio'] if r['ratio'] is not None else np.nan for r in rows])
            target=rows[0]['physical_target'];observable=rows[0]['observable_target']
            assignments=2 if key[1]=='coded' else 1;fields=len(rows)//assignments
            matrix=joint.reshape(fields,assignments)
            first=rows[0]
            label=f"{first['tag']}-{first['model']}-{first['world']}-{key[0]}-{key[1]}-{key[2]}"
            # Original arithmetic bootstrap and paired variant use identical ids.
            seed=int.from_bytes(hashlib.sha256((label+'|bootstrap').encode()).digest()[:8],'little')
            rng=np.random.default_rng(seed);ids=rng.integers(fields,size=(2000,fields))
            ci=np.quantile(matrix[ids].mean(axis=(1,2)),[.025,.975]).tolist()
            bias=float(np.mean(ratio[available])-target) if available.any() else None
            width=np.array([r['upper']-r['lower'] if r['usable'] else np.nan for r in rows])
            resolves=np.array([r['usable'] and r['upper']<1 for r in rows])
            quiet=np.array([r['D']['quiet'] for r in rows]);valid_d=np.array([r['D']['valid'] for r in rows])
            bad=available&(abs(ratio-target)>.01)
            gates=dict(ratio_bias=bias is not None and abs(bias)<=.01,
                       usable_fraction=usable.mean()>=.90,usable_coverage=joint.mean()>=.90,
                       nominal_coverage_in_bootstrap=ci[0]<=.95<=ci[1])
            summary={k:first[k] for k in ('tag','model','world','scope','method','q','variant','group_count','physical_target','observable_target','arithmetic_physical_target')}
            summary.update(records=len(rows),physical_fields=fields,assignments_per_field=assignments,
                primary=first['world'] in PRIMARY,available_fraction=float(available.mean()),
                mean_ratio=float(np.mean(ratio[available])) if available.any() else None,
                ratio_bias=bias,observable_ratio_bias=float(np.mean(ratio[available])-observable) if available.any() else None,
                original_arithmetic_target_bias=float(np.mean(ratio[available])-first['arithmetic_physical_target']) if available.any() else None,
                target_change=target-first['arithmetic_physical_target'],
                ratio_rmse=float(np.sqrt(np.mean((ratio[available]-target)**2))) if available.any() else None,
                ratio_quantiles=np.quantile(ratio[available],[0,.01,.5,.99,1]).tolist() if available.any() else [],
                usable_fraction=float(usable.mean()),unconditional_coverage=float(covers.mean()),
                observable_coverage=float(np.mean([r['observable_covers'] for r in rows])),
                usable_joint_coverage=float(joint.mean()),usable_joint_coverage_interval=ci,
                conditional_usable_coverage=float(covers[usable].mean()) if usable.any() else None,
                median_usable_interval_width=float(np.median(width[usable])) if usable.any() else None,
                usable_covers_and_excludes_zero_fraction=float((resolves&covers).mean()),
                D_invalid=int((~valid_d).sum()),bad_physical_errors=int(bad.sum()),bad_errors_quiet_D=int((bad&quiet).sum()),
                set_types=dict(Counter(r['kind'] for r in rows)),
                coverage_bootstrap_nominal_failure_direction='none' if gates['nominal_coverage_in_bootstrap'] else 'above_nominal' if ci[0]>.95 else 'below_nominal',
                qualified=all(gates.values()),failed_gates=[k for k,v in gates.items() if not v])
            summaries.append(summary)
    write(HERE/'summary.json',dict(complete=True,conditional_only=True,physical_validation=False,
          summaries=summaries,freeze_sha256=sha(HERE/'freeze.json')))
    fields=[k for k,v in summaries[0].items() if not isinstance(v,(dict,list))]
    with (HERE/'summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        writer.writerows({k:r[k] for k in fields} for r in summaries)
    print('Summary rows:',len(summaries),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','execute','summarize'])
    args=parser.parse_args();globals()[args.stage]()
