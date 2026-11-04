"""Independent frozen-formula recount; scientific artifacts are read only."""
import ast
import gzip
import json
import numpy as np
from scipy.stats import t
from methods import contains
from run import HERE, SOURCE, sha, write, assert_frozen


def frozen_formula():
    node=next(n for n in ast.parse((SOURCE/'validation.py').read_text()).body
              if isinstance(n,ast.FunctionDef) and n.name=='fieller')
    namespace=dict(np=np,student_t=t)
    exec(compile(ast.Module(body=[node],type_ignores=[]),'frozen-Fieller','exec'),namespace)
    return namespace['fieller']


def main():
    assert_frozen();formula=frozen_formula();cases=0;windows=0;arithmetic=0;worlds=0
    for name,digest in json.loads((HERE/'source-hashes.json').read_text()).items():
        assert sha(SOURCE/name)==digest,name
    for file in sorted((HERE/'worlds').glob('*.jsonl.gz')):
        label=file.name.removesuffix('.jsonl.gz');folder=SOURCE/'worlds'/label
        marker=json.loads(file.with_suffix('.json').read_text());assert sha(file)==marker['output_sha256']
        loaded_rep=None;arrays=None;worlds+=1
        with gzip.open(file,'rt') as stream:
            for line in stream:
                row=json.loads(line)
                if loaded_rep!=row['rep']:
                    if arrays is not None:arrays.close()
                    arrays=np.load(folder/f"field-{row['rep']:03d}.npz");loaded_rep=row['rep']
                    original=json.loads((folder/f"field-{row['rep']:03d}.json").read_text())
                    old={(r['scope'],r['method'],r['assignment'],r['q']):r for r in original['cases']}
                data=(arrays['coded'][row['assignment']] if row['method']=='coded' else arrays['known_source'])
                a=data[:row['group_count'],(0 if row['q']==.1 else 1)]
                prior=old[row['scope'],row['method'],row['assignment'],row['q']]
                assert row['D']==prior['D']
                if row['variant']=='arithmetic':
                    assert row['kind']==prior['kind'] and row['usable']==prior['usable']
                    assert row['covers']==prior['covers'];arithmetic+=1
                else:
                    ratios=[]
                    for w,observed in enumerate(row['windows']):
                        f=formula(a[:,w,0],a[:,w,1],1.,alpha=.025)
                        assert observed['kind']==f['kind'] and observed['usable']==f['usable']
                        if f['kind']=='bounded':
                            np.testing.assert_allclose([observed['intervals'][0]['lower'],observed['intervals'][0]['upper']],
                                                       [f['lower'],f['upper']],rtol=2e-8,atol=1e-10)
                        mean=a[:,w].mean(0);ratios.append(mean[0]/mean[1] if mean[1]>0 and mean[0]>0 else None)
                        windows+=1
                    expected=float(np.sqrt(ratios[0])*np.sqrt(ratios[1])) if all(r is not None for r in ratios) else None
                    if expected is None:assert row['ratio'] is None
                    else:assert np.isclose(row['ratio'],expected,rtol=1e-10,atol=1e-12)
                    assert row['covers']==contains(row['intervals'],row['physical_target'])
                cases+=1
        if arrays is not None:arrays.close()
    summary=json.loads((HERE/'summary.json').read_text())['summaries']
    old=json.loads((SOURCE/'summary.json').read_text())['records']
    index={(r['tag'],r['model'],r['world'],r['scope'],r['method'],r['q']):r for r in old}
    arithmetic_summaries=0
    for row in summary:
        if row['variant']!='arithmetic':continue
        prior=index[tuple(row[k] for k in ('tag','model','world','scope','method','q'))]
        assert row['qualified']==prior['qualified']
        for key in ('ratio_bias','usable_fraction','usable_joint_coverage','median_usable_interval_width'):
            if row[key] is None:assert prior[key] is None
            else:assert np.isclose(row[key],prior[key],rtol=2e-8,atol=1e-10),(key,row,prior)
        np.testing.assert_allclose(row['usable_joint_coverage_interval'],prior['usable_joint_coverage_interval'],atol=1e-12)
        arithmetic_summaries+=1
    assert worlds==126 and cases==84672 and arithmetic==42336 and windows==84672
    assert arithmetic_summaries==1008 and len(summary)==2016
    write(HERE/'verification.json',dict(all_pass=True,worlds=worlds,cases=cases,
          arithmetic_original_cases=arithmetic,geometric_window_formula_checks=windows,
          arithmetic_original_summaries=arithmetic_summaries,summary_rows=len(summary),
          independent_physical_validation=False,freeze_sha256=sha(HERE/'freeze.json')))
    print('Verified',cases,'cases',flush=True)


if __name__=='__main__':main()
