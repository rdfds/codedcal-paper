"""Replay exact code assignments and preserve the prior evidence chain."""
from pathlib import Path
import argparse,hashlib,itertools,json
import numpy as np
from core import coded_components,collapsed,witness_vector
from calibration_feasibility import SOURCE,EXPECTED

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
FROZEN=['optical-information-resolution','reference-uncertainty-diagnosis',
        'shared-bias-control-qualification','reference-archive-completion',
        'finite-code-covariance']


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert not args.output.exists();prior={}
    for name in FROZEN:
        parent=BASE.parent/name;manifest=json.loads((parent/'ARTIFACTS.json').read_text());files=manifest['files']
        if isinstance(files,list):files={entry['path']:entry for entry in files}
        for relative,entry in files.items():
            path=parent/relative
            assert path.stat().st_size==entry['bytes'],path
            assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'],path
        prior[name]=len(files)
    assert sum(prior[name] for name in FROZEN[:4])==344
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
    with np.load(SOURCE) as d:
        z=d['calibration_z'].copy();n=d['calibration_counts'].copy()
        original=d['prior_groups'];old=d['groups'];source_ids=d['test_sources']
    frequency_source=BASE.parent/'optical-information-resolution/commercial1-calibration-design-v2.json'
    assert json.loads(frequency_source.read_text())['primary_q']==.008
    result=json.loads((BASE/'calibration-feasibility/result.json').read_text())
    unit=result['common_count_scale'];z/=unit;n/=unit
    newfile=BASE/'calibration-feasibility/calibration-group-design.npz'
    with np.load(newfile) as d:
        new=d['groups'];v=d['witness'];np.testing.assert_array_equal(d['source_ids'],source_ids)
        np.testing.assert_array_equal(d['baseline_groups'],original);np.testing.assert_array_equal(d['prior_cubic_groups'],old)
    assert len(set(source_ids))==len(z)==279
    np.testing.assert_allclose(v,witness_vector(.96,result['training_theta'],result['training_scale']),atol=1e-15)
    info={}
    for label,groups in [('original_groups',original),('prior_cubic_calibration_groups',old),('new_quartic_calibration_groups',new)]:
        assert sorted(np.unique(groups))==list(range(139))
        assert sorted(np.bincount(groups))==[2]*138+[3]
        group_var=0.;means=[];tables=0
        for group in range(139):
            ix=np.flatnonzero(groups==group)
            codes=np.array(list(itertools.product(range(4),repeat=len(ix))))
            mask=codes[:,None,:]==np.arange(4)[None,:,None]
            zz=np.einsum('gks,st->gkt',mask,z[ix]);nn=np.einsum('gks,st->gkt',mask,n[ix])
            actual=coded_components(zz,nn)
            np.testing.assert_allclose(actual.mean(0),collapsed(z[ix],n[ix]).sum(0),rtol=1e-10,atol=1e-10)
            group_var+=float((actual@v).var());means.append(actual.mean(0));tables+=len(codes)
        expected=result['configurations'][label];methods={x['method']:x for x in expected['methods']}
        np.testing.assert_allclose(group_var,methods['ordinary_H4']['conditional_code_witness_variance'],rtol=1e-8,atol=1e-10)
        np.testing.assert_allclose(np.sum(means,axis=0),collapsed(z,n).sum(0),rtol=1e-11,atol=1e-9)
        info[label]=dict(groups=139,sources=279,enumerated_assignments=tables,
                         literal_contraction_witness_code_variance=group_var)
    sim=json.loads((BASE/'simulation-results.json').read_text())
    assert len(sim['results'])==10 and sum(x['repetitions'] for x in sim['results'])==12400
    assert not sim['physical_data_used'] and not sim['original_physical_comparison_passed']
    for row in sim['results']:
        for method in row['methods'].values():
            assert 0<=method['e']['calls']<=row['repetitions']
    out=dict(prior_manifests_verified=prior,preceding_four_frozen_files_unchanged=344,
             calibration_input_sha256=EXPECTED,frequency_source=str(frequency_source.relative_to(ROOT)),
             frequency_source_sha256=hashlib.sha256(frequency_source.read_bytes()).hexdigest(),
             calibration_arrays_decoded=['calibration_z','calibration_counts','prior_groups','groups','test_sources'],
             validation_arrays_decoded=False,exact_assignment_replay=info,synthetic_repetitions=12400,
             original_physical_criterion_unchanged=True,qualified_independent_physical_control=False)
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
