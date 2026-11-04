"""Independent evidence recount; does not import covariance or study drivers."""
from pathlib import Path
import argparse,hashlib,json,itertools,sys
import numpy as np
from scipy.stats import t as student_t
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.uncertainty import fieller_set
from codedcal.geometric_uncertainty import geometric_ratio_set


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def contains(parts,t):
    return any((p['lower'] is None or p['lower']<=t) and (p['upper'] is None or t<=p['upper']) for p in parts)



@lru_cache(None)
def critical(n,confidence):
    return float(student_t.ppf((1+confidence)/2,n-1)**2)


def independent_ratio(mean,cov,n,confidence):
    x,y=map(float,mean);z=critical(n,confidence)
    a=y*y-z*cov[1,1];b=-2*x*y+2*z*cov[0,1];c=x*x-z*cov[0,0]
    if a==0:
        if b==0:parts=[(-np.inf,np.inf)] if c<=0 else [];kind='all_real' if c<=0 else 'empty'
        else:
            root=-c/b;parts=[(-np.inf,root)] if b>0 else [(root,np.inf)];kind='half_line'
    else:
        disc=b*b-4*a*c
        if abs(disc)<=8*np.finfo(float).eps*(b*b+abs(4*a*c)):disc=0.
        if disc<0:parts=[(-np.inf,np.inf)] if a<0 else [];kind='all_real' if a<0 else 'empty'
        elif disc==0:
            root=-b/(2*a);parts=[(root,root)] if a>0 else [(-np.inf,np.inf)];kind='bounded' if a>0 else 'all_real'
        else:
            lo,hi=sorted(((-b-np.sqrt(disc))/(2*a),(-b+np.sqrt(disc))/(2*a)))
            parts=[(lo,hi)] if a>0 else [(-np.inf,lo),(hi,np.inf)];kind='bounded' if a>0 else 'disjoint'
    usable=kind=='bounded' and parts[0][0]>0 and y>0
    ratio=x/y if y!=0 else None
    return kind,parts,usable,ratio,y>0


def independent_case(mean,cov,n,variant,target,case):
    if variant=='arithmetic':
        transform=np.array([[.5,0,.5,0],[0,.5,0,.5]])
        kind,parts,usable,ratio,positive=independent_ratio(transform@mean,transform@cov@transform.T,n,.95)
        available=positive and ratio is not None and ratio>0
    else:
        windows=[independent_ratio(mean[2*i:2*i+2],cov[2*i:2*i+2,2*i:2*i+2],n,.975) for i in range(2)]
        available=all(w[4] and w[3] is not None and w[3]>0 for w in windows)
        ratio=float(np.sqrt(windows[0][3])*np.sqrt(windows[1][3])) if available else None
        usable=all(w[2] for w in windows)
        positive=[]
        for w in windows:positive.append([(max(lo,0.),hi) for lo,hi in w[1] if hi>0])
        candidates=sorted((np.sqrt(a*c),np.sqrt(b)*np.sqrt(d)) for a,b in positive[0] for c,d in positive[1])
        parts=[]
        for lo,hi in candidates:
            if parts and lo<=parts[-1][1]:parts[-1]=(parts[-1][0],max(parts[-1][1],hi))
            else:parts.append((lo,hi))
        kind='empty' if not parts else 'disjoint' if len(parts)>1 else 'unbounded_positive' if np.isinf(parts[0][1]) else 'bounded'
    covers=any(lo<=target<=hi for lo,hi in parts)
    contains_one=any(lo<=1<=hi for lo,hi in parts)
    assert bool(available)==case['available']
    assert bool(usable)==case['usable']
    assert kind==case['kind'],(kind,case['kind'])
    assert bool(covers)==case['covers']
    assert bool(usable and covers)==case['usable_covers']
    assert bool(usable and covers and not contains_one)==case['usable_covers_excludes_one']
    if available:np.testing.assert_allclose(ratio,case['ratio'],rtol=1e-12,atol=1e-12)
    assert len(parts)==len(case['intervals'])
    for own,stored in zip(parts,case['intervals']):
        for number,key in zip(own,('lower','upper')):
            if np.isinf(number):assert stored[key] is None
            else:np.testing.assert_allclose(number,stored[key],rtol=1e-8,atol=1e-9)


def main():
    p=argparse.ArgumentParser();p.add_argument('--prefix',required=True);p.add_argument('--finite',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    pref,finite=Path(args.prefix),Path(args.finite)
    assert json.loads((pref/'state.json').read_text())['stage']=='complete'
    assert json.loads((finite/'state.json').read_text())['stage']=='complete'
    source={};worlds=0;prefix_calls=0;raw_api=0;finite_cases=0;energy_recounts=0
    worst_psd=0.;variance_z=[]
    for file in sorted(pref.glob('*.json')):
        if file.name=='state.json':continue
        d=json.loads(file.read_text());worlds+=1
        source.update(d['source_hashes']);prefix_calls+=len(d['replay'])*3
        assert all(r['known']['bit_equal'] and all(x['bit_equal'] for x in r['coded']) for r in d['replay'])
        if not d['theorem_applicable']:
            assert d['world']=='code_offset';continue
        gs=d['group_moments'];n=16
        pair=np.mean([r.get('pair',np.zeros((8,8))) for r in gs],axis=0)
        triple=np.mean([r.get('triple',np.zeros((8,8))) for r in gs],axis=0)
        mean=np.mean([r['mean'] for r in gs],axis=0)
        np.testing.assert_allclose(pair,d['pair_per_group'],rtol=1e-13,atol=1e-7)
        np.testing.assert_allclose(triple,d['triple_per_group'],rtol=1e-13,atol=1e-7)
        cov=(pair/3+triple/6)/n
        np.testing.assert_allclose(cov,d['verification']['predicted_covariance'],rtol=1e-13,atol=1e-7)
        eig=np.linalg.eigvalsh(cov);worst_psd=min(worst_psd,float(eig.min()/max(np.linalg.norm(cov),1e-300)))
        mc=np.array(d['verification']['mc_second_moment']);se=np.array(d['verification']['variance_mc_se'])
        np.testing.assert_allclose((np.diag(mc)-np.diag(cov))/se,d['verification']['variance_z'],rtol=1e-10,atol=1e-10)
        variance_z.extend(d['verification']['variance_z'])
        f=json.loads((finite/file.name).read_text())
        arrays_path=finite/file.with_suffix('.npz').name
        assert sha(arrays_path)==f['sufficient_statistics_sha256']
        np.testing.assert_allclose(f['mean'],mean,rtol=1e-12,atol=1e-8)
        np.testing.assert_allclose(f['pair'],pair,rtol=1e-13,atol=1e-7)
        np.testing.assert_allclose(f['triple'],triple,rtol=1e-13,atol=1e-7)
        floor=np.array(f['physical_covariance']);pop=np.array(f['mean'])
        with np.load(arrays_path) as arr:
            for r in f['comparisons']:
                if not r['supported']:continue
                k=r['K'];g=256//k;iq=0 if r['q']==.1 else 1;sl=slice(4*iq,4*iq+4)
                means=arr[f'K{k}_means'];sample_cov=arr[f'K{k}_mean_covariances']
                ca,ra,cb,rb=pop[sl]
                grad=np.array([1/(ca+cb),-1/(ra+rb),1/(ca+cb),-1/(ra+rb)]) if r['variant']=='arithmetic' else .5*np.array([1/ca,-1/ra,1/cb,-1/rb])
                prediction=grad@(floor+pair/(k-1)+triple/((k-1)*(k-2)))[sl,sl]@grad/g
                np.testing.assert_allclose(prediction,r['predicted_influence_variance'],rtol=1e-11)
                residual=(means[:,sl]-pop[sl])@grad
                np.testing.assert_allclose(np.mean(residual**2),r['observed_influence_second_moment'],rtol=1e-12)
                cases=r['cases'];assert len(cases)==200;finite_cases+=len(cases);energy_recounts+=1
                for m,sc,c in zip(means,sample_cov,cases):
                    independent_case(m[sl],sc[sl,sl],g,r['variant'],r['target'],c)
                for key in ('available','usable','covers','usable_covers','usable_covers_excludes_one'):
                    assert r[key]==sum(c[key] for c in cases)
                err=np.array([c['error'] for c in cases if c['error'] is not None])
                if len(err):np.testing.assert_allclose(np.sqrt(np.mean(err**2)),r['finite_available_rmse'],rtol=1e-12)
                raw=arr[f'K{k}_first_group_statistics'][:,sl]
                np.testing.assert_allclose(raw.mean(0),means[0,sl],rtol=1e-12,atol=1e-10)
                np.testing.assert_allclose(np.cov(raw,rowvar=False)/g,sample_cov[0,sl,sl],rtol=1e-10,atol=1e-7)
                ref=fieller_set(raw[:,[0,2]].mean(1),raw[:,[1,3]].mean(1)) if r['variant']=='arithmetic' else geometric_ratio_set(raw.reshape(g,2,2))
                c=cases[0]
                assert ref['kind']==c['kind'] and ref['usable']==c['usable']
                assert contains(ref['intervals'],r['target'])==c['covers']
                assert len(ref['intervals'])==len(c['intervals'])
                for x,y in zip(ref['intervals'],c['intervals']):
                    for endpoint in ('lower','upper'):
                        if x[endpoint] is None:assert y[endpoint] is None
                        else:np.testing.assert_allclose(x[endpoint],y[endpoint],rtol=1e-8,atol=1e-10)
                raw_api+=1
    for path,digest in source.items():assert sha(ROOT/path)==digest,path
    assert worlds==126 and len(list(finite.glob('*-*.json')))==108
    report=dict(status='pass',worlds=worlds,applicable_worlds=108,excluded_code_offset_worlds=18,
                historical_prefix_vectors=prefix_calls,unique_source_hashes=len(source),
                finite_cases_recounted=finite_cases,variance_recounts=energy_recounts,
                independent_raw_shipped_interval_checks=raw_api,
                independent_sufficient_statistic_interval_checks=finite_cases,smallest_scaled_covariance_eigenvalue=worst_psd,
                mc_variance_z_abs_max=float(np.max(np.abs(variance_z))),
                mc_variance_z_abs_over_3=int(np.sum(np.abs(variance_z)>3)),
                note='MC studentized variance discrepancies retained; exact covariance tests are separate. No physical model/priority/coverage validation follows from integrity checks.')
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
