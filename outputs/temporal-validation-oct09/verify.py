"""Independent numerical recount of retained group observations and truth."""
import json
import hashlib
import subprocess
import sys
import scipy
from pathlib import Path
import numpy as np
from scipy.stats import t
import validation as V
from common import save, sha


def equal(got,expected,tolerance=1e-10):
    if expected is None:return not np.isfinite(got)
    return bool(np.isclose(got,expected,rtol=tolerance,atol=tolerance))


def main():
    summary=json.loads((V.P/'summary.json').read_text());assert summary['complete']
    manifest=json.loads((V.P/'artifact-manifest.json').read_text())
    frozen=json.loads((V.P/'freeze.json').read_text())
    assert frozen['legacy_components_sha256']==sha(V.ROOT/'outputs/stationarity-diagnostic-sep15/split_components.py')
    assert summary['code_sha256']==sha(V.P/'summarize.py')
    dependencies={}
    for rel in ['outputs/stationarity-diagnostic-sep15/split_components.py',
                'outputs/count-balanced-calibration-sep14/balanced_estimator.py',
                'outputs/noise-equalization-sep14/test_equalization.py',
                'outputs/reference-reconciliation-oct07/common.py',
                'outputs/reference-reconciliation-oct07/simulate_reference.py',
                'outputs/reference-reconciliation-oct07/reference_characterize.py',
                'outputs/oracle-crossval-sep15/oracle_pipeline.py']:
        frozen_bytes=subprocess.check_output(['git','show',f'44b98a8:{rel}'],cwd=V.ROOT)
        digest=hashlib.sha256(frozen_bytes).hexdigest()
        assert sha(V.ROOT/rel)==digest,rel
        dependencies[rel]=digest
    checks={};bad=[];cases=0;fields=0;catalogs=0;truths=0;information=[]
    for name,digest in manifest['artifacts'].items():
        assert sha(V.P/name)==digest,name
    for tag in V.TAGS:
        template,meta,cat=V.templates(tag)
        saved=np.load(V.P/f'catalog-{tag}.npz')
        for k,value in cat.items():assert np.array_equal(saved[k],value),(tag,k)
        original=np.load(V.REF/f'empirical-{tag}.npz')
        assert cat['ec'].sum()==len(original['events'])
        assert cat['lc'].sum()==original['events'][:,5].sum()
        assert abs(cat['factor'].mean()-1)<1e-12
        assert len(cat['pore'])==meta['inferred_sources'];catalogs+=1
        # Derive stationary pair kernel independently from covariances / empirical residual pool.
        kernels={}
        for model in V.MODELS:
            rows=[]
            for g in template:
                if g is None:rows.append(np.empty((0,2)));continue
                covproj=np.einsum('ai,sij,aj->sa',V.U,g['cov'],V.U)
                K=np.exp(-V.Q[None,:,None]**2*covproj[:,None,:])
                if model=='residual':
                    phi=np.zeros((len(g['n']),2,4),complex)
                    residual=g['res']*np.sqrt(g['n'][g['source']]/np.maximum(g['df'][g['source']],1))[:,None]
                    ph=np.exp(1j*(residual@V.U.T)[:,None,:]*V.Q[None,:,None])
                    np.add.at(phi,g['source'],ph)
                    phi/=g['n'][:,None,None]
                    K=np.where((g['df']>=4)[:,None,None],abs(phi)**2,K)
                rows.append(K.mean(-1))
            kernels[model]=rows
            weights=[];pore=[];factor=[]
            for g,K in zip(template,rows):
                if g is None:continue
                triple=sum(g['lc'][:,i]*g['lc'][:,j]*g['lc'][:,k]/2 for win in V.TRIPLES for i,j,k in win)
                weights.append(K*triple[:,None]);pore.append(np.full(len(g['n']),g['pore']));factor.append(g['factor'])
            weights=np.concatenate(weights);pore=np.concatenate(pore);factor=np.concatenate(factor)
            for iq,q in enumerate(V.Q):
                w=weights[:,iq];groupweight=np.bincount(pore,weights=w,minlength=meta['groups'])
                heterogeneous=(w*np.exp(-3*factor*q*q)).sum()/w.sum()
                information.append(dict(tag=tag,model=model,q=float(q),inferred_sources=meta['inferred_sources'],
                    triple_supported_sources=meta['source_triple_supported'],fallback_sources=meta['fallback_sources'],
                    effective_weighted_sources=float(w.sum()**2/(w*w).sum()),
                    effective_weighted_groups=float(groupweight.sum()**2/(groupweight*groupweight).sum()),
                    largest_source_weight_fraction=float(w.max()/w.sum()),
                    largest_group_weight_fraction=float(groupweight.max()/groupweight.sum()),
                    heterogeneous_source_mean_T=3.,heterogeneous_weighted_mean_T=float((w*3*factor).sum()/w.sum()),
                    heterogeneous_finite_q_apparent_T=float(-np.log(heterogeneous)/(q*q)),
                    fallback_kernel_weight_fraction=float(w[cat['df']<4].sum()/w.sum())))
            for world in V.WORLDS:
                stored=np.load(V.P/'worlds'/f'{tag}-{model}-{world}'/'truth.npz')
                independent=[]
                for g,K in zip(template,rows):
                    value=np.zeros((2,2,2,2))
                    if g is not None:
                        T=np.full(len(g['n']),0. if world=='zero' else 3.)
                        if world=='heterogeneous':T*=g['factor']
                        att=np.exp(-T[:,None]*V.Q[None,:]**2)
                        extra=(25.*g['factor'] if world=='source_noise_shift' else np.full(len(T),25.))
                        f=np.exp(-.5*extra[:,None]*V.Q[None,:]**2) if world in ('noise_shift','source_noise_shift') else np.ones_like(att)
                        for w,triples in enumerate(V.TRIPLES):
                            weights=sum(g['lc'][:,i]*g['lc'][:,j]*g['lc'][:,k]/2 for i,j,k in triples)
                            base=K*weights[:,None]
                            value[0,:,w,1]=base.sum(0);value[0,:,w,0]=(base*att).sum(0)
                            value[1,:,w,1]=(base*(f**2 if w==1 else 1)).sum(0)
                            value[1,:,w,0]=(base*f*att).sum(0)
                    independent.append(value)
                independent=np.stack(independent)
                assert np.allclose(stored['by_template'],independent,rtol=2e-12,atol=1e-8),(tag,model,world)
                assert np.allclose(stored['population'],independent.mean(0),rtol=2e-12,atol=1e-8)
                truths+=1
    # Replay the first sixteen groups of one frozen field per view/noise law.
    # This audits recorded seeds and observation provenance; law controls are separate.
    replayed=0
    for tag in V.TAGS:
        template,meta,_=V.templates(tag)
        for model in V.MODELS:
            name=f'{tag}-{model}-homogeneous';folder=V.P/'worlds'/name
            recorded=np.load(folder/'field-000.npz')
            rng=np.random.default_rng(V.seed(name,0,'physical'))
            code_rngs=[np.random.default_rng(V.seed(name,0,'codes',a)) for a in range(2)]
            ids=rng.integers(meta['groups'],size=max(meta['groups'],512))
            assert np.array_equal(ids,recorded['template_ids'])
            z=np.zeros((2,16,4,4,2,4),complex);counts=np.zeros((2,16,4,4));known=np.zeros((16,2,2,2))
            for j,idx in enumerate(ids[:16]):
                g=template[idx];ob=V.simulate_group(g,model,'homogeneous',rng)
                if ob is None:continue
                known[j]=V.diagonal(ob['z'],ob['n'])
                for a,crng in enumerate(code_rngs):
                    codes=crng.integers(4,size=len(g['n']))
                    z[a,j],counts[a,j]=V.aggregate_observations(ob['z'],ob['n'],codes)
            coded=np.stack([V.estimate(z[a],counts[a]) for a in range(2)])
            assert np.allclose(coded,recorded['coded'][:,:16],rtol=1e-13,atol=1e-8),(tag,model,'coded seed replay')
            assert np.array_equal(known,recorded['known_source'][:16]),(tag,model,'source seed replay')
            replayed+=1
    # Recompute means/covariances, confidence polynomial and D without estimate()/fieller().
    for folder in sorted((V.P/'worlds').iterdir()):
        state=json.loads((folder/'state.json').read_text())
        assert state['complete'] and state['expected']==state['completed']
        for entry in state['replicates']:
            rep=entry['rep'];a=np.load(folder/entry['filename']);fields+=1
            record=json.loads((folder/f'field-{rep:03d}.json').read_text())
            for r in record['cases']:
                array=a['coded'][r['assignment']] if r['method']=='coded' else a['known_source']
                iq=int(np.flatnonzero(V.Q==r['q'])[0]);group=array[:r['group_count'],iq]
                values=group.mean(1);x,y=values[:,0],values[:,1]
                mx,my=x.mean(),y.mean();N=len(x);crit=t.ppf(.975,N-1)**2
                A=my**2-crit*y.var(ddof=1)/N
                B=-2*mx*my+2*crit*np.cov(x,y,ddof=1)[0,1]/N
                C=mx**2-crit*x.var(ddof=1)/N
                for key,value in [('A',A),('B',B),('C',C),('numerator',mx),('reference',my),('ratio',mx/my if my>0 else np.nan)]:
                    if not equal(value,r[key],1e-8):bad.append((folder.name,rep,iq,key))
                disc=B*B-4*A*C
                if A>0: kind='bounded' if disc>=0 else 'empty'
                elif A<0: kind='disjoint' if disc>0 else 'all_real'
                else: kind='half_line' if B!=0 else ('all_real' if C<=0 else 'empty')
                if kind!=r['kind']:bad.append((folder.name,rep,iq,'set_type'))
                if disc>=0 and A!=0:
                    bounds=sorted(((-B-np.sqrt(disc))/(2*A),(-B+np.sqrt(disc))/(2*A)))
                    for key,value in zip(('lower','upper'),bounds):
                        if not equal(value,r[key],1e-8):bad.append((folder.name,rep,iq,key))
                    usable=kind=='bounded' and bounds[0]>0 and my>0
                    if bool(usable)!=r['usable']:bad.append((folder.name,rep,iq,'usability'))
                contrast=x-r['physical_target']*y
                covers=bool(contrast.mean()**2<=crit*contrast.var(ddof=1)/N+1e-12*max(mx**2,my**2))
                if covers!=r['covers']:bad.append((folder.name,rep,iq,'covers'))
                moments=group.mean(0);ra=moments[0,0]/moments[0,1];rb=moments[1,0]/moments[1,1]
                valid=moments[0,1]>0 and moments[1,1]>0 and ra>0 and rb>0
                if valid:
                    influence=(group[:,0,0]/moments[0,0]-group[:,0,1]/moments[0,1])-(group[:,1,0]/moments[1,0]-group[:,1,1]/moments[1,1])
                    se=influence.std(ddof=1)/np.sqrt(N)
                    z=np.log(ra/rb)/se if se>0 else (0. if ra==rb else np.inf)
                    if not equal(z,r['D']['z'],1e-8):bad.append((folder.name,rep,iq,'D_z'))
                    if bool(abs(z)<=t.ppf(.975,N-1))!=r['D']['quiet']:bad.append((folder.name,rep,iq,'D_quiet'))
                elif r['D']['valid']:bad.append((folder.name,rep,iq,'D_invalid'))
                cases+=1
    assert not bad,bad[:20]
    assert fields==3528 and cases==42336 and truths==126 and catalogs==9,(fields,cases,truths,catalogs)
    # Recount all summary fractions directly from field record manifests.
    for s in summary['records']:
        folder=V.P/'worlds'/f"{s['tag']}-{s['model']}-{s['world']}"
        rows=[]
        for path in sorted(folder.glob('field-*.json')):
            rows.extend(r for r in json.loads(path.read_text())['cases'] if r['method']==s['method'] and r['scope']==s['scope'] and r['q']==s['q'])
        assert len(rows)==s['records']
        assert equal(np.mean([r['usable'] for r in rows]),s['usable_fraction'])
        assert equal(np.mean([r['usable'] and r['covers'] for r in rows]),s['usable_joint_coverage'])
        assert equal(np.mean([r['covers'] for r in rows]),s['unconditional_coverage'])
    save(V.P/'population-information.json',dict(conditional_only=True,records=information))
    result=dict(all_pass=True,dependency_freeze_commit='44b98a8',dependencies_sha256=dependencies,
                runtime_versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),conditional_only=True,physical_validation=False,
                distinct_acquisitions=8,views=9,physical_fields=fields,coded_fields=7056,
                case_recounts=cases,truth_worlds=truths,source_catalogs=catalogs,
                verified_artifacts=len(manifest['artifacts']),summary_records=len(summary['records']),seed_replayed_fields=replayed,seed_replayed_groups=replayed*16,
                artifact_sha256={name:sha(V.P/name) for name in ('summary.json','summary.csv','artifact-manifest.json','controls.json','freeze.json','PROTOCOL.md','validation.py','summarize.py','verify.py','population-information.json','execution-state.json','execution.log')})
    save(V.P/'verification.json',result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
