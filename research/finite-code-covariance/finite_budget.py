"""Finite-budget verification under the frozen sixteen-state observation law."""
from pathlib import Path
import argparse,itertools,json,sys,time
import numpy as np
from scipy.stats import t as student_t
from covariance import four_panel,gradients
from run import ROOT,OLD,V,save,sha,replay_close
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.uncertainty import _quadratic_set
from codedcal.geometric_uncertainty import _positive_parts,_product_root,_contains


def ratio_from_moments(mean,cov,n,confidence=.95):
    x,y=mean
    scale=max(abs(x),abs(y),np.sqrt(np.abs(cov).max()),1e-300)
    x,y=mean/scale;v=cov/scale**2
    critical=student_t.ppf((1+confidence)/2,n-1)**2
    a=float(y*y-critical*v[1,1]);b=float(-2*x*y+2*critical*v[0,1]);c=float(x*x-critical*v[0,0])
    kind,parts=_quadratic_set(a,b,c)
    return dict(kind=kind,intervals=parts,usable=bool(kind=='bounded' and parts[0]['lower']>0 and y>0),
                ratio=float(x/y) if y!=0 else None,positive_reference=bool(y>0))


def summarize_case(mu,cov,n,variant,target):
    # cov is the estimated covariance of the mean, not of individual groups.
    if variant=='arithmetic':
        transform=np.array([[.5,0,.5,0],[0,.5,0,.5]])
        f=ratio_from_moments(transform@mu,transform@cov@transform.T,n)
        parts=f['intervals'];available=f['positive_reference'] and f['ratio'] is not None and f['ratio']>0
    else:
        windows=[ratio_from_moments(mu[2*w:2*w+2],cov[2*w:2*w+2,2*w:2*w+2],n,.975) for w in range(2)]
        parts=_product_root(*(_positive_parts(w['intervals']) for w in windows))
        available=all(w['positive_reference'] and w['ratio'] is not None and w['ratio']>0 for w in windows)
        f=dict(usable=all(w['usable'] for w in windows),ratio=float(np.exp(np.mean([np.log(w['ratio']) for w in windows]))) if available else None,
               kind=('empty' if not parts else 'disjoint' if len(parts)>1 else 'unbounded_positive' if parts[0]['upper'] is None else 'bounded'),
               window_kinds=[w['kind'] for w in windows])
    width=(parts[0]['upper']-parts[0]['lower']) if f['usable'] and len(parts)==1 and parts[0]['upper'] is not None else None
    return dict(available=available,ratio=f['ratio'] if available else None,kind=f['kind'],usable=f['usable'],
                covers=_contains(parts,target),width=width,
                usable_covers=bool(f['usable'] and _contains(parts,target)),
                usable_covers_excludes_one=bool(f['usable'] and _contains(parts,target) and not _contains(parts,1.)),
                error=(f['ratio']-target) if available else None,
                intervals=parts,window_kinds=f.get('window_kinds'))


def reconstruct(name,prepared):
    tag,model,world=name.rsplit('-',2);groups,meta,_=prepared
    rng=np.random.default_rng(V.seed(name,0,'physical'))
    chosen=rng.integers(meta['groups'],size=max(meta['groups'],512))
    crng=[np.random.default_rng(V.seed(name,0,'codes',a)) for a in range(2)]
    with np.load(OLD/'worlds'/name/'field-000.npz') as f:
        np.testing.assert_array_equal(chosen,f['template_ids']);diag=f['known_source'][:16].reshape(16,8).copy();coded=f['coded'][:,:16].copy()
    contrasts=[]
    for j,idx in enumerate(chosen[:16]):
        ob=V.simulate_group(groups[idx],model,world,rng)
        if ob is None:contrasts.append(None);continue
        c=four_panel(ob['z'],ob['n']);contrasts.append(c)
        for a in range(2):replay_close(c.evaluate(crng[a].integers(4,size=len(ob['n']))),coded[a,j].reshape(8))
        replay_close(c.moments()['mean'],diag[j])
    return contrasts,diag


def compute(name,prepared,prefix,out):
    contrasts,diag=reconstruct(name,prepared)
    mu=diag.mean(0);center=diag-mu
    floor=center.T@center/16
    pair=np.array(prefix['pair_per_group']);triple=np.array(prefix['triple_per_group'])
    result=dict(name=name,law='uniform on fixed sixteen observation groups, fresh iid codes',
                mean=mu,physical_covariance=floor,pair=pair,triple=triple,
                prefix_sources_sha256=prefix['source_hashes'],budget=256,replicates=200,comparisons=[])
    arrays={}
    for k in (4,8,16):
        n=256//k;rng=np.random.default_rng(V.seed('finite-budget-oct27',name,k))
        chosen=rng.integers(16,size=(200,n));observations=np.zeros((200,n,8))
        for j,c in enumerate(contrasts):
            if c is None:continue
            reps,slots=np.where(chosen==j)
            assignments=rng.integers(k,size=(len(reps),len(c.f)))
            for lo in range(0,len(reps),32):
                observations[reps[lo:lo+32],slots[lo:lo+32]]=c.evaluate(assignments[lo:lo+32],k)
        means=observations.mean(1);deviations=observations-means[:,None]
        covs=np.einsum('rgi,rgj->rij',deviations,deviations)/(n*(n-1))
        arrays[f'K{k}_means']=means;arrays[f'K{k}_mean_covariances']=covs
        # Retain one full experiment to independently check sufficient-statistic
        # interval evaluation against the shipped raw-group API.
        arrays[f'K{k}_first_group_statistics']=observations[0]
        predicted=(floor+pair/(k-1)+triple/((k-1)*(k-2)))/n
        observed=(means-mu).T@(means-mu)/200
        for iq,q in enumerate(V.Q):
            sl=slice(iq*4,iq*4+4);g=gradients(mu[sl])
            for variant in ('arithmetic','geometric'):
                if g is None:
                    result['comparisons'].append(dict(K=k,G=n,q=float(q),variant=variant,supported=False));continue
                a=mu[sl]
                target=(a[0]+a[2])/(a[1]+a[3]) if variant=='arithmetic' else np.sqrt((a[0]/a[1])*(a[2]/a[3]))
                cases=[summarize_case(m[sl],cov[sl,sl],n,variant,target) for m,cov in zip(means,covs)]
                v=g[variant]
                residual=(means[:,sl]-a)@v;squared=residual**2
                errors=np.array([r['error'] for r in cases if r['error'] is not None])
                widths=[r['width'] for r in cases if r['width'] is not None]
                result['comparisons'].append(dict(K=k,G=n,q=float(q),variant=variant,supported=True,target=float(target),
                    predicted_influence_variance=float(v@predicted[sl,sl]@v),
                    observed_influence_second_moment=float(squared.mean()),
                    influence_second_moment_mc_se=float(squared.std(ddof=1)/np.sqrt(200)),
                    available=sum(r['available'] for r in cases),usable=sum(r['usable'] for r in cases),
                    covers=sum(r['covers'] for r in cases),usable_covers=sum(r['usable_covers'] for r in cases),
                    usable_covers_excludes_one=sum(r['usable_covers_excludes_one'] for r in cases),
                    finite_available_bias=float(errors.mean()) if len(errors) else None,
                    finite_available_rmse=float(np.sqrt(np.mean(errors**2))) if len(errors) else None,
                    finite_available_abs_error_quantiles=np.quantile(np.abs(errors),[.5,.9,.99,1]) if len(errors) else None,
                    median_usable_width=float(np.median(widths)) if widths else None,cases=cases))
    np.savez_compressed(out/f'{name}.npz',**arrays)
    result['sufficient_statistics_sha256']=sha(out/f'{name}.npz')
    save(out/f'{name}.json',result)


def main():
    p=argparse.ArgumentParser();p.add_argument('--prefix',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    prefix=Path(args.prefix);out=Path(args.output)
    if out.exists():raise ValueError('new output directory required')
    out.mkdir(parents=True);completed=[];start=time.time()
    try:
        for tag in V.TAGS:
            prepared=V.templates(tag)
            for model,world in itertools.product(V.MODELS,tuple(w for w in V.WORLDS if w!='code_offset')):
                name=f'{tag}-{model}-{world}';tick=time.time()
                source=json.loads((prefix/f'{name}.json').read_text())
                compute(name,prepared,source,out);completed.append(name)
                save(out/'state.json',dict(stage='running',completed=completed,seconds=time.time()-start))
                print(name,'done',round(time.time()-tick,2),'seconds',flush=True)
        save(out/'state.json',dict(stage='complete',completed=completed,seconds=time.time()-start,
                                 implementation_sha256=sha(__file__)))
    except Exception as e:
        save(out/'state.json',dict(stage='failed',completed=completed,error=repr(e),seconds=time.time()-start));raise

if __name__=='__main__':main()
