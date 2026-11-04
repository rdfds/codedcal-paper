"""Assignment-free self/distinct temporal calibration; synthetic mechanism test.

Four random physical codes must be fixed per source and independent of its trace.
Source identities are used only in oracle audit arrays, never in code estimates.
"""
from pathlib import Path
import sys, itertools, json, time
import numpy as np
from scipy.optimize import minimize_scalar

P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent/'code-time-identification-sep14'))
from matched_code_moments import H, aggregate, diagonal_moments
Q = np.array([.10, .15])

def components(coded):
    f = np.concatenate([coded, coded.conj()], axis=-1)
    v = np.einsum('...cf,cb->...fb', f, H)
    _, _, t4 = diagonal_moments(f)
    def t22(a,b,c,d):
        mixed = sum(v[...,a,u]*v[...,c,u]*v[...,b,w]*v[...,d,w]
                    for u,w in itertools.permutations(range(3),2))/6
        return mixed-t4(a,b,c,d)
    return np.stack([t4(0,1,6,7).real, t4(0,5,2,7).real,
                     t22(0,1,6,7).real, t22(0,5,6,3).real], axis=-1)

def oracle(z):
    a,b,c,d = np.moveaxis(z,-1,0)
    def distinct(a,b,c,d):
        return ((a*c).sum(-1)*(b*d).sum(-1)-(a*b*c*d).sum(-1)).real
    return np.stack([(a*b*c.conj()*d.conj()).sum(-1).real,
                     (a*b.conj()*c*d.conj()).sum(-1).real,
                     distinct(a,b,c.conj(),d.conj()),
                     distinct(a,b.conj(),c.conj(),d)],axis=-1)

def exact_check():
    rng=np.random.default_rng(2609143100)
    z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    codes=np.array(list(itertools.product(range(4),repeat=4)))
    got=components(aggregate(np.broadcast_to(z,(len(codes),)+z.shape),codes)).mean(0)
    want=oracle(z)
    err=np.abs(got-want)
    assert err.max()<1e-10, (got,want)
    return dict(code_assignments=len(codes),max_absolute_error=float(err.max()),
                estimated=got.tolist(),source_oracle=want.tolist())

def fit(c,r):
    # c,r shape independent block,q. Signed parameters and invalid references
    # remain visible. Intervals are a first-order approximation, not certification.
    ratio=c.sum(0)/r.sum(0)
    obj=lambda t: np.sum((ratio-np.exp(-4*Q**2*t))**2)
    t=float(minimize_scalar(obj,bounds=(-30,60),method='bounded').x)
    d=-4*Q**2*np.exp(-4*Q**2*t)
    residual=ratio-np.exp(-4*Q**2*t)
    d2=16*Q**4*np.exp(-4*Q**2*t)
    curvature=np.sum(d*d-residual*d2)
    influence=((c-r*ratio)/r.mean(0))@d/curvature
    se=float(influence.std(ddof=1)/np.sqrt(len(c)))
    return dict(value=t,se=se,interval=[t-1.96*se,t+1.96*se],ratio=ratio.tolist(),
                frequency_estimates=(-np.log(np.where(ratio>0,ratio,np.nan))/(4*Q**2)).tolist(),
                reference_mean=r.mean(0).tolist(),
                boundary=bool(t < -29.9 or t > 59.9)),influence

def summarize(arr):
    # arr shape independent blocks,q,4 component sums
    total,ti=fit(arr[...,0],arr[...,1])
    common,ci=fit(arr[...,2],arr[...,3])
    relative=total['value']-common['value']
    se=float((ti-ci).std(ddof=1)/np.sqrt(len(ti)))
    return dict(total=total,common=common,
                independent=dict(value=relative,se=se,interval=[relative-1.96*se,relative+1.96*se]),
                relative_pair_variance=dict(value=2*relative,se=2*se,
                    interval=[2*relative-3.92*se,2*relative+3.92*se]))

def simulate(sources,gamma,tau,seed,R=32,G=1024):
    rng=np.random.default_rng(seed)
    # Keep individual event marginals identical across covariance worlds.
    code_blocks=[];oracle_blocks=[];group_sums=[]
    for rep in range(R):
        shape=(G,sources)
        direction=np.broadcast_to(np.linspace(-1,1,sources),shape)
        theta=35*direction;v=25+40*direction**2
        base=np.exp(rng.normal(np.log(8),.35,shape))
        lam=base[...,None]*(1+.85*direction[...,None]*np.array([1,1,-1,-1]))*np.array([.8,1.2,1.1,.9])
        n=rng.poisson(lam);maxn=n.max();mask=np.arange(maxn)<n[...,None]
        common=np.repeat(rng.normal(0,gamma,size=(G,1,2)),2,axis=-1)
        individual=np.repeat(rng.normal(0,tau,size=shape+(2,)),2,axis=-1)
        x=theta[...,None,None]+common[...,None]+individual[...,None]+np.sqrt(v-gamma**2-tau**2)[...,None,None]*rng.normal(size=shape+(4,maxn))
        weighted=mask*rng.integers(1,6,size=x.shape)
        codes=rng.integers(4,size=shape)
        ca=[];oa=[]
        for q in Q:
            z=(np.exp(1j*q*x)*weighted).sum(-1)
            ca.append(components(aggregate(z,codes)))
            oa.append(oracle(z))
        ca=np.stack(ca,axis=1);oa=np.stack(oa,axis=1)
        group_sums.append(ca);code_blocks.append(ca.sum(0));oracle_blocks.append(oa.sum(0))
    return dict(code=np.array(code_blocks),oracle=np.array(oracle_blocks),groups=np.concatenate(group_sums))

def budget_audit(groups):
    out={}
    for g in (256,1024,4096):
        # Non-overlapping data blocks; uncertainty within each comes from groups.
        fits=[summarize(x) for x in groups.reshape(-1,g,len(Q),4)]
        out[str(g)]={}
        for name in ('total','common','independent','relative_pair_variance'):
            vals=np.array([f[name]['value'] for f in fits])
            out[str(g)][name]=dict(median=float(np.median(vals)),
                percentile_5_95=np.quantile(vals,[.05,.95]).tolist(),
                estimates=vals.tolist())
        out[str(g)]['fit_boundary_count']=sum(f['total']['boundary'] or f['common']['boundary'] for f in fits)
    return out

def main():
    configs=[dict(sources=s,gamma=g,tau=t,seed=2609143101+j*4+k)
             for j,s in enumerate((8,32)) for k,(g,t) in enumerate(((0,0),(3,0),(0,3),(3,3)))]
    design=dict(q=Q.tolist(),replicates=32,groups_per_replicate=1024,
        configs=configs,assumptions=['uniform independent fixed source code',
        'Gaussian shared and independent window shifts independent of yields',
        'symmetric independent event errors between four panels'],
        primary='separate common and independent persistent variance at identical marginal event variance',
        limitations='Synthetic mechanism test; not real molecular recovery or novelty proof.')
    (P/'design.json').write_text(json.dumps(design,indent=2)+'\n')
    out=dict(exact_identity=exact_check(),results=[])
    start=time.time()
    for j,cfg in enumerate(configs):
        a=simulate(**cfg)
        np.savez_compressed(P/f'world_{j}.npz',**a,q=Q)
        result=dict(**cfg,truth=dict(total=cfg['gamma']**2+cfg['tau']**2,
                    common=cfg['gamma']**2,independent=cfg['tau']**2,relative_pair_variance=2*cfg['tau']**2),
                    code=summarize(a['code']),oracle=summarize(a['oracle']),budget=budget_audit(a['groups']))
        out['results'].append(result)
        (P/'results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({**cfg,'code':result['code'],'oracle':result['oracle']}),flush=True)
    print('DONE seconds',time.time()-start,flush=True)

if __name__=='__main__':main()
