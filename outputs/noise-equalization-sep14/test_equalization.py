"""Pilot: learn relative noise from code/quality strata, then lower moment order.

Quality is an observed covariate, not source identity. Its noise coefficient is
learned on independent groups. Misspecification is explicitly stress-tested.
"""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'relative-code-calibration-sep14'))
import experiment as old
Q=old.Q; H=old.H; LEVELS=np.arange(6)

def json_safe(x):
    if isinstance(x,dict):return {k:json_safe(v) for k,v in x.items()}
    if isinstance(x,list):return [json_safe(v) for v in x]
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def pair_components(z,n):
    # z,n: group,code,panel. Counts have no q dependence.
    s=np.einsum('gct,cb->gtb',z,H)
    k=np.einsum('gct,cb->gtb',n,H)
    def moments(pairs):
        phase=sum((s[:,a]*s[:,b].conj()).mean(-1).real for a,b in pairs)/len(pairs)
        count=sum((k[:,a]*k[:,b]).mean(-1) for a,b in pairs)/len(pairs)
        return phase,count
    cross,nc=moments([(0,2),(0,3),(1,2),(1,3)])
    within,nw=moments([(0,1),(2,3)])
    return np.stack([cross,nc,within,nw],-1)

def generate(sources,gamma,tau,unmodeled,seed,groups):
    rng=np.random.default_rng(seed);zs=[];ns=[]
    for start in range(0,groups,256):
        g=min(256,groups-start);shape=(g,sources)
        d=np.broadcast_to(np.linspace(-1,1,sources),shape)
        base=np.exp(rng.normal(np.log(8),.35,shape))
        lam=base[...,None]*(1+.85*d[...,None]*np.array([1,1,-1,-1]))*np.array([.8,1.2,1.1,.9])
        count=rng.poisson(lam);m=count.max();mask=np.arange(m)<count[...,None]
        shape4=shape+(4,m)
        # The observed quality varies between binding events at the same source.
        quality=np.clip(np.rint(4*d[...,None,None]**2+rng.uniform(-1,1.5,size=shape4)),0,5).astype(int)
        event_var=25-gamma**2-tau**2+10*quality+unmodeled*d[...,None,None]
        assert event_var.min()>0
        cg=np.repeat(rng.normal(0,gamma,size=(g,1,2)),2,-1)
        ui=np.repeat(rng.normal(0,tau,size=shape+(2,)),2,-1)
        x=35*d[...,None,None]+cg[...,None]+ui[...,None]+np.sqrt(event_var)*rng.normal(size=shape4)
        weight=mask*rng.integers(1,6,size=shape4)
        codes=rng.integers(4,size=shape)
        zb=np.empty((g,4,4,6,len(Q)),complex);nb=np.empty((g,4,4,6))
        for lev in LEVELS:
            w=weight*(quality==lev)
            nb[:,:,:,lev]=old.aggregate(w.sum(-1),codes)
            for iq,q in enumerate(Q):
                raw=(w*np.exp(1j*q*x)).sum(-1)
                zb[:,:,:,lev,iq]=old.aggregate(raw,codes)
        zs.append(zb);ns.append(nb)
    return np.concatenate(zs),np.concatenate(ns)

def calibration(z,n):
    a=np.empty((len(z),len(Q),6,2))
    for j in range(len(Q)):
        for k in LEVELS:
            m=pair_components(z[:,:,:,k,j],n[:,:,:,k])
            a[:,j,k]=m[:,[2,3]]
    blocks=a.reshape(-1,256,len(Q),6,2).sum(1)
    def fit(b):
        s=b.sum(0);ratio=s[...,0]/s[...,1]
        # Positive expected count mass determines included quality strata.
        valid=(s[...,1]>0)&np.isfinite(ratio)
        model=lambda p:np.exp(-Q[:,None]**2*(p[0]+p[1]*LEVELS))
        sol=least_squares(lambda p:(model(p)-ratio)[valid],x0=[15,10],bounds=([-30,-5],[100,30]))
        return sol.x
    pars=fit(blocks)
    rng=np.random.default_rng(2609144100)
    boot=np.array([fit(blocks[rng.integers(len(blocks),size=len(blocks))]) for _ in range(200)])
    return dict(intercept=float(pars[0]),beta=float(pars[1]),beta_se=float(boot[:,1].std(ddof=1)),
        beta_interval=np.quantile(boot[:,1],[.025,.975]).tolist(),
        note='Independent calibration groups; intercept is fitted but not used for equalization.'),blocks

def pair_fit(a):
    # Independent blocks,q,4 moments. Ratio of separately normalized kernels.
    m=a.mean(0)
    ratio=(m[:,0]/m[:,1])/(m[:,2]/m[:,3])
    fun=lambda t:np.sum((ratio-np.exp(-Q**2*t))**2)
    t=float(minimize_scalar(fun,bounds=(-30,60),method='bounded').x)
    d=-Q**2*np.exp(-Q**2*t)
    d2=Q**4*np.exp(-Q**2*t)
    hess=np.sum(d*d-(ratio-np.exp(-Q**2*t))*d2)
    jac=np.stack([m[:,3]/m[:,1]/m[:,2],-ratio/m[:,1],-ratio/m[:,2],m[:,0]/m[:,1]/m[:,2]],-1)
    inf=np.sum((a-m)*jac,axis=-1)@d/hess
    return t,inf,ratio

def features(z,n,beta):
    counts=n.sum(-1)
    answer=[]
    for j,q in enumerate(Q):
        amp=np.exp(-.5*q*q*beta*(5-LEVELS))
        answer.append(pair_components((z[...,j]*amp).sum(-1),counts))
    return np.stack(answer,1)

def summarize(pair,quartic,beta_se=0.,sensitivity=0.):
    t,ti,ratio=pair_fit(pair)
    c,ci=old.fit(quartic[...,2],quartic[...,3])
    val=t-c['value']
    se=float(np.sqrt(np.var(ti-ci,ddof=1)/len(pair)+(beta_se*sensitivity)**2))
    return dict(total=t,common=c['value'],independent=val,independent_se=se,
        interval=[val-1.96*se,val+1.96*se],ratio=ratio.tolist(),
        boundary=bool(t< -29.9 or t>59.9 or c['boundary']))

def main():
    configs=[dict(sources=32,gamma=g,tau=t,unmodeled=u,seed=2609144101+j)
        for j,(g,t,u) in enumerate([(0,0,0),(3,0,0),(0,3,0),(3,3,0),(0,0,8),(3,3,8)])]
    design=dict(configs=configs,training_groups=4096,test_groups=16384,q=Q.tolist(),
        fixed_quality_noise_slope_truth=10,quality='Observed event-varying level; no source identities supplied.',
        primary='Learn relative noise slope on independent groups; bounded noise equalization followed by pair-order temporal contrast.',
        stress='Unobserved source-specific variance 8*direction, correlated with changing yields.',
        status='Synthetic pilot with explicit extra observable; not a real-data or novelty claim.')
    (P/'design.json').write_text(json.dumps(design,indent=2)+'\n')
    out=[];start=time.time()
    for j,cfg in enumerate(configs):
        args={k:v for k,v in cfg.items() if k!='seed'}
        z,n=generate(**args,seed=cfg['seed'],groups=4096)
        cal,calblocks=calibration(z,n)
        assert cal['beta']>0, 'Learned negative noise slope: retain and investigate, do not silently clip.'
        z,n=generate(**args,seed=cfg['seed']+100,groups=16384)
        quartic=np.stack([old.components(z[...,iq].sum(-1)) for iq in range(len(Q))],1)
        qb=quartic.reshape(-1,1024,len(Q),4).sum(1)
        beta=cal['beta']
        methods={};saved=dict(quartic=quartic,calibration=calblocks,q=Q)
        plus=features(z,n,beta+.01).reshape(-1,1024,len(Q),4).sum(1)
        minus=features(z,n,beta-.01).reshape(-1,1024,len(Q),4).sum(1)
        sensitivity=(pair_fit(plus)[0]-pair_fit(minus)[0])/.02
        saved['beta_plus']=plus;saved['beta_minus']=minus
        for name,b in [('raw_pairs',0),('learned_equalization',beta),('known_relative_noise',10),('wrong_slope',5)]:
            f=features(z,n,b);blocks=f.reshape(-1,1024,len(Q),4).sum(1)
            bs=cal['beta_se'] if name=='learned_equalization' else 0
            stats=summarize(blocks,qb,bs,sensitivity)
            fits=[summarize(x,y) for x,y in zip(f.reshape(-1,1024,len(Q),4),quartic.reshape(-1,1024,len(Q),4))]
            vals=np.array([x['independent'] for x in fits])
            stats.update(rmse_1024=float(np.sqrt(np.mean((vals-cfg['tau']**2)**2))),
                block_estimates=vals.tolist(),boundary_count=sum(x['boundary'] for x in fits))
            methods[name]=stats;saved[name]=f
        oldfits=[old.summarize(x)['independent']['value'] for x in quartic.reshape(-1,1024,len(Q),4)]
        methods['original_quartic']=dict(pooled=old.summarize(qb),block_estimates=oldfits,
            rmse_1024=float(np.sqrt(np.mean((np.array(oldfits)-cfg['tau']**2)**2))))
        result=dict(**cfg,calibration=cal,calibration_sensitivity=sensitivity,methods=methods)
        out.append(result);np.savez_compressed(P/f'case_{j}.npz',**saved)
        (P/'results.json').write_text(json.dumps(json_safe(out),indent=2,allow_nan=False)+'\n')
        print(json.dumps({**cfg,'calibration':cal,'estimates':{k:round(v.get('independent',v.get('pooled',{}).get('independent',{}).get('value',0)),3) for k,v in methods.items()},'rmse1024':{k:v['rmse_1024'] for k,v in methods.items()}}),flush=True)
    print('DONE seconds',time.time()-start,flush=True)

if __name__=='__main__':main()
