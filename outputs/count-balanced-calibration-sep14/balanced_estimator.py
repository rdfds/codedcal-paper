"""Count-balanced cubic temporal contrast; source identities used only by generator.

Both sides of each contrast involve the same three panels and quality stratum.
Two panels supply phase, the other count. Swap phase/count roles to compare
within- and cross-window motion at exactly matched population source weights.
"""
from pathlib import Path
import itertools,json,sys,time
import numpy as np
from scipy.optimize import minimize_scalar
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'noise-equalization-sep14'))
from test_equalization import H,Q,old,pair_components,json_safe
PERMS=list(itertools.permutations(range(3)))
TRIPLES=[(0,1,2),(0,1,3),(2,3,0),(2,3,1)]

def cubic(z,n):
    # group,code,panel. No source or quality metadata input.
    a=np.einsum('gct,cb->gtb',z,H);b=np.einsum('gct,cb->gtb',n,H)
    def t3(x,y,z):return sum(x[:,p[0]]*y[:,p[1]]*z[:,p[2]] for p in PERMS).real/6
    cross=0.;within=0.
    for i,j,k in TRIPLES:
        within=within+t3(a[:,i],a[:,j].conj(),b[:,k])/4
        cross=cross+(t3(a[:,i],b[:,j],a[:,k].conj())+t3(b[:,i],a[:,j],a[:,k].conj()))/8
    return np.stack([cross,within],-1)

def fit_ratio(a):
    m=a.mean(0);ratio=m[:,0]/m[:,1]
    t=float(minimize_scalar(lambda t:np.sum((ratio-np.exp(-Q**2*t))**2),bounds=(-30,60),method='bounded').x)
    return t,dict(ratio=ratio.tolist(),reference=m[:,1].tolist(),boundary=bool(t< -29.9 or t>59.9))

def strat_fit(a):
    m=a.mean(0);c,nc,w,nw=np.moveaxis(m,-1,0)
    valid=(nc>0)&(nw>0);weight=np.sqrt(np.maximum(nc,0)*np.maximum(nw,0))
    cw=np.sum(np.where(valid,c/np.where(valid,nc,1)*weight,0),-1)
    ww=np.sum(np.where(valid,w/np.where(valid,nw,1)*weight,0),-1)
    t,det=fit_ratio(np.stack([cw,ww],-1)[None])
    det['excluded_strata']=int((~valid).sum())
    return t,det

def generate(seed,groups=8192,sources=32,gamma=3,tau=3,hidden=40,quality_shift=0,hidden_shift=0,quality_source=False,keep_raw=False):
    rng=np.random.default_rng(seed);out=[]
    for start in range(0,groups,256):
        g=min(256,groups-start);shape=(g,sources)
        d=np.broadcast_to(np.linspace(-1,1,sources),shape)
        base=np.exp(rng.normal(np.log(8),.35,shape))
        lam=base[...,None]*(1+.85*d[...,None]*np.array([1,1,-1,-1]))*np.array([.8,1.2,1.1,.9])
        count=rng.poisson(lam);m=count.max();shape4=shape+(4,m)
        valid=np.arange(m)<count[...,None]
        # h unrelated to source in the static counterexample; varying h across
        # windows tests conditional cubic correction when pooled cubic fails.
        h0=(np.rint(4*d[...,None,None]**2+rng.uniform(-1,1.5,size=shape4)).astype(int)
            if quality_source else rng.integers(0,6,size=shape4))
        h=np.clip(h0+quality_shift*np.array([0,0,1,1])[None,None,:,None],0,5)
        ev=25-gamma**2-tau**2+10*h+hidden*d[...,None,None]**2+hidden_shift*d[...,None,None]**2*np.array([0,0,1,1])[None,None,:,None]
        assert ev.min()>0
        cg=np.repeat(rng.normal(0,gamma,size=(g,1,2)),2,-1)
        ui=np.repeat(rng.normal(0,tau,size=shape+(2,)),2,-1)
        x=35*d[...,None,None]+cg[...,None]+ui[...,None]+np.sqrt(ev)*rng.normal(size=shape4)
        w=valid*rng.integers(1,6,size=shape4);colors=rng.integers(4,size=shape)
        z=np.zeros((g,4,4,6,len(Q)),complex);n=np.zeros((g,4,4,6))
        phases=np.exp(1j*x[...,None]*Q)
        for lev in range(6):
            ww=w*(h==lev);n[:,:,:,lev]=old.aggregate(ww.sum(-1),colors)
            for iq in range(len(Q)):z[:,:,:,lev,iq]=old.aggregate((ww*phases[...,iq]).sum(-1),colors)
        a={}
        a['quartic']=np.stack([old.components(z[...,j].sum(-1)) for j in range(len(Q))],1)
        a['cubic_pooled']=np.stack([cubic(z[...,j].sum(-1),n.sum(-1)) for j in range(len(Q))],1)
        a['cubic_quality']=np.stack([sum(cubic(z[:,:,:,h,j],n[:,:,:,h]) for h in range(6)) for j in range(len(Q))],1)
        a['pair_quality']=np.stack([np.stack([pair_components(z[:,:,:,h,j],n[:,:,:,h]) for h in range(6)],1) for j in range(len(Q))],1)
        if keep_raw:a.update(z=z,n=n)
        out.append(a)
    return {k:np.concatenate([a[k] for a in out]) for k in out[0]}

def summarize(data,tau):
    out={};q=data['quartic'];common=old.fit(q[...,2],q[...,3])[0]['value']
    for name in ['cubic_pooled','cubic_quality','pair_quality','quartic']:
        a=data[name]
        if name=='quartic':val=old.summarize(q)['independent']['value'];det={}
        else:
            total,det=(strat_fit(a) if name=='pair_quality' else fit_ratio(a));val=total-common
        vals=[];details=[]
        for i in range(0,len(q),1024):
            b=a[i:i+1024];qb=q[i:i+1024]
            if name=='quartic':v=old.summarize(qb)['independent']['value'];d={}
            else:
                t,d=(strat_fit(b) if name=='pair_quality' else fit_ratio(b))
                v=t-old.fit(qb[...,2],qb[...,3])[0]['value']
            vals.append(v);details.append(d)
        out[name]=dict(pooled=val,pooled_details=det,block_estimates=vals,
            rmse=float(np.sqrt(np.mean((np.array(vals)-tau*tau)**2))),block_details=details)
    return out

def exact():
    rng=np.random.default_rng(2609145100);raw=rng.normal(size=(3,4))+1j*rng.normal(size=(3,4));n=rng.uniform(.1,5,size=(3,4))
    got=[]
    for colors in itertools.product(range(4),repeat=3):
        zc=np.zeros((1,4,4),complex);nc=np.zeros((1,4,4))
        for i,c in enumerate(colors):zc[0,c]+=raw[i];nc[0,c]+=n[i]
        got.append(cubic(zc,nc)[0])
    c=0;w=0
    for i,j,k in TRIPLES:
        w+=np.sum(raw[:,i]*raw[:,j].conj()*n[:,k]).real/4
        c+=(np.sum(raw[:,i]*n[:,j]*raw[:,k].conj())+np.sum(n[:,i]*raw[:,j]*raw[:,k].conj())).real/8
    error=float(np.max(abs(np.mean(got,0)-[c,w])));assert error<1e-12
    return dict(color_assignments=64,max_error=error)

def main():
    configs=[dict(seed=2609145101+j,gamma=g,tau=t,hidden=hidden,quality_shift=qs,hidden_shift=hs)
        for j,(g,t,hidden,qs,hs) in enumerate([(0,0,40,0,0),(3,3,40,0,0),(3,3,40,2,0),(3,3,40,2,40)])]
    (P/'design.json').write_text(json.dumps(dict(configs=configs,groups=8192,sources=32,q=Q.tolist(),
        primary='Count-balanced cubic quality contrast cancels stationary hidden source noise.',
        limitation_test='The fourth world changes hidden source variance across windows.'),indent=2)+'\n')
    print('EXACT',exact(),flush=True);out=[];start=time.time()
    for j,cfg in enumerate(configs):
        a=generate(**cfg);np.savez_compressed(P/f'case_{j}.npz',**a)
        r=dict(**cfg,methods=summarize(a,cfg['tau']));out.append(r)
        (P/'results.json').write_text(json.dumps(json_safe(out),indent=2,allow_nan=False)+'\n')
        print(json.dumps(dict(**cfg,pooled={k:v['pooled'] for k,v in r['methods'].items()},rmse={k:v['rmse'] for k,v in r['methods'].items()})),flush=True)
    print('DONE',time.time()-start,flush=True)
if __name__=='__main__':main()
