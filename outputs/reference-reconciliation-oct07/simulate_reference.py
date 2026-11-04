"""Field-derived sensitivity of frozen clustering/QC, not physical truth."""
import argparse
import json
import sys
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from common import P, ROOT, save, sha
sys.path.insert(0,str(ROOT/'outputs/oracle-crossval-sep15'))
import oracle_pipeline as O
from reference_characterize import ANGLES, count_weight

def fast_clusters(pos,thresh=10.):
    """Same strict 2D single-linkage partition without a dense N-by-N array."""
    if len(pos)==0:return np.full(0,-1,int)
    pairs=cKDTree(pos[:,:2]).query_pairs(thresh,output_type='ndarray')
    if len(pairs):
        distance=np.linalg.norm(pos[pairs[:,0],:2]-pos[pairs[:,1],:2],axis=1)
        pairs=pairs[distance<thresh]
    graph=coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(len(pos),len(pos)))
    return connected_components(graph,directed=False,return_labels=True)[1]

def templates(tag):
    z=np.load(P/f'empirical-{tag}.npz')
    events=z['events'];offsets=z['offsets']; result=[];pooled=[]
    for lo,hi in zip(offsets[:-1],offsets[1:]):
        e=events[lo:hi];labels=e[:,6].astype(int);pan=e[:,3].astype(int)
        sources=[]
        for s in np.unique(labels):
            m=labels==s;xy=e[m,:2];pt=pan[m];res=xy.copy();df=0
            for t in np.unique(pt):
                k=pt==t;res[k]-=res[k].mean(0);df+=max(int(k.sum())-1,0)
            covariance=res.T@res/max(df,1)
            if df>=4:pooled.append(covariance)
            sources.append(dict(center=e[m,:3].mean(0),pan=pt,residual=res,cov=covariance,df=df,
                                var=e[m,4],multiplicity=e[m,5],z_residual=e[m,2]-e[m,2].mean()))
        result.append(sources)
    # Same global fallback for both noise models, without using oracle QC.
    global_cov=np.median(np.stack(pooled),axis=0)
    def positive(c):
        vals,vec=np.linalg.eigh(c);return (vec*np.maximum(vals,.01))@vec.T
    global_cov=positive(global_cov)
    for sources in result:
        for s in sources:
            s['cov']=positive(s['cov'] if s['df']>=4 else global_cov)
            s['chol']=np.linalg.cholesky(s['cov'])
    return result,json.loads(str(z['qc_json'])),dict(groups=len(result),global_fallback_cov=global_cov,
        fallback_sources=sum(s['df']<4 for g in result for s in g),
        inferred_sources=sum(len(g) for g in result), empirical_sha256=sha(P/f'empirical-{tag}.npz'))

def run_world(template,qc,T,model,seed,draws):
    rng=np.random.default_rng(seed);records=[];rejections={};impure=0;passes=0;candidate_events=0;passing_events=0
    by_group=[]
    for g in range(draws):
        sources=template[rng.integers(len(template))]
        xyz=[];pan=[];var=[];mult=[];labels=[];truth=[]
        for i,s in enumerate(sources):
            k=len(s['pan']);window_motion=rng.normal(0,np.sqrt(T),size=(2,2))
            if model=='gaussian':noise=rng.normal(size=(k,2))@s['chol'].T
            else:
                noise=np.empty((k,2))
                for t in np.unique(s['pan']):
                    m=s['pan']==t;pool=s['residual'][m];n=len(pool)
                    if n>=3:
                        noise[m]=pool[rng.integers(n,size=n)]*np.sqrt(n/(n-1))
                    else:noise[m]=rng.normal(size=(n,2))@s['chol'].T
            xy=s['center'][:2]+window_motion[s['pan']//2]+noise
            # Reuse z noise as a nuisance; it is not a newly validated 3D model.
            zz=s['center'][2]+s['z_residual'][rng.integers(k,size=k)]
            xyz.append(np.column_stack([xy,zz]));pan.extend(s['pan']);var.extend(s['var']);mult.extend(s['multiplicity'])
            labels.extend([i]*k);truth.append(window_motion)
        pos=np.concatenate(xyz);pan=np.asarray(pan,int);var=np.asarray(var);mult=np.asarray(mult);labels=np.asarray(labels)
        lab=fast_clusters(pos)
        ok,reasons=O.source_qc(pos,lab,pan,**qc)
        for name,v in reasons.items():rejections[name]=rejections.get(name,0)+v
        rows=O.oracle_stats(pos,var,lab,pan,ok,ANGLES)
        candidate_events+=len(pos);gn=np.zeros(3);gd=np.zeros(3)
        for row in rows:
            m=lab==row['source'];members=labels[m]
            unique,cnt=np.unique(members,return_counts=True);purity=cnt.max()/cnt.sum()
            impure+=len(unique)>1;passes+=1;passing_events+=int(m.sum())
            lw=np.bincount(pan[m],weights=mult[m],minlength=4)
            weights=np.array([1.,count_weight(row['counts']),count_weight(lw)])
            value=float(row['that'].mean());gn+=weights*value;gd+=weights
        by_group.append(np.r_[gn,gd])
    a=np.asarray(by_group);num=a[:,:3].sum(0);den=a[:,3:].sum(0)
    estimate=np.divide(num,den,out=np.full(3,np.nan),where=den>0)
    influence=a[:,:3]-a[:,3:]*estimate
    se=np.sqrt(draws)*influence.std(0,ddof=1)/den
    return dict(T=T,model=model,seed=seed,group_draws=draws,passing_sources=passes,
                estimate=estimate,se=se,rejections=rejections,
                impure_fraction=impure/passes if passes else None,
                event_coverage=passing_events/candidate_events,
                by_group=a)

def main(tag,draws):
    template,qc,metadata=templates(tag);out=[];t0=time.time()
    metadata.update(clustering='sparse strict 2D connected components; partition-equivalent to legacy',
                    simulation_source_sha256=sha(__file__))
    for model in ('gaussian','residual'):
        for T in (0.,1.,3.,5.):
            for seed in (2610071101,2610071102):
                record=run_world(template,qc,T,model,seed,draws)
                groups=record.pop('by_group')
                name=f'{tag}-{model}-T{int(T)}-s{seed}'
                np.savez_compressed(P/f'control-{name}.npz',by_group=groups)
                record['group_records']=f'control-{name}.npz';out.append(record)
                save(P/f'controls-{tag}.json',dict(tag=tag,metadata=metadata,qc=qc,records=out,
                    complete=len(out)==16,runtime_seconds=time.time()-t0))
                print(tag,model,T,seed,record['estimate'],record['passing_sources'],flush=True)
    print('complete',tag,'seconds',time.time()-t0,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tag',required=True);p.add_argument('--draws',type=int,default=400)
    a=p.parse_args();assert a.draws>=400;main(a.tag,a.draws)
