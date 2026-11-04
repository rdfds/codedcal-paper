"""Assignment-based oracle for persistent-displacement variance.

Clusters localizations to individual sources within each (pore, round),
links them into binding events, and computes per-source panel means. The
oracle contrast per source and coordinate is

    C   = (a1+a2)/2 - (b1+b2)/2      E[C^2] = 2T + (n1+n2+n3+n4)/4
    W_A = a1-a2, W_B = b1-b2         E[W_A^2] = n1+n2, E[W_B^2] = n3+n4
    That = [C^2 - (W_A^2+W_B^2)/4] / 2

with n_t the noise variance of the panel-t mean. That is unbiased for T
per source under within-window persistence, the same assumption the
count-balanced code estimator makes, but it uses source assignments the
code estimator never sees. Clustering quality is controlled and must be
validated on synthetic data before the real comparison is trusted.
"""
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
# Measured event-noise law for this field (drift-free, assignment-based,
# interleaved pair differences vs reported precision; see report):
# true event variance per axis = NOISE_A + NOISE_B * reported event variance.
NOISE_A=20.06;NOISE_B=1.71
PANELS=[(0,3115),(3135,6240),(6260,9365),(9385,12500)]

def link_events(frames,xy,lp,extra=None,gap=5,jump=25.):
    # localizations sorted by frame -> event index array.
    order=np.argsort(frames,kind='stable')
    f=frames[order];p=xy[order]
    new=np.ones(len(f),bool)
    if len(f)>1:
        d=np.linalg.norm(np.diff(p,axis=0),axis=1)
        new[1:]=(np.diff(f)>gap)|(d>jump)
    ev=np.cumsum(new)-1
    out=np.empty(len(f),int);out[order]=ev
    return out

def event_table(frames,xy,lp,ev):
    # precision-weighted event positions; event frame = mean frame.
    k=int(ev.max())+1;w=1/np.maximum(lp,1e-3)**2
    sw=np.bincount(ev,weights=w,minlength=k)
    pos=np.stack([np.bincount(ev,weights=w*xy[:,i],minlength=k)/sw for i in range(xy.shape[1])],-1)
    fr=np.bincount(ev,weights=frames.astype(float),minlength=k)/np.bincount(ev,minlength=k)
    n=np.bincount(ev,minlength=k)
    var=1/sw # nominal (reported-precision) variance of the event mean, per axis
    return pos,fr,n,var

def assign_panels(fr):
    pan=np.full(len(fr),-1)
    for t,(lo,hi) in enumerate(PANELS):pan[(fr>=lo)&(fr<hi)]=t
    return pan

def cluster_events(pos,thresh=10.,min_events=8):
    # Greedy agglomerative single-linkage on 2D event positions with a
    # distance cap; small n per pore-round makes O(n^2) fine.
    n=len(pos)
    if n==0:return np.full(0,-1)
    lab=np.arange(n)
    d=np.linalg.norm(pos[:,None,:2]-pos[None,:,:2],axis=-1)
    pairs=np.argwhere((d<thresh)&(np.arange(n)[:,None]<np.arange(n)[None,:]))
    order=np.argsort(d[pairs[:,0],pairs[:,1]]) if len(pairs) else []
    def find(i):
        while lab[i]!=i:lab[i]=lab[lab[i]];i=lab[i]
        return i
    for idx in order:
        i,j=pairs[idx];lab[find(i)]=find(j)
    out=np.array([find(i) for i in range(n)])
    _,out=np.unique(out,return_inverse=True)
    return out

def split_gain(r):
    # 2-means style merge detector on centered 2D positions: fraction of
    # within-cluster variance removed by the best split along the principal
    # axis. Merged ~12 nm pairs score high; single sources score low.
    if len(r)<6:return 0.
    u,s,_=np.linalg.svd(r-r.mean(0),full_matrices=False)
    proj=np.sort(u[:,0]*s[0]);tot=(proj**2).sum()
    if tot<=0:return 0.
    best=0.
    for k in range(2,len(proj)-1):
        a=proj[:k];b=proj[k:]
        rem=((a-a.mean())**2).sum()+((b-b.mean())**2).sum()
        best=max(best,1-rem/tot)
    return float(best)

def source_qc(pos,lab,pan,min_per_panel=2,min_events=12,max_rms=10.,min_sep=30.,
              max_split=.8,max_zsd=22.,neighbor_min_events=5):
    # Returns per-source dict rows that pass quality control, plus counts of
    # each rejection reason. Isolation is measured against neighboring
    # clusters with at least neighbor_min_events events (real competing
    # sources), not against single background events.
    # split_gain and z spread are always recorded; max_split filters if set.
    ok=[];reasons=dict(size=0,panel_coverage=0,rms=0,isolation=0,split=0,zsd=0)
    sizes=np.bincount(lab)
    ids=np.arange(len(sizes))
    centers=np.array([pos[lab==s,:2].mean(0) for s in ids])
    real=sizes>=neighbor_min_events
    for s in ids:
        m=lab==s
        if sizes[s]<min_events:reasons['size']+=1;continue
        counts=[int(((pan==t)&m).sum()) for t in range(4)]
        if min(counts)<min_per_panel:reasons['panel_coverage']+=1;continue
        # QC statistics must be invariant to the between-window displacement
        # being measured, or selection biases the oracle low: center per
        # window before computing rms and the merge detector.
        r=pos[m,:2].copy();win=pan[m]//2;nw=0
        for wv in (0,1):
            mw=win==wv
            if mw.any():r[mw]-=r[mw].mean(0);nw+=1
        dof=max(len(r)-nw,1)
        rms=float(np.sqrt((r**2).sum()/ (2*dof)))*np.sqrt(2)
        if rms>max_rms:reasons['rms']+=1;continue
        comp=real.copy();comp[s]=False
        sep=float(np.linalg.norm(centers[comp]-centers[s],axis=1).min()) if comp.any() else np.inf
        if sep<min_sep:reasons['isolation']+=1;continue
        gain=split_gain(r)
        if max_split is not None and gain>max_split:reasons['split']+=1;continue
        zs=float(pos[m,2].std(ddof=1)) if pos.shape[1]>2 else None
        if max_zsd is not None and zs is not None and zs>max_zsd:
            reasons['zsd']+=1;continue
        ok.append(dict(source=int(s),counts=counts,center=centers[s],sep=sep,
            rms=rms,split_gain=gain,z_sd=zs))
    return ok,reasons

def oracle_stats(pos,var,lab,pan,ok,angles):
    # Per passing source: precision-weighted panel means on projected
    # coordinates, contrasts, and count weights. Weighting by reported
    # per-event variance makes high-noise tail events nearly irrelevant, so
    # linkage-radius truncation cannot pull window means toward the pooled
    # center (a downward bias measured on synthetic data with unweighted
    # means). The W-based noise correction is valid for any fixed weights.
    # NOTE: no event-selection refinement is applied inside a cluster. A
    # per-window median cut was tried and rejected: selecting events around
    # a center shared by the two panels of a window correlates their mean
    # errors, inflating C^2 and deflating W^2 (+1.2 nm^2 at true T=0 on
    # synthetic). Static contaminant dilution is instead calibrated on
    # synthetic fields with the measured junk rate.
    # Weights follow the measured noise law (nearly uniform, because the
    # additive linkage floor dominates); the W-based correction is valid for
    # any fixed weights, so the law's exact constants affect efficiency only.
    w_all=(1/(NOISE_A+NOISE_B*np.clip(var,.25,None)) if var is not None
           else np.ones(len(pos)))
    w_rep=1/np.clip(var,.25,None) if var is not None else np.ones(len(pos))
    rows=[]
    for rec in ok:
        m=lab==rec['source'];pm=[];nm=[];rep=[];repm=[]
        for t in range(4):
            mt=m&(pan==t)
            proj=pos[mt,0][:,None]*np.cos(angles)+pos[mt,1][:,None]*np.sin(angles)
            w=w_all[mt];sw=w.sum()
            mean=(w[:,None]*proj).sum(0)/sw
            pm.append(mean)
            # empirical noise of the weighted panel mean from event scatter
            nm.append(((w[:,None]*(proj-mean)**2).sum(0)/sw)*(w**2).sum()/sw**2
                      *mt.sum()/max(mt.sum()-1,1))
            rep.append(1/w_rep[mt].sum()) # reported-precision variance
            repm.append((w**2).sum()/sw**2*0+1/sw) # model-law variance
        pm=np.array(pm);nm=np.array(nm);rep=np.array(rep);repm=np.array(repm)
        C=(pm[0]+pm[1])/2-(pm[2]+pm[3])/2
        WA=pm[0]-pm[1];WB=pm[2]-pm[3]
        that=(C**2-(WA**2+WB**2)/4)/2
        rows.append(dict(**rec,C=C,WA=WA,WB=WB,that=that,
            noise=nm.sum(0)/4,reported_var=rep,model_var=repm,panel_means=pm))
    return rows
