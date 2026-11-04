"""Replay frozen references and retain field-derived ingredients before QC."""
import argparse
import sys
import time
import h5py
import numpy as np
from common import CACHE, P, ROOT, cd20_registry, read, save, sha, input_path

sys.path.insert(0, str(ROOT / 'outputs/oracle-crossval-sep15'))
sys.path.insert(0, str(ROOT / 'outputs/stationarity-diagnostic-sep15'))
sys.path.insert(0, str(ROOT / 'outputs/replication-cd20-sep15'))
import oracle_pipeline as O
import cd20_common as E
from split_components import Q, components, split_summary
ANGLES = np.arange(4) * np.pi / 4
TRIPLES = (((0, 1, 2), (0, 1, 3)), ((2, 3, 0), (2, 3, 1)))

def assign_panels(frame, panels):
    result = np.full(len(frame), -1, int)
    for t, (lo, hi) in enumerate(panels):
        result[(frame >= lo) & (frame < hi)] = t
    return result

def count_weight(counts):
    n = np.asarray(counts, float)
    return sum(n[..., i] * n[..., j] * n[..., k] for win in TRIPLES for i, j, k in win)

def phase_sums(xy, pan):
    z = np.zeros((4, len(Q), len(ANGLES)), complex)
    n = np.bincount(pan[pan >= 0], minlength=4).astype(float)
    proj = xy[:, 0, None] * np.cos(ANGLES) + xy[:, 1, None] * np.sin(ANGLES)
    for t in range(4):
        z[t] = np.exp(1j * proj[pan == t, None, :] * Q[None, :, None]).sum(0)
    return z, n

def diagonal(z, n):
    """Source-diagonal terms with the same localizations/counts as cubic."""
    out = np.zeros((len(Q), 2, 4))
    for w, triples in enumerate(TRIPLES):
        for i, j, k in triples:
            within = (z[i] * z[j].conj()).real * n[k] / 2
            o1 = (z[i] * z[k].conj()).real * n[j] / 2
            o2 = (z[j] * z[k].conj()).real * n[i] / 2
            out[:, w] += np.stack([(o1 + o2) / 2, within, o1, o2], -1).mean(1)
    return out

def boot_means(that, weights, pore, draws=1000):
    ids, inv = np.unique(pore, return_inverse=True)
    numerator = np.bincount(inv, weights=that * weights)
    denominator = np.bincount(inv, weights=weights)
    rng = np.random.default_rng(26100741)
    indices = rng.integers(len(ids), size=(draws, len(ids)))
    return np.quantile(numerator[indices].sum(1) / denominator[indices].sum(1), [.025, .975])

def registry():
    rows = cd20_registry()
    npc = dict(tag='npc', prefix='Figure2,EDF06,EDF07/',
               center=str(ROOT/'outputs/npc-random-code-validation-sep13/groups.npz'),
               oracle=str(ROOT/'outputs/oracle-crossval-sep15/real_oracle_sources.npz'),
               panels={str(r): O.PANELS for r in range(1, 5)},
               frames={str(r): 12500 for r in range(1, 5)},
               old=read('outputs/oracle-crossval-sep15/real_oracle.json')['points']['A'])
    rows.insert(0, npc)
    f = read('outputs/independent-field-sep15/frozen_constants.json')
    panels = [(0,4200),(4250,8450),(8500,12700),(12750,17000)]
    sdc = dict(tag='sdc_corrected', prefix=None,
               center=str(ROOT/'outputs/independent-field-sep15/npc_centers.npz'),
               oracle=str(ROOT/'outputs/independent-field-sep15/oracle_sources.npz'),
               panels={str(r): panels for r in range(1,5)},
               frames={str(r): 17000 for r in range(1,5)},
               frozen=f, old=read('outputs/independent-field-sep15/field_results.json'))
    rows.append(sdc)
    rows.append(dict(sdc, tag='sdc_raw', oracle=None))
    return rows

def load_round(row, r, centers):
    path=input_path(row,r)
    if row['tag'].startswith('sdc_'):
        prior=read('outputs/heldout-nup96-sep30/sdc-recovery.json')[f'round{r}']
        assert sha(path)==prior['reduced_sha256'],path
        d=np.load(path)
        lpx, lpy=d['lpx'],d['lpy']
        lp=np.sqrt((lpx**2+lpy**2)/2)
        k=np.isfinite(lp)&(lp<=row['frozen']['lp_max_nm'])
        xy=np.column_stack([d['x'][k],d['y'][k]])-centers[d['pick'][k]]
        if row['tag']=='sdc_corrected':
            traj=np.load(ROOT/f'outputs/independent-field-sep15/drift_trace_round{r}.npy')
            tt=(np.arange(len(traj))+.5)*row['frozen']['drift_bin_frames']
            for c in range(2): xy[:,c]-=np.interp(d['frame'][k],tt,traj[:,c])
        return dict(frame=d['frame'][k],xy=xy,already_centered=True,lp=lp[k],lpx=lpx[k],lpy=lpy[k],
                    pick=d['pick'][k],z=None,path=str(path),sha256=sha(path),
                    source_sha256=prior['source_sha256'])
    with h5py.File(path) as f: a=f['locs'][:]
    lp=130*np.sqrt((a['lpx']**2+a['lpy']**2)/2)
    k=np.isfinite(lp)&(lp<=8);a=a[k];lp=lp[k]
    xy=np.column_stack([a['x'],a['y']])*130
    pick=a['group'].astype(int) if row['tag']=='npc' else E.assign(xy,centers[:,:2])
    m=pick>=0
    return dict(frame=a['frame'][m].astype(int),xy=xy[m],lp=lp[m],
                lpx=a['lpx'][m]*130,lpy=a['lpy'][m]*130,pick=pick[m],
                z=a['z'][m].astype(float) if 'z' in a.dtype.names else None,
                path=str(path),sha256=sha(path))

def characterize(row):
    tag=row['tag']; t0=time.time()
    g=np.load(row['center']); centers=g['centers'];ng=len(centers)
    oracle_rows=[]; ingredients=[]; group_index=[]; rejection={}; support=[]; hashes=[]
    selected_z=np.zeros((ng,4,4,len(Q),len(ANGLES)),complex)
    selected_n=np.zeros((ng,4,4))
    full_z=np.zeros_like(selected_z);full_n=np.zeros_like(selected_n)
    diag=np.zeros((ng,len(Q),2,4));loc_weights=[];event_weights=[]
    all_local=0;selected_local=0;all_events=0;selected_events=0
    for r in range(1,5):
        d=load_round(row,r,centers);panels=row['panels'][str(r)]
        hashes.append(dict(round=r,path=d['path'],sha256=d['sha256'],source_sha256=d.get('source_sha256')))
        assert int(d['frame'].max())<int(row['frames'][str(r)])
        order=np.argsort(d['pick'],kind='stable');bounds=np.searchsorted(d['pick'][order],np.arange(ng+1))
        group_before=len(ingredients)
        for p in range(ng):
            ii=order[bounds[p]:bounds[p+1]]
            xy=d['xy'][ii] if d.get('already_centered') else d['xy'][ii]-centers[p,:2]
            frame=d['frame'][ii];lp=d['lp'][ii]
            lp_pan=assign_panels(frame,panels)
            z,n=phase_sums(xy,lp_pan);full_z[p,r-1]+=z;full_n[p,r-1]+=n
            all_local+=int(n.sum())
            if len(ii)<20:continue
            ev=O.link_events(frame,xy,lp,jump=row.get('frozen',{}).get('link_jump_nm',25.))
            zcol=d['z'][ii] if d['z'] is not None else np.zeros(len(ii))
            pos,evf,nloc,var=O.event_table(frame,np.column_stack([xy,zcol]),lp,ev)
            epan=assign_panels(evf,panels);keep=epan>=0;ids=np.flatnonzero(keep)
            pos=pos[keep];epan=epan[keep];var=var[keep];evf=evf[keep];nloc=nloc[keep]
            if len(pos)<8:continue
            lab=O.cluster_events(pos,thresh=row.get('frozen',{}).get('cluster_thresh_nm',10.))
            qc=dict(row.get('frozen',{}).get('oracle_qc',{}))
            if d['z'] is None:qc['max_zsd']=None
            ok,reasons=O.source_qc(pos,lab,epan,**qc)
            for k,v in reasons.items():rejection[k]=rejection.get(k,0)+v
            # Retain ALL candidate event clusters, before source eligibility.
            # [x,y,z,panel,reported event variance,localization multiplicity,label,frame]
            table=np.column_stack([pos,epan,var,nloc,lab,evf])
            ingredients.append(table);group_index.append([r-1,p,int(d['z'] is not None)])
            all_events+=len(pos)
            recs=O.oracle_stats(pos,var,lab,epan,ok,ANGLES)
            for rec in recs:
                rec.update(pore=p,round=r-1);oracle_rows.append(rec)
                eids=ids[lab==rec['source']]
                mask=np.isin(ev,eids)&(lp_pan>=0)
                zs,ns=phase_sums(xy[mask],lp_pan[mask])
                selected_z[p,r-1]+=zs;selected_n[p,r-1]+=ns
                diag[p]+=diagonal(zs,ns)
                loc_weights.append(count_weight(ns));event_weights.append(count_weight(rec['counts']))
                selected_local+=int(ns.sum());selected_events+=int((lab==rec['source']).sum())
                support.append([r-1,p,*rec['counts'],*ns,
                                *[float(evf[(lab==rec['source'])&(epan==t)].mean()) for t in range(4)]])
        print(tag,'round',r,'groups',len(ingredients)-group_before,'sources',len(oracle_rows),flush=True)
    that=np.array([rec['that'].mean() for rec in oracle_rows]);pore=np.array([rec['pore'] for rec in oracle_rows])
    rnd=np.array([rec['round'] for rec in oracle_rows]);event_weights=np.array(event_weights);loc_weights=np.array(loc_weights)
    oldnp=np.load(row['oracle']) if row['oracle'] else None
    replay=None
    if oldnp is not None:
        replay=dict(sources_equal=len(that)==len(oldnp['that']),
                    pore_equal=np.array_equal(pore,oldnp['pore']),round_equal=np.array_equal(rnd,oldnp['round']))
        if len(that)==len(oldnp['that']):
            replay.update(max_abs_that=float(np.max(abs(that-oldnp['that']))),
                          event_count_weights_equal=np.array_equal(event_weights,oldnp['countw']))
        replay['passed']=bool(replay['sources_equal'] and replay['pore_equal'] and replay['round_equal']
                              and replay.get('max_abs_that',np.inf)<1e-9 and replay.get('event_count_weights_equal'))
    def code_comp(z,n):
        return sum(np.stack([components(z[:,:,:,j,k],n) for j in range(len(Q))],1)
                   for k in range(len(ANGLES)))/len(ANGLES)
    selected=code_comp(selected_z,selected_n);full=code_comp(full_z,full_n)
    weights=dict(unweighted=np.ones(len(that)),event_count=event_weights,localization_count=loc_weights)
    estimates={k:dict(raw=float(np.average(that,weights=w)),interval=boot_means(that,w,pore),
                      effective_sources=float(w.sum()**2/(w**2).sum())) for k,w in weights.items()}
    offsets=np.r_[0,np.cumsum([len(a) for a in ingredients])]
    np.savez_compressed(P/f'empirical-{tag}.npz',events=np.concatenate(ingredients),
                        offsets=offsets,group_index=group_index,has_z=d['z'] is not None,
                        qc_json=np.array(__import__('json').dumps(qc)))
    np.savez_compressed(P/f'reference-{tag}.npz',that=that,pore=pore,round=rnd,
                        event_countw=event_weights,localization_countw=loc_weights,support=np.array(support),
                        diagonal=diag,selected_code=selected,full_code=full,centers=centers)
    # Spatial-block bootstrap / Fieller checks are deliberately a separate analysis.
    summary=dict(tag=tag,distinct_acquisition=tag!='sdc_raw',panels=row['panels'],
                 frames=row['frames'],input_members=hashes,reference_replay=replay,
                 passing_sources=len(that),candidate_event_groups=len(ingredients),rejections=rejection,
                 estimates=estimates,event_coverage=selected_events/all_events,
                 localization_coverage=selected_local/all_local,
                 weight_rank_correlation=float(__import__('scipy').stats.spearmanr(event_weights,loc_weights).statistic),
                 per_round={str(r+1):dict(sources=int((rnd==r).sum()),
                    raw=float(np.mean(that[rnd==r])) if (rnd==r).any() else None) for r in range(4)},
                 full_code=split_summary(full),selected_code=split_summary(selected),
                 assigned_diagonal=split_summary(diag),runtime_seconds=time.time()-t0)
    save(P/f'characterization-{tag}.json',summary)
    print('finished',tag,'replay',replay,'estimates',estimates,flush=True)
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--tag',required=True);args=parser.parse_args()
    row=next(row for row in registry() if row['tag']==args.tag)
    characterize(row)
