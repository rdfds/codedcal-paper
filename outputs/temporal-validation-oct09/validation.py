"""Frozen finite-panel validation. Only aggregate observations enter estimate()."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import os
import sys
import time
import numpy as np
from scipy.stats import t as student_t

P = Path(__file__).resolve().parent
ROOT = P.parent.parent
REF = ROOT / 'outputs/reference-reconciliation-oct07'
sys.path.insert(0, str(REF))
from common import save, sha
from simulate_reference import fast_clusters
from reference_characterize import registry, Q, ANGLES, TRIPLES
from split_components import components
import oracle_pipeline as O

TAGS = ('npc', 'edf10rep1', 'edf10rep2', 'fig4untreated', 'rtx_rep1',
        'rtx_rep2', 'rtx_rep3', 'sdc_corrected', 'sdc_raw')
MODELS = ('gaussian', 'residual')
PRIMARY = ('zero', 'homogeneous', 'group_common', 'heterogeneous')
STRESS = ('noise_shift', 'source_noise_shift', 'code_offset')
WORLDS = PRIMARY + STRESS
U = np.column_stack([np.cos(ANGLES), np.sin(ANGLES)])
SETTINGS = dict(primary_fields=40, stress_fields=12, codes_per_field=2,
                standardized_groups=512, bootstrap_draws=2000,
                cluster_fields=8, cluster_groups=16,
                variance_floor=.01, fallback_min_df=4, T=3., added_noise=25.,
                code_offsets=[-15., -5., 5., 15.],
                heterogeneity='f_s=(0.5+trace(cov_s)/(2*median_axis_var)); normalized by arithmetic mean over ALL retained sources',
                source_noise_shift='additional window-B covariance=25*f_s*I',
                bias_limit=.01, usable_limit=.90, coverage_limit=.90,
                alpha=.05, provenance='conditional simulation from all retained pre-QC Source Data 69 candidate clusters')


def seed(*parts):
    return int.from_bytes(hashlib.sha256('|'.join(map(str, parts)).encode()).digest()[:8], 'little')


def positive(cov):
    vals, vecs = np.linalg.eigh(cov)
    return (vecs * np.maximum(vals, SETTINGS['variance_floor'])[..., None, :]) @ np.swapaxes(vecs, -1, -2)


def bincount_complex(index, weights, size):
    out = np.empty((size, *weights.shape[1:]), complex)
    flat = weights.reshape(len(index), -1)
    target = out.reshape(size, -1)
    for j in range(flat.shape[1]):
        target[:, j] = np.bincount(index, weights=flat[:, j].real, minlength=size)
        target[:, j] += 1j * np.bincount(index, weights=flat[:, j].imag, minlength=size)
    return out


def templates(tag):
    """Synthetic identities are unique within each original round/group."""
    d = np.load(REF / f'empirical-{tag}.npz')
    events, offsets, index = d['events'], d['offsets'], d['group_index']
    ng = len(np.load(REF / f'reference-{tag}.npz')['centers'])
    by_pore = [[] for _ in range(ng)]
    for a, (lo, hi) in enumerate(zip(offsets[:-1], offsets[1:])):
        e = events[lo:hi]
        by_pore[int(index[a, 1])].append((int(index[a, 0]), e))
    groups = []
    for pore, chunks in enumerate(by_pore):
        if not chunks:
            groups.append(None)
            continue
        parts, origins, offset = [], [], 0
        for rnd, original in chunks:
            e = original.copy()
            labs, inv = np.unique(e[:, 6].astype(int), return_inverse=True)
            e[:, 6] = inv + offset
            origins.extend((rnd, int(lab)) for lab in labs)
            offset += len(labs)
            parts.append(e)
        e = np.concatenate(parts)
        e = e[np.argsort(e[:, 6], kind='stable')]
        s = e[:, 6].astype(int); ns = offset; n = np.bincount(s, minlength=ns)
        pan = e[:, 3].astype(int); sp = 4*s + pan
        ec = np.bincount(sp, minlength=ns*4).reshape(ns, 4)
        lc = np.bincount(sp, weights=e[:, 5], minlength=ns*4).reshape(ns, 4)
        means = np.stack([np.bincount(sp, weights=e[:, j], minlength=ns*4) for j in range(2)], -1)
        means /= np.maximum(ec.reshape(-1, 1), 1)
        res = e[:, :2] - means[sp]
        df = n - (ec > 0).sum(1)
        cov = np.stack([np.bincount(s, weights=res[:, a]*res[:, b], minlength=ns)
                        for a, b in ((0, 0), (0, 1), (1, 0), (1, 1))], -1).reshape(ns, 2, 2)
        cov /= np.maximum(df, 1)[:, None, None]
        center = np.stack([np.bincount(s, weights=e[:, j], minlength=ns)/n for j in range(3)], -1)
        starts = np.r_[0, np.cumsum(n)[:-1]]
        groups.append(dict(source=s, pan=pan, ec=ec, lc=lc, n=n, df=df, cov=cov,
                           center=center, res=res, starts=starts, var=e[:, 4], mult=e[:, 5],
                           zres=e[:, 2]-center[s, 2], frame=e[:, 7], origins=np.array(origins), pore=pore))
    supported = np.concatenate([g['cov'][g['df'] >= 4] for g in groups if g is not None])
    assert len(supported), 'No fallback covariance can be estimated'
    fallback = positive(np.median(supported, axis=0))
    for g in groups:
        if g is None: continue
        g['cov'] = positive(np.where((g['df'] >= 4)[:, None, None], g['cov'], fallback))
        g['chol'] = np.linalg.cholesky(g['cov'])
    allvar = np.concatenate([np.trace(g['cov'], axis1=1, axis2=2)/2 for g in groups if g is not None])
    median = float(np.median(allvar))
    normalizer = float((.5+allvar/median).mean())
    for g in groups:
        if g is not None:
            g['factor'] = (.5+np.trace(g['cov'], axis1=1, axis2=2)/(2*median))/normalizer
    row = next(r for r in registry() if r['tag'] == tag)
    qc = json.loads(str(d['qc_json']))
    meta = dict(tag=tag, groups=ng, nonempty_groups=sum(g is not None for g in groups),
                inferred_sources=len(allvar), supported_covariance_sources=len(supported),
                fallback_sources=sum(int((g['df'] < 4).sum()) for g in groups if g is not None),
                events=len(events), localizations=int(events[:, 5].sum()),
                reported_event_variance_quantiles=np.quantile(events[:, 4], [0,.25,.5,.75,1]),
                multiplicity_quantiles=np.quantile(events[:, 5], [0,.25,.5,.75,1]),
                median_axis_variance=median, factor_normalizer=normalizer, fallback_cov=fallback,
                source_triple_supported=sum(int((g['ec'] > 0).sum(1).clip(0,4).astype(int).__ge__(3).sum()) for g in groups if g is not None),
                empirical_sha256=sha(REF/f'empirical-{tag}.npz'),
                reference_sha256=sha(REF/f'reference-{tag}.npz'),
                historical_characterization_sha256=sha(REF/f'characterization-{tag}.json'),
                panels=row['panels'], frames=row['frames'], qc=qc,
                cluster_radius=row.get('frozen', {}).get('cluster_thresh_nm', 10.))
    cat = {k: np.concatenate([g[k] for g in groups if g is not None]) for k in ('center','cov','df','ec','lc','n','factor','origins')}
    cat['pore'] = np.concatenate([np.full(len(g['n']), g['pore']) for g in groups if g is not None])
    # Catalog sufficient to audit the conditioned population; event residuals remain in hashed inputs.
    return groups, meta, cat


def kernel(g, model):
    v = np.einsum('ai,sij,aj->sa', U, g['cov'], U)
    gaussian = np.exp(-.5 * Q[None, :, None]**2 * v[:, None, :])
    if model == 'gaussian': return gaussian.astype(complex)
    res = g['res'] * np.sqrt(g['n'][g['source']]/np.maximum(g['df'][g['source']], 1))[:, None]
    phases = np.exp(1j * (res @ U.T)[:, None, :] * Q[None, :, None])
    out = bincount_complex(g['source'], phases, len(g['n']))/g['n'][:, None, None]
    return np.where((g['df'] >= 4)[:, None, None], out, gaussian)


def diagonal(z, n):
    """Known-source comparator. z: source,panel,q,angle; no code averaging."""
    out = np.zeros((len(Q), 2, 2))
    for w, triples in enumerate(TRIPLES):
        for i, j, k in triples:
            within = (z[:, i]*z[:, j].conj()).real * n[:, k, None, None]/2
            cross = ((z[:, i]*z[:, k].conj()).real*n[:, j, None, None]
                     + (z[:, j]*z[:, k].conj()).real*n[:, i, None, None])/4
            out[:, w, 0] += cross.sum(0).mean(-1)
            out[:, w, 1] += within.sum(0).mean(-1)
    return out


def truth_group(g, model, world):
    if g is None: return np.zeros((2, len(Q), 2, 2)), np.empty(0)
    if model not in g.setdefault('_kernels', {}): g['_kernels'][model] = kernel(g, model)
    phi = g['_kernels'][model]
    T = np.zeros(len(g['n'])) if world == 'zero' else np.full(len(g['n']), 3.)
    if world == 'heterogeneous': T *= g['factor']
    atten = np.exp(-T[:, None]*Q[None, :]**2)
    baseline = np.zeros((len(g['n']), len(Q), 2, 2))
    observable = np.zeros_like(baseline)
    panel_phi = np.repeat(phi[:, None], 4, axis=1)
    if world in ('noise_shift', 'source_noise_shift'):
        extra = 25.*(g['factor'] if world == 'source_noise_shift' else np.ones(len(T)))
        panel_phi[:, 2:] *= np.exp(-.5*extra[:, None, None]*Q[None, :, None]**2)[:, None]
    for w, triples in enumerate(TRIPLES):
        for i,j,k in triples:
            wt = g['lc'][:, i]*g['lc'][:, j]*g['lc'][:, k]
            base_w = (phi*phi.conj()).real.mean(-1)*wt[:, None]/2
            baseline[:, :, w, 1] += base_w
            baseline[:, :, w, 0] += base_w*atten
            a,b,c = panel_phi[:,i], panel_phi[:,j], panel_phi[:,k]
            observable[:, :, w, 1] += (a*b.conj()).real.mean(-1)*wt[:,None]/2
            observable[:, :, w, 0] += ((a*c.conj()).real+(b*c.conj()).real).mean(-1)*wt[:,None]*atten/4
    # Weight for sampled pure-source reference coverage: sum both windows and q=.10.
    mass = baseline[:, 0, :, 1].sum(-1)
    return np.stack([baseline.sum(0), observable.sum(0)]), mass


def aggregate_observations(z, n, code):
    """Source labels used ONLY to generate observed code aggregate arrays."""
    observed_z = np.zeros((4, 4, len(Q), len(ANGLES)), complex)
    observed_n = np.zeros((4, 4))
    np.add.at(observed_z, code, z); np.add.at(observed_n, code, n)
    return observed_z, observed_n


def estimate(observed_z, observed_n):
    """Inputs: group,code,panel,q,angle and group,code,panel; NO source labels."""
    return sum(np.stack([components(observed_z[:, :, :, j, a], observed_n)[..., :2]
                         for j in range(len(Q))], 1) for a in range(len(ANGLES)))/len(ANGLES)


def simulate_group(g, model, world, rng):
    if g is None: return None
    s, pan = g['source'], g['pan']; ns = len(g['n']); N = len(s)
    T = np.zeros(ns) if world == 'zero' else np.full(ns, 3.)
    if world == 'heterogeneous': T *= g['factor']
    if world == 'group_common':
        motion = rng.normal(size=(ns, 2, 2))*np.sqrt(1.5)
        motion += rng.normal(size=(1, 2, 2))*np.sqrt(1.5)
    else: motion = rng.normal(size=(ns, 2, 2))*np.sqrt(T)[:, None, None]
    normal = rng.normal(size=(N, 2))
    noise = np.einsum('eij,ej->ei', g['chol'][s], normal)
    if model == 'residual':
        draw = g['starts'][s] + (rng.random(N)*g['n'][s]).astype(int)
        residual = g['res'][draw]*np.sqrt(g['n'][s]/np.maximum(g['df'][s], 1))[:, None]
        noise = np.where((g['df'][s] >= 4)[:, None], residual, noise)
    if world in ('noise_shift', 'source_noise_shift'):
        extra = 25.*(g['factor'][s] if world == 'source_noise_shift' else np.ones(N))
        noise += rng.normal(size=(N, 2))*np.sqrt(extra)[:, None]*(pan >= 2)[:, None]
    xy = g['center'][s, :2]+motion[s, pan//2]+noise
    phases = np.exp(1j*(xy@U.T)[:, None, :]*Q[None,:,None])*g['mult'][:,None,None]
    z = bincount_complex(s*4+pan, phases, ns*4).reshape(ns,4,len(Q),len(ANGLES))
    return dict(z=z, n=g['lc'], xy=xy, pan=pan, source=s, g=g)


def clustering(ob, code, rng, qc, radius, mass):
    g = ob['g']; pos = np.column_stack([ob['xy'], g['center'][g['source'],2]+g['zres']])
    totals = dict(candidate_clusters=0, passing_clusters=0, impure_passing=0,
                  events=len(pos), localizations=float(g['mult'].sum()), passing_localizations=0.,
                  pure_kernel_mass=0., total_kernel_mass=float(mass.sum()), rejections={})
    pure_ids = set()
    for c in range(4):
        ids = np.flatnonzero(code[g['source']] == c)
        if not len(ids): continue
        lab = fast_clusters(pos[ids], radius)
        ok, reasons = O.source_qc(pos[ids], lab, g['pan'][ids], **qc)
        totals['candidate_clusters'] += int(lab.max()+1); totals['passing_clusters'] += len(ok)
        for key,v in reasons.items(): totals['rejections'][key] = totals['rejections'].get(key,0)+v
        for row in ok:
            m = ids[lab == row['source']]; unique = np.unique(g['source'][m])
            totals['passing_localizations'] += float(g['mult'][m].sum())
            totals['impure_passing'] += int(len(unique)>1)
            if len(unique)==1:
                # Fragmentation is possible: source weight credited only by captured panel triples.
                sid = int(unique[0]); counts = np.bincount(g['pan'][m], weights=g['mult'][m], minlength=4)
                captured = sum(counts[i]*counts[j]*counts[k] for win in TRIPLES for i,j,k in win)
                full = sum(g['lc'][sid,i]*g['lc'][sid,j]*g['lc'][sid,k] for win in TRIPLES for i,j,k in win)
                totals['pure_kernel_mass'] += float(mass[sid]*captured/full) if full>0 else 0.
    assert totals['pure_kernel_mass'] <= totals['total_kernel_mass']*(1+1e-10)
    return totals


def fieller(x, y, target, alpha=.05):
    n = len(x); mx,my = np.mean(x),np.mean(y)
    cov = np.cov(np.stack([x,y]), ddof=1)/n
    critical = student_t.ppf(1-alpha/2, n-1)**2
    A = my*my-critical*cov[1,1]
    B = -2*mx*my+2*critical*cov[0,1]
    C = mx*mx-critical*cov[0,0]
    disc = B*B-4*A*C
    lo,hi = np.nan,np.nan
    if A>0:
        kind = 'bounded' if disc>=0 else 'empty'
        if disc>=0: lo,hi = sorted(((-B-np.sqrt(disc))/(2*A),(-B+np.sqrt(disc))/(2*A)))
    elif A<0:
        kind = 'disjoint' if disc>0 else 'all_real'
        if disc>=0: lo,hi = sorted(((-B-np.sqrt(disc))/(2*A),(-B+np.sqrt(disc))/(2*A)))
    else:
        kind = 'half_line' if B!=0 else ('all_real' if C<=0 else 'empty')
    value = A*target**2+B*target+C
    covers = bool(value <= 1e-12*max(abs(A*target**2),abs(B*target),abs(C),np.finfo(float).tiny))
    usable = bool(kind == 'bounded' and lo>0 and my>0)
    return dict(kind=kind, lower=lo, upper=hi, A=A, B=B, C=C,
                covers=covers, usable=usable, ratio=mx/my if my>0 else np.nan,
                positive_reference=bool(my>0), reference=my, numerator=mx)


def diagnostic(comp):
    m = comp.mean(0); ra,rb = m[0,0]/m[0,1],m[1,0]/m[1,1]
    if not (m[0,1]>0 and m[1,1]>0 and ra>0 and rb>0):
        return dict(valid=False, quiet=False, z=np.nan)
    influence = (comp[:,0,0]-ra*comp[:,0,1])/(m[0,1]*ra)-(comp[:,1,0]-rb*comp[:,1,1])/(m[1,1]*rb)
    se = influence.std(ddof=1)/np.sqrt(len(comp))
    z = np.log(ra/rb)/se if se>0 else (0. if ra==rb else np.inf)
    return dict(valid=bool(np.isfinite(z)), quiet=bool(abs(z)<=student_t.ppf(.975,len(comp)-1)),z=z)


def field_records(coded, diag, targets, groups, rep):
    records = []
    for scope,n in (('observed',groups),('standardized',512)):
        for name, data in [('known_source',diag[None]),('coded',coded)]:
            for assignment, a in enumerate(data):
                for j,q in enumerate(Q):
                    pooled = a[:n,j].mean(1)
                    physical = targets[0,j,:,0].sum()/targets[0,j,:,1].sum()
                    observable = targets[1,j,:,0].sum()/targets[1,j,:,1].sum()
                    f = fieller(pooled[:,0],pooled[:,1],physical)
                    obs_cover = f['A']*observable**2+f['B']*observable+f['C']<=0
                    d = diagnostic(a[:n,j])
                    records.append(dict(rep=rep, scope=scope, group_count=n, method=name,
                                        assignment=assignment, q=float(q), physical_target=float(physical),
                                        observable_target=float(observable), observable_covers=bool(obs_cover),
                                        D=d, **f))
    return records


def atomic_npz(path, **arrays):
    temp = Path(str(path)+'.tmp')
    with temp.open('wb') as f: np.savez_compressed(f, **arrays)
    os.replace(temp,path)


def freeze():
    manifest = dict(settings=SETTINGS, tags=TAGS, models=MODELS, worlds=WORLDS,
                    protocol_sha256=sha(P/'PROTOCOL.md'), code_sha256=sha(__file__),
                    legacy_components_sha256=sha(ROOT/'outputs/stationarity-diagnostic-sep15/split_components.py'),
                    expected_worlds=126, physical_fields=3528, coded_fields=7056,
                    estimator_inputs=['observed_z[group,code,panel,q,angle]', 'observed_n[group,code,panel]'])
    save(P/'freeze.json',manifest)
    return manifest


def run(tag, model, world, prepared=None):
    name=f'{tag}-{model}-{world}'; folder=P/'worlds'/name; folder.mkdir(parents=True,exist_ok=True)
    groups,meta,cat=prepared if prepared is not None else templates(tag); ng=meta['groups']; N=max(ng,512)
    truth=[];mass=[]
    for g in groups:
        a,b=truth_group(g,model,world);truth.append(a);mass.append(b)
    truth=np.stack(truth); target=truth.mean(0)
    if not (P/f'catalog-{tag}.npz').exists():
        atomic_npz(P/f'catalog-{tag}.npz',**cat); save(P/f'catalog-{tag}.json',meta)
    truthfile=folder/'truth.npz'
    if not truthfile.exists(): atomic_npz(truthfile, by_template=truth, population=target)
    else: assert np.array_equal(np.load(truthfile)['population'],target)
    expected=40 if world in PRIMARY else 12
    statepath=folder/'state.json'
    state=json.loads(statepath.read_text()) if statepath.exists() else dict(tag=tag,model=model,world=world,completed=0,expected=expected,replicates=[],freeze_sha256=sha(P/'freeze.json'),complete=False)
    assert state['freeze_sha256']==sha(P/'freeze.json'), 'Code/settings changed: preserve previous results and repair explicitly'
    for rep in range(state['completed'],expected):
        started=time.time(); rng=np.random.default_rng(seed(name,rep,'physical'))
        code_rngs=[np.random.default_rng(seed(name,rep,'codes',a)) for a in range(2)]
        chosen=rng.integers(ng,size=N)
        observed_z=np.zeros((2,N,4,4,len(Q),len(ANGLES)),complex)
        observed_n=np.zeros((2,N,4,4)); diag=np.zeros((N,len(Q),2,2)); cluster=[]
        for j,idx in enumerate(chosen):
            g=groups[idx];ob=simulate_group(g,model,world,rng)
            if ob is None: continue
            diag[j]=diagonal(ob['z'],ob['n'])
            for a,crng in enumerate(code_rngs):
                code=crng.integers(4,size=len(g['n']))
                z=ob['z']
                if world=='code_offset':
                    shift=np.array(SETTINGS['code_offsets'])[code]
                    z=z*np.exp(1j*shift[:,None,None,None]*Q[None,None,:,None]*U[None,None,None,:,0])
                observed_z[a,j],observed_n[a,j]=aggregate_observations(z,ob['n'],code)
                if world=='homogeneous' and rep<8 and j<16:
                    c=clustering(ob,code,rng,meta['qc'],meta['cluster_radius'],mass[idx]);c.update(group=j,template=int(idx),assignment=a)
                    cluster.append(c)
        coded=np.stack([estimate(observed_z[a],observed_n[a]) for a in range(2)])
        filename=f'field-{rep:03d}.npz'
        atomic_npz(folder/filename,coded=coded,known_source=diag,template_ids=chosen)
        record=dict(rep=rep,physical_seed=seed(name,rep,'physical'),code_seeds=[seed(name,rep,'codes',a) for a in range(2)],
                    filename=filename,sha256=sha(folder/filename),seconds=time.time()-started,
                    cases=field_records(coded,diag,target,ng,rep),clustering=cluster)
        save(folder/f'field-{rep:03d}.json',record)
        state['replicates'].append(dict(rep=rep,filename=filename,sha256=record['sha256'],record_sha256=sha(folder/f'field-{rep:03d}.json')))
        state.update(completed=rep+1,complete=rep+1==expected)
        save(statepath,state)
        save(P/'execution-state.json',dict(pid=os.getpid(),world=name,completed=rep+1,expected=expected,stage='running',last_update_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        print(name,rep+1,'/',expected,'seconds',round(record['seconds'],2),flush=True)
    return state


def main():
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true');p.add_argument('--tag',choices=TAGS);p.add_argument('--model',choices=MODELS);p.add_argument('--world',choices=WORLDS);p.add_argument('--all',action='store_true');args=p.parse_args()
    if args.freeze: freeze();return
    frozen=json.loads((P/'freeze.json').read_text());assert frozen['code_sha256']==sha(__file__)
    from controls import main as controls
    controls()
    configurations=list(itertools.product(TAGS,MODELS,WORLDS)) if args.all else [(args.tag,args.model,args.world)]
    assert all(all(v is not None for v in cfg) for cfg in configurations)
    previous=None; prepared=None
    for tag,model,world in configurations:
        if tag!=previous: prepared=templates(tag); previous=tag
        run(tag,model,world,prepared)
    save(P/'execution-state.json',dict(pid=os.getpid(),stage='complete',worlds=len(configurations),finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))

if __name__=='__main__':main()
