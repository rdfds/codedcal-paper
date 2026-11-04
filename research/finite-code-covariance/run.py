"""Finite read-only replay and conditional design analysis; writes only --output.

Never calls legacy run(), main(), freeze(), clustering() or save(). Source
reconstruction verifies retained prefix observations before comparative use.
"""
from pathlib import Path
import argparse, csv, hashlib, itertools, json, os, sys, time
import numpy as np
from covariance import four_panel, gradients, cost_variance
ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'outputs/temporal-validation-oct09'
sys.path.insert(0,str(OLD))
import validation as V
KS=(3,4,6,8,12,16,24,32,48,64)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def save(path,value):
    def clean(v):
        if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [clean(x) for x in v]
        if isinstance(v,np.ndarray):return clean(v.tolist())
        if isinstance(v,np.generic):return clean(v.item())
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    tmp=Path(str(path)+'.tmp');tmp.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n');tmp.replace(path)


def replay_close(got,want):
    err=float(np.max(np.abs(got-want))) if np.size(got) else 0.
    scale=max(float(np.max(np.abs(want))),1.) if np.size(want) else 1.
    if err>2e-12*scale:raise AssertionError(('retained replay mismatch',err,scale))
    return dict(absolute=err,relative_scale=err/scale,bit_equal=bool(np.array_equal(got,want)))


def retained(folder,meta):
    state=json.loads((folder/'state.json').read_text());n=meta['groups']
    floor=[];coding=[];records=[];hashes={}
    pop=np.load(folder/'truth.npz')['population'][1].reshape(8)
    for rec in state['replicates']:
        path=folder/rec['filename'];rpath=path.with_suffix('.json')
        for p,key in ((path,'sha256'),(rpath,'record_sha256')):
            digest=sha(p)
            if digest!=rec[key]: raise AssertionError(('source hash mismatch',str(p)))
            hashes[str(p.relative_to(ROOT))]=digest
        with np.load(path) as data:
            diag=data['known_source'][:n].reshape(n,8)
            coded=data['coded'][:,:n].reshape(2,n,8)
            residual=diag-pop
            floor.append(residual.T@residual/n)
            err=coded-diag[None]
            coding.append(np.einsum('agi,agj->ij',err,err)/(2*n))
    return dict(population_mean=pop,physical_floor=np.mean(floor,axis=0),
                observed_code_covariance=np.mean(coding,axis=0),
                physical_floor_by_field=np.asarray(floor),
                code_covariance_by_field=np.asarray(coding),
                physical_fields=len(floor),paired_assignments=2,group_count=n,
                source_hashes=hashes)


def process(tag,model,world,prepared,out):
    name=f'{tag}-{model}-{world}';folder=OLD/'worlds'/name
    groups,meta,_=prepared
    with np.load(folder/'field-000.npz') as stored:
        chosen=stored['template_ids'].copy();diag=stored['known_source'][:16].copy();coded=stored['coded'][:,:16].copy()
    rng=np.random.default_rng(V.seed(name,0,'physical'))
    np.testing.assert_array_equal(rng.integers(meta['groups'],size=max(meta['groups'],512)),chosen)
    crng=[np.random.default_rng(V.seed(name,0,'codes',a)) for a in range(2)]
    check=[];moments=[];monte=np.zeros((512,8));mean=np.zeros(8);pair=np.zeros((8,8));triple=np.zeros((8,8))
    for j,idx in enumerate(chosen[:16]):
        g=groups[idx];ob=V.simulate_group(g,model,world,rng)
        if ob is None:
            check.append(dict(group=j,empty=True,known=replay_close(np.zeros_like(diag[j]),diag[j]),
                              coded=[replay_close(np.zeros_like(coded[a,j]),coded[a,j]) for a in range(2)]))
            moments.append(dict(group=j,empty=True,mean=np.zeros(8),pair=np.zeros((8,8)),triple=np.zeros((8,8))))
            continue
        codes=[r.integers(4,size=len(g['n'])) for r in crng]
        entry=dict(group=j,template=int(idx),sources=len(g['n']),known=replay_close(V.diagonal(ob['z'],ob['n']),diag[j]),coded=[])
        for a in range(2):
            z=ob['z']
            if world=='code_offset':
                shift=np.array(V.SETTINGS['code_offsets'])[codes[a]]
                z=z*np.exp(1j*shift[:,None,None,None]*V.Q[None,None,:,None]*V.U[None,None,None,:,0])
            oz,on=V.aggregate_observations(z,ob['n'],codes[a])
            computed=V.estimate(oz[None],on[None])[0]
            entry['coded'].append(replay_close(computed,coded[a,j]))
        if world!='code_offset':
            contrast=four_panel(ob['z'],ob['n'])
            for a in range(2):
                replay_close(contrast.evaluate(codes[a]),coded[a,j].reshape(8))
            m=contrast.moments()
            if m['triple_min_eigenvalue'] < -1e-10*m['triple_roundoff_scale']:
                raise ArithmeticError('triple covariance loses PSD beyond roundoff')
            replay_close(m['mean'],diag[j].reshape(8))
            mean+=m['mean']/16;pair+=m['pair_gram']/16;triple+=m['triple_gram']/16
            moments.append(dict(group=j,empty=False,mean=m['mean'],pair=m['pair_gram'],triple=m['triple_gram'],
                                min_eigenvalue_scaled=m['triple_min_eigenvalue']/m['triple_roundoff_scale']))
            code_rng=np.random.default_rng(V.seed('finite-code-covariance-oct27',name,j))
            assignments=code_rng.integers(4,size=(512,len(g['n'])))
            for lo in range(0,512,32):monte[lo:lo+32]+=contrast.evaluate(assignments[lo:lo+32])/16
        check.append(entry)
    base=retained(folder,meta)
    hashes=base.pop('source_hashes')
    hashes[str((OLD/'validation.py').relative_to(ROOT))]=sha(OLD/'validation.py')
    hashes[str((folder/'truth.npz').relative_to(ROOT))]=sha(folder/'truth.npz')
    hashes[str((folder/'state.json').relative_to(ROOT))]=sha(folder/'state.json')
    for p in (ROOT/'outputs/reference-reconciliation-oct07'/f'empirical-{tag}.npz',
              ROOT/'outputs/reference-reconciliation-oct07'/f'reference-{tag}.npz'):
        hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(name=name,tag=tag,model=model,world=world,theorem_applicable=world!='code_offset',
                prefix_groups=16,replay=check,retained=base,source_hashes=hashes,
                source_scope='first 16 groups of replicate 0; conditional exposed development observations')
    if world!='code_offset':
        pred=(pair/3+triple/6)/16
        centered=monte-mean
        mc=centered.T@centered/512
        squares=centered**2;se=squares.std(0,ddof=1)/np.sqrt(512)
        result['verification']=dict(predicted_mean=mean,predicted_covariance=pred,mc_second_moment=mc,
            mc_mean=monte.mean(0),mc_mean_se=monte.std(0,ddof=1)/np.sqrt(512),
            variance_mc_se=se,variance_z=(np.diag(mc)-np.diag(pred))/np.maximum(se,1e-300),
            relative_frobenius_error=np.linalg.norm(mc-pred)/max(np.linalg.norm(pred),1e-300),
            assignments=512,seed_rule='sha256 finite-code-covariance-oct27|world|group; independent assignments per group')
        result['prefix_mean']=mean;result['pair_per_group']=pair;result['triple_per_group']=triple
        result['group_moments']=moments
        designs=[]
        for iq,q in enumerate(V.Q):
            sl=slice(iq*4,iq*4+4)
            grad=gradients(base['population_mean'][sl])
            for variant in ('arithmetic','geometric'):
                if grad is None:
                    designs.append(dict(q=float(q),variant=variant,supported=False));continue
                v=grad[variant]
                a=float(v@base['physical_floor'][sl,sl]@v)
                p=float(v@pair[sl,sl]@v);t=float(v@triple[sl,sl]@v)
                if min(a,p,t)<-1e-9*max(a,p,t,1e-300):raise ArithmeticError('negative influence energy')
                # Tiny negative contraction roundoff is retained above; no clipping in design.
                row=dict(q=float(q),variant=variant,supported=True,A=a,P=p,Q=t,
                         coding_fraction_at_K4=(p/3+t/6)/(a+p/3+t/6),
                         physical_floor_source='all retained observed-scope known-source groups',
                         coding_source='fixed replicate-0 16-group prefix')
                scenarios={}
                for policy in ('fixed_groups','linear','quarter_marginal'):
                    grid=[]
                    for k in KS:
                        if policy=='quarter_marginal' and k<4:continue
                        cost=1 if policy=='fixed_groups' else (k if policy=='linear' else 1+.25*(k-4))
                        value=cost_variance(a,p,t,k,cost)
                        grid.append(dict(K=k,cost=cost,variance_cost=value))
                    ref=next(x['variance_cost'] for x in grid if x['K']==4)
                    for point in grid:point['relative_to_four']=point['variance_cost']/ref
                    best=min(grid,key=lambda x:x['variance_cost'])
                    scenarios[policy]=dict(best_K=best['K'],best_ratio=best['relative_to_four'],grid=grid)
                row['scenarios']=scenarios
                field_a=np.einsum('i,fij,j->f',v,base['physical_floor_by_field'][:,sl,sl],v)
                row['A_field_mean_se']=float(field_a.std(ddof=1)/np.sqrt(len(field_a)))
                designs.append(row)
        result['designs']=designs
    save(out/f'{name}.json',result)
    return dict(name=name,applicable=world!='code_offset',prefix_groups=16,
                exact_known=sum(x['known']['bit_equal'] for x in check),
                exact_coded=sum(y['bit_equal'] for x in check for y in x['coded']),
                relative_mc_error=result.get('verification',{}).get('relative_frobenius_error'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);args=parser.parse_args()
    out=Path(args.output).resolve()
    if out.exists(): raise ValueError('require a new output directory; preserve previous attempts')
    out.mkdir(parents=True)
    checks=[];start=time.time()
    try:
        for tag in V.TAGS:
            prepared=V.templates(tag)
            for model,world in itertools.product(V.MODELS,V.WORLDS):
                name=f'{tag}-{model}-{world}';tick=time.time()
                row=process(tag,model,world,prepared,out);checks.append(row)
                save(out/'state.json',dict(stage='running',pid=os.getpid(),completed=checks,seconds=time.time()-start))
                print(name,'done',round(time.time()-tick,2),'seconds',flush=True)
        save(out/'state.json',dict(stage='complete',completed=checks,seconds=time.time()-start,
            implementation_sha256=sha(__file__),covariance_sha256=sha(Path(__file__).with_name('covariance.py'))))
    except Exception as exc:
        save(out/'state.json',dict(stage='failed',pid=os.getpid(),completed=checks,error=repr(exc),seconds=time.time()-start))
        raise

if __name__=='__main__':main()
