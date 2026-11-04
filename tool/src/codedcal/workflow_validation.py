"""Budgeted computational qualification of image-analysis workflows.

Known-inverse references define computational consistency, never RNA truth or
physical accuracy. Adapter execution never receives the imposed maps/reference.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol
import hashlib
import importlib
import itertools
import json
import math
import time
import numpy as np
from scipy.ndimage import shift
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree


def spatial_match(a,b,radius,max_component=2000):
    """Three-dimensional cardinality-first match; retain unresolved components."""
    a=np.asarray(a,float).reshape(-1,3);b=np.asarray(b,float).reshape(-1,3)
    if not len(a) or not len(b):return np.empty((0,2),int),dict(unresolved_components=0,ambiguous_reference=0)
    neighbors=cKDTree(a).query_ball_tree(cKDTree(b),radius)
    reverse=[[] for _ in b]
    for i,js in enumerate(neighbors):
        for j in js:reverse[j].append(i)
    visited=set();pairs=[];unresolved=0
    for initial,js in enumerate(neighbors):
        if not js or initial in visited:continue
        left={initial};right=set();stack=[initial];visited.add(initial)
        while stack:
            i=stack.pop()
            for j in neighbors[i]:
                if j in right:continue
                right.add(j)
                for k in reverse[j]:
                    if k not in visited:visited.add(k);left.add(k);stack.append(k)
        aa=sorted(left);bb=sorted(right)
        if len(aa)+len(bb)>max_component:unresolved+=1;continue
        distance=np.linalg.norm(a[aa,None,:]-b[None,bb,:],axis=2)
        reward=(min(len(aa),len(bb))+1)*radius
        cost=np.zeros((len(aa)+len(bb),len(aa)+len(bb)))
        cost[:len(aa),:len(bb)]=np.where(distance<=radius,distance-reward,reward)
        rows,cols=linear_sum_assignment(cost)
        for i,j in zip(rows,cols):
            if i<len(aa) and j<len(bb) and distance[i,j]<=radius:pairs.append((aa[i],bb[j]))
    return np.array(pairs,int).reshape(-1,2),dict(unresolved_components=unresolved,ambiguous_reference=sum(len(js)>1 for js in neighbors))


@dataclass
class ImageCase:
    images: np.ndarray
    support: np.ndarray
    codebook: np.ndarray
    rounds: np.ndarray
    colors: np.ndarray
    class_ids: np.ndarray
    blank_classes: np.ndarray
    name: str = 'input'
    preparation: str = 'unspecified'
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.images=np.asarray(self.images);self.support=np.asarray(self.support)
        self.codebook=np.asarray(self.codebook,float)
        for value in [self.rounds,self.colors,self.class_ids,self.blank_classes]:
            a=np.asarray(value)
            if not np.isfinite(a).all() or np.any(a!=np.floor(a)):raise ValueError('Readout and class IDs must be integers')
        self.rounds=np.asarray(self.rounds,int);self.colors=np.asarray(self.colors,int)
        self.class_ids=np.asarray(self.class_ids,int);self.blank_classes=np.asarray(self.blank_classes,int)
        if self.images.ndim!=4 or self.images.shape!=self.support.shape:raise ValueError('Expected matching FZYX images/support')
        f=len(self.images)
        if self.codebook.ndim!=2 or self.codebook.shape[0]!=f or self.rounds.shape!=(f,) or self.colors.shape!=(f,):raise ValueError('Incompatible readout metadata')
        if not np.isfinite(self.images).all() or np.any(self.images<0) or not np.isfinite(self.codebook).all() or np.any(self.codebook<0):raise ValueError('Nonnegative finite observations required')
        if not np.isfinite(self.support).all() or np.any((self.support<0)|(self.support>1)):raise ValueError('Support must be in [0,1]')
        if np.any(self.codebook.sum(axis=0)<=0) or self.class_ids.shape!=(self.codebook.shape[1],) or np.any(self.class_ids<0):raise ValueError('Empty or unlabelled codeword')

    def replace_images(self, images, support):
        return ImageCase(images,support,self.codebook,self.rounds,self.colors,self.class_ids,self.blank_classes,self.name,self.preparation,self.metadata)

    def for_adapter(self):
        # Provenance paths to historical catalogs stay outside the execution API.
        metadata={k:v for k,v in self.metadata.items() if k not in {'historical_reference','reference_catalog','imposed_maps','challenge_truth'}}
        return ImageCase(self.images,self.support,self.codebook,self.rounds,self.colors,self.class_ids,self.blank_classes,self.name,self.preparation,metadata)

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as a:
            data={k:a[k] for k in ['images','support','codebook','rounds','colors','class_ids','blank_classes']}
            for k in ['name','preparation']:
                if k in a:data[k]=str(a[k].item())
            if 'metadata_json' in a:data['metadata']=json.loads(str(a['metadata_json'].item()))
        return cls(**data)


@dataclass
class Catalog:
    points: np.ndarray
    classes: np.ndarray
    def __post_init__(self):
        self.points=np.asarray(self.points,float).reshape(-1,3)
        ids=np.asarray(self.classes)
        if not np.isfinite(ids).all() or np.any(ids!=np.floor(ids)):raise ValueError('Catalog class IDs must be integers')
        self.classes=np.asarray(ids,int)
        if self.classes.shape!=(len(self.points),) or not np.isfinite(self.points).all() or np.any(self.classes<0):raise ValueError('Invalid ZYX catalog')
    def as_dict(self):return dict(points=self.points.tolist(),classes=self.classes.tolist())
    @classmethod
    def from_dict(cls,d):return cls(d['points'],d['classes'])


@dataclass
class WorkflowResult:
    catalog: Catalog
    corrections: np.ndarray
    coverage: float
    diagnostics: dict=field(default_factory=dict)
    def __post_init__(self):
        self.corrections=np.asarray(self.corrections,float)
        if not np.isfinite(self.coverage) or not 0<=self.coverage<=1:raise ValueError('Invalid support coverage')


class WorkflowAdapter(Protocol):
    def run(self, case: ImageCase, policy: str) -> WorkflowResult: ...


def warp(images,maps):
    """Nonnegative order-one resampling, fixed zero extension, no periodic wrap."""
    images=np.asarray(images);maps=np.asarray(maps,float)
    if maps.shape!=(len(images),3) or not np.isfinite(maps).all():raise ValueError('Expected F by ZYX shifts')
    return np.stack([shift(im,m,order=1,mode='constant',cval=0.,prefilter=False) for im,m in zip(images,maps)])


def canonical(maps):
    maps=np.asarray(maps,float)
    if maps.ndim!=2 or maps.shape[1]!=3 or not np.isfinite(maps).all():raise ValueError('Invalid corrections')
    return maps-maps[:1]


def apply_corrections(case,maps):
    """Fix the first-readout gauge before the decoder sees any resampled input."""
    maps=canonical(maps)
    return case.replace_images(warp(case.images,maps),warp(case.support.astype(np.float32),maps)),maps


def pad_for_challenges(case,margin=16):
    """Retain all shifted original pixels; padding is known zero extension.

    Calibration and decoding keep their original contexts. No original pixel is
    reconstructed or filled in. Integer inverse references are therefore exact
    over the entire padded array, not just an optimistically selected interior.
    """
    if 'challenge_padding_yx' in case.metadata:return case
    z,y,x=case.images.shape[1:];metadata=dict(case.metadata)
    defaults={'context':[[0,z],[16,y-16],[16,x-16]],
        'scoring_region':[[0,z],[24,y-24],[24,x-24]],
        'calibration_region':[[0,z],[0,y],[0,x]]}
    defaults['coverage_region']=metadata.get('coverage_region',metadata.get('context',defaults['context']))
    for key,value in defaults.items():
        region=np.asarray(metadata.get(key,value),int).copy();region[1:]+=margin;metadata[key]=region.tolist()
    metadata['challenge_padding_yx']=margin
    pad=((0,0),(0,0),(margin,margin),(margin,margin))
    return ImageCase(np.pad(case.images,pad),np.pad(case.support,pad),case.codebook,case.rounds,case.colors,
        case.class_ids,case.blank_classes,case.name,case.preparation,metadata)


def compare(reference,prediction,blank_classes=(),radius=3.):
    """Label-aware maximum-cardinality matching, plus assignment-free composition."""
    if radius<=0 or not np.isfinite(radius):raise ValueError('Positive matching radius required')
    negative=np.asarray(blank_classes,int)
    ar=~np.isin(reference.classes,negative);br=~np.isin(prediction.classes,negative)
    ap,ac=reference.points[ar],reference.classes[ar];bp,bc=prediction.points[br],prediction.classes[br]
    nr,npred=len(ac),len(bc)
    if nr==0:return dict(supported=False,reason='empty_nonblank_reference',reference_calls=0,predicted_calls=npred)
    correct=0;unresolved=0
    for word in np.union1d(ac,bc):
        pairs,info=spatial_match(ap[ac==word],bp[bc==word],radius);correct+=len(pairs);unresolved+=info['unresolved_components']
    keys=np.union1d(ac,bc);aa=np.array([np.sum(ac==k) for k in keys])/nr
    bb=np.array([np.sum(bc==k) for k in keys])/max(npred,1)
    tv=float(np.abs(aa-bb).sum()/2) if npred else 1.
    # A separate spatial-first matching sensitivity retains label ambiguity.
    pairs,diag=spatial_match(ap,bp,radius)
    spatial_correct=int(np.sum(ac[pairs[:,0]]==bc[pairs[:,1]])) if len(pairs) else 0
    return dict(supported=unresolved==0,identity_loss=1-2*correct/(nr+npred),composition_tv=tv,unresolved_identity_components=unresolved,
        spatial_first_identity_loss=None if diag['unresolved_components'] else 1-2*spatial_correct/(nr+npred),
        spatial_first_supported=diag['unresolved_components']==0,matching=diag,
        reference_calls=nr,predicted_calls=npred,blank_calls=int((~br).sum()),
        yield_ratio=npred/nr,empty_prediction=npred==0)


def assess(reference,result,blank_classes,minimum_reference=10,minimum_coverage=.99,radius=3.,tolerance=.05):
    nref=int(np.sum(~np.isin(reference.classes,blank_classes)))
    if nref<minimum_reference:
        return dict(status='unsupported',reason='insufficient_reference_calls',failed=None,reference_calls=nref)
    if result.coverage<minimum_coverage:
        return dict(status='workflow_failed',reason='insufficient_output_support',failed=True,coverage=float(result.coverage))
    d=compare(reference,result.catalog,blank_classes,radius)
    if not d['supported'] or d['reference_calls']<minimum_reference:
        return d|dict(status='unsupported',reason='unresolved_matching' if not d['supported'] else 'insufficient_reference_calls',failed=None)
    return dict(status='scored',failed=d['identity_loss']>tolerance or d['composition_tv']>tolerance,coverage=float(result.coverage),**d)


def _key(mask):
    m=np.asarray(mask,float)
    # Sign reversals have the same readout split; direction is specified separately.
    a=tuple(np.round(m,12));b=tuple(np.round(-m,12));return min(a,b)


def mask_bank(case,seed=110126):
    """A common empirical shift histogram across round/color/frame allocations."""
    f=len(case.images);rs=np.unique(case.rounds);cs=np.unique(case.colors)
    g=math.gcd(len(rs),len(cs))
    levels=np.array([-1.,1.]) if g%2==0 else (np.array([-1.,0.,1.]) if g%3==0 else None)
    if levels is None or f%len(levels):raise ValueError('No supported common round/color histogram')
    records=[];seen=set()
    def add(mask,family,label):
        mask=np.asarray(mask,float);key=_key(mask)
        if len(mask)!=f or not np.array_equal(np.sort(mask),np.repeat(levels,f//len(levels))):return
        if key not in seen:records.append(dict(mask=mask.tolist(),family=family,label=label));seen.add(key)
    for family,ids in [('round',case.rounds),('color',case.colors)]:
        group=np.unique(ids);signature=np.tile(levels,len(group)//len(levels))
        assignments=[]
        if len(group)<=6:assignments=sorted(set(itertools.permutations(signature)))
        else:
            assignments=[signature,np.sort(signature),np.roll(signature,1)]
            rng=np.random.default_rng(seed+(0 if family=='round' else 1))
            assignments += [rng.permutation(signature) for _ in range(12)]
        for k,values in enumerate(assignments):
            add(np.array([values[np.searchsorted(group,i)] for i in ids]),family,f'{family}-{k}')
    round_masks=[r for r in records if r['family']=='round'];color_masks=[r for r in records if r['family']=='color']
    if not round_masks or not color_masks:raise ValueError('Round/color challenge unavailable')
    for i,r in enumerate(round_masks[:4]):
        for j,c in enumerate(color_masks[:4]):
            product=np.array(r['mask'])*np.array(c['mask'])
            add(product,'interaction',f'interaction-{i}-{j}')
    signature=np.repeat(levels,f//len(levels));rng=np.random.default_rng(seed+2)
    for k in range(64):add(rng.permutation(signature),'frame',f'balanced-{k}')
    return records


def suites(case,budget=4,seed=110126):
    """Historical unconstrained strategies, retained for v1 reproduction only."""
    if budget!=4:raise ValueError('This specified design uses four challenge pairs')
    bank=mask_bank(case,seed);families={k:[i for i,r in enumerate(bank) if r['family']==k] for k in ['round','color','interaction','frame']}
    def fill(ids):
        ids=list(dict.fromkeys(ids))
        for i in families['round']+families['color']+families['frame']:
            if len(ids)>=budget:break
            if i not in ids:ids.append(i)
        if len(ids)!=budget:raise ValueError('Insufficient distinct probes')
        return ids
    standard=fill([families['round'][0]]+families['color'][:3])
    mixed=fill([families['round'][0],families['color'][0]]+families['interaction'][:1]+families['frame'][:1])
    random=fill(families['frame'][:4])
    # No reference catalog, injected held-out error or scored outcome enters this weight.
    b=case.codebook/np.maximum(case.codebook.sum(axis=0),1e-20)
    energy=np.sqrt(np.mean(np.square(case.images.astype(float)),axis=(1,2,3)))
    gram=b@b.T;gram*=np.sqrt(np.outer(energy,energy));ii,jj=np.triu_indices(len(b),1);weights=gram[ii,jj]
    if weights.sum()==0:raise ValueError('No supported code/readout coupling')
    features=np.array([np.square(np.array(r['mask'])[ii]-np.array(r['mask'])[jj])/4 for r in bank])
    chosen=[];covered=np.zeros(len(weights))
    for allowed in [families['round'],families['color'],list(range(len(bank))),list(range(len(bank)))]:
        candidates=[i for i in allowed if i not in chosen]
        selected=min(candidates,key=lambda i:(-float(np.maximum(features[i]-covered,0)@weights),i))
        chosen.append(selected);covered=np.maximum(covered,features[selected])
    return dict(bank=bank,strategies=dict(conventional=standard,mixed_conventional=mixed,
                balanced_random=random,readout_aware=chosen),budget_pairs=budget,seed=seed,
                selection_uses_reference_calls=False,selection_uses_heldout_outcomes=False)


def legacy_challenge_plan(case,seed=110126):
    selection=suites(case,seed=seed);bank=selection['bank']
    directions=[(0.,2.),(2.,0.),(2.,2.),(-2.,2.)];records={};strategy_ids={}
    def add(mask,direction,label,role):
        maps=np.zeros((len(case.images),3));maps[:,1:]=np.array(mask)[:,None]*np.array(direction)
        maps=canonical(maps)
        token=hashlib.sha256(maps.astype('<f8').tobytes()).hexdigest()
        if token not in records:records[token]=dict(id=token,label=label,role=role,maps=maps.tolist())
        return token
    for method,indices in selection['strategies'].items():
        strategy_ids[method]=[add(bank[k]['mask'],d,f"{bank[k]['label']}-probe{j}",'probe') for j,(k,d) in enumerate(zip(indices,directions))]
    groups={k:[r for r in bank if r['family']==k] for k in ['round','color','interaction','frame']}
    # Five held-out interventions; exact maps differ from every probe.
    specs=[(groups['round'][-1],(1.,0.)),(groups['color'][-1],(0.,3.)),
           ((groups['interaction'] or groups['frame'])[-1],(1.,1.)),
           (groups['frame'][-1],(-3.,3.)),(groups['frame'][-2],(1.,-1.))]
    heldout=[add(r['mask'],d,r['label']+'-heldout','heldout') for r,d in specs]
    probes={k for ids in strategy_ids.values() for k in ids}
    if probes.intersection(heldout):raise AssertionError('Held-out challenge overlaps a probe')
    fractional=[add(groups['frame'][-3]['mask'],(.5,.25),'fractional-frame','fractional'),
                add(groups['color'][0]['mask'],(-.5,.25),'fractional-color','fractional')]
    return dict(selection=selection,strategies=strategy_ids,heldout=heldout,fractional=fractional,
                challenges=list(records.values()),budget_pairs=4,logical_executions_per_workflow=47,
                scope='Five held-out digital challenges, not all possible errors or new physical observations')


def challenge_plan(case,seed=110126,strategy='mixed_conventional',include_holdout=False,include_fractional=False):
    from .workflow_design import plan_challenges
    return plan_challenges(case,seed,strategy,include_holdout,include_fractional)


def import_adapter(spec,config):
    module,name=spec.split(':',1)
    return getattr(importlib.import_module(module),name)(config)


def _json_safe(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return _json_safe(x.item())
    if isinstance(x,dict):return {str(k):_json_safe(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [_json_safe(v) for v in x]
    if isinstance(x,float) and not np.isfinite(x):return None
    return x


def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(_json_safe(value),indent=2,allow_nan=False)+'\n');tmp.replace(path)


def qualify_input(case,adapter,minimum_native=50,native_policy='identity'):
    """Native, no-effect and exact-inverse controls, before comparative outcomes.

    These qualify the declared computational input/target, not biological truth.
    Fractional interpolation is separately checked against the composed operator.
    Four logical calls are charged even when an exact repeat is cached.
    """
    if not isinstance(native_policy,str) or not native_policy.strip():raise ValueError('A nonempty native policy is required')
    case=pad_for_challenges(case)
    t=np.zeros((len(case.images),3));t[:,1]=np.arange(len(t))%2*2
    observed=case.replace_images(warp(case.images,t),warp(case.support,t))
    restored,_=apply_corrections(observed,-t)
    context=adapter.context(case) if hasattr(adapter,'context') else tuple(slice(8,-8) if n>16 else slice(None) for n in case.images.shape[1:])
    exact_images=np.array_equal(restored.images[(slice(None),)+context],case.images[(slice(None),)+context])
    exact_support=np.array_equal(restored.support[(slice(None),)+context],case.support[(slice(None),)+context])
    # Explicitly exercise the no-effect transformation, not just an alias.
    zero=case.replace_images(warp(case.images,np.zeros_like(t)),warp(case.support,np.zeros_like(t)))
    calls=[]
    for role,value in [('native',case),('no_effect',zero),('known_inverse',restored)]:
        start=time.perf_counter()
        try:
            r=adapter.run(value.for_adapter(),native_policy)
            if r.corrections.shape!=t.shape or not np.isfinite(r.corrections).all() or not np.allclose(r.corrections[0],0,atol=1e-10):raise ValueError('Invalid native workflow correction gauge')
            calls.append(dict(role=role,status='completed',catalog=r.catalog.as_dict(),corrections=r.corrections,coverage=r.coverage,diagnostics=r.diagnostics,wall_seconds=time.perf_counter()-start))
        except Exception as e:
            import traceback
            calls.append(dict(role=role,status='failed',error=type(e).__name__+': '+str(e),traceback=traceback.format_exc(),wall_seconds=time.perf_counter()-start))
    # Fourth call is the fractional composed oracle: its pixels need not equal native.
    fractional=t/4;obs=case.replace_images(warp(case.images,fractional),warp(case.support,fractional))
    composed,_=apply_corrections(obs,-fractional)
    analytic=np.array_equal(composed.images,warp(warp(case.images,fractional),canonical(-fractional)))
    start=time.perf_counter()
    try:
        r=adapter.run(composed.for_adapter(),native_policy)
        if r.corrections.shape!=t.shape or not np.isfinite(r.corrections).all() or not np.allclose(r.corrections[0],0,atol=1e-10):raise ValueError('Invalid native workflow correction gauge')
        calls.append(dict(role='fractional_oracle',status='completed',catalog=r.catalog.as_dict(),corrections=r.corrections,coverage=r.coverage,diagnostics=r.diagnostics,wall_seconds=time.perf_counter()-start))
    except Exception as e:calls.append(dict(role='fractional_oracle',status='failed',error=type(e).__name__+': '+str(e),wall_seconds=time.perf_counter()-start))
    native=calls[0];count=0 if native['status']!='completed' else int(np.sum(~np.isin(np.array(native['catalog']['classes']),case.blank_classes)))
    same=all(x['status']=='completed' and x['catalog']==native.get('catalog') for x in calls[1:3])
    passed=all(x['status']=='completed' and x['coverage']>=.99 for x in calls[:3]) and count>=minimum_native and exact_images and exact_support and same and analytic
    fractional=calls[3]
    fractional_count=0 if fractional['status']!='completed' else int(np.sum(~np.isin(fractional['catalog']['classes'],case.blank_classes)))
    return dict(passed=passed,name=case.name,preparation=case.preparation,native_policy=native_policy,blank_classes=case.blank_classes,native_nonblank_calls=count,minimum_native=minimum_native,
        integer_inverse_images_exact=exact_images,integer_inverse_support_exact=exact_support,no_effect_and_inverse_catalogs_exact=same,
        fractional_composed_operator_exact=analytic,calls=calls,logical_workflow_executions=4,
        fractional_input_supported=fractional['status']=='completed' and fractional['coverage']>=.99 and fractional_count>=10,
        fractional_reference_calls=fractional_count,minimum_reference=10,
        actual_workflow_executions=sum(x.get('diagnostics',{}).get('actual_workflow_executions',1) for x in calls),physical_accuracy_claim=False)


def run_challenge(case,adapter,policy,challenge,minimum_reference=10):
    """Two full workflow calls. The imposed maps are kept outside the adapter."""
    case=pad_for_challenges(case)
    t=np.array(challenge['maps'],float);observed=case.replace_images(warp(case.images,t),warp(case.support.astype(np.float32),t))
    restored,_=apply_corrections(observed,-t)
    records=[]
    for role,input_case in [('challenged',observed),('oracle_reference',restored)]:
        start=time.perf_counter()
        try:
            result=adapter.run(input_case.for_adapter(),policy)
            if result.corrections.shape!=t.shape or not np.isfinite(result.corrections).all():raise ValueError('Invalid workflow corrections')
            if not np.allclose(result.corrections[0],0,atol=1e-10):raise ValueError('Adapter did not canonicalize before decoding')
            row=dict(status='completed',catalog=result.catalog.as_dict(),corrections=result.corrections,
                coverage=float(result.coverage),diagnostics=result.diagnostics)
        except Exception as e:
            import traceback
            result=None;row=dict(status='failed',error=type(e).__name__+': '+str(e),traceback=traceback.format_exc())
        row.update(role=role,wall_seconds=time.perf_counter()-start,actual_workflow_executions=row.get('diagnostics',{}).get('actual_workflow_executions',1));records.append(row)
    a,b=records
    nref=0 if b['status']!='completed' else int(np.sum(~np.isin(b['catalog']['classes'],case.blank_classes)))
    if b['status']!='completed':assessment=dict(status='unsupported',reason='oracle_reference_execution_failed',failed=None)
    elif b['coverage']<.99:assessment=dict(status='unsupported',reason='oracle_reference_support_failed',failed=None)
    elif nref<minimum_reference:assessment=dict(status='unsupported',reason='insufficient_reference_calls',failed=None,reference_calls=nref)
    elif a['status']!='completed':assessment=dict(status='workflow_failed',reason='challenged_execution_failed',failed=True,reference_calls=nref)
    else:
        reference=Catalog.from_dict(b['catalog']);r=WorkflowResult(Catalog.from_dict(a['catalog']),np.array(a['corrections']),a['coverage'],a['diagnostics'])
        assessment=assess(reference,r,case.blank_classes,minimum_reference=minimum_reference)
        assessment['matching_sensitivity']={str(rad):compare(reference,r.catalog,case.blank_classes,rad) for rad in [1.5,5.]}
        residual=canonical(np.array(a['corrections'])+t-np.array(b['corrections']))
        assessment['known_probe_map_rms']=float(np.sqrt(np.mean(np.sum(residual**2,axis=1))))
    return dict(challenge=challenge,policy=policy,runs=records,assessment=assessment,
                logical_workflow_executions=2,actual_workflow_executions=sum(x['actual_workflow_executions'] for x in records),physical_accuracy_claim=False)


def run_suite(case,adapter,policies,out,plan=None,deadline_seconds=10800):
    """Run each policy against its own qualified native catalog.

    The deadline is checked between calls; adapters must bound an in-flight call.
    Every completed native/trial record survives an interruption.
    """
    out=Path(out)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Choose an empty trial directory; evidence is never overwritten')
    if not math.isfinite(deadline_seconds) or deadline_seconds<=0:raise ValueError('Positive finite deadline required')
    plan=plan or challenge_plan(case)
    from .workflow_design import require_design
    require_design(case,plan)
    out.mkdir(parents=True,exist_ok=True)
    policies=list(policies)
    if len(set(policies))!=len(policies) or not policies or any(not p or Path(p).name!=p or p in ('.','..') for p in policies):raise ValueError('Unique policy names without path separators required')
    write_json(out/'plan.json',plan);start=time.monotonic();rows=[];case=pad_for_challenges(case)
    state=dict(stage='running',policies=policies,completed_pairs=0,logical_workflow_dispatches=0,native_gates={})
    def save():
        state['elapsed_seconds']=time.monotonic()-start;write_json(out/'state.json',state)
    def deadline():
        if time.monotonic()-start>deadline_seconds:raise TimeoutError('Finite study deadline reached; partial evidence retained')
    save()
    try:
        for policy in policies:
            deadline();path=out/policy/'native.json'
            if path.exists():raise FileExistsError('Do not overwrite trial evidence: '+str(path))
            state.update(current_policy=policy,current_challenge='native');state['logical_workflow_dispatches']+=1;save();begin=time.perf_counter()
            try:
                r=adapter.run(case.for_adapter(),policy)
                if r.corrections.shape!=(len(case.images),3) or not np.isfinite(r.corrections).all() or not np.allclose(r.corrections[0],0,atol=1e-10):raise ValueError('Invalid native workflow correction gauge')
                native=dict(status='completed',catalog=r.catalog.as_dict(),corrections=r.corrections,coverage=r.coverage,diagnostics=r.diagnostics)
            except Exception as e:native=dict(status='failed',error=type(e).__name__+': '+str(e))
            count=0 if native['status']!='completed' else int(np.sum(~np.isin(native['catalog']['classes'],case.blank_classes)))
            native.update(policy=policy,wall_seconds=time.perf_counter()-begin,nonblank_calls=count,
                native_workflow_qualified=native['status']=='completed' and native['coverage']>=.99 and count>=50,
                actual_workflow_executions=native.get('diagnostics',{}).get('actual_workflow_executions',1),logical_workflow_executions=1)
            write_json(path,native);state['native_gates'][policy]=native['native_workflow_qualified'];save()
            if not native['native_workflow_qualified']:continue
            for challenge in plan['challenges']:
                deadline();path=out/policy/(challenge['id']+'.json')
                if path.exists():raise FileExistsError('Do not overwrite or silently reuse trial evidence: '+str(path))
                state['current_challenge']=challenge['id'];state['logical_workflow_dispatches']+=2;save()
                row=run_challenge(case,adapter,policy,challenge);oracle=row['runs'][1]
                if oracle['status']=='completed':
                    row['oracle_vs_native']=compare(Catalog.from_dict(native['catalog']),Catalog.from_dict(oracle['catalog']),case.blank_classes)
                    if challenge['role']!='fractional' and (oracle['catalog']!=native['catalog'] or not np.array_equal(oracle['corrections'],native['corrections'])):
                        row['assessment']=dict(status='unsupported',reason='integer_oracle_native_mismatch',failed=None)
                write_json(path,row);rows.append(row);state['completed_pairs']+=1;save()
        state['stage']='completed';save()
    except BaseException as e:
        state.update(stage='interrupted',error=type(e).__name__+': '+str(e));save();raise
    return rows
