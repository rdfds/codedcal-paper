"""Pre-execution checks of representative readout-error mechanisms.

This is ordinary linear projection and constrained coverage, not a new
identification theorem or a guarantee about an untested imaging error.
"""
from pathlib import Path
import hashlib
import numpy as np

FAMILIES=('round','color','interaction','frame')
STRATEGIES=('mixed_conventional','mechanism_weighted')


def _basis(case):
    return np.column_stack([np.eye(len(np.unique(ids)))[np.searchsorted(np.unique(ids),ids)]
                            for ids in (case.rounds,case.colors)])


def _residual(values,basis):
    return values-basis@np.linalg.lstsq(basis,values,rcond=None)[0]


def _valid(case,values,family,interaction=None):
    values=np.asarray(values,float)
    if values.shape!=(len(case.images),3) or not np.isfinite(values).all():return False
    centered=values-values.mean(axis=0);norm=float(np.linalg.norm(centered));tol=1e-9*max(1.,norm)
    if norm<=tol:return False
    def group_residual(ids):
        return np.linalg.norm(values-np.stack([values[ids==i].mean(axis=0) for i in ids]))
    rr=group_residual(case.rounds);cc=group_residual(case.colors)
    if family=='round':return rr<=tol and cc>tol
    if family=='color':return cc<=tol and rr>tol
    additive=_basis(case)
    if np.linalg.norm(_residual(values,additive))<=tol:return False
    if family=='interaction':return True
    if family=='frame' and interaction is not None:
        if np.asarray(interaction).shape!=values.shape or not np.isfinite(interaction).all():return False
        return np.linalg.norm(_residual(values,np.column_stack([additive,interaction])))>tol
    return False


def audit_plan(case,plan):
    """Check actual maps; family labels alone cannot satisfy the guard."""
    errors=[];strategies={};basis=_basis(case);rank=int(np.linalg.matrix_rank(basis))
    challenges=plan.get('challenges',[]);by_id={r.get('id'):r for r in challenges}
    if len(by_id)!=len(challenges):errors.append('duplicate_challenge_ids')
    if not plan.get('strategies'):errors.append('no_probe_strategy')
    for name,ids in plan.get('strategies',{}).items():
        problems=[];families={};maps=[];interaction=None
        if len(ids)!=4 or len(set(ids))!=4:problems.append('four_distinct_probe_pairs_required')
        for token in ids:
            record=by_id.get(token)
            if record is None:problems.append('missing_probe_record');continue
            family=record.get('mechanism_family')
            values=np.asarray(record.get('maps'),float)
            if family not in FAMILIES:problems.append('undeclared_mechanism_family');continue
            if family in families:problems.append('duplicate_mechanism_family:'+family)
            families[family]=values
        interaction=families.get('interaction')
        for family in FAMILIES:
            if family not in families:problems.append('missing_mechanism_family:'+family);continue
            values=families[family]
            if not _valid(case,values,family,interaction):problems.append('invalid_mechanism_maps:'+family)
            elif values.shape==(len(case.images),3):maps.append(values)
        sampled=max(0,int(np.linalg.matrix_rank(np.column_stack([basis]+maps)))-rank) if maps else 0
        strategies[name]=dict(supported=not problems,reasons=problems,
            representative_families=list(families),sampled_nonadditive_rank=sampled)
        errors.extend(name+':'+p for p in problems)
    return dict(supported=not errors,reasons=errors,strategies=strategies,
        readout_count=len(case.images),round_color_rank=rank,
        nonadditive_allocation_dimension=len(case.images)-rank,
        scope='One representative per declared family; residual rank does not establish coverage of all directions, magnitudes or nonlinear workflow failures.')


def _selection(case,bank,strategy):
    if strategy not in STRATEGIES:raise ValueError('Unknown challenge strategy: '+strategy)
    # The mixed default is independent of observed outcomes and signal weights.
    weights=features=None
    if strategy=='mechanism_weighted':
        b=case.codebook/np.maximum(case.codebook.sum(axis=0),1e-20)
        energy=np.sqrt(np.mean(np.square(case.images.astype(float)),axis=(1,2,3)))
        gram=(b@b.T)*np.sqrt(np.outer(energy,energy));ii,jj=np.triu_indices(len(b),1);weights=gram[ii,jj]
        if not np.isfinite(weights).all() or weights.sum()<=0:raise ValueError('No finite supported code/readout coupling for weighted selection')
        features=np.array([(np.asarray(r['mask'])[ii]-np.asarray(r['mask'])[jj])**2/4 for r in bank])
    chosen=[];covered=None if weights is None else np.zeros(len(weights));interaction=None
    for family in FAMILIES:
        candidates=[]
        for i,record in enumerate(bank):
            if record['family']!=family:continue
            values=np.column_stack([np.zeros(len(case.images)),record['mask'],np.zeros(len(case.images))])
            if _valid(case,values,family,interaction):candidates.append(i)
        if not candidates:raise ValueError('No supported representative for mechanism family: '+family)
        selected=candidates[0] if weights is None else min(candidates,key=lambda i:(-float(np.maximum(features[i]-covered,0)@weights),i))
        chosen.append(selected)
        if weights is not None:covered=np.maximum(covered,features[selected])
        if family=='interaction':interaction=np.column_stack([np.zeros(len(case.images)),bank[selected]['mask'],np.zeros(len(case.images))])
    return chosen


def plan_challenges(case,seed=110126,strategy='mixed_conventional',include_holdout=False,include_fractional=False):
    from .workflow_validation import mask_bank,canonical
    plan=dict(design_version=2,selected_strategy=strategy,budget_pairs=4,seed=seed,
        strategies={},heldout=[],fractional=[],challenges=[],
        development_status='Post-hoc design repair after the exposed negative comparison; no new comparative validation',
        scope='Preservation of native catalogs on specified digital challenges; no general safety, natural-error risk or physical-accuracy conclusion.')
    try:
        bank=mask_bank(case,seed);chosen=_selection(case,bank,strategy)
    except ValueError as error:
        plan.update(design_status='unsupported',design_audit=dict(supported=False,reasons=[str(error)]),logical_executions_per_workflow=0)
        return plan
    records={}
    def add(record,direction,role,label):
        maps=np.zeros((len(case.images),3));maps[:,1:]=np.asarray(record['mask'])[:,None]*np.asarray(direction);maps=canonical(maps)
        token=hashlib.sha256(maps.astype('<f8').tobytes()).hexdigest()
        if token in records:raise ValueError('Challenge allocations must be distinct')
        records[token]=dict(id=token,label=label,role=role,maps=maps.tolist(),mechanism_family=record['family'])
        return token
    directions=[(0.,2.),(2.,0.),(2.,2.),(-2.,2.)]
    plan['strategies'][strategy]=[add(bank[i],d,'probe',bank[i]['label']+'-probe'+str(j)) for j,(i,d) in enumerate(zip(chosen,directions))]
    families={family:[r for r in bank if r['family']==family] for family in FAMILIES}
    if include_holdout:
        specs=[(families['round'][-1],(1.,0.)),(families['color'][-1],(0.,3.)),
               (families['interaction'][-1],(1.,1.)),(families['frame'][-1],(-3.,3.)),(families['frame'][-2],(1.,-1.))]
        plan['heldout']=[add(r,d,'heldout',r['label']+'-heldout') for r,d in specs]
    if include_fractional:
        plan['fractional']=[add(families['frame'][-3],(.5,.25),'fractional','fractional-frame'),
                            add(families['color'][0],(-.5,.25),'fractional','fractional-color')]
    plan.update(challenges=list(records.values()),selection=dict(bank=bank,selected_indices=chosen,
        selection_uses_reference_calls=False,selection_uses_heldout_outcomes=False,
        mechanism_constraints=list(FAMILIES),weighted_variant_comparatively_validated=False),
        logical_executions_per_workflow=1+2*len(records))
    plan['design_audit']=audit_plan(case,plan);plan['design_status']='supported' if plan['design_audit']['supported'] else 'unsupported'
    return plan


def require_design(case,plan):
    audit=audit_plan(case,plan)
    if not audit['supported']:raise ValueError('Unsupported challenge design: '+'; '.join(audit['reasons']))
    return audit


def write_design_report(output,case,plan):
    from .workflow_validation import write_json
    output=Path(output);audit=plan.get('design_audit') or audit_plan(case,plan)
    write_json(output/'design.json',audit)
    lines=['Challenge design: '+('supported' if audit['supported'] else 'unsupported'),
           'Selected strategy: '+str(plan.get('selected_strategy')),
           'Design support is not evidence of workflow accuracy.']
    if 'nonadditive_allocation_dimension' in audit:
        lines.append('Nonadditive allocation dimension: '+str(audit['nonadditive_allocation_dimension']))
    for name,report in audit.get('strategies',{}).items():
        lines.append(name+': '+', '.join(report['representative_families'])+'; sampled nonadditive rank '+str(report['sampled_nonadditive_rank']))
    lines.extend(audit.get('reasons',[]));lines.append(plan['scope'])
    (output/'design.txt').write_text('\n'.join(lines)+'\n')
    return audit
