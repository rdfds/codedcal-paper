"""Exposed calibration-only variance/design feasibility; no validation scoring."""
from pathlib import Path
import argparse,hashlib,itertools,json,sys
import numpy as np
import networkx as nx
from core import collapsed,covariance,coherence,witness_vector,realize

BASE=Path(__file__).resolve().parent
SOURCE=BASE.parent/'optical-information-resolution/commercial1-group-plan/design.npz'
EXPECTED='f63a8de2a1bcf075e731e6e0462e6ab4296fd816da8a71d7a07e56079f12c2ec'

def summarize(z,n,groups,v):
    ng=int(groups.max()+1);diagonal=collapsed(z,n);target=diagonal.sum(0)
    group_means=[];code_var=0.;cv_var=0.;collision_var=0.;expected_log=0.;mass=0.
    bet=1.;g=.96;theta=np.arctan2(v[1],v[0]);scale=-2*v[2]/g
    for j in range(ng):
        ix=np.flatnonzero(groups==j);model=covariance(z[ix],n[ix],True)
        m=model['mean'];s=model['covariance'];w=float(m[:4]@v)
        group_means.append(m[:4]);vv=float(v@s[:4,:4]@v);code_var+=vv
        cov=float(v@s[:4,4]);cv_var+=max(0.,vv-cov*cov/s[4,4]) if s[4,4]>0 else vv
        p=np.prod([(4-k)/4 for k in range(len(ix))]);collision_var+=(1/p-1)*w*w
        cap=n[ix].sum(0).prod()*(1+g*(scale+1/scale)/2)
        values=[]
        for codes in itertools.product(range(4),repeat=len(ix)):
            score=realize(z[ix],n[ix],codes)@v
            factor=1+bet*score/cap
            if factor<=0:values.append(-np.inf)
            else:values.append(np.log(factor))
        expected_log+=float(np.mean(values));mass+=m[4]
    group_means=np.array(group_means);w=group_means@v
    between=float(np.sum((w-w.mean())**2))
    info=[]
    for name,var,extra in [('ordinary_H4',code_var,'aggregate code-wise moments only'),
                           ('count_mass_control',cv_var,'exact current source-count-product mass plus calibration-selected control coefficient'),
                           ('collision_aware',collision_var,'collision flag known independently from source-code design'),
                           ('source_resolved_oracle',0.,'source identities resolved on scored readouts')]:
        total=var+between
        info.append(dict(method=name,conditional_code_witness_variance=var,
            calibration_between_group_scatter=between,scatter_proxy=total,
            expected_signal_to_scatter=float(w.sum()/np.sqrt(total)) if total>0 else None,
            extra_information=extra))
    np.testing.assert_allclose(group_means.sum(0),target,rtol=1e-12,atol=1e-8)
    return dict(groups=ng,sources=len(z),minimum_four_panel_source_count_product=float(n.prod(1).min()),
        total_four_panel_count_mass=float(mass),finite_calibration_diagonal_coherence=coherence(target),
        calibration_witness_signal=float(w.sum()),methods=info,
        expected_conditional_code_log_e_for_fitted_training_witness=expected_log,
        log_nominal_e_threshold=float(np.log(20)),physical_sampling_coverage_established=False,
        note='Calibration-only design diagnostics. Fitted direction and scatter are not independent confidence or validation results.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert not a.output.exists();assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
    with np.load(SOURCE) as d:
        z=d['calibration_z'].copy();n=d['calibration_counts'].copy()
        prior=d['prior_groups'].copy();existing=d['groups'].copy();sources=d['test_sources'].copy()
    assert len(z)==279 and len(np.unique(existing))==len(np.unique(prior))==139
    # One common numerical scale, never source-specific reweighting.
    unit=float(np.median(n.prod(1))**.25);z/=unit;n/=unit
    diag=collapsed(z,n);tot=diag.sum(0);assert np.all(tot[2:]>0)
    theta=float(np.arctan2(tot[1],tot[0]));scale=float(np.sqrt(tot[3]/tot[2]));v=witness_vector(.96,theta,scale)
    original=summarize(z,n,prior,v);old_design=summarize(z,n,existing,v)
    triple_ids=[np.flatnonzero(existing==j) for j in np.unique(existing) if np.sum(existing==j)==3]
    assert len(triple_ids)==1;triple=triple_ids[0];remain=np.setdiff1d(np.arange(len(z)),triple)
    both=collapsed(z[:,None]+z[None],n[:,None]+n[None]);mixed=both-diag[:,None]-diag[None]
    per_source=diag@v;global_group_mean=per_source.sum()/139
    cost=(mixed@v)**2/3+(per_source[:,None]+per_source[None]-global_group_mean)**2
    graph=nx.Graph();graph.add_nodes_from(remain.tolist())
    for j,i in enumerate(remain):
        for k in remain[j+1:]:graph.add_edge(int(i),int(k),weight=float(cost[i,k]))
    pairs=nx.algorithms.matching.min_weight_matching(graph,maxcardinality=True,weight='weight')
    groups=np.full(len(z),-1,int)
    for j,pair in enumerate(sorted(tuple(sorted(p)) for p in pairs)):groups[list(pair)]=j
    groups[triple]=len(pairs);assert np.all(groups>=0) and len(pairs)==138
    new=summarize(z,n,groups,v)
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'calibration-group-design.npz',groups=groups,baseline_groups=prior,
                        prior_cubic_groups=existing,source_ids=sources,witness=v)
    out=dict(stage='exposed_calibration_only_feasibility',input=str(SOURCE.relative_to(BASE.parents[1])),input_sha256=EXPECTED,
        approved_npz_members=['calibration_z','calibration_counts','prior_groups','groups','test_sources'],
        validation_coordinates_or_outcomes_read=False,q=.008,g0=.96,alpha=.05,
        common_count_scale=unit,training_theta=theta,training_scale=scale,
        configurations=dict(original_groups=original,prior_cubic_calibration_groups=old_design,new_quartic_calibration_groups=new),
        all_sources_and_four_panel_count_weights_retained=True,ordinary_code_seeds_selected=False,
        source_support_assumptions_not_physically_certified=True,old_physical_failure_unchanged=True)
    (a.output/'result.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
