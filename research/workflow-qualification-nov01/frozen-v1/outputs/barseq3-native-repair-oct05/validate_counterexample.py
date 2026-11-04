"""Validate a preselected algebraic witness; no empirical prevalence estimate."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL','2')
os.environ.setdefault('TF_NUM_INTRAOP_THREADS','1')
os.environ.setdefault('TF_NUM_INTEROP_THREADS','1')
import sys
import numpy as np
from condition import HERE,ROOT,read,write,sha
sys.path.insert(0,str(ROOT.parent/'bardensr-upstream.nosync'))
import bardensr

spec=read(HERE/'counterexample-specification.json')
classes=read(HERE/'codebook-classes.json')
B=np.zeros((40,len(classes['sequences'])),float)
for j,word in enumerate(classes['sequences']):
    for r,c in enumerate(word):B[4*r+'GTCA'.index(c),j]=1
blank=np.array(classes['blank'],bool)
signals=np.array([spec['complete_point_signal'],spec['isolated_point_signal']],float).T
expected=(signals.T@B)/((np.linalg.norm(signals,axis=0)[:,None]+.05)*np.sqrt(10))
density=np.asarray(bardensr.spot_calling.estimate_density_singleshot(signals[:,None,None,:],B,.05)).reshape(2,-1)
np.testing.assert_allclose(density,expected,rtol=1e-6,atol=1e-7)
parents=spec['parent_indices'];target=spec['absent_target_index']
support=signals[:,0]>0
permitted=np.flatnonzero(~np.any(B[~support]>0,axis=0))
assert set(permitted)==set(parents)
assert np.linalg.matrix_rank(B[:,parents])==2
weights=np.linalg.lstsq(B[:,parents],signals[:,0],rcond=None)[0]
np.testing.assert_allclose(weights,[.5,.5],atol=1e-12)
np.testing.assert_allclose(B[:,parents]@weights,signals[:,0],atol=1e-12)
records={}
for i,mode in enumerate(['complete','one_camera_locus']):
    calls=np.flatnonzero(density[i]>.775)
    records[mode]=dict(called_indices=calls.tolist(),called_names=[classes['names'][j] for j in calls],
        blank_calls=int(blank[calls].sum()),absent_target_called=bool(target in calls),
        absent_target_evidence=float(density[i,target]),parent_evidence=density[i,parents].tolist(),
        known_false_nonblank_calls=int(sum(j not in parents and not blank[j] for j in calls)),
        known_true_parent_calls=int(sum(j in parents for j in calls)))
assert records['one_camera_locus']['absent_target_called']
out=dict(status='passed',cases=records,kernel_matches_algebra=True,
    complete_known_gain_nonnegative_solution_unique=True,complete_solution_weights=weights.tolist(),
    general_blank_support_bound_does_not_cover_all_blank_codes=not spec['all_blanks_below_correlation_cut_for_any_positive_camera_only_signal'],
    scope='Selected noiseless point-vector counterexample from the public codebook, not a spatial peak benchmark, empirical prevalence estimate, true-FDR calibration or prospective molecular validation.',
    fixed_gain_scope='Unique complete nonnegative mixture is an algebraic oracle; it is not an execution or success claim for the gain-learning regression kernel.',
    specification_sha256=sha(HERE/'counterexample-specification.json'),script_sha256=sha(__file__))
write(HERE/'counterexample-validation.json',out)
print(records,flush=True)
