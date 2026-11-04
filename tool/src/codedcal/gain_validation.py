"""Support and amplitude-unit checks for shared readout normalization."""
import numpy as np
from .label_evaluation import training_support


def gain_training_support(words,traces,min_on_pixels=10):
    words=np.asarray(words,float)
    if words.ndim!=2 or not np.isfinite(words).all() or np.any(words<0):
        raise ValueError('finite nonnegative codebook required')
    if not isinstance(min_on_pixels,int) or min_on_pixels<1:
        raise ValueError('positive integer support threshold required')
    report=training_support(words,traces)
    counts=np.array([len(t) for t in traces])
    on_pixels=counts@(words>0).astype(int)
    report.update(on_pixels_per_readout=on_pixels.tolist(),
                  min_on_pixels=min_on_pixels,
                  adequate=bool(report['connected_complete_support'] and np.all(on_pixels>=min_on_pixels)))
    return report


def anchored_gain_state(initial,relative,optimizer_succeeded,lower=.25,upper=4.):
    """Normalize the cosine objective's arbitrary global gain and reject extremes.

    Limits are adapter qualification, not evidence of biological correctness.
    Invalid proposals return the unchanged training percentile state.
    """
    initial=np.asarray(initial,float);relative=np.asarray(relative,float)
    if initial.ndim!=1 or not np.isfinite(initial).all() or np.any(initial<=0):
        raise ValueError('positive finite initial scales required')
    if relative.shape!=initial.shape or not 0<lower<1<upper:
        raise ValueError('gain shape or bounds invalid')
    valid=bool(optimizer_succeeded and np.isfinite(relative).all() and np.all(relative>0))
    anchored=None
    if valid:
        anchored=relative/np.exp(np.mean(np.log(relative)))
        valid=bool(np.all((anchored>=lower)&(anchored<=upper)))
    return (initial*anchored if valid else initial.copy()),dict(
        accepted=valid,anchored_relative=anchored.tolist() if anchored is not None else None,
        reason=None if valid else 'failed_optimizer_nonfinite_or_extreme_relative_gains',
        amplitude_gauge='geometric_mean_relative_gain_one',relative_bounds=[lower,upper])


def partial_gain_training_support(words,traces,min_on_pixels=10):
    """Support for a constrained partial fit; uncovered bits stay at a prior.

    All fitted bits must be connected through observed co-active codewords.
    Priors outside this subgraph are assumed, not measured calibration gains.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    report=gain_training_support(words,traces,min_on_pixels)
    words=np.asarray(words,float);active=np.array(report['on_pixels_per_readout'])>=min_on_pixels
    seen=np.array([len(x)>0 for x in traces]);on=words[seen][:,active]>0
    adjacency=on.astype(int).T@on.astype(int)
    components=int(connected_components(csr_matrix(adjacency),directed=False)[0]) if active.any() else 0
    report.update(fitted_readouts=np.flatnonzero(active).tolist(),prior_readouts=np.flatnonzero(~active).tolist(),
        fitted_incidence_components=components,partial_fit_adequate=bool(active.sum()>=3 and components==1),
        prior_readouts_are_not_identified=True)
    return report
