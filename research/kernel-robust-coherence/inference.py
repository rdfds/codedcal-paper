"""Prespecified finite-sample witnesses and explicitly approximate diagnostics."""
import numpy as np
from scipy.special import logsumexp
from scipy.stats import t as student_t
from core import witness_vector


def _inputs(components,capacities):
    x=np.asarray(components,float);b=np.asarray(capacities,float)
    if x.ndim!=2 or x.shape[1]!=4 or b.shape!=(len(x),) or len(x)<2:raise ValueError('at least two groups with four components required')
    if not np.isfinite(x).all() or not np.isfinite(b).all() or np.any(b<0):raise ValueError('finite capacities required')
    if np.any(np.max(abs(x),axis=1)>b+1e-8*np.maximum(1,b)):raise ValueError('component exceeds supplied deterministic range')
    if not np.any(b>0):raise ValueError('no count support')
    return x,b


def hoeffding_rectangle(components,capacities,alpha=.05):
    x,b=_inputs(components,capacities)
    if not 0<alpha<1:raise ValueError('alpha')
    scale=b.sum();mu=x.sum(0)/scale
    radius=float(np.sqrt(2*np.sum((b/scale)**2)*np.log(8/alpha)))
    lo=mu-radius;hi=mu+radius
    # Requiring positive lower references is an explicit numerical-support gate.
    usable=bool(np.all(lo[2:]>0))
    dist=np.maximum(np.maximum(lo[:2],-hi[:2]),0)
    lower=float(np.linalg.norm(dist)/np.sqrt(hi[2]*hi[3])) if usable else 0.
    return dict(method='simultaneous finite bounded-group Hoeffding rectangle',alpha=alpha,
        component_mean=mu.tolist(),component_lower=lo.tolist(),component_upper=hi.tolist(),
        radius=radius,positive_reference_support=usable,coherence_lower=lower,
        note='Valid only for independent bounded groups under the population model; not physical assumption validation.')


DEFAULT_PHASES=tuple(np.arange(16)*np.pi/8)
DEFAULT_SCALES=(.5,2**-.5,1.,2**.5,2.)


def witness_test(components,capacities,g=.96,alpha=.05,
                 phases=DEFAULT_PHASES,scales=DEFAULT_SCALES,bets=(.25,.5,1.)):
    x,b=_inputs(components,capacities)
    if not 0<alpha<1 or not 0<=g<=1:raise ValueError('threshold/alpha')
    phases=tuple(phases);scales=tuple(scales);bets=tuple(bets)
    if not phases or not scales or not bets or any(not 0<l<=1 for l in bets):raise ValueError('fixed nonempty mixture required')
    logs=[]
    for theta in phases:
        for t in scales:
            v=witness_vector(g,theta,t)
            L=b*(1+g*(t+1/t)/2)
            w=x@v
            u=np.divide(w,L,out=np.zeros_like(w),where=L>0)
            if np.any(abs(u)>1+1e-10):raise ValueError('witness range violation')
            for lam in bets:
                factors=1+lam*u
                if np.any(factors<0):raise ValueError('negative betting factor')
                with np.errstate(divide='ignore'):logs.append(float(np.log(factors).sum()))
    log_e=float(logsumexp(logs)-np.log(len(logs)))
    return dict(method='fixed finite mixture of bounded-group mean-null e-values',g=g,alpha=alpha,
        log_evalue=log_e,log_rejection_threshold=float(np.log(1/alpha)),
        reject_attenuation_at_most_g=bool(log_e>=np.log(1/alpha)),mixture_members=len(logs),
        usable_groups=int((b>0).sum()),
        note='Model-conditional threshold test. Independent group null means and deterministic capacities required; no same-data witness selection.')


def t_witness_diagnostic(components,g=.96,theta=0.,scale=1.,alpha=.05):
    """Approximate, not finite-sample certified, not a physical validation claim."""
    x=np.asarray(components,float)
    if x.ndim!=2 or x.shape[1]!=4 or len(x)<2 or not np.isfinite(x).all():raise ValueError('components')
    w=x@witness_vector(g,theta,scale);se=float(w.std(ddof=1)/np.sqrt(len(w)));mean=float(w.mean())
    lower=mean-float(student_t.ppf(1-alpha,len(w)-1))*se
    return dict(method='approximate one-sided student-t linear-witness diagnostic',g=g,
        theta=theta,scale=scale,alpha=alpha,witness_mean=mean,witness_SE=se,witness_lower=lower,
        nominal_rejection=bool(lower>0),finite_sample_guarantee=False,
        note='If theta/scale were fitted on these reads this is a training diagnostic only.')
