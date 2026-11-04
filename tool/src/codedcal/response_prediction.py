"""Training-only readout responses and paired spatial-block decisions.

These are predictive checks, not molecular truth or a temporal-budget estimator.
A response with zero code amplitude carries only background evidence. Spatial
bootstrap intervals are descriptive within an acquisition, not safety bounds.
"""
import numpy as np
from scipy.optimize import nnls


def word_background_features(stack, words, heldout):
    stack=np.asarray(stack,float);words=np.asarray(words,float)
    if stack.ndim!=3 or words.ndim!=2 or words.shape[1]!=len(stack) or not 0<=heldout<len(stack):
        raise ValueError('inconsistent prediction dimensions')
    if not np.isfinite(stack).all() or not np.isfinite(words).all() or np.any(stack<0) or np.any(words<0):
        raise ValueError('finite nonnegative prediction inputs required')
    fit=np.setdiff1d(np.arange(len(stack)),[heldout]);b=words[:,fit]
    n=len(fit);s=b.sum(1);v=(b*b).sum(1);det=n*v-s*s
    if np.any(det<=1e-12):raise ValueError('word and background not separable')
    x=stack[fit].reshape(n,-1).T;dot=x@b.T;total=x.sum(1)[:,None];xx=(x*x).sum(1)[:,None]
    amp=(n*dot-total*s)/det;bg=(total*v-dot*s)/det
    loss=xx-2*amp*dot-2*bg*total+amp*amp*v+2*amp*bg*s+n*bg*bg
    loss[(amp<0)|(bg<0)]=np.inf
    wa=np.maximum(dot/v,0);wl=xx-2*wa*dot+wa*wa*v;use=wl<loss
    loss=np.where(use,wl,loss);amp=np.where(use,wa,amp);bg=np.where(use,0,bg)
    ba=np.maximum(total/n,0);bl=xx-2*ba*total+n*ba*ba;use=bl<loss
    loss=np.where(use,bl,loss);amp=np.where(use,0,amp);bg=np.where(use,ba,bg)
    chosen=loss.argmin(1);i=np.arange(len(x));amplitude=amp[i,chosen];background=bg[i,chosen];on=words[chosen,heldout]
    return np.column_stack([amplitude*on,amplitude*(1-on),background,np.ones(len(x))]),chosen


def fit_responses(training, words, heldout, stride=2):
    if not isinstance(stride,int) or stride<1:raise ValueError('positive integer stride required')
    small=np.asarray(training)[:,::stride,::stride];result={}
    for bit in heldout:
        features,_=word_background_features(small,words,bit);target=small[bit].ravel();coef,residual=nnls(features,target)
        result[str(bit)]=dict(coefficients=coef.tolist(),training_sse=float(residual**2),training_energy=float(target@target),
                             code_sensitive=bool(abs(coef[0]-coef[1])>1e-8))
    return dict(responses=result,training_shape=list(small.shape),feature_order=['on_word_amplitude','off_word_amplitude','shared_background','constant'])


def predict_response(stack,words,bit,state):
    features,choices=word_background_features(stack,words,bit)
    coefficients=np.asarray(state['responses'][str(bit)]['coefficients'],float)
    if coefficients.shape!=(4,) or not np.isfinite(coefficients).all() or np.any(coefficients<0):raise ValueError('invalid response state')
    return (features@coefficients).reshape(np.asarray(stack).shape[1:])


def choose_paired_blocks(losses,baseline='iterative',alpha=.05,resamples=4000,seed=104,minimum_blocks=8):
    names=list(losses);arrays={k:np.asarray(v,float) for k,v in losses.items()}
    if baseline not in arrays or not 0<alpha<1 or resamples<100:raise ValueError('invalid comparison setup')
    n=len(arrays[baseline])
    if n<minimum_blocks or any(v.shape!=(n,) or not np.isfinite(v).all() or np.any(v<0) for v in arrays.values()):
        return dict(choice=baseline,abstained=True,reason='unsupported_common_spatial_blocks')
    rng=np.random.default_rng(seed);indices=rng.integers(0,n,size=(resamples,n));eligible=[baseline];rows=[]
    for name in names:
        if name==baseline:continue
        difference=arrays[baseline]-arrays[name];distribution=difference[indices].mean(1)
        lower=float(np.quantile(distribution,alpha/(len(names)-1)));gain=float(difference.mean())
        passes=lower>0;rows.append(dict(arm=name,mean_loss_reduction=gain,simultaneous_lower_reduction=lower,passes=bool(passes)))
        if passes:eligible.append(name)
    best=min(eligible,key=lambda name:(arrays[name].mean(),names.index(name)))
    return dict(choice=best,abstained=best==baseline,reason='paired_spatial_block_prediction',rows=rows,blocks=n,alpha=alpha,
                interpretation='Descriptive spatial bootstrap within exposed acquisition; not calibrated molecular accuracy or independent confirmation')
