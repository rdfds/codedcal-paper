"""Training-only quality filtering using codebook blanks.

Area, measured intensity and fitted coefficient follow the JSIT calling-stage
features. This adapter uses finite quantile bins, fixed pseudocounts and rejects
unseen/out-of-range bins rather than assigning them zero error. The estimate
is a blank-derived proxy, not molecular truth or a calibrated FDR bound.
"""
import numpy as np


def fit_blank_filter(features,word_indices,blank_mask,bins=4,pseudocount=.5):
    x=np.asarray(features,float);w=np.asarray(word_indices);blank=np.asarray(blank_mask,bool)
    if x.ndim!=2 or x.shape[1]!=3 or w.ndim!=1 or not np.issubdtype(w.dtype,np.integer) or len(w)!=len(x) or np.any(x<0) or not np.isfinite(x).all() or not isinstance(bins,int) or bins<2 or not np.isfinite(pseudocount) or pseudocount<=0:raise ValueError('invalid quality inputs')
    if blank.ndim!=1 or not blank.any() or blank.all() or np.any(w<0) or np.any(w>=len(blank)):raise ValueError('invalid codebook labels')
    if not len(x):return dict(qualified=False,reason='no_quality_training_candidates')
    log=np.log1p(x);edges=[np.unique(np.quantile(log[:,j],np.linspace(0,1,bins+1))) for j in range(3)]
    if any(len(e)<3 for e in edges):return dict(qualified=False,reason='degenerate_quality_dimensions')
    shape=tuple(len(e)-1 for e in edges);index=[]
    for j,e in enumerate(edges):index.append(np.clip(np.searchsorted(e,log[:,j],side='right')-1,0,len(e)-2))
    key=np.ravel_multi_index(tuple(index),shape);n=np.prod(shape);g=np.bincount(key[~blank[w]],minlength=n);b=np.bincount(key[blank[w]],minlength=n)
    proxy=(b+pseudocount)/(g+pseudocount)*int((~blank).sum())/int(blank.sum())
    return dict(qualified=True,edges=[e.tolist() for e in edges],shape=list(shape),gene_counts=g.tolist(),blank_counts=b.tolist(),proxy=proxy.tolist(),pseudocount=pseudocount,
       inputs='area, measured intensity norm, mean fitted source coefficient; no outcome/reference labels')


def blank_filter_scores(features,state):
    x=np.asarray(features,float)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all() or np.any(x<0):raise ValueError('invalid features')
    if not state['qualified']:return np.full(len(x),np.inf)
    log=np.log1p(x);index=[];valid=np.ones(len(x),bool)
    for j,e in enumerate(state['edges']):
        e=np.array(e);valid&=(log[:,j]>=e[0])&(log[:,j]<=e[-1]);index.append(np.clip(np.searchsorted(e,log[:,j],side='right')-1,0,len(e)-2))
    key=np.ravel_multi_index(tuple(index),tuple(state['shape']));scores=np.array(state['proxy'])[key];counts=np.array(state['gene_counts'])+np.array(state['blank_counts']);scores[(~valid)|(counts[key]==0)]=np.inf
    return scores
