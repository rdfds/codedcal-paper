"""Measured-vector checks on externally fitted source positions and labels.

This never generates labels or re-estimates gains. Shared photometric checks
can be applied to different spatial estimators; they do not establish RNA truth.
"""
import numpy as np
from scipy.ndimage import map_coordinates


def spectral_consistency(stack, positions, labels, codebook,
                         max_distance=.6, minimum_norm=2.):
    stack=np.asarray(stack,float);xy=np.asarray(positions,float)
    labels=np.asarray(labels);book=np.asarray(codebook,float)
    if stack.ndim!=3 or book.ndim!=2 or book.shape[1]!=len(stack):
        raise ValueError('inconsistent stack/codebook dimensions')
    if labels.ndim!=1 or not np.issubdtype(labels.dtype,np.integer) or xy.shape!=(len(labels),2):
        raise ValueError('source positions and integer labels must match')
    if any(not np.isfinite(x).all() for x in [stack,xy,book]) or np.any(stack<0) or np.any(book<0):
        raise ValueError('finite nonnegative photometric inputs required')
    if np.any(labels<0) or np.any(labels>=len(book)) or np.any(np.linalg.norm(book,axis=1)==0):
        raise ValueError('invalid source labels or unsupported codebook')
    if np.any(xy<0) or np.any(xy[:,0]>stack.shape[2]-1) or np.any(xy[:,1]>stack.shape[1]-1):
        raise ValueError('source outside measured support')
    if not np.isfinite(max_distance) or not 0<=max_distance<=2 or not np.isfinite(minimum_norm) or minimum_norm<0:
        raise ValueError('invalid photometric thresholds')
    traces=np.array([map_coordinates(channel,xy[:,::-1].T,order=1) for channel in stack]).T
    norm=np.linalg.norm(traces,axis=1)
    unit=traces/np.maximum(norm[:,None],np.finfo(float).tiny)
    word=book[labels];word=word/np.linalg.norm(word,axis=1)[:,None]
    distance=np.linalg.norm(unit-word,axis=1)
    keep=(norm>0)&(norm>=minimum_norm)&(distance<=max_distance)
    return keep,dict(distance=distance,norm=norm)
