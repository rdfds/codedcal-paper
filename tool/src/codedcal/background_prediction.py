"""Readout-fold prediction with an explicit shared-background nuisance."""
import numpy as np


def single_word_background_prediction(stack,codebook,heldout):
    """NNLS word amplitude plus common background from non-held readouts.

    Unlike a cosine-only fit, this predicts the full specified intensity model.
    It cannot distinguish correlated background from optical truth and does not
    claim source identity. All choices use fitting readouts only.
    """
    from .readout_validation import conditional_reconstruction
    stack=np.asarray(stack,float);book=np.asarray(codebook,float);heldout=np.asarray(heldout,int)
    _,diagnostic=conditional_reconstruction(stack,book,heldout)
    training=np.asarray(diagnostic['training_readouts']);b=book[:,training];x=stack[training].reshape(len(training),-1).T
    n=len(training);s=b.sum(1);v=(b*b).sum(1);det=n*v-s*s
    valid=det>1e-10
    if not valid.any():raise ValueError('all word amplitudes confounded with background')
    dot=x@b.T;total=x.sum(1)[:,None];xx=(x*x).sum(1)[:,None]
    denom=np.where(valid,det,1.)[None]
    amplitude=(n*dot-total*s[None])/denom
    background=(total*v[None]-dot*s[None])/denom
    loss=xx-2*amplitude*dot-2*background*total+amplitude**2*v[None]+2*amplitude*background*s[None]+n*background**2
    loss[(amplitude<0)|(background<0)|(~valid[None])]=np.inf
    axis_amplitude=np.maximum(dot/np.maximum(v[None],1e-12),0)
    axis_loss=xx-2*axis_amplitude*dot+axis_amplitude**2*v[None]
    only_word=axis_loss<loss;loss=np.where(only_word,axis_loss,loss)
    amplitude=np.where(only_word,axis_amplitude,amplitude);background=np.where(only_word,0.,background)
    bg=np.maximum(total/n,0);bg_loss=xx-2*bg*total+n*bg**2
    only_background=bg_loss<loss;loss=np.where(only_background,bg_loss,loss)
    amplitude=np.where(only_background,0.,amplitude);background=np.where(only_background,bg,background)
    loss[:,~valid]=np.inf
    chosen=np.argmin(loss,axis=1);pixel=np.arange(len(x))
    predicted=amplitude[pixel,chosen,None]*book[chosen][:,heldout]+background[pixel,chosen,None]
    matrix=b/np.maximum(np.sqrt(v[:,None]),1e-12)
    _,groups,counts=np.unique(matrix,axis=0,return_inverse=True,return_counts=True)
    # Affine collinearity with the background can create further collisions.
    centered=b-b.mean(1)[:,None];length=np.linalg.norm(centered,axis=1)
    centered/=np.maximum(length[:,None],1e-12)
    _,cg,cc=np.unique(np.round(centered,12),axis=0,return_inverse=True,return_counts=True)
    supported=valid[chosen]&(counts[groups[chosen]]==1)&(cc[cg[chosen]]==1)
    diag=dict(training_readouts=training.tolist(),heldout_readouts=heldout.tolist(),
              identified_pixel_mask=supported.reshape(stack.shape[1:]),
              nuisance='per-pixel nonnegative common background, fit on training readouts only',
              background_confounded_words=int((~valid).sum()))
    return predicted.T.reshape((len(heldout),)+stack.shape[1:]),diag
