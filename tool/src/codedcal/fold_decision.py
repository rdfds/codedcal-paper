"""Observable calibration decisions with singleton readout folds.

Jointly hiding several bits can collapse otherwise distinct codewords. Each
singleton fold retains more barcode information; arbitrary source mixtures
are still not uniquely identified by this single-word/background model.
"""
import numpy as np
from scipy.ndimage import binary_dilation
from .background_prediction import single_word_background_prediction
from .readout_validation import transported_prediction_error,choose_calibration


def fold_calibration_decision(stacks,words,raw_target,maps,origin,heldout_readouts,relative_margin=.1,border=8):
    raw=np.asarray(raw_target,float);words=np.asarray(words,float)
    heldout=list(heldout_readouts)
    if not heldout or len(set(heldout))!=len(heldout) or any(int(k)!=k or k<0 or k>=len(raw) for k in heldout):
        raise ValueError('distinct valid heldout indices required')
    if set(stacks)!=set(maps) or 'none' not in stacks:raise ValueError('matching candidates and none baseline required')
    if border<1 or 2*border>=min(raw.shape[1:]):raise ValueError('invalid prediction border')
    errors={name:dict(squared_error=0.,observed_energy=0.,selected_pixels=0,supported_pixels=0) for name in stacks};folds=[]
    for bit in heldout:
        fit=np.setdiff1d(np.arange(len(raw)),[bit])
        mask=binary_dilation(np.linalg.norm(raw[fit],axis=0)>.5,iterations=6)
        mask[:border]=False;mask[-border:]=False;mask[:,:border]=False;mask[:,-border:]=False
        local={}
        for name,stack in stacks.items():
            prediction,diagnostic=single_word_background_prediction(stack,words,[bit])
            support=diagnostic.pop('identified_pixel_mask')
            error=transported_prediction_error(prediction,raw,[bit],maps[name],origin,mask,support)
            local[name]=dict(**error,model=diagnostic)
            for key in ['squared_error','observed_energy','selected_pixels','supported_pixels']:errors[name][key]+=error[key]
        folds.append(dict(heldout_readout=int(bit),candidates=local))
    # Any incomplete fold excludes the candidate; combining folds cannot dilute it.
    for name,error in errors.items():
        if any(f['candidates'][name]['supported_pixels']!=f['candidates'][name]['selected_pixels'] for f in folds):
            error['supported_pixels']=0;error['reason']='incomplete_singleton_fold'
        else:error['reason']=None
    decision=choose_calibration(errors,relative_margin=relative_margin,priority=list(stacks))
    return dict(decision=decision,errors=errors,folds=folds,
        input_contract='observed validation intensities, codebook, candidate maps, frozen disjoint gains; no injected dose or call-reference labels',
        interpretation='specified single-word/background predictive model; no source-mixture uniqueness or temporal-budget validation')
