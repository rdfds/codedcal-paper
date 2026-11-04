"""Objective diagnosis on tracked objects in unchanged source-conditioned inputs."""
import gc,json,time
import numpy as np
from scipy.ndimage import label,center_of_mass
from skimage.morphology import local_maxima
import adapter as A
import local_trace as L
import tissue as T
import fullfield_tissue as F
import probe as P
DEST=A.HERE/'source-anchor-score'
ORDER=['none','common_iterative','template_2d_multistart','template_3d_multistart','author_multistart','source_round_3d']

def anchors(images,support):
    out=[]
    for c in range(4):
        peaks=local_maxima(images[c],connectivity=3,allow_borders=True)&(images[c]>.2);tags,n=label(peaks,structure=np.ones((3,3,3)))
        for center in center_of_mass(peaks,tags,np.arange(1,n+1)):
            z,y,x=np.floor(np.array(center)+.5).astype(int)
            if 6<=z<24 and 32<=y<96 and 32<=x<96 and np.all(support[:,z-1:z+2,y-2:y+3,x-2:x+3]):out.append((c,x,y,z))
    return np.asarray(out,int).reshape(-1,4)

def score(images,support,points,maps,book,flows=None):
    rounds=len(images)//4;normbook=book/np.sqrt((book*book).sum(axis=0,keepdims=True));losses=[];missing=0;duplicates=0
    offsets=np.array([(x,y,z) for z in [-1,0,1] for y in [-2,-1,0,1,2] for x in [-2,-1,0,1,2]],int)
    for c,x,y,z in points:
        target=np.array([x,y,z])+maps[c][::-1];grid=target[None]+offsets;values=[];valid=True
        for f in range(len(images)):
            pull=grid.astype(float)
            if flows is not None:
                if np.any(grid<0) or np.any(grid>=np.array(images.shape[1:])[::-1]):valid=False;break
                pull+=flows[f//4][grid[:,2],grid[:,1],grid[:,0]]
            source=np.rint(pull-maps[f][::-1]).astype(int);unique=np.unique(source,axis=0);duplicates+=len(source)-len(unique)
            if np.any(unique<0) or np.any(unique>=np.array(images.shape[1:])[::-1]) or not np.all(support[f,unique[:,2],unique[:,1],unique[:,0]]):valid=False;break
            values.append(float(images[f,unique[:,2],unique[:,1],unique[:,0]].sum(dtype=np.float64)))
        if not valid:losses.append(1.);missing+=1;continue
        values=np.asarray(values).reshape(rounds,4);energy=np.sqrt((values*values).sum(axis=1,keepdims=True));normalized=values/(energy+1e-6);signal=float((normalized*normalized).sum())
        if signal<rounds-1e-4:losses.append(1.);missing+=1;continue
        projection=normbook[:,book[c]>0].T@normalized.ravel();losses.append(float(np.clip(1-projection.max()**2/signal,0,1)))
    return dict(loss=float(np.mean(losses)) if losses else None,object_losses=losses,anchors=len(points),unsupported_or_missing=missing,duplicate_inverse_voxels=duplicates)

def flows(observed,alpha,record):
    maps=np.asarray(record['maps_zyx'],int);out=T.move3d(observed,maps);fixed=out[:4].max(axis=0);result=[]
    for r in range(4):
        field,diag=L.fit(fixed,out[r*4:r*4+4].max(axis=0));original=record['local'][r]
        if 'field_sha256' in original:assert diag['field_sha256']==original['field_sha256']
        else:assert original['identity'] and not np.any(field)
        result.append(field)
    assert not np.any(result[0]);return maps,result

def main():
    T.guard();assert json.loads((DEST/'controls.json').read_text())['passed'];codes,_=T.book()
    with np.load(T.CACHE/'fullfield-conditioned-crop-4.npz') as data:base=data['images'];support=data['support']
    decisions=[]
    for kind in ['zero','common','color','outclass']:
        geom=json.loads((F.DEST/f'{kind}-geometry.json').read_text());imposed=T.offsets(kind);observed=P.warp(base,imposed);alpha=P.warp(support,imposed);points=anchors(observed,alpha);scores={}
        for policy in ORDER+['global_group']:
            if policy=='source_round_3d':maps,field=flows(observed,alpha,geom['diagnostics'][policy])
            else:maps=np.asarray(geom['maps_zyx'][policy],int);field=None
            scores[policy]=score(observed,alpha,points,maps,codes,field)
        valid=[p for p in ORDER if scores[p]['loss'] is not None];best=min(valid,key=lambda p:scores[p]['loss']);eligible=[p for p in valid if scores[p]['loss']<=scores[best]['loss']+.005];selected=min(eligible,key=ORDER.index)
        if scores['global_group']['loss'] is not None and scores[best]['loss']-scores['global_group']['loss']>=.005:selected='global_group'
        if len(points)<10:selected=None
        delta=np.asarray(scores[best]['object_losses'])-np.asarray(scores['global_group']['object_losses'])
        if len(delta):
            rng=np.random.default_rng(61015);replicate=delta[rng.integers(0,len(delta),(1000,len(delta)))].mean(axis=1);uncertainty=dict(mean=float(delta.mean()),paired_object_bootstrap95=np.quantile(replicate,[.025,.975]).tolist(),iid_object_assumption_not_verified=True)
        else:uncertainty=None
        metrics={}
        for dec in A.DECODERS:
            if selected is not None:
                record=json.loads((F.DEST/'outcomes'/f'{kind}__{selected}__{dec}.json').read_text());metrics[dec]=record['metrics']['3']['identity_f1']
            else:metrics[dec]=None
        out=dict(kind=kind,selected=selected,best_conventional_loss_policy=best,scores=scores,abstained_insufficient_objects=selected is None,group_advantage_uncertainty=uncertainty,selected_F1_for_evaluation_only=metrics,development=True,independent_confirmation=False,score_script_sha256=A.sha(__file__));A.write(DEST/f'{kind}.json',out);decisions.append(out);print('SOURCE ANCHOR SCORE',kind,selected,{p:s['loss'] for p,s in scores.items()},metrics,flush=True);gc.collect()
    A.write(DEST/'state.json',dict(stage='complete',cases=4,development=True,independent_confirmation=False,general_recommendation_established=False,decisions=decisions))
if __name__=='__main__':main()
