"""Identity/composition sensitivity and real kernel ambiguity controls."""
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import local_trace as L
import quality_catalog as Q

def main():
    L.guard();codes,negative,names,words=A.book();rng=np.random.default_rng(61009);genes=rng.choice(998,60,replace=False)
    centers=np.array([(y,x) for y in range(29,100,10) for x in range(29,100,10)])[:60];depth=np.array([6,14,22])[np.arange(60)%3]
    images=np.zeros((24,30,128,128),np.float32)
    for gene,(y,x),z in zip(genes,centers,depth):
        for r,c in enumerate(words[gene]):images[r*4+c,z,y,x]=1
    images=gaussian_filter(images,(0,.7,1,1));images/=images.max();support=np.ones_like(images,np.uint8)
    reference={(int(x),int(y),int(z)):int(j) for (y,x),z,j in zip(centers,depth,genes)};records=[]
    for decoder in A.DECODERS:
        xyz,calls,scores,eligible,runtime=Q.scored(images,support,codes,decoder)
        correct=sum(reference.get(tuple(p))==int(j) for p,j,e in zip(xyz,calls,eligible) if e)
        assert correct==60 and eligible.sum()==60 and not negative[calls[eligible]].any()
        full,selected=Q.curve(scores,calls,eligible,negative);assert selected is not None and selected['coding']==60
        records.append(dict(decoder=decoder,identity_f1=1.,coding=60,negative=0,runtime=runtime))
    # Count-preserving label corruption must defeat identity even while gene
    # composition, yield and negative controls remain unchanged.
    corrupted=np.roll(genes,1);identity=float(np.mean(genes==corrupted));tv=float(np.abs(np.bincount(genes,minlength=998)-np.bincount(corrupted,minlength=998)).sum()/120)
    assert identity==0 and tv==0 and len(corrupted)==60 and not negative[corrupted].any()
    shifts=np.tile([0,3,-2],(24,1));observed=np.stack([A.integer_warp(im,t) for im,t in zip(images,shifts)]);oracle=np.stack([A.integer_warp(im,-t) for im,t in zip(observed,shifts)])
    assert np.array_equal(oracle[:,:,16:112,16:112],images[:,:,16:112,16:112])
    ambiguous=np.zeros_like(images)
    for frame in range(24):ambiguous[frame,14,64,64]=1
    ambiguous=gaussian_filter(ambiguous,(0,.7,1,1));ambiguous/=ambiguous.max()
    ambiguity=[]
    for decoder in A.DECODERS:
        _,_,_,eligible,_=Q.scored(ambiguous,support,codes,decoder);assert not eligible.any(),decoder
        ambiguity.append(dict(decoder=decoder,accepted=0,reason='symmetric equal-color ambiguity'))
    A.write(Q.DEST/'controls.json',dict(passed=True,synthetic=records,ambiguity=ambiguity,seed=61009,count_preserving_corruption=dict(identity_f1=identity,composition_tv=tv,yield_unchanged=True,negative_proxy_unchanged=True),geometric_oracle_exact=True,supported_object_target=True,development=True,independent_confirmation=False,script_sha256=A.sha(__file__),quality_script_sha256=A.sha(Q.__file__)))
    print('PASSED three quality kernels, ambiguity abstentions and endpoint/oracle controls',flush=True)

if __name__=='__main__':main()
