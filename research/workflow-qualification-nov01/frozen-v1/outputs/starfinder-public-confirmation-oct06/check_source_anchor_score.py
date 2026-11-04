"""Equal-object angular score gains/support and coordinate controls."""
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import tissue as T
import probe as P
import source_anchor_score as O

def main():
    codes,words=T.book();rng=np.random.default_rng(61014);image=np.zeros((16,30,128,128),np.float32);points=[]
    centers=[(x,y,z) for z in [8,16,22] for y in range(37,93,9) for x in range(37,93,9)][:60];genes=rng.choice(64,60,replace=False)
    for j,(x,y,z) in zip(genes,centers):image[:,z,y,x]+=codes[:,j];points.append((words[j,0],x,y,z))
    image=gaussian_filter(image,(0,.7,1,1));image/=image.max();support=np.ones_like(image,np.uint8);points=np.array(points);zero=np.zeros((16,3),int)
    native=O.score(image,support,points,zero,codes);assert native['loss']<1e-6
    gains=np.repeat([1,.2,5,2],4);changed=image*gains[:,None,None,None];gain=O.score(changed,support,points,zero,codes);assert abs(gain['loss']-native['loss'])<1e-7
    imposed=T.offsets('color');observed=P.warp(image,imposed);alpha=P.warp(support,imposed);rawpoints=points.copy();rawpoints[:,1:3]+=imposed[points[:,0]][:,::-1];inverse=np.c_[np.zeros(16,int),-imposed];correct=O.score(observed,alpha,rawpoints,inverse,codes);wrong=O.score(observed,alpha,rawpoints,zero,codes);assert correct['loss']<1e-6 and wrong['loss']>correct['loss']+.05
    masked=np.zeros_like(support);failed=O.score(image,masked,points,zero,codes);assert failed['loss']==1 and failed['unsupported_or_missing']==60
    flow=[np.zeros(image.shape[1:]+(3,)) for _ in range(4)];identity=O.score(image,support,points,zero,codes,flow);assert identity['object_losses']==native['object_losses']
    mislabeled=points.copy();mislabeled[:,0]=(mislabeled[:,0]+1)%4;bad=O.score(image,support,mislabeled,zero,codes);assert bad['loss']>native['loss']+.05
    A.write(O.DEST/'controls.json',dict(passed=True,native_loss=native['loss'],gain_invariant_loss=gain['loss'],corrected_loss=correct['loss'],wrong_loss=wrong['loss'],missing_loss=failed['loss'],zero_flow_identical=True,incorrect_source_channel_loss=bad['loss'],seed=61014,script_sha256=A.sha(O.__file__),control_script_sha256=A.sha(__file__),protocol_sha256=A.sha(A.HERE/'SOURCE_ANCHOR_PROTOCOL.md'),development=True))
    print('PASS tracked source coordinates, round gains, inverse sampling and support controls',flush=True)
if __name__=='__main__':main()
