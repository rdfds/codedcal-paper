import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import tissue as T
import volumetric_model as M

def main():
    codes,words=T.book();rng=np.random.default_rng(61016);image=np.zeros((16,8,24,24),np.float32)
    for j in rng.choice(64,30,replace=False):
        z,y,x=rng.integers([2,6,6],[6,18,18]);image[:,z,y,x]+=codes[:,j]
    image=gaussian_filter(image,(0,.5,.7,.7));m=np.tile([[0,1,0],[0,0,1],[0,-1,0],[0,0,-1]],(4,1));obs=T.move3d(image,m);inv=-m;native=M.loss(image,codes,np.zeros((16,3),int));correct=M.loss(obs,codes,inv);wrong=M.loss(obs,codes,np.zeros((16,3),int));assert np.allclose(T.move3d(obs,inv),image,atol=1e-5);assert abs(correct-native)<1e-5 and wrong>correct+.01
    model=A.R.qualify.bardensr.registration.find_translations_using_model(obs.astype(float),codes,maximum_wiggle=(3,10,10),niter=3,initial_guess=inv.astype(float));assert model.corrections.shape==(16,3) and np.isfinite(model.losses).all()
    A.write(M.DEST/'controls.json',dict(passed=True,native_loss=native,corrected_loss=correct,wrong_loss=wrong,actual_author_3d_finite=True,seed=61016,script_sha256=A.sha(M.__file__),control_sha256=A.sha(__file__),protocol_sha256=A.sha(A.HERE/'VOLUMETRIC_MODEL_PROTOCOL.md')));print('PASS full-volume author API, XYZ inverse and loss controls',flush=True)
if __name__=='__main__':main()
