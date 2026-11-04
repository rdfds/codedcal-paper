"""Focused unmodified full-volume barcode comparator on exposed tissue."""
import gc,json,subprocess,time,traceback
import numpy as np
import tissue as T
import adapter as A
import probe as P
from bardensr.registration import lowrankregistration as LR
DEST=A.HERE/'volumetric-model'; PAD=np.array([3,10,10])
def loss(image,codes,m):return float(LR.calc_loss(image.astype(float),codes,-PAD-m,zero_padding=tuple(PAD)))
def main():
    T.guard()
    for name in ['volumetric_model.py','VOLUMETRIC_MODEL_PROTOCOL.md','volumetric-model/controls.json']:
        path=A.HERE/name;assert subprocess.check_output(['git','show','HEAD:'+str(path.relative_to(A.ROOT))],cwd=A.ROOT)==path.read_bytes()
    ctrl=json.loads((DEST/'controls.json').read_text());assert ctrl['passed'] and ctrl['script_sha256']==A.sha(__file__)
    codes,_=T.book();data=[]
    for k in [0,4]:
        path=T.CACHE/f'fullfield-conditioned-crop-{k}.npz';assert A.sha(path)==json.loads((A.HERE/'fullfield-tissue'/f'crop-{k}.json').read_text())['array_sha256']
        with np.load(path) as d:data.append((d['images'],d['support']))
    imposed=T.offsets('color');cal=P.warp(data[0][0],imposed);obs=P.warp(data[1][0],imposed);alpha=P.warp(data[1][1],imposed)
    geom=json.loads((A.HERE/'fullfield-tissue/color-geometry.json').read_text());seeds=[('zero',np.zeros((16,3),int))]+[(p,np.asarray(geom['maps_zyx'][p],int)) for p in ['template_3d_multistart','author_multistart']];candidates=[];models=[];start=time.monotonic()
    for name,m in seeds:
        initial=loss(cal,codes,m);candidates.append((initial,float(np.square(m).sum()),name+'_initial',m));tick=time.monotonic();print('VOLUMETRIC START',name,initial,flush=True)
        try:
            model=A.R.qualify.bardensr.registration.find_translations_using_model(cal.astype(float),codes,maximum_wiggle=tuple(PAD),niter=50,initial_guess=-m.astype(float));maps=np.rint(-model.corrections).astype(int);ls=np.asarray(model.losses);assert np.isfinite(ls).all() and np.isfinite(maps).all();value=loss(cal,codes,maps);candidates.append((value,float(np.square(maps).sum()),name+'_refined',maps));models.append(dict(name=name,status='completed',initial_loss=initial,rounded_loss=value,losses=ls.tolist(),maps=maps.tolist(),fit_s=time.monotonic()-tick));print('VOLUMETRIC DONE',name,value,maps.tolist(),flush=True)
        except Exception:models.append(dict(name=name,status='failed',traceback=traceback.format_exc(),fit_s=time.monotonic()-tick))
        A.write(DEST/'fit-progress.json',dict(models=models,planned=3,completed=len(models),development=True));gc.collect()
    value,_,selected,maps=min(candidates,key=lambda t:t[:3]);geometry=dict(selected=selected,loss=value,maps=maps.tolist(),models=models,fit_s=time.monotonic()-start,candidates=[dict(name=n,loss=l,maps=m.tolist()) for l,_,n,m in candidates],development=True,independent_confirmation=False);A.write(DEST/'geometry.json',geometry)
    out=T.move3d(obs,maps);support=T.move3d(alpha,maps);rows=[]
    for dec in A.DECODERS:
        ref=json.loads((A.HERE/'fullfield-tissue/native'/f'4__{dec}.json').read_text())['catalog'];tick=time.monotonic();pred=T.catalog(out,support,codes,dec);metrics={str(r):P.metrics(ref,pred,r) for r in [1.5,3,5]};row=dict(kind='color',policy='author_full_volume',decoder=dec,status='completed',selected=selected,prediction=pred,metrics=metrics,decode_s=time.monotonic()-tick,policy_geometry_s=geometry['fit_s'],development=True,independent_confirmation=False);A.write(DEST/f'color__{dec}.json',row);rows.append(row);print('VOLUMETRIC SCORE',dec,metrics['3']['identity_f1'],flush=True)
    A.write(DEST/'state.json',dict(stage='complete',completed=3,planned=3,failed_starts=sum(m['status']=='failed' for m in models),development=True,independent_confirmation=False))
if __name__=='__main__':main()
