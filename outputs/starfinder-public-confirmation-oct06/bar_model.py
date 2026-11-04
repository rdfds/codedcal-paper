"""Strengthened conventional barcode calibration; no group-search initialization."""
import gc,json,subprocess,time,traceback
import numpy as np
import adapter as A
import initialization as I
import recovery as R
import global_repair as G
from bardensr.registration import lowrankregistration as LR
DEST=A.HERE/'bar-model'

def main():
    A.guard()
    for name in ['bar_model.py','BAR_MODEL_PROTOCOL.md']:
        p=A.HERE/name;assert subprocess.check_output(['git','show','HEAD:'+str(p.relative_to(A.ROOT))],cwd=A.ROOT)==p.read_bytes()
    book,blank,labels=R.qualify.codebook();names=list(dict.fromkeys(labels));gene_ids=np.array([names.index(n) for n in labels]);outcomes=[]
    for k in [0,1]:
        images,support=R.inputs(k);refs={d:json.loads((R.HERE/'recovery'/f'field-{k}'/f'reference-{d}.json').read_text()) for d in R.DECODERS}
        for kind in ['zero','common','color']:
            imposed=R.offsets(kind);observed=R.warp(images,imposed);alpha=R.warp(support,imposed);start=time.monotonic();template,diag=I.fit(observed);seeds=[('zero',np.zeros((40,2),int)),('common',G.common_fit(observed,10))]+[(r['start'],np.asarray(r['maps'],int)) for r in diag['starts']];seen=set();candidates=[];models=[]
            for name,guess in seeds:
                token=guess.tobytes()
                if token in seen:continue
                seen.add(token);initial_loss=float(LR.calc_loss(observed.astype(float),book,-20-guess,zero_padding=20));candidates.append((initial_loss,float(np.square(guess).sum()),name+'_initial',guess.copy()))
                try:
                    model=R.qualify.bardensr.registration.find_translations_using_model(observed[:,None].astype(float),book,maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(40),-guess].astype(float));maps=np.rint(-np.asarray(model.corrections)[:,1:]).astype(int);losses=np.asarray(model.losses);assert np.isfinite(losses).all();loss=float(LR.calc_loss(observed.astype(float),book,-20-maps,zero_padding=20));candidates.append((loss,float(np.square(maps).sum()),name+'_refined',maps));models.append(dict(start=name,status='completed',initial_loss=initial_loss,rounded_loss=loss,losses=losses.tolist(),maps=maps.tolist()))
                except Exception:models.append(dict(start=name,status='failed',traceback=traceback.format_exc()))
            loss,_,selected,maps=min(candidates,key=lambda t:t[:3]);geometry=dict(field=k,kind=kind,selected=selected,coding_loss=loss,maps=maps.tolist(),models=models,candidates=[dict(loss=l,name=n,maps=m.tolist()) for l,_,n,m in candidates],fit_s=time.monotonic()-start,unique_starts=len(seen));A.write(DEST/f'{k}_{kind}-geometry.json',geometry)
            corrected=R.warp(observed,maps);oa=R.warp(alpha,maps);exact=np.array_equal(corrected[:,16:112,16:112],images[:,16:112,16:112]);fraction=float(np.all(oa[:,24:104,24:104]>=1-1e-6,axis=0).mean())
            for dec in R.DECODERS:
                ref=refs[dec]
                if exact:xy=np.array(ref['xy']).reshape(-1,2);words=np.array(ref['words'],int);runtime=dict(exact_native_input_call_reuse=True,kernel_s=0.)
                else:xy,words,runtime=R.decode(corrected,book,dec)
                metrics={str(r):R.metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,gene_ids,names,r) for r in [1.5,3,5]};row=dict(field=k,kind=kind,decoder=dec,status='completed',selected=selected,xy=xy.tolist(),words=words.tolist(),metrics=metrics,quality=R.quality(words,blank,fraction),runtime=runtime,geometry_sha256=A.sha(DEST/f'{k}_{kind}-geometry.json'),script_sha256=A.sha(__file__),development=True,independent_confirmation=False);A.write(DEST/f'{k}_{kind}_{dec}.json',row);outcomes.append(row);print('STRONG BAR MODEL',k,kind,dec,selected,row['metrics']['3']['gene_class_identity_f1'],flush=True);gc.collect()
            A.write(DEST/'state.json',dict(stage='complete' if len(outcomes)==12 else 'partial',planned=12,completed=len(outcomes),development=True,independent_confirmation=False))
if __name__=='__main__':main()
