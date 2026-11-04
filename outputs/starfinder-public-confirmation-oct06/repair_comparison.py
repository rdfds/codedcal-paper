"""Unchanged comparison parameters on a declared exposed source-adapter repair."""
import os,subprocess
from scipy.ndimage import label
from skimage.morphology import local_maxima
from tissue import *
import tissue as T
import fullfield_tissue as F
DEST=F.DEST
CACHE=T.CACHE/'fullfield-comparison'

def select_plane(images,support):
    assert images.shape==(16,30,128,128) and support.shape==images.shape
    evidence=np.zeros(30,int)
    for c in range(4):
        a=images[c]
        peaks=local_maxima(a,connectivity=3,allow_borders=True)&(a>.2)
        # Count connected plateaus, not their voxel area. Lower z wins ties.
        labels,n=label(peaks,structure=np.ones((3,3,3)))
        from scipy.ndimage import center_of_mass
        for center in center_of_mass(peaks,labels,np.arange(1,n+1)):
            z,y,x=np.floor(np.asarray(center)+.5).astype(int)
            if 24<=y<104 and 24<=x<104:evidence[z]+=1
    safe=np.array([z for z in range(4,26) if np.all(support[:,z-1:z+2,24:104,24:104])])
    assert len(safe),'No supported candidate plane'
    selected=int(min(safe,key=lambda z:(-int(evidence[z]),int(z))))
    return selected,dict(selected_z=selected,reference_spot_counts_by_z=evidence.tolist(),safe_planes=safe.tolist(),selected_reference_spots=int(evidence[selected]),selection_uses_round_one_only=True,selection_uses_gene_identity_or_method_outcomes=False)


def guard():
    T.guard();spec=json.loads((DEST/'comparison-freeze.json').read_text())
    for rel,digest in spec['files'].items():
        assert A.sha(A.ROOT/rel)==digest,rel
        assert subprocess.check_output(['git','show','HEAD:'+rel],cwd=A.ROOT)==(A.ROOT/rel).read_bytes(),rel
    assert json.loads((DEST/'state.json').read_text())['stage']=='eligible'

def execute():
    guard();state=json.loads((DEST/'native-state.json').read_text());assert state['stage']=='eligible'
    codes,_=book();data=[]
    for k in [0,4]:
        with np.load(CACHE/f'conditioned-crop-{k}.npz') as d:data.append((d['images'],d['support']))
    references={dec:json.loads((DEST/'native'/f'4__{dec}.json').read_text())['catalog'] for dec in A.DECODERS};outcomes=[]
    for kind in ['zero','common','color','outclass']:
        imposed=offsets(kind);cal=P.warp(data[0][0],imposed);ca=P.warp(data[0][1],imposed);obs=P.warp(data[1][0],imposed);oa=P.warp(data[1][1],imposed)
        z,plane=select_plane(cal,ca);start=time.monotonic();maps={'none':np.zeros((16,3),int)};diagnostics={};failures={};timing={}
        tick=time.monotonic();common=G.common_fit(cal[:,z],4);timing['common_iterative_s']=time.monotonic()-tick;tick=time.monotonic();template,diag=I.fit(cal[:,z]);timing['template_2d_multistart_s']=time.monotonic()-tick;maps['common_iterative']=np.c_[np.zeros(16,int),common];maps['template_2d_multistart']=np.c_[np.zeros(16,int),template];diagnostics['template_2d_multistart']=diag
        tick=time.monotonic();models=[]
        for name,initial in [('zero',np.zeros((16,2),int)),('common',common),('template',template)]:
            try:
                model=A.R.qualify.bardensr.registration.find_translations_using_model(cal[:,z,None].astype(float),codes,maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(16),-initial].astype(float));m=-np.asarray(model.corrections);losses=np.asarray(model.losses);assert np.isfinite(m).all() and np.isfinite(losses).all();models.append((float(losses[-1]),np.rint(m).astype(int),name))
            except Exception:failures['author_'+name]=traceback.format_exc()
        timing['author_multistart_s']=time.monotonic()-tick
        if models:loss,m,name=min(models,key=lambda t:t[0]);maps['author_multistart']=m;diagnostics['author_multistart']=dict(start=name,loss=loss)
        tick=time.monotonic()
        try:maps['template_3d_multistart']=maps3d(cal)
        except Exception:failures['template_3d_multistart']=traceback.format_exc()
        timing['template_3d_multistart_s']=time.monotonic()-tick;tick=time.monotonic()
        try:m,diag=G.fit(cal[:,z],ca[:,z],codes);maps['global_group']=np.c_[np.zeros(16,int),m];diagnostics['global_group']=diag
        except Exception:failures['global_group']=traceback.format_exc()
        timing['global_group_s']=time.monotonic()-tick
        maps['oracle']=np.c_[np.zeros(16,int),-imposed]
        aligned={p:(move3d(obs,m),move3d(oa,m)) for p,m in maps.items()}
        tick=time.monotonic()
        try:a,b,diag=source_prepare(obs,oa);aligned['source_round_3d']=(a,b);diagnostics['source_round_3d']=diag
        except Exception:failures['source_round_3d']=traceback.format_exc()
        timing['source_round_3d_s']=time.monotonic()-tick
        vz,vplane=select_plane(obs,oa);losses={p:float(LR.calc_loss(a[:,vz].astype(float),codes,np.full((16,2),-20),zero_padding=20)) for p,(a,b) in aligned.items() if p not in ['none','oracle']}
        conventional=[p for p in losses if p!='global_group'];best=min(conventional,key=lambda p:losses[p]);selected='global_group' if 'global_group' in losses and losses[best]-losses['global_group']>=.005 else best
        geometry=dict(kind=kind,calibration_plane=plane,validation_plane=vplane,maps_zyx={p:m.tolist() for p,m in maps.items()},coding_losses=losses,policy_geometry_runtime=timing,selected=selected,failures=failures,diagnostics=diagnostics,geometry_s=time.monotonic()-start,imposed_for_evaluation_only=imposed.tolist());A.write(DEST/f'{kind}-geometry.json',geometry);print('TISSUE FIT',kind,selected,losses,flush=True)
        for policy in POLICIES:
            for decoder in A.DECODERS:
                try:
                    a,b=aligned[policy];tick=time.monotonic()
                    if policy=='oracle':
                        points,_=V.candidates(a,b);original,_=V.candidates(*data[1]);assert np.array_equal(points,original) and np.array_equal(traces(a,points),traces(data[1][0],original));pred=references[decoder]|dict(exact_reference_kernel_reused=True)
                    else:pred=catalog(a,b,codes,decoder)
                    metrics={str(r):P.metrics(references[decoder],pred,r) for r in [1.5,3,5]};row=dict(kind=kind,policy=policy,decoder=decoder,status='completed',prediction=pred,metrics=metrics,decode_s=time.monotonic()-tick,policy_geometry_s=timing.get(policy+'_s',0.))
                except Exception:row=dict(kind=kind,policy=policy,decoder=decoder,status='failed',traceback=traceback.format_exc())
                outcomes.append(row);A.write(DEST/'outcomes'/f'{kind}__{policy}__{decoder}.json',row);print('TISSUE SCORE',kind,policy,decoder,row['status'],row.get('metrics',{}).get('3',{}).get('identity_f1'),flush=True);gc.collect()
        A.write(DEST/'comparison-state.json',dict(stage='complete' if len(outcomes)==96 else 'partial',planned=96,completed=sum(r['status']=='completed' for r in outcomes),failed=sum(r['status']=='failed' for r in outcomes),source='exposed tissue source-geometry repair',workflow_frozen_before_pixels=False,development=True,independent_confirmation=False,negative_proxy_unavailable=True,native_RNA_truth_unavailable=True))

if __name__=='__main__':
    CACHE.mkdir(exist_ok=True)
    for k in range(5):
        source=T.CACHE/f'fullfield-conditioned-crop-{k}.npz';target=CACHE/f'conditioned-crop-{k}.npz'
        assert A.sha(source)==json.loads((DEST/f'crop-{k}.json').read_text())['array_sha256']
        if not target.exists():os.link(source,target)
        assert A.sha(target)==A.sha(source)
    A.write(DEST/'native-state.json',json.loads((DEST/'state.json').read_text()))
    execute()
