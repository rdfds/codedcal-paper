"""Frozen four-round public tissue adapter/comparison, conditional catalog target."""
import argparse,gc,hashlib,io,json,time,traceback,zlib
from pathlib import Path
import numpy as np
import tifffile
import adapter as A
import active_plane as AP
import local_trace as L
import volume_native as V
import initialization as I
import probe as P
import global_repair as G
from bardensr.registration import lowrankregistration as LR

DEST=A.HERE/'tissue';CACHE=A.ROOT.parent/'starfinder-tissue-inputs-oct06.nosync'
CENTERS=[(768,768),(768,2304),(1536,1536),(2304,768),(2304,2304)]
CUTS={'source_spot':.019306540489196777,'correlation':0.,'regression':0.}
POLICIES=['none','common_iterative','template_2d_multistart','author_multistart','template_3d_multistart','source_round_3d','global_group','oracle']

def book():
    source=json.loads((A.AUDIT/'starfinder-color-codebooks.json').read_text())['source_books'][0]
    words=np.array([[int(c)-1 for c in r['color_word']] for r in source['items']],int)
    assert words.shape==(64,4) and len(set(map(tuple,words)))==64
    codes=np.zeros((4,4,64),np.float32)
    for j,w in enumerate(words):codes[np.arange(4),w,j]=1
    assert np.all(np.bitwise_xor.reduce(words,axis=1)==0)
    return codes.reshape(16,64),words

def traces(images,points):
    raw=np.stack([images[:,z-1:z+2,y-2:y+3,x-2:x+3].sum(axis=(1,2,3),dtype=np.float64) for x,y,z in points],axis=1) if len(points) else np.empty((16,0))
    cube=raw.T.reshape(-1,4,4);cube/=np.sqrt((cube*cube).sum(axis=2,keepdims=True))+1e-6
    return cube.reshape(-1,16).T.astype(np.float32)

def catalog(images,support,codes,decoder):
    start=time.monotonic();points,rejected=V.candidates(images,support);coverage=len(points)/(len(points)+rejected) if len(points)+rejected else 0.
    values=traces(images,points)
    if not len(points):return dict(xyz=[],words=[],coding=0,negative=None,candidate_objects=0,unsupported_candidates=rejected,abstained=0,coverage=coverage,runtime=dict(elapsed_s=time.monotonic()-start))
    if decoder=='source_spot':
        cube=values.T.reshape(-1,4,4);ordered=np.sort(cube,axis=2);score=(ordered[:,:,-1]-ordered[:,:,-2]).min(axis=1)
        lookup={tuple(w):j for j,w in enumerate(codes.reshape(4,4,-1).argmax(axis=1).T)}
        calls=np.array([lookup.get(tuple(w),-1) for w in cube.argmax(axis=2)],int);eligible=(calls>=0)&(score>0)
    else:
        data=values[:,None,None,:]
        if decoder=='correlation':density=A.R.qualify.bardensr.spot_calling.estimate_density_singleshot(data,codes,.05)
        else:density,_=A.R.qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),codes,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
        density=np.asarray(density,np.float32).reshape(len(points),64);assert np.isfinite(density).all()
        ordered=np.sort(density,axis=1);first=ordered[:,-1];second=ordered[:,-2];calls=density.argmax(axis=1)
        score=first-second if decoder=='correlation' else (first-second)/np.maximum(density.sum(axis=1),1e-20)
        eligible=(score>0)&(first>(.775 if decoder=='correlation' else .2*density.max()))
    keep=eligible&(score>=CUTS[decoder])
    return dict(xyz=points[keep].tolist(),words=calls[keep].tolist(),coding=int(keep.sum()),negative=None,candidate_objects=len(points),unsupported_candidates=rejected,abstained=int((~keep).sum()),coverage=coverage,runtime=dict(elapsed_s=time.monotonic()-start),negative_proxy_unavailable=True)

def maps3d(images,roundwise=False):
    rounds=len(images)//4
    merged=images.reshape(rounds,4,*images.shape[1:]).max(axis=1);ref=A.spectrum(merged[0]);common=np.stack([A.register(ref,a) for a in merged]);common=np.repeat(np.rint(common-common.mean(axis=0)).astype(int),4,axis=0)
    if roundwise:return common
    ref=A.spectrum(images[0]);first=np.stack([A.register(ref,a) for a in images]);first=np.rint(first-first.mean(axis=0)).astype(int)
    choices=[]
    for initial in [np.zeros((len(images),3),int),common,first]:
        def bounded(m):return np.max(np.abs(m[:,0]))<=3 and np.max(np.abs(m[:,1:]))<=10
        if not bounded(initial):continue
        maps=initial.copy()
        for step in range(10):
            aligned=move3d(images,maps);ref=A.spectrum(aligned.mean(axis=0));new=maps+np.stack([A.register(ref,a) for a in aligned]);new=np.rint(new-new.mean(axis=0)).astype(int)
            if not bounded(new) or np.array_equal(new,maps):break
            maps=new
        aligned=move3d(images,maps)[:,4:26,16:112,16:112].astype(float);loss=float(np.square(aligned-aligned.mean(axis=0)).sum()/max(np.square(aligned).sum(),1e-20));choices.append((loss,float(np.square(maps).sum()),maps))
    assert choices;return min(choices,key=lambda t:t[:2])[2]

def move3d(images,maps):return np.stack([A.integer_warp(frame,shift) for frame,shift in zip(images,maps)])

def source_prepare(images,support):
    start=time.monotonic();maps=maps3d(images,True);out=move3d(images,maps);alpha=move3d(support,maps);fixed=out[:4].max(axis=0);diag=[]
    for r in range(len(images)//4):
        flow,record=L.fit(fixed,out[r*4:r*4+4].max(axis=0))
        for c in range(4):f=r*4+c;out[f]=L.warp(out[f],flow);alpha[f]=L.warp(alpha[f],flow)>=1-1e-6
        diag.append(record)
    return out,alpha,dict(maps_zyx=maps.tolist(),local=diag,geometry_s=time.monotonic()-start)

def offsets(kind):
    v=np.array([4,2])
    if kind=='zero':return np.zeros((16,2),int)
    if kind=='common':return np.repeat(np.r_[np.tile(v,(2,1)),np.tile(-v,(2,1))],4,axis=0)
    if kind=='color':return np.tile(np.r_[np.tile(v,(2,1)),np.tile(-v,(2,1))],(4,1))
    return np.tile(np.array([[3,1],[-2,4],[1,-3],[-2,-2]]),(4,1))

def guard():
    L.guard();spec=json.loads((DEST/'freeze.json').read_text())
    for rel,digest in spec['files'].items():assert A.sha(A.ROOT/rel)==digest,rel
    assert json.loads((DEST/'controls.json').read_text())['passed']
    return spec

def prepare():
    guard();meta=json.loads((A.HERE/'tissue-retrieval.json').read_text());assert meta['stage']=='complete' and meta['completed']==17
    state=dict(stage='preparing',source='tissue-2D',field='tile_1',raw_intensity_pixels_examined=True,independent_comparative_freeze_sha256=A.sha(DEST/'freeze.json'),profiles_completed=0);A.write(DEST/'exposure.json',state)
    crops=np.zeros((5,16,30,128,128),np.float32);records=[];reference=None
    for r in range(4):
        for c,channel in enumerate(A.CHANNELS):
            source=next(t for t in meta['records'] if f'/round{r+1}/' in t['member'] and t['member'].endswith('_'+channel+'.tif'))
            path=Path(source['local_compressed_path']);assert A.sha(path)==source['compressed_sha256'];raw=zlib.decompress(path.read_bytes(),-15);assert hashlib.sha256(raw).hexdigest()==source['sha256'] and zlib.crc32(raw)==int(source['crc32'],16)
            image=tifffile.imread(io.BytesIO(raw));assert image.shape==(30,3072,3072) and image.dtype==np.uint8
            image=np.rot90(image,k=-1,axes=(1,2));counts=A.histogram(image);minmax=A.minmax_lut(counts)
            if reference is None:reference=A.mapped_hist(counts,minmax)
            match=A.matching_lut(A.mapped_hist(counts,minmax),reference);lut=match[minmax]
            for k,(y,x) in enumerate(CENTERS):crops[k,r*4+c]=lut[image[:,y-64:y+64,x-64:x+64]].astype(np.float32)/255
            records.append(dict(round=r,color=c,source_sha256=source['sha256'],source_compressed_sha256=source['compressed_sha256'],histogram=counts.tolist(),minmax_lut=minmax.tolist(),histogram_lut=match.tolist()))
            state['profiles_completed']=len(records);A.write(DEST/'exposure.json',state);print('TISSUE PROFILE',r,c,flush=True);del raw,image;gc.collect()
    A.write(DEST/'profiles.json',dict(records=records,source_retrieval_sha256=A.sha(A.HERE/'tissue-retrieval.json')))
    for k in range(5):
        images,support,diagnostics=source_prepare(crops[k],np.ones_like(crops[k],np.uint8));target=CACHE/f'conditioned-crop-{k}.npz';np.savez_compressed(target,images=images,support=support)
        A.write(DEST/f'crop-{k}-preparation.json',dict(crop=k,center_yx=CENTERS[k],diagnostics=diagnostics,array_sha256=A.sha(target),freeze_sha256=A.sha(DEST/'freeze.json')));print('TISSUE PREPARED',k,flush=True)
    state['stage']='prepared';A.write(DEST/'exposure.json',state)

def native():
    guard();codes,_=book();rows=[]
    for k in range(5):
        path=CACHE/f'conditioned-crop-{k}.npz';assert A.sha(path)==json.loads((DEST/f'crop-{k}-preparation.json').read_text())['array_sha256']
        with np.load(path) as d:images=d['images'];support=d['support']
        for decoder in A.DECODERS:
            r=dict(crop=k,decoder=decoder,status='completed',catalog=catalog(images,support,codes,decoder),input_sha256=A.sha(path));rows.append(r);A.write(DEST/'native'/f'{k}__{decoder}.json',r);print('TISSUE NATIVE',k,decoder,r['catalog']['coding'],r['catalog']['coverage'],flush=True)
    groups={}
    for dec in A.DECODERS:
        part=[r for r in rows if r['decoder']==dec];groups[dec]=dict(coding=sum(r['catalog']['coding'] for r in part),minimum_coverage=min(r['catalog']['coverage'] for r in part),designated_counts=[next(r['catalog']['coding'] for r in part if r['crop']==k) for k in [0,4]])
        g=groups[dec];g['eligible']=g['coding']>=50 and g['minimum_coverage']>=.90 and min(g['designated_counts'])>=10
    A.write(DEST/'native-state.json',dict(stage='eligible' if all(g['eligible'] for g in groups.values()) else 'failed_eligibility',groups=groups,completed=15,planned=15,negative_proxy_unavailable=True,conditional_target=True,original_whole_interface_gates_not_passed=True))

def execute():
    guard();state=json.loads((DEST/'native-state.json').read_text());assert state['stage']=='eligible'
    codes,_=book();data=[]
    for k in [0,4]:
        with np.load(CACHE/f'conditioned-crop-{k}.npz') as d:data.append((d['images'],d['support']))
    references={dec:json.loads((DEST/'native'/f'4__{dec}.json').read_text())['catalog'] for dec in A.DECODERS};outcomes=[]
    for kind in ['zero','common','color','outclass']:
        imposed=offsets(kind);cal=P.warp(data[0][0],imposed);ca=P.warp(data[0][1],imposed);obs=P.warp(data[1][0],imposed);oa=P.warp(data[1][1],imposed)
        z,plane=AP.select(cal,ca);start=time.monotonic();maps={'none':np.zeros((16,3),int)};diagnostics={};failures={};timing={}
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
        vz,vplane=AP.select(obs,oa);losses={p:float(LR.calc_loss(a[:,vz].astype(float),codes,np.full((16,2),-20),zero_padding=20)) for p,(a,b) in aligned.items() if p not in ['none','oracle']}
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
        A.write(DEST/'comparison-state.json',dict(stage='complete' if len(outcomes)==96 else 'partial',planned=96,completed=sum(r['status']=='completed' for r in outcomes),failed=sum(r['status']=='failed' for r in outcomes),source='distinct tissue experiment',workflow_frozen_before_pixels=True,negative_proxy_unavailable=True,native_RNA_truth_unavailable=True))

if __name__=='__main__':
    cli=argparse.ArgumentParser();cli.add_argument('mode',choices=['prepare','native','execute']);args=cli.parse_args();globals()[args.mode]()
