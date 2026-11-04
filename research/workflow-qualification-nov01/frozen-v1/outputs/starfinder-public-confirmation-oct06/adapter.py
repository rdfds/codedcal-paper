"""3D source-family qualification; central-plane endpoint, no outcome selection."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '2')
os.environ.setdefault('TF_NUM_INTRAOP_THREADS', '2')
os.environ.setdefault('TF_NUM_INTEROP_THREADS', '1')
import argparse
import gc
import hashlib
import importlib.metadata as metadata
import itertools
import json
from pathlib import Path
import subprocess
import sys
import time
import traceback
import numpy as np
from scipy import fft
from scipy.ndimage import center_of_mass, label, maximum_filter
from skimage.morphology import local_maxima
import tifffile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
AUDIT = HERE.parent / 'barseq-public-substitute-oct06'
CACHE = ROOT.parent / 'starfinder-inputs-oct06.nosync'
CENTERS = [(374,374),(374,1122),(748,748),(1122,374),(1122,1122)]
CHANNELS = ['ch00','ch02','ch01','ch03']
DECODERS = ['source_spot','correlation','regression']
sys.path.insert(0, str(HERE.parent / 'barseq3-native-repair-oct05'))
import recovery as R

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block := stream.read(1024**2): digest.update(block)
    return digest.hexdigest()

def write(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n'); temp.replace(path)

def book():
    source = json.loads((AUDIT / 'starfinder-color-codebooks.json').read_text())['source_books'][1]
    words = [tuple(int(c)-1 for c in row['color_word']) for row in source['items']]
    absent = [w for w in itertools.product(range(4), repeat=6)
              if np.bitwise_xor.reduce(w)==0 and w not in words]
    assert len(words)==998 and len(set(words))==998 and len(absent)==26
    codes = np.zeros((6,4,1024), np.float32)
    for j, word in enumerate(words + absent): codes[np.arange(6),word,j] = 1
    names = [row['name'] for row in source['items']] + ['no_gene_'+''.join(str(c+1) for c in w) for w in absent]
    return codes.reshape(24,1024), np.r_[np.zeros(998,bool),np.ones(26,bool)], names, words+absent

def histogram(a):
    counts = np.zeros(256, np.int64)
    for plane in a: counts += np.bincount(plane.ravel(), minlength=256)
    assert counts.sum()==a.size
    return counts

def minmax_lut(counts):
    occupied = np.flatnonzero(counts); assert len(occupied)>1
    lo, hi = int(occupied[0]), int(occupied[-1])
    return np.floor(np.clip((np.arange(256)-lo)/(hi-lo),0,1)*255+.5).astype(np.uint8)

def mapped_hist(counts, lut):
    out = np.zeros(256,np.int64); np.add.at(out,lut,counts); return out

def matching_lut(counts, reference):
    bins = np.floor(np.arange(256)*63/255+.5).astype(int)
    ref = np.zeros(64); np.add.at(ref,bins,reference)
    a = np.cumsum(counts)/counts.sum(); b = np.cumsum(ref)/ref.sum()
    indices = np.argmin(np.abs(a[:,None]-b[None,:]),axis=1)
    return np.floor(indices*255/63+.5).astype(np.uint8)

def integer_warp(a, shift):
    shift = np.asarray(shift,int); assert a.ndim==len(shift)
    out = np.roll(a, tuple(shift), axis=tuple(range(a.ndim)))
    for axis, amount in enumerate(shift):
        slices = [slice(None)]*a.ndim
        if amount>0: slices[axis]=slice(0,amount)
        elif amount<0: slices[axis]=slice(amount,None)
        else: continue
        out[tuple(slices)]=0
    return out

def spectrum(a): return fft.rfftn(a.astype(np.float32), workers=2)

def register(reference_fft, moving):
    moving_fft = spectrum(moving)
    np.conjugate(moving_fft,out=moving_fft); moving_fft *= reference_fft
    correlation = fft.irfftn(moving_fft, s=moving.shape, workers=2)
    del moving_fft
    np.abs(correlation,out=correlation)
    peak = np.array(np.unravel_index(np.argmax(correlation),moving.shape))
    shape = np.array(moving.shape); peak = np.where(peak>=(shape+1)//2,peak-shape,peak)
    return peak.astype(int)

def source_decode(images, codes):
    # images are frame,z,y,x. Each channel's 26-connected plateaus are detected
    # independently, as source SpotFindingMax3D does. Positive plateau centroids
    # round half up, MATLAB's int16 conversion convention.
    t0=time.monotonic(); rounds=images.shape[0]//4; assert rounds==6
    words = np.argmax(codes.reshape(6,4,-1),axis=1).T
    lookup = {tuple(w):j for j,w in enumerate(words)}
    candidates=[]
    for c in range(4):
        a=images[c]
        # Reference-round full-volume min/max + LUT has maximum 1. A cropped
        # maximum would lower the author threshold in dim ROIs.
        region=local_maxima(a,connectivity=3,allow_borders=True)&(a>.2)
        plateaus,n=label(region,structure=np.ones((3,3,3)))
        for centroid in center_of_mass(region,plateaus,np.arange(1,n+1)):
            center=np.floor(np.asarray(centroid)+.5).astype(int)
            z,y,x=center
            if z==14 and 24<=y<104 and 24<=x<104: candidates.append(center)
    candidates=np.asarray(candidates,int).reshape(-1,3)
    if not len(candidates): return np.empty((0,2)),np.empty(0,int),dict(kernel_s=time.monotonic()-t0,candidates=0)
    values=np.empty((len(candidates),rounds,4),np.float32)
    for r in range(rounds):
        for c in range(4):
            # Exact clipped 5x5x3 sums at supported interior coordinates.
            a=images[r*4+c]
            values[:,r,c]=[a[z-1:z+2,y-2:y+3,x-2:x+3].sum(dtype=np.float64) for z,y,x in candidates]
    maxima=values.max(axis=2); unique=(values==maxima[:,:,None]).sum(axis=2)==1
    readwords=values.argmax(axis=2); xy=[]; calls=[]
    for point,word,good in zip(candidates,readwords,unique.all(axis=1)):
        if good and tuple(word) in lookup: xy.append([int(point[2]),int(point[1])]); calls.append(lookup[tuple(word)])
    order=np.lexsort((np.asarray(xy).reshape(-1,2)[:,0],np.asarray(xy).reshape(-1,2)[:,1],calls)) if calls else []
    return np.asarray(xy,float).reshape(-1,2)[order],np.asarray(calls,int)[order],dict(kernel_s=time.monotonic()-t0,candidates=len(candidates),declared_qscore_not_applied=True)

def decode(images, codes, decoder):
    if decoder=='source_spot': return source_decode(images,codes)
    start=time.monotonic(); data=images[:,13:16,16:112,16:112]
    if decoder=='correlation': density=R.qualify.bardensr.spot_calling.estimate_density_singleshot(data,codes,.05)
    else: density,_=R.qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),codes,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
    density=np.asarray(density,np.float32); assert np.isfinite(density).all()
    kernel_s=time.monotonic()-start; cut=float(.775 if decoder=='correlation' else .2*density.max())
    maxima=maximum_filter(density,size=(3,3,3,1),mode='constant',cval=np.finfo(np.float32).min)
    z,y,x,j=np.where((density==maxima)&(density>cut)); x+=16; y+=16
    keep=(z==1)&(x>=24)&(x<104)&(y>=24)&(y<104); x,y,j=x[keep],y[keep],j[keep]
    order=np.lexsort((x,y,j))
    return np.c_[x[order],y[order]].astype(float),j[order],dict(kernel_s=kernel_s,elapsed_s=time.monotonic()-start,absolute_cut=cut,density_max=float(density.max()),density_sha256=hashlib.sha256(density.tobytes()).hexdigest())

def guard():
    frozen=json.loads((HERE/'freeze.json').read_text())
    for rel,digest in frozen['files'].items():
        assert sha(ROOT/rel)==digest,rel
        assert subprocess.check_output(['git','show','HEAD:'+rel],cwd=ROOT)==(ROOT/rel).read_bytes(),rel
    for name,version in frozen['versions'].items(): assert metadata.version(name)==version,name
    for rel,digest in frozen['upstream_files'].items():assert sha(ROOT.parent/rel)==digest,rel
    assert json.loads((HERE/'controls.json').read_text())['passed']
    return frozen

def prepare():
    total_start=time.monotonic()
    guard(); dest=CACHE/'conditioned-Position351'; dest.mkdir(exist_ok=True)
    marker=dest/'preparation.json'
    if marker.exists():
        out=json.loads(marker.read_text()); assert out['freeze_sha256']==sha(HERE/'freeze.json')
        assert all(sha(dest/r['file'])==r['sha256'] for r in out['outputs']); return out
    records=json.loads((AUDIT/'starfinder-culture-retrieval.json').read_text())['records']
    write(HERE/'exposure-state.json',dict(field='Position351',stage='development_full_volume_profiles',intensity_pixels_examined=True,profiles_completed=0,independent_study_frozen=False,independent_confirmation=False))
    inputs=[]; profiles={}; luts={}
    def find(r,c): return next(t for t in records if '/round'+str(r+1)+'/' in t['member'] and t['member'].endswith('_'+CHANNELS[c]+'.tif'))
    def load(record):
        path=Path(record['local_path']); assert sha(path)==record['sha256']
        a=tifffile.imread(path); assert a.shape==(30,1496,1496) and a.dtype==np.uint8
        return np.rot90(a,k=-1,axes=(1,2))
    for r in range(6):
        for c in range(4):
            record=find(r,c); a=load(record); counts=histogram(a); lut=minmax_lut(counts)
            profiles[r,c]=dict(source_sha256=record['sha256'],counts=counts.tolist(),minmax_lut=lut.tolist())
            luts[r,c]=lut; inputs.append(dict(member=record['member'],sha256=record['sha256']))
            write(HERE/'exposure-state.json',dict(field='Position351',stage='development_full_volume_profiles',intensity_pixels_examined=True,profiles_completed=len(inputs),independent_study_frozen=False,independent_confirmation=False))
            del a; gc.collect()
    reference_hist=mapped_hist(np.array(profiles[0,0]['counts']),luts[0,0])
    for key,lut in list(luts.items()):
        match=matching_lut(mapped_hist(np.array(profiles[key]['counts']),lut),reference_hist)
        profiles[key]['histogram_lut']=match.tolist(); luts[key]=match[lut]
    def merged(r):
        out=np.zeros((30,1496,1496),np.uint8)
        for c in range(4):
            a=load(find(r,c)); np.maximum(out,luts[r,c][a],out=out); del a
        return out
    t0=time.monotonic(); ref=merged(0); ref_fft=spectrum(ref); del ref; gc.collect()
    shifts=[]; outputs=[[] for _ in CENTERS]; supports=[[] for _ in CENTERS]; timings=[]
    for r in range(6):
        start=time.monotonic()
        if r==0: maps=np.zeros(3,int)
        else:
            moving=merged(r); maps=register(ref_fft,moving); del moving; gc.collect()
        assert np.max(np.abs(maps[1:]))<=100 and abs(maps[0])<=10,'Global support shift too large'
        shifts.append(maps.tolist())
        for c in range(4):
            a=load(find(r,c)); a=integer_warp(luts[r,c][a],maps)
            for k,(y,x) in enumerate(CENTERS):
                outputs[k].append(a[:,y-64:y+64,x-64:x+64].astype(np.float32)/255)
                zgrid,ygrid,xgrid=np.ogrid[:30,y-64:y+64,x-64:x+64]
                support=((zgrid-maps[0]>=0)&(zgrid-maps[0]<30)&(ygrid-maps[1]>=0)&(ygrid-maps[1]<1496)&(xgrid-maps[2]>=0)&(xgrid-maps[2]<1496))
                supports[k].append(support.astype(np.uint8))
            del a;gc.collect()
        timings.append(time.monotonic()-start);print('PREPARED',r,maps.tolist(),timings[-1],flush=True)
    del ref_fft; paths=[]
    for k in range(5):
        path=dest/f'crop-{k}.npz'; np.savez_compressed(path,images=np.stack(outputs[k]),support=np.stack(supports[k]))
        paths.append(dict(file=path.name,sha256=sha(path),center_yx=list(CENTERS[k])))
    result=dict(field='Position351',intensity_pixels_examined=True,inputs=inputs,profiles={f'{r}-{c}':p for (r,c),p in profiles.items()},shifts_zyx=shifts,outputs=paths,round_processing_s=timings,geometry_and_output_s=time.monotonic()-t0,total_preparation_s=time.monotonic()-total_start,freeze_sha256=sha(HERE/'freeze.json'))
    write(marker,result); write(HERE/'preparation.json',result); return result

def native():
    preparation=prepare(); codes,negative,names,_=book(); results=[]
    for k in range(5):
        path=CACHE/'conditioned-Position351'/f'crop-{k}.npz'
        with np.load(path) as data: images=data['images'];support=data['support']
        fraction=float(np.all(support[:,13:16,24:104,24:104],axis=(0,1)).mean())
        for decoder in DECODERS:
            target=HERE/'native'/f'crop-{k}__{decoder}.json'
            if target.exists():
                result=json.loads(target.read_text());assert result['freeze_sha256']==sha(HERE/'freeze.json')
            else:
                try:
                    xy,calls,runtime=decode(images,codes,decoder)
                    result=dict(status='completed',crop=k,decoder=decoder,xy=xy.tolist(),words=calls.tolist(),coding=int((~negative[calls]).sum()),negative=int(negative[calls].sum()),support=fraction,runtime=runtime)
                except Exception: result=dict(status='failed',crop=k,decoder=decoder,traceback=traceback.format_exc())
                result.update(input_sha256=sha(path),freeze_sha256=sha(HERE/'freeze.json'));write(target,result)
            results.append(result);print('NATIVE',k,decoder,result['status'],result.get('coding'),result.get('negative'),flush=True)
            gc.collect()
    grouped={}
    for decoder in DECODERS:
        rows=[r for r in results if r['decoder']==decoder];n=sum(r.get('coding',0) for r in rows);b=sum(r.get('negative',0) for r in rows)
        complete=all(r['status']=='completed' for r in rows);proxy=b/n*998/26 if n else None;support=min(r.get('support',0) for r in rows)
        grouped[decoder]=dict(completed=complete,coding=n,negative=b,scaled_proxy=proxy,support=support,qualified=complete and n>=50 and proxy<=.1 and support>=.99)
    write(HERE/'state.json',dict(stage='native_qualified' if all(r['qualified'] for r in grouped.values()) else 'completed_native_gate_failure',planned=15,completed=sum(r['status']=='completed' for r in results),failed=sum(r['status']=='failed' for r in results),groups=grouped,adapter_qualified=all(r['qualified'] for r in grouped.values()),independent_comparison_completed=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--native',action='store_true');parser.add_argument('--prepare',action='store_true');args=parser.parse_args()
    try:
        if args.native:native()
        elif args.prepare:prepare()
    except Exception:
        write(HERE/'execution-failure.json',dict(status='failed',traceback=traceback.format_exc(),adapter_sha256=sha(__file__),independent_comparison_completed=False));raise
