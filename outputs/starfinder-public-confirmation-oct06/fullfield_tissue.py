"""Source full-field spatial information repair; original failure remains intact."""
import gc,hashlib,io,json,time,zlib
from pathlib import Path
import numpy as np
import tifffile
import adapter as A
import local_trace as L
import tissue as T
DEST=A.HERE/'fullfield-tissue';CACHE=T.CACHE

def source_map(fixed,moving):
    coarse=A.register(A.spectrum(fixed[:,::2,::2]),moving[:,::2,::2]);coarse=coarse*np.array([1,2,2]);shape=np.array(fixed.shape);choices=[]
    for dz in [-1,0,1]:
        for dy in [-1,0,1]:
            for dx in [-1,0,1]:
                m=coarse+np.array([dz,dy,dx]);lo=np.maximum(m,0);hi=np.minimum(shape,shape+m)
                if np.any(hi<=lo):continue
                f=tuple(slice(int(a),int(b),1 if k==0 else 2) for k,(a,b) in enumerate(zip(lo,hi)));q=tuple(slice(int(a-c),int(b-c),1 if k==0 else 2) for k,(a,b,c) in enumerate(zip(lo,hi,m)))
                score=float(np.multiply(fixed[f],moving[q],dtype=np.float32).sum(dtype=np.float64));choices.append((-score,int(np.square(m).sum()),tuple(map(int,m))))
    result=np.array(min(choices) [2],int)
    return result,dict(coarse_zyx=coarse.tolist(),selected_zyx=result.tolist(),candidate_scores=[dict(score=-score,shift_zyx=list(m)) for score,_,m in choices],XY_subsample=2,Z_subsample=1)

def load(record):
    path=Path(record['local_compressed_path']);assert A.sha(path)==record['compressed_sha256'];raw=zlib.decompress(path.read_bytes(),-15);assert hashlib.sha256(raw).hexdigest()==record['sha256'] and zlib.crc32(raw)==int(record['crc32'],16)
    image=tifffile.imread(io.BytesIO(raw));assert image.shape==(30,3072,3072) and image.dtype==np.uint8
    return np.rot90(image,k=-1,axes=(1,2))

def main():
    T.guard();controls=json.loads((DEST/'controls.json').read_text());assert controls['passed'] and controls['script_sha256']==A.sha(__file__)
    records=json.loads((A.HERE/'tissue-retrieval.json').read_text())['records'];profiles=json.loads((T.DEST/'profiles.json').read_text())['records'];outputs=np.zeros((5,16,30,128,128),np.float32);alpha=np.zeros_like(outputs,np.uint8);maps=[];diagnostics=[];reference=None
    def source(r,c):return next(x for x in records if f'/round{r+1}/' in x['member'] and x['member'].endswith('_'+A.CHANNELS[c]+'.tif'))
    def lut(r,c):
        profile=next(x for x in profiles if x['round']==r and x['color']==c);return np.asarray(profile['histogram_lut'],np.uint8)[np.asarray(profile['minmax_lut'],np.uint8)]
    for r in range(4):
        start=time.monotonic();merged=np.zeros((30,3072,3072),np.uint8)
        for c in range(4):
            image=load(source(r,c));np.maximum(merged,lut(r,c)[image],out=merged);del image;gc.collect()
        if reference is None:reference=merged.copy();shift=np.zeros(3,int);diag=dict(identity=True)
        else:shift,diag=source_map(reference,merged)
        del merged;gc.collect();maps.append(shift.tolist());diagnostics.append(diag)
        assert abs(shift[0])<26 and max(abs(shift[1:]))<384,'Full-field geometry unsupported by fixed scoring windows'
        for c in range(4):
            image=load(source(r,c));table=lut(r,c)
            for k,(y,x) in enumerate(T.CENTERS):
                y0=y-64-shift[1];x0=x-64-shift[2];patch=table[image[:,y0:y0+128,x0:x0+128]];assert patch.shape==(30,128,128)
                outputs[k,r*4+c]=A.integer_warp(patch,(shift[0],0,0)).astype(np.float32)/255
                z=np.arange(30)-shift[0];alpha[k,r*4+c]=((z>=0)&(z<30))[:,None,None]
            del image;gc.collect()
        print('FULLFIELD ROUND',r,shift.tolist(),time.monotonic()-start,flush=True);A.write(DEST/'geometry.json',dict(maps_zyx=maps,diagnostics=diagnostics,source_profiles_sha256=A.sha(T.DEST/'profiles.json'),repair_development=True))
    del reference;gc.collect()
    for k in range(5):
        images=outputs[k].copy();support=alpha[k].copy();fixed=images[:4].max(axis=0);local=[]
        for r in range(4):
            flow,diag=L.fit(fixed,images[r*4:r*4+4].max(axis=0));local.append(diag)
            for c in range(4):f=r*4+c;images[f]=L.warp(images[f],flow);support[f]=L.warp(support[f],flow)>=1-1e-6
        target=CACHE/f'fullfield-conditioned-crop-{k}.npz';np.savez_compressed(target,images=images,support=support)
        A.write(DEST/f'crop-{k}.json',dict(crop=k,array_sha256=A.sha(target),local=local,source_geometry_sha256=A.sha(DEST/'geometry.json')));print('FULLFIELD CROP',k,flush=True)
    codes,_=T.book();rows=[]
    for k in range(5):
        with np.load(CACHE/f'fullfield-conditioned-crop-{k}.npz') as d:images=d['images'];support=d['support']
        for dec in A.DECODERS:
            row=dict(crop=k,decoder=dec,status='completed',catalog=T.catalog(images,support,codes,dec));rows.append(row);A.write(DEST/'native'/f'{k}__{dec}.json',row);print('FULLFIELD NATIVE',k,dec,row['catalog']['coding'],row['catalog']['coverage'],flush=True)
    groups={}
    for dec in A.DECODERS:
        part=[r for r in rows if r['decoder']==dec];g=dict(coding=sum(r['catalog']['coding'] for r in part),minimum_coverage=min(r['catalog']['coverage'] for r in part),designated_counts=[next(r['catalog']['coding'] for r in part if r['crop']==k) for k in [0,4]]);g['eligible']=g['coding']>=50 and g['minimum_coverage']>=.90 and min(g['designated_counts'])>=10;groups[dec]=g
    A.write(DEST/'state.json',dict(stage='eligible' if all(g['eligible'] for g in groups.values()) else 'failed_eligibility',groups=groups,completed=15,development=True,original_frozen_tissue_test_failed=True,independent_confirmation=False,negative_proxy_unavailable=True))
if __name__=='__main__':main()
