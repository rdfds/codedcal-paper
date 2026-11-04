"""Stream the audited published conditioning; no decode or outcome tuning."""
import argparse
import ast
import gc
import hashlib
import json
from pathlib import Path
import time
import traceback
import types
import cv2
import h5py
import numpy as np
import SimpleITK as sitk
from scipy.ndimage import convolve1d
from skimage.feature import blob_log

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CACHE = ROOT.parent / 'barseq3-inputs-oct05.nosync'
SOURCE = CACHE / 'source-audit/barseq3-Snakefile_nopatch.txt'
RAW = CACHE / 'adult-cb-full.hdf5'
MIX = np.array([[1,.15,0,0],[.5,1,0,0],[0,0,1,.5],[0,0,.2,1]]).T
sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(2)
cv2.setNumThreads(2)
cv2.setRNGSeed(503)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    return h.hexdigest()

def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(p)

def read(p):return json.loads(Path(p).read_text())

def bg(a,sigma):
    # Exact last-axis-only behavior of the pinned background_subtraction loop.
    xs=np.arange(-3*sigma,3*sigma+1,dtype=np.float64)
    kernel=np.exp(-.5*xs**2/sigma**2);kernel/=kernel.sum()
    kernel=kernel.astype(a.dtype)
    return np.maximum(a-convolve1d(a,kernel,axis=-1,mode='constant',cval=0),0)

def source_functions():
    # Execute only inspected function ASTs, never the Snakemake file or globals.
    text=SOURCE.read_text();ns=dict(np=np,cv2=cv2,blob_log=blob_log,
        use_cv2_blob_detector=False,path='',
        bardensr=types.SimpleNamespace(preprocessing=types.SimpleNamespace(
            background_subtraction=lambda a,s:bg(a,s[-1]))))
    captured=[]
    ns['save_hdf5']=lambda arrays,p:captured.extend(arrays)
    for name in ['blob_keypoints','channel_alignment_get_transforms_opencv_sm']:
        section=text[text.index('def '+name+'('):]
        section=section.split('\ndef ',1)[0]
        node=ast.parse(section).body[0]
        exec(compile(ast.Module(body=[node],type_ignores=[]),str(SOURCE),'exec'),ns)
    return ns['channel_alignment_get_transforms_opencv_sm'],captured

def features(a):
    b=bg((a/a.max()).astype(float),5).sum(axis=0)
    b=np.clip(b,0,np.percentile(b,99))
    return (b/b.max()*255).astype(np.uint8)

def descriptors(im):
    blobs=blob_log(im,max_sigma=10,threshold=.1)
    pts=[cv2.KeyPoint(x=float(b[1]),y=float(b[0]),size=10) for b in blobs]
    return cv2.SIFT_create().compute(im,pts)

def affine(moving,fixed_desc):
    k,d=descriptors(moving);kf,df=fixed_desc
    if d is None or df is None:raise ValueError('No source round descriptors')
    matches=sorted(cv2.BFMatcher(cv2.NORM_L2,crossCheck=True).match(d,df),key=lambda m:m.distance)
    matches=matches[:int(len(matches)*.8)]
    src=np.float32([k[m.queryIdx].pt for m in matches]).reshape(-1,1,2)
    dst=np.float32([kf[m.trainIdx].pt for m in matches]).reshape(-1,1,2)
    H,inliers=cv2.estimateAffinePartial2D(src,dst,method=cv2.RANSAC,ransacReprojThreshold=1.)
    if H is None or not np.isfinite(H).all():raise ValueError('Source partial-affine fit failed')
    return H.astype(np.float32),dict(retained_matches=len(matches),inliers=int(inliers.sum()))

def camera_warp(a,H):
    a=a.copy()
    for c in [1,3]:
        # Published assignment casts interpolated float32 back to raw uint16.
        a[c]=cv2.warpPerspective(a[c].astype(np.float32),H,(2304,2304))
    return a

def processed(index,r,H,A,support=False):
    if support:a=np.ones((4,2304,2304),np.float32)
    else:
        with h5py.File(RAW) as f:a=f[f'list/{index}'][r]
    a=camera_warp(a,H)
    for c in range(4):a[c]=cv2.warpAffine(a[c],A,(2304,2304))
    return a if support else np.clip(np.einsum('cxy,cd->dxy',a,np.linalg.inv(MIX)),0,None)

def controls():
    rng=np.random.default_rng(503)
    data=rng.uniform(size=(4,1,43,47))
    np.savez_compressed(CACHE/'native-background-control.npz',a64=data,a32=data.astype(np.float32),
        b64s5=bg(data,5),b32s5=bg(data.astype(np.float32),5),
        b64s20=bg(data,20),b32s20=bg(data.astype(np.float32),20))
    a=rng.uniform(size=(45,49));image=sitk.GetImageFromArray(a)
    field=np.zeros((45,49,2));field[...,0]=.35;field[...,1]=-.45
    transform=sitk.DisplacementFieldTransform(sitk.GetImageFromArray(field,isVector=True))
    full=sitk.GetArrayFromImage(sitk.Resample(image,image,transform,sitk.sitkLinear,0.,image.GetPixelID()))
    roi=sitk.Image([17,15],sitk.sitkFloat64);roi.SetOrigin((11.,13.))
    small=sitk.GetArrayFromImage(sitk.Resample(image,roi,transform,sitk.sitkLinear,0.,image.GetPixelID()))
    assert np.array_equal(full[13:28,11:28],small)
    # Known SIFT/affine translation validates map sign and native fit convention.
    im=np.zeros((192,192),np.uint8)
    for x,y,v in rng.integers(20,170,size=(70,3)):cv2.circle(im,(int(x),int(y)),2,int(v),-1)
    shift=np.array([[1,0,3],[0,1,-2]],np.float32)
    moved=cv2.warpAffine(im,shift,(192,192))
    A,diag=affine(moved,descriptors(im))
    assert np.max(np.abs(A-np.array([[1,0,-3],[0,1,2]])))<.05
    write(HERE/'geometry-controls.json',dict(condition_sha256=sha(__file__),
        crop_resample_exact=True,known_partial_affine_max_error=float(np.max(np.abs(A-np.array([[1,0,-3],[0,1,2]])))),
        native_fit=diag,background_control_sha256=sha(CACHE/'native-background-control.npz')))

def geometry(index):
    assert read(HERE/'background-controls.json')['passed']
    assert read(HERE/'full-input-verification.json')['whole_source_verified']
    dest=HERE/f'field-{index}';dest.mkdir(exist_ok=True)
    completed=dest/'geometry.json'
    if completed.exists():return read(completed)
    fn,capture=source_functions();common=np.ones((2304,2304),bool)
    with h5py.File(RAW) as f:
        old=CACHE/f'native-context-source-{index}.npz'
        if old.exists():
            rawold=np.load(old)['images']
            assert np.array_equal(rawold,f[f'list/{index}'][:10,:,960:1344,960:1344].reshape(40,384,384))
        fixed_desc=None
        for r in range(10):
            record=dest/f'round-{r}-geometry.json';raw=f[f'list/{index}'][r]
            if record.exists():
                out=read(record);assert out['condition_sha256']==sha(__file__)
                H=np.array(out['camera']);A=np.array(out['affine'],np.float32)
                camera=camera_warp(raw,H)
                if r==0:fixed_desc=descriptors(features(camera))
            else:
                start=time.monotonic();capture.clear();fn([raw]);H=capture[0]
                if H is None or not np.isfinite(H).all():raise ValueError('Native camera homography failed')
                camera=camera_warp(raw,H);feature=features(camera)
                if r==0:fixed_desc=descriptors(feature)
                A,diag=affine(feature,fixed_desc)
                out=dict(round=r,camera=H.tolist(),affine=A.tolist(),diagnostics=diag,
                    elapsed_s=time.monotonic()-start,raw_sha256=hashlib.sha256(raw.tobytes()).hexdigest(),condition_sha256=sha(__file__))
                write(record,out)
            del camera,raw
            a=processed(index,r,H,A)
            common&=np.all(a!=0,axis=0)
            del a;gc.collect()
            print('geometry',index,r,flush=True)
    yy,xx=np.where(common)
    if not len(yy):raise ValueError('No native common overlap')
    box=[int(yy.min()+20),int(yy.max()-20),int(xx.min()+20),int(xx.max()-20)]
    origin=1088 if index==0 else 1408
    assert box[0]<=origin and box[1]>=origin+128 and box[2]<=origin and box[3]>=origin+128,'Frozen ROI outside native crop'
    out=dict(source_index=index,bbox_y0_y1_x0_x1=box,roi_origin_yx=[origin,origin],
        condition_sha256=sha(__file__),source_sha256=sha(SOURCE),full_input_verification_sha256=sha(HERE/'full-input-verification.json'),
        previous_selective_context_exact=True)
    write(completed,out);return out

def crop(a,box):return a[:,box[0]:box[1],box[2]:box[3]]

def normalized(index,r,geo,profile):
    t=read(HERE/f'field-{index}/round-{r}-geometry.json')
    a=crop(processed(index,r,np.array(t['camera']),np.array(t['affine'],np.float32)),geo['bbox_y0_y1_x0_x1'])
    low=np.array(profile['min1'][r])[:,None,None];high=np.array(profile['max1'][r])[:,None,None]
    assert np.all(high>0)
    a=bg(np.clip((a-low)/high,0,1).astype(np.float32),20)
    if 'min2' in profile:
        low=np.array(profile['min2'][r])[:,None,None];high=np.array(profile['max2'][r])[:,None,None]
        assert np.all(high>0);a=np.clip((a-low)/high,0,1)
    return a

def profile_fit():
    geo=geometry(0);profile=dict(source_index=0,min1=[],max1=[],min2=[],max2=[],
        precision=.005,condition_sha256=sha(__file__),geometry_sha256=sha(HERE/'field-0/geometry.json'),
        deviation='Only source-zero normalization field available; not published five-field median')
    for r in range(10):
        t=read(HERE/f'field-0/round-{r}-geometry.json')
        a=crop(processed(0,r,np.array(t['camera']),np.array(t['affine'],np.float32)),geo['bbox_y0_y1_x0_x1'])
        lo,hi=np.percentile(a,[.005,99.995],axis=(1,2))
        profile['min1'].append(lo.tolist());profile['max1'].append(hi.tolist());del a
        print('normalization stage1',r,flush=True)
    partial={k:v for k,v in profile.items() if k not in ['min2','max2']}
    for r in range(10):
        a=normalized(0,r,geo,partial);lo,hi=np.percentile(a,[.005,99.995],axis=(1,2))
        profile['min2'].append(lo.tolist());profile['max2'].append(hi.tolist());del a
        print('normalization stage2',r,flush=True)
    write(HERE/'normalization-profile.json',profile)

def condition(index):
    profile=read(HERE/'normalization-profile.json')
    assert profile['condition_sha256']==sha(__file__)
    if index==1:
        assert read(HERE/'field-0/qualification.json')['qualified'],'Source zero not qualified'
        import subprocess
        assert subprocess.check_output(['git','show','HEAD:outputs/barseq3-native-repair-oct05/normalization-profile.json'],cwd=ROOT)==(HERE/'normalization-profile.json').read_bytes()
    geo=geometry(index);box=geo['bbox_y0_y1_x0_x1'];origin=geo['roi_origin_yx'][0]
    y=origin-box[0];x=origin-box[2];dest=HERE/f'field-{index}'
    fixed=sitk.GetImageFromArray(normalized(index,0,geo,profile).sum(axis=0))
    reference=sitk.Image([128,128],sitk.sitkFloat64);reference.SetOrigin((float(x),float(y)))
    for r in range(10):
        output=CACHE/f'native-repair-source-{index}-round-{r}.npz'
        if output.exists():continue
        start=time.monotonic();a=normalized(index,r,geo,profile)
        if r==0:
            transform=sitk.Transform(2,sitk.sitkIdentity);stats=dict(identity=True)
        else:
            demons=sitk.DemonsRegistrationFilter();demons.SetNumberOfIterations(150);demons.SetStandardDeviations(2.)
            displacement=demons.Execute(fixed,sitk.GetImageFromArray(a.sum(axis=0)))
            field=sitk.GetArrayViewFromImage(displacement)
            stats=dict(iterations=demons.GetElapsedIterations(),rms_change=demons.GetRMSChange(),
                full_field_sha256=hashlib.sha256(field.tobytes()).hexdigest(),roi_displacement=field[y:y+128,x:x+128].tolist())
            del field
            transform=sitk.DisplacementFieldTransform(displacement)
        aligned=[]
        for c in range(4):
            im=sitk.GetImageFromArray(a[c]);aligned.append(sitk.GetArrayFromImage(sitk.Resample(im,reference,transform,sitk.sitkLinear,0.,im.GetPixelID())))
        t=read(dest/f'round-{r}-geometry.json')
        support=crop(processed(index,r,np.array(t['camera']),np.array(t['affine'],np.float32),True),box)
        supports=[]
        for c in range(4):
            im=sitk.GetImageFromArray(support[c]);supports.append(sitk.GetArrayFromImage(sitk.Resample(im,reference,transform,sitk.sitkLinear,0.,im.GetPixelID())))
        np.savez_compressed(output,images=np.stack(aligned),support=np.stack(supports))
        stats.update(elapsed_s=time.monotonic()-start,output_sha256=sha(output),condition_sha256=sha(__file__),normalization_profile_sha256=sha(HERE/'normalization-profile.json'))
        write(dest/f'round-{r}-demons.json',stats)
        print('conditioned',index,r,stats['elapsed_s'],flush=True)
        del a,aligned,support,supports,transform;gc.collect()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['controls','geometry','profile','condition']);p.add_argument('--index',type=int,default=0,choices=[0,1]);args=p.parse_args()
    try:
        if args.stage=='controls':controls()
        elif args.stage=='geometry':geometry(args.index)
        elif args.stage=='profile':profile_fit()
        else:condition(args.index)
    except Exception:
        write(HERE/f'{args.stage}-{args.index}-failure.json',dict(condition_sha256=sha(__file__),traceback=traceback.format_exc()))
        raise
