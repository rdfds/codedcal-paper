from pathlib import Path
import json,sys,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(R/'tool/src'),str(R.parent/'serval-upstream.nosync')]
from codedcal.workflow_validation import write_json
from serval.image import ImageStack
from serval.transform import FiducialAlignmentImageTransform,HighPassImageTransform,DeconvoleImageTransform,LowPassImageTransform
from serval.decode.utils import get_imgs_hist,get_init_scaling_factors
import skimage.transform
C=R.parent/'workflow-qualification-fresh.nosync';entry=json.loads((H/'fresh-input.json').read_text());old=np.load(entry['path']);s=json.loads((C/'state.json').read_text());raw={};fids=[]
for r in s['records']:
 p=Path(r['extraction']);assert hashlib.sha256(p.read_bytes()).hexdigest()==r['extraction_sha256']
 with np.load(p) as a:
  fids.append(a['fiducial_frame6'])
  for k in a.files:
   if k.startswith('bit'):raw[int(k[3:])]=a[k]
raw=np.stack([raw[i] for i in range(1,17)]);rr=old['rounds'];cc=old['colors'];native_fid=FiducialAlignmentImageTransform.get_fiducial_transforms(ImageStack(np.stack(fids)),sigma=3,upsample_factor=100)
transform=FiducialAlignmentImageTransform({i:int(r) for i,r in enumerate(rr)},{0:native_fid})
transforms=[transform,HighPassImageTransform(sigma=2),DeconvoleImageTransform(filter_size=9,num_iters=10,sigma=1.2),LowPassImageTransform(sigma=1.)]
stack=ImageStack(raw,fov=0,z=7);hist=get_imgs_hist(stack,transforms);scales=get_init_scaling_factors(hist,q=.9)
processed=stack
for t in transforms:processed=t.transform(processed,preserve_dtype=False)
images=(processed.imgs/scales[:,None,None])[:,None].astype(np.float32)
support=np.stack([skimage.transform.warp(np.ones(raw[0].shape,np.float32),native_fid[int(r)],preserve_range=True) for r in rr])[:,None].astype(np.float32)
metadata=json.loads(str(old['metadata_json'].item()));metadata['preprocessing']='actual_source_fiducial_highpass2_RL10_sigma1.2_lowpass1_hist90_initial_scale';metadata['native_adapter_repair']='one_pinned_source_conformant_repair_no_comparative_outcomes'
p=C/'source-conditioned-input.npz';assert not p.exists();data={k:old[k] for k in old.files if k not in ['images','support','metadata_json']};np.savez_compressed(p,images=images,support=support,metadata_json=json.dumps(metadata),**data)
entry.update(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),metadata=metadata,source_conditioning_scales=scales,
 source_fiducial_matrices=[t.params for t in native_fid],source_conditioning_sha256={str(p.relative_to(R.parent/'serval-upstream.nosync')):hashlib.sha256(p.read_bytes()).hexdigest() for p in (R.parent/'serval-upstream.nosync/serval').rglob('*.py')})
write_json(H/'fresh-input-v2.json',entry);print('PREPARED',images.shape,scales.tolist(),flush=True)
