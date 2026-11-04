"""Native full-plane conditioning, streamed before fixed ROI extraction."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(R/'tool/src'),str(R.parent/'serval-upstream.nosync')]
from codedcal.workflow_validation import write_json
from serval.image import ImageStack
from serval.transform import FiducialAlignmentImageTransform,HighPassImageTransform,DeconvoleImageTransform,LowPassImageTransform
from serval.decode.utils import get_init_scaling_factors
import skimage.transform


def stream_bit(a,bit,transforms,preserve_dtype):
    dtype=a.dtype
    for t in transforms:
        a=t._transform(a,bit,0,7)
        if preserve_dtype:a=a.astype(dtype)
    return a


def main(size=256):
 assert size in (256,512)
 low=(1608-size)//2;high=low+size
 C=R.parent/'workflow-qualification-fullfield.nosync';state=json.loads((C/'state.json').read_text());assert state['stage']=='ready_for_native_conditioning'
 rows=[r for r in state['records'] if r.get('stage')=='extracted'];assert len(rows)==8
 original=json.loads((H/'fresh-input.json').read_text());source=np.load(original['path']);old=json.loads((H/'fresh-collection-receipts.json').read_text())['records'];mapping={};fids=[]
 for row in rows:
  p=Path(row['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']
  previous=next(r for r in old if r['round']==row['round'])
  with np.load(p) as a,np.load(H/'fresh-public-patches'/f'round-{row["round"]:02}.npz') as b:
   np.testing.assert_array_equal(a['fiducial_frame6'][676:932,676:932],b['fiducial_frame6']);fids.append(a['fiducial_frame6'])
   for entry in previous['mapping']:
    key='frame'+str(entry['frame_index']);np.testing.assert_array_equal(a[key][676:932,676:932],b['bit'+str(entry['bit'])])
    mapping[entry['bit']]=(p,key)
 native_fid=FiducialAlignmentImageTransform.get_fiducial_transforms(ImageStack(np.stack(fids)),sigma=3,upsample_factor=100)
 del fids
 rounds=source['rounds'];transforms=[FiducialAlignmentImageTransform({i:int(r) for i,r in enumerate(rounds)},{0:native_fid}),HighPassImageTransform(sigma=2),DeconvoleImageTransform(filter_size=9,num_iters=10,sigma=1.2),LowPassImageTransform(sigma=1.)]
 # Source histogram path quantizes each transform; prediction path preserves float.
 hist=[];patches=[];supports=[]
 for bit in range(16):
  p,key=mapping[bit+1]
  with np.load(p) as a:plane=a[key]
  quantized=stream_bit(plane,bit,transforms,True)
  hist.append(np.histogram(quantized,bins=np.arange(np.iinfo(plane.dtype).max))[0]);del quantized
  conditioned=stream_bit(plane,bit,transforms,False)
  patches.append(conditioned[low:high,low:high].astype(np.float32).copy());del conditioned
  alpha=skimage.transform.warp(np.ones(plane.shape,np.float32),native_fid[int(rounds[bit])],preserve_range=True)
  supports.append(alpha[low:high,low:high].copy());print('CONDITIONED',bit+1,flush=True)
 scales=get_init_scaling_factors(np.array(hist),q=.9);images=(np.stack(patches)/scales[:,None,None])[:,None].astype(np.float32);support=np.stack(supports)[:,None].astype(np.float32)
 metadata=json.loads(str(source['metadata_json'].item()));metadata['preprocessing']='actual_fullfield_source_fiducial_HP2_RL10_sigma1.2_LP1_hist90_then_fixed_crop';metadata['native_adapter_repair']='full_native_context_fixed_ROI_and_criteria'
 for key in ['context','scoring_region','coverage_region']:
  margin=16 if key=='context' else 24;metadata[key]=[[0,1],[margin,size-margin],[margin,size-margin]]
 metadata['patch_size']=size;metadata['patch_origin_yx']=[low,low]
 p=C/('input.npz' if size==256 else 'input-512.npz');assert not p.exists();data={k:source[k] for k in source.files if k not in ['images','support','metadata_json']};np.savez_compressed(p,images=images,support=support,metadata_json=json.dumps(metadata),**data)
 original['source_conditioning_sha256']={str(v.relative_to(R.parent/'serval-upstream.nosync')):hashlib.sha256(v.read_bytes()).hexdigest() for v in (R.parent/'serval-upstream.nosync/serval').rglob('*.py')}
 original.update(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),metadata=metadata,source_conditioning_scales=scales,source_fiducial_matrices=[t.params for t in native_fid],raw_fixed_patch_cells_exactly_reproduced=True,source_fullfield_state_sha256=hashlib.sha256((C/'state.json').read_bytes()).hexdigest())
 write_json(H/('fresh-input-fullfield.json' if size==256 else 'fresh-input-fullfield-512.json'),original);print('PREPARED',images.shape,scales.tolist(),flush=True)


if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--size',type=int,choices=[256,512],default=256);args=parser.parse_args();main(args.size)
