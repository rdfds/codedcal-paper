from pathlib import Path
import sys,json,csv,io,hashlib
import numpy as np
from scipy.ndimage import gaussian_filter
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path.insert(0,str(R/'tool/src'))
from codedcal.workflow_bardensr import bounded_translation
from codedcal.workflow_validation import ImageCase,warp,write_json
C=R.parent/'workflow-qualification-fresh.nosync';s=json.loads((C/'state.json').read_text());assert s['stage']=='raw_input_retained'
assert len(s['records'])==9 and all(r['publisher_md5_verified'] for r in s['records'])
bookrows=list(csv.DictReader(io.StringIO((H/'fresh-originals/12555510.bin').read_text())))
org=list(csv.DictReader(io.StringIO((H/'fresh-originals/12555511.bin').read_text())));org=sorted([r for r in org if r.get('readoutName','').startswith('RS')],key=lambda r:int(r['bitNumber']))
labels=[r['gene_name'] for r in bookrows];ids=np.array([labels.index(n) for n in labels]);blanks=np.array(sorted({ids[i] for i,n in enumerate(labels) if n.lower().startswith(('blank','neg','unused'))}),int)
book=np.array([[float(row[r['readoutName']]) for row in bookrows] for r in org])
assert book.shape[0]==16 and np.all(book.sum(axis=0)==4)
fids=[];raw={}
for r in s['records']:
    p=Path(r['extraction']);assert hashlib.sha256(p.read_bytes()).hexdigest()==r['extraction_sha256']
    with np.load(p) as a:
        fids.append(a['fiducial_frame6'].astype(np.float32))
        for key in a.files:
            if key.startswith('bit'):raw[int(key[3:])]=a[key].astype(np.float32)
maps=np.stack([bounded_translation(fids[0]-np.mean(fids[0]),a-np.mean(a),(32,32)) for a in fids]);maps-=maps[:1]
images=np.stack([raw[i] for i in range(1,17)])[:,None];rr=np.array([int(r['imagingRound'])-1 for r in org]);cc=np.array([int(r['color']) for r in org]);corrections=np.c_[np.zeros(16),maps[rr]]
corrected=warp(images,corrections);support=warp(np.ones_like(images),corrections)
processed=[];scales=[]
for a in corrected[:,0]:
    b=gaussian_filter(np.maximum(a-gaussian_filter(a,3),0),.7);scale=float(np.percentile(b[b>0],99.9)) if np.any(b>0) else 0.;assert scale>0
    processed.append(b/scale);scales.append(scale)
images=np.stack(processed)[:,None].astype(np.float32)
metadata=dict(context=[[0,1],[16,240],[16,240]],scoring_region=[[0,1],[24,232],[24,232]],coverage_region=[[0,1],[24,232],[24,232]],
    target='image_localization_catalog',public_source='doi:10.7910/DVN/NE0WOP',field='000',z_index=7,preprocessing='pinned_fiducial_round_correction_and_bandpass',microns_per_pixel_reported=.024)
case=ImageCase(images,support,book,rr,cc,ids,blanks,'4t1-culture-000','4t1-culture-NE0WOP',metadata)
p=C/'input.npz';assert not p.exists();np.savez_compressed(p,images=images,support=support,codebook=book,rounds=rr,colors=cc,class_ids=ids,blank_classes=blanks,name=case.name,preparation=case.preparation,metadata_json=json.dumps(metadata))
write_json(H/'fresh-input.json',dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),shape=images.shape,class_names=labels,blank_classes=blanks,metadata=metadata,
    fiducial_round_maps=maps,scales=scales,original_records=[{k:r[k] for k in ['name','id','md5','actual_md5','sha256','bytes_received','extraction_sha256','mapping']} for r in s['records']],comparative_outcomes_examined=False))
print('PREPARED',images.shape,book.shape,'blanks',len(blanks),'support',np.all(support[:,:,24:232,24:232]>=1-1e-6,axis=0).mean(),flush=True)
