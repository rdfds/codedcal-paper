"""Stock native qualification; no threshold search and no recovery inference."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL','2')
os.environ.setdefault('TF_NUM_INTRAOP_THREADS','2')
os.environ.setdefault('TF_NUM_INTEROP_THREADS','1')
import argparse
import csv
import gc
import hashlib
import sys
import time
import traceback
from pathlib import Path
import numpy as np
from scipy.ndimage import maximum_filter
from condition import CACHE,HERE,ROOT,read,write,sha
sys.path.insert(0,str(ROOT.parent/'bardensr-upstream.nosync'))
import bardensr

def codebook():
    rows=list(csv.DictReader((HERE.parent/'independent-comparison-oct05/barseq3-adult-codebook.csv').open()))
    groups={}
    for r in rows:groups.setdefault(r['Barcode'],set()).add(r['Name'])
    names=['|'.join(sorted(v)) for v in groups.values()]
    blank=np.array([all('no_gene' in n.lower() for n in v) for v in groups.values()])
    book=np.zeros((10,4,len(groups)),float)
    for j,s in enumerate(groups):
        for r,c in enumerate(s):book[r,'GTCA'.index(c),j]=1
    write(HERE/'codebook-classes.json',dict(names=names,sequences=list(groups),blank=blank.tolist(),
        codebook_sha256=sha(HERE.parent/'independent-comparison-oct05/barseq3-adult-codebook.csv')))
    return book.reshape(40,-1),blank,names

def execute(index,decoder):
    dest=HERE/f'field-{index}';record=dest/f'{decoder}-qualification.json'
    assert read(HERE/'background-controls.json')['condition_sha256']==sha(HERE/'condition.py')
    if index==1:assert read(HERE/'field-0/qualification.json')['qualified']
    book,blank,names=codebook();images=[];supports=[]
    for r in range(10):
        path=CACHE/f'native-repair-source-{index}-round-{r}.npz'
        marker=read(dest/f'round-{r}-demons.json')
        assert marker['output_sha256']==sha(path) and marker['condition_sha256']==sha(HERE/'condition.py')
        with np.load(path) as data:images.extend(data['images']);supports.extend(data['support'])
    data=np.stack(images)[:,None];support=np.stack(supports)[:,24:104,24:104]
    support_fraction=float(np.all(support>=1-1e-6,axis=0).mean());del images,supports,support;gc.collect()
    threshold=.775 if decoder=='correlation' else .20
    try:
        start=time.monotonic()
        if decoder=='correlation':density=bardensr.spot_calling.estimate_density_singleshot(data,book,.05)
        else:density,_=bardensr.spot_calling.estimate_density_iterative(data.astype(float),book,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
        density=np.asarray(density,np.float32);assert np.isfinite(density).all();kernel_s=time.monotonic()-start
        cut=np.float32(threshold*float(density.max()) if decoder=='regression' else threshold)
        maxima=maximum_filter(density,size=(1,3,3,1),mode='constant',cval=np.finfo(np.float32).min)
        z,y,x,j=np.where((density==maxima)&(maxima>cut));inside=(x>=24)&(x<104)&(y>=24)&(y<104)
        x,y,j=x[inside],y[inside],j[inside]
        order=np.lexsort((x,y,j));x,y,j=x[order],y[order],j[order]
        nonblank=int((~blank[j]).sum());nblank=int(blank[j].sum());proxy=nblank/max(nonblank,1)*int((~blank).sum())/int(blank.sum())
        out=dict(status='completed',source_index=index,decoder=decoder,nonblank_calls=nonblank,blank_calls=nblank,
            blank_proxy=proxy,support_fraction=support_fraction,qualified=nonblank>=50 and proxy<=.10 and support_fraction>=.99,
            xy=np.c_[x,y].tolist(),words=j.tolist(),threshold=threshold,absolute_threshold=float(cut),
            density_max=float(density.max()),density_sha256=hashlib.sha256(density.tobytes()).hexdigest(),
            ambiguous_calls=int((j==names.index('Lamc1|Rex2')).sum()),kernel_s=kernel_s,
            qualify_sha256=sha(__file__),condition_sha256=sha(HERE/'condition.py'),
            normalization_profile_sha256=sha(HERE/'normalization-profile.json'),
            conclusion='Native qualification only; blank proxy is not molecular truth')
    except Exception:out=dict(status='failed',qualified=False,traceback=traceback.format_exc(),source_index=index,decoder=decoder,qualify_sha256=sha(__file__))
    write(record,out);print(index,decoder,out.get('nonblank_calls'),out.get('blank_proxy'),out['qualified'],flush=True)
    records=[dest/f'{d}-qualification.json' for d in ['correlation','regression']]
    if all(p.exists() for p in records):write(dest/'qualification.json',dict(qualified=all(read(p)['qualified'] for p in records),decoders={p.name.split('-')[0]:read(p) for p in records}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,choices=[0,1],default=0);p.add_argument('--decoder',choices=['correlation','regression'],required=True);a=p.parse_args();execute(a.index,a.decoder)
