"""Read public remote HDF5 schema only, with an enforced transport cap.

Never calls Dataset.__getitem__: image/annotation values remain unopened.
All fetched ranges must be exact 206 responses; no fallback whole downloads.
"""
import io
import json
from pathlib import Path
import requests
import h5py

HERE=Path(__file__).resolve().parent
meta=json.loads((HERE/'barseq3-zenodo.json').read_text())['body']
file=next(f for f in meta['files'] if f['key'].startswith('adult_cb_raw'))

class Ranges(io.RawIOBase):
    def __init__(self,url,size,cap=8*1024**2,block=32768):
        self.url=url;self.size=size;self.cap=cap;self.block=block
        self.position=0;self.blocks={};self.bytes=0;self.ledger=[]
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset
        if self.position<0:raise ValueError('negative seek')
        return self.position
    def readinto(self,b):
        data=self.read(len(b));b[:len(data)]=data;return len(data)
    def read(self,n=-1):
        if n<0:raise ValueError('unbounded read forbidden')
        if self.position>=self.size:return b''
        n=min(n,self.size-self.position)
        output=[]
        end=self.position+n
        while self.position<end:
            base=self.position//self.block*self.block
            if base not in self.blocks:
                last=min(base+self.block,self.size)-1
                if self.bytes+last-base+1>self.cap:raise ValueError('metadata transport cap exceeded')
                r=requests.get(self.url,params={'schema_range':base},
                    headers={'Range':f'bytes={base}-{last}'},timeout=30)
                assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {base}-{last}/{self.size}'
                assert len(r.content)==last-base+1
                self.blocks[base]=r.content;self.bytes+=len(r.content)
                self.ledger.append(dict(start=base,end=last,bytes=len(r.content)))
            stop=min(end,base+len(self.blocks[base]))
            output.append(self.blocks[base][self.position-base:stop-base]);self.position=stop
        return b''.join(output)

remote=Ranges(file['links']['self'],file['size'])
schema={}
with h5py.File(remote,'r') as f:
    attrs={k:str(v) for k,v in f.attrs.items()}
    def visit(name,obj):
        if isinstance(obj,h5py.Dataset):
            schema[name]=dict(shape=list(obj.shape),dtype=str(obj.dtype),
                chunks=list(obj.chunks) if obj.chunks else None,compression=obj.compression)
    f.visititems(visit)
result=dict(file=file,attrs=attrs,datasets=schema,transport_bytes=remote.bytes,
            range_ledger=remote.ledger,exposure='HDF5 schema only; no dataset values read',
            date_caveat='API created June 2026; publication_date says March 2028. Publication status not inferred.')
(HERE/'barseq3-remote-schema.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['datasets','transport_bytes','exposure']},indent=2))
