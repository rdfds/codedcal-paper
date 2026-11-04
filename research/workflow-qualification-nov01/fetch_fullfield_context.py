"""Strict finite native-context re-read of already whole-checksummed TIFFs."""
from pathlib import Path
import io,os,json,time,hashlib,fcntl,shutil,urllib.request,urllib.error
import numpy as np
import tifffile
H=Path(__file__).resolve().parent;R=H.parents[1];C=R.parent/'workflow-qualification-fullfield.nosync';C.mkdir(exist_ok=True)
lock=(C/'collector.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
state=C/'state.json'
prior=json.loads(state.read_text())
assert prior['stage']=='stopped_with_error' and prior['requests']==1 and prior['transport_bytes']==0 and prior['records'][-1]['status']==412
assert (H/'fullfield-context-precondition-stop.json').read_bytes()==state.read_bytes()
assert shutil.disk_usage(C).free>=536870912+125829120
state.rename(C/'state-precondition-stop.json')
pins=json.loads((H/'fresh-pixel-pins.json').read_text())[:8];old=json.loads((H/'fresh-collection-receipts.json').read_text())['records']
for p in pins:assert next(r for r in old if r['id']==p['id'])['publisher_md5_verified']
rows=prior['records'];started=time.time();last=prior['started_epoch'];transport=prior['transport_bytes'];requests=prior['requests']
for p in pins:p['etag']=next(r for r in old if r['id']==p['id'])['headers']['ETag']

def save(stage,error=None):
 d=dict(pid=os.getpid(),command=str(Path(__file__).resolve()),stage=stage,requests=requests,transport_bytes=transport,records=rows,started_epoch=started,error=error,comparative_outcomes_examined=False)
 p=state.with_suffix('.tmp');p.write_text(json.dumps(d,indent=2)+'\n');p.replace(state)

def fetch(p,start,end,reason):
 global last,transport,requests
 assert 0<=start<=end<p['bytes'] and end-start+1<=16777216
 assert requests<32 and transport+end-start+1<=134217728 and time.time()-started<1800
 assert shutil.disk_usage(C).free>=536870912+16777216
 if last:time.sleep(max(0,10-(time.time()-last)))
 row=dict(id=p['id'],url=p['url'],start=start,end=end,reason=reason,expected_etag=p['etag'],stage='pinned');rows.append(row);save('running')
 requests+=1;last=time.time();body=b''
 try:
  req=urllib.request.Request(p['url'],headers={'Range':f'bytes={start}-{end}','User-Agent':'codedcal-native-context/0.1','If-Match':p['etag']})
  with urllib.request.urlopen(req,timeout=40) as r:
   row.update(status=r.status,headers={k:v for k,v in r.headers.items() if k.lower() in ['content-range','content-length','etag','content-type','last-modified','retry-after']})
   assert r.status==206 and r.headers.get('Content-Range')==f'bytes {start}-{end}/{p["bytes"]}'
   assert not r.headers.get('Retry-After') and r.headers.get('ETag','')==p['etag']
   body=r.read(end-start+2);transport+=len(body)
   assert len(body)==end-start+1
  row.update(stage='validated',bytes=len(body),sha256=hashlib.sha256(body).hexdigest());save('running');return body
 except urllib.error.HTTPError as e:
  row.update(status=e.code,retry_after=e.headers.get('Retry-After'),error=str(e));save('stopped_with_error',str(e));raise
 except BaseException as e:row.update(error=type(e).__name__+': '+str(e),bytes_received=len(body));save('stopped_with_error',row['error']);raise

class Sparse(io.RawIOBase):
 def __init__(self,n,parts):self.n=n;self.parts=parts;self.pos=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else (self.pos+offset if whence==1 else self.n+offset)
  if not 0<=self.pos<=self.n:raise ValueError('Seek out of bounds')
  return self.pos
 def read(self,n=-1):
  if n<0:n=self.n-self.pos
  if self.pos+n>self.n:n=self.n-self.pos
  for start,body in self.parts:
   if start<=self.pos and self.pos+n<=start+len(body):
    out=body[self.pos-start:self.pos-start+n];self.pos+=n;return out
  raise ValueError(f'Unpinned TIFF metadata read at{self.pos} size{n}')

save('running')
try:
 for p in pins:
  if p['round']==1:
   prefix=H/'fresh-headers/12555727.bin';receipt=json.loads((H/'fresh-headers/receipts.json').read_text())['records'][-1]
   assert hashlib.sha256(prefix.read_bytes()).hexdigest()==receipt['sha256'] and receipt['headers']['ETag']==p['etag']
   header=prefix.read_bytes()[:512];rows.append(dict(stage='verified_header_cache_reuse',id=p['id'],bytes=512,source_sha256=receipt['sha256'],etag=p['etag']))
  else:header=fetch(p,0,511,'native TIFF header')
  tailstart=p['bytes']-16384;tail=fetch(p,tailstart,p['bytes']-1,'all native IFD metadata')
  (C/f'{p["id"]}-header.bin').write_bytes(header);(C/f'{p["id"]}-tail.bin').write_bytes(tail)
  with tifffile.TiffFile(Sparse(p['bytes'],[(0,header),(tailstart,tail)])) as t:
   assert len(t.pages)==45;info={}
   for f in [6,22,23]:
    page=t.pages[f];assert page.shape==(1608,1608) and page.dtype==np.dtype('uint16') and int(page.compression)==1
    assert len(page.dataoffsets)==len(page.databytecounts)==1 and page.databytecounts[0]==5171328
    info[f]=(int(page.dataoffsets[0]),int(page.databytecounts[0]))
   assert info[23][0]==sum(info[22]);endian=t.byteorder
  (C/f'{p["id"]}-ifds.json').write_text(json.dumps(dict(info=info,endian=endian),indent=2)+'\n')
  f=fetch(p,info[6][0],sum(info[6])-1,'fixed native fiducial frame6')
  rna=fetch(p,info[22][0],sum(info[23])-1,'fixed RNA z7 frames22/23')
  dtype=np.dtype(endian+'u2');arrays={}
  arrays['fiducial_frame6']=np.frombuffer(f,dtype).reshape(1608,1608).T[::-1,::-1].copy()
  for k,frame in enumerate([22,23]):arrays['frame'+str(frame)]=np.frombuffer(rna[k*5171328:(k+1)*5171328],dtype).reshape(1608,1608).T[::-1,::-1].copy()
  path=C/f'round-{p["round"]:02}.npz';np.savez_compressed(path,**arrays)
  rows.append(dict(stage='extracted',id=p['id'],round=p['round'],path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_whole_md5_verified_previously=p['md5']))
  save('running');print(p['round'],'context retained',requests,transport,flush=True)
 save('ready_for_native_conditioning')
except BaseException as e:save('stopped_with_error',type(e).__name__+': '+str(e));raise
