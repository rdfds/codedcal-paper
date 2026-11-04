"""Offline native-context proof from receipts, metadata and retained integer cells."""
from pathlib import Path
import json,hashlib,struct
import numpy as np
H=Path(__file__).resolve().parent;C=H.parents[1].parent/'workflow-qualification-fullfield.nosync'
s=json.loads((H/'fullfield-context-receipts.json').read_text());assert s['stage']=='ready_for_native_conditioning' and s['requests']<=32 and s['transport_bytes']<=134217728
originals=json.loads((H/'fresh-collection-receipts.json').read_text())['records'];report=[]
for e in [r for r in s['records'] if r.get('stage')=='extracted']:
 fid=e['id'];old=next(r for r in originals if r['id']==fid);assert old['publisher_md5_verified'] and old['actual_md5']==old['md5']==e['source_whole_md5_verified_previously']
 header=(H/'fullfield-metadata'/f'{fid}-header.bin').read_bytes();tail=(H/'fullfield-metadata'/f'{fid}-tail.bin').read_bytes();size=old['bytes'];parts=[(0,header),(size-len(tail),tail)]
 def read(offset,n):
  for start,body in parts:
   if start<=offset and offset+n<=start+len(body):return body[offset-start:offset-start+n]
  raise ValueError('Metadata pointer outside retained coverage')
 endian='<' if header[:2]==b'II' else ('>' if header[:2]==b'MM' else None);assert endian and struct.unpack(endian+'H',header[2:4])[0]==42
 def u16(offset):return struct.unpack(endian+'H',read(offset,2))[0]
 def u32(offset):return struct.unpack(endian+'I',read(offset,4))[0]
 offset=u32(4);seen=set();pages=[];tag_ids=set()
 while offset:
  assert offset not in seen;seen.add(offset);n=u16(offset);tags={}
  for i in range(n):
   pos=offset+2+12*i;tag,kind,count=struct.unpack(endian+'HHI',read(pos,8));tag_ids.add(tag)
   assert tag in {256,257,258,259,262,270,273,277,278,279,282,283,296,305,339},('Unexpected tag',tag)
   width={1:1,2:1,3:2,4:4,5:8}.get(kind);assert width is not None
   blob=read(pos+8,4)[:width*count] if width*count<=4 else read(u32(pos+8),width*count)
   if kind in (3,4) and count==1:tags[tag]=struct.unpack(endian+('H' if kind==3 else 'I'),blob)[0]
  assert tags[256]==tags[257]==1608 and tags[258]==16 and tags[259]==1 and tags[277]==1 and tags[279]==5171328
  pages.append(tags);offset=u32(offset+2+12*n)
 assert len(pages)==45
 path=Path(e['path']);assert hashlib.sha256(path.read_bytes()).hexdigest()==e['sha256'];dtype=np.dtype(endian+'u2')
 def wire(a):return a[::-1,::-1].T.astype(dtype,copy=False).tobytes(order='C')
 with np.load(path) as a:
  bodies={'fixed native fiducial frame6':wire(a['fiducial_frame6']),'fixed RNA z7 frames22/23':wire(a['frame22'])+wire(a['frame23'])}
 for record in [r for r in s['records'] if r.get('id')==fid and r.get('stage')=='validated']:
  assert record['status']==206 and record['headers']['ETag']==old['headers']['ETag']
  assert record['headers']['Content-Range']==f'bytes {record["start"]}-{record["end"]}/{size}'
  if record['reason']=='native TIFF header':body=header
  elif record['reason']=='all native IFD metadata':body=tail
  else:body=bodies[record['reason']]
  assert len(body)==record['end']-record['start']+1==record['bytes'] and hashlib.sha256(body).hexdigest()==record['sha256']
  if 'fiducial' in record['reason']:assert record['start']==pages[6][273]
  if 'RNA' in record['reason']:assert record['start']==pages[22][273] and pages[23][273]==pages[22][273]+5171328
 report.append(dict(file_id=fid,all45_ifds_checked=True,tag_ids=sorted(tag_ids),retained_native_cells_reproduce_range_hashes=True,source_whole_md5_verified_at_collection=True))
assert len(report)==8
out=dict(passed=True,records=report,requests_including_failed_conditional_request=s['requests'],received_response_body_bytes=s['transport_bytes'],whole_TIFF_hashes_recomputed_offline=False,metadata_and_selected_pixel_ranges_replayed=True,physical_validation_added=False)
(H/'native-context-verification.json').write_text(json.dumps(out,indent=2)+'\n');print('VERIFIED',len(report),'native contexts; all selected-wire hashes reproduced')
