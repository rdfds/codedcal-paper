"""Losslessly retain verified source TIFFs before removing expanded duplicates."""
import hashlib,json,zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent;CACHE=ROOT.parent/'starfinder-inputs-oct06.nosync'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024**2):h.update(b)
 return h.hexdigest()
def main():
 source=ROOT/'outputs/barseq-public-substitute-oct06/starfinder-culture-retrieval.json';manifest=json.loads(source.read_text());target=CACHE/'Position351-source-lossless.zip';partial=target.with_suffix('.partial')
 assert len(manifest['records'])==25
 if not target.exists():
  with zipfile.ZipFile(partial,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
   for r in manifest['records']:
    path=Path(r['local_path']);assert path.stat().st_size==r['bytes'] and sha(path)==r['sha256'];z.write(path,arcname=r['member']);print('PACKED',r['member'],flush=True)
  partial.replace(target)
 verified=[]
 with zipfile.ZipFile(target) as z:
  assert sorted(z.namelist())==sorted(r['member'] for r in manifest['records'])
  for r in manifest['records']:
   item=z.getinfo(r['member']);assert item.file_size==r['bytes'] and item.CRC==int(r['crc32'],16)
   h=hashlib.sha256();total=0
   with z.open(item) as f:
    while b:=f.read(1024**2):h.update(b);total+=len(b)
   assert total==r['bytes'] and h.hexdigest()==r['sha256'];verified.append(dict(original_path=r['local_path'],zip_member=r['member'],bytes=total,sha256=h.hexdigest(),crc32_verified=True))
 record=dict(stage='fully_verified_before_duplicate_removal',archive_path=str(target.resolve()),archive_bytes=target.stat().st_size,archive_sha256=sha(target),original_manifest_sha256=sha(source),records=verified,all_original_source_bytes_preserved=True,original_archive_md5_not_verified=True)
 out=HERE/'culture-packed-storage.json';out.write_text(json.dumps(record,indent=2)+'\n')
 for r in manifest['records']:
  path=Path(r['local_path'])
  if path.exists():assert sha(path)==r['sha256'];path.unlink()
 record['stage']='complete';record['expanded_duplicates_removed']=25;record['net_bytes_freed']=sum(r['bytes'] for r in manifest['records'])-target.stat().st_size;out.write_text(json.dumps(record,indent=2)+'\n');print('FREED',record['net_bytes_freed'],'bytes; all25original payloads losslessly retained',flush=True)
if __name__=='__main__':main()
