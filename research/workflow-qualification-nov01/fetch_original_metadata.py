from pathlib import Path
import hashlib,json,shutil,time,urllib.request,urllib.error
H=Path(__file__).resolve().parent;O=H/'fresh-originals';O.mkdir(exist_ok=True);assert not (O/'receipts.json').exists()
entries=json.loads((H/'public-metadata/metadata-1.txt').read_text())['data']['latestVersion']['files'];rows=[]
for i,fid in enumerate([12555510,12555511,12555507]):
    if i:time.sleep(10)
    assert shutil.disk_usage(O).free>=1073741824
    e=next(e for e in entries if e['dataFile']['id']==fid);assert not e['restricted'];f=e['dataFile']
    url=f'https://dataverse.harvard.edu/api/access/datafile/{fid}?format=original';r=dict(url=url,file_id=fid,expected_md5=f['checksum']['value']);body=b''
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'codedcal-public-qualification/0.1'}),timeout=40) as a:
            r.update(status=a.status,headers={k:v for k,v in a.headers.items() if k.lower() in ['content-length','content-type','etag','content-disposition','retry-after']});body=a.read(524289);assert len(body)<=524288
            r['success']=a.status==200
    except urllib.error.HTTPError as e:r.update(status=e.code,error=str(e),retry_after=e.headers.get('Retry-After'),success=False)
    except Exception as e:r.update(error=type(e).__name__+': '+str(e),success=False)
    (O/(str(fid)+'.bin')).write_bytes(body);r.update(bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),md5=hashlib.md5(body).hexdigest(),publisher_md5_matches=hashlib.md5(body).hexdigest()==f['checksum']['value']);rows.append(r)
    (O/'receipts.json').write_text(json.dumps({'records':rows},indent=2)+'\n');print(fid,r.get('status'),r['bytes'],r['publisher_md5_matches'],flush=True)
    if r.get('status')==429 or r.get('retry_after') or r.get('headers',{}).get('Retry-After'):break
