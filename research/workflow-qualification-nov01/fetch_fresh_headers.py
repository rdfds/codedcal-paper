from pathlib import Path
import hashlib,json,shutil,time,urllib.request,urllib.error
HERE=Path(__file__).resolve().parent;OUT=HERE/'fresh-headers';OUT.mkdir(exist_ok=True)
assert not (OUT/'receipts.json').exists()
entries=json.loads((HERE/'public-metadata/metadata-1.txt').read_text())['data']['latestVersion']['files']
receipts=[]
for i,fid in enumerate([12555510,12555511,12555509,12555507,12555727]):
    if i:time.sleep(10)
    assert shutil.disk_usage(OUT).free>=1073741824
    e=next(x for x in entries if x['dataFile']['id']==fid);assert e['restricted'] is False;f=e['dataFile']
    url='https://dataverse.harvard.edu/api/access/datafile/'+str(fid);headers={'User-Agent':'codedcal-public-qualification/0.1'}
    if i==4:headers['Range']='bytes=0-65535'
    row=dict(url=url,file_id=fid,publisher=f,range_requested=i==4)
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=40) as r:
            row.update(status=r.status,headers={k:v for k,v in r.headers.items() if k.lower() in ['content-length','content-range','content-type','content-disposition','etag','last-modified','retry-after']})
            cap=65536 if i==4 else 524288;body=r.read(cap+1)
            assert len(body)<=cap,'Body exceeds cap'
            if i==4:assert r.status==206 and r.headers.get('Content-Range')=='bytes 0-65535/232717320' and len(body)==65536
            else:assert r.status==200
            row['success']=True
    except urllib.error.HTTPError as x:row.update(status=x.code,error=str(x),retry_after=x.headers.get('Retry-After'),success=False);body=x.read(8192)
    except Exception as x:row.update(status=row.get('status'),error=type(x).__name__+': '+str(x),success=False);body=locals().get('body',b'')[:524288]
    p=OUT/(str(fid)+'.bin');p.write_bytes(body);row.update(file=p.name,bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),publisher_md5_matches=hashlib.md5(body).hexdigest()==f['checksum']['value'] if i<4 else False)
    receipts.append(row);(OUT/'receipts.json').write_text(json.dumps(dict(records=receipts,full_tiff_hash_verified=False,pixel_outcomes_examined=False),indent=2)+'\n');print(fid,row['status'],len(body),row['success'],row['publisher_md5_matches'],flush=True)
    if row['status']==429 or row.get('retry_after') or row.get('headers',{}).get('Retry-After'):break
