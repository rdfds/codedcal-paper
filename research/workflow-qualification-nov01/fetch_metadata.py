"""Finite, pre-pinned public metadata reader; no images or automatic retries."""
from pathlib import Path
import hashlib,json,shutil,time,urllib.request,urllib.error
BASE=Path(__file__).resolve().parent;plan=json.loads((BASE/'PUBLIC_INPUT_PLAN.json').read_text())
out=BASE/'public-metadata';out.mkdir(exist_ok=True)
assert not (out/'receipts.json').exists(),'Do not overwrite an access attempt'
rows=[]
for i,url in enumerate(plan['urls']):
    if i:time.sleep(plan['minimum_spacing_seconds'])
    assert shutil.disk_usage(BASE).free>=plan['required_disk_reserve_bytes']
    row=dict(url=url,cap=plan['max_body_bytes_per_request'],index=i)
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'codedcal-public-metadata/0.1'}),timeout=plan['timeout_seconds']) as r:
            row.update(status=r.status,content_type=r.headers.get('Content-Type'),retry_after=r.headers.get('Retry-After'),final_url=r.url)
            body=r.read(plan['max_body_bytes_per_request']+1)
            if len(body)>plan['max_body_bytes_per_request']:raise ValueError('Metadata body exceeded pinned cap')
    except urllib.error.HTTPError as e:
        row.update(status=e.code,retry_after=e.headers.get('Retry-After'),error=str(e));body=e.read(8192)
    except Exception as e:row.update(status=None,error=type(e).__name__+': '+str(e));body=b''
    p=out/f'metadata-{i}.txt';p.write_bytes(body);row.update(bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),file=p.name)
    rows.append(row);(out/'receipts.json').write_text(json.dumps(dict(records=rows,network_calls=len(rows),raw_images_read=False),indent=2)+'\n')
    print(i,row['status'],len(body),row.get('error',''),flush=True)
    if row['status']==429 or row.get('retry_after'):
        (out/'STOP.json').write_text(json.dumps(row,indent=2)+'\n');break
