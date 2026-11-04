"""Finite member transport and header verification; no intensity evaluation."""
import gc,hashlib,importlib.util,io,json,shutil,time,zipfile,zlib
from pathlib import Path
import requests,tifffile
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent;AUDIT=ROOT/'outputs/barseq-public-substitute-oct06';CACHE=ROOT.parent/'starfinder-tissue-inputs-oct06.nosync'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while b:=f.read(1024**2):h.update(b)
    return h.hexdigest()
def write(path,obj):
    path.parent.mkdir(exist_ok=True,parents=True);temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2)+'\n');temp.replace(path)
def main():
    source=ROOT/'outputs/public-opportunities-sep27/inspect_archives.py';spec=importlib.util.spec_from_file_location('reader',source);reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    meta=json.loads((AUDIT/'11176779-metadata.json').read_text());archive=next(f for f in meta['files'] if f['key']=='sample-dataset.zip');remote=reader.Remote(archive['links']['self'],archive['size'])
    CACHE.mkdir(exist_ok=True);path=HERE/'tissue-retrieval.json';state=json.loads(path.read_text()) if path.exists() else dict(source='tissue-2D',field='tile_1',planned=17,completed=0,records=[],raw_intensity_pixels_examined=False,full_archive_md5_verified=False,archive_url=remote.url,archive_bytes=remote.size,script_sha256=sha(Path(__file__)),protocol_sha256=sha(HERE/'TISSUE_RETRIEVAL_PROTOCOL.md'))
    assert state['script_sha256']==sha(Path(__file__))
    with zipfile.ZipFile(remote) as z:
        entries=[x for x in z.infolist() if '/tissue-2D/round' in x.filename and '/tile_1/' in x.filename and x.filename.endswith('.tif')];assert len(entries)==17
        for item in entries:
            name=Path(item.filename).parts[2]+'__'+Path(item.filename).name;target=CACHE/(name+'.deflate');old=next((r for r in state['records'] if r['member']==item.filename),None)
            if old:assert sha(target)==old['compressed_sha256'];continue
            assert item.compress_type==8 and item.compress_size<100_000_000 and item.file_size==283498208
            assert shutil.disk_usage(CACHE).free>item.compress_size+item.file_size+1024**3
            remote.seek(item.header_offset);header=remote.read(30);assert header[:4]==b'PK\x03\x04';n=int.from_bytes(header[26:28],'little');extra=int.from_bytes(header[28:30],'little');start=item.header_offset+30+n+extra;end=start+item.compress_size-1
            url=remote.url+f'?audit_range={start}-{end}'
            begin=time.monotonic()
            with requests.get(url,headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(30,60)) as resp:
                assert resp.status_code==206 and resp.headers.get('Content-Range')==f'bytes {start}-{end}/{remote.size}'
                temp=target.with_suffix('.partial');h=hashlib.sha256();count=0
                with temp.open('wb') as f:
                    for chunk in resp.iter_content(4*1024**2):f.write(chunk);h.update(chunk);count+=len(chunk)
            assert count==item.compress_size;payload=zlib.decompress(temp.read_bytes(),-15);assert len(payload)==item.file_size and zlib.crc32(payload)==item.CRC
            with tifffile.TiffFile(io.BytesIO(payload)) as t:pages=[dict(shape=list(page.shape),dtype=str(page.dtype)) for page in t.pages]
            assert len(pages)==30 and all(p==dict(shape=[3072,3072],dtype='uint8') for p in pages)
            record=dict(member=item.filename,local_compressed_path=str(target.resolve()),compression='raw DEFLATE original ZIP member',compressed_bytes=count,compressed_sha256=h.hexdigest(),bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest(),crc32=f'{item.CRC:08x}',zip_crc32_verified=True,header_offset=item.header_offset,data_start=start,pages=30,page_shape=[3072,3072],dtype='uint8',pixels_decoded=False,retrieval_s=time.monotonic()-begin)
            temp.replace(target);state['records'].append(record);state.update(completed=len(state['records']),stage='complete' if len(state['records'])==17 else 'partial');write(path,state);print('TISSUE MEMBER',len(state['records']),item.filename,count,flush=True);del payload;gc.collect()
    print('COMPLETE',state['completed'],'source members; intensity pixels unopened',flush=True)
if __name__=='__main__':main()
