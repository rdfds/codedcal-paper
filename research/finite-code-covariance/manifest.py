"""Hash closure for this standalone theory extension; excludes manifest itself."""
from pathlib import Path
import argparse,hashlib,json
P=Path(__file__).resolve().parent

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
target=P/'ARTIFACTS.json'
if args.verify:
    d=json.loads(target.read_text());names=set()
    for r in d['files']:
        path=Path(r['path'])
        if path.is_absolute() or '..' in path.parts or path.as_posix()!=r['path']:raise ValueError('unsafe artifact path')
        f=P/path
        if f.is_symlink() or not f.is_file() or f.stat().st_size!=r['bytes'] or sha(f)!=r['sha256']:raise ValueError('artifact mismatch: '+r['path'])
        if r['path'] in names:raise ValueError('duplicate path')
        names.add(r['path'])
    assert names=={f.relative_to(P).as_posix() for f in P.rglob('*') if f.is_file() and f!=target}
    print(json.dumps(dict(passed=True,files=len(names),bytes=sum(r['bytes'] for r in d['files']))))
else:
    rows=[]
    for f in sorted(P.rglob('*')):
        if f.is_symlink():raise ValueError('linked artifact')
        if f.is_file() and f!=target:rows.append(dict(path=f.relative_to(P).as_posix(),bytes=f.stat().st_size,sha256=sha(f)))
    target.write_text(json.dumps(dict(format=1,files=rows),indent=2)+'\n')
    print(json.dumps(dict(written=True,files=len(rows))))
