"""Read-only integrity of this independent research snapshot; standard library only."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'provenance/files.json'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()

def paths():
 import os
 for base,dirs,files in os.walk(ROOT):
  dirs[:]=[d for d in dirs if d not in {'.git','__pycache__','.pytest_cache','build','dist'} and not d.startswith('.venv') and not d.endswith('.egg-info')]
  for name in files:
   p=Path(base)/name
   if p!=MANIFEST and not name.endswith(('.pyc','.pyo')) and name!='.DS_Store':yield p

def verify():
 expected=json.loads(MANIFEST.read_text())['files']
 names=set()
 for row in expected:
  rel=PurePosixPath(row['path'])
  if rel.is_absolute() or '..' in rel.parts or str(rel) in names:raise ValueError('Unsafe/repeated evidence path')
  names.add(str(rel));p=ROOT.joinpath(*rel.parts)
  if p.is_symlink() or any(q.is_symlink() for q in p.parents if q!=ROOT.parent):raise ValueError('Linked evidence path')
  if p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:raise ValueError('Altered or missing file: '+str(rel))
 if names!={p.relative_to(ROOT).as_posix() for p in paths()}:raise ValueError('Unsealed files or incomplete inventory')
 if (ROOT/'.git').is_file() or (ROOT/'.git/objects/info/alternates').exists():raise ValueError('External Git dependency')
 imported=json.loads((ROOT/'provenance/import.json').read_text())['files']
 for r in imported:
  if sha(ROOT/r['path'])!=r['sha256']:raise ValueError('Import ledger mismatch: '+r['path'])
 return dict(passed=True,files=len(expected),bytes=sum(r['bytes'] for r in expected),
             imported_files=len(imported),imported_byte_identical=sum(r['byte_identical'] for r in imported),
             external_history_required=False,downloads=0,new_outcome_study=False)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--seal',action='store_true');a=p.parse_args()
 if a.seal:
  rows=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths())]
  with MANIFEST.open('x') as f:json.dump(dict(files=rows,scope='Current research snapshot; manifest itself and runtime caches excluded.'),f,indent=2);f.write('\n')
 print(json.dumps(verify(),indent=2))
