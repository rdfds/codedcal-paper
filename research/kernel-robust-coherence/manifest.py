"""Seal/check this isolated research extension, excluding only the manifest itself."""
from pathlib import Path
import argparse,hashlib,json

BASE=Path(__file__).resolve().parent
TARGET=BASE/'ARTIFACTS.json'


def entries():
    return {str(p.relative_to(BASE)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(BASE.rglob('*')) if p.is_file() and p!=TARGET and '__pycache__' not in p.parts}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    actual=entries()
    if a.verify:
        stored=json.loads(TARGET.read_text())['files']
        if actual!=stored:
            changed=[k for k in set(actual)|set(stored) if actual.get(k)!=stored.get(k)]
            raise ValueError(f'artifact mismatch: {sorted(changed)}')
        print(f'Verified {len(actual)} stable extension files; integrity does not establish physical validation.')
    else:
        if TARGET.exists():raise FileExistsError('sealed manifest must not be overwritten')
        TARGET.write_text(json.dumps(dict(scope='Conditional research extension; no physical claim upgrade',files=actual),indent=2)+'\n')
        print(f'Sealed {len(actual)} stable extension files.')

if __name__=='__main__':main()
