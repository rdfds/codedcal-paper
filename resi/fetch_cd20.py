"""Fetch RESI round localizations (CD20 views and the Nup96 NPC view) from the public
archive, one member at a time.

Each member is retrieved with an HTTP range request, wrapped as a single-member RAR5
archive and extracted with unrar, which checks the member CRC. Only the fields used
downstream are kept, in nm, in a compact npz per view and round:

    x, y (float32, nm), frame (int32), lp (float32, nm; sqrt((lpx^2 + lpy^2) / 2))

RESI molecule tables (one per round) are kept as x, y (nm), lp (nm), n.

Usage: python fetch_cd20.py CACHE_DIR [view ...]
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
PIXEL_NM = 130.0
END_MARKER = bytes.fromhex('1d77565103050400')


def fetch(url, lo, hi, dest):
    """Bytes [lo, hi) of url into dest."""
    subprocess.run(['curl', '-sSfL', '--retry', '5', '-r', f'{lo}-{hi - 1}', '-o', str(dest), url],
                   check=True)
    assert dest.stat().st_size == hi - lo, (dest, dest.stat().st_size, hi - lo)


def extract(url, head, member, work):
    seg = work / 'segment.bin'
    fetch(url, member['header_offset'], member['end_offset'], seg)
    archive = work / 'member.rar'
    with archive.open('wb') as f:
        f.write(head)
        with seg.open('rb') as g:
            shutil.copyfileobj(g, f, 1 << 22)
        f.write(END_MARKER)
    seg.unlink()
    out = work / 'x'
    out.mkdir(exist_ok=True)
    r = subprocess.run(['unrar', 'x', '-y', '-o+', str(archive), member['name'], str(out) + '/'],
                       capture_output=True, text=True)
    archive.unlink()
    path = out / member['name']
    # The synthetic end-of-archive block makes unrar report a trailing-archive error;
    # the member itself must extract with a passing CRC ("OK") at the full size.
    lines = [l for l in r.stdout.splitlines() if l.startswith('Extracting') and member['name'] in l]
    assert path.exists() and lines and lines[-1].rstrip().endswith('OK'), (
        member['name'], r.stdout[-500:], r.stderr[-500:])
    assert path.stat().st_size == member['unpacked']
    return path


def compact_locs(path):
    with h5py.File(path, 'r') as f:
        a = f['locs'][:]
    lp = PIXEL_NM * np.sqrt((a['lpx'].astype(float) ** 2 + a['lpy'].astype(float) ** 2) / 2)
    out = dict(x=(a['x'] * PIXEL_NM).astype(np.float32), y=(a['y'] * PIXEL_NM).astype(np.float32),
               frame=a['frame'].astype(np.int32), lp=lp.astype(np.float32))
    if 'group' in a.dtype.names:          # pick (e.g. one nuclear pore) when present
        out['group'] = a['group'].astype(np.int32)
    return out


def compact_resi(path):
    with h5py.File(path, 'r') as f:
        a = f['locs'][:]
    out = dict(x=a['x'] * PIXEL_NM, y=a['y'] * PIXEL_NM,
               lp=PIXEL_NM * np.sqrt((a['lpx'] ** 2 + a['lpy'] ** 2) / 2))
    for key in ('n', 'group'):
        if key in a.dtype.names:
            out[key] = a[key]
    return out


def main(cache, views):
    manifest = json.loads((HERE / 'cd20_members.json').read_text())
    url = manifest['source'] + '?download=1'
    cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=cache) as tmp:
        work = Path(tmp)
        head_path = work / 'head.bin'
        fetch(url, 0, manifest['main_header_end'], head_path)
        head = head_path.read_bytes()
        assert head[:8] == b'Rar!\x1a\x07\x01\x00'
        for view in views:
            members = manifest['views'][view]
            if 'info' in members and not (cache / f'{view}-information.txt').exists():
                path = extract(url, head, members['info'], work)
                (cache / f'{view}-information.txt').write_bytes(path.read_bytes())
                path.unlink()
            for kind, convert in (('resi', compact_resi), ('locs', compact_locs)):
                for r, member in sorted(members[kind].items()):
                    dest = cache / f'{view}-{kind}-R{r}.npz'
                    if dest.exists():
                        print('cached', dest.name, flush=True)
                        continue
                    path = extract(url, head, member, work)
                    np.savez_compressed(dest, **convert(path))
                    path.unlink()
                    print('wrote', dest.name, round(dest.stat().st_size / 1e6, 1), 'MB', flush=True)


if __name__ == '__main__':
    manifest_views = list(json.loads((HERE / 'cd20_members.json').read_text())['views'])
    main(Path(sys.argv[1]), sys.argv[2:] or manifest_views)
