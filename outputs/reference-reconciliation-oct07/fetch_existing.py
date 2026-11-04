"""Recover just the 24 previously analyzed CD20 picked members, exactly."""
import fcntl
import importlib.util
import shutil
import subprocess
from common import CACHE, P, ROOT, cd20_registry, historical_members, read, save, sha

def main():
    CACHE.mkdir(exist_ok=True)
    spec = importlib.util.spec_from_file_location(
        'verified_ranges', ROOT / 'outputs/shaw-sensitivity-sep26/fetch_inputs.py')
    ranges = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ranges)
    idx = read('outputs/random-code-reference-search-sep13/resi2023_inventory.json')
    lookup = {e['name']: e for e in idx['files']}
    old = historical_members()
    manifest_path = P / 'input-manifest.json'
    manifest = read(str(manifest_path.relative_to(ROOT))) if manifest_path.exists() else {}
    head = CACHE / 'head.bin'
    ranges.fetch(idx['url'], 0, 29, head)
    unrar = shutil.which('unrar')
    assert unrar, 'unrar required for independent extraction / CRC'
    for row in cd20_registry():
        tag = row['tag']
        for r in range(1, 5):
            name = row['prefix'] + f'RESI/RESI rounds/RESI_R{r}_apicked.hdf5'
            entry = lookup[name]
            prior = old[name]
            expected = next(v['sha256'] for k, v in prior.items()
                            if k in ('unrar', 'unar', 'bsdtar') and 'sha256' in v)
            dest = CACHE / 'raw' / name
            if dest.exists() and dest.stat().st_size == entry['size'] and sha(dest) == expected:
                print('verified cached', tag, r, flush=True)
                continue
            segment = CACHE / f'{tag}-{r}.range'
            lo, hi = entry['header_offset'], entry['data_offset'] + entry['compressed']
            print('fetch', tag, r, round((hi-lo)/1e6, 1), 'MB', flush=True)
            ranges.fetch(idx['url'], lo, hi, segment)
            archive = CACHE / 'member.rar'
            with archive.open('wb') as f:
                f.write(head.read_bytes())
                with segment.open('rb') as g:
                    shutil.copyfileobj(g, f, 1 << 22)
                f.write(bytes.fromhex('1d77565103050400'))
            raw = CACHE / 'raw'; raw.mkdir(exist_ok=True)
            duplicate = CACHE / 'duplicate'; duplicate.mkdir(exist_ok=True)
            b = subprocess.run(['/usr/bin/bsdtar', '-x', '-f', str(archive), '-C', str(duplicate), name], capture_output=True)
            u = subprocess.run([unrar, 'x', '-y', '-o+', str(archive), name, str(raw) + '/'], capture_output=True)
            t = subprocess.run([unrar, 't', str(archive)], capture_output=True, text=True)
            copy = duplicate / name
            assert dest.exists() and copy.exists(), (name, b.stderr, u.stderr)
            assert dest.stat().st_size == entry['size'] == copy.stat().st_size
            assert sha(dest) == expected == sha(copy), name
            assert any(line.strip().startswith('Testing') and name in line and line.strip().endswith('OK')
                       for line in t.stdout.splitlines()), (name, t.stdout)
            manifest[name] = dict(field=tag, round=r, sha256=expected, bytes=entry['size'],
                                  archive_url=idx['url'], range=[lo, hi],
                                  duplicate_sha256=sha(copy), crc=True,
                                  path=str(dest), original_hash_matched=True)
            save(manifest_path, manifest)
            copy.unlink(); archive.unlink(); segment.unlink()
            print('retained verified', tag, r, flush=True)

if __name__ == '__main__':
    CACHE.mkdir(exist_ok=True)
    with (CACHE / 'fetch.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        main()
