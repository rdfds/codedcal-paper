"""Read single members of a remote zip archive with HTTP range requests (via curl)."""
import io
import subprocess
import zipfile


def _curl(args):
    return subprocess.run(['curl', '-sSfL', '--retry', '5'] + args, check=True, capture_output=True).stdout


class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only view of a remote file served with Range support."""

    def __init__(self, url):
        # some hosts refuse HEAD; a one-byte ranged GET reports the total size
        head = _curl(['-r', '0-0', '-D', '-', '-o', '/dev/null', url]).decode(errors='replace').lower()
        sizes = [int(l.rsplit('/', 1)[1]) for l in head.splitlines() if l.startswith('content-range:')]
        self.url, self.size, self.pos = url, sizes[-1], 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=io.SEEK_SET):
        self.pos = {io.SEEK_SET: offset, io.SEEK_CUR: self.pos + offset, io.SEEK_END: self.size + offset}[whence]
        return self.pos

    def readinto(self, b):
        if self.pos >= self.size:
            return 0
        end = min(self.pos + len(b), self.size) - 1
        data = _curl(['-r', f'{self.pos}-{end}', self.url])
        b[:len(data)] = data
        self.pos += len(data)
        return len(data)


def open_remote_zip(url, buffer=1 << 22):
    return zipfile.ZipFile(io.BufferedReader(HTTPRangeFile(url), buffer_size=buffer))
