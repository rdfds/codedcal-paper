"""Running the published retained check must not mutate frozen evidence."""
from pathlib import Path
import subprocess
import sys


def test_default_verification_preserves_all_in_tree_file_bytes(tmp_path):
    p = Path(__file__).parent
    before = {f.name: f.read_bytes() for f in p.iterdir() if f.is_file()}
    result = subprocess.run([sys.executable, str(p/'verify.py')], cwd=tmp_path,
        capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    assert '"raw_inputs_freshly_verified": false' in result.stdout
    after = {f.name: f.read_bytes() for f in p.iterdir() if f.is_file()}
    assert before == after
