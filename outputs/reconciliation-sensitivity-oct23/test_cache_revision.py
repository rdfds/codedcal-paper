"""The execution repair changes array loading, not scientific calculations."""
from pathlib import Path
import numpy as np


def test_performance_revision_only_changes_loading_and_metadata():
    p = Path(__file__).parent
    old = (p/'sensitivity-initial.py').read_text()
    old = old.replace("    empirical = np.load(OLD/f'empirical-{tag}.npz')\n",
        "    empirical = np.load(OLD/f'empirical-{tag}.npz')\n"
        "    # NpzFile __getitem__ decompresses the whole member on every call. Keep\n"
        "    # one in-memory array instead of rereading millions of events per group.\n"
        "    empirical_events = empirical['events']\n")
    old = old.replace("old_events = empirical['events'][lo:hi]", "old_events = empirical_events[lo:hi]")
    old = old.replace("    result = dict(tag=tag, distinct_acquisition=tag != 'sdc_raw', arms=results,\n",
        "    result = dict(tag=tag, distinct_acquisition=tag != 'sdc_raw', arms=results,\n"
        "        implementation_sha256=sha(__file__),\n")
    assert old == (p/'sensitivity.py').read_text()


def test_cached_and_reloaded_npz_member_slices_are_identical(tmp_path):
    rng = np.random.default_rng(26102318)
    source = rng.normal(size=(400, 8)); target = tmp_path/'events.npz'
    np.savez_compressed(target, events=source)
    with np.load(target) as f:
        cached = f['events']
        for lo, hi in ((0, 13), (13, 200), (200, 400), (300, 300)):
            assert np.array_equal(cached[lo:hi], f['events'][lo:hi])
