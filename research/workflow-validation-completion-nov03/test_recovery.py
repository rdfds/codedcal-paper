"""Offline safety/control tests; no native workflow executions or fits."""
import copy
import json
from pathlib import Path
import unittest
from recover import H, P, Stop, missing_pairs, restore_cache, validate_dispatch


class RecoveryControls(unittest.TestCase):
    def setUp(self):
        self.f = json.loads((P/'freeze.json').read_text())
        self.pairs = missing_pairs(self.f, P/'run')
        self.calls = [json.loads(p.read_text()) for p in sorted((P/'run/calls').glob('*.json'))]
        self.rule = dict(maximum_new_dispatches=6,
            missing_pairs=[dict(policy=p, challenge=c) for p, c in self.pairs],
            retry=dict(input_fingerprint=self.calls[-1]['input_fingerprint']))

    def test_exact_missing_frozen_probes_only(self):
        self.assertEqual(len(self.pairs), 3)
        self.assertTrue(all(p=='group_search' and ch['role']=='probe' for p,ch in self.pairs))
        self.assertEqual([p+'/'+ch['id'] for p,ch in self.pairs], json.loads((P/'summary.json').read_text())['missing'])

    def test_interrupted_output_never_becomes_cache(self):
        cache = restore_cache(self.calls)
        self.assertNotIn(self.calls[-1]['input_fingerprint'], cache)
        self.assertGreater(len(cache), 0)

    def test_conflicting_completed_cache_is_rejected(self):
        row = next(c for c in self.calls if c['stage']=='completed' and c['policy']=='group_search')
        changed = copy.deepcopy(row)
        changed['coverage'] = row['coverage']-.01
        with self.assertRaises(ValueError):
            restore_cache([row, changed])

    def test_completed_duplicates_allow_timing_difference_only(self):
        row = next(c for c in self.calls if c['stage']=='completed' and c['policy']=='group_search')
        other = copy.deepcopy(row)
        other['diagnostics']['wall_seconds'] = 123456.
        self.assertEqual(len(restore_cache([row, other])), 1)

    def test_changed_retry_input_is_rejected(self):
        with self.assertRaises(Stop):
            validate_dispatch(self.rule, 0, self.pairs[0][1]['id'], 'wrong-input', False)

    def test_uncached_reference_and_extra_call_are_rejected(self):
        with self.assertRaises(Stop):
            validate_dispatch(self.rule, 1, self.pairs[0][1]['id'], 'missing-cache', False)
        with self.assertRaises(Stop):
            validate_dispatch(self.rule, 6, self.pairs[-1][1]['id'], 'extra', True)

    def test_only_original_six_call_sequence_is_allowed(self):
        for i in range(6):
            key = self.rule['retry']['input_fingerprint'] if i==0 else str(i)
            validate_dispatch(self.rule, i, self.pairs[i//2][1]['id'], key, bool(i%2))
        with self.assertRaises(Stop):
            validate_dispatch(self.rule, 2, self.pairs[0][1]['id'], 'wrong-order', False)


if __name__ == '__main__':
    unittest.main()
