"""Tests for codedcal.certificate against cases whose answer is known."""
import itertools
from pathlib import Path
import numpy as np
import pytest

from codedcal import certificate as C


def test_rm24_is_blind():
    rec = C.certify('RM(2,4)', C.rm24_shell())
    assert rec['verdict'] == 'C'


def test_random_per_source_is_conditional_not_empirically_verified():
    rec = C.certify('rand', np.eye(4, dtype=int), random_per_source=True)
    assert rec['verdict'] == 'CONDITIONAL'
    assert not rec['assumptions_verified']
    assert not rec['budget_contrasts_certified']


def test_no_shared_pairs_is_class_N():
    # weight-1 words: no plane pair shares a target
    assert C.certify('onehot', np.eye(6, dtype=int))['verdict'] == 'N'


def test_partition_layout_has_invariants():
    # seqFISH-style one-channel-per-round layout, every code of 4 rounds x 3
    # channels: structurally rank-deficient (certificate-sep18)
    codes = list(itertools.product(range(3), repeat=4))
    rec = C.certify('partition', C.partition_code(codes, 4, 3))
    assert rec['informative_invariants'] > 0


def test_exact_rank_matches_float_rank():
    rng = np.random.default_rng(7)
    M = (rng.random((60, 40)) < 0.15).astype(int)
    r = np.linalg.matrix_rank(M.astype(float))
    assert C.rank_mod_p_fast(M, C.Q1) == C.rank_mod_p_fast(M, C.Q2) == r


def test_weak_projection_is_not_reported_as_exact_blindness():
    root = Path(__file__).resolve().parents[2]
    _, W = C.load_codebook_csv(root / 'outputs/code-design-sep16/codebook_VA142.csv')
    rec = C.certify('VA142', W, replicated_planes=[0])
    assert rec['nuisance_free_dimension'] == 6
    assert rec['verdict'] == 'C'
    assert rec['node_frac']['max'] > 0
    assert rec['exact_blindness'] is None
    assert rec['classification_basis'] == 'below practical projection threshold'
    assert not rec['replication_assessed']


def test_blank_policy_and_binary_metadata_ids(tmp_path):
    path = tmp_path / 'book.csv'
    path.write_text('name,id,b0,b1,b2\nGeneA,0,1,1,0\nGeneB,1,1,0,1\nBlank1,0,0,1,1\n')
    names, W = C.load_codebook_csv(path)
    assert names == ['GeneA', 'GeneB']
    assert W.shape == (2, 3)  # id is metadata, never an extra bit
    _, all_codes = C.load_codebook_csv(path, include_blanks=True)
    assert all_codes.shape == (3, 3)
    assert C.certify('genes', W)['covered_pairs'] == 2
    assert C.certify('all', all_codes)['covered_pairs'] == 3


@pytest.mark.parametrize('text', ['name,id,b0\nA,0,2\n', 'name,id,b0\nA,0\n',
                                  'name,id,b0\nBlank,0,1\n'])
def test_malformed_or_empty_codebook_fails(text, tmp_path):
    path = tmp_path / 'book.csv'
    path.write_text(text)
    with pytest.raises(ValueError):
        C.load_codebook_csv(path)
