"""Codebook incidence and sensitivity diagnostics (revised 25 Sep 2026).

Given a codebook, report -- before any data are taken -- which components
of the cross-round registration error that design's own data can identify.
The theory is borrowed, not ours (novelty-blindcal-sep18): coded imaging
is a closure / blind-calibration system in which the CODEBOOK plays the
role of the array. Per-plane displacement enters pairwise (edge)
statistics; target-level nuisance enters through the edge-codeword
incidence N (N[e,g] = 1 iff both planes of edge e lie in codeword g's
support). The invariants of that nuisance are counted exactly as closure
invariants are: dim = covered edges - rank(N)  [Theorem B, corrected in
identifiability-theorem-sep16 Addendum 4; convention- and weight-
independent]. Per Lemma D, per-plane quantities are gauge-variant, so a
single-snapshot design can at best identify budget-plus-blur CONTRASTS.
Replication must be evaluated in an augmented observation model: its
presence alone does not certify separation for every plane.

Verdict classes:
  CONDITIONAL  random-source design declared; physical/model assumptions unchecked
  B  legacy class: median unweighted node projection >= 0.05
  C  legacy class: no nuisance-free invariants OR projection below 0.05;
     inspect exact_blindness and classification_basis to distinguish these
  N  no code redundancy (no plane pair shares a target): cross-plane
     information comes only from a shared reference channel
The projection is a diagnostic in the historical unweighted convention,
not an abundance-aware certificate of budget-contrast estimability.
Modular full row rank proves no invariants for the stated model. Deficient
ranks are cross-checked with SVD, not proved by an exact upper-bound witness.

Provenance: outputs/certificate-sep18/certificate.py (incidence, certify,
partition_code, rm24_shell) and fastrank.py (rank_mod_p_fast, Q1, Q2),
with revised claim semantics and blank-codeword handling.
"""
import csv
import itertools
import re

import numpy as np

Q1, Q2 = 2147483629, 2147483587          # primes < 2**31


def rank_mod_p_fast(M, p):
    A = np.asarray(M, dtype=np.int64) % p
    if A.shape[0] > A.shape[1]:
        A = A.T.copy()                    # eliminate along the short side
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0:
            continue
        piv = r + nz[0]
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        below = A[r + 1:, c].copy()
        idx = np.flatnonzero(below)
        if idx.size:
            A[r + 1 + idx] = (A[r + 1 + idx] -
                              (below[idx, None] * A[r][None, :]) % p) % p
        r += 1
        if r == rows:
            break
    return r


NODE_FRAC_MIN = 0.05       # historical practical threshold, not an exact zero
BLANK = re.compile(r'blank|empty|notarget|negprb|negative|systemcontrol|'
                   r'falsecode|^neg_|unassigned', re.I)


def incidence(W):
    W = np.asarray(W, int)
    n = W.shape[1]
    E = list(itertools.combinations(range(n), 2))
    N = np.array([[int(W[g, i] and W[g, j]) for g in range(W.shape[0])]
                  for (i, j) in E], dtype=int)
    return N, E


def certify(label, W, random_per_source=False, replicated_planes=(),
            notes=''):
    W = np.asarray(W)
    if W.ndim != 2 or min(W.shape) == 0 or not np.isin(W, [0, 1]).all():
        raise ValueError('codebook must be a nonempty binary target-by-plane matrix')
    W = W.astype(int)
    G, n = W.shape
    rec = dict(label=label, n_planes=int(n), n_targets=int(G),
               codeword_weights=sorted(set(int(x) for x in W.sum(1))),
               random_per_source=bool(random_per_source),
               replicated_planes=list(replicated_planes), notes=notes,
               replication_assessed=False, budget_contrasts_certified=False,
               model='covered-pair incidence with unrestricted target nuisance',
               exact_blindness=None)
    if random_per_source:
        rec.update(verdict='CONDITIONAL', assumptions_verified=False,
                   classification_basis='user-declared physical randomization',
                   reason='random-source identity is conditional on exogeneity, '
                   'panel stability and the temporal observation model; a flag '
                   'does not verify them or certify budget/blur separation')
        return rec
    N, E = incidence(W)
    cov = N.sum(1) > 0
    rec.update(pairs=len(E), covered_pairs=int(cov.sum()))
    if cov.sum() == 0:
        rec.update(verdict='N', informative_invariants=0,
                   nuisance_free_dimension=0,
                   classification_basis='no covered plane pairs',
                   reason='no plane pair shares a target; cross-plane '
                   'registration can only come from a reference channel')
        return rec
    Nk = N[cov]
    r1, r2 = rank_mod_p_fast(Nk, Q1), rank_mod_p_fast(Nk, Q2)
    s = np.linalg.svd(Nk.astype(float), compute_uv=False)
    rs = int((s > max(Nk.shape) * np.finfo(float).eps * s[0]).sum())
    if not (r1 == r2 == rs):
        rec.update(verdict='UNCERTIFIED', rank_mod_p=[r1, r2], rank_svd=rs)
        return rec
    k = int(cov.sum()) - rs
    rec.update(rank_N=rs, informative_invariants=k,
               nuisance_free_dimension=k,
               rank_mod_p=[r1, r2], rank_svd=rs,
               rank_evidence=('modular full row rank proves exact blindness'
                              if k == 0 else 'modular lower bounds agree with SVD; '
                              'no exact upper-bound witness'),
               exact_blindness=True if k == 0 else None,
               corollaryC_bound=max(0, int(cov.sum()) - G))
    frac = np.zeros(n)
    if k > 0:
        U, _, _ = np.linalg.svd(Nk.astype(float), full_matrices=True)
        K = U[:, rs:]                         # invariant edge functionals
        Ek = [e for e, c in zip(E, cov) if c]
        node = np.array([[1.0 if i in e else 0.0 for e in Ek]
                         for i in range(n)]).T
        frac = np.linalg.norm(K @ (K.T @ node), axis=0) / \
            np.maximum(np.linalg.norm(node, axis=0), 1e-12)
    rec['node_frac'] = dict(min=round(float(frac.min()), 4),
                            median=round(float(np.median(frac)), 4),
                            max=round(float(frac.max()), 4))
    rec['node_frac_threshold'] = NODE_FRAC_MIN
    rec['projection_convention'] = 'historical unweighted incidence and node signals'
    if k > 0 and np.median(frac) >= NODE_FRAC_MIN:
        rec['verdict'] = 'B'
        rec['classification_basis'] = 'above practical projection threshold'
        rec['reason'] = (f'{k} nuisance-free invariants; median unweighted node '
                         'projection meets the 0.05 sensitivity threshold; '
                         'specific budget contrasts are not certified')
    else:
        rec['verdict'] = 'C'
        rec['classification_basis'] = ('below practical projection threshold'
                                       if k > 0 else 'full row rank')
        rec['reason'] = (f'{k} nuisance-free invariants; median unweighted node '
                         'projection is below 0.05, which does not establish '
                         'exact non-identifiability' if k > 0 else
                         'no nuisance-free invariants: full row rank on covered '
                         'pairs under the stated target-nuisance model')
    if replicated_planes:
        rec['reason'] += '; replication declared but not modeled by this diagnostic'
    return rec


def partition_code(codes, n_rounds, n_channels):
    """Round x channel one-hot layout: codes[g][r] = channel of target g
    in round r (ISS / seqFISH-style); returns the binary plane matrix."""
    W = np.zeros((len(codes), n_rounds * n_channels), int)
    for g, cw in enumerate(codes):
        for r, c in enumerate(cw):
            W[g, r * n_channels + c] = 1
    return W


def rm24_shell():
    pts = np.array([[(q >> k) & 1 for k in range(4)] for q in range(16)])
    basis = ([np.ones(16, int)] + [pts[:, i] for i in range(4)] +
             [pts[:, i] * pts[:, j] for i in range(4)
              for j in range(i + 1, 4)])
    words = (np.array([[(m >> k) & 1 for k in range(11)]
                       for m in range(2048)]) @ np.array(basis)) % 2
    return words[words.sum(1) == 4]


def load_codebook_csv(path, include_blanks=False):
    """MERFISH-style codebook CSV (name, id, then one 0/1 column per bit):
    returns (names, W). By default the historical survey's name-based
    blank/control rule is used: unoccupied codes add no source nuisance.
    Use include_blanks=True only when those rows represent actual targets.
    Named id and barcodeType columns are metadata wherever they occur,
    even when all ids happen to be 0 or 1.
    Native MERlin space-separated barcode files are not this CSV format.
    """
    with open(path, encoding='utf-8-sig', newline='') as handle:
        rows = [r for r in csv.reader(handle) if r]
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError('expected CSV header name,[id,]bit... and target rows')
    head, body = rows[0], rows[1:]
    if any(len(r) != len(head) for r in body):
        raise ValueError('codebook rows must have the same number of columns')
    bitcols = [j for j, h in enumerate(head[1:], 1)
               if h.strip().lower() not in ('id', 'barcodetype')]
    if not bitcols or any(r[j].strip() not in ('0', '1') for r in body for j in bitcols):
        raise ValueError('non-metadata columns must contain binary bits')
    selected = [r for r in body if include_blanks or not BLANK.search(r[0])]
    if not selected:
        raise ValueError('no target rows remain after blank/control exclusion')
    names = [r[0] for r in selected]
    W = np.array([[int(r[j]) for j in bitcols] for r in selected], int)
    return names, W
