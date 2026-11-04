"""Bounded observed-data sensitivity; all historical evidence is read-only."""
import argparse
import inspect
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.stats import t as student

P = Path(__file__).resolve().parent
ROOT = P.parents[1]
OLD = ROOT/'outputs/reference-reconciliation-oct07'
sys.path.insert(0, str(OLD))
sys.path.insert(0, str(ROOT/'tool/src'))
from common import save, sha
from reference_characterize import (registry, load_round, assign_panels, phase_sums,
    diagonal, count_weight, ANGLES, Q, O, components, split_summary)
from simulate_reference import fast_clusters
from codedcal.uncertainty import attenuation_uncertainty, fieller_set

PROFILES = ('lenient', 'original', 'stringent')
RADII = (.8, 1., 1.2)
WEIGHTS = ('unweighted', 'event_count', 'localization_count')


def settings(row):
    default = {k: v.default for k, v in inspect.signature(O.source_qc).parameters.items()
               if v.default is not inspect.Parameter.empty}
    default.update(row.get('frozen', {}).get('oracle_qc', {}))
    if row['tag'] != 'npc':
        default['max_zsd'] = None
    result = []
    for radius in RADII:
        for profile in PROFILES:
            qc = dict(default)
            if profile != 'original':
                loose = profile == 'lenient'
                qc['min_events'] = int(np.ceil(default['min_events']*(.8 if loose else 1.2)))
                qc['min_per_panel'] = max(1, default['min_per_panel']+(-1 if loose else 1))
                for key in ('max_rms', 'max_zsd'):
                    if default[key] is not None:
                        qc[key] = default[key]*(1.2 if loose else .8)
                qc['min_sep'] = default['min_sep']*(.8 if loose else 1.2)
                if default['max_split'] is not None:
                    qc['max_split'] = min(1., default['max_split']+(.05 if loose else -.05))
            result.append(dict(name=f'r{radius:g}-{profile}', radius_factor=radius,
                radius_nm=radius*row.get('frozen', {}).get('cluster_thresh_nm', 10.),
                profile=profile, qc=qc))
    return result


def event_phase_sums(xy, panel, event, event_count):
    """Individual-row panels remain distinct from mean-event eligibility panels."""
    n = np.zeros((event_count, 4))
    z = np.zeros((event_count, 4, len(Q), len(ANGLES)), complex)
    projection = xy[:, 0, None]*np.cos(ANGLES)+xy[:, 1, None]*np.sin(ANGLES)
    for t in range(4):
        m = panel == t
        n[:, t] = np.bincount(event[m], minlength=event_count)
        for j, q in enumerate(Q):
            phase = np.exp(1j*q*projection[m])
            for k in range(len(ANGLES)):
                z[:, t, j, k] = (np.bincount(event[m], weights=phase[:, k].real,
                    minlength=event_count)+1j*np.bincount(event[m], weights=phase[:, k].imag,
                    minlength=event_count))
    return z, n


def code_components(z, n):
    return sum(np.stack([components(z[:, :, :, j, k], n)
        for j in range(len(Q))], 1) for k in range(len(ANGLES)))/len(ANGLES)


def tiles(centers):
    xy = centers[:, :2]; lo = xy.min(0); extent = np.maximum(xy.max(0)-lo, 1.)
    cell = np.minimum(((xy-lo)/extent*[4, 3]).astype(int), [3, 2])
    return cell[:, 0]*3+cell[:, 1]


def blocks(a, tile):
    return np.stack([a[tile == k].sum(0) for k in np.unique(tile)])


def signed_t(a):
    """Last two axes are window/component; preserve nonpositive support failures."""
    num = a[..., 0].mean(-1); den = a[..., 1].mean(-1)
    ratio = np.divide(num, den, out=np.full_like(num, np.nan), where=den != 0)
    valid = (ratio > 0) & (den > 0) & np.isfinite(ratio)
    return np.where(valid, -np.log(np.where(valid, ratio, 1.))/Q**2, np.nan)


def describe(comp, tile):
    b = blocks(comp, tile)
    uncertainty = attenuation_uncertainty(b[..., :2].mean(-2), Q)
    summary = split_summary(comp) if np.any(comp) else None
    loo = signed_t(b.sum(0)[None]-b)
    return dict(fit=summary, uncertainty=uncertainty, point_T=signed_t(b.sum(0)),
        leave_one_block_out_T=loo)


def paired(left, right, tile):
    a, b = blocks(left, tile), blocks(right, tile)
    rng = np.random.default_rng(26102361)
    ii = rng.integers(len(a), size=(2000, len(a)))
    draw = signed_t(a[ii].sum(1))-signed_t(b[ii].sum(1))
    point = signed_t(a.sum(0))-signed_t(b.sum(0))
    result = []
    for j, q in enumerate(Q):
        good = np.isfinite(draw[:, j])
        result.append(dict(q=q, difference=point[j], valid_draws=int(good.sum()),
            total_draws=len(good), conditional_pointwise_interval=(
                np.quantile(draw[good, j], [.025, .975]) if good.mean() >= .95 else None)))
    return result


def analyze(row):
    tag = row['tag']; started = time.time(); arms = settings(row)
    original = np.load(OLD/f'reference-{tag}.npz'); centers = original['centers']; ng = len(centers)
    empirical = np.load(OLD/f'empirical-{tag}.npz')
    lookup = {(int(r), int(p)): (int(lo), int(hi)) for (r, p, _), lo, hi in
        zip(empirical['group_index'], empirical['offsets'][:-1], empirical['offsets'][1:])}
    selected_z = np.zeros((9, ng, 4, 4, len(Q), len(ANGLES)), complex)
    selected_n = np.zeros((9, ng, 4, 4))
    diag = np.zeros((9, ng, len(Q), 2, 4))
    quadratic = np.zeros((9, ng, 3, 2)); weight_squares = np.zeros((9, 3))
    counts = np.zeros(9, int); selected_events = np.zeros(9, int)
    selected_local = np.zeros(9, int); overlap = np.zeros((9, 2), int)
    rejects = [dict() for _ in arms]; all_events = 0; all_local = 0
    inputs = []; baseline_source_records = []
    baseline_idx = next(i for i, a in enumerate(arms) if a['name'] == 'r1-original')
    for r in range(1, 5):
        d = load_round(row, r, centers)
        inputs.append(dict(round=r, path=d['path'], sha256=d['sha256']))
        expected = json.loads((OLD/f'characterization-{tag}.json').read_text())['input_members'][r-1]
        assert d['sha256'] == expected['sha256'], (tag, r, 'raw hash')
        order = np.argsort(d['pick'], kind='stable')
        bounds = np.searchsorted(d['pick'][order], np.arange(ng+1))
        for p in range(ng):
            ii = order[bounds[p]:bounds[p+1]]
            xy = d['xy'][ii] if d.get('already_centered') else d['xy'][ii]-centers[p, :2]
            frame, lp = d['frame'][ii], d['lp'][ii]
            panel = assign_panels(frame, row['panels'][str(r)])
            all_local += int((panel >= 0).sum())
            if (r-1, p) not in lookup:
                continue
            ev = O.link_events(frame, xy, lp, jump=row.get('frozen', {}).get('link_jump_nm', 25.))
            zz = d['z'][ii] if d['z'] is not None else np.zeros(len(ii))
            pos, evf, nloc, var = O.event_table(frame, np.column_stack([xy, zz]), lp, ev)
            epan = assign_panels(evf, row['panels'][str(r)]); keep = epan >= 0
            ez, en = event_phase_sums(xy, panel, ev, len(pos))
            pos, var, epan, evf, nloc = pos[keep], var[keep], epan[keep], evf[keep], nloc[keep]
            ez, en = ez[keep], en[keep]
            lo, hi = lookup[(r-1, p)]; old_events = empirical['events'][lo:hi]
            assert np.allclose(np.column_stack([pos, epan, var, nloc, evf]),
                old_events[:, [0, 1, 2, 3, 4, 5, 7]], rtol=0, atol=1e-9), (tag, r, p, 'event replay')
            all_events += len(pos); masks = []
            labels = {radius: fast_clusters(pos, radius*row.get('frozen', {}).get('cluster_thresh_nm', 10.))
                      for radius in RADII}
            for i, arm in enumerate(arms):
                lab = labels[arm['radius_factor']]
                ok, reasons = O.source_qc(pos, lab, epan, **arm['qc'])
                for key, value in reasons.items():
                    rejects[i][key] = rejects[i].get(key, 0)+value
                records = O.oracle_stats(pos, var, lab, epan, ok, ANGLES)
                mask = np.isin(lab, [rec['source'] for rec in records]); masks.append(mask)
                selected_z[i, p, r-1] += ez[mask].sum(0)
                selected_n[i, p, r-1] += en[mask].sum(0)
                selected_events[i] += int(mask.sum()); selected_local[i] += int(en[mask].sum())
                for rec in records:
                    m = lab == rec['source']; ns = en[m].sum(0); zs = ez[m].sum(0)
                    diag[i, p] += diagonal(zs, ns)
                    weights = np.array([1., count_weight(rec['counts']), count_weight(ns)])
                    quadratic[i, p, :, 0] += weights*float(rec['that'].mean())
                    quadratic[i, p, :, 1] += weights
                    weight_squares[i] += weights**2
                    counts[i] += 1
                    if i == baseline_idx:
                        baseline_source_records.append((r-1, p, float(rec['that'].mean()), *weights[1:]))
            base_mask = masks[baseline_idx]
            for i, mask in enumerate(masks):
                overlap[i, 0] += int(en[mask & base_mask].sum())
                overlap[i, 1] += int(en[mask | base_mask].sum())
        print(tag, 'round', r, 'sources', counts.tolist(), flush=True)
    selected = np.stack([code_components(z, n) for z, n in zip(selected_z, selected_n)])
    baseline_errors = {k: float(np.max(abs(v-original[k]))) for k, v in (
        ('selected_code', selected[baseline_idx]), ('diagonal', diag[baseline_idx]))}
    for key, value in [('selected_code', selected[baseline_idx]), ('diagonal', diag[baseline_idx])]:
        assert np.allclose(value, original[key], rtol=1e-10, atol=1e-7), (tag, key, baseline_errors)
    base_records = np.array(baseline_source_records)
    # Sparse graph labels can reorder source rows: compare by group/round and moments.
    expected_rows = np.column_stack([original['round'], original['pore'], original['that'],
                                    original['event_countw'], original['localization_countw']])
    sort = lambda a: a[np.lexsort((a[:, 4], a[:, 3], a[:, 2], a[:, 1], a[:, 0]))]
    assert np.allclose(sort(base_records), sort(expected_rows), rtol=1e-10, atol=1e-8), (tag, 'sources')
    tile = tiles(centers); full = original['full_code']; full_desc = describe(full, tile)
    results = []
    for i, arm in enumerate(arms):
        sd, dd = describe(selected[i], tile), describe(diag[i], tile)
        sf = sd['fit']['total_arithmetic'] if sd['fit'] else None
        df = dd['fit']['total_arithmetic'] if dd['fit'] else None
        old_gap = abs(full_desc['fit']['total_arithmetic']-df) if df is not None else None
        new_gap = abs(sf-df) if sf is not None and df is not None else None
        weight_reports = {}
        for j, name in enumerate(WEIGHTS):
            b = blocks(quadratic[i, :, j], tile)
            point = b[:, 0].sum()/b[:, 1].sum() if b[:, 1].sum() else None
            weight_reports[name] = dict(point=point, uncertainty=fieller_set(b[:, 0], b[:, 1]),
                effective_sources=(b[:, 1].sum()**2/weight_squares[i, j] if weight_squares[i, j] else 0))
        results.append(dict(**arm, sources=int(counts[i]), rejection_counts=rejects[i],
            event_coverage=float(selected_events[i]/all_events),
            localization_coverage=float(selected_local[i]/all_local),
            original_row_jaccard=float(overlap[i, 0]/overlap[i, 1]) if overlap[i, 1] else None,
            selected_code=sd, diagonal=dd, quadratic=weight_reports,
            full_to_diagonal_fit_gap=old_gap, selected_to_diagonal_fit_gap=new_gap,
            fit_gap_narrows=(new_gap < old_gap if new_gap is not None else None),
            paired_selected_difference=paired(selected[i], diag[i], tile)))
    np.savez_compressed(P/f'groups-{tag}.npz', selected_code=selected, diagonal=diag,
        quadratic=quadratic, full_code=full, tile=tile, source_counts=counts,
        weight_squares=weight_squares, selected_localizations=selected_local,
        selected_events=selected_events, overlap=overlap, all_localizations=all_local,
        all_events=all_events)
    result = dict(tag=tag, distinct_acquisition=tag != 'sdc_raw', arms=results,
        full_code=full_desc, baseline_replay=dict(passed=True, errors=baseline_errors,
            source_rows=len(base_records)), input_members=inputs,
        old_reference_sha256=sha(OLD/f'reference-{tag}.npz'),
        old_empirical_sha256=sha(OLD/f'empirical-{tag}.npz'),
        groups_sha256=sha(P/f'groups-{tag}.npz'), runtime_seconds=time.time()-started)
    save(P/f'result-{tag}.json', result)
    print('finished', tag, 'seconds', round(result['runtime_seconds'], 1), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--tag', required=True)
    args = parser.parse_args()
    if (P/f'result-{args.tag}.json').exists():
        parser.error('Completed result exists; preserve it rather than overwrite')
    analyze(next(r for r in registry() if r['tag'] == args.tag))
