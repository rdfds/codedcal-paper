"""codedcal command line interface.

Commands:
  codedcal demo-chen2015          one-command reproduction: download the
                                  public Chen 2015 MERFISH example (17 MB),
                                  run the fiducial-free trace vs the bead
                                  oracle, print the comparison table.
  codedcal fish-trace <npz>       fiducial-free trace + color offsets + QC
                                  on a reduced per-plane spot/bead table.
  codedcal smlm-budget <table>    attenuation confidence sets + window-split
                                  D diagnostic from an SMLM event table.
  codedcal magnification <npz>    assignment-free per-colour magnification
                                  (pre-decode) with its quality gate;
                                  --round1 uses the first round only.
  codedcal certificate <csv>      design-time certificate: which cross-round
                                  error components a codebook can identify.
"""
import argparse
import hashlib
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

from . import __version__, certificate, chromatic, fish, retrofit, smlm

CHEN_URL = 'https://zhuanglab.s3.amazonaws.com/merfish/MERFISH_Examples.zip'
# The canonical S3 object intermittently returns 403; the Internet Archive
# holds a snapshot of the very same URL, tried second.
CHEN_MIRRORS = [
    CHEN_URL,
    'https://web.archive.org/web/20240813230954id_/' + CHEN_URL,
]
CHEN_SIZE_MB = 17
# SHA-256 of MERFISH_Examples.zip as validated with this package (a mismatch
# is reported as a warning, not an error, in case the object is re-uploaded).
CHEN_SHA256 = '0117048b8471b16f8aff8461ce1eaeaccff0bc7d361e9f7e1fbea5819f86751f'


def _json_safe(x):
    # From outputs/noise-equalization-sep14/test_equalization.py json_safe
    # (non-finite floats -> null rather than invalid JSON).
    if isinstance(x, dict):
        return {k: _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_safe(v) for v in x]
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def _write_json(path, obj):
    Path(path).write_text(json.dumps(_json_safe(obj), indent=2, default=float,
                                     allow_nan=False) + '\n')
    print(f'wrote {path}')


# ---------------------------------------------------------------- demo ----

def _find_example_dir(root):
    root = Path(root)
    if (root / 'STORM_00_01_c1_alist.bin').exists():
        return root
    for p in sorted(root.rglob('example_data')):
        if (p / 'STORM_00_01_c1_alist.bin').exists():
            return p
    return None


def cmd_demo(args):
    cache = Path(args.dir).expanduser()
    ex = _find_example_dir(cache)
    if ex is None:
        cache.mkdir(parents=True, exist_ok=True)
        zpath = cache / 'MERFISH_Examples.zip'
        if not zpath.exists():
            err = None
            for url in CHEN_MIRRORS:
                print(f'downloading {url} (~{CHEN_SIZE_MB} MB) -> {zpath}')
                try:
                    sha = hashlib.sha256()
                    with urllib.request.urlopen(url, timeout=300) as r, \
                            open(zpath, 'wb') as f:
                        while True:
                            chunk = r.read(1 << 20)
                            if not chunk:
                                break
                            f.write(chunk)
                            sha.update(chunk)
                            print('.', end='', flush=True)
                    print(' done')
                    print(f'sha256 {sha.hexdigest()}')
                    if CHEN_SHA256 and sha.hexdigest() != CHEN_SHA256:
                        print(f'WARNING: sha256 differs from the validated '
                              f'release ({CHEN_SHA256}); the upstream object '
                              'may have changed', file=sys.stderr)
                    err = None
                    break
                except Exception as e:  # noqa: BLE001
                    print(f'  failed: {e}', file=sys.stderr)
                    err = e
                    zpath.unlink(missing_ok=True)
            if err is not None:
                print('ERROR: could not download MERFISH_Examples.zip from '
                      'any mirror; download it manually and place it in '
                      f'{cache}', file=sys.stderr)
                return 1
        print(f'extracting {zpath}')
        with zipfile.ZipFile(zpath) as z:
            z.extractall(cache)
        ex = _find_example_dir(cache)
        if ex is None:
            print('ERROR: example_data with STORM_*_c1_alist.bin not found in the archive',
                  file=sys.stderr)
            return 1
    print(f'data: {ex}')
    out = fish.chen2015_trace(ex)
    print()
    print(out['dataset'])
    print('fiducial-free, assignment-free signal trace vs the dataset\'s own '
          'fiducial-bead oracle')
    print(f"bead pairs: >= {out['bead_pair_inliers']['min']} inliers/pair, "
          f"residual sd {out['bead_pair_inliers']['residual_sd_px']:.3f} px; "
          f"bleedthrough guard dropped {sum(out['bleedthrough_guard_dropped'])} signal spots")
    print()
    print('round  n_sig  n_bead    bead trace (px)     signal trace (px)      diff (px)')
    for r in range(len(out['n_signal'])):
        b = out['trace_bead_px'][r]; s = out['trace_signal_px'][r]; d = out['trace_diff_px'][r]
        print(f"{r:5d} {out['n_signal'][r]:6d} {out['n_bead'][r]:7d}   "
              f"{b[0]:7.3f} {b[1]:7.3f}    {s[0]:7.3f} {s[1]:7.3f}    "
              f"{d[0]:7.3f} {d[1]:7.3f}")
    print()
    print(f"RMS difference: {out['rms_diff_px']:.3f} px   "
          f"(max |diff| {out['max_abs_diff_px']:.3f} px; "
          f"corr x {out['corr_x']:.3f}, y {out['corr_y']:.3f})")
    print(f"pre-registered criterion {out['criterion']}: "
          f"{'PASS' if out['passed'] else 'FAIL'}   [expected ~0.17 px; "
          f"{out['px_nm_note']}]")
    if args.json:
        _write_json(args.json, out)
    return 0 if out['passed'] else 2


# ---------------------------------------------------------- fish-trace ----

def cmd_fish_trace(args):
    d = np.load(args.npz)
    if 'spots' not in d:
        print(f"ERROR: {args.npz} has no 'spots' array (expected keys "
              f"'spots' and optionally 'beads'; found {list(d.keys())})",
              file=sys.stderr)
        return 1
    spots = d['spots']
    beads = d['beads'] if 'beads' in d else np.empty((0, 3))
    excl = [tuple(int(v) for v in p.split(',')) for p in args.exclude_plane]
    out = fish.fov_trace(beads, spots, nr=args.rounds, nc=args.colors,
                         excluded_planes=excl, min_spots=args.min_spots,
                         bins=args.bins, bw=args.bin_width,
                         search_px=args.search,
                         subtract_mean_histogram=args.subtract_mean_histogram)
    print(f'{args.npz}: {len(spots)} spots, {len(beads)} beads; '
          f'{args.rounds} rounds x {args.colors} colors')
    print(f"planes kept: {out['n_planes_kept']} "
          f"(>= {args.min_spots} spots, excluding {out['excluded_planes']}); "
          f"outlier planes dropped: {out['n_outlier_planes']}")
    if out['rounds_missing']:
        print(f"INCOMPLETE TRACE: no estimate for round(s) "
              f"{out['rounds_missing']} — QC cannot pass")
    print(f"QC gate: {'PASS' if out['qc_pass'] else 'FAIL'} "
          f"({out['n_planes_usable']} usable planes, median pair weight "
          f"{out['median_pair_weight']:.3g}, "
          f"normalized {out['norm_weight']:.3g})")
    if out['n_pairs_peak_outside_window']:
        print(f"SHIFT SEARCH TOO NARROW: {out['n_pairs_peak_outside_window']} "
              f"pair(s) have their correlation peak beyond the "
              f"{out['search_px']:.0f} px search window, so their shifts "
              f"cannot be measured and QC FAILS. Re-run with "
              f"--search larger than the motion you expect (the largest "
              f"in-window estimate here is "
              f"{out['max_abs_pair_shift_px']:.1f} px). A capped search does "
              f"NOT announce itself: every pair is wrong the same way, so the "
              f"trace stays self-consistent and the residual budget stays "
              f"small.")
    print(f"pair/node consistency: max residual "
          f"{out['max_pair_residual_px']:.3f} px "
          f"(small residuals do NOT by themselves imply a correct trace)")
    print()
    print('round   signal trace (px)      bead trace (px)        diff (px)')
    for r in range(args.rounds):
        s = out['trace_signal_px'][r]; b = out['trace_bead_px'][r]; f = out['trace_diff_px'][r]
        print(f'{r:5d}   {s[0]:8.3f} {s[1]:8.3f}   {b[0]:8.3f} {b[1]:8.3f}   '
              f'{f[0]:8.3f} {f[1]:8.3f}')
    rms = out['trace_rms_vs_bead_px']
    print(f'\ntrace RMS vs bead oracle: '
          + (f'{rms:.3f} px' if np.isfinite(rms) else 'n/a (no beads)'))
    print('per-color rigid offsets rel. color 0 (px, jackknife SE over rounds):')
    for c in range(args.colors):
        o = out['color_offsets_px'][c]; se = out['color_offsets_se_px'][c]
        print(f'  color {c}: ({o[0]:.3f}, {o[1]:.3f}) +- ({se[0]:.3f}, {se[1]:.3f})')
    if args.json:
        _write_json(args.json, out)
    # Deliberate exit-status contract: 0 = ran and QC PASS, 2 = ran but QC
    # FAIL (trace printed for inspection, not trustworthy), 1 = input error.
    return 0 if out['qc_pass'] else 2


# --------------------------------------------------------- smlm-budget ----

EVENT_COLUMNS = """\
Expected event-table columns (CSV with header, or .npz):
  group   integer grouping unit (e.g. picked structure / origami / FOV cell);
          statistics are averaged over independent groups (not localizations)
  code    physical code of the event's source, integer 0..3
  panel   temporal panel, integer 0..3: panels 0,1 = window A (share one
          persistent displacement), panels 2,3 = window B
  x       event coordinate along ONE axis (units: 1/q; nm for default q)
  weight  optional positive event weight (photons/quality); default 1
.npz alternative: either an 'events' (n,4|5) array with the same column
order, or precomputed 'z' (groups,4,4,nq) complex and 'n' (groups,4,4)
aggregates plus optional 'q' (nq)."""


def _load_events(path):
    p = Path(path)
    if p.suffix == '.npz':
        d = np.load(p)
        if 'z' in d and 'n' in d:
            return ('zn', d['z'], d['n'], d['q'] if 'q' in d else None)
        if 'events' in d:
            ev = np.asarray(d['events'], float)
            if ev.ndim != 2 or ev.shape[1] not in (4, 5):
                raise ValueError("'events' must be (n, 4) or (n, 5): "
                                 "group, code, panel, x[, weight]")
            w = ev[:, 4] if ev.shape[1] == 5 else None
            return ('events', (ev[:, 0], ev[:, 1], ev[:, 2], ev[:, 3], w), None, None)
        raise ValueError(f"{path}: expected 'events' or 'z'+'n' arrays; "
                         f"found {list(d.keys())}")
    a = np.atleast_1d(np.genfromtxt(p, delimiter=',', names=True))
    names = a.dtype.names or ()
    for col in ('group', 'code', 'panel', 'x'):
        if col not in names:
            raise ValueError(f'{path}: missing column {col!r} (have {names})\n'
                             + EVENT_COLUMNS)
    w = a['weight'] if 'weight' in names else None
    return ('events', (a['group'], a['code'], a['panel'], a['x'], w), None, None)


def _format_confidence_set(intervals):
    if intervals is None:
        return 'unavailable'
    if not intervals:
        return 'empty'
    pieces = []
    for interval in intervals:
        lo, hi = interval['lower'], interval['upper']
        pieces.append(('[' if interval['lower_closed'] else '(')
                      + ('-inf' if lo is None else f'{lo:.6g}') + ', '
                      + ('+inf' if hi is None else f'{hi:.6g}')
                      + (']' if interval['upper_closed'] else ')'))
    return ' union '.join(pieces)


def _temporal_report(args):
    from .reporting import budget as reporting_budget
    kind, ev, n, qfile = _load_events(args.table)
    q = np.asarray(args.q, float) if args.q else None
    if kind == 'zn':
        z = ev
        if q is None and qfile is not None:
            q = np.asarray(qfile)
        if q is None:
            q = smlm.Q
        nev = None
    else:
        group, code, panel, x, w = ev
        if q is None:
            q = smlm.Q
        z, n, _ = smlm.aggregate_events(group, code, panel, x, w, q)
        nev = len(x)
    return reporting_budget(z, n, q, confidence_level=args.confidence_level), nev


def cmd_smlm_budget(args):
    try:
        out, nev = _temporal_report(args)
    except (ValueError, OSError) as error:
        print(f'smlm-budget: {error}', file=sys.stderr)
        return 2
    print(f"{args.table}: {out['n_groups']} groups"
          + (f', {nev} events' if nev is not None else '')
          + f"; q = {' '.join('%g' % v for v in out['q'])} (1/coordinate unit)")
    variants = ('arithmetic', 'geometric') if args.uncertainty_variant == 'both' else (args.uncertainty_variant,)
    reports = {'arithmetic': out['uncertainty'], 'geometric': out['geometric_uncertainty']}
    for variant in variants:
        uncertainty = reports[variant]
        if variant == 'arithmetic':
            print(f'\nPointwise approximate {100*args.confidence_level:g}% Fieller sets '
                  'for arithmetic matched-count attenuation:')
        else:
            print(f'\nPointwise conservative {100*args.confidence_level:g}% geometric Fieller product sets '
                  f"(two {100*uncertainty['window_confidence_level']:g}% window sets):")
            print('  Joint support uses Bonferroni; windows need not be independent. '
                  'This target can differ from arithmetic under heterogeneous motion.')
        for row in uncertainty['per_frequency']:
            status = 'finite positive set' if row['usable'] else f"UNRESOLVED ({row['reason']})"
            print(f"  q={row['q']:g}: {row['kind']}; {status}")
            ratio = 'undefined' if row['ratio'] is None else f"{row['ratio']:.6g}"
            print(f"    ratio={ratio}; set={_format_confidence_set(row['intervals'])}")
            print('    Gaussian-equivalent T set (conditional, coordinate units^2): '
                  + _format_confidence_set(row['gaussian_equivalent_T_intervals']))
            for name, window in zip('AB', row.get('windows', [])):
                print(f"    window {name}: {window['kind']}; "
                      + ('usable' if window['usable'] else f"UNRESOLVED ({window['reason']})")
                      + f"; full ratio set={_format_confidence_set(window['intervals'])}")
            if row['usable']:
                print('    ratio=1 excluded: ' + ('yes' if row['unit_attenuation_excluded'] else 'no'))
            if row['attenuation_model_compatible'] is False:
                print('    MODEL TENSION: set has no support in 0 < attenuation <= 1; '
                      'negative window estimates are retained.')
            concentration = row.get('group_concentration', {}).get('reference', {})
            share = concentration.get('largest_absolute_fraction')
            if share is not None:
                print(f'    largest absolute reference-group contribution: {share:.1%} '
                      '(descriptive; not an effective independent sample size)')
    sp = out['split']
    if args.legacy_fits:
        print('\nLegacy bounded fits (historical reproduction; not confidence bounds):')
        print(f"  T (pooled cubic fit):        {out['T_pooled']:8.3f}"
              + ('   [BOUNDARY HIT]' if out['pooled']['boundary'] else ''))
        print(f"  T (window-split arithmetic): {sp['total_arithmetic']:8.3f}"
              + ('   [BOUNDARY HIT]' if sp['boundary']['arithmetic'] else ''))
        print(f"  T (geometric, drift-robust): {sp['total_geometric']:8.3f}"
              + ('   [BOUNDARY HIT]' if sp['boundary']['geometric'] else ''))
    print('\nWindow-split D diagnostic (necessary model check; a quiet D '
          'does not establish stationarity):')
    print('     q     ratio_A   ratio_B      D      z(logD)   drift est. (units^2)')
    for i, qq in enumerate(out['q']):
        z_i = sp['D_z'][i]
        print(f"  {qq:6.3f}  {sp['ratio_A'][i]:8.4f}  {sp['ratio_B'][i]:8.4f}  "
              f"{sp['D'][i]:7.4f}  {z_i if z_i is None else format(z_i, '7.2f')}   "
              f"{sp['drift_estimate_nm2'][i]:10.3f}")
    if sp['negative_ratio']:
        print('  WARNING: a window ratio is non-positive; log-based D entries '
              'are undefined (kept as null, not clipped)')
    print('\nPhysical assumptions are NOT verified. Groups must be independent; '
          'localized events within a group are not independent replicates.\n'
          'Coverage is approximate and pointwise, and may fail with concentrated '
          'contributions. Finite sets do not imply useful resolution.\n'
          'The Gaussian variance interpretation requires the stated motion model; '
          'heterogeneous source motion changes the weighted target.\n'
          'Full ratio sets are retained; T sets use only their positive-ratio part. '
          'Historical fit keys remain in JSON; use --legacy-fits to print them.')
    print('Exit status checks requested uncertainty variants only: ' + ', '.join(variants)
          + '. No automatic variant selection or physical qualification is performed.')
    out['reported_uncertainty_variants'] = list(variants)
    out['reported_uncertainty_usable'] = all(reports[v]['all_frequencies_usable'] for v in variants)
    if args.json:
        _write_json(args.json, out)
    return 0 if out['reported_uncertainty_usable'] else 2


def cmd_smlm_phase_budget(args):
    """Explicit opt-in path; never substitutes for the historical estimators."""
    from .phase_attenuation import phase_budget
    try:
        kind, data, n, qfile = _load_events(args.table)
        q = np.asarray(args.q, float) if args.q else qfile
        q = smlm.Q if q is None else np.asarray(q)
        if kind == 'zn':
            z = data
        else:
            group, code, panel, x, weight = data
            z, n, _ = smlm.aggregate_events(group, code, panel, x, weight, q)
        out = phase_budget(z, n, q, args.confidence_level)
    except (ValueError, OSError) as error:
        print(f'smlm-phase-budget: {error}', file=sys.stderr)
        return 2
    report = out['phase_uncertainty']
    print(f"{args.table}: {out['n_groups']} groups; opt-in phase-preserving attenuation magnitude")
    print(f"Pointwise {100*args.confidence_level:g}% Fieller Bonferroni sets "
          f"(four {100*report['component_confidence_level']:g}% component sets).")
    for row in report['per_frequency']:
        status = 'finite positive magnitude support' if row['usable'] else f"UNRESOLVED ({row['reason']})"
        print(f"  q={row['q']:g}: {row['kind']}; {status}")
        ratio = 'undefined' if row['ratio'] is None else f"{row['ratio']:.6g}"
        print(f"    magnitude={ratio}; set={_format_confidence_set(row['intervals'])}")
        print('    Gaussian-equivalent half-increment-variance set (conditional, units^2): '
              + _format_confidence_set(row['gaussian_equivalent_T_intervals']))
        if row['gaussian_transform_reason']:
            print('    Gaussian transform unavailable: '+row['gaussian_transform_reason'])
        for name, window in zip('AB',row['windows']):
            for axis in ('real','imag'):
                part = window[axis]
                print(f"    window {name} {axis}: {part['kind']}; "
                      + _format_confidence_set(part['intervals']))
        if row['attenuation_model_compatible'] is False:
            print('    MODEL TENSION: no support in 0 < magnitude <= 1; estimates are not clipped.')
    print('Physical assumptions are NOT verified. This opt-in extension cancels a nonzero '
          'source-uniform complex kernel change under the declared common-motion model.\n'
          'Negative real or imaginary components are legitimate. Mean displacement and '
          'kernel phase are confounded; neither is recovered. Source-dependent changes, '
          'code dependence and shared errors remain unresolved.\n'
          'Groups must be independent; within-group events are not replicates. Internal '
          'groups cannot estimate across-acquisition shared-motion variation.\n'
          'Coverage is approximate and pointwise; the four-component rectangle can lose '
          'resolution. The sample magnitude is biased. No automatic variant selection '
          'or independent physical validation is supplied.\n'
          'Exit status checks finite positive attenuation support, not physical accuracy '
          'or a finite Gaussian variance interval. Existing smlm-budget results are unchanged.')
    if args.json:
        _write_json(args.json,out)
    return 0 if report['all_frequencies_usable'] else 2


# --------------------------------------------------------------- parser ---


# ------------------------------------------------------- magnification ----

def cmd_magnification(args):
    d = np.load(args.npz)
    spots = np.asarray(d['spots'], float)
    colours = [c.strip() for c in args.colours.split(',')] if args.colours else \
        [str(int(c)) for c in sorted(set(spots[:, 1].astype(int)))]
    planes, k = {}, 0
    for r in sorted(set(spots[:, 0].astype(int))):
        for ci in sorted(set(spots[spots[:, 0] == r, 1].astype(int))):
            planes[k] = (r, ci); k += 1
    if args.round1:
        first = min(r for r, _ in planes.values()) if args.first_round is None else args.first_round
        est = chromatic.round1_qc(spots, planes, colours, first_round=first,
                                  field_px=args.field_px)
        ok, why = est['gate_pass'], est['gate_reason']
    else:
        est = chromatic.estimate_scales(spots, planes, colours, field_px=args.field_px,
                                        within_round_only=True)
        ok, why = chromatic.gate(est)
    out = dict(mode='round1' if args.round1 else 'all rounds', gate_pass=bool(ok), gate_reason=why,
               scales=est.get('scales'), scale_se=est.get('scale_se'),
               scale_se_gauge=est.get('scale_se_gauge', colours[0]),
               uncertainty_status='uncalibrated residual-fit diagnostic; not a confidence interval',
               gate_interpretation='frozen empirical eligibility rule; not an accuracy guarantee',
               same_colour_control=est.get('same_colour_control_median'),
               residual_rms_permille=est.get('residual_rms_permille'), ntile=est.get('ntile'))
    if est.get('scales'):
        ref = colours[-1]
        rel = chromatic.relative_to(est['scales'], ref)
        out['scales_relative_to_' + ref] = rel
        print(f"per-colour scale relative to {ref}: " +
              ', '.join(f'{c} {rel[c]:.5f} ({(rel[c] - 1) * 1e3:+.2f} per mille)' for c in colours))
    print(f"gate: {'PASS' if ok else 'FAIL'} ({why})")
    print('Accuracy interval: unavailable. The gate is an empirical eligibility rule; '
          'the fitted SE is not calibrated uncertainty.')
    if args.json:
        _write_json(args.json, out)
    return 0 if ok else 2


def cmd_certificate(args):
    names, W = certificate.load_codebook_csv(args.csv, include_blanks=args.include_blanks)
    rec = certificate.certify(Path(args.csv).name, W,
                              random_per_source=args.random_per_source)
    rec['blank_policy'] = 'included' if args.include_blanks else 'excluded by name'
    print(f"{rec['label']}: {rec['n_targets']} targets x {rec['n_planes']} planes -> "
          f"class {rec['verdict']}: {rec.get('reason', '')}")
    if 'informative_invariants' in rec:
        print(f"  covered plane pairs {rec.get('covered_pairs')}, rank {rec.get('rank_N')}, "
              f"nuisance-free dimension {rec['informative_invariants']}")
    if args.json:
        _write_json(args.json, rec)
    return 0


def cmd_retrofit(args):
    """What one extra per-molecule randomized readout would buy this codebook."""
    names, W = certificate.load_codebook_csv(args.csv, include_blanks=args.include_blanks)
    rep = retrofit.report(W, n_spots=tuple(args.spots), seed=args.seed)
    rep['blank_policy'] = 'included' if args.include_blanks else 'excluded by name'
    idm = rep['identified']
    print(f"{args.csv}: {rep['n_genes']} genes x {rep['n_bits']} bits")
    print(f"  mixture cancellation (max |residual|): "
          f"{rep['mixture_cancellation_max_abs']:.2e}  (exact by construction)")
    print(f"  budget map rank {idm['rank']} of {idm['n_budgets']}; "
          f"all budget-plus-blur contrasts identified: {idm['all_contrasts_identified']}")
    if idm['all_contrasts_identified']:
        print('  only the global budget-plus-blur level is unrecoverable')
    else:
        print('  additional contrast directions are unidentifiable; do not invert as full budgets')
    print('  displacement/blur separation under independent plane blur: '
          + str(rep['blur_separation']['all_targets_identified']))
    for k, v in rep['noise'].items():
        print(f"  {k:<9} contrast noise sd {v['sd']:8.2f} nm^2  "
              f"(falls as 1/sqrt(N); not a floor)")
    print(f"  estimand: {rep['estimand']}")
    print(f"  {rep['status']}")
    if args.json:
        _write_json(args.json, {k: v for k, v in rep.items()
                                if k != 'identified'} |
                    {'identified': {kk: vv for kk, vv in idm.items()
                                    if kk != 'null_vector'}})
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog='codedcal',
        description='Calibration error budgets for coded localization '
                    'measurements (count-balanced SMLM contrasts; '
                    'fiducial-free coded-FISH registration).')
    ap.add_argument('--version', action='version', version=f'codedcal {__version__}')
    sub = ap.add_subparsers(dest='cmd', required=True)

    d = sub.add_parser('demo-chen2015',
                       help='download the public Chen 2015 MERFISH example '
                            f'({CHEN_SIZE_MB} MB) and reproduce the '
                            'fiducial-free trace vs bead oracle (~0.17 px RMS)',
                       description=f'Downloads {CHEN_URL} (~{CHEN_SIZE_MB} MB, '
                                   'cached), runs the fiducial-free trace on the '
                                   '16-round example FOV and compares it with the '
                                   'fiducial-bead oracle. Expected: ~0.17 px RMS '
                                   'agreement (pre-registered criterion 0.2 px).')
    d.add_argument('--dir', default=str(Path.home() / '.cache' / 'codedcal' / 'chen2015'),
                   help='download/cache directory (default: %(default)s)')
    d.add_argument('--json', help='also write the full result dict to this JSON file')
    d.set_defaults(func=cmd_demo)

    f = sub.add_parser('fish-trace',
                       help='fiducial-free trace + color offsets + QC on a '
                            'reduced spot/bead npz',
                       description="Input npz: 'spots' (m,>=4) rows "
                                   "[round, color_index, x, y, ...] and optional "
                                   "'beads' (n,>=3) rows [round, x, y, ...] "
                                   '(fiducials, used only as the comparison '
                                   'oracle). Coordinates in px on a '
                                   '[0, bins*bin_width)^2 field. Defaults match '
                                   'the TMS/CZ-Biohub MERFISH liver extraction '
                                   '(8 rounds x 3 colors, 2048 px field).')
    f.add_argument('npz')
    f.add_argument('--rounds', type=int, default=8, help='imaging rounds (default 8)')
    f.add_argument('--colors', type=int, default=3, help='color channels (default 3)')
    f.add_argument('--exclude-plane', action='append', metavar='R,C',
                   default=None,
                   help='exclude plane round,color from the trace (repeatable; '
                        'default: 6,1 and 7,0 — the VZG116 codebook-unused bits). '
                        "pass --exclude-plane '' to exclude nothing")
    f.add_argument('--min-spots', type=int, default=200,
                   help='minimum spots for a plane to enter the trace graph '
                        '(default 200)')
    f.add_argument('--bins', type=int, default=1024, help='histogram bins (default 1024)')
    f.add_argument('--bin-width', type=float, default=2.0,
                   help='histogram bin width, px (default 2.0; field = bins*width)')
    f.add_argument('--search', type=float, default=30.0, metavar='PX',
                   help='positive coarse shift-search radius in pixels (default 30); '
                        'use e.g. 150 for larger inter-round motion')
    f.add_argument('--subtract-mean-histogram', action='store_true',
                   help='reproduce the historical frozen variant (per-color '
                        'mean-histogram subtraction). OFF by default: the '
                        'sensitivity audit and blind clean-room showed it '
                        'manufactures ghost correlations; NOSUB is the '
                        'validated primary')
    f.add_argument('--json', help='also write the full result dict to this JSON file')
    f.set_defaults(func=cmd_fish_trace)

    s = sub.add_parser('smlm-budget',
                       help='attenuation confidence sets + window-split D diagnostic from '
                            'an SMLM event table (CSV or npz)',
                       description=EVENT_COLUMNS,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    s.add_argument('table', help='CSV (header: group,code,panel,x[,weight]) or .npz')
    s.add_argument('--q', type=float, nargs='+', default=None,
                   help='spatial frequencies, 1/coordinate-unit '
                        '(default 0.10 0.15, i.e. nm^-1 for nm coordinates)')
    s.add_argument('--json', help='also write the full result dict to this JSON file')
    s.add_argument('--confidence-level', type=float, default=.95,
                   help='pointwise approximate confidence level, strictly between 0 and 1 (default .95)')
    s.add_argument('--uncertainty-variant', choices=['arithmetic', 'geometric', 'both'], default='arithmetic',
                   help='reports determining printed confidence sets and exit status; JSON always includes both '
                        '(default arithmetic, preserving prior exit semantics; no automatic selection)')
    s.add_argument('--legacy-fits', action='store_true',
                   help='also print historical bounded T fits; these are not confidence bounds')
    s.set_defaults(func=cmd_smlm_budget)

    p = sub.add_parser('smlm-phase-budget',
                       help='opt-in phase-preserving attenuation magnitude with full confidence sets',
                       description=EVENT_COLUMNS,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('table', help='CSV or complex aggregate NPZ; same layout as smlm-budget')
    p.add_argument('--q', type=float, nargs='+', default=None,
                   help='positive spatial frequencies in inverse coordinate units')
    p.add_argument('--confidence-level',type=float,default=.95,
                   help='pointwise approximate joint confidence level (default .95)')
    p.add_argument('--json',help='also write full component and magnitude confidence sets')
    p.set_defaults(func=cmd_smlm_phase_budget)

    m = sub.add_parser('magnification',
                       help='assignment-free per-colour magnification with its '
                            'quality gate (--round1: first round only, pre-decode)',
                       description="Input npz: 'spots' (m,>=4) rows [round, "
                                   "colour_index, x, y, ...]. Uses within-round "
                                   'cross-colour plane pairs only (immune to '
                                   'inter-round jumps). Scales are reported '
                                   'relative to the LAST colour. For --round1, '
                                   'detect spots on each image alone.')
    m.add_argument('npz')
    m.add_argument('--colours', help='comma-separated colour labels in colour-index '
                                     'order (e.g. 750,650,560)')
    m.add_argument('--round1', action='store_true',
                   help='estimate from the first round only (gate without the '
                        'same-colour control, which needs two rounds)')
    m.add_argument('--first-round', type=int, default=None)
    m.add_argument('--field-px', type=float, default=2048.)
    m.add_argument('--json', help='also write the result dict to this JSON file')
    m.set_defaults(func=cmd_magnification)

    r = sub.add_parser('retrofit',
                       help='what one randomized readout would buy a codebook '
                            '(mixture cancellation, identified contrasts, noise)',
                       description='Input: MERFISH-style codebook CSV. Reports '
                                   'the abundance-weighted contrast that cancels '
                                   'the per-gene mixture identically, which budget '
                                   'contrasts become identified, and the finite-N '
                                   'randomization noise. SIMULATION ONLY: no '
                                   'instrument has run a per-molecule randomized '
                                   'readout.')
    r.add_argument('csv')
    r.add_argument('--include-blanks', action='store_true',
                   help='include named blank/control codewords as targets; default excludes them')
    r.add_argument('--spots', type=int, nargs='+', default=[100, 1000],
                   help='spots per gene for the noise estimate (default 100 1000)')
    r.add_argument('--seed', type=int, default=0)
    r.add_argument('--json', help='write the full report as JSON')
    r.set_defaults(func=cmd_retrofit)

    c = sub.add_parser('certificate',
                       help='design-time certificate of a codebook (which '
                            'cross-round error components its data can identify)',
                       description='Input: MERFISH-style codebook CSV (name, id, '
                                   'one 0/1 column per bit). B/C are legacy '
                                   'sensitivity classes; C is not always exact '
                                   'non-identifiability. Random-source claims '
                                   'are conditional on assumptions not checked here.')
    c.add_argument('csv')
    c.add_argument('--include-blanks', action='store_true',
                   help='include named blank/control codewords as targets; default excludes them')
    c.add_argument('--random-per-source', action='store_true',
                   help='declare physical per-source randomization (conditional; not verified)')
    c.add_argument('--json', help='also write the certificate record to this JSON file')
    c.set_defaults(func=cmd_certificate)

    args = ap.parse_args(argv)
    if hasattr(args, 'confidence_level') and (not np.isfinite(args.confidence_level)
                                              or not 0 < args.confidence_level < 1):
        ap.error('--confidence-level must be finite and between 0 and 1')
    if getattr(args, 'cmd', None) == 'smlm-budget' and args.q is not None:
        if not np.isfinite(args.q).all() or np.any(np.asarray(args.q) <= 0):
            ap.error('--q must contain finite positive spatial frequencies')
    if hasattr(args, 'search') and (not np.isfinite(args.search) or args.search <= 0):
        ap.error('--search must be finite and positive')
    if getattr(args, 'exclude_plane', 'missing') is None:
        args.exclude_plane = ['6,1', '7,0']
    elif getattr(args, 'exclude_plane', None) == ['']:
        args.exclude_plane = []
    return args.func(args)


if __name__ == '__main__':
    raise SystemExit(main())
