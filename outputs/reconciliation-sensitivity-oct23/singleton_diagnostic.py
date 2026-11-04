"""Post-hoc diagnostic: singleton-group code/diagonal agreement is algebraic."""
from sensitivity import *
from report import supported_fit_comparison


def run():
    results = []
    for row in registry():
        tag = row['tag']; d = np.load(P/f'groups-{tag}.npz')
        primary = json.loads((P/f'result-{tag}.json').read_text()); arms = []
        for i, arm in enumerate(primary['arms']):
            count = d['quadratic'][i, :, 0, 1]
            single = count == 1; multi = count > 1
            selected = d['selected_code'][i]; diagonal_comp = d['diagonal'][i]
            error = np.max(abs(selected[single]-diagonal_comp[single])) if single.any() else 0.
            assert np.allclose(selected[single], diagonal_comp[single], rtol=1e-10, atol=1e-7), (tag, arm['name'], error)
            mass = np.abs(diagonal_comp[..., 1]).sum((0, 2))
            single_mass = np.abs(diagonal_comp[single, ..., 1]).sum((0, 2))
            fraction = np.divide(single_mass, mass, out=np.full_like(mass, np.nan), where=mass > 0)
            compared = dict(selected_code=describe(selected*multi[:, None, None, None], d['tile']),
                diagonal=describe(diagonal_comp*multi[:, None, None, None], d['tile']))
            full = describe(d['full_code']*multi[:, None, None, None], d['tile'])
            supported = bool(multi.any() and supported_fit_comparison(compared, full))
            if supported:
                df = compared['diagonal']['fit']['total_arithmetic']
                old_gap = abs(full['fit']['total_arithmetic']-df)
                new_gap = abs(compared['selected_code']['fit']['total_arithmetic']-df)
            else:
                old_gap = new_gap = None
            arms.append(dict(name=arm['name'], singleton_groups=int(single.sum()),
                multi_source_groups=int(multi.sum()), singleton_identity_max_abs_error=float(error),
                singleton_reference_absolute_mass_fraction=fraction, supported=supported,
                full_to_diagonal_fit_gap=old_gap, selected_to_diagonal_fit_gap=new_gap,
                gap_narrows=new_gap < old_gap if supported else None,
                selected_code=compared['selected_code'], diagonal=compared['diagonal'], full_code=full))
        results.append(dict(tag=tag, arms=arms))
    count = dict(narrowed=sum(a['gap_narrows'] is True for r in results for a in r['arms']),
        worsened=sum(a['gap_narrows'] is False for r in results for a in r['arms']),
        unavailable=sum(a['gap_narrows'] is None for r in results for a in r['arms']))
    save(P/'singleton-diagnostic.json', dict(post_hoc=True, primary_grid_unchanged=True,
        scope='Remove deterministic singleton agreement; all full/selected/diagonal contrasts restricted to the same multi-source coarse groups. Inferred labels remain unvalidated.',
        totals=count, views=results))
    print(json.dumps(count, indent=2))


if __name__ == '__main__':
    run()
