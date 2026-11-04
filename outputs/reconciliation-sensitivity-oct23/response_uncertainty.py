"""Propagate retained observation and conditional response-fit uncertainty."""
from sensitivity import *


def rectangle_image(observation, intercept, slope):
    if observation is None or slope[0] <= 0:
        return dict(kind='unsupported_or_unbounded', interval=None,
            reason='nonfinite_observation_set' if observation is None else 'slope_not_strictly_positive')
    corners = [(v-a)/b for v in observation for a in intercept for b in slope]
    return dict(kind='bounded', interval=[min(corners), max(corners)], reason=None)


def fit_response(rows, index):
    x = np.array([[1., r['T']] for r in rows]); y = np.array([r['estimate'][index] for r in rows])
    se = np.maximum([r['se'][index] for r in rows], .001)
    xx = x/se[:, None]; beta = np.linalg.lstsq(xx, y/se, rcond=None)[0]
    covariance_y = np.zeros((len(rows), len(rows)))
    for seed in sorted({r['seed'] for r in rows}):
        indices = [i for i, r in enumerate(rows) if r['seed'] == seed]
        influence = []
        for i in indices:
            group = np.load(OLD/rows[i]['group_records'])['by_group']
            influence.append((group[:, index]-y[i]*group[:, index+3])/group[:, index+3].mean())
        influence = np.stack(influence, 1)
        covariance_y[np.ix_(indices, indices)] = np.cov(influence.T, ddof=1)/len(influence)
    operator = np.linalg.solve(xx.T@xx, x.T/se**2)
    covariance = operator@covariance_y@operator.T
    critical = student.ppf(.99375, len(influence)-1)
    coefficient_intervals = beta[:, None]+np.array([-1, 1])*critical*np.sqrt(np.diag(covariance))[:, None]
    return beta, covariance, coefficient_intervals


def run():
    old = json.loads((OLD/'comparison.json').read_text())
    results = []
    for row in registry():
        tag = row['tag']; path = OLD/f'controls-{tag}.json'
        if not path.exists():
            results.append(dict(tag=tag, status='no_matched_preparation_response_controls',
                calibrated_interval=None)); continue
        controls = json.loads(path.read_text()); prior = next(a for a in old['acquisitions'] if a['tag'] == tag)
        observed = np.load(P/f'groups-{tag}.npz'); original = 4
        reports = {}
        for j, name in enumerate(WEIGHTS):
            b = blocks(observed['quadratic'][original, :, j], observed['tile'])
            observation = fieller_set(b[:, 0], b[:, 1], .975)
            interval = ([observation['intervals'][0]['lower'], observation['intervals'][0]['upper']]
                if observation['kind'] == 'bounded' else None)
            models = {}
            for model in ('gaussian', 'residual'):
                rows = [r for r in controls['records'] if r['model'] == model]
                beta, covariance, coefficient_intervals = fit_response(rows, j)
                old_fit = prior['response_sensitivity'][name]['models'][model]
                assert np.allclose(beta, [old_fit['intercept'], old_fit['slope']], rtol=0, atol=1e-10)
                image = rectangle_image(interval, coefficient_intervals[0], coefficient_intervals[1])
                models[model] = dict(intercept=beta[0], slope=beta[1], covariance=covariance,
                    coefficient_intervals=coefficient_intervals,
                    point=(observation['ratio']-beta[0])/beta[1], propagated=image,
                    observation_only_interval=([(v-beta[0])/beta[1] for v in interval]
                        if interval is not None and beta[1] > 0 else None))
            images = [v['propagated'] for v in models.values()]
            envelope = ([min(v['interval'][0] for v in images), max(v['interval'][1] for v in images)]
                if all(v['kind'] == 'bounded' for v in images) else None)
            reports[name] = dict(observation=observation, models=models,
                noise_model_sensitivity_envelope=envelope,
                original_response_qualified=prior['response_sensitivity'][name]['qualified_conditional_response'])
        results.append(dict(tag=tag, status='conditional_propagation_completed', weights=reports,
            original_controls_sha256=sha(path)))
    save(P/'response-uncertainty.json', dict(results=results,
        confidence='Approximate Bonferroni rectangle, nominal joint >=95%; observation 97.5%, each response coefficient 98.75%.',
        conditional_on='Inferred templates, original QC and each specified noise law; template/physical assignment uncertainty not covered.',
        noise_envelope_is_confidence_interval=False))


if __name__ == '__main__':
    run()
