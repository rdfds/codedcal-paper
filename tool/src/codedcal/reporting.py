"""Numerical admission checks around the unchanged historical temporal budget.

The 0.5.0 numerical modules remain available for exact reproduction. This
reporting entry point rejects nonrepresentable Gaussian transforms rather than
letting JSON nulls masquerade as statistical infinity or underflow as zero.
It changes neither the attenuation estimator nor its uncertainty construction.
"""
import numpy as np

from . import smlm


def budget(z, n, q=None, confidence_level=.95):
    """Return the historical budget only when its current report is representable.

    Coordinates and inverse-coordinate frequencies must use consistent units.
    Unsupported floating-point scales raise ValueError; rescaling coordinates
    and frequencies together can restore representability. A successful check
    still says nothing about physical assumptions or independent-group validity.
    """
    raw_q = smlm.Q if q is None else np.asarray(q)
    if np.iscomplexobj(raw_q) or np.iscomplexobj(n):
        raise ValueError('frequencies and counts must be real')
    frequencies = np.atleast_1d(np.asarray(raw_q, float))
    if (frequencies.ndim != 1 or not len(frequencies)
            or not np.isfinite(frequencies).all() or np.any(frequencies <= 0)):
        raise ValueError('q must contain finite positive spatial frequencies')
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        squared = frequencies*frequencies
    if not np.isfinite(squared).all() or np.any(squared < np.finfo(float).tiny):
        raise ValueError('frequency squared is outside normal floating-point support; '
                         'rescale coordinates and inverse-coordinate frequencies together')
    with np.errstate(over='ignore', under='ignore', invalid='ignore', divide='ignore'):
        result = smlm.budget(z, n, frequencies, confidence_level)
    for name in ('uncertainty', 'geometric_uncertainty'):
        for row in result[name]['per_frequency']:
            if row['ratio'] is not None and not np.isfinite(row['ratio']):
                raise ValueError('nonrepresentable attenuation point estimate; '
                                 'rescale aggregate weights consistently')
            for part in row['intervals'] or []:
                if any(part[key] is not None and not np.isfinite(part[key])
                       for key in ('lower', 'upper')):
                    raise ValueError('nonrepresentable attenuation confidence bound; '
                                     'rescale aggregate weights consistently')
            for part in row['gaussian_equivalent_T_intervals'] or []:
                if any(part[key] is not None and not np.isfinite(part[key])
                       for key in ('lower', 'upper')):
                    raise ValueError('nonrepresentable Gaussian confidence transform; '
                                     'rescale coordinates and frequencies together')
            # A finite positive non-unit endpoint cannot transform to exact zero.
            for part in row['intervals'] or []:
                for key in ('lower', 'upper'):
                    value = part[key]
                    if value is not None and value > 0 and value != 1:
                        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
                            transformed = -np.log(value)/row['q']**2
                        if transformed == 0 or not np.isfinite(transformed):
                            raise ValueError('nonrepresentable Gaussian confidence transform; '
                                             'rescale coordinates and frequencies together')
    result['numerical_reporting_guard'] = dict(
        passed=True, historical_estimator_changed=False,
        physical_assumptions_verified=False,
        scope='finite confidence endpoints and representable Gaussian transforms; '
              'not an exhaustive floating-point or input-validity certificate')
    return result
