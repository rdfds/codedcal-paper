"""Amplitude information conditional on the other selected source components.

The diagonal of a template's Gram alone ignores competing components. This
fixed-dictionary calculation includes them. Unknown source positions, missing
candidate words and empirical noise-model error remain outside its scope.
"""
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from .source_debias import spectral_gram


def conditional_amplitude_information(templates, codes, jitter_relative=1e-10):
    """Return active-dictionary profile curvature, with disclosed regularization.

    The tiny fixed ridge permits a numerical solve of singular dictionaries;
    redundant columns have vanishing information as this ridge tends to zero.
    Its curvature is an approximation, not a calibrated confidence interval.
    """
    gram = spectral_gram(templates,np.asarray(codes,float)).toarray()
    if (not np.isfinite(gram).all() or not np.isfinite(jitter_relative)
            or not 0 < jitter_relative <= 1e-6):
        raise ValueError('finite dictionary and small positive numerical ridge required')
    if not len(gram):
        return np.zeros(0),dict(columns=0,numerical_ridge=0.)
    diagonal = np.diag(gram).copy()
    ridge = jitter_relative*max(float(diagonal.max()),1e-12)
    L = cholesky(gram+ridge*np.eye(len(gram)),lower=True)
    inverse = solve_triangular(L,np.eye(len(gram)),lower=True)
    variance = np.sum(inverse**2,axis=0)
    information = np.minimum(1/variance,np.maximum(diagonal,0))
    return information,dict(columns=len(gram),numerical_ridge=ridge,
                            minimum_information_fraction=float(np.min(information/np.maximum(diagonal,1e-20))),
                            scope='Other selected components, fixed positions; not all unknown alternatives or physical confidence')
