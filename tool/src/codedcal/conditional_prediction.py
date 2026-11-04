"""Calibration-region conditional image forecasts, without molecular claims.

Models are fit on a separate spatial calibration region. A forecast uses only
other readouts at the scoring pixels. Ridge accounts for correlated image
structure; it does not identify an RNA or a temporal calibration budget.
"""
import numpy as np


def fit_readout_forecast(calibration, *, ridge=0.1, source_predictions=None):
    """Fit every singleton conditional mean using separate calibration pixels.

    The penalty is ridge times the number of pixels after centering/scaling
    each feature. An optional fitting-only source forecast adds one feature per
    held bit. The source forecast must itself exclude that bit's target cells.
    """
    data = np.asarray(calibration, float)
    source = None if source_predictions is None else np.asarray(source_predictions, float)
    if (data.ndim != 2 or data.shape[0] < 2 or data.shape[1] < 2
            or not np.isfinite(data).all() or not np.isfinite(ridge) or ridge <= 0
            or source is not None and (source.shape != data.shape or not np.isfinite(source).all())):
        raise ValueError('finite pixels-by-readouts calibration and positive ridge required')
    models = []
    for held in range(data.shape[1]):
        fitting = np.delete(np.arange(data.shape[1]), held)
        X = data[:, fitting]
        if source is not None:
            X = np.column_stack([X, source[:, held]])
        center = X.mean(axis=0); scale = X.std(axis=0)
        usable = scale > 1e-12
        scale = np.where(usable, scale, 1.)
        Z = (X-center)/scale
        yc = float(data[:, held].mean())
        beta = np.linalg.solve(Z.T@Z+ridge*len(Z)*np.eye(Z.shape[1]), Z.T@(data[:, held]-yc))
        beta[~usable] = 0.
        models.append(dict(heldout=held, fitting_readouts=fitting.tolist(),
                           feature_center=center.tolist(), feature_scale=scale.tolist(),
                           coefficients=beta.tolist(), target_center=yc,
                           source_feature=source is not None))
    return dict(readouts=data.shape[1], pixels=len(data), ridge=float(ridge), models=models,
                interpretation='Conditional image forecast; not source identity, optical truth or a workflow recommendation')


def predict_readout(model, fitting_values, source_prediction=None):
    """Forecast one readout. API never accepts its held target values."""
    X = np.asarray(fitting_values, float)
    if X.ndim != 2 or X.shape[1] != len(model['fitting_readouts']) or not np.isfinite(X).all():
        raise ValueError('finite fitting-only pixel features required')
    if model['source_feature']:
        source = np.asarray(source_prediction, float)
        if source.shape != (len(X),) or not np.isfinite(source).all():
            raise ValueError('matching fitting-only source forecast required')
        X = np.column_stack([X, source])
    elif source_prediction is not None:
        raise ValueError('unexpected source feature')
    center = np.asarray(model['feature_center']); scale = np.asarray(model['feature_scale'])
    beta = np.asarray(model['coefficients'])
    if center.shape != (X.shape[1],) or scale.shape != center.shape or beta.shape != center.shape:
        raise ValueError('model feature dimensions disagree')
    return (X-center)/scale@beta+model['target_center']


def robust_readout_covariance(calibration, clipping=3., shrinkage=.1):
    """Winsorized empirical covariance for a separate calibration region.

    Signals may remain in this covariance; it is not certified camera noise.
    Diagonal shrinkage makes the estimated spectral metric positive definite.
    """
    x = np.asarray(calibration, float)
    if (x.ndim != 2 or len(x) < 2 or x.shape[1] < 2 or not np.isfinite(x).all()
            or not np.isfinite(clipping) or clipping <= 0
            or not np.isfinite(shrinkage) or not 0 < shrinkage <= 1):
        raise ValueError('finite calibration, clipping and positive shrinkage required')
    center = np.median(x, axis=0)
    mad = 1.4826*np.median(abs(x-center), axis=0)
    if np.any(mad <= 0):
        raise ValueError('every readout needs nonzero calibration dispersion')
    z = np.clip(x-center, -clipping*mad, clipping*mad)
    z -= z.mean(axis=0)
    C = z.T@z/len(z)
    C = (1-shrinkage)*C+shrinkage*np.diag(np.diag(C))
    np.linalg.cholesky(C)
    return C
