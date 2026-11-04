"""Experimental compact-source prediction, not a native decoding pipeline.

Spatial proposals, candidate words, sparse selection and amplitude refitting
use fitting readouts only. Held measurements score predictions after fitting.
Multiple words may occupy one site. A fixed global signed PSF eliminates tile
boundaries; this does not establish that the PSF describes a public instrument.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter, label, center_of_mass
from scipy.sparse import csr_matrix, vstack
from scipy.optimize import nnls
from .source_debias import debias_amplitudes, spectral_gram
from .profile_information import conditional_amplitude_information


def compact_templates(shape, positions, raw_sigma=1., foreground_sigma=.6,
                      background_sigma=np.sqrt(4.36), radius=9):
    """Globally centered Gaussian-difference PSFs in sparse C-order rows."""
    h, w = shape; rows = []
    if (min(h, w) < 1 or np.asarray(positions).shape != (len(positions), 2)
            or not np.isfinite(positions).all() or raw_sigma <= 0
            or foreground_sigma < 0 or background_sigma <= foreground_sigma
            or radius < 1):
        raise ValueError('valid global PSF geometry required')
    small = np.arange(-radius, radius+1)
    yy, xx = np.meshgrid(small, small, indexing='ij')
    for x, y in positions:
        ix, iy = int(np.rint(x)), int(np.rint(y))
        raw = np.exp(-((xx+ix-x)**2+(yy+iy-y)**2)/(2*raw_sigma**2))
        raw /= raw.sum()
        psf = (gaussian_filter(raw, foreground_sigma, mode='constant')
               - gaussian_filter(raw, background_sigma, mode='constant'))
        xs, ys = xx+ix, yy+iy
        valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
        rows.append(csr_matrix((psf[valid], (np.zeros(valid.sum(), int),
                                            (ys[valid]*w+xs[valid]).astype(int))),
                               shape=(1, h*w)))
    return vstack(rows, format='csr') if rows else csr_matrix((0, h*w))


def select_amplitudes(templates, codes, target, penalties, iterations=2000,
                      tolerance=1e-6):
    """Nonnegative L1 selection in a sparse fixed-source dictionary."""
    gram = spectral_gram(templates, codes)
    b = np.sum((templates @ target)*codes, axis=1)
    L = float(np.max(np.asarray(abs(gram).sum(axis=1))))
    x = np.zeros(len(codes)); z = x.copy(); t = 1.; relative = None
    for it in range(1, iterations+1):
        previous = x
        x = np.maximum(z-(gram @ z-b+penalties)/L, 0.)
        old = t; t = .5*(1+np.sqrt(1+4*t*t))
        z = x+(old-1)/t*(x-previous)
        if it % 25 == 0 or it == iterations:
            fixed = np.maximum(x-(gram @ x-b+penalties)/L, 0.)
            relative = float(np.linalg.norm(x-fixed)/max(np.linalg.norm(x), 1e-12))
            if relative <= tolerance:
                break
    return x, dict(converged=bool(relative <= tolerance), iterations=it,
                   fixed_point_relative=relative)


def fit_compact_fold(stack, codebook, heldout, noise, *, gains=None, top_words=12,
                     peak_snr=4., selection_sigma=3., source_snr=4., border=9,
                     max_relative_error=.6, raw_sigma=1., nuisance_codes=None,
                     calibrated_readouts=None, candidate_rounds=3,
                     noise_covariance=None, noise_center=None, spatial_noise_factor=1.,
                     profile_sources=False, balanced_candidates=False):
    """Predict withheld images from jointly fitted compact multiword sources.

    Noise and gains must be fixed using a separate calibration region. An
    empirical noise floor is not a calibrated probability or discovery rate.
    Unidentifiable fitting words and unconverged fits cannot validate a source.
    Optional covariance uses generalized least squares and a fitting-only
    conditional nuisance forecast. It must be frozen on a separate region.
    It is an empirical model, not certified noise or molecular error control.
    """
    stack = np.asarray(stack, float); book = np.asarray(codebook, float)
    held = np.asarray(heldout); noise = np.asarray(noise, float)
    gains = np.ones(len(stack)) if gains is None else np.asarray(gains, float)
    nuisance = (np.empty((0, len(stack))) if nuisance_codes is None
                else np.asarray(nuisance_codes, float))
    calibrated = (np.ones(len(stack), bool) if calibrated_readouts is None
                  else np.asarray(calibrated_readouts, bool))
    if (stack.ndim != 3 or book.ndim != 2 or book.shape[1] != len(stack)
            or noise.shape != (len(stack),) or gains.shape != noise.shape
            or not np.isfinite(stack).all() or not np.isfinite(book).all()
            or not np.isfinite(noise).all() or not np.isfinite(gains).all()
            or np.any(book < 0) or np.any(noise <= 0) or np.any(gains < 0)
            or nuisance.ndim != 2 or nuisance.shape[1] != len(stack)
            or not np.isfinite(nuisance).all() or np.any(nuisance < 0)
            or calibrated.shape != (len(stack),)
            or held.ndim != 1 or not np.issubdtype(held.dtype, np.integer)
            or not len(held) or len(np.unique(held)) != len(held)
            or np.any(held < 0) or np.any(held >= len(stack)) or len(held) == len(stack)
            or not isinstance(top_words, int) or top_words < 1
            or not isinstance(candidate_rounds, int) or candidate_rounds < 1
            or peak_snr <= 0 or selection_sigma <= 0 or source_snr <= 0
            or max_relative_error < 0 or border < 1
            or not np.isfinite(spatial_noise_factor) or spatial_noise_factor < 1):
        raise ValueError('finite supported compact-fold inputs required')
    fitting = np.setdiff1d(np.arange(len(stack)), held)
    response = np.vstack([book, nuisance])*gains
    center = np.zeros(len(stack)) if noise_center is None else np.asarray(noise_center, float)
    if center.shape != (len(stack),) or not np.isfinite(center).all():
        raise ValueError('finite fixed calibration center required')
    conditioning = np.zeros((len(fitting),len(held)))
    if noise_covariance is None:
        whitening = np.diag(1/noise[fitting])
        codes = response[:,fitting]/noise[fitting]
    else:
        covariance = np.asarray(noise_covariance,float)
        if (covariance.shape != (len(stack),len(stack)) or not np.isfinite(covariance).all()
                or not np.allclose(covariance,covariance.T,rtol=1e-12,atol=1e-12)):
            raise ValueError('symmetric fixed calibration covariance required')
        try:
            np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as error:
            raise ValueError('positive definite calibration covariance required') from error
        C = covariance[np.ix_(fitting,fitting)]
        whitening = np.linalg.inv(np.linalg.cholesky(C).T)
        conditioning = np.linalg.solve(C,covariance[np.ix_(fitting,held)])
        codes = response[:,fitting]@whitening
    norms = np.linalg.norm(codes, axis=1); supported = norms > 0
    available = np.flatnonzero(supported[:len(book)])
    if not len(available):
        raise ValueError('no fitting code support')
    unit = codes[available]/norms[available, None]
    _, groups, counts = np.unique(np.round(unit, 12), axis=0, return_inverse=True, return_counts=True)
    identified = np.zeros(len(response), bool); identified[available] = counts[groups] == 1
    measured_fit = np.moveaxis(stack[fitting]-center[fitting,None,None],0,-1).reshape(-1,len(fitting))
    whitened = ((stack[fitting]-center[fitting,None,None])/noise[fitting,None,None]
                if noise_covariance is None else
                np.moveaxis((measured_fit@whitening).reshape(*stack.shape[1:],len(fitting)),-1,0))
    # Physical nonnegative bright-pixel proposals are invariant to rotations of
    # the GLS whitening basis; a positive-part norm in that basis is not.
    proposal = (stack[fitting]-center[fitting,None,None])/noise[fitting,None,None]
    image = np.linalg.norm(np.maximum(proposal, 0), axis=0)
    peaks = (image == maximum_filter(image, size=3)) & (image >= peak_snr*np.sqrt(len(fitting)))
    peaks[:border] = False; peaks[-border:] = False
    peaks[:, :border] = False; peaks[:, -border:] = False
    components, n = label(peaks)
    xy = np.asarray([(x, y) for y, x in center_of_mass(image, components, np.arange(1, n+1))]).reshape(-1, 2)
    templates = compact_templates(stack.shape[1:], xy, raw_sigma=raw_sigma)
    target = np.moveaxis(whitened, 0, -1).reshape(-1, len(fitting))
    baseline = center[held]+measured_fit@conditioning
    prediction = baseline.reshape(*stack.shape[1:],len(held))
    if not len(xy):
        return dict(source_rows=[], nuisance_rows=[], prediction=prediction, fitting_readouts=fitting.tolist(),
                    heldout_readouts=held.tolist(), spatial_sites=0,
                    candidate_columns=0,fitted_nuisance_sources=0,
                    selection=dict(converged=True), debias=dict(converged=True),
                    correlated_nuisance=noise_covariance is not None,
                    noise_center=center.tolist(),conditioning=conditioning.tolist(),
                    spatial_noise_factor=float(spatial_noise_factor),profile_information=[],
                    balanced_candidates=bool(balanced_candidates),
                    interpretation='Compact fitting-only experimental dictionary, fixed photometry; not truth/FDR/native JSIT')
    projection = templates @ target
    candidate_spaces = [(unit,codes,projection)]
    if balanced_candidates:
        # Likelihood-weighted cosine can hide a real word behind strong
        # responses in other channels. Also search the calibrated support
        # pattern with equal readout weights; joint GLS still uses the actual
        # responses/noise. This expands candidates, never certifies a label.
        channel = np.max(response[:,fitting],axis=0)
        pattern = (response[:,fitting] > 0).astype(float)
        equalized = (projection@np.linalg.inv(whitening))/np.maximum(channel,1e-12)
        pnorm = np.linalg.norm(pattern[available],axis=1)
        candidate_spaces.append((pattern[available]/pnorm[:,None],pattern,equalized))
    k = min(top_words, len(available))
    nuisance_ids = np.flatnonzero(supported[len(book):])+len(book)
    sites = []; words = []
    # One initial cosine shortlist can omit a weak co-located word behind
    # aliases of the dominant word. Expand with fitting-only orthogonal
    # residuals. This changes candidate discovery, not the held scoring rule.
    for site in range(len(projection)):
        discovered = set()
        for search_unit,search_codes,search_trace in candidate_spaces:
            trace=search_trace[site];residual_trace = trace.copy();pursuit = []
            for _ in range(candidate_rounds):
                order = np.argsort(-(search_unit @ residual_trace), kind='stable')[:k]
                discovered.update(available[order].tolist())
                best = int(available[order[0]])
                if best in pursuit:
                    break
                pursuit.append(best)
                coefficient, _ = nnls(search_codes[pursuit].T, trace)
                residual_trace = trace-coefficient @ search_codes[pursuit]
        candidates = sorted(discovered)+nuisance_ids.tolist()
        sites.extend([site]*len(candidates)); words.extend(candidates)
    sites = np.array(sites, int); words = np.array(words, int)
    T = templates[sites]; C = codes[words]
    energies = np.asarray(T.multiply(T).sum(axis=1)).ravel()
    penalties = selection_sigma*spatial_noise_factor*np.sqrt(energies)*norms[words]
    amplitudes, selection = select_amplitudes(T, C, target, penalties)
    snr = amplitudes*np.sqrt(energies)*norms[words]/spatial_noise_factor
    active = snr >= source_snr
    T = T[active]; C = C[active]; words = words[active]; sites = sites[active]
    amplitudes, debias = debias_amplitudes(T, C, target, amplitudes[active])
    profile_history = []; profile_snr = None
    if profile_sources:
        # Own-column SNR counts evidence shared with competing components as
        # though it belonged entirely to this source. Use active-dictionary
        # conditional information, then refit on fitting readouts only.
        for _ in range(3):
            information, info = conditional_amplitude_information(T,C)
            profile_snr = amplitudes*np.sqrt(information)/spatial_noise_factor
            keep = profile_snr >= source_snr
            profile_history.append(dict(**info,retained=int(keep.sum()),removed=int((~keep).sum())))
            if keep.all():
                break
            T,C,words,sites = T[keep],C[keep],words[keep],sites[keep]
            amplitudes,debias = debias_amplitudes(T,C,target,amplitudes[keep])
        information,info = conditional_amplitude_information(T,C)
        profile_snr = amplitudes*np.sqrt(information)/spatial_noise_factor
    # Positive dependent same-site codes or an exact nonnegative alternative
    # cannot support uniquely identified sources. This is a local exact check,
    # not a certificate against near aliases or other unknown positions.
    ambiguous = np.zeros(len(words), bool)
    alias_codes = codes if noise_covariance is None else response[:,fitting]
    for site in np.unique(sites):
        ids = np.flatnonzero((sites == site) & (amplitudes > 1e-12))
        active_codes = alias_codes[words[ids]]
        alias = np.linalg.matrix_rank(active_codes, tol=1e-8) < len(ids)
        mixture = amplitudes[ids] @ active_codes
        # Nonnegative alternatives cannot contain an on-bit where the target
        # is exactly zero. Exclude one active gene at a time: alternatives may
        # legitimately retain other active genes or replace a single gene
        # with a mixture. Excluding the entire active set misses those aliases.
        within_support = np.all(alias_codes[:, mixture <= 1e-12] <= 1e-12, axis=1)
        for word in words[ids]:
            if alias or word >= len(book):
                continue
            other = np.flatnonzero(supported & within_support & (np.arange(len(response)) != word))
            if len(other):
                try:
                    _, error = nnls(alias_codes[other].T, mixture, maxiter=max(100, 20*len(other)))
                    alias = error <= 1e-8*max(np.linalg.norm(mixture), 1e-12)
                except RuntimeError:
                    alias = True
        ambiguous[ids] = alias
    codes_held = response[words][:, held]
    conditional_codes = codes_held-response[words][:,fitting]@conditioning
    prediction = (baseline+np.asarray(T.T @ (amplitudes[:,None]*conditional_codes))).reshape(*stack.shape[1:],len(held))
    residual = np.moveaxis(stack[held], 0, -1).reshape(-1, len(held))-prediction.reshape(-1, len(held))
    projections = T @ residual; energies = np.asarray(T.multiply(T).sum(axis=1)).ravel()
    rows = []; nuisance_rows = []
    for i, word in enumerate(words):
        if word >= len(book):
            nuisance_rows.append(dict(column_index=i, nuisance_index=int(word-len(book)),
                                      position=xy[sites[i]].tolist(), peak=float(amplitudes[i]),
                                      prediction_trace=(amplitudes[i]*conditional_codes[i]).tolist(),
                                      expected_trace=(amplitudes[i]*codes_held[i]).tolist()))
            continue
        expected = amplitudes[i]*codes_held[i]
        if energies[i] <= 0:
            continue
        observed = projections[i]/energies[i]+expected
        on = np.linalg.norm(expected) > 1e-12
        error = float(np.linalg.norm(observed-expected)/np.linalg.norm(expected)) if on else None
        deletion = amplitudes[i]*conditional_codes[i]
        gain = float(energies[i]*(deletion @ deletion)+2*(projections[i] @ deletion))
        response_supported = bool(np.all(calibrated[book[word] > 0]))
        rows.append(dict(column_index=i, word=int(word), position=xy[sites[i]].tolist(), peak=float(amplitudes[i]),
                         fitting_identified=bool(identified[word]), expected_trace=expected.tolist(),
                         prediction_trace=(amplitudes[i]*conditional_codes[i]).tolist(),
                         profile_snr=float(profile_snr[i]) if profile_snr is not None else None,
                         profile_supported=bool(profile_snr is None or profile_snr[i] >= source_snr),
                         heldout_on_readouts=int((book[word, held] > 0).sum()),
                         ambiguous=bool(ambiguous[i]), response_supported=response_supported,
                         measured_trace=observed.tolist(), relative_error=error,
                         deletion_sse_gain=gain, accepted=bool(identified[word] and not ambiguous[i]
                         and (profile_snr is None or profile_snr[i] >= source_snr)
                         and response_supported and selection['converged']
                         and debias['converged'] and on and gain > 0 and error <= max_relative_error)))
    return dict(source_rows=rows, nuisance_rows=nuisance_rows, prediction=prediction, fitting_readouts=fitting.tolist(),
                heldout_readouts=held.tolist(), spatial_sites=len(xy), candidate_columns=len(sites),
                selection=selection, debias=debias,
                fitted_nuisance_sources=int((words >= len(book)).sum()),
                correlated_nuisance=noise_covariance is not None,
                noise_center=center.tolist(), conditioning=conditioning.tolist(),
                spatial_noise_factor=float(spatial_noise_factor),
                profile_information=profile_history,
                balanced_candidates=bool(balanced_candidates),
                interpretation='Compact fitting-only experimental dictionary, fixed photometry; not truth/FDR/native JSIT')
