"""codedcal: calibration error budgets for coded localization measurements.

Two estimator families, both assignment-free (no source identities, no
fiducial input to the estimate itself):

- ``codedcal.smlm``: count-balanced cubic temporal contrast for coded SMLM
  event tables, with approximate Fieller attenuation confidence sets and a
  conditional Gaussian variance interpretation, plus the window-split D
  necessary-model diagnostic.
- ``codedcal.fish``: fiducial-free cross-round registration trace for coded
  FISH (MERFISH-class) per-plane spot tables, with per-color offsets, a QC
  gate, and comparison against a fiducial-bead oracle.
- ``codedcal.chromatic``: assignment-free per-colour magnification from spot
  positions alone (pre-decode; round 1 suffices), its quality gate, and the
  inter-colour displacement field.
- ``codedcal.certificate``: design-time certificate of which cross-round
  error components a codebook's own data can identify.

Historical point-estimator numerics are copied from the research implementations.
Confidence-set handling extends Source Data 70/71. Arithmetic Fieller sets and
conservative geometric joint sets are documented in ``codedcal.uncertainty``
and ``codedcal.geometric_uncertainty``. They can target different functionals;
no automatic variant selector or physical-model verification is supplied.
The opt-in ``codedcal.phase_attenuation`` extension retains complex cross
moments and cancels a nonzero source-uniform complex kernel change under its
stated common-motion model; it does not change the historical estimators.
"""
__version__ = "0.5.1"

from . import certificate, chromatic, fish, io, smlm  # noqa: F401
