"""Unsupported fits cannot count as successful reconciliation."""
from copy import deepcopy
from report import supported_fit_comparison


def test_unsupported_or_boundary_fit_is_not_a_narrowing_success():
    d = dict(fit=dict(boundary=dict(arithmetic=False)), point_T=[1., -2.])
    a = dict(selected_code=deepcopy(d), diagonal=deepcopy(d))
    assert supported_fit_comparison(a, d)
    a['selected_code']['point_T'][0] = None
    assert not supported_fit_comparison(a, d)
    a['selected_code'] = deepcopy(d)
    a['selected_code']['fit']['boundary']['arithmetic'] = True
    assert not supported_fit_comparison(a, d)
    a['selected_code']['fit'] = None
    assert not supported_fit_comparison(a, d)
