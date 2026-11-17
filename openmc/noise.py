"""Post-processing helpers for frequency-domain neutron noise tallies.

In noise mode (:attr:`openmc.Settings.noise`) particles carry complex weights
and tallies with a :class:`openmc.ComplexPartFilter` hold the real and
imaginary parts of each score as two bins. The statepoint additionally stores
the per-realization sum of Re x Im, which gives the covariance of the two
means and hence first-order uncertainties on amplitude and phase.
"""

import numpy as np

import openmc


def _complex_axis(tally):
    """Return (stride, n_filter_bins) for the ComplexPartFilter of a tally."""
    stride = 1
    found = None
    for f in reversed(tally.filters):
        if isinstance(f, openmc.ComplexPartFilter):
            found = stride
            break
        stride *= f.num_bins
    if found is None:
        raise ValueError(f'Tally {tally.id} has no ComplexPartFilter')
    n_bins = int(np.prod([f.num_bins for f in tally.filters]))
    return found, n_bins


def complex_tally_stats(tally):
    """Means, uncertainties and Re-Im covariance of a complex-part tally.

    Parameters
    ----------
    tally : openmc.Tally
        Tally read from a statepoint that has a :class:`openmc.ComplexPartFilter`

    Returns
    -------
    dict of numpy.ndarray
        Arrays of shape (n_filter_bins // 2, n_nuclides, n_scores), with the
        filter bins in tally order with the complex-part filter removed. Keys:
        ``re``, ``im``, ``re_std``, ``im_std`` (standard deviations of the
        means), ``cov`` (covariance of the two means), ``corr`` (correlation
        coefficient), ``value`` (complex mean), ``amplitude``,
        ``amplitude_std``, ``phase`` [rad] and ``phase_std`` [rad]. The
        amplitude and phase uncertainties are first-order propagations:

        .. math::

            \\sigma_A^2 = (x^2\\sigma_x^2 + y^2\\sigma_y^2 + 2xy\\,c)/A^2,
            \\quad
            \\sigma_\\varphi^2 = (y^2\\sigma_x^2 + x^2\\sigma_y^2 - 2xy\\,c)/A^4

    """
    stride, n_bins = _complex_axis(tally)
    idx = np.arange(n_bins)
    re_idx = idx[(idx // stride) % 2 == 0]
    im_idx = re_idx + stride

    n = tally.num_realizations
    s = np.asarray(tally.sum)
    s2 = np.asarray(tally.sum_sq)
    if tally.sum_cross is None:
        raise ValueError(f'Tally {tally.id} has no cross-product data')
    sc = np.asarray(tally.sum_cross)

    x, y = s[re_idx] / n, s[im_idx] / n
    if n > 1:
        vx = np.clip((s2[re_idx] / n - x**2) / (n - 1), 0.0, None)
        vy = np.clip((s2[im_idx] / n - y**2) / (n - 1), 0.0, None)
        cov = (sc[re_idx] / n - x * y) / (n - 1)
    else:
        vx = np.full_like(x, np.inf)
        vy = np.full_like(y, np.inf)
        cov = np.zeros_like(x)

    sx, sy = np.sqrt(vx), np.sqrt(vy)
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = cov / (sx * sy)
        A = np.hypot(x, y)
        varA = (x**2 * vx + y**2 * vy + 2 * x * y * cov) / A**2
        varp = (y**2 * vx + x**2 * vy - 2 * x * y * cov) / A**4
    return {
        're': x, 'im': y, 're_std': sx, 'im_std': sy, 'cov': cov,
        'corr': corr, 'value': x + 1j * y, 'amplitude': A,
        'amplitude_std': np.sqrt(np.clip(varA, 0.0, None)),
        'phase': np.arctan2(y, x),
        'phase_std': np.sqrt(np.clip(varp, 0.0, None)),
    }
