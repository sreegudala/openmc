"""Shared helpers for the frequency-domain neutron-noise reference solvers.

Conventions (see OpenMC_neutron_noise_transport_design.md, Section 1.1):

* Fourier transform g(w) = int g(t) exp(-i w t) dt, so d/dt -> +i w and the
  operator gains +i w / v.
* All fission production is divided by the static k_eff of the *same* operator.
* Delayed production of family k carries lambda_k / (lambda_k + i w).
* Group-to-group scattering matrices are stored OpenMC-style as
  ``scatter[g_in, g_out]``; the emission into group g is
  ``sum_g' scatter[g', g] * phi[g']``.

All arrays are numpy ``complex128`` where a frequency is involved.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Material:
    """Multigroup material with prompt/delayed fission data.

    Parameters
    ----------
    total : (G,) total cross section [1/cm]
    scatter : (G, G) scatter[g_in, g_out] [1/cm]
    inv_v : (G,) inverse speed [s/cm]
    nu_fission_p : (G,) prompt nu-fission per incident group [1/cm]
    nu_fission_d : (K, G) delayed nu-fission per family and incident group
    chi_p : (G,) prompt spectrum
    chi_d : (K, G) delayed spectra
    lam : (K,) decay constants [1/s]
    """
    total: np.ndarray
    scatter: np.ndarray
    inv_v: np.ndarray
    nu_fission_p: np.ndarray = None
    nu_fission_d: np.ndarray = None
    chi_p: np.ndarray = None
    chi_d: np.ndarray = None
    lam: np.ndarray = field(default_factory=lambda: np.zeros(0))

    def __post_init__(self):
        self.total = np.atleast_1d(np.asarray(self.total, float))
        G = self.total.size
        self.scatter = np.asarray(self.scatter, float).reshape(G, G)
        self.inv_v = np.atleast_1d(np.asarray(self.inv_v, float))
        self.lam = np.atleast_1d(np.asarray(self.lam, float))
        K = self.lam.size
        if self.nu_fission_p is None:
            self.nu_fission_p = np.zeros(G)
        self.nu_fission_p = np.atleast_1d(np.asarray(self.nu_fission_p, float))
        if self.nu_fission_d is None:
            self.nu_fission_d = np.zeros((K, G))
        self.nu_fission_d = np.asarray(self.nu_fission_d, float).reshape(K, G)
        if self.chi_p is None:
            self.chi_p = np.eye(G)[0]
        self.chi_p = np.atleast_1d(np.asarray(self.chi_p, float))
        if self.chi_d is None:
            self.chi_d = np.tile(self.chi_p, (K, 1))
        self.chi_d = np.asarray(self.chi_d, float).reshape(K, G)

    @property
    def G(self):
        return self.total.size

    @property
    def nu_fission(self):
        return self.nu_fission_p + self.nu_fission_d.sum(axis=0)

    @property
    def absorption(self):
        return self.total - self.scatter.sum(axis=1)

    @classmethod
    def from_beta(cls, total, scatter, inv_v, nu_fission=None, chi=None,
                  beta=(), lam=(), chi_d=None):
        """Build from total nu-fission and family fractions beta_k
        (independent of incident group)."""
        total = np.atleast_1d(np.asarray(total, float))
        G = total.size
        beta = np.atleast_1d(np.asarray(beta, float))
        nuf = np.zeros(G) if nu_fission is None else np.atleast_1d(
            np.asarray(nu_fission, float))
        nfp = nuf * (1.0 - beta.sum())
        nfd = np.outer(beta, nuf)
        if chi is None:
            chi = np.eye(G)[0]
        return cls(total, scatter, inv_v, nfp, nfd, chi, chi_d, lam)

    def scaled(self, factor):
        """Return a copy with every reaction scaled by `factor` (density)."""
        return Material(self.total * factor, self.scatter * factor,
                        self.inv_v, self.nu_fission_p * factor,
                        self.nu_fission_d * factor, self.chi_p, self.chi_d,
                        self.lam)


def delayed_factor(lam, omega):
    """lambda_k / (lambda_k + i omega), shape (K,)."""
    lam = np.atleast_1d(np.asarray(lam, float))
    return lam / (lam + 1j * omega)


def chi_omega(mat, omega):
    """Frequency-dependent production matrix P[g_out, g_in] such that the
    fission emission is P @ phi (before division by k)."""
    P = np.outer(mat.chi_p, mat.nu_fission_p).astype(complex)
    f = delayed_factor(mat.lam, omega)
    for k in range(mat.lam.size):
        P += f[k] * np.outer(mat.chi_d[k], mat.nu_fission_d[k])
    return P


def scatter_out_in(mat):
    """Scatter matrix arranged as [g_out, g_in]."""
    return mat.scatter.T


def stats_report(name, value, reference, rtol):
    rel = abs(value - reference) / abs(reference)
    ok = rel <= rtol
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: value={value:.6e} "
          f"ref={reference:.6e} rel={rel:.2e} (tol {rtol:.0e})")
    return ok
