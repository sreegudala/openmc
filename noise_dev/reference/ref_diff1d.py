"""One-group bare-slab diffusion noise reference (design doc Section 3.3).

    -D dphi'' + (Sigma_a + i w/v - Upsilon(w) nuSigma_f / k) dphi = dS

on |x| < a (a = extrapolated half-width, dphi(+-a) = 0).  With
B^2(w) = (Upsilon nuSf / k - Sigma_a - i w/v) / D and the Green's function
G'' + B^2 G = delta(x - x'),

    G(x, x', w) = -sin B(a + x<) sin B(a - x>) / (B sin 2Ba),
    dphi(x)     = -(1/D) int G(x, x') dS(x') dx'.

k is the slab's own diffusion eigenvalue nuSf / (Sigma_a + D (pi/2a)^2)
unless given.  For an absorber source dS = -dSigma_a phi0 this is the
SSM-report form  nabla^2 dphi + B^2 dphi = dSigma_a phi0 / D.
"""

import argparse

import numpy as np

from ref_g0 import upsilon


def greens(x, xp, B, a):
    xl = np.minimum(x, xp)
    xg = np.maximum(x, xp)
    return -np.sin(B * (a + xl)) * np.sin(B * (a - xg)) / (B * np.sin(2 * B * a))


def b_squared(omega, D, sa, nuSf, inv_v, beta, lam, k):
    ups = upsilon(omega, beta, lam)
    return (ups * nuSf / k - sa - 1j * omega * inv_v) / D


def slab_noise(x, omega, a, D, sa, nuSf, inv_v, beta, lam, src_region,
               strength, k=None, source_fn=None, n_quad=400):
    """dphi at points x for a uniform source `strength` on src_region, or
    for a general source density source_fn(x') on src_region."""
    if k is None:
        k = nuSf / (sa + D * (np.pi / (2 * a))**2)
    B = np.sqrt(b_squared(omega, D, sa, nuSf, inv_v, beta, lam, k) + 0j)
    t, wq = np.polynomial.legendre.leggauss(n_quad)
    x1, x2 = src_region
    xp = 0.5 * (x2 - x1) * t + 0.5 * (x1 + x2)
    wq = 0.5 * (x2 - x1) * wq
    s = strength * np.ones_like(xp) if source_fn is None else source_fn(xp)
    x = np.atleast_1d(x)
    Gm = greens(x[:, None], xp[None, :], B, a)
    return -(Gm * (wq * s)[None, :]).sum(axis=1) / D


def self_test():
    """Check the closed form against a finite-difference solve."""
    ok = True
    print("ref_diff1d self-test")
    a, D, sa, iv = 50.0, 1 / 3, 0.01, 1 / 2.2e5
    nuSf = sa + D * (np.pi / (2 * a))**2
    for f in [0.1, 1, 100]:
        w = 2 * np.pi * f
        n = 4000
        h = 2 * a / n
        x = -a + h * (np.arange(n - 1) + 1)
        k = 1.0
        ups = upsilon(w, [0.0065], [0.0784])
        diag = 2 * D / h**2 + sa + 1j * w * iv - ups * nuSf / k
        A = (np.diag(np.full(n - 1, diag)) + np.diag(np.full(n - 2, -D / h**2), 1)
             + np.diag(np.full(n - 2, -D / h**2), -1))
        s = ((x >= 14) & (x <= 15)).astype(complex)
        s[np.isclose(x, 14) | np.isclose(x, 15)] = 0.5  # trapezoid ends
        fd = np.linalg.solve(A, s)
        cf = slab_noise(x, w, a, D, sa, nuSf, iv, [0.0065], [0.0784],
                        (14, 15), 1.0, k=1.0)
        mask = np.abs(x - 14.5) > 2
        err = (np.abs(fd - cf) / np.abs(cf))[mask].max()
        good = err < 2e-3
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] closed form vs FD f={f} Hz: "
              f"max rel err {err:.2e}")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--test", action="store_true")
    args = p.parse_args()
    if args.test:
        raise SystemExit(0 if self_test() else 1)


if __name__ == "__main__":
    main()
