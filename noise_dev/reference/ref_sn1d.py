"""1-D slab multigroup S_N reference for frequency-domain neutron noise.

* Gauss-Legendre quadrature, diamond difference, isotropic scattering.
* Vacuum or reflective boundaries on either side.
* Static eigenproblem, then the complex fixed-source noise problem
      [mu d/dx + Sigma_t + i w/v] dpsi = (S + P(w)/k) dphi / 2 + q
  with the fission source scaled by the solver's *own* k.
* Prescribed isotropic sources, arbitrary angular sources, or the
  perturbation source -dSigma_t psi0 + dS^T phi0/2 + P(w) dnuSf phi0/(2k)
  built from the solver's own static angular flux.
* Adjoint (bilinear transpose) mode: group matrices transposed, same +i w.

The linear system is solved densely.  Because streaming does not couple
groups, the transport sweep is built once per group as a response matrix
for unit isotropic cell sources and unit incoming boundary fluxes; the
boundary angular fluxes are explicit unknowns so reflective boundaries are
exact (no boundary iteration).

Angular index convention: m < M/2 are mu = +mu_p[m], m >= M/2 are
mu = -mu_p[m - M/2].  Scalar flux phi = sum_m w_m psi_m with sum w = 2.
"""

import argparse
import time

import numpy as np
import scipy.linalg as sla
from scipy.sparse.linalg import LinearOperator, eigs

from noise_common import Material, chi_omega, stats_report


class Slab:
    """Slab geometry made of regions.

    Parameters
    ----------
    regions : list of (width [cm], Material, n_cells)
    x0 : float, left edge coordinate
    bc : (left, right), each 'vacuum' or 'reflective'
    n_angles : even number of Gauss-Legendre points
    """

    def __init__(self, regions, x0=0.0, bc=("vacuum", "vacuum"), n_angles=16):
        widths, mats, cells = [], [], []
        self.materials = []
        for width, mat, n in regions:
            # materials are matched by identity (Material holds arrays, so
            # == is elementwise)
            idx = [i for i, m in enumerate(self.materials) if m is mat]
            if not idx:
                self.materials.append(mat)
                idx = [len(self.materials) - 1]
            widths += [width / n] * n
            cells += [idx[0]] * n
        self.h = np.array(widths)
        self.mat_index = np.array(cells)
        self.edges = x0 + np.concatenate([[0.0], np.cumsum(self.h)])
        self.x = 0.5 * (self.edges[1:] + self.edges[:-1])
        self.N = self.h.size
        self.G = self.materials[0].G
        self.bc = bc
        mu, w = np.polynomial.legendre.leggauss(n_angles)
        self.mu_p = mu[n_angles // 2:]
        self.w_p = w[n_angles // 2:]
        self.M2 = n_angles // 2
        self.mu = np.concatenate([self.mu_p, -self.mu_p])
        self.w = np.concatenate([self.w_p, self.w_p])

    # ----------------------------------------------------------- materials
    def cell_array(self, attr):
        return np.array([getattr(self.materials[m], attr)
                         for m in self.mat_index])

    def sigma(self, omega):
        st = self.cell_array("total").astype(complex)
        iv = self.cell_array("inv_v")
        return st + 1j * omega * iv  # (N, G)

    def coupling(self, omega, adjoint=False):
        """Per-cell scattering S[i, g_out, g_in] and production
        P[i, g_out, g_in] (not divided by k)."""
        S = np.empty((self.N, self.G, self.G))
        P = np.empty((self.N, self.G, self.G), complex)
        cache = {}
        for j, m in enumerate(self.materials):
            s = m.scatter.T  # [g_out, g_in]
            p = chi_omega(m, omega)
            cache[j] = (s.T, p.T) if adjoint else (s, p)
        for i, m in enumerate(self.mat_index):
            S[i], P[i] = cache[m]
        return S, P

    # --------------------------------------------------------------- sweep
    def sweep(self, sig, q_ang, psiL, psiR, keep_angular=False):
        """Transport sweep for one group, vectorised over right-hand sides.

        sig : (N,) complex total
        q_ang : (nrhs, N, M) angular emission density per unit mu
        psiL : (nrhs, M2) incoming at left (mu > 0)
        psiR : (nrhs, M2) incoming at right (mu < 0)
        Returns phi (nrhs, N), psiR_out (nrhs, M2), psiL_out (nrhs, M2)
        and optionally cell-average psi (nrhs, N, M).
        """
        nrhs = q_ang.shape[0]
        N, M2 = self.N, self.M2
        phi = np.zeros((nrhs, N), complex)
        psi_avg = np.zeros((nrhs, N, 2 * M2), complex) if keep_angular else None
        inv_mu = 1.0 / self.mu_p
        psi = psiL.astype(complex).copy()
        for i in range(N):
            t = sig[i] * self.h[i] * inv_mu
            c = self.h[i] * inv_mu
            out = (psi * (1 - 0.5 * t) + q_ang[:, i, :M2] * c) / (1 + 0.5 * t)
            avg = 0.5 * (psi + out)
            phi[:, i] += avg @ self.w_p
            if keep_angular:
                psi_avg[:, i, :M2] = avg
            psi = out
        psiR_out = psi
        psi = psiR.astype(complex).copy()
        for i in range(N - 1, -1, -1):
            t = sig[i] * self.h[i] * inv_mu
            c = self.h[i] * inv_mu
            out = (psi * (1 - 0.5 * t) + q_ang[:, i, M2:] * c) / (1 + 0.5 * t)
            avg = 0.5 * (psi + out)
            phi[:, i] += avg @ self.w_p
            if keep_angular:
                psi_avg[:, i, M2:] = avg
            psi = out
        return phi, psiR_out, psi, psi_avg

    def _bc_map(self, psiR_out, psiL_out):
        """New incoming boundary fluxes [left_in, right_in] from outgoing."""
        left = psiL_out if self.bc[0] == "reflective" else 0 * psiL_out
        right = psiR_out if self.bc[1] == "reflective" else 0 * psiR_out
        return np.concatenate([left, right], axis=1)

    def response(self, sig_g):
        """Response matrix for one group: columns are unit isotropic cell
        sources (N) then unit incoming boundary fluxes (2 M2); rows are
        [phi (N), new incoming boundary (2 M2)]."""
        N, M2 = self.N, self.M2
        nb = 2 * M2
        n = N + nb
        q = np.zeros((n, N, 2 * M2))
        q[np.arange(N), np.arange(N), :] = 0.5  # isotropic unit density
        psiL = np.zeros((n, M2))
        psiR = np.zeros((n, M2))
        psiL[N + np.arange(M2), np.arange(M2)] = 1.0
        psiR[N + M2 + np.arange(M2), np.arange(M2)] = 1.0
        phi, pRo, pLo, _ = self.sweep(sig_g, q, psiL, psiR)
        b = self._bc_map(pRo, pLo)
        return np.concatenate([phi, b], axis=1).T  # (n, n)

    # ------------------------------------------------------------- systems
    def _blocks(self, omega):
        sig = self.sigma(omega)
        return [self.response(sig[:, g]) for g in range(self.G)]

    def operators(self, omega, adjoint=False):
        """Return (Ks, Kf) with x = Ks x + Kf x / k + b."""
        N, G = self.N, self.G
        nb = 2 * self.M2
        ng = N + nb
        T = self._blocks(omega)
        S, P = self.coupling(omega, adjoint)
        Ks = np.zeros((G * ng, G * ng), complex)
        Kf = np.zeros((G * ng, G * ng), complex)
        for go in range(G):
            r = slice(go * ng, (go + 1) * ng)
            Tq = T[go][:, :N]
            Ks[r, go * ng + N:(go + 1) * ng] = T[go][:, N:]
            for gi in range(G):
                c = slice(gi * ng, gi * ng + N)
                Ks[r, c] += Tq * S[:, go, gi][None, :]
                Kf[r, c] += Tq * P[:, go, gi][None, :]
        return Ks, Kf

    def unpack(self, x):
        ng = self.N + 2 * self.M2
        x = x.reshape(self.G, ng)
        return x[:, :self.N].T.copy(), x[:, self.N:].T.copy()  # (N,G), (nb,G)

    def eigenvalue(self, adjoint=False):
        Ks, Kf = self.operators(0.0, adjoint)
        Ks, Kf = Ks.real, Kf.real
        A = np.eye(Ks.shape[0]) - Ks
        lu = sla.lu_factor(A)
        op = LinearOperator(A.shape, matvec=lambda v: sla.lu_solve(lu, Kf @ v),
                            dtype=float)
        vals, vecs = eigs(op, k=1, which="LM", tol=1e-13, maxiter=10000)
        k = vals[0].real
        x = vecs[:, 0].real
        phi, b = self.unpack(x)
        s = np.sign(phi.sum())
        phi, b = phi * s, b * s
        self.k = k
        self.phi0 = phi
        self.b0 = b
        self.adjoint_static = adjoint
        return k, phi

    def static_angular(self):
        """Cell-average angular static flux psi0 (N, M, G)."""
        S, P = self.coupling(0.0, self.adjoint_static)
        q = (np.einsum("igh,ih->ig", S, self.phi0)
             + np.einsum("igh,ih->ig", P.real, self.phi0) / self.k)
        sig = self.sigma(0.0)
        psi = np.zeros((self.N, 2 * self.M2, self.G))
        for g in range(self.G):
            qa = np.repeat(0.5 * q[None, :, g, None], 2 * self.M2, axis=2)
            _, _, _, pa = self.sweep(sig[:, g], qa, self.b0[None, :self.M2, g],
                                     self.b0[None, self.M2:, g],
                                     keep_angular=True)
            psi[:, :, g] = pa[0].real
        return psi

    def solve(self, omega, q_ang, k=None, adjoint=False, return_angular=False):
        """Solve the noise problem for angular source q_ang (N, M, G) per
        unit mu.  Returns dphi (N, G) and optionally dpsi (N, M, G)."""
        if k is None:
            k = self.k
        N, G, M2 = self.N, self.G, self.M2
        ng = N + 2 * M2
        Ks, Kf = self.operators(omega, adjoint)
        sig = self.sigma(omega)
        rhs = np.zeros(G * ng, complex)
        zero = np.zeros((1, M2))
        for g in range(G):
            phi, pRo, pLo, _ = self.sweep(sig[:, g], q_ang[None, :, :, g],
                                          zero, zero)
            rhs[g * ng:g * ng + N] = phi[0]
            rhs[g * ng + N:(g + 1) * ng] = self._bc_map(pRo, pLo)[0]
        A = np.eye(G * ng) - Ks - Kf / k
        x = np.linalg.solve(A, rhs)
        dphi, b = self.unpack(x)
        if not return_angular:
            return dphi
        S, P = self.coupling(omega, adjoint)
        q = (np.einsum("igh,ih->ig", S, dphi)
             + np.einsum("igh,ih->ig", P, dphi) / k)
        dpsi = np.zeros((N, 2 * M2, G), complex)
        for g in range(G):
            qa = q_ang[None, :, :, g] + 0.5 * q[None, :, g, None]
            _, _, _, pa = self.sweep(sig[:, g], qa, b[None, :M2, g],
                                     b[None, M2:, g], keep_angular=True)
            dpsi[:, :, g] = pa[0]
        if adjoint:
            # Solver direction mu corresponds to adjoint direction -mu
            dpsi = np.concatenate([dpsi[:, M2:], dpsi[:, :M2]], axis=1)
        return dphi, dpsi

    # ------------------------------------------------------------- sources
    def isotropic_source(self, region, spectrum):
        """Isotropic source density (per cm^3) `spectrum` (G,) on cells
        whose centres lie in region=(xlo, xhi); returned per unit mu.  Cells
        cut by the region edges get the overlapping fraction."""
        xlo, xhi = region
        lo = np.maximum(self.edges[:-1], xlo)
        hi = np.minimum(self.edges[1:], xhi)
        frac = np.clip(hi - lo, 0, None) / self.h
        q = 0.5 * frac[:, None, None] * np.asarray(spectrum, complex)[None, None, :]
        return np.repeat(q, 2 * self.M2, axis=1)

    def perturbation_source(self, omega, region, d_total=None, d_scatter=None,
                            d_nu_fission_scale=None):
        """-dSigma_t psi0 + dS^T phi0 / 2 + P(w) dnuSf phi0 / (2k) on cells
        in region (overlap-weighted).  d_total (G,), d_scatter (G, G)
        [g_in, g_out], d_nu_fission_scale scalar relative change."""
        xlo, xhi = region
        lo = np.maximum(self.edges[:-1], xlo)
        hi = np.minimum(self.edges[1:], xhi)
        frac = np.clip(hi - lo, 0, None) / self.h
        psi0 = self.static_angular()
        G = self.G
        q = np.zeros((self.N, 2 * self.M2, G), complex)
        if d_total is not None:
            q -= np.asarray(d_total)[None, None, :] * psi0
        iso = np.zeros((self.N, G), complex)
        if d_scatter is not None:
            iso += self.phi0 @ np.asarray(d_scatter)
        if d_nu_fission_scale is not None:
            _, P = self.coupling(omega)
            iso += d_nu_fission_scale * np.einsum("igh,ih->ig", P,
                                                  self.phi0) / self.k
        q += 0.5 * iso[:, None, :]
        return q * frac[:, None, None]

    def region_average(self, values, region):
        lo = np.maximum(self.edges[:-1], region[0])
        hi = np.minimum(self.edges[1:], region[1])
        fr = np.clip(hi - lo, 0, None)
        return np.tensordot(fr, values, axes=(0, 0)) / fr.sum()


def critical_nu_fission(make_slab, lo, hi, tol=1e-12):
    """Find a scalar s in [lo, hi] with k(make_slab(s)) = 1 by secant."""
    def f(s):
        sl = make_slab(s)
        k, _ = sl.eigenvalue()
        return k - 1.0
    a, b = lo, hi
    fa, fb = f(a), f(b)
    for _ in range(60):
        c = b - fb * (b - a) / (fb - fa)
        fc = f(c)
        if abs(fc) < tol:
            return c
        a, fa, b, fb = b, fb, c, fc
    return c


# --------------------------------------------------------------------- tests
def _test_infinite(ok):
    from ref_inf import cortex_fuel, ferrero_case, solve as inf_solve, static_k
    for name, mat, spec in [("Ferrero 1g", ferrero_case(), [1.0]),
                            ("CORTEX fuel 2g", cortex_fuel(), [0.3, 1.0])]:
        sl = Slab([(10.0, mat, 20)], bc=("reflective", "reflective"),
                  n_angles=8)
        k, phi0 = sl.eigenvalue()
        kinf, _ = static_k(mat)
        ok &= stats_report(f"{name} k vs k_inf", k, kinf, 1e-10)
        q = sl.isotropic_source((0, 10), spec)
        for f in [0.1, 1, 10, 100, 1e4]:
            w = 2 * np.pi * f
            for adj in (False, True):
                d = sl.solve(w, q, adjoint=adj)
                ref = inf_solve(mat, w, spec, kinf, adjoint=adj)
                err = np.abs(d - ref[None, :]).max() / np.abs(ref).max()
                good = err < 1e-6
                ok &= good
                print(f"  [{'PASS' if good else 'FAIL'}] {name} f={f:g} "
                      f"{'adj' if adj else 'fwd'}: max rel err {err:.2e}")
    return ok


def one_group_4a(nuSf=None, n_cells=60, n_angles=16, width=60.0):
    """Stage 4a slab material: Sigma_t=1, Sigma_s=0.9, v=2.2e5 cm/s,
    beta=0.0065, lambda=0.0784."""
    def mat(s):
        return Material.from_beta([1.0], [[0.9]], [1 / 2.2e5], [s], [1.0],
                                  beta=[0.0065], lam=[0.0784])
    if nuSf is None:
        nuSf = critical_nu_fission(
            lambda s: Slab([(width, mat(s), 60)], x0=-width / 2,
                           n_angles=16), 0.1, 0.105)
    return Slab([(width, mat(nuSf), n_cells)], x0=-width / 2,
                n_angles=n_angles), nuSf


def _test_convergence(ok):
    """Mesh halving and S32 vs S64 on the Stage 4a problem."""
    t0 = time.time()
    _, nuSf = one_group_4a()
    print(f"  Stage 4a critical nuSf (60 cells, S16) = {nuSf:.10f}")
    results = {}
    for nc, na in [(240, 32), (480, 32), (240, 64)]:
        sl, _ = one_group_4a(nuSf, nc, na)
        k, _ = sl.eigenvalue()
        q = sl.isotropic_source((14, 15), [1.0])
        res = {}
        for f in [0.1, 1, 10, 100]:
            d = sl.solve(2 * np.pi * f, q)[:, 0]
            # average onto 1 cm bins
            res[f] = d.reshape(60, -1).mean(axis=1)
        results[(nc, na)] = (k, res)
        print(f"    cells={nc} S{na}: k={k:.8f} ({time.time()-t0:.0f} s)")
    for label, other in [("mesh halving", (480, 32)), ("S32 vs S64", (240, 64))]:
        for f in [0.1, 1, 10, 100]:
            a = results[(240, 32)][1][f]
            b = results[other][1][f]
            err = np.abs(a - b).max() / np.abs(b).max()
            good = err < 1e-3
            ok &= good
            print(f"  [{'PASS' if good else 'FAIL'}] {label} f={f:g} Hz: "
                  f"max |diff|/max = {err:.2e}")
    return ok


def _test_diffusion(ok):
    from ref_diff1d import slab_noise
    a = 50.0
    D = 1 / 3.0
    sa = 0.01
    ext = 0.7104 / 1.0
    Be = np.pi / (2 * (a + ext))
    nuSf_d = sa + D * Be**2  # diffusion-critical
    mat = Material.from_beta([1.0], [[1 - sa]], [1 / 2.2e5], [nuSf_d], [1.0],
                             beta=[0.0065], lam=[0.0784])
    sl = Slab([(2 * a, mat, 400)], x0=-a, n_angles=32)
    k, _ = sl.eigenvalue()
    print(f"  thick slab: transport k = {k:.6f}")
    src = (14.0, 15.0)
    q = sl.isotropic_source(src, [1.0])
    for f in [1.0, 100.0]:
        w = 2 * np.pi * f
        d = sl.solve(w, q)[:, 0]
        dd = slab_noise(sl.x, w, a + ext, D, sa, nuSf_d, 1 / 2.2e5,
                        [0.0065], [0.0784], src, 1.0, k=None)
        mask = (np.abs(sl.x - 14.5) > 5) & (np.abs(sl.x) < a - 5)
        err = (np.abs(d - dd) / np.abs(dd))[mask].max()
        good = err < 0.02
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] SN vs diffusion f={f:g} Hz:"
              f" max rel err {err:.2e} (away from source and boundary)")
    return ok


def self_test(which):
    ok = True
    print("ref_sn1d self-test")
    if "inf" in which:
        ok = _test_infinite(ok)
    if "diff" in which:
        ok = _test_diffusion(ok)
    if "conv" in which:
        ok = _test_convergence(ok)
    return bool(ok)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--test", nargs="*", default=None,
                   help="run self-tests: inf diff conv (default all)")
    args = p.parse_args()
    if args.test is not None:
        which = args.test or ["inf", "diff", "conv"]
        raise SystemExit(0 if self_test(which) else 1)


if __name__ == "__main__":
    main()
