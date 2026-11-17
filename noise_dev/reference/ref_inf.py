"""Multigroup infinite-medium noise reference (design doc Section 3.3).

Forward:  [diag(Sigma_t + i w/v) - S - P(w)/k] dphi = dS
Adjoint:  the plain (bilinear) transpose of the same matrix at the same +i w.

k is the dominant eigenvalue of the same matrices at w = 0 unless given.
"""

import argparse

import numpy as np

from noise_common import Material, chi_omega, scatter_out_in, stats_report
from ref_g0 import g0, upsilon


def static_k(mat):
    A = np.diag(mat.total) - scatter_out_in(mat)
    P = chi_omega(mat, 0.0).real
    M = np.linalg.solve(A, P)
    w, v = np.linalg.eig(M)
    i = np.argmax(w.real)
    phi = np.abs(v[:, i].real)
    return w[i].real, phi


def operator(mat, omega, k):
    return (np.diag(mat.total + 1j * omega * mat.inv_v)
            - scatter_out_in(mat) - chi_omega(mat, omega) / k)


def solve(mat, omega, source, k=None, adjoint=False):
    """Return dphi for an isotropic source density `source` (G,)."""
    if k is None:
        k, _ = static_k(mat)
    A = operator(mat, omega, k)
    if adjoint:
        A = A.T
    return np.linalg.solve(A, np.asarray(source, complex))


def perturbation_source(mat, omega, k, phi0, d_total=None, d_scatter=None,
                        d_nu_fission_scale=None):
    """dS = -dSigma_t phi0 + dS^T phi0 + P(w) * scale * phi0 / k
    (scalar-flux form; exact in an infinite medium with isotropic flux)."""
    G = mat.G
    s = np.zeros(G, complex)
    if d_total is not None:
        s -= np.asarray(d_total) * phi0
    if d_scatter is not None:
        s += np.asarray(d_scatter).T @ phi0
    if d_nu_fission_scale is not None:
        s += chi_omega(mat, omega) @ (d_nu_fission_scale * phi0) / k
    return s


def ferrero_case():
    """MGMC case from the Ferrero thesis: Sigma_t=1, Sigma_s=0, Sigma_f=0.4,
    nu=2.5 so k=1, beta=7e-3, lambda=0.08, v=3e5 cm/s."""
    return Material.from_beta(total=[1.0], scatter=[[0.0]], inv_v=[1 / 3e5],
                              nu_fission=[2.5 * 0.4], chi=[1.0],
                              beta=[7e-3], lam=[0.08])


def cortex_fuel():
    """CORTEX 2-group fuel data (design doc Section 3.4, from Abeille)."""
    nu_p, nu_d = 2.576819862, 0.013860138
    sf = np.array([0.0057671, 0.10622])
    return Material(total=[0.3779, 0.55064],
                    scatter=[[0.35128029, 0.00086471], [0.0, 0.39276]],
                    inv_v=[1 / 1.82304e7, 1 / 4.13067e5],
                    nu_fission_p=nu_p * sf, nu_fission_d=[nu_d * sf],
                    chi_p=[1, 0], chi_d=[[1, 0]], lam=[0.0851])


def cortex_water():
    return Material(total=[0.25411, 1.2182],
                    scatter=[[0.22519143, 0.028124], [0.0, 1.188884]],
                    inv_v=[1 / 1.82304e7, 1 / 4.13067e5])


def self_test():
    ok = True
    print("ref_inf self-test")
    # Gate: Ferrero case at 1e4 Hz with fission amplitude pi*eps_f.
    mat = ferrero_case()
    k, phi0 = static_k(mat)
    phi0 = phi0 / phi0[0]
    ok &= stats_report("Ferrero k", k, 1.0, 1e-12)
    w = 2 * np.pi * 1e4
    src = perturbation_source(mat, w, k, phi0,
                              d_nu_fission_scale=np.pi * 5e-3)
    d = solve(mat, w, src, k)[0]
    print(f"  dphi/phi0 = {d.real:.4e} {d.imag:+.4e} i")
    ok &= stats_report("Ferrero Re", d.real, 2.486e-3, 2e-4)
    ok &= stats_report("Ferrero Im", d.imag, -7.439e-2, 2e-4)
    # Analytic check against G0 Upsilon (one group, one family)
    Lam = 1 / (3e5 * 1.0)
    ref = g0(w, [7e-3], [0.08], Lam) * upsilon(w, [7e-3], [0.08]) * np.pi * 5e-3
    ok &= stats_report("Ferrero vs G0*Ups*pi*eps", abs(d - ref) + 1, 1.0, 1e-10)
    # Published reference (rounded): 2.49e-3 - 7.44e-2 i
    ok &= stats_report("Ferrero vs published Re", d.real, 2.49e-3, 3e-3)
    ok &= stats_report("Ferrero vs published Im", d.imag, -7.44e-2, 3e-3)
    # CORTEX fuel k_inf and FEMFFUSION nuSigma_f_1 consistency
    fuel = cortex_fuel()
    kf, _ = static_k(fuel)
    print(f"  CORTEX fuel k_inf = {kf:.6f}")
    ok &= stats_report("CORTEX nuSf_1", fuel.nu_fission[0], 0.01494071, 1e-6)
    beta = fuel.nu_fission_d.sum() / fuel.nu_fission.sum()
    ok &= stats_report("CORTEX beta", beta, 0.00535, 1e-3)
    # Reciprocity of the transposed solve: <s2, A^-1 s1> = <A^-T s2, s1>
    rng = np.random.default_rng(1)
    s1, s2 = rng.random(2), rng.random(2)
    for f in [0.1, 1, 100]:
        w = 2 * np.pi * f
        lhs = s2 @ solve(fuel, w, s1, kf)
        rhs = solve(fuel, w, s2, kf, adjoint=True) @ s1
        ok &= stats_report(f"reciprocity {f} Hz",
                           abs(lhs - rhs) / abs(lhs) + 1, 1, 1e-12)
    # sesquilinear adjoint at w equals conj of bilinear adjoint at -w... and
    # forward at -w is conj of forward at w
    w = 2 * np.pi
    a = solve(fuel, w, s1, kf)
    b = solve(fuel, -w, s1, kf)
    ok &= stats_report("conj symmetry", abs(a - b.conj()).max() + 1, 1, 1e-12)
    return bool(ok)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--test", action="store_true")
    args = p.parse_args()
    if args.test:
        raise SystemExit(0 if self_test() else 1)


if __name__ == "__main__":
    main()
