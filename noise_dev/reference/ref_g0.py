"""Zero-power transfer function G0(w) and collision-based cut-off frequencies.

G0(w) = 1 / ( i w ( Lambda + sum_k beta_k / (i w + lambda_k) ) )
(e^{-i w t} convention; Hansen & Hansson 2025 Eqs. 2.68-2.69.)

Cut-offs for the collision-based scheme (Fauvel et al. 2025):
    kappa = sqrt(1 + eta^2) - eta
    w_CC  ~ beta / (kappa Lambda)          cancellation needed above
    w_D   = 1 / (kappa Lambda)             divergent above
    w_VLF ~ 2 kappa lambda^2 Lambda / beta (one family; lambda_eff otherwise)
Extended generation radius: K~cb = Upsilon~ / (1 - kappa w Lambda), with
Upsilon~ = 1 - beta + sum_k beta_k lambda_k / |lambda_k + i w|.
"""

import argparse

import numpy as np

# Keepin U-235 thermal six-family decay constants [1/s].  The TRIGA beta_k
# come from the user's IFP run (results/kinetics/kinetics.json); the
# lambda_k were not tallied there, so Keepin values are used.
KEEPIN_U235_LAMBDA = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01])
TRIGA_BETA_K = np.array([0.0002644705355382282, 0.0013058963004798163,
                         0.001225756664711052, 0.00273594574264271,
                         0.0011547702801722568, 0.0004637971873227474])
TRIGA_LAMBDA = 4.735711965107854e-05


def g0(omega, beta_k, lam_k, Lambda):
    omega = np.asarray(omega, float)
    beta_k = np.atleast_1d(beta_k)
    lam_k = np.atleast_1d(lam_k)
    iw = 1j * omega[..., None]
    s = (beta_k / (iw + lam_k)).sum(axis=-1)
    return 1.0 / (1j * omega * (Lambda + s))


def upsilon(omega, beta_k, lam_k):
    beta_k = np.atleast_1d(beta_k)
    lam_k = np.atleast_1d(lam_k)
    iw = 1j * np.asarray(omega, float)[..., None]
    return 1.0 - beta_k.sum() + (lam_k * beta_k / (lam_k + iw)).sum(axis=-1)


def upsilon_abs(omega, beta_k, lam_k):
    beta_k = np.atleast_1d(beta_k)
    lam_k = np.atleast_1d(lam_k)
    w = np.asarray(omega, float)[..., None]
    return 1.0 - beta_k.sum() + (lam_k * beta_k / np.hypot(lam_k, w)).sum(-1)


def kappa(eta):
    return np.sqrt(1.0 + eta**2) - eta


def cutoffs(eta, beta, Lambda, lam_eff):
    k = kappa(eta)
    return {"kappa": k,
            "omega_CC": beta / (k * Lambda),
            "omega_D": 1.0 / (k * Lambda),
            "omega_VLF": 2 * k * lam_eff**2 * Lambda / beta}


def ktilde_cb(omega, eta, beta_k, lam_k, Lambda):
    """Extended (absolute-weight) generation radius, collision-based."""
    return upsilon_abs(omega, beta_k, lam_k) / (1 - kappa(eta) * omega * Lambda)


def ctilde_cb(omega, eta, c_inf, tau):
    """Extended within-generation radius, collision-based.
    tau = 1/(v Sigma_t)."""
    wt = omega * tau
    return (c_inf + np.sqrt(1 + eta**2) * wt) / (1 + eta * wt)


def lambda_eff(beta_k, lam_k):
    beta_k = np.atleast_1d(beta_k)
    return beta_k.sum() / (beta_k / np.atleast_1d(lam_k)).sum()


def self_test():
    ok = True
    print("ref_g0 self-test")
    # Design doc Section 1.5 table: one family lambda=0.08, beta=715 pcm,
    # Lambda=47.4 us.
    table = {0.1: (140.9, -7.5), 1: (139.7, -3.1), 10: (129.1, -22.7),
             24: (98.9, -45.0), 100: (32.6, -76.5)}
    for f, (amp, ph) in table.items():
        g = g0(2 * np.pi * f, [715e-5], [0.08], 47.4e-6)
        a, p = abs(g), np.degrees(np.angle(g))
        good = abs(a - amp) < 0.06 and abs(p - ph) < 0.06
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] f={f:5} Hz |G0|={a:7.2f} "
              f"(doc {amp}) arg={p:6.2f} (doc {ph})")
    # Cut-off table, Section 2.3
    for eta, fcc in [(0.5, 38.8), (1, 58.0), (2, 101.7), (3, 147.9)]:
        c = cutoffs(eta, 715e-5, 47.4e-6, 0.08)
        f = c["omega_CC"] / (2 * np.pi)
        good = abs(f - fcc) < 0.1
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] eta={eta}: f_CC={f:.1f} Hz "
              f"(doc {fcc})")
    kt = ktilde_cb(2 * np.pi * 100, 1.0, [715e-5], [0.08], 47.4e-6)
    good = abs(kt - 1.0053) < 2e-4
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] K~cb(100 Hz, eta=1)={kt:.5f} "
          "(doc 1.0053)")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--test", action="store_true", help="run self-test")
    p.add_argument("--triga", action="store_true",
                   help="print G0 and cut-offs for the TRIGA kinetics")
    args = p.parse_args()
    if args.test:
        raise SystemExit(0 if self_test() else 1)
    if args.triga:
        b, lam, L = TRIGA_BETA_K, KEEPIN_U235_LAMBDA, TRIGA_LAMBDA
        print(f"beta={b.sum():.6f} Lambda={L:.4e} s lambda_eff="
              f"{lambda_eff(b, lam):.4f}/s")
        for f in [0.1, 0.3, 1, 3, 10, 30, 100]:
            g = g0(2 * np.pi * f, b, lam, L)
            print(f"  f={f:6} Hz |G0|={abs(g):8.2f} arg={np.degrees(np.angle(g)):7.2f}"
                  f" K~cb(eta=1)={ktilde_cb(2*np.pi*f, 1, b, lam, L):.5f}"
                  f" K~cb(eta=3)={ktilde_cb(2*np.pi*f, 3, b, lam, L):.5f}")
        for eta in [0.5, 1, 2, 3, 5]:
            c = cutoffs(eta, b.sum(), L, lambda_eff(b, lam))
            print(f"  eta={eta}: f_CC={c['omega_CC']/2/np.pi:7.1f} Hz "
                  f"f_D={c['omega_D']/2/np.pi:9.1f} Hz")


if __name__ == "__main__":
    main()
