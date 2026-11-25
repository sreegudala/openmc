"""8c n = 1 at the vibrating pin: half-pin detectors (thermal), Abeille map
sum, OpenMC map sums at each cancellation setting, and the adjoint
contraction with each forward run's static moments. Run from noise_dev/."""
import sys
sys.path[:0] = ["stage8", "stage6", "common", "reference"]
sys.argv = sys.argv[:1]
import numpy as np  # noqa: E402
import run_stage8c as S  # noqa: E402

c = slice(68, 70)
s = np.load("stage6/stage6c_abeille_static.npz")
st = np.transpose(s["avg"], (0, 2, 1))
B = np.transpose(np.load("stage8/stage8c_abeille_n1.npz")["vals"],
                 (0, 1, 3, 2)) / st[1, c, c].mean()


def ms(v):
    return v.mean(), np.hypot(v.real.std(ddof=1),
                              v.imag.std(ddof=1)) / np.sqrt(len(v))


for det in ("pin_top", "pin_bottom"):
    (ix0, ix1), (iy0, iy1) = S.det_bins(det)
    adj = np.load(f"stage8/stage8c_adjoint_n1_{det}_probe.npz")
    print(f"== {det}")
    m, s_ = ms(B[:, 1, iy0:iy1, ix0:ix1].sum(axis=(1, 2)))
    print(f"   Abeille map sum {m:.4f} +- {s_:.4f}")
    for tag in ("c0_diag", "c170_diag", "c340_diag", "c170_probe"):
        A, a, norm = S.normalised_openmc(S.HERE /
                                         f"stage8c_openmc_n1_{tag}.npz")
        m, s_ = ms(A[:, 1, iy0:iy1, ix0:ix1].sum(axis=(1, 2)))
        mc, sc = ms(S.contraction(a, adj, 1, det).sum(1) / norm)
        print(f"   OpenMC {tag}: map sum {m:.4f} +- {s_:.4f}; contraction "
              f"{mc:.4f} +- {sc:.4f}")
