"""1-D vibrating fuel plate (two opposite interfaces) in the water gap of
the 8a slab: geometry and S_N reference. Run from noise_dev/."""
import sys
sys.path[:0] = ["stage8", "common", "reference", "stage4"]; sys.argv = sys.argv[:1]
import numpy as np
import run_stage8 as S
from noise_common import chi_omega
from ref_sn1d import Slab
mix, wat = S.materials()
EPS, HP, XC = 0.2, 0.3657, 12.0
F1, F2 = XC - HP, XC + HP          # faces: (water | mix) at F1, (mix | water) at F2
REG = [(60.0, mix, 240), (F1 - EPS - 10.0, wat, 6), (EPS, wat, 20), (EPS, mix, 20),
       (2*HP - 2*EPS, mix, 20), (EPS, mix, 20), (EPS, wat, 20), (14.0 - F2 - EPS, wat, 6), (36.0, mix, 144)]
FACES = [(F1, wat, mix), (F2, mix, wat)]        # (x0, negative side, positive side)

def sn(n, n_angles=32):
    sl = Slab(REG, x0=-50.0, n_angles=n_angles); k, _ = sl.eigenvalue(); psi0 = sl.static_angular()
    w = 2*np.pi*n; q = np.zeros_like(psi0, dtype=complex)
    for x0, L, R in FACES:
        c = S.c_avg(n, sl.edges[:-1] - x0 + S.X0, sl.edges[1:] - x0 + S.X0, EPS)   # c_avg is centred at S.X0
        dT = L.total - R.total; dS = L.scatter - R.scatter; dP = chi_omega(L, w) - chi_omega(R, w)
        q += -c[:, None, None]*dT[None, None, :]*psi0
        q += 0.5*(c[:, None]*(sl.phi0 @ dS + np.einsum("gh,ih->ig", dP, sl.phi0)/k))[:, None, :]
    d = sl.solve(w, q)
    return sl, k, d

if __name__ == "__main__":
    for n in (1, 2):
        sl, k, d = sn(n)
        r = d[:, 1] / sl.phi0[:, 1]
        near = (sl.x > XC - 1.5) & (sl.x < XC + 1.5)
        D = (sl.x > 20) & (sl.x < 25)
        print(f"n={n} k={k:.6f}: thermal dphi/phi0 near plate min/max Im {r[near].imag.min():+.4f} {r[near].imag.max():+.4f}; "
              f"far (x=0) {r[np.argmin(abs(sl.x))]:.4f}; detector [20,25] {(d[D,1]*sl.h[D]).sum()/(sl.phi0[D,1]*sl.h[D]).sum():.4f}")
