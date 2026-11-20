# Copied from ce_mg_triga (TRIGA proposal repository):
#   step5e_pool_pt.py lines 243-269  pt_terms
# Source commit: 540ab07832bafeea48927861e14ea34d4531ddd7 (2026-11-20)
# Copied: 2026-11-20
#
# Copied verbatim. Do not edit here; extend in adjoint_noise_lib.py.
# Changes go into ce_mg_triga first and are then copied over again (refresh
# the fixture with it; see test_pt_lib.py and fixtures/pt_lib/).
#
# Only the import block below is new. pt_terms is the first-order
# perturbation contraction of a forward and an adjoint run, P0 and P1 parts
# (design doc, Stage 7 amendment, test 7c). It casts its inputs to float:
# for the complex adjoint of the noise problem use
# adjoint_noise_lib.pt_terms_complex, which calls it on the real and
# imaginary parts.

import numpy as np


def pt_terms(F, A, w, dSt, dN0, dN1):
    """First-order change of a detector response, P0 and P1 parts. Pure (brief v20 G2).

    F, A    (V, G, 4) volume-integrated forward and adjoint moments per voxel, library
            order (highest energy first); moment 0 is the scalar flux, 1..3 the Y1m (current)
            components. A is the adjoint run's tally times its source strength S.
    w       (V,) weight per voxel: f/V for a product of voxel averages, 1/(f V) for
            water-only moments (Phi_w Phi+_w / V_w), 0 to skip.
    dSt     (G,) or (V, G) change in the total cross section
    dN0     (G, G) change in the P0 nu-scatter matrix g -> g' (or (V, G, G))
    dN1     (G, G) change in its first Legendre moment (or (V, G, G))
    Returns (P0, P1):
      P0 = -sum_v w [ sum_g dSt F0 A0 - sum_gg' F0_g dN0(g->g') A0_g' ]
      P1 = +3 sum_v w sum_m [ sum_g dSt F_m A_m - sum_gg' F_m,g dN1(g->g') A_m,g' ]
    The +3 is the Legendre factor (2l+1) with the sign flip of the odd adjoint moment
    (psi+(Omega) = 4 pi psi_MC(-Omega), see the module docstring).
    """
    F, A, w = np.asarray(F, float), np.asarray(A, float), np.asarray(w, float)
    dSt = np.asarray(dSt, float)
    dSt = np.broadcast_to(dSt, F.shape[:2]) if dSt.ndim == 1 else dSt
    gg = "gh" if np.ndim(dN0) == 2 else "vgh"
    F0, A0 = F[:, :, 0] * w[:, None], A[:, :, 0]
    F1, A1 = F[:, :, 1:] * w[:, None, None], A[:, :, 1:]
    p0 = -(np.sum(F0 * A0 * dSt) - np.einsum(f"vg,{gg},vh->", F0, dN0, A0))
    p1 = 3 * (np.sum(F1 * A1 * dSt[:, :, None]) - np.einsum(f"vgm,{gg},vhm->", F1, dN1, A1))
    return float(p0), float(p1)

