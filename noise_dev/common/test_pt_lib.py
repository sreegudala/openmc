"""Sync test: the copied pt_lib.pt_terms reproduces ce_mg_triga's, and the
complex wrapper pt_terms_complex is exact.

    python -m pytest noise_dev/common/test_pt_lib.py

"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import pt_lib  # noqa: E402
from adjoint_noise_lib import pt_terms_complex  # noqa: E402

FIXTURE = HERE / "fixtures" / "pt_lib" / "sync.npz"
KEYS = ("F", "A", "w", "dSt", "dN0", "dN1")


def test_copy_matches_triga():
    d = np.load(FIXTURE)
    n = len({k.split("_")[0] for k in d.files})
    assert n == 2
    for i in range(n):
        args = {k: d[f"{i}_{k}"] for k in KEYS}
        assert np.array_equal(np.array(pt_lib.pt_terms(**args)),
                              d[f"{i}_out"]), i


def _direct(F, A, w, dSt, dN0, dN1):
    """The same contraction written out in complex arithmetic."""
    dSt = np.broadcast_to(dSt, F.shape[:2]) if np.ndim(dSt) == 1 else dSt
    gg = "gh" if np.ndim(dN0) == 2 else "vgh"
    F0, A0 = F[:, :, 0] * w[:, None], A[:, :, 0]
    F1, A1 = F[:, :, 1:] * w[:, None, None], A[:, :, 1:]
    p0 = -(np.sum(F0 * A0 * dSt) - np.einsum(f"vg,{gg},vh->", F0, dN0, A0))
    p1 = 3 * (np.sum(F1 * A1 * dSt[:, :, None])
              - np.einsum(f"vgm,{gg},vhm->", F1, dN1, A1))
    return p0, p1


def test_complex_wrapper():
    rng = np.random.default_rng(7)
    V, G = 3, 2
    c = lambda *s: rng.random(s) - 0.5 + 1j * (rng.random(s) - 0.5)
    F = rng.random((V, G, 4))
    w = rng.random(V)
    for per_voxel in (False, True):
        s = (V, G) if per_voxel else (G,)
        args = (F, c(V, G, 4), w, c(*s), c(*s, G), c(*s, G))
        got = pt_terms_complex(*args)
        ref = _direct(*args)
        np.testing.assert_allclose(got, ref, rtol=1e-13, atol=0)
        # real inputs give the copied function's result
        real = tuple(np.real(a) for a in args)
        np.testing.assert_allclose(pt_terms_complex(*real),
                                   pt_lib.pt_terms(*real), rtol=1e-15)
