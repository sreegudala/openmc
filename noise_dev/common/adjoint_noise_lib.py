"""Adjoint library for the full noise operator (brief N2).

The scattering part comes from the copied, frozen ce_mg_triga transposition
(:func:`adjoint_lib.make_adjoint_library`); this module adds the transposed
fission and delayed-neutron data that the noise operator needs.

Fission form: **matrix**. With no ``chi`` dataset, OpenMC's MG reader takes
``prompt-nu-fission`` [g_in, g_out] and ``delayed-nu-fission``
[family, g_in, g_out] as production matrices (``XsData::
fission_matrix_no_beta_from_hdf5``). It samples the family from the per-g_in
row sums and the outgoing group from the normalised row. The noise fission
path (``noise::create_fission_sites``) reads only the total production
``nu_fission(g_in)``, i.e. the row sum of the transposed matrices, and
``sample_fission_energy``. Without delayed families a single ``nu-fission``
[g_in, g_out] matrix is written.

Transposition. The forward production matrix stored by OpenMC is
F[g_in, g_out] = chi_{g_out}(g_in) nuSigma_f(g_in) per family (prompt and
each delayed family k). The adjoint library stores F'[g_in, g_out] =
F[g_out, g_in], so

* adjoint production in group g = sum_g' F[g', g] = chi_g sum_g' nuSigma_f,g'
  for a rank-1 family: it is zero in groups with no fission spectrum;
* the outgoing (adjoint) group is sampled proportional to nuSigma_f,g';
* lambda_k / (lambda_k + i omega) is a scalar per family and is applied by
  the code as in forward mode; ``decay-rate`` is copied unchanged;
* ``inverse-velocity`` and ``total`` are unchanged (kept by the copied
  transposition).

Keys in the transposed library:

* ``prompt-nu-fission``, ``delayed-nu-fission`` (or ``nu-fission``): the
  transposed production matrices above. Rows with zero production make
  OpenMC's normalised chi row 0/0 (NaN); such rows are never sampled
  because their production is zero.
* ``absorption``: N1's balance quantity, max(0, total - in-scatter). It is
  not a reaction rate and does not contain fission: in the noise collision
  the continuing weight is (total - absorption)/total~ times the row
  multiplicity, i.e. in-scatter/total~, and fission sites are a separate
  branch with expected number nu_fission/(total~ k).
* ``fission``, ``kappa-fission``, ``chi``, ``beta`` and the vector forms are
  absent (removed by the copied code). Scores such as ``fission`` or
  ``absorption`` are therefore meaningless in adjoint runs; tally flux only.

Only isotropic representations are supported.
"""

import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adjoint_lib  # noqa: E402
import pt_lib  # noqa: E402


def _normalised(x, axis=-1):
    x = np.asarray(x, float)
    return x / x.sum(axis=axis, keepdims=True)


def _beta_by_group(beta, K, G):
    beta = np.asarray(beta, float)
    if beta.shape == (K,):
        return np.repeat(beta[:, None], G, axis=1)
    if beta.shape == (K, G):
        return beta
    raise ValueError(f"beta of shape {beta.shape} not supported")


def production_matrices(data, G, K):
    """Forward production matrices as OpenMC's MG reader builds them.

    Parameters
    ----------
    data : h5py.Group
        One temperature group of an XSdata entry (isotropic)
    G, K : int
        Numbers of energy groups and of delayed families (library level)

    Returns
    -------
    Mp : numpy.ndarray, shape (G, G)
        Prompt production [g_in, g_out]
    Md : numpy.ndarray, shape (K, G, G)
        Delayed production per family [k, g_in, g_out]

    """
    Md = np.zeros((K, G, G))
    if "chi" in data or "chi-prompt" in data:
        if K == 0:
            Mp = np.outer(data["nu-fission"][()], _normalised(data["chi"][()]))
        elif "beta" in data:
            chi = _normalised(data["chi"][()])
            nuf = data["nu-fission"][()]
            b = _beta_by_group(data["beta"][()], K, G)
            Mp = np.outer(nuf * (1.0 - b.sum(0)), chi)
            for k in range(K):
                Md[k] = np.outer(b[k] * nuf, chi)
        else:
            chi_p = _normalised(data["chi-prompt"][()])
            chi_d = _normalised(np.reshape(data["chi-delayed"][()], (K, G)))
            pnf = data["prompt-nu-fission"][()]
            dnf = np.reshape(data["delayed-nu-fission"][()], (K, G))
            Mp = np.outer(pnf, chi_p)
            for k in range(K):
                Md[k] = np.outer(dnf[k], chi_d[k])
    else:
        if K == 0:
            Mp = np.reshape(data["nu-fission"][()], (G, G)).copy()
        elif "beta" in data:
            M = np.reshape(data["nu-fission"][()], (G, G))
            b = _beta_by_group(data["beta"][()], K, G)
            Mp = (1.0 - b.sum(0))[:, None] * M
            for k in range(K):
                Md[k] = b[k][:, None] * M
        else:
            Mp = np.reshape(data["prompt-nu-fission"][()], (G, G)).copy()
            Md = np.reshape(data["delayed-nu-fission"][()], (K, G, G)).copy()
    return Mp, Md


def _check_isotropic(group):
    rep = group.attrs.get("representation", b"isotropic")
    rep = rep.decode() if isinstance(rep, bytes) else rep
    if rep != "isotropic":
        raise ValueError(f"{group.name}: representation {rep} not supported")


def make_adjoint_noise_library(src, dst):
    """Transposed MGXS library for the adjoint noise operator.

    Parameters
    ----------
    src : path-like
        Forward MGXS library
    dst : path-like
        Output path

    Returns
    -------
    dict
        The report of :func:`adjoint_lib.make_adjoint_library` plus
        ``fission_form`` ("matrix") and ``fissionable`` (material names)

    """
    fission = {}
    with h5py.File(src, "r") as f:
        G = int(f.attrs["energy_groups"])
        K = int(f.attrs.get("delayed_groups", 0))
        for name, group in f.items():
            if not isinstance(group, h5py.Group):
                continue
            _check_isotropic(group)
            if not group.attrs.get("fissionable", False):
                continue
            for key, data in group.items():
                if key == "kTs" or not isinstance(data, h5py.Group):
                    continue
                Mp, Md = production_matrices(data, G, K)
                lam = data["decay-rate"][()] if K > 0 else None
                fission[name, key] = (Mp, Md, lam)

    report = adjoint_lib.make_adjoint_library(src, dst)

    with h5py.File(dst, "r+") as f:
        for (name, key), (Mp, Md, lam) in fission.items():
            data = f[name][key]
            if K == 0:
                data.create_dataset("nu-fission", data=Mp.T)
            else:
                data.create_dataset("prompt-nu-fission", data=Mp.T)
                data.create_dataset("delayed-nu-fission",
                                    data=Md.transpose(0, 2, 1))
                data.create_dataset("decay-rate", data=lam)
            f[name].attrs["fissionable"] = True
    report["fission_form"] = "matrix"
    report["fissionable"] = sorted({n for n, _ in fission})
    return report


def effective_operator(path, material, temperature=None):
    """Infinite-medium operator that OpenMC transports for one material.

    Scattering is what the MG collision actually does: scatter branch
    (total - absorption)/total, outgoing group from the scatter matrix
    divided by the multiplicity, weight times the multiplicity.

    Returns
    -------
    dict
        ``total`` (G,), ``scatter`` (G, G) [g_in, g_out] expected transfers,
        ``Mp`` (G, G), ``Md`` (K, G, G) [g_in, g_out], ``decay_rate`` (K,),
        ``inv_v`` (G,)

    """
    with h5py.File(path, "r") as f:
        G = int(f.attrs["energy_groups"])
        K = int(f.attrs.get("delayed_groups", 0))
        group = f[material]
        nb, fmt = adjoint_lib._nbins(group)
        key = temperature or next(k for k in group if k != "kTs")
        data = group[key]
        N, M = adjoint_lib._dense(data["scatter_data"], G, nb)
        N0 = N[..., 0] if fmt == "legendre" else N.sum(2)
        total = data["total"][()]
        absorption = data["absorption"][()]
        plain = np.where(M > 0, N0 / np.where(M > 0, M, 1), 0.0)
        row = plain.sum(1)
        p_out = plain / np.where(row > 0, row, 1)[:, None]
        scatter = (total - absorption)[:, None] * p_out * M
        if group.attrs.get("fissionable", False):
            Mp, Md = production_matrices(data, G, K)
            lam = data["decay-rate"][()] if K > 0 else np.zeros(0)
        else:
            Mp, Md, lam = np.zeros((G, G)), np.zeros((K, G, G)), np.zeros(K)
        inv_v = data["inverse-velocity"][()]
    return {"total": total, "absorption": absorption, "scatter": scatter,
            "p_out": p_out, "multiplicity": M, "Mp": Mp, "Md": Md,
            "decay_rate": lam, "inv_v": inv_v}


def second_moment_radius(op):
    """Spectral radius of the per-collision second-moment operator.

    Scattering only, at omega = 0, with implicit capture: a collision in
    group g moves the particle to g' with probability p(g, g') and multiplies
    its weight by f(g, g') = (1 - Sigma_a/Sigma_t)_g M(g, g'). The expected
    squared weight then evolves with Q = p f^2. With rho(Q) >= 1 the second
    moment of a history's score is infinite in an infinite medium (no
    leakage, no splitting): sample means converge slowly and erratically and
    replica variances are meaningless. Forward libraries (M = 1) always have
    rho(Q) < 1; transposed ones with multiplicity > 1 need not.
    """
    f = (1.0 - op["absorption"] / op["total"])[:, None] * op["multiplicity"]
    Q = op["p_out"] * f**2
    return float(np.max(np.abs(np.linalg.eigvals(Q))))


def infinite_k(op):
    """Dominant eigenvalue of the infinite-medium static operator."""
    F = op["Mp"] + op["Md"].sum(0)
    A = np.diag(op["total"]) - op["scatter"].T
    w = np.linalg.eigvals(np.linalg.solve(A, F.T))
    return float(w[np.argmax(w.real)].real)


def pt_terms_complex(F, A, w, dSt, dN0, dN1):
    """pt_lib.pt_terms (copied from ce_mg_triga, P0 and P1 parts of the
    first-order detector change) for a complex adjoint A and complex
    perturbation (dSt, dN0, dN1); F (static forward moments) and w are real.

    pt_terms casts to float, so it is called on real and imaginary parts:
    the contraction is bilinear in A and in the perturbation, hence
    P(F, Ar + i Ai, Xr + i Xi) = P(Ar, Xr) - P(Ai, Xi)
    + i (P(Ai, Xr) + P(Ar, Xi)).

    Returns
    -------
    (P0, P1) : complex
    """
    F = np.asarray(F)
    if np.iscomplexobj(F) and np.any(np.imag(F) != 0):
        raise ValueError("F (static forward moments) must be real")
    F = np.real(F)
    w = np.real(np.asarray(w))
    A = np.asarray(A, complex)
    X = [np.asarray(x, complex) for x in (dSt, dN0, dN1)]
    Xr = [x.real for x in X]
    Xi = [x.imag for x in X]

    def P(a, x):
        return np.array(pt_lib.pt_terms(F, a, w, *x))

    out = (P(A.real, Xr) - P(A.imag, Xi)
           + 1j * (P(A.imag, Xr) + P(A.real, Xi)))
    return complex(out[0]), complex(out[1])
