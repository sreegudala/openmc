# Copied from ce_mg_triga (TRIGA proposal repository):
#   common.py         lines 448-484  FISSION_DATASETS, remove_fission
#   step3a_adjoint.py lines 69-206   _nbins, _dense, _write_sparse,
#                                    transpose_library, make_adjoint_library
# Source commit: 540ab07832bafeea48927861e14ea34d4531ddd7 (2026-11-20)
# Copied: 2026-11-22 (refresh; previous copy e018a83, 2026-11-20)
#
# Copied verbatim. Do not edit here; extend in adjoint_noise_lib.py.
# Changes go into ce_mg_triga first (logged in its CHANGELOG_transpose.md)
# and are then copied over again; refresh the fixtures with it (see
# test_adjoint_lib.py and fixtures/adjoint_lib/).
#
# Only the import block below is new. transpose_library is the pure
# forward-to-adjoint transformation of one material's scattering data, as
# used for the TRIGA source-to-detector matrix G; make_adjoint_library is
# the HDF5 wrapper that applies it to every material and removes fission:
#   - nu-scatter transposed, g -> g' becomes g' -> g, each transfer keeping
#     its mu histogram (or Legendre moments);
#   - absorption' = max(0, total - in-scatter), in-scatter_g =
#     sum_g'' nuSs(g'' -> g);
#   - multiplicity c_g = in-scatter_g / (total_g - absorption'_g), constant
#     along the row (multiplicity_matrix rewritten, scatter_matrix holds
#     nu-scatter);
#   - total unchanged;
#   - all FISSION_DATASETS removed and fissionable = False.
# The 2026-11-20 refresh split transpose_library out of make_adjoint_library
# with unchanged arithmetic (ce_mg_triga CHANGELOG_transpose.md); both sync
# fixtures reproduce exactly.
# Note: the transposed absorption is a balance quantity, not a reaction rate;
# it can be zeroed where in-scatter exceeds total.
"""Adjoint (transposed) MGXS library, copied verbatim from ce_mg_triga."""
import os
from pathlib import Path

import h5py
import numpy as np


# ---------------------------------------------------------------------------
FISSION_DATASETS = ("nu-fission", "chi", "fission", "kappa-fission",
                    "prompt-nu-fission", "delayed-nu-fission", "chi-prompt", "chi-delayed",
                    "beta", "decay-rate")


def remove_fission(src, dst):
    """Copy an MGXS library with fission turned off: no fission neutrons are emitted.

    Absorption and total keep their fission contribution, so neutrons absorbed by
    fission still disappear. This is the transport operator without fission used
    for the source-to-detector matrix G (brief, task 3) and for fixed-source
    FW-CADIS weight windows. Returns the number of materials changed.
    """
    import os
    import shutil
    import h5py

    # Build a private temporary copy and rename it into place, so concurrent jobs
    # never read (or write) a half-finished file.
    dst = Path(dst)
    tmp = dst.with_name(f".{dst.name}.tmp.{os.getpid()}")
    shutil.copyfile(src, tmp)
    changed = 0
    with h5py.File(tmp, "r+") as f:
        for group in f.values():
            if not isinstance(group, h5py.Group) or not group.attrs.get("fissionable", False):
                continue
            group.attrs["fissionable"] = False
            for key, data in group.items():
                if isinstance(data, h5py.Group):
                    for name in FISSION_DATASETS:
                        if name in data:
                            del data[name]
            changed += 1
    os.replace(tmp, dst)
    return changed


def _nbins(group):
    fmt = group.attrs["scatter_format"].decode() if isinstance(
        group.attrs["scatter_format"], bytes) else group.attrs["scatter_format"]
    order = int(group.attrs["order"])
    if fmt == "legendre":
        return order + 1, fmt
    if fmt == "histogram":
        return order, fmt
    raise ValueError(f"scatter_format {fmt} not supported")


def _dense(sd, G, nb):
    gmin, gmax = sd["g_min"][()], sd["g_max"][()]
    flat = sd["scatter_matrix"][()].reshape(-1, nb)
    mult = sd["multiplicity_matrix"][()] if "multiplicity_matrix" in sd else np.ones(len(flat))
    N = np.zeros((G, G, nb))
    M = np.ones((G, G))
    i = 0
    for g, (a, b) in enumerate(zip(gmin, gmax)):
        k = b - a + 1
        N[g, a - 1:b] = flat[i:i + k]
        M[g, a - 1:b] = mult[i:i + k]
        i += k
    return N, M


def _write_sparse(sd, N, M):
    G, _, nb = N.shape
    gmin, gmax, flat, mult = [], [], [], []
    for g in range(G):
        nz = np.flatnonzero(np.abs(N[g]).sum(1) > 0)
        a, b = (nz[0], nz[-1]) if len(nz) else (g, g)
        gmin.append(a + 1)
        gmax.append(b + 1)
        flat.append(N[g, a:b + 1].ravel())
        mult.append(M[g, a:b + 1])
    for key in ("g_min", "g_max", "scatter_matrix", "multiplicity_matrix"):
        if key in sd:
            del sd[key]
    sd.create_dataset("g_min", data=np.array(gmin, dtype=np.int64))
    sd.create_dataset("g_max", data=np.array(gmax, dtype=np.int64))
    sd.create_dataset("scatter_matrix", data=np.concatenate(flat))
    sd.create_dataset("multiplicity_matrix", data=np.concatenate(mult))


def transpose_library(N, total, absorption=None, multiplicity=None, scatter_format="histogram"):
    """Adjoint cross sections of one material at one temperature. Pure: arrays in, arrays out.

    Copied verbatim by the noise project (mg-complex-weight-noise-solver); keep it free of
    paths and I/O, and log every change in CHANGELOG_transpose.md (brief v20).

    N[g, g', b]   forward nu-scatter g -> g', groups highest energy first, last axis the
                  angular bins (histogram bins of mu, or Legendre moments P0..PL)
    total[g]      total cross section (kept)
    absorption    forward absorption (only for the balance diagnostic)
    multiplicity  forward multiplicity matrix (only for the balance diagnostic; default 1)

    The adjoint transfer g -> g' is the forward g' -> g, with the same angular bins (the
    kernel depends only on mu). OpenMC MG scatters with probability 1 - Sa/St and
    multiplies the weight by the multiplicity, so with in-scatter_g = sum_g'' nuSs(g'' -> g):
        Sa'_g = max(0, St_g - in-scatter_g),  c_g = in-scatter_g / (St_g - Sa'_g)  (1 if 0/0)
    gives exactly in-scatter_g / St_g secondaries per collision; c_g >= 1 and is constant
    across row g of the adjoint multiplicity matrix. The stored matrix stays nu-scatter.

    Returns (N_adj[g, g', b], absorption_adj[g], multiplicity_adj[g, g'], info) with info =
    {"balance_error": max |(Sa + sum nuSs/M)/St - 1| of the forward data (nan if absorption
    is None), "groups_absorption_zeroed": int, "max_multiplicity": float}.
    """
    N = np.asarray(N, float)
    total = np.asarray(total, float)
    G = len(total)
    if N.shape[:2] != (G, G):
        raise ValueError(f"scatter matrix {N.shape} does not match {G} groups")
    N0 = N[..., 0] if scatter_format == "legendre" else N.sum(2)
    info = {"balance_error": float("nan")}
    if absorption is not None:
        M = np.ones((G, G)) if multiplicity is None else np.asarray(multiplicity, float)
        scat = np.where(M > 0, N0 / np.where(M > 0, M, 1), 0).sum(1)
        ok = (total > 0) & (scat > 0)
        info["balance_error"] = (float(np.abs((np.asarray(absorption) + scat)[ok] / total[ok] - 1).max())
                                 if ok.any() else 0.0)
    NA = N.transpose(1, 0, 2).copy()
    inscatter = N0.sum(0)                                 # sum_g'' nuSs(g'' -> g)
    abs_adj = np.maximum(0.0, total - inscatter)
    branch = total - abs_adj
    c = np.where(branch > 0, inscatter / np.where(branch > 0, branch, 1), 1.0)
    info["groups_absorption_zeroed"] = int(np.sum((inscatter > total) & (total > 0)))
    info["max_multiplicity"] = float(c.max())
    MA = np.repeat(c[:, None], G, axis=1)
    return NA, abs_adj, MA, info


def make_adjoint_library(src, dst):
    """Adjoint, fission-free copy of a forward MGXS library (see module docstring); the
    file I/O around transpose_library."""
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.name}.tmp.{os.getpid()}")
    remove_fission(src, tmp)
    report = {"materials": 0, "max_balance_error": 0.0, "max_multiplicity": 1.0,
              "groups_absorption_zeroed": 0}
    with h5py.File(tmp, "r+") as f:
        G = int(f.attrs["energy_groups"])
        for name, group in f.items():
            if not isinstance(group, h5py.Group):
                continue
            nb, fmt = _nbins(group)
            for key, data in group.items():
                if key == "kTs" or not isinstance(data, h5py.Group):
                    continue
                sd = data["scatter_data"]
                N, M = _dense(sd, G, nb)
                NA, abs_adj, MA, info = transpose_library(N, data["total"][()], data["absorption"][()], M, fmt)
                if info["balance_error"] > report["max_balance_error"]:
                    report["max_balance_error"] = info["balance_error"]
                    report["worst_balance_material"] = name
                report["groups_absorption_zeroed"] += info["groups_absorption_zeroed"]
                report["max_multiplicity"] = max(report["max_multiplicity"], info["max_multiplicity"])
                _write_sparse(sd, NA, MA)
                del data["absorption"]
                data.create_dataset("absorption", data=abs_adj)
            report["materials"] += 1
    os.replace(tmp, dst)
    # absorption + sum(nu-scatter / multiplicity) should equal total. Legendre
    # libraries balance to round-off; the 32-bin histogram library is off by ~5e-4
    # (its mu-binned tally is not exactly the P0 scatter tally). A wrong reading of
    # the bins (densities instead of bin integrals) would be off by O(1), so only
    # that is fatal.
    if report["max_balance_error"] > 1e-4:
        print(f"note: forward balance error {report['max_balance_error']:.2e} "
              f"({report.get('worst_balance_material')}); absorption' = total - in-scatter absorbs it")
    if report["max_balance_error"] > 1e-2:
        tmp.unlink(missing_ok=True)
        dst.unlink(missing_ok=True)
        raise RuntimeError(f"{src}: forward balance error {report['max_balance_error']:.2e}; "
                           "scatter-matrix convention not understood, adjoint library not written")
    return report

