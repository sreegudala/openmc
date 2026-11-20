# Copied from ce_mg_triga (TRIGA proposal repository):
#   common.py         lines 415-450  FISSION_DATASETS, remove_fission
#   step3a_adjoint.py lines 68-171   _nbins, _dense, _write_sparse,
#                                    make_adjoint_library
# Source commit: e018a837d4859a68e1ec6ab9417f09b667941727 (2026-11-19)
# Copied: 2026-11-20
#
# Copied verbatim. Do not edit here; extend in adjoint_noise_lib.py.
# Changes go into ce_mg_triga first and are then copied over again (refresh
# the fixture in tests/ with it; see test_adjoint_lib.py).
#
# Only the import block below is new. The copied code is the forward-to-
# adjoint transformation of a fission-free MGXS library, as used for the
# TRIGA source-to-detector matrix G:
#   - nu-scatter transposed, g -> g' becomes g' -> g, each transfer keeping
#     its mu histogram (or Legendre moments);
#   - absorption' = max(0, total - in-scatter), in-scatter_g =
#     sum_g'' nuSs(g'' -> g);
#   - multiplicity c_g = in-scatter_g / (total_g - absorption'_g), constant
#     along the row (multiplicity_matrix rewritten, scatter_matrix holds
#     nu-scatter);
#   - total unchanged;
#   - all FISSION_DATASETS removed and fissionable = False.
# Checked against the code on 2026-11-20: matches the brief's description.
# Note: the transposed absorption is a balance quantity, not a reaction rate;
# it can be zeroed where in-scatter exceeds total.
"""Adjoint (transposed) MGXS library, copied verbatim from ce_mg_triga."""
import os
import shutil
from pathlib import Path

import h5py
import numpy as np


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


# ---------------------------------------------------------------- library
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


def make_adjoint_library(src, dst):
    """Adjoint, fission-free copy of a forward MGXS library (see module docstring)."""
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
                N0 = N[..., 0] if fmt == "legendre" else N.sum(2)
                total, absorption = data["total"][()], data["absorption"][()]
                # forward consistency: absorption + scatter (nu-scatter / multiplicity) = total
                scat = np.where(M > 0, N0 / np.where(M > 0, M, 1), 0).sum(1)
                ok = (total > 0) & (scat > 0)
                if ok.any():
                    err = float(np.abs((absorption + scat)[ok] / total[ok] - 1).max())
                    if err > report["max_balance_error"]:
                        report["max_balance_error"] = err
                        report["worst_balance_material"] = name
                NA = N.transpose(1, 0, 2).copy()
                inscatter = N0.sum(0)                         # sum_g'' nuSs(g'' -> g)
                abs_adj = np.maximum(0.0, total - inscatter)
                branch = total - abs_adj
                c = np.where(branch > 0, inscatter / np.where(branch > 0, branch, 1), 1.0)
                report["groups_absorption_zeroed"] += int(np.sum((inscatter > total) & (total > 0)))
                report["max_multiplicity"] = max(report["max_multiplicity"], float(c.max()))
                MA = np.repeat(c[:, None], G, axis=1)
                # stored matrix is nu-scatter; OpenMC divides by multiplicity to sample
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

