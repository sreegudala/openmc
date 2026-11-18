"""Replica runner and the design document's standard pass criterion (STD).

STD (design doc Section 5): for each tallied quantity z = (q - q_ref)/sigma for
Re and Im separately, with sigma the *replica* standard deviation of the mean.
Pass if |z| < 3 for at least 99 % of entries and, with >= 30 entries, the z
histogram has mean within +-0.3 and standard deviation within 0.8-1.2. In
addition amplitude must agree within amp_tol (relative) and phase within
phase_tol (degrees), evaluated only where the replica sigma is below one
third of that tolerance.

STD-C (std_correlated) replaces the histogram test for strongly correlated
tallies such as mesh bins with an integral check and a Hotelling T^2 test;
see its docstring.
"""

import os
import shlex
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import openmc
import scipy.stats
import openmc.noise


def _env_int(name, default):
    value = os.environ.get(name)
    return int(value) if value else default


# Execution settings, overridable from the environment (e.g. in a SLURM job):
#   NOISE_THREADS      OpenMP threads per OpenMC process (default 12)
#   NOISE_JOBS         replicas run concurrently (default 1)
#   NOISE_OPENMC_CMD   command used to launch OpenMC (default "openmc")
THREADS = _env_int("NOISE_THREADS", 12)
JOBS = _env_int("NOISE_JOBS", 1)
OPENMC_CMD = shlex.split(os.environ.get("NOISE_OPENMC_CMD", "openmc"))


def run_openmc(cwd, threads=None):
    """Run OpenMC in `cwd` (which holds model.xml); return the wall time.

    Output goes to cwd/openmc.out. Thread binding variables are removed so
    that concurrent processes on one node are not pinned to the same cores.
    """
    cwd = Path(cwd)
    env = dict(os.environ)
    # OPENMC_CHAIN_FILE: not needed in multigroup mode, and every process
    # would parse it (slow, and has failed on a network filesystem).
    for var in ("OMP_PROC_BIND", "OMP_PLACES", "OPENMC_CHAIN_FILE"):
        env.pop(var, None)
    cmd = OPENMC_CMD + ["-s", str(threads or THREADS)]
    t0 = time.time()
    with open(cwd / "openmc.out", "w") as fh:
        res = subprocess.run(cmd, cwd=cwd, stdout=fh,
                             stderr=subprocess.STDOUT, env=env)
    elapsed = time.time() - t0
    if res.returncode != 0:
        tail = (cwd / "openmc.out").read_text().splitlines()[-20:]
        raise RuntimeError(f"OpenMC failed in {cwd} (code {res.returncode}):\n"
                           + "\n".join(tail))
    return elapsed


def last_statepoint(d):
    sps = sorted(Path(d).glob("statepoint.*.h5"),
                 key=lambda q: int(q.name.split(".")[1]))
    if not sps:
        raise FileNotFoundError(f"no statepoint in {d}")
    return sps[-1]


def run_replicas(make_model, workdir, n_rep=20, threads=None, base_seed=1,
                 tally_names=("tl",), keep=False, extra=None, jobs=None):
    """Run `n_rep` independent replicas (different seeds).

    make_model() must return a fresh openmc.Model. Up to `jobs` replicas
    (default NOISE_JOBS) run concurrently with `threads` OpenMP threads each
    (default NOISE_THREADS). Returns dict with complex replica means per
    tally name, shape (R, ...), the batch standard deviations of the last
    replica, the wall time of each replica and the `extra` results.
    """
    workdir = Path(workdir)
    dirs = []
    for r in range(n_rep):
        d = workdir / f"rep{r:02d}"
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
        model = make_model()
        model.settings.seed = base_seed + 7919 * r
        model.export_to_model_xml(d / "model.xml")
        dirs.append(d)
    n_jobs = max(1, min(jobs or JOBS, n_rep))
    with ThreadPoolExecutor(max_workers=n_jobs) as pool:
        times = list(pool.map(lambda d: run_openmc(d, threads), dirs))

    out = {name: [] for name in tally_names}
    batch_std = {}
    extras = []
    for r, d in enumerate(dirs):
        with openmc.StatePoint(last_statepoint(d)) as sp:
            for name in tally_names:
                st = openmc.noise.complex_tally_stats(sp.get_tally(name=name))
                out[name].append(st["value"].squeeze())
                batch_std[name] = (st["re_std"].squeeze(),
                                   st["im_std"].squeeze())
            if extra is not None:
                extras.append(extra(sp, d))
        if not keep and r < n_rep - 1:
            shutil.rmtree(d)
    res = {name: np.array(v) for name, v in out.items()}
    return {"values": res, "batch_std": batch_std, "times": np.array(times),
            "extras": extras}


def summarize(values):
    """Grand mean and replica standard deviation of the mean (Re, Im)."""
    R = values.shape[0]
    mean = values.mean(axis=0)
    s_re = values.real.std(axis=0, ddof=1) / np.sqrt(R)
    s_im = values.imag.std(axis=0, ddof=1) / np.sqrt(R)
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.array([np.corrcoef(values.real[:, i], values.imag[:, i])[0, 1]
                         for i in range(values.shape[1])]) \
            if values.ndim == 2 else np.corrcoef(values.real, values.imag)[0, 1]
    return mean, s_re, s_im, corr


def std_criterion(mean, s_re, s_im, ref, amp_tol=0.02, phase_tol_deg=1.0,
                  mask=None):
    """Apply STD. Returns (ok, info dict)."""
    mean = np.atleast_1d(mean)
    ref = np.atleast_1d(ref).astype(complex)
    s_re = np.atleast_1d(s_re)
    s_im = np.atleast_1d(s_im)
    if mask is None:
        mask = np.ones(mean.shape, bool)
    z_re = (mean.real - ref.real) / s_re
    z_im = (mean.imag - ref.imag) / s_im
    z = np.concatenate([z_re[mask], z_im[mask]])
    z = z[np.isfinite(z)]
    n = z.size
    n_bad = int(np.sum(np.abs(z) >= 3))
    ok_z = n_bad <= np.floor(0.01 * n)
    ok_hist = True
    if n >= 30:
        ok_hist = abs(z.mean()) <= 0.3 and 0.8 <= z.std(ddof=1) <= 1.2
    # amplitude / phase where precise enough
    A = np.abs(mean)
    Aref = np.abs(ref)
    sA = np.sqrt(((mean.real * s_re)**2 + (mean.imag * s_im)**2)) / A
    sphi = np.degrees(np.sqrt((mean.imag * s_re)**2
                              + (mean.real * s_im)**2) / A**2)
    rel_amp = np.abs(A / Aref - 1)
    dphi = np.degrees(np.angle(mean / ref))
    amp_eval = mask & (sA / Aref < amp_tol / 3)
    ph_eval = mask & (sphi < phase_tol_deg / 3)
    ok_amp = bool(np.all(rel_amp[amp_eval] <= amp_tol))
    ok_ph = bool(np.all(np.abs(dphi[ph_eval]) <= phase_tol_deg))
    info = {"n": n, "n_bad": n_bad, "z_mean": z.mean() if n else np.nan,
            "z_std": z.std(ddof=1) if n > 1 else np.nan,
            "max_abs_z": np.abs(z).max() if n else np.nan,
            "amp_eval": int(amp_eval.sum()), "ph_eval": int(ph_eval.sum()),
            "max_rel_amp": rel_amp[amp_eval].max() if amp_eval.any() else np.nan,
            "max_dphi": np.abs(dphi[ph_eval]).max() if ph_eval.any() else np.nan,
            "ok_z": ok_z, "ok_hist": ok_hist, "ok_amp": ok_amp, "ok_ph": ok_ph,
            "z": z}
    return ok_z and ok_hist and ok_amp and ok_ph, info


def hotelling_pca(vals, ref, mask=None, frac=0.99):
    """Hotelling T^2 test of mean == ref for a correlated complex vector.

    Re and Im of the gated entries are stacked and standardised by their
    replica standard deviations; the test uses the k leading principal
    components of the replica covariance that explain `frac` of the
    variance (k <= R - 3). Under the null, (R - k) / (k (R - 1)) T^2 follows
    F(k, R - k).

    Returns
    -------
    k, T2, p_value
    """
    vals = np.asarray(vals)
    ref = np.asarray(ref).astype(complex)
    if mask is None:
        mask = np.ones(ref.shape, bool)
    v = vals[:, mask]
    r = ref[mask]
    X = np.concatenate([v.real, v.imag], axis=1)
    x0 = np.concatenate([r.real, r.imag])
    R = X.shape[0]
    sd = X.std(axis=0, ddof=1)
    keep = sd > 0
    Xs = (X[:, keep] - X[:, keep].mean(axis=0)) / sd[keep]
    ds = (X[:, keep].mean(axis=0) - x0[keep]) / sd[keep]
    _, S, Vt = np.linalg.svd(Xs, full_matrices=False)
    ev = S**2 / (R - 1)
    k = int(np.searchsorted(np.cumsum(ev) / ev.sum(), frac) + 1)
    k = max(1, min(k, R - 3))
    proj = Vt[:k] @ ds
    T2 = float(R * np.sum(proj**2 / ev[:k]))
    F = (R - k) / (k * (R - 1)) * T2
    return k, T2, float(scipy.stats.f.sf(F, k, R - k))


def std_correlated(vals, ref, mask=None, amp_tol=0.02, phase_tol_deg=1.0,
                   p_min=0.0027):
    """STD-C: standard criterion for strongly correlated tallies (mesh bins).

    Gated: (i) |z| < 3 for >= 99 % of raw entries and the amplitude/phase
    limits, as in STD; (ii) |z| < 3 for Re and Im of the integral sum(v)
    over the gated entries; (iii) Hotelling T^2 on the leading principal
    components (hotelling_pca) with p > p_min (0.0027, the two-sided 3 sigma
    level). Logged only: the raw z histogram and the z histogram of the
    per-replica shape v / sum(v). See the design document, "Amendment
    (2026-11-18)".

    Parameters
    ----------
    vals : (R, N) complex replica values
    ref : (N,) complex reference
    mask : (N,) bool, entries to gate (default all)

    Returns
    -------
    ok : bool
    info : dict with "raw" and "shape" (std_criterion infos), "z_int",
        "hotelling" (k, T2, p) and the flags "ok_raw", "ok_int", "ok_hot"
    """
    vals = np.asarray(vals)
    ref = np.asarray(ref).astype(complex)
    if mask is None:
        mask = np.ones(ref.shape, bool)
    R = vals.shape[0]
    mean, s_re, s_im, _ = summarize(vals)
    _, raw = std_criterion(mean, s_re, s_im, ref, amp_tol, phase_tol_deg,
                           mask=mask)
    ok_raw = bool(raw["ok_z"] and raw["ok_amp"] and raw["ok_ph"])

    v = vals[:, mask]
    r = ref[mask]
    T = v.sum(axis=1)
    rT = r.sum()
    z_int = ((T.real.mean() - rT.real) / (T.real.std(ddof=1) / np.sqrt(R)),
             (T.imag.mean() - rT.imag) / (T.imag.std(ddof=1) / np.sqrt(R)))
    ok_int = bool(abs(z_int[0]) < 3 and abs(z_int[1]) < 3)

    hot = hotelling_pca(vals, ref, mask)
    ok_hot = hot[2] > p_min

    ms, ss_re, ss_im, _ = summarize(v / T[:, None])
    _, shape = std_criterion(ms, ss_re, ss_im, r / rT, amp_tol=np.inf,
                             phase_tol_deg=np.inf)
    info = {"raw": raw, "shape": shape, "z_int": z_int, "hotelling": hot,
            "ok_raw": ok_raw, "ok_int": ok_int, "ok_hot": ok_hot,
            "corr_neighbour": float(np.median(np.diag(
                np.corrcoef(v.real.T), 1))) if v.shape[1] > 1 else np.nan}
    return ok_raw and ok_int and ok_hot, info


def fmt_info_correlated(info):
    raw, sh = info["raw"], info["shape"]
    k, T2, pv = info["hotelling"]
    flag = lambda ok: "ok" if ok else "FAIL"
    return (f"bins [{flag(info['ok_raw'])}] n={raw['n']} bad={raw['n_bad']} "
            f"max|z|={raw['max_abs_z']:.2f}; amp {raw['amp_eval']} evaluated, "
            f"max rel err {raw['max_rel_amp']:.2e}; phase {raw['ph_eval']} "
            f"evaluated, max {raw['max_dphi']:.3f} deg | integral "
            f"[{flag(info['ok_int'])}] z Re={info['z_int'][0]:+.2f} "
            f"Im={info['z_int'][1]:+.2f} | Hotelling [{flag(info['ok_hot'])}] "
            f"k={k} T2={T2:.1f} p={pv:.3f} | logged: raw hist "
            f"mean={raw['z_mean']:+.2f} std={raw['z_std']:.2f}, shape hist "
            f"mean={sh['z_mean']:+.2f} std={sh['z_std']:.2f}, neighbour corr "
            f"{info['corr_neighbour']:.2f}")


def fmt_info(info):
    return (f"n={info['n']} bad={info['n_bad']} z_mean={info['z_mean']:+.2f} "
            f"z_std={info['z_std']:.2f} max|z|={info['max_abs_z']:.2f}; "
            f"amp: {info['amp_eval']} evaluated, max rel err "
            f"{info['max_rel_amp']:.2e}; phase: {info['ph_eval']} evaluated, "
            f"max {info['max_dphi']:.3f} deg")


def fom(values, times):
    """FOM = 1/(sigma_rel^2 T) for Re and Im of the grand mean, with T the
    total wall time of all replicas."""
    mean, s_re, s_im, _ = summarize(values)
    T = times.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        return (1 / ((s_re / np.abs(mean.real))**2 * T),
                1 / ((s_im / np.abs(mean.imag))**2 * T))
