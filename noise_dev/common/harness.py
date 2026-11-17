"""Replica runner and the design document's standard pass criterion (STD).

STD (design doc Section 5): for each tallied quantity z = (q - q_ref)/sigma for
Re and Im separately, with sigma the *replica* standard deviation of the mean.
Pass if |z| < 3 for at least 99 % of entries and, with >= 30 entries, the z
histogram has mean within +-0.3 and standard deviation within 0.8-1.2. In
addition amplitude must agree within amp_tol (relative) and phase within
phase_tol (degrees), evaluated only where the replica sigma is below one
third of that tolerance.
"""

import shutil
import time
from pathlib import Path

import numpy as np
import openmc
import openmc.noise


def run_replicas(make_model, workdir, n_rep=20, threads=12, base_seed=1,
                 tally_names=("tl",), keep=False, extra=None):
    """Run `n_rep` independent replicas (different seeds).

    make_model() must return a fresh openmc.Model. Returns dict with
    complex replica means per tally name, shape (R, ...), the batch standard
    deviations of the last replica, wall times and the last statepoint path.
    """
    workdir = Path(workdir)
    out = {name: [] for name in tally_names}
    batch_std = {}
    times = []
    extras = []
    for r in range(n_rep):
        d = workdir / f"rep{r:02d}"
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
        model = make_model()
        model.settings.seed = base_seed + 7919 * r
        t0 = time.time()
        sp_path = model.run(cwd=d, threads=threads, output=False)
        times.append(time.time() - t0)
        with openmc.StatePoint(sp_path) as sp:
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
