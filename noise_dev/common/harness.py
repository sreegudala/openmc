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
                  mask=None, significance=None):
    """Apply STD. Returns (ok, info dict).

    significance : None (default) or k. None: every evaluated entry must lie
        within the amplitude and phase tolerances (the rule up to Stage 9a).
        k (amendment 2026-11-26b, k = 3 from Stage 9b on): an evaluated entry
        fails only if its difference exceeds the tolerance by more than k of
        its own standard deviation, i.e. the tolerance is exceeded at k
        sigma. The evaluation thresholds (sigma < tol / 3) are unchanged.
    """
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
    k = 0.0 if significance is None else float(significance)
    amp_excess = rel_amp - k * sA / Aref - amp_tol
    ph_excess = np.abs(dphi) - k * sphi - phase_tol_deg
    ok_amp = bool(np.all(amp_excess[amp_eval] <= 0))
    ok_ph = bool(np.all(ph_excess[ph_eval] <= 0))
    info = {"n": n, "n_bad": n_bad, "z_mean": z.mean() if n else np.nan,
            "z_std": z.std(ddof=1) if n > 1 else np.nan,
            "max_abs_z": np.abs(z).max() if n else np.nan,
            "amp_eval": int(amp_eval.sum()), "ph_eval": int(ph_eval.sum()),
            "max_rel_amp": rel_amp[amp_eval].max() if amp_eval.any() else np.nan,
            "max_dphi": np.abs(dphi[ph_eval]).max() if ph_eval.any() else np.nan,
            "n_amp_fail": int(np.sum(amp_excess[amp_eval] > 0)),
            "n_ph_fail": int(np.sum(ph_excess[ph_eval] > 0)),
            "significance": significance,
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


def _z(diff, sigma):
    """diff / sigma, with 0 for an exact match without variance (e.g. the
    Im part of a real tally against a real reference)."""
    if sigma == 0 and diff == 0:
        return 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        return diff / sigma


def std_correlated(vals, ref, mask=None, amp_tol=0.02, phase_tol_deg=1.0,
                   p_min=0.0027, significance=None):
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
    significance : passed to std_criterion (None: every evaluated entry
        within tolerance; 3: amendment 2026-11-26b, from Stage 9b on)

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
                           mask=mask, significance=significance)
    ok_raw = bool(raw["ok_z"] and raw["ok_amp"] and raw["ok_ph"])

    v = vals[:, mask]
    r = ref[mask]
    T = v.sum(axis=1)
    rT = r.sum()
    z_int = (_z(T.real.mean() - rT.real, T.real.std(ddof=1) / np.sqrt(R)),
             _z(T.imag.mean() - rT.imag, T.imag.std(ddof=1) / np.sqrt(R)))
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


def paired_correlated(on, off, mask=None, p_min=0.0027):
    """STD-C for comparing two replica sets of the same problem (e.g.
    weight cancellation on vs off) run with the same seeds.

    Replica r of both sets shares its seed, so the pairs are positively
    correlated and the per-replica differences d_r = on_r - off_r are
    i.i.d.; the test is on d (paired), which is valid whether or not the
    pairs are correlated. Gated: paired |t| < 3 for >= 99 % of entries (Re
    and Im), |t| < 3 for the integral of d, and Hotelling T^2 of mean(d) = 0
    on the leading principal components with p > p_min.

    Returns
    -------
    ok : bool
    info : dict with "t" (per-entry paired t, Re then Im), "n_bad",
        "t_int", "hotelling", "ok_bins", "ok_int", "ok_hot", "corr_pair"
        (median Re correlation between on and off replicas)
    """
    on = np.asarray(on)
    off = np.asarray(off)
    if mask is None:
        mask = np.ones(on.shape[1], bool)
    R = on.shape[0]
    d = (on - off)[:, mask]
    se_re = d.real.std(axis=0, ddof=1) / np.sqrt(R)
    se_im = d.imag.std(axis=0, ddof=1) / np.sqrt(R)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.concatenate([d.real.mean(0) / se_re, d.imag.mean(0) / se_im])
    t = t[np.isfinite(t)]
    n_bad = int(np.sum(np.abs(t) >= 3))
    ok_bins = n_bad <= np.floor(0.01 * t.size)
    D = d.sum(axis=1)
    t_int = (D.real.mean() / (D.real.std(ddof=1) / np.sqrt(R)),
             D.imag.mean() / (D.imag.std(ddof=1) / np.sqrt(R)))
    ok_int = bool(abs(t_int[0]) < 3 and abs(t_int[1]) < 3)
    hot = hotelling_pca(d, np.zeros(d.shape[1]))
    ok_hot = hot[2] > p_min
    corr = np.median([np.corrcoef(on[:, i].real, off[:, i].real)[0, 1]
                      for i in np.where(mask)[0]])
    info = {"t": t, "n_bad": n_bad, "t_int": t_int, "hotelling": hot,
            "ok_bins": ok_bins, "ok_int": ok_int, "ok_hot": ok_hot,
            "corr_pair": float(corr)}
    return bool(ok_bins and ok_int and ok_hot), info


def welch_p(xa, xb, rel_sys_b=0.0):
    """Two-sided Welch p-values of mean(xa) = mean(xb) along axis 0 (real
    arrays), with an optional relative systematic on b treated as known.
    Entries that agree exactly with no variance in either set (e.g. the Im
    part of two real tallies) get t = 0 and p = 1."""
    Ra, Rb = xa.shape[0], xb.shape[0]
    va = xa.var(axis=0, ddof=1) / Ra
    vb = xb.var(axis=0, ddof=1) / Rb
    vs = (rel_sys_b * xb.mean(axis=0))**2
    diff = xa.mean(axis=0) - xb.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = diff / np.sqrt(va + vb + vs)
        df = (va + vb)**2 / (va**2 / (Ra - 1) + vb**2 / (Rb - 1))
        p = 2 * scipy.stats.t.sf(np.abs(t), df)
    exact = (va + vb + vs == 0) & (diff == 0)
    return np.where(exact, 0.0, t), np.where(exact, 1.0, p)


def two_sample_correlated(a, b, mask=None, b_rel_sys=0.0, p_min=0.0027):
    """STD-C between two independent Monte Carlo replica sets (e.g. OpenMC
    against an Abeille reference), for strongly correlated tallies.

    Gated: (i) Welch p > p_min (the two-sided 3 sigma level at the Welch
    degrees of freedom) for >= 99 % of entries (Re and Im), where b also
    carries a relative systematic b_rel_sys (e.g. the uncertainty of a
    normalisation applied to every replica of b); (ii) the same for Re and
    Im of the integral over the gated entries; (iii) a two-sample Hotelling
    T^2 of mean_a = mean_b on the leading principal components (99 % of the
    pooled within-set variance, k <= Ra + Rb - 3), with
    (Ra + Rb - k - 1) / (k (Ra + Rb - 2)) T^2 ~ F(k, Ra + Rb - k - 1) and
    pass if p > p_min. b_rel_sys is not included in (iii).

    On fields dominated by a few smooth modes, (i) has a false-alarm rate
    well above p_min because exceedances cluster (about 4 % for 20 + 20
    replicas of a six-mode synthetic field); (ii) and (iii) hold their
    nominal rates.

    Returns
    -------
    ok : bool
    info : dict with "z" (Welch t, Re then Im), "n_bad", "z_int",
        "hotelling", "ok_bins", "ok_int", "ok_hot"
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if mask is None:
        mask = np.ones(a.shape[1], bool)
    Ra, Rb = a.shape[0], b.shape[0]
    va, vb = a[:, mask], b[:, mask]

    tr, pr = welch_p(va.real, vb.real, b_rel_sys)
    ti, pi = welch_p(va.imag, vb.imag, b_rel_sys)
    t = np.concatenate([tr, ti])
    pv = np.concatenate([pr, pi])
    fin = np.isfinite(t)
    z, pv = t[fin], pv[fin]
    n_bad = int(np.sum(pv <= p_min))
    ok_bins = n_bad <= np.floor(0.01 * z.size)

    Ta, Tb = va.sum(1), vb.sum(1)
    zi_r, pi_r = welch_p(Ta.real[:, None], Tb.real[:, None], b_rel_sys)
    zi_i, pi_i = welch_p(Ta.imag[:, None], Tb.imag[:, None], b_rel_sys)
    z_int = (float(zi_r[0]), float(zi_i[0]))
    ok_int = bool(pi_r[0] > p_min and pi_i[0] > p_min)

    Xa = np.concatenate([va.real, va.imag], axis=1)
    Xb = np.concatenate([vb.real, vb.imag], axis=1)
    C = np.vstack([Xa - Xa.mean(0), Xb - Xb.mean(0)])
    n = Ra + Rb
    sd = np.sqrt((C**2).sum(0) / (n - 2))
    keep = sd > 0
    Cs = C[:, keep] / sd[keep]
    d = (Xa.mean(0) - Xb.mean(0))[keep] / sd[keep]
    _, S, Vt = np.linalg.svd(Cs, full_matrices=False)
    ev = S**2 / (n - 2)
    k = int(np.searchsorted(np.cumsum(ev) / ev.sum(), 0.99) + 1)
    k = max(1, min(k, n - 3))
    proj = Vt[:k] @ d
    T2 = float(Ra * Rb / n * np.sum(proj**2 / ev[:k]))
    F = (n - k - 1) / (k * (n - 2)) * T2
    p = float(scipy.stats.f.sf(F, k, n - k - 1))
    ok_hot = p > p_min
    info = {"z": z, "n_bad": n_bad, "z_int": z_int, "hotelling": (k, T2, p),
            "ok_bins": bool(ok_bins), "ok_int": ok_int, "ok_hot": ok_hot}
    return bool(ok_bins and ok_int and ok_hot), info


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
