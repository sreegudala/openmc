"""Stage 2 gates: the i w / v term (copy reaction), no fission.

Infinite medium = 20 cm cube with reflective faces, uniform isotropic real
unit source. Reference for the volume-integrated flux per source particle:
the infinite-medium solve (ref_inf), = 1/(Sigma_a + i w/v) in one group.

2a  one group, Sigma_t=1, Sigma_s=0.5, v=2.2e5, x = w/(v Sigma_a) in
    {0.01, 0.1, 0.5, 1, 5}, eta in {0.5, 1, 5}; plus F7 at eta=1.
2b  two groups, CORTEX water, source in group 1, x defined with group 2.
2c  collision vs track-length estimators agree within statistics.

Pass: STD with 1 % amplitude and 0.5 deg phase; results at different eta
mutually consistent (|z| < 3).

Cases whose extended within-generation radius c~ (absolute-weight game) is
>= 1 are divergent by Fauvel et al.; one is run to confirm the guard fires.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "reference"))

from build import box_source, make_library, slab_model  # noqa
from harness import fmt_info, fom, run_replicas, std_criterion, summarize  # noqa
from noise_common import Material  # noqa
from ref_inf import cortex_water, solve as inf_solve  # noqa

WORK = HERE / "work"
LOG = []
F7 = [0.1, 0.3, 1, 3, 10, 30, 100]


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def one_group():
    return Material(total=[1.0], scatter=[[0.5]], inv_v=[1 / 2.2e5])


def ctilde(mat, omega, eta):
    """Spectral radius of the absolute-weight collision kernel (infinite
    medium, collision-based scheme)."""
    sw = eta * omega * mat.inv_v
    st = mat.total + sw
    M = mat.scatter.T / st[None, :]
    M = M + np.diag(sw / st * np.sqrt(1 + eta**2) / eta)
    return np.abs(np.linalg.eigvals(M)).max()


def cube_model(G, omega, eta, particles, batches, lib):
    f = omega / (2 * np.pi)
    model, _ = slab_model([(20.0, "m")], bc=("reflective", "reflective"),
                          x0=-10.0, transverse=20.0, library=lib)
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.noise = {"frequency": f, "eta": eta}
    model.settings.source = [box_source(-10, 10, 0, G, transverse=20.0)]
    filters = [openmc.EnergyFilter(np.arange(G + 1.0)),
               openmc.ComplexPartFilter()]
    tallies = []
    for name, est in [("tl", "tracklength"), ("col", "collision")]:
        t = openmc.Tally(name=name)
        t.filters = filters
        t.scores = ["flux"]
        t.estimator = est
        tallies.append(t)
    model.tallies = openmc.Tallies(tallies)
    return model


def reference(mat, omega):
    src = np.zeros(mat.G)
    src[0] = 1.0
    # EnergyFilter bins are in ascending energy: reverse to OpenMC order
    return inf_solve(mat, omega, src, k=1.0)[::-1]


def run_case(label, mat, lib, omega, eta, particles, batches, n_rep, results):
    ct = ctilde(mat, omega, eta)
    x_info = f"f={omega/2/np.pi:.4g} Hz eta={eta} c~={ct:.4f}"
    if ct >= 1.0:
        log(f"  {label}: {x_info} -> predicted divergent (c~ >= 1), not run")
        results.append({"label": label, "eta": eta, "omega": omega,
                        "ctilde": ct, "status": "divergent"})
        return None
    t0 = time.time()
    out = run_replicas(lambda: cube_model(mat.G, omega, eta, particles,
                                          batches, lib),
                       WORK / label, n_rep=n_rep, tally_names=("tl", "col"))
    ref = reference(mat, omega)
    entry = {"label": label, "eta": eta, "omega": omega, "ctilde": ct,
             "status": "run", "time": time.time() - t0}
    for est in ("tl", "col"):
        vals = out["values"][est].reshape(n_rep, -1)
        mean, sre, sim, corr = summarize(vals)
        ok, info = std_criterion(mean, sre, sim, ref, amp_tol=0.01,
                                 phase_tol_deg=0.5)
        fre, fim = fom(vals, out["times"])
        entry[est] = {"mean": mean, "s_re": sre, "s_im": sim, "ok": ok,
                      "z": info["z"], "corr": corr, "fom_re": fre,
                      "fom_im": fim, "vals": vals}
        log(f"  {label} {est}: {x_info} {'PASS' if ok else 'FAIL'} "
            f"{fmt_info(info)}; corr(Re,Im)={np.round(corr, 2)} "
            f"FOM Re={np.round(fre, 1)} Im={np.round(fim, 1)} "
            f"T={out['times'].sum():.0f}s")
    # 2c: estimators agree (paired replicas)
    d = entry["tl"]["vals"] - entry["col"]["vals"]
    dm, dre, dim, _ = summarize(d)
    zc = np.concatenate([dm.real / dre, dm.imag / dim])
    entry["z_est"] = zc
    log(f"  {label} 2c tl-col: max|z|={np.abs(zc).max():.2f}")
    results.append(entry)
    return entry


def guard_check(mat, lib, omega, eta):
    """Run a predicted-divergent case briefly and confirm the guard fires."""
    d = WORK / "guard"
    d.mkdir(parents=True, exist_ok=True)
    model = cube_model(mat.G, omega, eta, 1000, 2, lib)
    model.settings.max_secondaries = 2000
    try:
        model.run(cwd=d, threads=12, output=False)
    except RuntimeError as e:
        fired = "max_secondaries" in str(e)
        return fired, str(e).strip().splitlines()[-1][:160]
    return False, "run completed without the guard firing"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=10000)
    p.add_argument("--batches", type=int, default=10)
    p.add_argument("--only", default="2a,2b,f7")
    args = p.parse_args()
    WORK.mkdir(exist_ok=True)
    only = args.only.split(",")
    xs = [0.01, 0.1, 0.5, 1.0, 5.0]
    etas = [0.5, 1.0, 5.0]
    results = []
    gate = []

    def consistency(entries, name):
        """Mutual consistency across eta at fixed x (tl estimator)."""
        ok = True
        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                a, b = entries[i]["tl"], entries[j]["tl"]
                z = np.concatenate([
                    (a["mean"].real - b["mean"].real) /
                    np.hypot(a["s_re"], b["s_re"]),
                    (a["mean"].imag - b["mean"].imag) /
                    np.hypot(a["s_im"], b["s_im"])])
                good = np.all(np.abs(z) < 3)
                ok &= good
                log(f"    eta {entries[i]['eta']} vs {entries[j]['eta']} "
                    f"({name}): max|z|={np.abs(z).max():.2f} "
                    f"{'ok' if good else 'INCONSISTENT'}")
        return ok

    if "2a" in only:
        log("Test 2a: one group infinite medium")
        mat = one_group()
        lib = WORK / "lib_1g.h5"
        make_library({"m": mat}, lib)
        for x in xs:
            omega = x * 2.2e5 * 0.5
            entries = [e for e in (run_case(f"2a_x{x}_eta{eta}", mat, lib,
                                            omega, eta, args.particles,
                                            args.batches, args.reps, results)
                                   for eta in etas) if e is not None]
            for e in entries:
                gate.append((e["label"], e["tl"]["ok"] and e["col"]["ok"]))
            if len(entries) > 1:
                gate.append((f"2a x={x} eta consistency",
                             consistency(entries, f"x={x}")))
        fired, msg = guard_check(mat, lib, 5.0 * 2.2e5 * 0.5, 1.0)
        log(f"  guard on divergent case x=5 eta=1: "
            f"{'fired' if fired else 'DID NOT FIRE'} ({msg})")
        gate.append(("2a guard fires on divergent case", fired))
    if "f7" in only:
        log("Test 2a-F7: one group, F7 frequencies, eta=1")
        mat = one_group()
        lib = WORK / "lib_1g.h5"
        make_library({"m": mat}, lib)
        for f in F7:
            e = run_case(f"2a_F7_{f}Hz", mat, lib, 2 * np.pi * f, 1.0,
                         args.particles, args.batches, args.reps, results)
            gate.append((e["label"], e["tl"]["ok"] and e["col"]["ok"]))
    if "2b" in only:
        log("Test 2b: two-group infinite medium (CORTEX water)")
        mat = cortex_water()
        lib = WORK / "lib_water.h5"
        make_library({"m": mat}, lib)
        for x in xs:
            omega = x * mat.absorption[1] / mat.inv_v[1]
            entries = [e for e in (run_case(f"2b_x{x}_eta{eta}", mat, lib,
                                            omega, eta, args.particles,
                                            args.batches, args.reps, results)
                                   for eta in etas) if e is not None]
            for e in entries:
                gate.append((e["label"], e["tl"]["ok"] and e["col"]["ok"]))
            if len(entries) > 1:
                gate.append((f"2b x={x} eta consistency",
                             consistency(entries, f"x={x}")))

    # 2c over all run cases: pooled z between estimators
    zc = np.concatenate([r["z_est"] for r in results if r["status"] == "run"])
    n_bad = int(np.sum(np.abs(zc) >= 3))
    ok2c = n_bad <= np.floor(0.01 * zc.size) and (
        zc.size < 30 or (abs(zc.mean()) <= 0.3 and 0.8 <= zc.std(ddof=1) <= 1.2))
    log(f"Test 2c pooled tl vs col: n={zc.size} bad={n_bad} "
        f"z_mean={zc.mean():+.2f} z_std={zc.std(ddof=1):.2f}")
    gate.append(("2c estimators agree (pooled)", ok2c))
    # pooled z histogram of all STD entries
    zall = np.concatenate([np.concatenate([r["tl"]["z"], r["col"]["z"]])
                           for r in results if r["status"] == "run"])
    log(f"Pooled z (all cases, both estimators): n={zall.size} "
        f"mean={zall.mean():+.2f} std={zall.std(ddof=1):.2f} "
        f"|z|>=3: {int(np.sum(np.abs(zall) >= 3))}")

    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 2 ({args.only}): {len(gate) - n_fail}/{len(gate)} gate "
        "checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_")
    with open(HERE / f"stage2_log_{tag}.md", "w") as f:
        f.write("# Stage 2 log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    np.save(HERE / f"stage2_z_{tag}.npy", zall)
    summary = [{k: (v if not isinstance(v, (np.ndarray, dict)) else None)
                for k, v in r.items()} for r in results]
    with open(HERE / f"stage2_cases_{tag}.json", "w") as f:
        json.dump(summary, f, indent=1, default=float)
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
