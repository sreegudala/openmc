"""Stage 3 gates: fission with 1/k, delayed neutrons, fission-generation loop.

Infinite medium (20 cm reflective cube), uniform real unit source.

3a  Sigma_t=1, Sigma_s=0, Sigma_f=0.4, nu=2.5 (k=1), beta=7e-3, lambda=0.08,
    v=3e5; F7 at eta=1 (reference G0 S/Sigma_a), plus 1e4 Hz with source
    pi eps_f Upsilon Sigma_a and eta=20 (reference 2.49e-3 - 7.44e-2 i).
3b  TRIGA-like: Sigma_t=1, Sigma_s=0.9, nuSigma_f=Sigma_a=0.1,
    v=2.1097e5 (Lambda=47.4 us), six families (IFP beta_k, Keepin lambda_k);
    F7 with eta=1 up to 30 Hz and eta=3 at 100 Hz.
3c  3b at 100 Hz with eta=1: measured generation ratio of sum|w| within 1 %
    of K~cb = Upsilon~/(1 - kappa w Lambda) (> 1) and the growth guard fires
    (sum|w| > growth_limit x its first-generation value; the doc's
    "50 consecutive increases" rule fires by chance in near-critical
    convergent games, see stage3 log).
3d  CORTEX fuel, two groups, keff = k_inf, F7, source in group 2.

Pass: STD for 3a, 3b, 3d; measured generation ratio within 1 % of K~cb in
convergent cases; generations to extinction of order 1/beta.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "reference"))

from build import box_source, make_library, slab_model  # noqa
from harness import fmt_info, fom, run_replicas, std_criterion, summarize  # noqa
from noise_common import Material  # noqa
from ref_g0 import (KEEPIN_U235_LAMBDA, TRIGA_BETA_K, g0, ktilde_cb,  # noqa
                    upsilon)
from ref_inf import cortex_fuel, ferrero_case, solve as inf_solve, static_k  # noqa

WORK = HERE / "work"
LOG = []
F7 = [0.1, 0.3, 1, 3, 10, 30, 100]


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def triga_1g():
    return Material.from_beta([1.0], [[0.9]], [1 / 2.1097e5], [0.1], [1.0],
                              beta=TRIGA_BETA_K, lam=KEEPIN_U235_LAMBDA)


def cube_model(G, f, eta, keff, particles, batches, lib, src_group=0,
               strength=1.0, phase=0.0, max_gen=100000, noise_extra=None):
    model, _ = slab_model([(20.0, "m")], bc=("reflective", "reflective"),
                          x0=-10.0, transverse=20.0, library=lib)
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.noise = {"frequency": f, "eta": eta, "keff": keff,
                            "max_generations": max_gen, **(noise_extra or {})}
    model.settings.source = [box_source(-10, 10, src_group, G,
                                        transverse=20.0, strength=strength,
                                        phase=phase)]
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


def read_noise(sp_path):
    with h5py.File(sp_path, "r") as f:
        g = f["noise"]
        return {k: g[k][()] for k in g}


def gen_ratio(S, nb, start=5, min_batches=None):
    """Fit log(S_g) vs g where all batches still reach g (or at least half)."""
    S = np.asarray(S)
    nb = np.asarray(nb)
    nmax = nb[0]
    # use generations that every batch still reaches, so the sum over
    # batches is an unbiased estimate of E[sum|w|_g] without truncation
    ok = nb >= nmax if min_batches is None else nb >= min_batches
    idx = np.arange(S.size)
    sel = ok & (idx >= start)
    if sel.sum() < 5:
        return np.nan, 0
    i = idx[sel]
    coef = np.polyfit(i, np.log(S[sel]), 1)
    return np.exp(coef[0]), sel.sum()


def run_case(label, mat, lib, f, eta, keff, ref, args, results, src_group=0,
             strength=1.0, phase=0.0, kt=None):
    extra = lambda sp, d: read_noise(d / f"statepoint.{args.batches}.h5")
    out = run_replicas(lambda: cube_model(mat.G, f, eta, keff, args.particles,
                                          args.batches, lib, src_group,
                                          strength, phase),
                       WORK / label, n_rep=args.reps, tally_names=("tl", "col"),
                       extra=extra)
    entry = {"label": label, "f": f, "eta": eta}
    ok_all = True
    for est in ("tl", "col"):
        vals = out["values"][est].reshape(args.reps, -1)
        mean, sre, sim, corr = summarize(vals)
        ok, info = std_criterion(mean, sre, sim, ref)
        fre, fim = fom(vals, out["times"])
        ok_all &= ok
        entry[est] = {"mean": mean, "s_re": sre, "s_im": sim, "z": info["z"],
                      "ok": ok}
        log(f"  {label} {est}: {'PASS' if ok else 'FAIL'} {fmt_info(info)}; "
            f"mean={np.round(mean, 6)} ref={np.round(ref, 6)} "
            f"corr={np.round(corr, 2)} FOM Re={np.round(fre, 2)} "
            f"Im={np.round(fim, 2)}")
    # generation statistics pooled over replicas
    S = None
    for ex in out["extras"]:
        s = ex["gen_sum_abs_wgt"]
        n = ex["gen_n_batches"]
        if S is None:
            S, NB, NS = s.copy(), n.copy(), ex["gen_n_sites"].copy()
        else:
            m = max(S.size, s.size)
            S = np.pad(S, (0, m - S.size)) + np.pad(s, (0, m - s.size))
            NB = np.pad(NB, (0, m - NB.size)) + np.pad(n, (0, m - n.size))
            NS = np.pad(NS, (0, m - NS.size)) + np.pad(
                ex["gen_n_sites"], (0, m - ex["gen_n_sites"].size))
    gens = np.concatenate([ex["generations_per_batch"] for ex in out["extras"]])
    ratio, npts = gen_ratio(S, NB)
    entry.update({"gens_mean": gens.mean(), "gens_max": int(gens.max()),
                  "ratio": ratio, "ktilde": kt, "time": out["times"].sum()})
    rmsg = ""
    if kt is not None and np.isfinite(ratio):
        rel = abs(ratio / kt - 1)
        entry["ratio_ok"] = rel < 0.01
        rmsg = (f"ratio={ratio:.5f} (fit over {npts} gens) K~cb={kt:.5f} "
                f"rel diff {rel:.2e} {'ok' if rel < 0.01 else 'OFF'}")
    log(f"  {label} generations/batch mean={gens.mean():.0f} max={gens.max()} "
        f"{rmsg} T={out['times'].sum():.0f}s")
    entry["ok"] = ok_all
    entry["S"] = S
    entry["NB"] = NB
    entry["NS"] = NS
    results.append(entry)
    return entry


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=5000)
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--only", default="3a,3b,3c,3d")
    args = p.parse_args()
    only = args.only.split(",")
    WORK.mkdir(exist_ok=True)
    results, gate = [], []

    if "3a" in only:
        log("Test 3a: published point data (one family)")
        mat = ferrero_case()
        lib = WORK / "lib_3a.h5"
        make_library({"m": mat}, lib)
        k, _ = static_k(mat)
        Lam = 1 / (3e5 * 1.0)
        for f in F7:
            w = 2 * np.pi * f
            ref = inf_solve(mat, w, [1.0], k)
            kt = ktilde_cb(w, 1.0, [7e-3], [0.08], Lam)
            e = run_case(f"3a_{f}Hz", mat, lib, f, 1.0, k, ref, args, results,
                         kt=kt)
            gate.append((e["label"], e["ok"]))
            gate.append((e["label"] + " gen ratio", e.get("ratio_ok", False)))
        # 1e4 Hz point, eta = 20
        f = 1e4
        w = 2 * np.pi * f
        amp = np.pi * 5e-3 * upsilon(w, [7e-3], [0.08]) * 1.0
        ref = inf_solve(mat, w, [amp], k)
        log(f"  1e4 Hz reference {ref[0]:.4e} (published 2.49e-3 - 7.44e-2 i)")
        e = run_case("3a_1e4Hz", mat, lib, f, 20.0, k, ref, args, results,
                     strength=abs(amp), phase=np.angle(amp),
                     kt=ktilde_cb(w, 20.0, [7e-3], [0.08], Lam))
        gate.append((e["label"], e["ok"]))

    if "3b" in only:
        log("Test 3b: TRIGA-like six-family kinetics")
        mat = triga_1g()
        lib = WORK / "lib_3b.h5"
        make_library({"m": mat}, lib)
        Lam = 1 / (2.1097e5 * 0.1)
        for f in F7:
            eta = 3.0 if f >= 100 else 1.0
            w = 2 * np.pi * f
            ref = inf_solve(mat, w, [1.0], 1.0)
            g0ref = g0(w, TRIGA_BETA_K, KEEPIN_U235_LAMBDA, Lam) / 0.1
            assert abs(ref[0] / g0ref - 1) < 1e-10
            kt = ktilde_cb(w, eta, TRIGA_BETA_K, KEEPIN_U235_LAMBDA, Lam)
            e = run_case(f"3b_{f}Hz_eta{eta}", mat, lib, f, eta, 1.0, ref,
                         args, results, kt=kt)
            gate.append((e["label"], e["ok"]))
            gate.append((e["label"] + " gen ratio", e.get("ratio_ok", False)))

    if "3c" in only:
        log("Test 3c: predicted failure, TRIGA-like at 100 Hz with eta=1")
        mat = triga_1g()
        lib = WORK / "lib_3b.h5"
        make_library({"m": mat}, lib)
        Lam = 1 / (2.1097e5 * 0.1)
        w = 2 * np.pi * 100
        kt = ktilde_cb(w, 1.0, TRIGA_BETA_K, KEEPIN_U235_LAMBDA, Lam)
        d = WORK / "3c"
        d.mkdir(exist_ok=True)
        model = cube_model(1, 100.0, 1.0, 1.0, args.particles * 4, 1, lib)
        model.settings.verbosity = 7
        t0 = time.time()
        import re
        import subprocess
        model.export_to_model_xml(d / "model.xml")
        res = subprocess.run(["openmc", "-s", "12"], cwd=d, capture_output=True,
                             text=True)
        text = res.stdout + res.stderr
        fired = res.returncode != 0 and "grew by a factor" in text
        msg = [l for l in text.splitlines() if "ERROR" in l or "grew by" in l]
        msg = msg[-1] if msg else f"return code {res.returncode}"
        S = []
        for line in res.stdout.splitlines():
            m = re.search(r"Noise fission generation\s+(\d+)\s+sites:\s+(\d+)"
                          r"\s+sum\|w\|:\s+(\S+)", line)
            if m:
                S.append(float(m.group(3)))
        S = np.array(S)
        ratio, npts = gen_ratio(S, np.full(S.size, 2), start=20)
        rel = abs(ratio / kt - 1) if np.isfinite(ratio) else np.nan
        log(f"  3c guard {'fired' if fired else 'DID NOT FIRE'}: {msg[:200]}")
        log(f"  3c measured ratio {ratio:.5f} over {npts} gens vs K~cb "
            f"{kt:.5f}: rel diff {rel:.2e}; T={time.time()-t0:.0f}s")
        np.save(HERE / "stage3c_sumw.npy", S)
        gate.append(("3c guard fires", fired))
        gate.append(("3c ratio > 1 and within 1 % of K~cb",
                     bool(ratio > 1 and rel < 0.01)))

    if "3d" in only:
        log("Test 3d: two-group CORTEX fuel, keff = k_inf")
        mat = cortex_fuel()
        lib = WORK / "lib_3d.h5"
        make_library({"m": mat}, lib)
        k, _ = static_k(mat)
        log(f"  k_inf = {k:.10f}")
        for f in F7:
            w = 2 * np.pi * f
            ref = inf_solve(mat, w, [0.0, 1.0], k)[::-1]
            e = run_case(f"3d_{f}Hz", mat, lib, f, 1.0, k, ref, args, results,
                         src_group=1)
            gate.append((e["label"], e["ok"]))

    zall = np.concatenate([np.concatenate([r["tl"]["z"], r["col"]["z"]])
                           for r in results])
    if zall.size:
        log(f"Pooled z: n={zall.size} mean={zall.mean():+.2f} "
            f"std={zall.std(ddof=1):.2f} |z|>=3: {int(np.sum(np.abs(zall) >= 3))}")
    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 3 ({args.only}): {len(gate) - n_fail}/{len(gate)} gate checks "
        "passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_")
    with open(HERE / f"stage3_log_{tag}.md", "w") as fh:
        fh.write("# Stage 3 log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    np.savez(HERE / f"stage3_gen_{tag}.npz",
             **{r["label"] + "_S": r["S"] for r in results},
             **{r["label"] + "_NB": r["NB"] for r in results})
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
