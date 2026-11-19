"""Stage 6 gates: noise source for an absorber of varying strength.

The noise source -dSigma psi_0 is sampled at the collisions of one static
generation per batch (after nskip - 1 plain static generations), and the
static flux phi_0 is tallied in the same generation through a
NoisePhaseFilter, so dphi / phi_0 needs no external normalisation.

6a  case 3b (TRIGA-like infinite medium) with a uniform
    dSigma_a / Sigma_a = 1e-3 (fission unchanged). Reference:
    dphi / phi_0 = -G_0 dSigma_a / Sigma_a, computed exactly with the
    infinite-medium noise operator (ref_inf.solve), F7.
6b  slab 4a with dSigma_a / Sigma_a = 0.01 in x in [14, 15] cm. Reference:
    ref_sn1d with source -dSigma_a psi_0 from its own static angular flux;
    dphi / phi_0 on 60 1-cm bins, 0.1, 1, 10, 100 Hz. nskip sensitivity
    (3 against 20) at 1 Hz is logged.
(6c, CORTEX against Abeille, is in run_stage6c.py.)

Pass: 6a STD; 6b STD-C (design doc amendment 2026-11-18); the sampled
source integrates to -dSigma_a x (static collision rate in the region) /
Sigma_t = -dSigma_t phi_0 (region) within statistics (|z| < 3).
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "reference"))
sys.path.insert(0, str(HERE.parent / "stage3"))
sys.path.insert(0, str(HERE.parent / "stage4"))

from build import box_source, make_library, slab_model, x_mesh  # noqa
from harness import (fmt_info, fmt_info_correlated, last_statepoint,  # noqa
                     run_replicas, std_correlated, std_criterion, summarize)
import run_stage3 as s3  # noqa
import run_stage4 as s4  # noqa
from ref_g0 import KEEPIN_U235_LAMBDA, TRIGA_BETA_K, g0  # noqa
from ref_inf import solve as inf_solve  # noqa

WORK = HERE / "work"
LOG = []
F7 = [0.1, 0.3, 1, 3, 10, 30, 100]


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def read_source(sp, d):
    """Per-replica batch mean of the sampled source (per static particle)."""
    with h5py.File(last_statepoint(d), "r") as f:
        g = f["noise"]
        n = f["n_particles"][()]
        ss = g["source_sum"][()].reshape(-1, g["n_perturbations"][()], 2)
        return {"source": (ss[..., 0] + 1j * ss[..., 1]).mean(0) / n,
                "static_k": g["static_k"][()],
                "n_src": g["source_n_sites"][()].mean(),
                "gens": g["generations_per_batch"][()].mean()}


def source_check(label, src, expected):
    """z of the sampled-source integral against -dSigma_t phi_0 (region),
    replica by replica (both from the same static generation)."""
    d = src - expected
    R = d.size
    out = []
    for part, x in (("Re", d.real), ("Im", d.imag)):
        s = x.std(ddof=1) / np.sqrt(R)
        out.append((part, x.mean() / s if s > 0 else 0.0))
    ok = all(abs(z) < 3 for _, z in out)
    log(f"  {label} source integral: mean {src.mean():.6e} vs -dSigma_t "
        f"phi_0 {expected.mean():.6e}; z " +
        ", ".join(f"{p} {z:+.2f}" for p, z in out) +
        f" -> {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles-6a", type=int, default=2000,
                   help="static particles per generation for 6a (about 18 "
                   "noise source sites each)")
    p.add_argument("--particles-6b", type=int, default=60000,
                   help="static particles per generation for 6b (about 0.37 "
                   "noise source sites each)")
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--nskip", type=int, default=20)
    p.add_argument("--static-inactive", type=int, default=200,
                   help="static generations before the first batch (slab 4a "
                   "has dominance ratio ~0.975)")
    p.add_argument("--only", default="6a,6b,nskip")
    p.add_argument("--freqs", default=None,
                   help="comma-separated subset of the 6a/6b frequencies")
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--tag", default="")
    args = p.parse_args()
    only = args.only.split(",")
    sfx = f"_{args.tag}" if args.tag else ""
    global WORK
    if args.tag:
        WORK = HERE / f"work_{args.tag}"
    WORK.mkdir(exist_ok=True)
    sel = None if args.freqs is None else [float(x) for x in
                                           args.freqs.split(",")]
    gate = []
    saved = {}

    # ---------------------------------------------------------------- 6a
    if "6a" in only:
        log("Test 6a: infinite medium (case 3b), uniform dSigma_a/Sigma_a "
            "= 1e-3")
        mat = s3.triga_1g()
        lib = WORK / "lib_3b.h5"
        make_library({"m": mat}, lib)
        sig_t = mat.total[0]
        sig_a = sig_t - mat.scatter[0].sum()
        eps = 1e-3
        d_total = eps * sig_a / sig_t
        Lam = 1 / (2.1097e5 * 0.1)
        for f in [x for x in F7 if sel is None or x in sel]:
            eta = 3.0 if f >= 100 else 1.0
            w = 2 * np.pi * f
            ref = inf_solve(mat, w, [-d_total * sig_t], 1.0)[0]
            ref_g0 = -g0(w, TRIGA_BETA_K, KEEPIN_U235_LAMBDA, Lam) * eps

            def make():
                model, _ = slab_model([(20.0, "m")],
                                      bc=("reflective", "reflective"),
                                      x0=-10.0, transverse=20.0, library=lib)
                cell = list(model.geometry.root_universe.cells.values())[0]
                model.settings.particles = args.particles_6a
                model.settings.batches = args.batches
                model.settings.noise = {
                    "frequency": f, "eta": eta, "keff": 1.0,
                    "nskip": args.nskip,
                    "static_inactive": args.static_inactive,
                    "perturbations": [{"cells": [cell], "total": d_total}]}
                model.settings.source = [box_source(-10, 10, 0, 1,
                                                    transverse=20.0)]
                t = openmc.Tally(name="tl")
                t.filters = [openmc.NoisePhaseFilter(),
                             openmc.ComplexPartFilter()]
                t.scores = ["flux"]
                model.tallies = openmc.Tallies([t])
                return model
            out = run_replicas(make, WORK / f"6a_{f}", n_rep=args.reps,
                               base_seed=args.base_seed, extra=read_source)
            v = out["values"]["tl"].reshape(args.reps, 2)
            phi0 = v[:, 0].real
            ratio = v[:, 1] / phi0
            mean, sre, sim, corr = summarize(ratio[:, None])
            ok, info = std_criterion(mean, sre, sim, np.array([ref]))
            log(f"  6a {f} Hz (eta {eta}): {'PASS' if ok else 'FAIL'} "
                f"{fmt_info(info)}; dphi/phi0 = {mean[0]:.6f} ref "
                f"{ref:.6f} (-G0 eps {ref_g0:.6f})")
            src = np.array([e["source"][0] for e in out["extras"]])
            ok_s = source_check(f"6a {f} Hz", src, -d_total * sig_t * phi0)
            ex = out["extras"]
            sk = np.concatenate([e["static_k"][args.static_inactive:]
                                 for e in ex])
            log(f"    static k {sk.mean():.5f} +- "
                f"{sk.std(ddof=1) / np.sqrt(sk.size):.5f} (per generation, "
                f"correlated); source sites/batch "
                f"{np.mean([e['n_src'] for e in ex]):.0f}; generations/batch "
                f"{np.mean([e['gens'] for e in ex]):.0f}; T="
                f"{out['times'].sum():.0f}s")
            gate.append((f"6a {f} Hz", ok))
            gate.append((f"6a {f} Hz source integral", ok_s))
            saved[f"6a_{f}_vals"] = v
            saved[f"6a_{f}_src"] = src
            saved[f"6a_{f}_times"] = out["times"]

    # ---------------------------------------------------------------- 6b
    def slab_case(label, f, nskip, freq_ref):
        mat = s4.mat_4a()
        lib = WORK / "lib_4a.h5"
        if not lib.exists():
            make_library({"m": mat}, lib)
        eta = 3.0 if f >= 100 else 1.0
        sl, k_sn, ref_cells, d_total, phi0_cells = freq_ref

        def make():
            model, edges = slab_model([(44.0, "m"), (1.0, "m"), (15.0, "m")],
                                      x0=-30.0, library=lib)
            cell = list(model.geometry.root_universe.cells.values())[1]
            model.settings.particles = args.particles_6b
            model.settings.batches = args.batches
            model.settings.noise = {
                "frequency": f, "eta": eta, "keff": k_sn, "nskip": nskip,
                "static_inactive": args.static_inactive,
                "perturbations": [{"cells": [cell], "total": d_total}]}
            model.settings.source = [box_source(-30, 30, 0, 1)]
            t = openmc.Tally(name="tl")
            t.filters = [openmc.MeshFilter(x_mesh(-30, 30, 60)),
                         openmc.NoisePhaseFilter(),
                         openmc.ComplexPartFilter()]
            t.scores = ["flux"]
            model.tallies = openmc.Tallies([t])
            return model
        out = run_replicas(make, WORK / f"{label}", n_rep=args.reps,
                           base_seed=args.base_seed, extra=read_source)
        v = out["values"]["tl"].reshape(args.reps, 60, 2)
        phi0 = v[:, :, 0].real
        ratio = v[:, :, 1] / phi0
        edges = np.arange(-30, 31, 1.0)
        ref = (s4.bin_to_mesh(sl, ref_cells, edges)[:, 0]
               / s4.bin_to_mesh(sl, phi0_cells, edges)[:, 0])
        mask = np.abs(ref) > 0.1 * np.abs(ref).max()
        ok, info = std_correlated(ratio, ref, mask)
        log(f"  {label} {f} Hz (eta {eta}, nskip {nskip}): "
            f"{'PASS' if ok else 'FAIL'} (STD-C) {fmt_info_correlated(info)}")
        m = ratio.mean(0)
        i14 = 44  # bin [14, 15]
        log(f"    dphi/phi0 at [14,15]: {m[i14]:.5e} ref {ref[i14]:.5e}; "
            f"at [0,1]: {m[30]:.5e} ref {ref[30]:.5e}")
        src = np.array([e["source"][0] for e in out["extras"]])
        ok_s = source_check(f"{label} {f} Hz", src,
                            -d_total * phi0[:, i14])  # Sigma_t = 1
        ex = out["extras"]
        sk = np.concatenate([e["static_k"][args.static_inactive:] for e in ex])
        log(f"    static k {sk.mean():.5f} (S_N {k_sn:.5f}); source "
            f"sites/batch {np.mean([e['n_src'] for e in ex]):.0f}; "
            f"generations/batch {np.mean([e['gens'] for e in ex]):.0f}; T="
            f"{out['times'].sum():.0f}s")
        saved[f"{label}_{f}_vals"] = v
        saved[f"{label}_{f}_src"] = src
        saved[f"{label}_{f}_times"] = out["times"]
        saved[f"{label}_{f}_ref"] = ref
        return ok, ok_s

    if "6b" in only or "nskip" in only:
        mat = s4.mat_4a()
        sl = s4.Slab([(60.0, mat, 480)], x0=-30.0, n_angles=64)
        k_sn, _ = sl.eigenvalue()
        sig_a = mat.total[0] - mat.scatter[0].sum()
        d_total = 0.01 * sig_a / mat.total[0]
        log(f"Test 6b: slab 4a, dSigma_a/Sigma_a = 0.01 in [14, 15] "
            f"(dSigma_t/Sigma_t = {d_total:.4g}); S_N k = {k_sn:.8f}")
        refs = {}
        for f in [0.1, 1, 10, 100]:
            q = sl.perturbation_source(2 * np.pi * f, (14, 15),
                                       d_total=[d_total * mat.total[0]])
            refs[f] = (sl, k_sn, sl.solve(2 * np.pi * f, q), d_total,
                       sl.phi0)
        if "6b" in only:
            for f in [x for x in [0.1, 1, 10, 100] if sel is None or x in sel]:
                ok, ok_s = slab_case("6b", f, args.nskip, refs[f])
                gate.append((f"6b {f} Hz", ok))
                gate.append((f"6b {f} Hz source integral", ok_s))
        if "nskip" in only:
            log("nskip sensitivity (logged, not gated): 6b at 1 Hz, nskip 3")
            slab_case("6b_nskip3", 1, 3, refs[1])

    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 6 ({args.only}): {len(gate) - n_fail}/{len(gate)} gate "
        "checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_") + sfx
    np.savez(HERE / f"stage6_vals_{tag}.npz", **saved)
    with open(HERE / f"stage6_log_{tag}.md", "w") as fh:
        fh.write("# Stage 6 log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
