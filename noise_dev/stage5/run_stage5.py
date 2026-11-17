"""Stage 5 gates: weight cancellation (mesh averaging + bank roulette).

5a  case 3b (TRIGA-like infinite medium) with cancellation on a 4x4x4 mesh:
    STD against the reference; log the generation ratio vs K and K~.
5b  bias: case 4a at 1 Hz, cancellation off vs on with h = 2, 1, 0.5 cm,
    bin-wise t = (on - off)/sqrt(s_on^2 + s_off^2) for Re and Im, bins with
    > 20 % relative error excluded; repeated with the mesh origin shifted h/2.
5c  case 3c (100 Hz, eta=1) now terminates and meets STD; the 1e4 Hz point
    of 3a at eta=1 meets STD.
5d  FOM on vs off for 4a at 0.1 and 1 Hz.
5e  5b repeated with a 32-bin histogram scattering library (linearly
    anisotropic, mu_bar = 0.3) against the same library with cancellation off.

Pass: 5a, 5c STD; 5b, 5e t-histograms standard normal at the two finest
meshes and no trend with origin shift (failure at h = 2 cm recorded as the
bias scale, not a defect); 5d FOM ratio > 1; sum|w| plateaus/decays.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "reference"))
sys.path.insert(0, str(HERE.parent / "stage3"))
sys.path.insert(0, str(HERE.parent / "stage4"))

from build import make_library  # noqa
from harness import fmt_info, fom, run_replicas, std_criterion, summarize  # noqa
import run_stage3 as s3  # noqa
import run_stage4 as s4  # noqa
from ref_g0 import KEEPIN_U235_LAMBDA, TRIGA_BETA_K, ktilde_cb, upsilon  # noqa
from ref_inf import ferrero_case, solve as inf_solve  # noqa

WORK = HERE / "work"
LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def cube_mesh(n):
    m = openmc.RegularMesh()
    m.lower_left = (-10, -10, -10)
    m.upper_right = (10, 10, 10)
    m.dimension = (n, n, n)
    return m


def slab_mesh(h, shift=0.0, a=30.0):
    n = int(round(2 * a / h)) + (1 if shift else 0)
    m = openmc.RegularMesh()
    m.lower_left = (-a - shift, -5, -5)
    m.upper_right = (-a - shift + n * h, 5, 5)
    m.dimension = (n, 1, 1)
    return m


def physical_K(omega, beta_k, lam_k, Lam):
    """Physical generation radius (flight-based expression of Fauvel et al.)
    |Upsilon / (1 + i w Lambda)|."""
    return abs(upsilon(omega, beta_k, lam_k) / (1 + 1j * omega * Lam))


def t_hist(name, on, off):
    """Bin-wise t statistics between two replica sets (complex)."""
    m_on, sre_on, sim_on, _ = summarize(on)
    m_off, sre_off, sim_off, _ = summarize(off)
    rel = np.maximum(np.hypot(sre_off, sim_off), np.hypot(sre_on, sim_on)) \
        / np.abs(m_off)
    mask = rel < 0.2
    t = np.concatenate([
        ((m_on.real - m_off.real) / np.hypot(sre_on, sre_off))[mask],
        ((m_on.imag - m_off.imag) / np.hypot(sim_on, sim_off))[mask]])
    ok = (abs(t.mean()) <= 0.3 and 0.8 <= t.std(ddof=1) <= 1.2
          and np.sum(np.abs(t) >= 3) <= np.floor(0.01 * t.size)) \
        if t.size >= 30 else bool(np.all(np.abs(t) < 3))
    log(f"  {name}: n={t.size} mean={t.mean():+.2f} std={t.std(ddof=1):.2f} "
        f"max|t|={np.abs(t).max():.2f} {'standard-normal' if ok else 'NOT standard-normal'}")
    return ok, t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=5000)
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--slab-particles", type=int, default=2000)
    p.add_argument("--slab-batches", type=int, default=10)
    p.add_argument("--only", default="5a,5c,5b,5d,5e")
    args = p.parse_args()
    only = args.only.split(",")
    WORK.mkdir(exist_ok=True)
    s3.WORK = WORK
    s4.WORK = WORK
    gate = []

    class A:
        reps = args.reps
        particles = args.particles
        batches = args.batches

    # ------------------------------------------------------------ 5a
    if "5a" in only:
        log("Test 5a: TRIGA-like infinite medium with cancellation (4x4x4)")
        mat = s3.triga_1g()
        lib = WORK / "lib_3b.h5"
        make_library({"m": mat}, lib)
        Lam = 1 / (2.1097e5 * 0.1)
        mesh = cube_mesh(4)
        orig = s3.cube_model
        for f in [0.1, 1, 10, 100]:
            eta = 3.0 if f >= 100 else 1.0
            w = 2 * np.pi * f
            ref = inf_solve(mat, w, [1.0], 1.0)
            s3.cube_model = lambda *a, **k: orig(
                *a, noise_extra={"cancellation_mesh": mesh}, **k)
            kt = ktilde_cb(w, eta, TRIGA_BETA_K, KEEPIN_U235_LAMBDA, Lam)
            e = s3.run_case(f"5a_{f}Hz", mat, lib, f, eta, 1.0, ref, A, [],
                            kt=None)
            s3.cube_model = orig
            Kp = physical_K(w, TRIGA_BETA_K, KEEPIN_U235_LAMBDA, Lam)
            log(f"    measured ratio {e['ratio']:.5f}; physical K {Kp:.5f}; "
                f"K~cb {kt:.5f}")
            gate.append((e["label"], e["ok"]))

    # ------------------------------------------------------------ 5c
    if "5c" in only:
        log("Test 5c: conditionally convergent cases with cancellation")
        mat = s3.triga_1g()
        lib = WORK / "lib_3b.h5"
        make_library({"m": mat}, lib)
        mesh = cube_mesh(4)
        orig = s3.cube_model
        s3.cube_model = lambda *a, **k: orig(
            *a, noise_extra={"cancellation_mesh": mesh}, **k)
        w = 2 * np.pi * 100
        ref = inf_solve(mat, w, [1.0], 1.0)
        e = s3.run_case("5c_3c_100Hz_eta1", mat, lib, 100.0, 1.0, 1.0, ref,
                        A, [])
        gate.append((e["label"], e["ok"]))
        # 1e4 Hz point of 3a at eta = 1
        mat = ferrero_case()
        lib = WORK / "lib_3a.h5"
        make_library({"m": mat}, lib)
        f = 1e4
        w = 2 * np.pi * f
        amp = np.pi * 5e-3 * upsilon(w, [7e-3], [0.08])
        ref = inf_solve(mat, w, [amp], 1.0)

        class B:
            reps = args.reps
            particles = 500000
            batches = 4
        e = s3.run_case("5c_3a_1e4Hz_eta1", mat, lib, f, 1.0, 1.0, ref, B, [],
                        strength=abs(amp), phase=np.angle(amp))
        gate.append((e["label"], e["ok"]))
        s3.cube_model = orig

    # ------------------------------------------------------- 5b, 5d, 5e
    def slab_runs(tag, lib, f, h_list, shifts=(0.0,)):
        regions = [(60.0, "m")]
        k = s4.Slab([(60.0, s4.mat_4a(), 480)], x0=-30.0,
                    n_angles=64).eigenvalue()[0]
        res = {}
        mk = lambda extra: (lambda: s4.noise_model(
            lib, regions, -30.0, 1, f, 1.0, k, (14, 15), 0, 60,
            args.slab_particles, args.slab_batches, noise_extra=extra))
        out = run_replicas(mk(None), WORK / f"{tag}_off_{f}", n_rep=args.reps)
        res["off"] = (out["values"]["tl"].reshape(args.reps, -1), out["times"])
        for h in h_list:
            for sh in shifts:
                extra = {"cancellation_mesh": slab_mesh(h, sh * h)}
                out = run_replicas(mk(extra), WORK / f"{tag}_on_{f}_{h}_{sh}",
                                   n_rep=args.reps)
                res[(h, sh)] = (out["values"]["tl"].reshape(args.reps, -1),
                                out["times"])
        return res

    if "5b" in only or "5d" in only:
        lib = WORK / "lib_4a.h5"
        make_library({"m": s4.mat_4a()}, lib)
    if "5b" in only:
        log("Test 5b: bias of mesh averaging, slab 4a at 1 Hz")
        res = slab_runs("5b", lib, 1.0, [2.0, 1.0, 0.5], shifts=(0.0, 0.5))
        for (h, sh), (vals, _) in [(k, v) for k, v in res.items() if k != "off"]:
            ok, _ = t_hist(f"h={h} shift={sh}h", vals, res["off"][0])
            if h < 2.0:
                gate.append((f"5b h={h} shift={sh}", ok))
            else:
                log(f"    (h=2 cm recorded as bias scale, not gated)")
    if "5d" in only:
        log("Test 5d: FOM with and without cancellation, slab 4a")
        for f in [0.1, 1.0]:
            res = slab_runs("5d", lib, f, [1.0])
            vals_off, t_off = res["off"]
            vals_on, t_on = res[(1.0, 0.0)]
            mask = np.abs(vals_off.mean(0)) > 0.1 * np.abs(vals_off.mean(0)).max()
            fr_off, fi_off = fom(vals_off, t_off)
            fr_on, fi_on = fom(vals_on, t_on)
            ratio_re = np.median(fr_on[mask] / fr_off[mask])
            ratio_im = np.median(fi_on[mask] / fi_off[mask])
            log(f"  {f} Hz: median FOM ratio on/off Re={ratio_re:.2f} "
                f"Im={ratio_im:.2f} (T off {t_off.sum():.0f}s, on "
                f"{t_on.sum():.0f}s)")
            gate.append((f"5d {f} Hz FOM ratio > 1",
                         ratio_re > 1 and ratio_im > 1))
    if "5e" in only:
        log("Test 5e: histogram scattering library (32 bins, mu_bar=0.3)")
        libh = WORK / "lib_4a_hist.h5"
        make_library({"m": s4.mat_4a()}, libh, histogram=(0.3, 32))
        res = slab_runs("5e", libh, 1.0, [1.0, 0.5])
        for (h, sh), (vals, _) in [(k, v) for k, v in res.items() if k != "off"]:
            ok, _ = t_hist(f"hist h={h}", vals, res["off"][0])
            gate.append((f"5e h={h}", ok))

    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 5 ({args.only}): {len(gate) - n_fail}/{len(gate)} gate "
        "checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_")
    with open(HERE / f"stage5_log_{tag}.md", "w") as fh:
        fh.write("# Stage 5 log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
