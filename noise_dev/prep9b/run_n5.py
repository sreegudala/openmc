"""N5: checks of the shared TRIGA library on the TRIGA model (brief
noise_stage7_brief.md N5; design doc Stage 9b preparation; prep9b/n5_log.md).

Modes:

  kadj      R multigroup eigenvalue replicas of the TRIGA model on the
            transposed library (N1 + N2, built here from the forward one),
            with a flat weight window (the transposed library's second
            moment is unbounded without one: multiplicity up to 3.08)
            -> n5_kadj_{form}.npz
  analyze   transposed k against ce_mg_triga's forward k of the same form
            -> n5_log_kadj.md
  g0        R forward noise replicas of the TRIGA model at one frequency:
            10 % absorption increase in the central thimble water (cell
            3789, r = 1.69 cm, full height), relative noise of the total
            fission rate dF/F0 -> n5_g0_{f}Hz.npz
  g0-analyze  point-kinetic limit: dF/F0 = G0(omega) drho with one real
            drho fitted over all frequencies, G0 from the IFP beta_eff,k and
            Lambda_eff (kinetics.json) and the library's tallied lambda_k;
            also with the library's physical beta_k -> n5_log_g0.md

The forward multigroup k is ce_mg_triga's own (user decision 2026-11-26):
results/mgxs/kinetics_lib.md section 9 (commit e26de6a), mean of 5 seeds,
100k x (100 + 200) per seed: vector 1.00137 +- 0.00009, matrix 1.00152 +-
0.00009; CE 1.0008 +- 0.0005.

Locations (environment, else defaults):
  TRIGA_DIR      ce_mg_triga checkout (model_zoned/model.xml); default
                 <parent of parent of $OPENMC_DEV_SRC>/ce_mg_triga (LS6),
                 else ~/TRIGA_proposal_work/ce_mg_triga
  TRIGA_LIB_DIR  the two libraries; default $TRIGA_DIR/results/mgxs_ang,
                 else ~/dev/triga_lib
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))

from harness import run_replicas  # noqa: E402
import adjoint_noise_lib as anl  # noqa: E402

sys.path.insert(0, str(HERE.parent / "reference"))
import ref_g0  # noqa: E402
import h5py  # noqa: E402

LIBS = {"vector": "casmo70_hist32_kin.h5",
        "matrix": "casmo70_hist32_kin_matrix.h5"}
K_FORWARD = {"vector": (1.00137, 0.00009), "matrix": (1.00152, 0.00009)}
K_CE = (1.0008, 0.0005)
# ce_mg_triga results/mgxs/kinetics_lib.md sections 4-5: core values
LAMBDA_K = np.array([0.01334, 0.03272, 0.12075, 0.30291, 0.85035, 2.85514])
BETA_PHYS_K = np.array([22.5, 116.6, 111.3, 249.7, 103.1, 43.1]) * 1e-5
# ce_mg_triga results/kinetics/kinetics.json (IFP, CE)
BETA_EFF_K = ref_g0.TRIGA_BETA_K
LAMBDA_EFF = ref_g0.TRIGA_LAMBDA
CENTRAL_CELL = 3789                     # central thimble water (Water 2295)
CENTRAL_MATERIAL = "Water_2295"
EPS_ABS = 0.1
F7 = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
# cancellation mesh over the fuel box (zones.json fuel segments)
FUEL_BOX = ((-24.44, -25.77, -39.25), (24.44, 25.77, 19.05))
LOG = []


def log(msg=""):
    print(msg, flush=True)
    LOG.append(msg)


def triga_dir():
    d = os.environ.get("TRIGA_DIR")
    if d:
        return Path(d)
    src = os.environ.get("OPENMC_DEV_SRC")
    if src:
        cand = Path(src).resolve().parent.parent / "ce_mg_triga"
        if (cand / "model_zoned" / "model.xml").exists():
            return cand
    return Path.home() / "TRIGA_proposal_work" / "ce_mg_triga"


def lib_dir():
    d = os.environ.get("TRIGA_LIB_DIR")
    if d:
        return Path(d)
    cand = triga_dir() / "results" / "mgxs_ang"
    if (cand / LIBS["vector"]).exists():
        return cand
    return Path.home() / "dev" / "triga_lib"


def adjoint_library(form, work):
    src = lib_dir() / LIBS[form]
    dst = (work / f"adj_{form}.h5").resolve()
    if not dst.exists() or dst.stat().st_mtime < src.stat().st_mtime:
        rep = anl.make_adjoint_noise_library(src, dst)
        log(f"transposed {src.name} -> {dst.name}: {rep['materials']} "
            f"materials, max multiplicity {rep['max_multiplicity']:.3f}")
    return dst


def central_amplitudes(library):
    """'total' relative amplitudes for a pure absorption change
    dSigma_a = EPS_ABS Sigma_a of the central thimble water."""
    with h5py.File(library, "r") as f:
        g = f[CENTRAL_MATERIAL]
        t = g[next(k for k in g if k != "kTs")]
        return EPS_ABS * t["absorption"][()] / t["total"][()]


def noise_model(library, frequency, particles, batches, nskip,
                static_inactive, cancellation=(25, 25, 8)):
    model = mg_model(library, particles, 0, batches, ww=None)
    s = model.settings
    s.run_mode = "fixed source"
    s.batches = batches
    s.inactive = 0
    cell = model.geometry.get_all_cells()[CENTRAL_CELL]
    macro = getattr(cell.fill, "_macroscopic", None)
    if macro != CENTRAL_MATERIAL:
        raise SystemExit(f"cell {CENTRAL_CELL} holds {macro}, expected "
                         f"{CENTRAL_MATERIAL}")
    eta = 1.0 if frequency <= 30 else 3.0
    noise = {"frequency": frequency, "eta": eta, "keff": K_FORWARD["vector"][0],
             "nskip": nskip, "static_inactive": static_inactive,
             "perturbations": [{"cells": [cell], "total": [
                 complex(a) for a in central_amplitudes(library)]}]}
    if cancellation:
        cm = openmc.RegularMesh()
        cm.lower_left, cm.upper_right = FUEL_BOX
        cm.dimension = cancellation
        noise["cancellation_mesh"] = cm
    s.noise = noise
    t = openmc.Tally(name="fission")
    t.filters = [openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
    t.scores = ["fission"]
    model.tallies = openmc.Tallies([t])
    return model


def run_g0(args):
    f = args.freq
    work = HERE / f"work_n5_g0_{f:g}Hz"
    work.mkdir(exist_ok=True)
    lib = (lib_dir() / LIBS["vector"]).resolve()
    make = lambda: noise_model(lib, f, args.g0_particles, args.batches,
                               args.nskip, args.static_inactive)
    out = run_replicas(make, work / "reps", n_rep=args.reps,
                       tally_names=("fission",), base_seed=args.base_seed)
    v = out["values"]["fission"].reshape(args.reps, 2)
    path = HERE / f"n5_g0_{f:g}Hz.npz"
    np.savez(path, static=v[:, 0].real, noise=v[:, 1], frequency=f,
             times=out["times"])
    r = v[:, 1] / v[:, 0].real
    log(f"G0 run {f:g} Hz: dF/F0 = {r.mean():.4e} (sigma Re "
        f"{r.real.std(ddof=1) / np.sqrt(len(r)):.1e}, Im "
        f"{r.imag.std(ddof=1) / np.sqrt(len(r)):.1e}), T="
        f"{out['times'].sum():.0f}s -> {path.name}")


def g0_analyze(args):
    log("N5 point-kinetic limit: total fission rate dF/F0 against "
        "G0(omega) drho, drho real, fitted over all frequencies")
    log(f"  IFP beta_eff {BETA_EFF_K.sum() * 1e5:.1f} pcm, Lambda_eff "
        f"{LAMBDA_EFF * 1e6:.2f} us; library physical beta "
        f"{BETA_PHYS_K.sum() * 1e5:.1f} pcm; tallied lambda_k {LAMBDA_K}")
    lam_eff = ref_g0.lambda_eff(BETA_EFF_K, LAMBDA_K)
    rows = []
    for f in F7:
        path = HERE / f"n5_g0_{f:g}Hz.npz"
        if not path.exists():
            continue
        d = np.load(path)
        r = d["noise"] / d["static"]
        R = r.size
        sre = r.real.std(ddof=1) / np.sqrt(R)
        sim = r.imag.std(ddof=1) / np.sqrt(R)
        rows.append((f, r.mean(), sre, sim, R, d["times"].sum()))
    if not rows:
        log("  no n5_g0_*Hz.npz")
        return
    w = 2 * np.pi * np.array([x[0] for x in rows])
    r = np.array([x[1] for x in rows])
    sre = np.array([x[2] for x in rows])
    sim = np.array([x[3] for x in rows])
    for label, beta in (("IFP beta_eff", BETA_EFF_K),
                        ("library physical beta", BETA_PHYS_K)):
        G = ref_g0.g0(w, beta, LAMBDA_K, LAMBDA_EFF)
        # least squares for real drho: residual r - G drho, weights 1/sigma^2
        wr, wi = 1 / sre**2, 1 / sim**2
        num = np.sum(wr * G.real * r.real + wi * G.imag * r.imag)
        den = np.sum(wr * G.real**2 + wi * G.imag**2)
        drho = num / den
        sdrho = 1 / np.sqrt(den)
        res = r - G * drho
        chi2 = np.sum(wr * res.real**2 + wi * res.imag**2)
        dof = 2 * len(r) - 1
        log(f"\n  {label}: drho = {drho * 1e5:+.2f} +- {sdrho * 1e5:.2f} pcm"
            f"; chi2 {chi2:.1f} / {dof}")
        for (f, m, a, b, R, T), g in zip(rows, G):
            p = g * drho
            sA = np.hypot(m.real * a, m.imag * b) / abs(m)
            sph = np.degrees(np.hypot(m.imag * a, m.real * b) / abs(m)**2)
            log(f"    {f:6g} Hz: dF/F0 {abs(m):.4e} at "
                f"{np.degrees(np.angle(m)):+8.2f} deg; G0 drho {abs(p):.4e}"
                f" at {np.degrees(np.angle(p)):+8.2f} deg; amp ratio - 1 "
                f"{abs(m) / abs(p) - 1:+.4f} (sigma {sA / abs(p):.4f}), "
                f"dphase {np.degrees(np.angle(m / p)):+.2f} deg (sigma "
                f"{sph:.2f}); R {R}, T {T:.0f}s")
    # diagnostics of the IFP residual: a beta/Lambda scale (degenerate: only
    # their ratio is seen once drho is free) and a reactivity offset rho0,
    # dF/F0 = drho / (1/G0 - rho0), the signature of noise.keff differing
    # from the model's own k
    from scipy.optimize import least_squares
    diag = {
        "drho, Lambda scale": (
            lambda p: ref_g0.g0(w, BETA_EFF_K, LAMBDA_K, LAMBDA_EFF * p[1])
            * p[0], [-3e-4, 1.0]),
        "drho, rho0": (
            lambda p: p[0] / (1 / ref_g0.g0(w, BETA_EFF_K, LAMBDA_K,
                                            LAMBDA_EFF) - p[1]),
            [-3e-4, 0.0]),
    }
    log("\n  diagnostic fits (IFP beta_eff):")
    for name, (fn, p0) in diag.items():
        fit = least_squares(lambda p: np.concatenate(
            [(fn(p).real - r.real) / sre, (fn(p).imag - r.imag) / sim]),
            p0, x_scale="jac")
        sig = np.sqrt(np.diag(np.linalg.inv(fit.jac.T @ fit.jac)))
        q = r / fn(fit.x)
        log(f"    {name}: p = " + ", ".join(
            f"{v:.5g} +- {e:.2g}" for v, e in zip(fit.x, sig))
            + f"; chi2 {2 * fit.cost:.1f} / {2 * len(r) - len(p0)}")
        log("      amp ratio - 1 [%]: " + " ".join(
            f"{100 * (abs(x) - 1):+.1f}" for x in q) + "; dphase [deg]: "
            + " ".join(f"{np.degrees(np.angle(x)):+.2f}" for x in q))
    log("\n  collision-scheme cut-offs (ref_g0.cutoffs, IFP values):")
    for eta in (1.0, 3.0):
        log(f"    eta {eta:g}: " + str(ref_g0.cutoffs(
            eta, BETA_EFF_K.sum(), LAMBDA_EFF, lam_eff)))
    (HERE / "n5_log_g0.md").write_text(
        "# N5 point-kinetic limit\n\n```\n" + "\n".join(LOG) + "\n```\n")


def mg_model(library, particles, inactive, active, ww=(0.5, 2.0)):
    """ce_mg_triga's zoned model in multigroup form (as its
    step3a_adjoint.mg_model), eigenvalue mode, with a flat weight window."""
    model = openmc.Model.from_model_xml(triga_dir() / "model_zoned"
                                        / "model.xml")
    model.tallies = openmc.Tallies([])
    model.convert_to_multigroup(method="material_wise", groups="CASMO-70",
                                mgxs_path=str(library),
                                overwrite_mgxs_library=False)
    model.materials.cross_sections = str(library)
    s = model.settings
    s.run_mode = "eigenvalue"
    s.particles, s.inactive, s.batches = particles, inactive, inactive + active
    if ww is not None:
        bb = model.geometry.bounding_box
        mesh = openmc.RegularMesh()
        mesh.lower_left, mesh.upper_right = bb.lower_left, bb.upper_right
        mesh.dimension = (1, 1, 1)
        s.weight_windows = [openmc.WeightWindows(
            mesh, lower_ww_bounds=[ww[0]], upper_ww_bounds=[ww[1]],
            particle_type="neutron")]
        s.weight_windows_on = True
        s.max_history_splits = 1000
    return model


def kadj_tag(form, forward=False, seed=1):
    """n5_kadj_{form}[_fwd][_s{seed}]: the N5 runs have no suffix."""
    return (f"{form}" + ("_fwd" if forward else "")
            + ("" if seed == 1 else f"_s{seed}"))


def run_kadj(args):
    tag = kadj_tag(args.form, args.kadj_forward, args.base_seed)
    work = HERE / f"work_n5_kadj_{tag}"
    work.mkdir(exist_ok=True)
    # --kadj-forward: the forward library under the same settings (flat
    # window, seeds), for a like-for-like difference (Stage 9b pilot)
    lib = ((lib_dir() / LIBS[args.form]).resolve() if args.kadj_forward
           else adjoint_library(args.form, work))
    make = lambda: mg_model(lib, args.particles, args.inactive, args.active)

    def extra(sp, d):
        return sp.keff.nominal_value, sp.keff.std_dev

    out = run_replicas(make, work / "reps", n_rep=args.reps, tally_names=(),
                       base_seed=args.base_seed, extra=extra)
    k = np.array(out["extras"])
    path = HERE / f"n5_kadj_{tag}.npz"
    np.savez(path, k=k[:, 0], k_std=k[:, 1], times=out["times"])
    which = "forward" if args.kadj_forward else "transposed"
    log(f"{which} ({args.form}) k: {k[:, 0].mean():.5f} +- "
        f"{k[:, 0].std(ddof=1) / np.sqrt(len(k)):.5f} ({len(k)} replicas, "
        f"T={out['times'].sum():.0f}s) -> {path.name}")


def analyze(args):
    log("N5: transposed library k on the TRIGA model against ce_mg_triga's "
        "forward k (kinetics_lib.md section 9)")
    log(f"  CE k {K_CE[0]:.5f} +- {K_CE[1]:.5f}")
    for form in LIBS:
        path = HERE / f"n5_kadj_{form}.npz"
        if not path.exists():
            log(f"  {form}: no {path.name}")
            continue
        d = np.load(path)
        R = d["k"].size
        ka = d["k"].mean()
        sa = max(d["k"].std(ddof=1), d["k_std"].mean()) / np.sqrt(R)
        kf, sf = K_FORWARD[form]
        diff = (ka - kf) * 1e5
        sd = np.hypot(sa, sf) * 1e5
        z = diff / sd
        log(f"  {form}: transposed {ka:.5f} +- {sa:.5f} ({R} replicas; "
            f"spread {d['k'].std(ddof=1) * 1e5:.0f} pcm, mean within-run "
            f"sigma {d['k_std'].mean() * 1e5:.0f} pcm), forward {kf:.5f} +- "
            f"{sf:.5f}: transposed - forward {diff:+.0f} +- {sd:.0f} pcm "
            f"({z:+.1f} sigma) -> {'consistent' if abs(z) < 3 else 'DIFFER'}"
            f"; against CE {(ka - K_CE[0]) * 1e5:+.0f} pcm")
    # like-for-like reruns (Stage 9b pilot): transposed and forward
    # libraries under the same settings, another seed
    def stat(d):
        R = d["k"].size
        return d["k"].mean(), max(d["k"].std(ddof=1),
                                  d["k_std"].mean()) / np.sqrt(R), R
    for fp in sorted(HERE.glob("n5_kadj_*_fwd_s*.npz")):
        form, seed = fp.stem.split("_")[2], fp.stem.split("_s")[-1]
        ap = HERE / f"n5_kadj_{kadj_tag(form, False, int(seed))}.npz"
        if not ap.exists():
            log(f"  {form} seed {seed}: no {ap.name}")
            continue
        (ka, sa, Ra), (kf, sf, Rf) = stat(np.load(ap)), stat(np.load(fp))
        diff, sd = (ka - kf) * 1e5, np.hypot(sa, sf) * 1e5
        log(f"  {form}, seed {seed}, same settings: transposed {ka:.5f} +- "
            f"{sa:.5f} ({Ra}), forward {kf:.5f} +- {sf:.5f} ({Rf}): "
            f"transposed - forward {diff:+.0f} +- {sd:.0f} pcm "
            f"({diff / sd:+.1f} sigma)")
    (HERE / "n5_log_kadj.md").write_text(
        "# N5 transposed k\n\n```\n" + "\n".join(LOG) + "\n```\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["kadj", "analyze", "g0", "g0-analyze"])
    p.add_argument("--freq", type=float, default=1.0)
    p.add_argument("--g0-particles", type=int, default=10000)
    p.add_argument("--batches", type=int, default=25)
    p.add_argument("--nskip", type=int, default=5)
    p.add_argument("--static-inactive", type=int, default=100)
    p.add_argument("--form", choices=list(LIBS), default="vector")
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--particles", type=int, default=200000)
    p.add_argument("--inactive", type=int, default=100)
    p.add_argument("--active", type=int, default=400)
    p.add_argument("--kadj-forward", action="store_true",
                   help="kadj with the forward library (same settings)")
    args = p.parse_args()
    log(f"TRIGA_DIR {triga_dir()}, TRIGA_LIB_DIR {lib_dir()}")
    if args.mode == "kadj":
        run_kadj(args)
    elif args.mode == "analyze":
        analyze(args)
    elif args.mode == "g0":
        run_g0(args)
    else:
        g0_analyze(args)


if __name__ == "__main__":
    main()
