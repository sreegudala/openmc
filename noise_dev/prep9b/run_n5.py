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

LIBS = {"vector": "casmo70_hist32_kin.h5",
        "matrix": "casmo70_hist32_kin_matrix.h5"}
K_FORWARD = {"vector": (1.00137, 0.00009), "matrix": (1.00152, 0.00009)}
K_CE = (1.0008, 0.0005)
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


def run_kadj(args):
    work = HERE / f"work_n5_kadj_{args.form}"
    work.mkdir(exist_ok=True)
    lib = adjoint_library(args.form, work)
    make = lambda: mg_model(lib, args.particles, args.inactive, args.active)

    def extra(sp, d):
        return sp.keff.nominal_value, sp.keff.std_dev

    out = run_replicas(make, work / "reps", n_rep=args.reps, tally_names=(),
                       base_seed=args.base_seed, extra=extra)
    k = np.array(out["extras"])
    path = HERE / f"n5_kadj_{args.form}.npz"
    np.savez(path, k=k[:, 0], k_std=k[:, 1], times=out["times"])
    log(f"transposed ({args.form}) k: {k[:, 0].mean():.5f} +- "
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
    (HERE / "n5_log_kadj.md").write_text(
        "# N5 transposed k\n\n```\n" + "\n".join(LOG) + "\n```\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["kadj", "analyze"])
    p.add_argument("--form", choices=list(LIBS), default="vector")
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--particles", type=int, default=200000)
    p.add_argument("--inactive", type=int, default=100)
    p.add_argument("--active", type=int, default=400)
    args = p.parse_args()
    log(f"TRIGA_DIR {triga_dir()}, TRIGA_LIB_DIR {lib_dir()}")
    if args.mode == "kadj":
        run_kadj(args)
    else:
        analyze(args)


if __name__ == "__main__":
    main()
