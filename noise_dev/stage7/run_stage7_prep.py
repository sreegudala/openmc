"""Stage 7 preparation (brief noise_stage7_brief.md, N3a and N3b).

Infinite medium (20 cm reflective cube, as in Stage 3), 3-group material
whose transposed library has multiplicity 2.4 in the thermal group
(absorption' zeroed there by the N1 rebalance) and a finite second moment
(per-collision second-moment radius 0.33, see
adjoint_noise_lib.second_moment_radius). The adjoint source R is a
uniform isotropic unit source in the thermal group. Each case runs OpenMC on
the transposed library at +omega and compares with ``ref_inf.solve(...,
adjoint=True)``, i.e. L(omega)^T chi = R (bilinear convention, no
conjugation).

N3a  fission-free: the copied ce_mg_triga transposition as is
     (adjoint_lib.make_adjoint_library); 100 Hz (plateau) and 3 kHz (break)
     with eta 1, 10 kHz (above) with eta 3: at eta 1 the thermal copy loop
     has gain x/(1+x) sqrt(2) = 1.10 > 1 (x = omega/(v Sigma_t) = 3.6), the
     Stage 2 collision-based cut-off, the same for forward and adjoint.
N3b  fission and six delayed families: adjoint_noise_lib (transposed
     production matrices), noise keff = k_inf of the material; 1 Hz
     (plateau), 30 Hz (break), 100 Hz (above, eta 3), eta 1 otherwise.

N3f  forward control, gated like N3b: the same material and frequencies
     run forward (forward library, thermal unit source, reference
     ref_inf.solve forward). Gives the forward/adjoint FOM ratio on one
     material.
N3x  informational, not gated: the first candidate material, whose
     transpose (multiplicity 2.4, thermal self-transfer probability 0.35)
     has second-moment radius 2.05, i.e. infinite variance. Fission-free,
     100 Hz. Records what a transposed library does without splitting or
     weight windows.

Pass: STD (design doc Section 5) for both estimators, as in Stages 2-3.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "common"))
sys.path.insert(0, str(ROOT / "reference"))
sys.path.insert(0, str(ROOT / "stage3"))

import adjoint_lib  # noqa: E402
import run_stage3 as s3  # noqa: E402
from adjoint_noise_lib import (effective_operator, infinite_k,  # noqa: E402
                               make_adjoint_noise_library,
                               second_moment_radius)
from build import make_library  # noqa: E402
from noise_common import Material  # noqa: E402
from ref_g0 import KEEPIN_U235_LAMBDA, TRIGA_BETA_K  # noqa: E402
from ref_inf import solve as inf_solve, static_k  # noqa: E402

WORK = HERE / "work"
LOG = []
R = np.array([0.0, 0.0, 1.0])     # adjoint source: thermal group


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def n3_material():
    """3-group material; its transpose has multiplicity 2.425 in group 3
    and second-moment radius 0.33. Thermal speed 2.2e4 cm/s gives
    Lambda = 37 us and a point-kinetic break near 30 Hz."""
    S = np.array([[0.05, 0.20, 0.04],
                  [0.00, 0.20, 1.80],
                  [0.00, 0.02, 0.10]])
    return Material.from_beta(
        total=[0.30, 2.20, 0.80], scatter=S,
        inv_v=1.0 / np.array([1e9, 1e7, 2.2e4]),
        nu_fission=[0.005, 0.02, 1.2], chi=[0.97, 0.03, 0.0],
        beta=TRIGA_BETA_K, lam=KEEPIN_U235_LAMBDA)


def n3x_material():
    """First candidate: transpose has multiplicity 2.4 in group 3 and
    second-moment radius 2.05 (infinite variance)."""
    S = np.array([[0.05, 0.20, 0.04],
                  [0.00, 0.20, 1.20],
                  [0.00, 0.02, 0.68]])
    return Material.from_beta(
        total=[0.30, 1.50, 0.80], scatter=S,
        inv_v=1.0 / np.array([1e9, 1e7, 2.2e5]),
        nu_fission=[0.005, 0.02, 0.12], chi=[0.97, 0.03, 0.0],
        beta=TRIGA_BETA_K, lam=KEEPIN_U235_LAMBDA)


def log_radius(label, fwd, adj):
    log(f"  {label} second-moment radius: forward "
        f"{second_moment_radius(effective_operator(fwd, 'm')):.3f}, "
        f"transposed {second_moment_radius(effective_operator(adj, 'm')):.3f}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=20000)
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--only", default="n3a,n3b,n3x")
    p.add_argument("--tag", default="")
    args = p.parse_args()
    only = args.only.split(",")
    global WORK
    if args.tag:
        WORK = HERE / f"work_{args.tag}"
    WORK.mkdir(parents=True, exist_ok=True)
    s3.WORK = WORK
    s3.LOG = LOG

    mat = n3_material()
    fwd = WORK / "lib_n3_forward.h5"
    make_library({"m": mat}, fwd)
    results, gate = [], []

    if "n3a" in only:
        adj = WORK / "lib_n3a_adjoint.h5"
        rep = adjoint_lib.make_adjoint_library(fwd, adj)
        log(f"Test N3a: fission-free transposed library (copied step3a), "
            f"report {rep}")
        log_radius("N3a", fwd, adj)
        m0 = Material(mat.total, mat.scatter, mat.inv_v)
        for f, eta in ((100.0, 1.0), (3000.0, 1.0), (10000.0, 3.0)):
            ref = inf_solve(m0, 2 * np.pi * f, R, k=1.0, adjoint=True)[::-1]
            e = s3.run_case(f"n3a_{f:g}Hz", mat, adj, f, eta, 1.0, ref, args,
                            results, src_group=2)
            gate.append((e["label"], e["ok"]))

    if "n3b" in only:
        adj = WORK / "lib_n3b_adjoint.h5"
        rep = make_adjoint_noise_library(fwd, adj)
        k, _ = static_k(mat)
        k_adj = infinite_k(effective_operator(adj, "m"))
        log(f"Test N3b: transposed library with fission and delayed data, "
            f"report {rep}")
        log(f"  k_inf forward (ref_inf) {k:.12f}, transposed library "
            f"{k_adj:.12f}")
        log_radius("N3b", fwd, adj)
        for f, eta in ((1.0, 1.0), (30.0, 1.0), (100.0, 3.0)):
            ref = inf_solve(mat, 2 * np.pi * f, R, k=k, adjoint=True)[::-1]
            e = s3.run_case(f"n3b_{f:g}Hz", mat, adj, f, eta, k, ref, args,
                            results, src_group=2)
            gate.append((e["label"], e["ok"]))

    if "n3f" in only:
        k, _ = static_k(mat)
        log("Test N3f: forward control on the N3 material (forward library)")
        for f, eta in ((1.0, 1.0), (30.0, 1.0), (100.0, 3.0)):
            ref = inf_solve(mat, 2 * np.pi * f, R, k=k)[::-1]
            e = s3.run_case(f"n3f_{f:g}Hz", mat, fwd, f, eta, k, ref, args,
                            results, src_group=2)
            gate.append((e["label"], e["ok"]))

    if "n3x" in only:
        mx = n3x_material()
        fwd_x = WORK / "lib_n3x_forward.h5"
        adj_x = WORK / "lib_n3x_adjoint.h5"
        make_library({"m": mx}, fwd_x)
        rep = adjoint_lib.make_adjoint_library(fwd_x, adj_x)
        log(f"Test N3x (informational, not gated): infinite-variance "
            f"transpose, report {rep}")
        log_radius("N3x", fwd_x, adj_x)
        m0 = Material(mx.total, mx.scatter, mx.inv_v)
        ref = inf_solve(m0, 2 * np.pi * 100.0, R, k=1.0, adjoint=True)[::-1]
        s3.run_case("n3x_100Hz", mx, adj_x, 100.0, 1.0, 1.0, ref, args, [],
                    src_group=2)

    zall = np.concatenate([np.concatenate([r["tl"]["z"], r["col"]["z"]])
                           for r in results]) if results else np.zeros(0)
    if zall.size:
        log(f"Pooled z: n={zall.size} mean={zall.mean():+.2f} "
            f"std={zall.std(ddof=1):.2f} |z|>=3: "
            f"{int(np.sum(np.abs(zall) >= 3))}")
    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 7 prep ({args.only}): {len(gate) - n_fail}/{len(gate)} "
        "checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_") + (f"_{args.tag}" if args.tag else "")
    with open(HERE / f"stage7_prep_log_{tag}.md", "w") as fh:
        fh.write("# Stage 7 preparation log (N3a, N3b)\n\n```\n"
                 + "\n".join(LOG) + "\n```\n")
    np.savez(HERE / f"stage7_prep_vals_{tag}.npz",
             **{r["label"] + "_" + est: r[est]["vals"] for r in results
                for est in ("tl", "col")})
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
