"""Stage 9a: C5G7 2-D noise (FEMFFUSION's examples/2D_C5G7_NN/NN1.2, a 5 %
capture oscillation in one MOX pin at 1 Hz), OpenMC against FEMFFUSION SP3
(design doc Stage 9; amendment 2026-11-25).

Modes (run separately; each writes an npz in this directory):

  femffusion  FEMFFUSION on its own example: the shipped configurations
              ref1_FE2 (the reference: mesh_ref1, one radial refinement, FE
              degree 2, SP3) and ref0_FE1, plus ref1_FE2 with SP1 and SP5,
              all with the full noise operator (-matrixfree_type_A
              full_allocated; amendment 2026-11-26), and the shipped default
              operator (non_diagonal) on ref1_FE2 SP3 and SP5, logged
              -> stage9a_femffusion.npz
  eigen       R OpenMC k-eigenvalue replicas -> stage9a_eigen{tag}.npz
              (k per replica; the noise runs divide fission by their mean)
  openmc      R OpenMC noise replicas -> stage9a_openmc{tag}.npz
  analyze     gate -> stage9a_log{tag}.md

Compared quantity: the relative noise dphi_g / phi_0,g per pin cell (51 x
51) and group (7), independent of either code's flux normalisation. OpenMC
takes phi_0 from the static sampling generation of each batch
(NoisePhaseFilter), FEMFFUSION from its static SP3 solution.

Gate (amendment 2026-11-25): on every pin cell outside the perturbed
assembly and every group where the OpenMC replica mean is precise enough
(amplitude sigma < 5/3 %, phase sigma < 2/3 deg, as in harness
std_criterion), amplitude within 5 % and phase within 2 deg of FEMFFUSION
ref1_FE2 SP3 with the full noise operator; at least half of those cells evaluated in every group; OpenMC's
k known to 3 pcm. Logged: the perturbed assembly, core sums, the static
shape, FEMFFUSION's discretisation (ref0_FE1) and SPN-order (SP1, SP5)
spreads, the two cancellation meshes against each other, and the same
comparison against the shipped default (non_diagonal) noise operator.

Amendment 2026-11-26: FEMFFUSION's default non_diagonal SPN noise operator
omits the cross-moment Marshak terms marshack_coeff[mi][mj] (mi != mj) at
vacuum boundaries (its off-diagonal moment blocks are volume-only mass
operators), which its static SPN operator includes. The noise and static
boundary conditions then differ: dphi/phi_0 falls towards the vacuum
boundary, more at SP5 than SP3 and not at SP1. full_allocated assembles the
same equations with those terms.

FEMFFUSION is GPL-3.0 and lives outside the repo ($FEMFFUSION_DIR, built by
tacc/build_femffusion.sh); its example inputs are read in place.
"""

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE))

from harness import run_replicas, std_criterion, summarize, \
    two_sample_correlated  # noqa
import c5g7  # noqa

LOG = []
N, G = c5g7.N, c5g7.G
FULL = "-matrixfree_type_A full_allocated"
FEMFF_CONFIGS = {             # name: (shipped prm, N_SPN, extra options)
    "ref1_FE2_SP3": ("c5g7_ref1_FE2_DSP3.prm", 3, FULL),
    "ref0_FE1_SP3": ("c5g7_ref0_FE1_DSP3.prm", 3, FULL),
    "ref1_FE2_SP1": ("c5g7_ref1_FE2_DSP3.prm", 1, FULL),
    "ref1_FE2_SP5": ("c5g7_ref1_FE2_DSP3.prm", 5, FULL),
    "ref1_FE2_SP3_nondiag": ("c5g7_ref1_FE2_DSP3.prm", 3, ""),
    "ref1_FE2_SP5_nondiag": ("c5g7_ref1_FE2_DSP3.prm", 5, ""),
}
SHIPPED = "ref1_FE2_SP3_nondiag"
REFERENCE = "ref1_FE2_SP3"
PERT_ASSEMBLY = (slice(0, 17), slice(17, 34))     # rows (y), columns (x)
PERT_PIN = (15, 18)
AMP_TOL, PHASE_TOL = 0.05, 2.0
K_SIGMA_MAX = 3e-5


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


# ---------------------------------------------------------------- FEMFFUSION

def parse_femffusion_out(path):
    lines = Path(path).read_text().splitlines()

    def block(title, cplx):
        i = lines.index(title)
        rows = lines[i + 1:i + 1 + N]
        conv = complex if cplx else float
        return np.array([[conv(t) for t in r.split()] for r in rows])

    i = lines.index("The Eigenvalues are: ")
    k = float(lines[i + 1].split()[0])
    static = np.array([block(f"Group {g + 1} flux", False) for g in range(G)])
    noise = np.array([block(f"Flux Noise Group {g + 1}", True)
                      for g in range(G)])
    return k, static, noise


def run_femffusion(args):
    root = c5g7.femffusion_dir()
    exe = root / "femffusion.exe"
    if not exe.exists():
        raise SystemExit(f"{exe} not found (tacc/build_femffusion.sh)")
    work = (HERE / "work_9a_femffusion").resolve()
    work.mkdir(exist_ok=True)
    configs = args.configs.split(",") if args.configs else list(FEMFF_CONFIGS)
    for name in configs:
        prm_name, nspn, opts = FEMFF_CONFIGS[name]
        text = (c5g7.example_dir() / prm_name).read_text()
        text = re.sub(r"(set N_SPN\s*=\s*)\d+", rf"\g<1>{nspn}", text)
        text = re.sub(r"(set Print_Grid_Flag\s*=\s*)\w+", r"\g<1>false", text)
        text = re.sub(r"(set Output_Filename\s*=\s*)\S+",
                      rf"\g<1>{work / name}.out", text)
        prm = work / f"{name}.prm"
        prm.write_text(text)
        cmd = f"OMP_NUM_THREADS=1 ./femffusion.exe -f {prm} {opts}"
        if (root / "modules.sh").exists():     # LS6
            cmd = f"source modules.sh && {cmd}"
        t0 = time.time()
        with open(work / f"{name}.log", "w") as fh:
            subprocess.run(["bash", "-c", cmd], cwd=root, stdout=fh,
                           stderr=subprocess.STDOUT, check=True)
        dt = time.time() - t0
        save_femffusion(name, work / f"{name}.out", dt)


def save_femffusion(name, out_file, dt=0.0):
    """Parse one FEMFFUSION output into stage9a_femffusion.npz (saved after
    every configuration, so a later failure keeps the earlier ones)."""
    k, static, noise = parse_femffusion_out(out_file)
    log(f"FEMFFUSION {name}: k = {k:.5f}, {dt:.0f} s")
    path = HERE / "stage9a_femffusion.npz"
    old = dict(np.load(path)) if path.exists() else {}
    old[f"{name}_k"] = np.array(k)
    old[f"{name}_static"] = static
    old[f"{name}_noise"] = noise
    old[f"{name}_time"] = np.array(dt)
    np.savez(path, **old)
    log(f"  -> {path.name}")


# ---------------------------------------------------------------- OpenMC

def library(work):
    lib = (work / "lib_c5g7.h5").resolve()
    c5g7.library(lib)
    return lib


def run_eigen(args):
    work = HERE / f"work_9a_eigen{args.sfx}"
    work.mkdir(exist_ok=True)
    lib = library(work)
    make = lambda: c5g7.model(lib, args.eigen_particles,
                              args.eigen_inactive + args.eigen_active,
                              run_mode="eigenvalue",
                              inactive=args.eigen_inactive)

    def extra(sp, d):      # real tally: not through complex_tally_stats
        k = sp.keff
        flux = sp.get_tally(name="flux").mean.ravel()
        return k.nominal_value, k.std_dev, flux

    out = run_replicas(make, work / "reps", n_rep=args.reps,
                       tally_names=(), base_seed=args.base_seed,
                       extra=extra)
    ks = np.array([e[:2] for e in out["extras"]])
    flux = c5g7.to_grid(np.array([e[2] for e in out["extras"]]))
    path = HERE / f"stage9a_eigen{args.sfx}.npz"
    np.savez(path, k=ks[:, 0], k_std=ks[:, 1], static=flux,
             times=out["times"])
    kmean, ksd = ks[:, 0].mean(), ks[:, 0].std(ddof=1) / np.sqrt(len(ks))
    log(f"OpenMC eigen: {args.reps} replicas x {args.eigen_particles} x "
        f"{args.eigen_active} active: k = {kmean:.6f} +- {ksd:.6f} "
        f"(replica spread; mean within-run sigma "
        f"{ks[:, 1].mean() / np.sqrt(len(ks)):.6f}), T={out['times'].sum():.0f}s"
        f" -> {path.name}")


def eigen_k(tag):
    path = HERE / f"stage9a_eigen{tag}.npz"
    e = np.load(path)
    return e["k"].mean(), e["k"].std(ddof=1) / np.sqrt(e["k"].size)


def run_openmc(args):
    keff = args.keff
    if keff is None:
        keff, ksd = eigen_k(args.eigen_sfx)
        log(f"keff from stage9a_eigen{args.eigen_sfx}.npz: {keff:.6f} +- "
            f"{ksd:.6f}")
    work = HERE / f"work_9a_openmc{args.sfx}"
    work.mkdir(exist_ok=True)
    lib = library(work)
    make = lambda: c5g7.model(lib, args.particles, args.batches, keff=keff,
                              nskip=args.nskip,
                              static_inactive=args.static_inactive,
                              cancellation=args.cancellation)
    out = run_replicas(make, work / "reps", n_rep=args.reps,
                       tally_names=("flux",), base_seed=args.base_seed)
    v = out["values"]["flux"].reshape(args.reps, N * N * G, 2)
    path = HERE / f"stage9a_openmc{args.sfx}.npz"
    np.savez(path, noise=c5g7.to_grid(v[..., 1]),
             static=c5g7.to_grid(v[..., 0].real), keff=np.array(keff),
             cancellation=np.array(args.cancellation or 0),
             times=out["times"])
    log(f"OpenMC noise: {args.reps} replicas x {args.batches} batches x "
        f"{args.particles}, cancellation {args.cancellation}, keff {keff:.6f}"
        f", T={out['times'].sum():.0f}s -> {path.name}")


# ---------------------------------------------------------------- analysis

def outside_mask():
    m = np.ones((N, N), bool)
    m[PERT_ASSEMBLY] = False
    return m


def fmt(z):
    return f"{abs(z):.5f} at {np.degrees(np.angle(z)):+7.2f} deg"


def analyze(args):
    F = np.load(HERE / "stage9a_femffusion.npz")
    ref_rel = F[f"{REFERENCE}_noise"] / F[f"{REFERENCE}_static"]
    kref = float(F[f"{REFERENCE}_k"])
    kmean, ksd = eigen_k(args.eigen_sfx)
    log(f"Stage 9a: C5G7 noise, OpenMC against FEMFFUSION {REFERENCE}")
    log(f"  k: OpenMC {kmean:.6f} +- {ksd:.6f} (C5G7 transport reference "
        f"1.18655); FEMFFUSION " + ", ".join(
            f"{n} {float(F[f'{n}_k']):.5f}" for n in FEMFF_CONFIGS
            if f"{n}_k" in F))
    gate = [("OpenMC k sigma <= 3 pcm", ksd <= K_SIGMA_MAX)]
    log(f"  k sigma {ksd * 1e5:.2f} pcm -> "
        f"{'PASS' if ksd <= K_SIGMA_MAX else 'FAIL'} (3 pcm)")

    out = outside_mask()
    fuel = np.zeros((N, N), bool)
    fuel[:34, :34] = True
    runs = {}
    for tag in args.openmc_tags.split(","):
        sfx = f"_{tag}" if tag else ""
        a = np.load(HERE / f"stage9a_openmc{sfx}.npz")
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = a["noise"] / a["static"]                # (R, G, N, N)
        # a cell without static flux in some replica is not evaluated
        runs[tag or "default"] = (a, rel)

    # FEMFFUSION's own spreads (logged)
    for name in FEMFF_CONFIGS:
        if name == REFERENCE or f"{name}_noise" not in F:
            continue
        r = F[f"{name}_noise"] / F[f"{name}_static"]
        msg = []
        for region, reg in (("fuel", out & fuel), ("reflector", out & ~fuel)):
            m = reg[None] & np.ones((G, 1, 1), bool)
            da = np.abs(np.abs(r[m]) / np.abs(ref_rel[m]) - 1)
            dp = np.abs(np.degrees(np.angle(r[m] / ref_rel[m])))
            msg.append(f"{region} max |amp ratio - 1| {da.max():.4f}, max "
                       f"|dphase| {dp.max():.2f} deg")
        log(f"  FEMFFUSION {name} vs {REFERENCE} (outside the perturbed "
            f"assembly): " + "; ".join(msg))

    for label, (a, rel) in runs.items():
        R = rel.shape[0]
        log(f"\n  OpenMC run '{label}': {R} replicas, cancellation "
            f"{int(a['cancellation'])}, keff {float(a['keff']):.6f}, "
            f"T={a['times'].sum():.0f}s")
        if abs(float(a["keff"]) - kmean) > 1e-7:
            log("    WARNING: run keff differs from the eigen mean")
        rule3 = []
        mean, s_re, s_im, _ = summarize(rel.reshape(R, -1))
        mean = mean.reshape(G, N, N)
        s_re = s_re.reshape(G, N, N)
        s_im = s_im.reshape(G, N, N)
        for g in range(G):
            with np.errstate(divide="ignore", invalid="ignore"):
                ok, info = std_criterion(mean[g], s_re[g], s_im[g],
                                         ref_rel[g], AMP_TOL, PHASE_TOL,
                                         mask=out & np.isfinite(mean[g]))
            n_cells = int(out.sum())
            enough = min(info["amp_eval"], info["ph_eval"]) >= n_cells / 2
            okg = info["ok_amp"] and info["ok_ph"] and enough
            gate.append((f"9a '{label}' group {g + 1}", okg))
            with np.errstate(divide="ignore", invalid="ignore"):
                _, info3 = std_criterion(mean[g], s_re[g], s_im[g],
                                         ref_rel[g], AMP_TOL, PHASE_TOL,
                                         mask=out & np.isfinite(mean[g]),
                                         significance=3)
            rule3.append(f"g{g + 1} "
                         f"{'pass' if info3['ok_amp'] and info3['ok_ph'] else 'FAIL'}"
                         f" ({info3['n_amp_fail']} amp, {info3['n_ph_fail']} "
                         f"phase cells)")
            fin = out & np.isfinite(mean[g])
            A = np.abs(mean[g][fin]) / np.abs(ref_rel[g][fin]) - 1
            P = np.degrees(np.angle(mean[g][fin] / ref_rel[g][fin]))
            log(f"    g{g + 1}: {'PASS' if okg else 'FAIL'}: evaluated amp "
                f"{info['amp_eval']}/{n_cells}, phase {info['ph_eval']}/"
                f"{n_cells}; max |amp ratio - 1| {info['max_rel_amp']:.4f}, "
                f"max |dphase| {info['max_dphi']:.2f} deg; all cells: amp "
                f"ratio - 1 median {np.median(A):+.4f} [5-95 % "
                f"{np.percentile(A, 5):+.4f}, {np.percentile(A, 95):+.4f}],"
                f" dphase median {np.median(P):+.2f} deg; z hist (logged) "
                f"mean {info['z_mean']:+.2f} std {info['z_std']:.2f}")
        log("    3-sigma tolerance rule (amendment 2026-11-26b, logged; 9a is "
            "graded by the rule above): " + ", ".join(rule3))
        # logged: the same check against the shipped default operator
        if f"{SHIPPED}_noise" in F:
            ship = F[f"{SHIPPED}_noise"] / F[f"{SHIPPED}_static"]
            res = []
            for g in range(G):
                with np.errstate(divide="ignore", invalid="ignore"):
                    _, info = std_criterion(mean[g], s_re[g], s_im[g],
                                            ship[g], AMP_TOL, PHASE_TOL,
                                            mask=out & np.isfinite(mean[g]))
                ok_s = info["ok_amp"] and info["ok_ph"]
                res.append(f"g{g + 1} {'pass' if ok_s else 'FAIL'} (max amp "
                           f"{info['max_rel_amp']:.3f})")
            log(f"    against {SHIPPED} (logged): " + ", ".join(res))
        # logged: perturbed pin, perturbed assembly, core sums
        for g in (0, 3, 6):
            r, c = PERT_PIN
            log(f"    g{g + 1} perturbed pin: OpenMC {fmt(mean[g, r, c])}, "
                f"FEMFFUSION {fmt(ref_rel[g, r, c])}")
        pa = np.zeros((N, N), bool)
        pa[PERT_ASSEMBLY] = True
        for g in (0, 6):
            A = np.abs(mean[g][pa]) / np.abs(ref_rel[g][pa]) - 1
            P = np.degrees(np.angle(mean[g][pa] / ref_rel[g][pa]))
            log(f"    g{g + 1} perturbed assembly (logged): amp ratio - 1 "
                f"[{A.min():+.4f}, {A.max():+.4f}], dphase [{P.min():+.2f}, "
                f"{P.max():+.2f}] deg")
        core = a["noise"].sum(axis=(2, 3)) / a["static"].sum(axis=(2, 3))
        cref = F[f"{REFERENCE}_noise"].sum(axis=(1, 2)) \
            / F[f"{REFERENCE}_static"].sum(axis=(1, 2))
        cm = core.mean(0)
        cs = np.hypot(core.real.std(0, ddof=1), core.imag.std(0, ddof=1)) \
            / np.sqrt(R)
        for g in range(G):
            log(f"    g{g + 1} core sum dphi/phi0: OpenMC {fmt(cm[g])} "
                f"(sigma {cs[g] / abs(cm[g]):.4f}), FEMFFUSION {fmt(cref[g])}"
                f"; ratio - 1 {abs(cm[g]) / abs(cref[g]) - 1:+.4f}, dphase "
                f"{np.degrees(np.angle(cm[g] / cref[g])):+.2f} deg")
        # static shape (logged)
        s0 = a["static"].mean(0)
        f0 = F[f"{REFERENCE}_static"]
        for g in (0, 6):
            sa = s0[g] / s0[g][fuel].sum()
            sb = f0[g] / f0[g][fuel].sum()
            q = np.abs(sa[fuel] / sb[fuel] - 1)
            log(f"    g{g + 1} static shape (fuel, sum-normalised): "
                f"|OpenMC/FEMFFUSION - 1| median {np.median(q):.4f}, max "
                f"{q.max():.4f}")

    # cancellation meshes against each other (logged)
    labels = list(runs)
    if len(labels) == 2:
        (_, r1), (_, r2) = runs[labels[0]], runs[labels[1]]
        for g in (0, 6):
            m = out.ravel()
            ok, info = two_sample_correlated(
                r1[:, g].reshape(r1.shape[0], -1),
                r2[:, g].reshape(r2.shape[0], -1), m)
            k, T2, p = info["hotelling"]
            log(f"  '{labels[0]}' vs '{labels[1]}' g{g + 1} (logged): "
                f"{'consistent' if ok else 'INCONSISTENT'}, max|z| "
                f"{np.abs(info['z']).max():.2f}, Hotelling p {p:.3f}")

    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 9a: {len(gate) - n_fail}/{len(gate)} gate checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    with open(HERE / f"stage9a_log{args.sfx}.md", "w") as fh:
        fh.write("# Stage 9a log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    return n_fail


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["femffusion", "eigen", "openmc",
                                    "analyze"])
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--tag", default="")
    p.add_argument("--configs", default="",
                   help="femffusion: comma-separated subset of "
                   + ",".join(FEMFF_CONFIGS))
    # eigen
    p.add_argument("--eigen-particles", type=int, default=200000)
    p.add_argument("--eigen-inactive", type=int, default=200)
    p.add_argument("--eigen-active", type=int, default=500)
    p.add_argument("--eigen-tag", default="",
                   help="openmc/analyze: tag of the eigen npz")
    # noise
    p.add_argument("--keff", type=float, default=None,
                   help="default: mean of stage9a_eigen{eigen-tag}.npz")
    p.add_argument("--particles", type=int, default=20000)
    p.add_argument("--batches", type=int, default=1000)
    p.add_argument("--nskip", type=int, default=20)
    p.add_argument("--static-inactive", type=int, default=200)
    p.add_argument("--cancellation", type=int, default=51,
                   help="n for an n x n cancellation mesh (0: off)")
    p.add_argument("--openmc-tags", default="c51,c102",
                   help="analyze: tags of the OpenMC npz files")
    args = p.parse_args()
    args.sfx = f"_{args.tag}" if args.tag else ""
    args.eigen_sfx = f"_{args.eigen_tag}" if args.eigen_tag else ""
    if args.mode == "femffusion":
        run_femffusion(args)
    elif args.mode == "eigen":
        run_eigen(args)
    elif args.mode == "openmc":
        run_openmc(args)
    else:
        raise SystemExit(1 if analyze(args) else 0)


if __name__ == "__main__":
    main()
