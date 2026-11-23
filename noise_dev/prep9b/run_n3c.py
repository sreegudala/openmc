"""Stage 9b preparation, N3c: weight windows with complex weights (design
doc, amendment 2026-11-22).

In noise mode the windows are absolute (in units of the noise source
weight scale), compare |w|, and split with the phase kept; inside the
window mesh they replace the noise roulette.

N3c-a  Deep penetration. Fission-free 2-group water slab (CORTEX water,
       inverse velocities of the CORTEX fuel), 60 cm = 15 fast / 73 thermal
       mean free paths, vacuum at both faces, reflective sides. Unit
       isotropic fast source in [0, 2] cm; detector = thermal flux in
       [58, 60] cm. 1 kHz (detector phase -32 deg) and 10 kHz (-101 deg: Re
       is negative at the detector and changes sign along the slab). eta 1
       at 1 kHz, 5 at 10 kHz: in water (thermal c = 0.976) the collision
       game's within-generation |w| ratio c~ = (c + sqrt(1 + eta^2) x) /
       (1 + eta x), x = omega / (v Sigma_t), is 1.025 at 10 kHz and eta 1
       (divergent, the Stage 2 cut-off) and 0.993 at eta 5.
       Windows: CADIS lower bounds C / |chi| on 1 cm x 2 groups, chi the
       ref_sn1d transposed solution for the detector at the same omega, C
       placing a unit source particle at its window centre (|chi| averaged
       over the source cells); upper/lower 5, max_split 10.
       Gate: STD of the windowed detector response against ref_sn1d; the
       two-sample test between windows on and off; FOM ratio on/off with a
       95 % bootstrap lower bound above 1. The mesh field against S_N
       (STD-C) is logged.
N3c-b  Adjoint efficiency. The 7a infinite-medium adjoint cases (N3
       material, CORTEX fuel) at 1, 30 and 100 Hz, with group windows
       C / |psi_g|, psi = the forward solution for the adjoint source R (the
       importance for <chi, R>), C placing the source group's unit weight
       at its window centre. Runs: adjoint without and with windows, and
       the forward control (forward library, same source group, no
       windows). Gate: STD against the ref_inf adjoint with windows on.
       Logged: FOM ratio with/without windows, and adjoint-with-windows
       against forward.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for sub in ("common", "reference", "stage3", "stage7"):
    sys.path.insert(0, str(ROOT / sub))

from adjoint_noise_lib import make_adjoint_noise_library  # noqa: E402
from build import box_source, make_library, slab_model, x_mesh  # noqa: E402
from harness import (fmt_info, fmt_info_correlated, run_replicas,  # noqa: E402
                     std_correlated, std_criterion, summarize,
                     two_sample_correlated)
from noise_common import Material  # noqa: E402
from ref_inf import (cortex_fuel, cortex_water, solve as inf_solve,  # noqa
                     static_k)
from ref_sn1d import Slab  # noqa: E402
import run_stage3 as s3  # noqa: E402
from run_stage7_prep import n3_material  # noqa: E402

WORK = HERE / "work"
LOG = []
SAVED = {}

L_SLAB, SRC, DET = 60.0, (0.0, 2.0), (58.0, 60.0)
NB = 60
FAST, THERMAL = 0, 1          # library order
E_FAST, E_THERMAL = 1, 0      # EnergyFilter / window energy bins (ascending)
RATIO, MAX_SPLIT = 5.0, 10
OFF_SEED = 500009             # seed offset of the windows-off set


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def fmt_c(z):
    return f"{z.real:+.5e}{z.imag:+.5e}j"


def fom_mod(vals, times):
    """FOM of the modulus of a complex scalar: 1 / (rel sigma^2 x T)."""
    m = vals.mean()
    s = np.hypot(vals.real.std(ddof=1), vals.imag.std(ddof=1)) \
        / np.sqrt(vals.size)
    if m == 0 or s == 0:       # no replica scored (or all equal)
        return (0.0 if m == 0 else np.inf), np.inf if m == 0 else 0.0
    return 1.0 / ((s / abs(m))**2 * times.sum()), s / abs(m)


def fom_ratio_ci(a, ta, b, tb, n_boot=4000, seed=1):
    """FOM(a)/FOM(b) and its 95 % bootstrap interval (replicas resampled
    with their times)."""
    rng = np.random.default_rng(seed)
    r0 = fom_mod(a, ta)[0] / fom_mod(b, tb)[0]
    rs = []
    for _ in range(n_boot):
        i = rng.integers(0, a.size, a.size)
        j = rng.integers(0, b.size, b.size)
        fb = fom_mod(b[j], tb[j])[0]
        rs.append(np.inf if fb == 0 else fom_mod(a[i], ta[i])[0] / fb)
    rs = np.array(rs)
    rs = rs[~np.isnan(rs)]
    return (r0, np.percentile(rs, 2.5, method="lower"),
            np.percentile(rs, 97.5, method="higher"))


def water():
    w = cortex_water()
    return Material(w.total, w.scatter, cortex_fuel().inv_v)


# ---------------------------------------------------------------- N3c-a
def slab_windows(sl, f, mesh):
    """CADIS lower bounds C / |chi| on the 1 cm x 2-group mesh."""
    q = sl.isotropic_source(DET, [0.0, 1.0])
    chi = sl.solve(2 * np.pi * f, q, k=1.0, adjoint=True)       # (N, G)
    edges = np.linspace(0, L_SLAB, NB + 1)
    chi_b = np.array([np.abs(chi[(sl.x > lo) & (sl.x < hi)]).mean(axis=0)
                      for lo, hi in zip(edges[:-1], edges[1:])])  # (NB, G)
    src = (sl.x > SRC[0]) & (sl.x < SRC[1])
    centre = 2.0 / (1.0 + RATIO)           # lower bound at the window centre
    C = np.abs(chi[src, FAST]).mean() * centre
    lower = C / chi_b[:, ::-1]             # energy bins ascending
    return openmc.WeightWindows(
        mesh, lower_ww_bounds=lower.reshape(NB, 1, 1, 2),
        upper_bound_ratio=RATIO, energy_bounds=[0.0, 1.0, 2.0],
        max_split=MAX_SPLIT), lower


def eta_a(f):
    return 5.0 if f >= 5000 else 1.0


def slab_model_n3c(lib, f, particles, batches, ww=None):
    model, edges = slab_model([(L_SLAB, "w")], x0=0.0, library=lib)
    s = model.settings
    s.particles = particles
    s.batches = batches
    s.noise = {"frequency": f, "eta": eta_a(f), "keff": 1.0}
    s.source = [box_source(SRC[0], SRC[1], FAST, 2)]
    if ww is not None:
        s.weight_windows = ww
    t = openmc.Tally(name="tl")
    t.filters = [openmc.MeshFilter(x_mesh(0.0, L_SLAB, NB)),
                 openmc.EnergyFilter(np.arange(3.0)),
                 openmc.ComplexPartFilter()]
    t.scores = ["flux"]
    model.tallies = openmc.Tallies([t])
    return model


def test_a(args, gate):
    w = water()
    lib = WORK / "lib_water.h5"
    make_library({"w": w}, lib)
    sl = Slab([(L_SLAB, w, int(L_SLAB * 4))], x0=0.0, n_angles=64)
    edges = np.linspace(0, L_SLAB, NB + 1)
    q = sl.isotropic_source(SRC, [1.0, 0.0])
    log(f"N3c-a: water slab {L_SLAB:.0f} cm ({L_SLAB * w.total[0]:.0f} fast /"
        f" {L_SLAB * w.total[1]:.0f} thermal mfp), fast source {SRC}, "
        f"thermal detector {DET}")
    for f in [float(x) for x in args.freqs_a.split(",")]:
        ref_f = sl.solve(2 * np.pi * f, q, k=1.0)
        ref_mesh = np.array([ref_f[(sl.x > lo) & (sl.x < hi)].mean(axis=0)
                             for lo, hi in zip(edges[:-1], edges[1:])])
        # ref_sn1d is per unit source density; OpenMC emits one particle
        # over the 2 cm x 10 x 10 source box and the tally integrates over
        # 1 cm x 10 x 10 bins: tally = (V_bin / V_src) x bin-average flux
        ref = ref_mesh[:, ::-1].reshape(-1) * (L_SLAB / NB) \
            / (SRC[1] - SRC[0])                         # [x][E ascending]
        # detector = [58, 60]: the last two 1 cm bins, thermal
        det_idx = [2 * (NB - 2) + E_THERMAL, 2 * (NB - 1) + E_THERMAL]
        ref_d = ref[det_idx].sum()
        mesh = x_mesh(0.0, L_SLAB, NB)
        ww, lower = slab_windows(sl, f, mesh)
        log(f"  {f:g} Hz (eta {eta_a(f):g}): S_N detector response "
            f"{fmt_c(ref_d)} "
            f"(|.| {abs(ref_d):.3e}, phase {np.degrees(np.angle(ref_d)):+.1f}"
            f" deg); window lower bounds {lower.min():.2e} .. "
            f"{lower.max():.2e}; Re of the S_N thermal field changes sign "
            f"{int(np.sum(np.diff(np.sign(ref_mesh[:, THERMAL].real)) != 0))}"
            f" times")
        outs = {}
        for tag, wwx, seed in (("on", ww, args.base_seed),
                               ("off", None, args.base_seed + OFF_SEED)):
            outs[tag] = run_replicas(
                lambda: slab_model_n3c(lib, f, args.particles, args.batches,
                                       wwx),
                WORK / f"a_{tag}_{f:g}", n_rep=args.reps, base_seed=seed)
        v = {t: o["values"]["tl"].reshape(args.reps, -1)
             for t, o in outs.items()}
        d = {t: x[:, det_idx].sum(axis=1) for t, x in v.items()}
        T = {t: o["times"] for t, o in outs.items()}
        # STD of the windowed detector response against S_N
        m, sre, sim, _ = summarize(d["on"][:, None])
        ok1, info = std_criterion(m, sre, sim, np.array([ref_d]))
        log(f"    detector, windows on vs S_N: {'PASS' if ok1 else 'FAIL'} "
            f"{fmt_info(info)}; mean {fmt_c(m[0])}")
        gate.append((f"N3c-a {f:g} Hz windows vs S_N", ok1))
        # two-sample on/off
        n_off = int(np.sum(np.abs(d["off"]) > 0))
        ok2, inf2 = two_sample_correlated(d["on"][:, None],
                                          d["off"][:, None])
        log(f"    detector, windows on vs off: {'PASS' if ok2 else 'FAIL'} "
            f"Welch t Re={inf2['z'][0]:+.2f} Im={inf2['z'][1]:+.2f}, "
            f"Hotelling p={inf2['hotelling'][2]:.3f}; off mean "
            f"{fmt_c(d['off'].mean())} ({n_off}/{args.reps} replicas score)")
        gate.append((f"N3c-a {f:g} Hz windows on vs off", ok2))
        # FOM
        r, lo, hi = fom_ratio_ci(d["on"], T["on"], d["off"], T["off"])
        f_on, s_on = fom_mod(d["on"], T["on"])
        f_off, s_off = fom_mod(d["off"], T["off"])
        ok3 = lo > 1.0
        log(f"    FOM (detector modulus): on {f_on:.3g} (sigma {s_on:.4f}, "
            f"T {T['on'].sum():.0f}s), off {f_off:.3g} (sigma {s_off:.4f}, "
            f"T {T['off'].sum():.0f}s); ratio {r:.3g} [95 % {lo:.3g}, "
            f"{hi:.3g}] -> {'PASS' if ok3 else 'FAIL'}")
        gate.append((f"N3c-a {f:g} Hz FOM ratio > 1", ok3))
        # field, logged
        mask = np.abs(ref) > 1e-3 * np.abs(ref).max()
        okf, inff = std_correlated(v["on"], ref, mask)
        log(f"    field, windows on vs S_N (logged): {'PASS' if okf else 'FAIL'}"
            f" (STD-C) {fmt_info_correlated(inff)}")
        for t in ("on", "off"):
            SAVED[f"a_{f:g}_{t}"] = v[t]
            SAVED[f"a_{f:g}_{t}_times"] = T[t]
        SAVED[f"a_{f:g}_ref"] = ref


# ---------------------------------------------------------------- N3c-b
def cube_windows(psi, src_group):
    """Group windows C / |psi_g| on one mesh cell covering the cube;
    psi in library order, src_group's unit weight at its window centre."""
    G = psi.size
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10.0, -10.0, -10.0)
    mesh.upper_right = (10.0, 10.0, 10.0)
    mesh.dimension = (1, 1, 1)
    centre = 2.0 / (1.0 + RATIO)
    lower = np.abs(psi[src_group]) * centre / np.abs(psi)
    return openmc.WeightWindows(
        mesh, lower_ww_bounds=lower[::-1].reshape(1, 1, 1, G),
        upper_bound_ratio=RATIO, energy_bounds=np.arange(G + 1.0),
        max_split=MAX_SPLIT), lower


def run_cube(label, mat, lib, f, eta, k, particles, args, ww=None, seed=0):
    G = mat.G

    def make():
        m = s3.cube_model(G, f, eta, k, particles, args.batches, lib, G - 1)
        if ww is not None:
            m.settings.weight_windows = ww
        return m
    out = run_replicas(make, WORK / label, n_rep=args.reps,
                       tally_names=("tl", "col"),
                       base_seed=args.base_seed + seed,
                       extra=lambda sp, d: s3.read_noise(
                           d / f"statepoint.{args.batches}.h5"))
    vals = out["values"]["tl"].reshape(args.reps, -1)
    gens = np.concatenate([e["generations_per_batch"]
                           for e in out["extras"]]).mean()
    return vals, out["times"], gens


def test_b(args, gate):
    log("N3c-b: infinite-medium adjoint with group windows (7a cases)")
    parts = dict(kv.split("=") for kv in args.particles_b.split(","))
    for name, mat in (("n3", n3_material()), ("cortex", cortex_fuel())):
        fwd = WORK / f"lib_b_{name}_forward.h5"
        adj = WORK / f"lib_b_{name}_adjoint.h5"
        make_library({"m": mat}, fwd)
        make_adjoint_noise_library(fwd, adj)
        k, _ = static_k(mat)
        G = mat.G
        R = np.eye(G)[G - 1]
        n_adj = int(parts[name])
        n_fwd = int(parts[name + "_fwd"])
        for f in [float(x) for x in args.freqs_b.split(",")]:
            w = 2 * np.pi * f
            eta = 5.0 if f >= 100 else 1.0
            psi = inf_solve(mat, w, R, k=k)                    # importance
            ww, lower = cube_windows(psi, G - 1)
            ref_a = inf_solve(mat, w, R, k=k, adjoint=True)[::-1]
            ref_f = psi[::-1]
            lab = f"b_{name}_{f:g}"
            va_on, ta_on, g_on = run_cube(lab + "_adj_on", mat, adj, f, eta,
                                          k, n_adj, args, ww)
            va_off, ta_off, g_off = run_cube(lab + "_adj_off", mat, adj, f,
                                             eta, k, n_adj, args, None,
                                             OFF_SEED)
            vf, tf, g_f = run_cube(lab + "_fwd", mat, fwd, f, eta, k, n_fwd,
                                   args, None, 2 * OFF_SEED)
            m, sre, sim, _ = summarize(va_on)
            ok, info = std_criterion(m, sre, sim, ref_a)
            log(f"  {name} {f:g} Hz (eta {eta:g}; window lower bounds "
                f"{np.array2string(lower, precision=3)} fastest first): "
                f"adjoint with windows vs ref: {'PASS' if ok else 'FAIL'} "
                f"{fmt_info(info)}")
            gate.append((f"N3c-b {name} {f:g} Hz adjoint with windows", ok))
            m2, sre2, sim2, _ = summarize(va_off)
            _, info2 = std_criterion(m2, sre2, sim2, ref_a)
            m3, sre3, sim3, _ = summarize(vf)
            _, info3 = std_criterion(m3, sre3, sim3, ref_f)
            log(f"    adjoint without windows vs ref (logged): "
                f"{fmt_info(info2)}")
            log(f"    forward control vs ref (logged): {fmt_info(info3)}")
            # FOM on the source-group response (the adjoint <chi, R>; the
            # forward <phi, R>, i.e. both the thermal bin, tally index 0)
            rows = []
            for tag, v, t, g in (("adj on", va_on, ta_on, g_on),
                                 ("adj off", va_off, ta_off, g_off),
                                 ("fwd", vf, tf, g_f)):
                fo, s = fom_mod(v[:, 0], t)
                rows.append((tag, fo, s, t.sum(), g))
            log("    FOM (source-group response): " + "; ".join(
                f"{tag} {fo:.3g} (sigma {s:.4f}, T {T:.0f}s, gens/batch "
                f"{g:.0f})" for tag, fo, s, T, g in rows))
            r1 = fom_ratio_ci(va_on[:, 0], ta_on, va_off[:, 0], ta_off)
            r2 = fom_ratio_ci(va_off[:, 0], ta_off, vf[:, 0], tf, seed=2)
            r3 = fom_ratio_ci(va_on[:, 0], ta_on, vf[:, 0], tf, seed=3)
            log(f"    FOM ratios [95 %]: adj on/off {r1[0]:.3g} "
                f"[{r1[1]:.3g}, {r1[2]:.3g}]; adj off/fwd {r2[0]:.3g} "
                f"[{r2[1]:.3g}, {r2[2]:.3g}]; adj on/fwd {r3[0]:.3g} "
                f"[{r3[1]:.3g}, {r3[2]:.3g}]")
            for tag, v, t in (("adj_on", va_on, ta_on),
                              ("adj_off", va_off, ta_off), ("fwd", vf, tf)):
                SAVED[f"{lab}_{tag}"] = v
                SAVED[f"{lab}_{tag}_times"] = t
            SAVED[f"{lab}_ref_adj"] = ref_a
            SAVED[f"{lab}_ref_fwd"] = ref_f


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--only", default="a,b")
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=20000,
                   help="particles per batch, N3c-a")
    p.add_argument("--particles-b",
                   default="n3=20000,n3_fwd=20000,cortex=200000,"
                   "cortex_fwd=20000")
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--freqs-a", default="1000,10000")
    p.add_argument("--freqs-b", default="1,30,100")
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--tag", default="")
    args = p.parse_args()
    only = args.only.split(",")
    global WORK
    sfx = f"_{args.tag}" if args.tag else ""
    if args.tag:
        WORK = HERE / f"work_{args.tag}"
    WORK.mkdir(parents=True, exist_ok=True)
    s3.WORK = WORK
    s3.LOG = LOG
    gate = []
    if "a" in only:
        test_a(args, gate)
    if "b" in only:
        test_b(args, gate)
    n_fail = sum(not ok for _, ok in gate)
    log(f"\nN3c ({args.only}): {len(gate) - n_fail}/{len(gate)} gate checks "
        "passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_") + sfx
    with open(HERE / f"n3c_log_{tag}.md", "w") as fh:
        fh.write("# N3c log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    np.savez(HERE / f"n3c_vals_{tag}.npz", **SAVED)
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
