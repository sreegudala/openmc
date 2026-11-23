"""Stage 8c: vibrating pin in the CORTEX assembly against Abeille (design
doc Stage 8, amendments 2026-11-22b and 2026-11-23).

The pin at (+2.52, -2.52) vibrates along y with eps = 0.2 cm and
omega_0 = 2 pi (Abeille's noise_vibration.yaml; the paper's x-vibrating pin
at (-2.52, -2.52) is the same problem rotated by 90 degrees). Two swept
layers, one per y face: top (fuel below, water above) and bottom (water
below, fuel above). Harmonics n = 1, 2 at n Hz.

Modes
  openmc     forward runs: 138 x 138 x 2-group flux map, the central
             water-hole detector (thermal flux in the 1.26 cm lattice
             element at the centre), and the static moments of phi_0 on
             NL bins across each layer (for the contraction). Cancellation
             meshes 170 and 340 (--cancellation).
  adjoint    transposed library, unit thermal source in the central water
             hole (R = 1, source strength V_det), moments of chi on the
             same layer bins.
  abeille    Abeille noise_vibration.yaml with approximate cancellation on
             the shipped 170 x 170 cancelator (the published runs use
             cancellation; the shipped input has it off), 138 x 138 real
             and imaginary flux, harmonic n.
  analyze    Gate:
             - two-sample STD-C against Abeille per group on cells above
               10 % of the maximum modulus (both normalised by the static
               thermal flux at the assembly centre, as in 6c);
             - OpenMC at the two cancellation meshes mutually consistent
               (two-sample STD-C);
             - the detector response against the adjoint contraction
               <chi, dS> (P0 + P1 + fission; paired test).
             The probe mode logs the forward sigma for sizing
             (amendment 2026-11-22b).
"""

import argparse
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for sub in ("common", "reference", "stage6"):
    sys.path.insert(0, str(ROOT / sub))

import cortex  # noqa: E402
import harness  # noqa: E402
from adjoint_noise_lib import (make_adjoint_noise_library,  # noqa: E402
                               pt_fission_complex, pt_terms_complex)
from harness import (paired_correlated, run_replicas,  # noqa: E402
                     two_sample_correlated)
from noise_common import chi_omega  # noqa: E402
from ref_inf import cortex_fuel, cortex_water  # noqa: E402
import run_stage6c as s6c  # noqa: E402

LOG = []
NL = 8                                   # bins across each layer
EPS = 0.2
F0 = 1.0
HP = cortex.HALF_PIN
DET_HALF = cortex.PITCH / 2              # central water hole
V_DET = (2 * DET_HALF)**2 * 4.0          # cm^3 (z from -2 to 2)
ADJ_SEED_OFFSET = 1000003


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def fmt_c(z):
    return f"{z.real:+.5e}{z.imag:+.5e}j"


def layer_mesh(y0, eps=EPS):
    xc = cortex.PIN_CENTRE[0]
    m = openmc.RegularMesh()
    m.lower_left = (xc - HP, y0 - eps, -2.0)
    m.upper_right = (xc + HP, y0 + eps, 2.0)
    m.dimension = (1, NL, 1)
    return m


def det_mesh():
    m = openmc.RegularMesh()
    m.lower_left = (-DET_HALF, -DET_HALF, -2.0)
    m.upper_right = (DET_HALF, DET_HALF, 2.0)
    m.dimension = (1, 1, 1)
    return m


def add_tallies(model, static):
    """Detector and layer-moment tallies (static phase bin if static)."""
    energy = openmc.EnergyFilter([0.0, 1.0, 2.0])
    phase = [openmc.NoisePhaseFilter()] if static else []
    cp = [openmc.ComplexPartFilter()]
    t = openmc.Tally(name="det")
    t.filters = [openmc.MeshFilter(det_mesh()), energy] + phase + cp
    t.scores = ["flux"]
    tallies = list(model.tallies) + [t]
    for name, y0, _, _ in cortex.vibration_layers(EPS):
        t = openmc.Tally(name=f"pt_{name}")
        t.filters = [openmc.MeshFilter(layer_mesh(y0)), energy,
                     openmc.SphericalHarmonicsFilter(1)] + phase + cp
        t.scores = ["flux"]
        tallies.append(t)
    model.tallies = openmc.Tallies(tallies)
    return model


def forward(lib, n, args, cancellation):
    m = cortex.model(lib, args.particles, args.batches, frequency=n * F0,
                     nskip=args.nskip, static_inactive=args.static_inactive,
                     cancellation=cancellation,
                     vibration={"eps": EPS, "harmonic": n})
    return add_tallies(m, static=True)


def adjoint(lib_adj, n, args):
    m = cortex.model(lib_adj, args.particles_adj, args.batches_adj,
                     frequency=n * F0, library_only_noise=True)
    m.settings.source = [openmc.IndependentSource(
        space=openmc.stats.Box((-DET_HALF, -DET_HALF, -2.0),
                               (DET_HALF, DET_HALF, 2.0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([0.5], [1.0]))]   # thermal group
    return add_tallies(m, static=False)


TALLY_NAMES = ("flux", "det", "pt_top", "pt_bottom")


def run_openmc(args):
    work = HERE / f"work_8c{args.sfx}"
    work.mkdir(exist_ok=True)
    lib = (work / "lib_cortex.h5").resolve()
    cortex.library(lib)
    for n in args.harm:
        for canc in args.cancellation:
            out = run_replicas(lambda: forward(lib, n, args, canc),
                               work / f"fwd_n{n}_c{canc}", n_rep=args.reps,
                               tally_names=TALLY_NAMES,
                               base_seed=args.base_seed)
            v = out["values"]
            R = args.reps
            flux = v["flux"].reshape(R, 138, 138, 2, 2)
            path = HERE / f"stage8c_openmc_n{n}_c{canc}{args.sfx}.npz"
            np.savez(path, noise=flux[..., 1], static=flux[..., 0].real,
                     det=v["det"].reshape(R, 2, 2),
                     pt_top=v["pt_top"].reshape(R, NL, 2, 4, 2),
                     pt_bottom=v["pt_bottom"].reshape(R, NL, 2, 4, 2),
                     times=out["times"])
            log(f"OpenMC n={n} cancellation {canc}: {R} x {args.batches} "
                f"batches, T={out['times'].sum():.0f}s -> {path.name}")
            probe(path, n, canc)


def probe(path, n, canc):
    """Forward sigma of the gated quantities (amendment 2026-11-22b)."""
    a = np.load(path)
    R = a["noise"].shape[0]
    c = slice(68, 70)
    norm = a["static"][:, c, c, 1].mean(axis=(1, 2))       # thermal centre
    A = a["noise"] / norm[:, None, None, None]
    for g, name in ((1, "thermal"), (0, "fast")):
        x = A[..., g].reshape(R, -1)
        m = np.abs(x.mean(0))
        mask = m > 0.1 * m.max()
        s = np.hypot(x.real.std(0, ddof=1), x.imag.std(0, ddof=1)) \
            / np.sqrt(R)
        rel = s[mask] / m[mask]
        log(f"  probe n={n} c={canc} {name}: {mask.sum()} cells above 10 %;"
            f" rel sigma of the mean: median {np.median(rel):.3f}, 90th "
            f"pct {np.percentile(rel, 90):.3f}")
    d = a["det"][:, 0, 1] / norm                    # thermal (E bin 0), noise
    s = np.hypot(d.real.std(ddof=1), d.imag.std(ddof=1)) / np.sqrt(R)
    log(f"  probe n={n} c={canc} detector: {fmt_c(d.mean())}, rel sigma "
        f"{s / abs(d.mean()):.3f}; T per replica "
        f"{a['times'].mean():.0f}s")


def run_adjoint(args):
    work = HERE / f"work_8c{args.sfx}"
    work.mkdir(exist_ok=True)
    lib = (work / "lib_cortex.h5").resolve()
    lib_adj = (work / "lib_cortex_adjoint.h5").resolve()
    cortex.library(lib)
    make_adjoint_noise_library(lib, lib_adj)
    for n in args.harm:
        out = run_replicas(lambda: adjoint(lib_adj, n, args),
                           work / f"adj_n{n}", n_rep=args.reps,
                           tally_names=TALLY_NAMES,
                           base_seed=args.base_seed + ADJ_SEED_OFFSET)
        v = out["values"]
        R = args.reps
        path = HERE / f"stage8c_adjoint_n{n}{args.sfx}.npz"
        np.savez(path, flux=v["flux"].reshape(R, 138, 138, 2),
                 pt_top=v["pt_top"].reshape(R, NL, 2, 4),
                 pt_bottom=v["pt_bottom"].reshape(R, NL, 2, 4),
                 times=out["times"])
        log(f"Adjoint n={n}: {R} x {args.batches_adj} batches, "
            f"T={out['times'].sum():.0f}s -> {path.name}")


# ---------------------------------------------------------------- Abeille
def abeille_input(n, nbatches, nignored, seed, max_minutes, particles):
    text = (s6c.abeille_dir() / "input_files" /
            "noise_vibration.yaml").read_text()
    geometry = text[:text.index("tally-filters:")]
    tallies = "".join(s6c.FLUX_TALLY.format(name=q)
                      for q in ("real-flux", "imag-flux"))
    w = 2 * np.pi * n * F0
    sim = f"""simulation:
  mode: noise
  nparticles: {particles}
  nbatches: {nbatches}
  nignored: {nignored}
  nskip: 3
  noise-cancellation: true
  inner-generations: true
  normalize-noise-source: true
  noise-angular-frequency: {w:.8f}
  keff: 0.99916

  cancelator:
    type: approximate
    shape: [170, 170, 1]
    low: [-10.71, -10.71, -2.00]
    hi: [10.71, 10.71, 2.00]
    beta: average-g

""" + s6c.SOURCE + """
  noise-sources:
    - type: flat-vibration
      low: [2.1543, -2.3543, -2.00]
      hi: [2.8857, -1.9543, 2.00]
      angular-frequency: 6.28318531 # 2 pi
      direction: Y
      positive-material: 2 # Water
      negative-material: 1 # Fuel

    - type: flat-vibration
      low: [2.1543, -3.0857, -2.00]
      hi: [2.8857, -2.6857, 2.00]
      angular-frequency: 6.28318531 # 2 pi
      direction: Y
      positive-material: 1 # Fuel
      negative-material: 2 # Water
"""
    return (geometry + s6c.TALLIES.format(tallies=tallies, seed=seed,
                                          max_minutes=max_minutes) + sim)


def run_abeille(args):
    exe = s6c.abeille_exe()
    for n in args.harm:
        work = HERE / f"work_8c_abeille_n{n}{args.sfx}"
        work.mkdir(exist_ok=True)
        dirs = []
        for r in range(args.reps):
            d = work / f"rep{r:02d}"
            d.mkdir(exist_ok=True)
            seed = 19073486328125 + 2 * 7919 * (r + args.base_seed)
            (d / "input.yaml").write_text(abeille_input(
                n, args.abeille_batches, args.abeille_ignored, seed,
                args.max_minutes, args.abeille_particles))
            dirs.append(d)

        def one(d):
            t0 = time.time()
            with open(d / "abeille.log", "w") as fh:
                subprocess.run([str(exe), "--input", "input.yaml",
                                "--threads", str(harness.THREADS),
                                "--output", "out.h5"], cwd=d, stdout=fh,
                               stderr=subprocess.STDOUT, check=True)
            return time.time() - t0
        with ThreadPoolExecutor(max_workers=min(harness.JOBS,
                                                args.reps)) as ex:
            times = list(ex.map(one, dirs))
        vals, nb = [], []
        for d in dirs:
            with h5py.File(d / "out.h5", "r") as f:
                r = f["results"]
                vals.append(r["real-flux/avg"][()]
                            + 1j * r["imag-flux/avg"][()])
            nb.append((d / "abeille.log").read_text()
                      .count("Simulated noise batch"))
        out = HERE / f"stage8c_abeille_n{n}{args.sfx}.npz"
        np.savez(out, vals=np.array(vals), batches=np.array(nb),
                 times=np.array(times))
        log(f"Abeille n={n}: {args.reps} replicas, batches {nb} -> "
            f"{out.name}")


# ---------------------------------------------------------------- analyze
def normalised_openmc(path):
    a = np.load(path)
    c = slice(68, 70)
    norm = a["static"][:, c, c, 1].mean(axis=(1, 2))
    # [R, iy, ix, E] -> [R, g, iy, ix] with g 0 = fast
    A = a["noise"] / norm[:, None, None, None]
    return np.stack([A[..., 1], A[..., 0]], axis=1), a, norm


def contraction(a, adj, n, norm):
    """<chi, dS> per replica over both layers (P0 + P1 + fission)."""
    fuel, water = cortex_fuel(), cortex_water()
    mats = {"fuel": fuel, "water": water}
    w = 2 * np.pi * n * F0
    out = np.zeros((a["det"].shape[0], 3), complex)
    for name, y0, neg, pos in cortex.vibration_layers(EPS):
        L, Rm = mats[neg], mats[pos]
        e = np.linspace(y0 - EPS, y0 + EPS, NL + 1)
        gx, gw = np.polynomial.legendre.leggauss(16)
        y = e[:-1, None] + 0.5 * (gx[None] + 1) * np.diff(e)[:, None]
        u = (y - y0) / EPS
        cn = np.where(np.abs(u) < 1, 2.0 / n * np.sin(
            n * np.arccos(np.clip(u, -1, 1))), 0) * np.exp(-1j * n * np.pi
                                                            / 2)
        c = (cn * gw[None]).sum(1) / 2.0
        vol = 2 * HP * (2 * EPS / NL) * 4.0
        wts = np.full(NL, 1.0 / vol)
        dSt = c[:, None] * (L.total - Rm.total)[None]
        dN0 = c[:, None, None] * (L.scatter - Rm.scatter)[None]
        dP = c[:, None, None] * ((chi_omega(L, w) - chi_omega(Rm, w))
                                 / cortex.KEFF)[None]
        F = a[f"pt_{name}"][:, :, ::-1, :, 0].real       # library order
        A = V_DET * adj[f"pt_{name}"][:, :, ::-1, :]
        for r in range(F.shape[0]):
            p0, p1 = pt_terms_complex(F[r], A[r], wts, dSt, dN0,
                                      np.zeros_like(dN0))
            out[r] += (p0, p1, pt_fission_complex(F[r], A[r], wts, dP))
    return out


def analyze(args):
    s = np.load(HERE.parent / "stage6" / "stage6c_abeille_static.npz")
    st = np.transpose(s["avg"], (0, 2, 1))
    st_sd = np.transpose(s["std"], (0, 2, 1))
    c = slice(68, 70)
    norm_b = st[1, c, c].mean()
    norm_b_rel = np.sqrt((st_sd[1, c, c]**2).sum()) / 4 / norm_b
    gate = []
    for n in args.harm:
        log(f"Stage 8c n={n}")
        A = {}
        for canc in args.cancellation:
            A[canc], a, norm = normalised_openmc(
                HERE / f"stage8c_openmc_n{n}_c{canc}{args.sfx}.npz")
        pb = HERE / f"stage8c_abeille_n{n}{args.sfx}.npz"
        if pb.exists():
            b = np.load(pb)
            B = np.transpose(b["vals"], (0, 1, 3, 2)) / norm_b
            log(f"  Abeille batches {b['batches'].tolist()}")
            for canc in args.cancellation:
                for g, name in ((0, "fast"), (1, "thermal")):
                    a_g = A[canc][:, g].reshape(A[canc].shape[0], -1)
                    b_g = B[:, g].reshape(B.shape[0], -1)
                    mod = np.abs(b_g.mean(0))
                    mask = mod > 0.1 * mod.max()
                    ok, info = two_sample_correlated(a_g, b_g, mask,
                                                     b_rel_sys=norm_b_rel)
                    rel = np.abs(a_g.mean(0)[mask]) / mod[mask] - 1
                    k, T2, p = info["hotelling"]
                    log(f"  c={canc} {name} vs Abeille: "
                        f"{'PASS' if ok else 'FAIL'} n={info['z'].size} "
                        f"bad={info['n_bad']} max|t|="
                        f"{np.abs(info['z']).max():.2f}; integral t "
                        f"Re={info['z_int'][0]:+.2f} "
                        f"Im={info['z_int'][1]:+.2f}; Hotelling k={k} "
                        f"p={p:.3f}; modulus OpenMC/Abeille - 1 median "
                        f"{np.median(rel):+.3f}")
                    gate.append((f"8c n={n} c={canc} {name} vs Abeille", ok))
        else:
            log(f"  (no Abeille result {pb.name})")
        if len(args.cancellation) == 2:
            c1, c2 = args.cancellation
            for g, name in ((0, "fast"), (1, "thermal")):
                x1 = A[c1][:, g].reshape(A[c1].shape[0], -1)
                x2 = A[c2][:, g].reshape(A[c2].shape[0], -1)
                mod = np.abs(x1.mean(0))
                mask = mod > 0.1 * mod.max()
                ok, info = two_sample_correlated(x1, x2, mask)
                log(f"  {name} cancellation {c1} vs {c2}: "
                    f"{'PASS' if ok else 'FAIL'} integral t "
                    f"Re={info['z_int'][0]:+.2f} Im={info['z_int'][1]:+.2f};"
                    f" Hotelling p={info['hotelling'][2]:.3f}")
                gate.append((f"8c n={n} {name} cancellation meshes", ok))
        pa = HERE / f"stage8c_adjoint_n{n}{args.sfx}.npz"
        if pa.exists():
            adj = np.load(pa)
            for canc in args.cancellation:
                a = np.load(HERE / f"stage8c_openmc_n{n}_c{canc}"
                            f"{args.sfx}.npz")
                dd = a["det"][:, 0, 1]               # thermal (E bin 0), noise
                P = contraction(a, adj, n, None)
                cv = P.sum(axis=1)
                ok, info = paired_correlated(dd[:, None], cv[:, None])
                log(f"  c={canc} detector vs contraction: "
                    f"{'PASS' if ok else 'FAIL'} paired t "
                    f"Re={info['t'][0]:+.2f} Im={info['t'][1]:+.2f}; "
                    f"forward {fmt_c(dd.mean())} contraction "
                    f"{fmt_c(cv.mean())} (P0 {fmt_c(P[:, 0].mean())}, P1 "
                    f"{fmt_c(P[:, 1].mean())}, fission "
                    f"{fmt_c(P[:, 2].mean())})")
                gate.append((f"8c n={n} c={canc} second estimator", ok))
    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 8c: {len(gate) - n_fail}/{len(gate)} gate checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    return n_fail


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["openmc", "adjoint", "abeille",
                                    "analyze"])
    p.add_argument("--harmonics", default="1,2")
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--tag", default="")
    p.add_argument("--particles", type=int, default=20000)
    p.add_argument("--batches", type=int, default=600)
    p.add_argument("--nskip", type=int, default=20)
    p.add_argument("--static-inactive", type=int, default=200)
    p.add_argument("--cancellation", default="170,340")
    p.add_argument("--particles-adj", type=int, default=20000)
    p.add_argument("--batches-adj", type=int, default=100)
    p.add_argument("--abeille-batches", type=int, default=600)
    p.add_argument("--abeille-ignored", type=int, default=100)
    p.add_argument("--abeille-particles", type=int, default=100000)
    p.add_argument("--max-minutes", type=float, default=1260.0)
    args = p.parse_args()
    args.sfx = f"_{args.tag}" if args.tag else ""
    args.harm = [int(x) for x in args.harmonics.split(",")]
    args.cancellation = [int(x) for x in args.cancellation.split(",")]
    n_fail = 0
    if args.mode == "openmc":
        run_openmc(args)
    elif args.mode == "adjoint":
        run_adjoint(args)
    elif args.mode == "abeille":
        run_abeille(args)
    else:
        n_fail = analyze(args)
    tag = args.mode + args.sfx
    with open(HERE / f"stage8c_log_{tag}.md", "w") as fh:
        fh.write("# Stage 8c log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
