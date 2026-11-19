"""Stage 6c: CORTEX 2-D UOX assembly, cross-section oscillation at 1 Hz,
OpenMC against Abeille (design doc Section 3.4 and Stage 6).

Modes (run separately; each writes an npz in this directory):

  abeille         R independent Abeille noise runs of noise_oscillation.yaml
                  (seeds differ) -> stage6c_abeille{tag}.npz
  abeille-static  one Abeille k-eigenvalue run with a flux tally on the same
                  138 x 138 mesh, for the static normalisation
                  -> stage6c_abeille_static.npz
  openmc          R OpenMC replicas -> stage6c_openmc{tag}.npz
  analyze         gate: stage6c_log{tag}.md

Both noise fields are divided by their own static thermal flux at the
assembly centre (mean of the four central bins of the 138 x 138 mesh; the
centre is a water hole). OpenMC gets phi_0 from the static sampling
generation (NoisePhaseFilter); Abeille from the k-eigenvalue run. Both are
per particle of one static generation. OpenMC is given Abeille's complex
amplitudes (-i pi epsilon_x), so no convention factor is needed.

Gate (user decision 2026-11-19): two-sample STD-C (harness
two_sample_correlated) on cells above 10 % of the maximum modulus, per
group; modulus within 2 % for the assembly-integrated noise per group and
the perturbed-pin noise per group. Per-cell modulus differences are
logged.

Abeille is GPL-3.0 and lives outside the repo ($ABEILLE_DIR, built by
tacc/build_abeille.sh); its input is read from its own input_files.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE))

import harness  # noqa
from harness import run_replicas, two_sample_correlated  # noqa
import cortex  # noqa

LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def abeille_dir():
    d = os.environ.get("ABEILLE_DIR")
    if d is None:
        src = os.environ.get("OPENMC_DEV_SRC")
        d = (Path(src).parent / "abeille" if src
             else Path.home() / "dev/external/abeille")
    return Path(d)


def abeille_exe():
    """bin/abeille (tacc/build_abeille.sh), else build/abeille (local)."""
    for sub in ("bin", "build"):
        exe = abeille_dir() / sub / "abeille"
        if exe.exists():
            return exe
    raise FileNotFoundError(f"no Abeille executable under {abeille_dir()}")


def abeille_geometry():
    """Materials, surfaces, cells and universes of the shipped
    noise_oscillation.yaml (everything before its tally filters)."""
    text = (abeille_dir() / "input_files" /
            "noise_oscillation.yaml").read_text()
    return text[:text.index("tally-filters:")]


TALLIES = """tally-filters:
  position-filters:
    - id: 2
      type: regular-cartesian-mesh
      low: [-10.79, -10.79, -2.00]
      high: [10.79, 10.79, 2.00]
      shape: [138, 138, 1]
  energy-filters:
    - id: 1
      energy-bounds: [0., 1., 2.]

tallies:
{tallies}
settings:
  ngroups: 2
  energy-mode: multi-group
  energy-bounds: [0., 1., 2.]
  seed: {seed}
  max-run-time: {max_minutes}

"""

FLUX_TALLY = """  - name: {name}
    quantity: {name}
    estimator: track-length
    position-filter: 2
    energy-filter: 1
"""

SOURCE = """  sources:
    - spatial: {type: box, low: [-10.71, -10.71, -2.00], hi: [10.71, 10.71, 2.00]}
      direction: {type: isotropic}
      energy: {type: mono-energetic, energy: 0.5}
      weight: 1.
"""

# Run controls of the shipped noise input, with batches, inactive
# generations and particles as parameters
NOISE_SIM = """simulation:
  mode: noise
  nparticles: {particles}
  nbatches: {nbatches}
  nignored: {nignored}
  nskip: 3
  noise-cancellation: true
  inner-generations: true
  normalize-noise-source: true
  noise-angular-frequency: 6.28318531
  keff: 0.99916

  cancelator:
    type: approximate
    shape: [170, 170, 170]
    low: [-10.71, -10.71, -2.00]
    hi: [10.71, 10.71, 2.00]

""" + SOURCE + """
  noise-sources:
    - type: square-oscillation
      low: [2.1543, -2.8857, -2.00]
      hi: [2.8857, -2.1543, 2.00]
      angular-frequency: 6.28318531
      epsilon-total:
        energy: [0., 2.]
        value: [0.041, 0.041]
      epsilon-fission:
        energy: [0., 2.]
        value: [0.021, 0.021]
      epsilon-scatter:
        energy: [0., 2.]
        value: [0.034, 0.034]
"""

STATIC_SIM = """simulation:
  mode: k-eigenvalue
  nparticles: {particles}
  ngenerations: {ngen}
  nignored: 50
""" + SOURCE


def abeille_input(nbatches, nignored, seed, max_minutes, particles):
    """Abeille noise input: shipped geometry and noise source, 138 x 138
    real/imag flux tallies, given run controls."""
    tallies = "".join(FLUX_TALLY.format(name=n)
                      for n in ("real-flux", "imag-flux"))
    return (abeille_geometry()
            + TALLIES.format(tallies=tallies, seed=seed,
                             max_minutes=max_minutes)
            + NOISE_SIM.replace("{particles}", str(particles))
            .replace("{nbatches}", str(nbatches))
            .replace("{nignored}", str(nignored)))


def run_abeille(args):
    exe = abeille_exe()
    work = HERE / f"work_6c_abeille{args.sfx}"
    work.mkdir(exist_ok=True)
    dirs = []
    for r in range(args.reps):
        d = work / f"rep{r:02d}"
        d.mkdir(exist_ok=True)
        doc = abeille_input(args.abeille_batches, args.abeille_ignored,
                            19073486328125 + 2 * 7919 * (r + args.base_seed),
                            args.max_minutes, args.abeille_particles)
        (d / "input.yaml").write_text(doc)
        dirs.append(d)

    def one(d):
        t0 = time.time()
        with open(d / "abeille.log", "w") as fh:
            subprocess.run([str(exe), "--input", "input.yaml", "--threads",
                            str(harness.THREADS), "--output", "out.h5"],
                           cwd=d, stdout=fh, stderr=subprocess.STDOUT,
                           check=True)
        return time.time() - t0
    with ThreadPoolExecutor(max_workers=min(harness.JOBS, args.reps)) as ex:
        times = list(ex.map(one, dirs))
    vals, nb = [], []
    for d in dirs:
        with h5py.File(d / "out.h5", "r") as f:
            r = f["results"]
            vals.append(r["real-flux/avg"][()] + 1j * r["imag-flux/avg"][()])
        log_text = (d / "abeille.log").read_text()
        nb.append(log_text.count("Simulated noise batch"))
    out = HERE / f"stage6c_abeille{args.sfx}.npz"
    np.savez(out, vals=np.array(vals), batches=np.array(nb),
             times=np.array(times))
    log(f"Abeille: {args.reps} replicas, batches {nb}, times "
        f"{np.round(times)} s -> {out.name}")


def run_abeille_static(args):
    exe = abeille_exe()
    d = HERE / "work_6c_abeille_static"
    d.mkdir(exist_ok=True)
    text = (abeille_geometry()
            + TALLIES.format(tallies=FLUX_TALLY.format(name="flux"),
                             seed=19073486328125, max_minutes=1.0e6)
            + STATIC_SIM.replace("{particles}", str(args.abeille_particles))
            .replace("{ngen}", str(args.static_generations + 50)))
    (d / "input.yaml").write_text(text)
    with open(d / "abeille.log", "w") as fh:
        subprocess.run([str(exe), "--input", "input.yaml", "--threads",
                        str(args.static_threads), "--output", "out.h5"],
                       cwd=d, stdout=fh, stderr=subprocess.STDOUT, check=True)
    save_abeille_static(d / "out.h5")


def save_abeille_static(path):
    with h5py.File(path, "r") as f:
        r = f["results"]
        avg, std = r["flux/avg"][()], r["flux/std"][()]
        k = (r.attrs["kcol-avg"], r.attrs["kcol-std"]) if "kcol-avg" in \
            r.attrs else (np.nan, np.nan)
    out = HERE / "stage6c_abeille_static.npz"
    np.savez(out, avg=avg, std=std, k=np.array(k))
    log(f"Abeille static: k = {k[0]:.5f} +- {k[1]:.5f} -> {out.name}")


def run_openmc(args):
    work = HERE / f"work_6c_openmc{args.sfx}"
    work.mkdir(exist_ok=True)
    lib = (work / "lib_cortex.h5").resolve()
    cortex.library(lib)
    make = lambda: cortex.model(lib, args.particles, args.batches,
                                nskip=args.nskip,
                                static_inactive=args.static_inactive,
                                cancellation=args.cancellation)
    out = run_replicas(make, work / "reps", n_rep=args.reps,
                       tally_names=("flux",), base_seed=args.base_seed)
    v = out["values"]["flux"].reshape(args.reps, 138, 138, 2, 2)
    path = HERE / f"stage6c_openmc{args.sfx}.npz"
    np.savez(path, noise=v[..., 1], static=v[..., 0].real,
             times=out["times"])
    log(f"OpenMC: {args.reps} replicas x {args.batches} batches, T="
        f"{out['times'].sum():.0f}s -> {path.name}")


def analyze(args):
    a = np.load(HERE / f"stage6c_openmc{args.sfx}.npz")
    b = np.load(HERE / f"stage6c_abeille{args.sfx}.npz")
    s = np.load(HERE / "stage6c_abeille_static.npz")
    c = slice(68, 70)   # four central bins
    omc = a["noise"]
    omc0 = a["static"]
    # OpenMC tally bins: mesh index i = ix + 138 iy (x fastest); reshape
    # (138, 138) in C order gives [iy, ix]. Energy bin 0 = thermal.
    omc = np.stack([omc[..., 1], omc[..., 0]], axis=1)  # (R, g, iy, ix)
    omc0 = np.stack([omc0[..., 1], omc0[..., 0]], axis=1)
    norm_a = omc0[:, 1, c, c].mean(axis=(1, 2))           # thermal centre
    A = omc / norm_a[:, None, None, None]
    # Abeille: (R, g, ix, iy) (tally shape [E, x, y]); g 0 = fast
    ab = np.transpose(b["vals"], (0, 1, 3, 2))             # (R, g, iy, ix)
    st = np.transpose(s["avg"], (0, 2, 1))
    st_sd = np.transpose(s["std"], (0, 2, 1))
    norm_b = st[1, c, c].mean()
    norm_b_rel = np.sqrt((st_sd[1, c, c]**2).sum()) / 4 / norm_b
    B = ab / norm_b
    log(f"Stage 6c: OpenMC {A.shape[0]} replicas, Abeille {B.shape[0]} "
        f"replicas (batches {b['batches'].tolist()}); Abeille static "
        f"normalisation rel sigma {norm_b_rel:.2e}; Abeille static k "
        f"{s['k'][0]:.5f} +- {s['k'][1]:.5f}")
    # static shape check (informational)
    sa = (omc0 / norm_a[:, None, None, None]).mean(0)
    sb = st / norm_b
    m = sb > 0.1 * sb.max()
    log(f"  static phi_0 / phi_0(centre, thermal): max |OpenMC / Abeille - "
        f"1| = {np.abs(sa[m] / sb[m] - 1).max():.3e} (fast and thermal, "
        f"OpenMC replica mean)")

    gate = []
    for g, name in ((0, "fast"), (1, "thermal")):
        a_g = A[:, g].reshape(A.shape[0], -1)
        b_g = B[:, g].reshape(B.shape[0], -1)
        mod = np.abs(b_g.mean(0))
        mask = mod > 0.1 * mod.max()
        ok, info = two_sample_correlated(a_g, b_g, mask, b_rel_sys=norm_b_rel)
        k, T2, p = info["hotelling"]
        rel = np.abs(a_g.mean(0)[mask]) / mod[mask] - 1
        log(f"  {name}: {'PASS' if ok else 'FAIL'} two-sample STD-C: "
            f"n={info['z'].size} bad={info['n_bad']} max|z|="
            f"{np.abs(info['z']).max():.2f}, z hist mean "
            f"{info['z'].mean():+.2f} std {info['z'].std(ddof=1):.2f}; "
            f"integral z Re={info['z_int'][0]:+.2f} Im={info['z_int'][1]:+.2f}"
            f"; Hotelling k={k} T2={T2:.1f} p={p:.3f}; per-cell modulus "
            f"OpenMC/Abeille - 1: median {np.median(rel):+.3f}, 5-95 % "
            f"[{np.percentile(rel, 5):+.3f}, {np.percentile(rel, 95):+.3f}]")
        gate.append((f"6c {name} STD-C", ok))
        # coarse quantities: assembly integral and perturbed pin
        iy = slice(int(np.ceil((-2.8857 + 10.79) / (21.58 / 138))),
                   int(np.floor((-2.1543 + 10.79) / (21.58 / 138))))
        ix = slice(int(np.ceil((2.1543 + 10.79) / (21.58 / 138))),
                   int(np.floor((2.8857 + 10.79) / (21.58 / 138))))
        for qname, fa, fb in (
                ("assembly integral", A[:, g].sum(axis=(1, 2)),
                 B[:, g].sum(axis=(1, 2))),
                ("perturbed pin", A[:, g, iy, ix].sum(axis=(1, 2)),
                 B[:, g, iy, ix].sum(axis=(1, 2)))):
            ma, mb = fa.mean(), fb.mean()
            sa_ = np.hypot(fa.real.std(ddof=1), fa.imag.std(ddof=1)) \
                / np.sqrt(fa.size) / abs(ma)
            sb_ = np.hypot(np.hypot(fb.real.std(ddof=1), fb.imag.std(ddof=1))
                           / np.sqrt(fb.size) / abs(mb), norm_b_rel)
            r = abs(ma) / abs(mb) - 1
            dph = np.degrees(np.angle(ma / mb))
            ok2 = abs(r) < 0.02
            log(f"    {name} {qname}: OpenMC {ma:.5e} Abeille {mb:.5e}; "
                f"|modulus ratio - 1| = {abs(r):.4f} (sigma OpenMC {sa_:.4f}, "
                f"Abeille {sb_:.4f}), phase diff {dph:+.2f} deg -> "
                f"{'PASS' if ok2 else 'FAIL'} (2 %)")
            gate.append((f"6c {name} {qname} modulus 2 %", ok2))
    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 6c: {len(gate) - n_fail}/{len(gate)} gate checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    with open(HERE / f"stage6c_log{args.sfx}.md", "w") as fh:
        fh.write("# Stage 6c log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    return n_fail


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["abeille", "abeille-static", "openmc",
                                    "analyze"])
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--tag", default="")
    # OpenMC
    p.add_argument("--particles", type=int, default=20000)
    p.add_argument("--batches", type=int, default=600)
    p.add_argument("--nskip", type=int, default=20)
    p.add_argument("--static-inactive", type=int, default=200)
    p.add_argument("--cancellation", type=int, default=68,
                   help="n for an n x n cancellation mesh (0: off)")
    # Abeille
    p.add_argument("--abeille-batches", type=int, default=600)
    p.add_argument("--abeille-ignored", type=int, default=100)
    p.add_argument("--abeille-particles", type=int, default=100000)
    p.add_argument("--max-minutes", type=float, default=1260.0,
                   help="Abeille max-run-time per replica (stops cleanly)")
    p.add_argument("--static-generations", type=int, default=1050)
    p.add_argument("--static-threads", type=int, default=12)
    args = p.parse_args()
    args.sfx = f"_{args.tag}" if args.tag else ""
    if args.mode == "abeille":
        run_abeille(args)
    elif args.mode == "abeille-static":
        run_abeille_static(args)
    elif args.mode == "openmc":
        run_openmc(args)
    else:
        raise SystemExit(1 if analyze(args) else 0)


if __name__ == "__main__":
    main()
