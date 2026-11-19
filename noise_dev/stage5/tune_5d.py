"""Stage 5d tuning study: what controls the FOM gain of weight cancellation.

Slab 4a (as in 5b/5d), cancellation off against on for a set of settings.
The FOM ratio is estimated from batch-to-batch variances pooled over
replicas: batches are independent (each is a full fission-generation game
from fresh source particles), so R replicas x B batches give R(B - 1)
degrees of freedom, against R - 1 for replica-to-replica variances. A
cross-check of the two variance estimates is logged for every setting.

For each setting and frequency the log gives the median over mesh bins of
FOM_on / FOM_off for Re and Im with a 95 % bootstrap interval (resampling
replicas), the time ratio, and cancellation diagnostics written by OpenMC
(sum|w| after / before averaging, sites per occupied bin, fraction of sites
removed by the bank roulette, generations per batch).
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "stage4"))

from build import make_library  # noqa
from harness import last_statepoint, run_replicas  # noqa
import run_stage4 as s4  # noqa

LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def slab_mesh(h, a=30.0):
    n = int(round(2 * a / h))
    m = openmc.RegularMesh()
    m.lower_left = (-a, -5, -5)
    m.upper_right = (a, 5, 5)
    m.dimension = (n, 1, 1)
    return m


# name -> (cancellation settings or None, particles multiplier)
CONFIGS = {
    "off": (None, 1),
    "h0.5": ({"h": 0.5}, 1),
    "h1": ({"h": 1.0}, 1),
    "h2": ({"h": 2.0}, 1),
    "h4": ({"h": 4.0}, 1),
    "h1_noroul": ({"h": 1.0, "cancellation_roulette": False}, 1),
    "off_x4": (None, 4),
    "h1_x4": ({"h": 1.0}, 4),
}


def read_extra(sp, d):
    """Per-replica batch statistics of the tally and noise diagnostics."""
    st = openmc.noise.complex_tally_stats(sp.get_tally(name="tl"))
    n = sp.get_tally(name="tl").num_realizations
    out = {"var_re": (st["re_std"].squeeze() ** 2) * n,
           "var_im": (st["im_std"].squeeze() ** 2) * n, "n": n}
    with h5py.File(last_statepoint(d), "r") as f:
        g = f["noise"]
        out["gens"] = g["generations_per_batch"][()].mean()
        for k in ("cancel_n_rendezvous", "cancel_n_sites", "cancel_n_occupied",
                  "cancel_abs_before", "cancel_abs_after", "cancel_n_removed"):
            out[k] = g[k][()] if k in g else np.nan
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=60)
    p.add_argument("--particles", type=int, default=20000)
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--freqs", default="1.0,0.1")
    p.add_argument("--configs", default="off,h0.5,h1,h2,h4,h1_noroul",
                   help="comma-separated names from CONFIGS; the first must "
                   "be 'off' (the reference)")
    p.add_argument("--configs-1hz-extra", default="off_x4,h1_x4",
                   help="additional configs run at 1 Hz only")
    p.add_argument("--base-seed", type=int, default=300001)
    p.add_argument("--tag", default="")
    args = p.parse_args()
    sfx = f"_{args.tag}" if args.tag else ""
    work = HERE / f"work_tune5d{sfx}"
    work.mkdir(exist_ok=True)
    lib = work / "lib_4a.h5"
    make_library({"m": s4.mat_4a()}, lib)
    regions = [(60.0, "m")]
    k = s4.Slab([(60.0, s4.mat_4a(), 480)], x0=-30.0,
                n_angles=64).eigenvalue()[0]
    saved = {}
    log(f"Stage 5d tuning: slab 4a, {args.reps} replicas x {args.batches} "
        f"batches, {args.particles} particles (x multiplier), seeds "
        f"{args.base_seed} + 7919 r")

    for f in [float(x) for x in args.freqs.split(",")]:
        names = args.configs.split(",")
        if f == 1.0 and args.configs_1hz_extra:
            names += args.configs_1hz_extra.split(",")
        log(f"\n{f} Hz")
        res = {}
        for name in names:
            canc, mult = CONFIGS[name]
            extra = None
            if canc is not None:
                extra = {"cancellation_mesh": slab_mesh(canc["h"])}
                extra.update({kk: v for kk, v in canc.items() if kk != "h"})
            mk = lambda: s4.noise_model(
                lib, regions, -30.0, 1, f, 1.0, k, (14, 15), 0, 60,
                args.particles * mult, args.batches, noise_extra=extra)
            out = run_replicas(mk, work / f"{name}_{f}", n_rep=args.reps,
                               base_seed=args.base_seed, extra=read_extra)
            ex = out["extras"]
            r = {"vals": out["values"]["tl"].reshape(args.reps, -1),
                 "t": out["times"],
                 "var_re": np.array([e["var_re"] for e in ex]),
                 "var_im": np.array([e["var_im"] for e in ex]),
                 "n": np.array([e["n"] for e in ex]),
                 "gens": np.array([e["gens"] for e in ex])}
            for kk in ("cancel_n_rendezvous", "cancel_n_sites",
                       "cancel_n_occupied", "cancel_abs_before",
                       "cancel_abs_after", "cancel_n_removed"):
                r[kk] = np.array([e[kk] for e in ex])
            res[name] = r
            for kk, v in r.items():
                saved[f"{name}_{f}_{kk}"] = v
            # replica-based vs batch-based variance of a replica mean
            vals = r["vals"]
            m = vals.mean(0)
            mask = np.abs(m) > 0.1 * np.abs(m).max()
            chk_re = np.median((vals.real.var(0, ddof=1)
                                / (r["var_re"].mean(0) / r["n"].mean()))[mask])
            chk_im = np.median((vals.imag.var(0, ddof=1)
                                / (r["var_im"].mean(0) / r["n"].mean()))[mask])
            diag = ""
            if canc is not None:
                diag = (f"; |sum w| after/before averaging "
                        f"{r['cancel_abs_after'].sum() / r['cancel_abs_before'].sum():.3f}, "
                        f"sites per occupied bin "
                        f"{r['cancel_n_sites'].sum() / r['cancel_n_occupied'].sum():.0f}, "
                        f"removed by roulette "
                        f"{r['cancel_n_removed'].sum() / r['cancel_n_sites'].sum():.2e}")
            log(f"  {name}: T={r['t'].sum():.0f}s gens/batch "
                f"{r['gens'].mean():.0f}; var check replica/batch Re "
                f"{chk_re:.2f} Im {chk_im:.2f}{diag}")

        # FOM ratios against the matching 'off' run
        rng = np.random.default_rng(12345)
        R = args.reps
        for name in names:
            if name.startswith("off"):
                continue
            ref = "off_x4" if name.endswith("_x4") else "off"
            a, b = res[ref], res[name]
            m = a["vals"].mean(0)
            mask = np.abs(m) > 0.1 * np.abs(m).max()

            def ratio(idx):
                tr = a["t"][idx].sum() / b["t"][idx].sum()
                return (np.median((a["var_re"][idx].mean(0)
                                   / b["var_re"][idx].mean(0))[mask]) * tr,
                        np.median((a["var_im"][idx].mean(0)
                                   / b["var_im"][idx].mean(0))[mask]) * tr,
                        tr)
            est = ratio(np.arange(R))
            bs = np.array([ratio(rng.integers(0, R, R)) for _ in range(2000)])
            lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
            log(f"  FOM {name}/{ref}: Re {est[0]:.2f} [{lo[0]:.2f}, "
                f"{hi[0]:.2f}]  Im {est[1]:.2f} [{lo[1]:.2f}, {hi[1]:.2f}]  "
                f"time off/on {est[2]:.3f}")

    np.savez(HERE / f"tune5d_vals{sfx}.npz", **saved)
    with open(HERE / f"tune5d_log{sfx}.md", "w") as fh:
        fh.write("# Stage 5d tuning log\n\n```\n" + "\n".join(LOG) + "\n```\n")


if __name__ == "__main__":
    main()
