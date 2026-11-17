"""Stage 4 gates: space dependence in a slab.

Slab |x| < a, vacuum at +-a, reflective transverse faces, isotropic
scattering.

4a  one group, 2a = 60 cm, Sigma_t=1, Sigma_s=0.9, nuSigma_f set so k ~ 1,
    v=2.2e5, beta=0.0065, lambda=0.0784; unit isotropic source in [14, 15];
    0.1, 1, 10, 100 Hz; Re/Im flux on a 1 cm mesh vs ref_sn1d (and
    ref_diff1d, recorded only).
4b  two groups (CORTEX fuel data in a slab sized to k ~ 1 with ref_sn1d),
    1 and 10 Hz, isotropic group-2 source in one 5 cm cell, vs ref_sn1d.
4c  point-kinetic limit for 4a at 1 Hz: G0_dyn from the OpenMC noise flux and
    the static adjoint of ref_sn1d, compared with G0 (2 %, 1 deg).

k_eff: the noise operator divides fission by the static k of the same
transport problem. OpenMC's eigenvalue for this slab is checked against the
converged S_N value within statistics, and the noise run uses the S_N value
(precision ~1e-7, far below an affordable MC sigma_k). Sensitivity is logged
by rerunning 4a at 0.1 Hz with k +- 2 sigma_k (OpenMC sigma_k).

Pass: STD on mesh cells where |dphi| > 10 % of its maximum; 4c amplitude
within 2 % and phase within 1 deg.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "reference"))

from build import box_source, make_library, slab_model, x_mesh  # noqa
from harness import fmt_info, fom, run_replicas, std_criterion, summarize  # noqa
from noise_common import Material  # noqa
from ref_diff1d import slab_noise  # noqa
from ref_g0 import g0  # noqa
from ref_inf import cortex_fuel  # noqa
from ref_sn1d import Slab, critical_nu_fission  # noqa

WORK = HERE / "work"
LOG = []
NUSF_4A = 0.1008719060  # critical at 60 cells / S16 (ref_sn1d self-test)


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def mat_4a(nusf=NUSF_4A):
    return Material.from_beta([1.0], [[0.9]], [1 / 2.2e5], [nusf], [1.0],
                              beta=[0.0065], lam=[0.0784])


def sn_reference(slab_regions, x0, G, freqs, src_region, spectrum,
                 n_angles=64, adjoint_static=False):
    sl = Slab(slab_regions, x0=x0, n_angles=n_angles)
    k, phi0 = sl.eigenvalue()
    q = sl.isotropic_source(src_region, spectrum)
    refs = {f: sl.solve(2 * np.pi * f, q) for f in freqs}
    return sl, k, refs


def bin_to_mesh(sl, values, edges):
    """Average cell values onto coarser bins with given edges (aligned)."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (sl.x > lo) & (sl.x < hi)
        out.append(values[m].mean(axis=0))
    return np.array(out)


def openmc_keff(lib, regions, x0, particles=200000, batches=300, inactive=50):
    model, _ = slab_model(regions, x0=x0, library=lib)
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.inactive = inactive
    lo, hi = x0, x0 + sum(w for w, _ in regions)
    model.settings.source = [box_source(lo, hi, 0, 1 if lib.name.startswith(
        "lib_4a") else 2, transverse=10.0)]
    d = WORK / ("keff_" + lib.stem)
    d.mkdir(parents=True, exist_ok=True)
    sp = model.run(cwd=d, threads=12, output=False)
    with openmc.StatePoint(sp) as s:
        return s.keff.n, s.keff.s


def noise_model(lib, regions, x0, G, f, eta, keff, src_region, src_group,
                n_mesh, particles, batches, noise_extra=None):
    model, edges = slab_model(regions, x0=x0, library=lib)
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.noise = {"frequency": f, "eta": eta, "keff": keff,
                            **(noise_extra or {})}
    model.settings.source = [box_source(src_region[0], src_region[1],
                                        src_group, G)]
    mesh = x_mesh(edges[0], edges[-1], n_mesh)
    filters = [openmc.MeshFilter(mesh),
               openmc.EnergyFilter(np.arange(G + 1.0)),
               openmc.ComplexPartFilter()]
    t = openmc.Tally(name="tl")
    t.filters = filters
    t.scores = ["flux"]
    model.tallies = openmc.Tallies([t])
    return model


def check(label, vals, ref, mask_frac=0.1):
    mean, sre, sim, corr = summarize(vals)
    mask = np.abs(ref) > mask_frac * np.abs(ref).max()
    ok, info = std_criterion(mean, sre, sim, ref, mask=mask)
    log(f"  {label}: {'PASS' if ok else 'FAIL'} {fmt_info(info)}")
    return ok, mean, sre, sim, info


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=20000)
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--only", default="4a,4b,4c,ksens")
    args = p.parse_args()
    only = args.only.split(",")
    WORK.mkdir(exist_ok=True)
    gate = []
    zs = []

    # ---------------------------------------------------------------- 4a
    freqs = [0.1, 1, 10, 100]
    regions = [(60.0, "m")]
    mat = mat_4a()
    lib = WORK / "lib_4a.h5"
    make_library({"m": mat}, lib)
    sl, k_sn, refs = sn_reference([(60.0, mat, 480)], -30.0, 1, freqs,
                                  (14, 15), [1.0])
    edges = np.arange(-30, 31, 1.0)
    log(f"Stage 4a: S_N (480 cells, S64) k = {k_sn:.8f}")
    k_mc, sk_mc = (None, None)
    if "4a" in only or "4c" in only or "ksens" in only:
        k_mc, sk_mc = openmc_keff(lib, regions, -30.0)
        z = (k_mc - k_sn) / sk_mc
        log(f"  OpenMC eigenvalue k = {k_mc:.6f} +- {sk_mc:.6f}; "
            f"(k_MC - k_SN)/sigma = {z:+.2f}")
        gate.append(("4a OpenMC k consistent with S_N k", abs(z) < 3))
    results_4a = {}
    if "4a" in only or "4c" in only:
        for f in freqs:
            eta = 1.0
            ref = bin_to_mesh(sl, refs[f], edges)[:, 0]
            out = run_replicas(lambda: noise_model(
                lib, regions, -30.0, 1, f, eta, k_sn, (14, 15), 0, 60,
                args.particles, args.batches), WORK / f"4a_{f}",
                n_rep=args.reps)
            vals = out["values"]["tl"].reshape(args.reps, -1)
            ok, mean, sre, sim, info = check(f"4a {f} Hz", vals, ref)
            zs.append(info["z"])
            results_4a[f] = (mean, sre, sim)
            gate.append((f"4a {f} Hz", ok))
            # diffusion, informational
            D = 1 / 3.0
            dd = slab_noise(0.5 * (edges[1:] + edges[:-1]), 2 * np.pi * f,
                            30 + 0.7104, D, 0.1, NUSF_4A, 1 / 2.2e5, [0.0065],
                            [0.0784], (14, 15), 1.0)
            m = np.abs(ref) > 0.1 * np.abs(ref).max()
            log(f"    diffusion vs S_N (informational): max |rel diff| "
                f"{np.abs(dd[m] / ref[m] - 1).max():.3f}; T="
                f"{out['times'].sum():.0f}s")
            np.savez(HERE / f"stage4a_{f}Hz.npz", mean=mean, s_re=sre,
                     s_im=sim, ref=ref, diff=dd, x=0.5 * (edges[1:] + edges[:-1]))

    if "ksens" in only and sk_mc is not None:
        f = 0.1
        ref = bin_to_mesh(sl, refs[f], edges)[:, 0]
        for sign in (-2, 2):
            kk = k_sn + sign * sk_mc
            out = run_replicas(lambda: noise_model(
                lib, regions, -30.0, 1, f, 1.0, kk, (14, 15), 0, 60,
                args.particles, args.batches), WORK / f"4a_ksens{sign}",
                n_rep=max(args.reps // 4, 3))
            mean = out["values"]["tl"].reshape(-1, 60).mean(0)
            m = np.abs(ref) > 0.1 * np.abs(ref).max()
            log(f"  k sensitivity 0.1 Hz, k = k_SN {sign:+d} sigma_k "
                f"({kk:.6f}): mean |dphi| ratio to reference "
                f"{np.mean(np.abs(mean[m]) / np.abs(ref[m])):.4f}, mean "
                f"phase shift {np.degrees(np.mean(np.angle(mean[m] / ref[m]))):+.2f} deg")

    # ---------------------------------------------------------------- 4c
    if "4c" in only and 1 in results_4a:
        mean, sre, sim = results_4a[1]
        # static adjoint (one group: = forward); normalised quantities
        sla, ka, phiA = Slab([(60.0, mat, 480)], x0=-30.0, n_angles=64), None, None
        ka, phiA = sla.eigenvalue(adjoint=True)
        phiA_b = bin_to_mesh(sla, phiA, edges)[:, 0]
        phi0_b = bin_to_mesh(sl, sl.phi0, edges)[:, 0]
        # delta rho: for an external source S, G0_dyn = <phi+ dphi/v>/<phi+ phi0/v>
        # relates to the source reactivity <phi+ S>/<phi+ nuSf phi0>.
        src = ((edges[:-1] >= 14) & (edges[1:] <= 15)).astype(float)
        drho = (phiA_b @ src) / (phiA_b @ (NUSF_4A * phi0_b))
        Lam = 1 / (2.2e5 * NUSF_4A) * 1.0
        G0 = g0(2 * np.pi * 1, [0.0065], [0.0784], (phiA_b @ phi0_b) / 2.2e5 /
                (phiA_b @ (NUSF_4A * phi0_b)))
        Gdyn = (phiA_b @ mean) / (phiA_b @ phi0_b) / drho
        rel = abs(abs(Gdyn) / abs(G0) - 1)
        dph = np.degrees(np.angle(Gdyn / G0))
        ok = rel < 0.02 and abs(dph) < 1.0
        log(f"Test 4c: G0_dyn = {Gdyn:.3f} (|.|={abs(Gdyn):.2f}, "
            f"{np.degrees(np.angle(Gdyn)):.2f} deg) vs G0 = {G0:.3f} "
            f"(|.|={abs(G0):.2f}, {np.degrees(np.angle(G0)):.2f} deg): "
            f"amp diff {rel:.3%}, phase diff {dph:+.2f} deg "
            f"{'PASS' if ok else 'FAIL'}")
        gate.append(("4c point-kinetic limit", ok))

    # ---------------------------------------------------------------- 4b
    if "4b" in only:
        fuel = cortex_fuel()
        # Size a bare fuel slab... fuel k_inf = 0.62, so add fuel/water:
        # use a homogenized fuel+water mixture scaled to k ~ 1 in a slab.
        log("Stage 4b: two-group slab")
        water = Material(total=[0.25411, 1.2182],
                         scatter=[[0.22519143, 0.028124], [0.0, 1.188884]],
                         inv_v=fuel.inv_v)
        # homogenized cell: fuel volume fraction from the CORTEX pin cell
        vf = 0.7314**2 / 1.26**2
        mix = Material(total=vf * fuel.total + (1 - vf) * water.total,
                       scatter=vf * fuel.scatter + (1 - vf) * water.scatter,
                       inv_v=fuel.inv_v,
                       nu_fission_p=vf * fuel.nu_fission_p,
                       nu_fission_d=vf * fuel.nu_fission_d,
                       chi_p=fuel.chi_p, chi_d=fuel.chi_d, lam=fuel.lam)

        def slab_k(width):
            s = Slab([(width, mix, int(width * 4))], x0=-width / 2, n_angles=16)
            return s.eigenvalue()[0]
        kinf = Slab([(10, mix, 10)], bc=("reflective", "reflective"),
                     n_angles=8).eigenvalue()[0]
        log(f"  homogenized fuel/water k_inf = {kinf:.5f}")
        # find width with k = 1 by bisection
        lo, hi = 20.0, 400.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if slab_k(mid) < 1:
                lo = mid
            else:
                hi = mid
        width = round(0.5 * (lo + hi) / 10) * 10  # multiple of 10 cm
        log(f"  critical width ~ {0.5*(lo+hi):.2f} cm; using {width:.0f} cm")
        lib2 = WORK / "lib_4b.h5"
        make_library({"mix": mix}, lib2)
        reg2 = [(width, "mix")]
        x0 = -width / 2
        src = (-10.0, -5.0)  # one 5 cm mesh cell
        sl2, k2, refs2 = sn_reference([(width, mix, int(width * 4))], x0, 2,
                                      [1, 10], src, [0.0, 1.0])
        log(f"  S_N k = {k2:.8f}")
        k_mc2, sk2 = openmc_keff(lib2, reg2, x0, particles=100000,
                                 batches=200)
        z = (k_mc2 - k2) / sk2
        log(f"  OpenMC k = {k_mc2:.6f} +- {sk2:.6f}; z = {z:+.2f}")
        gate.append(("4b OpenMC k consistent with S_N k", abs(z) < 3))
        nb = int(width / 5)
        e2 = np.linspace(x0, -x0, nb + 1)
        for f in [1, 10]:
            ref = 5.0 * bin_to_mesh(sl2, refs2[f], e2)[:, ::-1].reshape(-1)
            out = run_replicas(lambda: noise_model(
                lib2, reg2, x0, 2, f, 1.0, k2, src, 1, nb, args.particles,
                args.batches), WORK / f"4b_{f}", n_rep=args.reps)
            vals = out["values"]["tl"].reshape(args.reps, -1)
            ok, mean, sre, sim, info = check(f"4b {f} Hz", vals, ref)
            zs.append(info["z"])
            gate.append((f"4b {f} Hz", ok))
            np.savez(HERE / f"stage4b_{f}Hz.npz", mean=mean, s_re=sre,
                     s_im=sim, ref=ref, edges=e2)

    if zs:
        z = np.concatenate(zs)
        log(f"Pooled z: n={z.size} mean={z.mean():+.2f} std={z.std(ddof=1):.2f}")
    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 4 ({args.only}): {len(gate) - n_fail}/{len(gate)} gate "
        "checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_")
    with open(HERE / f"stage4_log_{tag}.md", "w") as fh:
        fh.write("# Stage 4 log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
