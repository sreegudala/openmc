"""Stage 8 gates on a slab: vibrating fuel-water interface (design doc
Stage 8, amendments 2026-11-22b and 2026-11-23).

Slab: the Stage 4b/7 mixture ("fuel", 100 cm, vacuum faces, reflective
sides) with a 4 cm water gap at [10, 14] cm. The fuel-water interface at
x0 = 10 cm vibrates as x0 + eps sin(omega_0 t), omega_0 = 2 pi (1 Hz);
fuel is on the negative side. Harmonic n runs at n Hz with the sampled
vibration source (perturbation 'vibration'), weight cancellation on a 1 cm
mesh, noise keff = the S_N k of the slab with its gap. The layer cells
[x0 - eps, x0] (fuel) and [x0, x0 + eps] (water) carry the source.

probe  Forward sigma of every gated quantity at eps = 0.2, n = 1 and 2
       (amendment 2026-11-22b: measure first, then size the gate runs).
8a     eps = 0.2, n = 1 and 2. Reference ref_sn1d with the source
       -dL psi_0 built from the same c_n on 0.01 cm cells across the layer
       (40 cells across 2 eps).
       Gated:
       - dphi/phi_0 on 1 cm bins, both groups (STD-C, bins above 10 % of
         the maximum);
       - the detector response, thermal flux in D = [20, 25] cm, as
         dd/phi_0(D) (STD);
       - the second estimator: forward dd against the adjoint contraction
         <chi, dS>, paired test as in 7c.
       The contraction uses chi from an adjoint run (transposed library,
       thermal source in D, R = 1, source strength V_D = 500 cm^3) and dS
       from tallies on 8 bins across the layer: the static moments of
       phi_0 (Y00, Y1m) from the forward run and Delta = fuel - water
       times the bin-averaged c_n, contracted with pt_lib.pt_terms
       (P0 + P1) plus the P0 fission term
       (adjoint_noise_lib.pt_fission_complex).
8b     n = 1 at eps = 0.05 and 0.1 (with 8a's 0.2): dd/phi_0(D)/eps
       consistent between eps (two-sample test against eps = 0.2), and
       each against S_N (STD). S_N shows n = 1 linear in eps to 1e-4; n = 2
       is O(eps^2) (int c_2 du = 0), so 8b uses n = 1 only.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for sub in ("common", "reference", "stage4"):
    sys.path.insert(0, str(ROOT / sub))

from adjoint_noise_lib import (make_adjoint_noise_library,  # noqa: E402
                               pt_fission_complex, pt_terms_complex)
from build import box_source, make_library, slab_model, x_mesh  # noqa: E402
from harness import (fmt_info, fmt_info_correlated, paired_correlated,  # noqa
                     run_replicas, std_correlated, std_criterion, summarize,
                     two_sample_correlated)
from noise_common import Material, chi_omega  # noqa: E402
from ref_inf import cortex_fuel, cortex_water  # noqa: E402
from ref_sn1d import Slab  # noqa: E402
import run_stage4 as s4  # noqa: E402

WORK = HERE / "work"
LOG = []
SAVED = {}

X0, GAP_END, X_LEFT, X_RIGHT = 10.0, 14.0, -50.0, 50.0
DET = (20.0, 25.0)
NB = 100                       # 1 cm mesh
I_DET = slice(70, 75)          # mesh bins of DET
NL = 8                         # bins across the layer
V_DET = 500.0                  # cm^3 (transverse 10 x 10)
F0 = 1.0                       # Hz, omega_0 / 2 pi
FAST, THERMAL = 0, 1
E_FAST, E_THERMAL = 1, 0
ADJ_SEED_OFFSET = 1000003


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def fmt_c(z):
    return f"{z.real:+.5e}{z.imag:+.5e}j"


def materials():
    water = cortex_water()
    return s4.mix_4b(), Material(water.total, water.scatter,
                                 cortex_fuel().inv_v)


def regions(eps):
    return [(50 + X0 - eps, "mix"), (eps, "mix"), (eps, "water"),
            (GAP_END - X0 - eps, "water"), (X_RIGHT - GAP_END, "mix")]


def c_n(n, u):
    u = np.asarray(u, float)
    s = np.where(np.abs(u) < 1,
                 2.0 / n * np.sin(n * np.arccos(np.clip(u, -1, 1))), 0.0)
    return s * np.exp(-1j * n * np.pi / 2)


def c_avg(n, lo, hi, eps, npts=16):
    """Average of c_n over [lo, hi] (arrays), Gauss-Legendre."""
    gx, gw = np.polynomial.legendre.leggauss(npts)
    lo, hi = np.asarray(lo)[:, None], np.asarray(hi)[:, None]
    x = lo + 0.5 * (gx[None, :] + 1) * (hi - lo)
    return (c_n(n, (x - X0) / eps) * gw[None, :]).sum(axis=1) / 2.0


# ------------------------------------------------------------------ S_N
class SNRef:
    def __init__(self, eps, n_angles=64):
        mix, wat = materials()
        self.mix, self.wat, self.eps = mix, wat, eps
        a, b = 50 + X0 - eps, GAP_END - X0 - eps
        # 0.01 cm cells across the layer (40 across 2 eps at eps = 0.2)
        nl = max(20, int(round(eps / 0.01)))
        self.sl = Slab([(a, mix, int(round(4 * a))), (eps, mix, nl),
                        (eps, wat, nl), (b, wat, int(round(4 * b))),
                        (X_RIGHT - GAP_END, mix, 4 * int(X_RIGHT - GAP_END))],
                       x0=X_LEFT, n_angles=n_angles)
        self.k, _ = self.sl.eigenvalue()
        self.psi0 = self.sl.static_angular()
        self.edges = np.linspace(X_LEFT, X_RIGHT, NB + 1)

    def source(self, n):
        sl = self.sl
        w = 2 * np.pi * n * F0
        c = c_avg(n, sl.edges[:-1], sl.edges[1:], self.eps)
        dT = self.mix.total - self.wat.total
        dS = self.mix.scatter - self.wat.scatter
        dP = chi_omega(self.mix, w) - chi_omega(self.wat, w)
        q = -c[:, None, None] * dT[None, None, :] * self.psi0
        iso = c[:, None] * (sl.phi0 @ dS
                            + np.einsum("gh,ih->ig", dP, sl.phi0) / sl.k)
        return q + 0.5 * iso[:, None, :]

    def ratios(self, n):
        """dphi/phi_0 on the 1 cm mesh ([x][E ascending]) and the detector
        ratio dd/phi_0(D, thermal)."""
        sl = self.sl
        d = sl.solve(2 * np.pi * n * F0, self.source(n))
        num = s4.bin_to_mesh(sl, d, self.edges)
        den = s4.bin_to_mesh(sl, sl.phi0, self.edges)
        field = (num / den)[:, ::-1].reshape(-1)
        D = (sl.x > DET[0]) & (sl.x < DET[1])
        det = (d[D, THERMAL] * sl.h[D]).sum() / (sl.phi0[D, THERMAL]
                                                  * sl.h[D]).sum()
        return field, det


# ---------------------------------------------------------------- models
def layer_mesh(eps):
    return x_mesh(X0 - eps, X0 + eps, NL)


def forward_model(lib, n, eps, k, args):
    model, edges = slab_model(regions(eps), x0=X_LEFT, library=lib)
    cells = list(model.geometry.root_universe.cells.values())
    mats = {m.name: m for m in model.materials}
    s = model.settings
    s.particles = args.particles
    s.batches = args.batches
    s.source = [box_source(X_LEFT, X_RIGHT, FAST, 2)]
    s.noise = {
        "frequency": n * F0, "eta": 1.0, "keff": k, "nskip": args.nskip,
        "static_inactive": args.static_inactive,
        "perturbations": [{"cells": cells[1:3], "vibration": {
            "axis": "x", "position": X0, "amplitude": eps, "harmonic": n,
            "negative_material": mats["mix"],
            "positive_material": mats["water"]}}]}
    if args.cancellation:
        s.noise["cancellation_mesh"] = x_mesh(X_LEFT, X_RIGHT, NB)
    energy = openmc.EnergyFilter(np.arange(3.0))
    t = openmc.Tally(name="tl")
    t.filters = [openmc.MeshFilter(x_mesh(X_LEFT, X_RIGHT, NB)), energy,
                 openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
    t.scores = ["flux"]
    t2 = openmc.Tally(name="pt")
    t2.filters = [openmc.MeshFilter(layer_mesh(eps)), energy,
                  openmc.SphericalHarmonicsFilter(1),
                  openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
    t2.scores = ["flux"]
    model.tallies = openmc.Tallies([t, t2])
    return model


def adjoint_model(lib_adj, n, eps, k, args):
    model, edges = slab_model(regions(eps), x0=X_LEFT, library=lib_adj)
    s = model.settings
    s.particles = args.particles_adj
    s.batches = args.batches
    s.source = [box_source(DET[0], DET[1], THERMAL, 2)]
    s.noise = {"frequency": n * F0, "eta": 1.0, "keff": k}
    energy = openmc.EnergyFilter(np.arange(3.0))
    t2 = openmc.Tally(name="pt")
    t2.filters = [openmc.MeshFilter(layer_mesh(eps)), energy,
                  openmc.SphericalHarmonicsFilter(1),
                  openmc.ComplexPartFilter()]
    t2.scores = ["flux"]
    model.tallies = openmc.Tallies([t2])
    return model


def read_gens(sp, d):
    import h5py
    from harness import last_statepoint
    with h5py.File(last_statepoint(d), "r") as f:
        g = f["noise"]
        return {k: g[k][()] for k in ("generations_per_batch",)
                if k in g}


def gens(out):
    g = [e["generations_per_batch"].mean() for e in out["extras"]
         if "generations_per_batch" in e]
    return f"{np.mean(g):.0f}" if g else "-"


def run_forward(lib, n, eps, k, args, label):
    out = run_replicas(lambda: forward_model(lib, n, eps, k, args),
                       WORK / label, n_rep=args.reps,
                       base_seed=args.base_seed, tally_names=("tl", "pt"),
                       extra=read_gens)
    R = args.reps
    tl = out["values"]["tl"].reshape(R, NB, 2, 2)      # [x][E][phase]
    pt = out["values"]["pt"].reshape(R, NL, 2, 4, 2)   # [b][E][Y][phase]
    phi0 = tl[..., 0].real
    with np.errstate(divide="ignore", invalid="ignore"):
        field = (tl[..., 1] / phi0).reshape(R, -1)
    dd = tl[:, I_DET, E_THERMAL, 1].sum(axis=1)
    phiD = phi0[:, I_DET, E_THERMAL].sum(axis=1)
    return {"field": field, "dd": dd, "phiD": phiD, "det": dd / phiD,
            "F": pt[:, :, ::-1, :, 0].real, "times": out["times"],
            "gens": gens(out), "tl": tl}


def rel_sigma(v):
    m = v.mean(axis=0)
    s = np.hypot(v.real.std(axis=0, ddof=1), v.imag.std(axis=0, ddof=1)) \
        / np.sqrt(v.shape[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        return s / np.abs(m)


def check_forward(tag, fw, field_ref, det_ref, gate, gated=True):
    mask = np.abs(field_ref) > 0.1 * np.abs(field_ref).max()
    ok1, info = std_correlated(fw["field"], field_ref, mask)
    rs = rel_sigma(fw["field"])[mask]
    log(f"  {tag} field dphi/phi0 vs S_N: {'PASS' if ok1 else 'FAIL'} "
        f"(STD-C) {fmt_info_correlated(info)}; rel sigma per bin median "
        f"{np.median(rs):.4f} max {rs.max():.4f}")
    m, sre, sim, _ = summarize(fw["det"][:, None])
    ok2, info2 = std_criterion(m, sre, sim, np.array([det_ref]))
    log(f"  {tag} detector dd/phi0(D) vs S_N: {'PASS' if ok2 else 'FAIL'} "
        f"{fmt_info(info2)}; MC {fmt_c(m[0])} S_N {fmt_c(det_ref)}; rel "
        f"sigma {rel_sigma(fw['det'][:, None])[0]:.4f}; T "
        f"{fw['times'].sum():.0f}s; generations/batch {fw['gens']}")
    if gated:
        gate.append((f"{tag} field", ok1))
        gate.append((f"{tag} detector", ok2))


def contraction(fw, adj, n, eps, k):
    """Per replica, <chi, dS> on the layer bins: P0 + P1 (pt_lib) + fission
    P0 (forward replica r with adjoint replica r)."""
    mix, wat = materials()
    w = 2 * np.pi * n * F0
    e = np.linspace(X0 - eps, X0 + eps, NL + 1)
    c = c_avg(n, e[:-1], e[1:], eps)                     # (NL,)
    vol = (2 * eps / NL) * 100.0
    wts = np.full(NL, 1.0 / vol)
    dT = mix.total - wat.total
    dS = mix.scatter - wat.scatter
    dP = (chi_omega(mix, w) - chi_omega(wat, w)) / k     # [g_out, g_in]
    dSt = c[:, None] * dT[None, :]
    dN0 = c[:, None, None] * dS[None, :, :]
    dN1 = np.zeros_like(dN0)
    dPv = c[:, None, None] * dP[None, :, :]
    out = []
    for r in range(fw["F"].shape[0]):
        F, A = fw["F"][r], adj[r]
        p0, p1 = pt_terms_complex(F, A, wts, dSt, dN0, dN1)
        pf = pt_fission_complex(F, A, wts, dPv)
        out.append((p0, p1, pf))
    return np.array(out)                                 # (R, 3)


def run_adjoint(lib_adj, n, eps, k, args, label):
    out = run_replicas(lambda: adjoint_model(lib_adj, n, eps, k, args),
                       WORK / label, n_rep=args.reps,
                       base_seed=args.base_seed + ADJ_SEED_OFFSET,
                       tally_names=("pt",), extra=read_gens)
    A = out["values"]["pt"].reshape(args.reps, NL, 2, 4)[:, :, ::-1, :]
    return V_DET * A, out["times"], gens(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--only", default="8a,8b")
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=60000,
                   help="static particles per generation (forward)")
    p.add_argument("--particles-adj", type=int, default=20000)
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--nskip", type=int, default=20)
    p.add_argument("--static-inactive", type=int, default=300)
    p.add_argument("--no-cancellation", dest="cancellation",
                   action="store_false")
    p.add_argument("--harmonics", default="1,2")
    p.add_argument("--eps-8b", default="0.05,0.1")
    p.add_argument("--base-seed", type=int, default=1)
    p.add_argument("--tag", default="")
    args = p.parse_args()
    only = args.only.split(",")
    global WORK
    sfx = f"_{args.tag}" if args.tag else ""
    if args.tag:
        WORK = HERE / f"work_{args.tag}"
    WORK.mkdir(parents=True, exist_ok=True)
    gate = []
    mix, wat = materials()
    lib = WORK / "lib_8.h5"
    lib_adj = WORK / "lib_8_adjoint.h5"
    make_library({"mix": mix, "water": wat}, lib)
    make_adjoint_noise_library(lib, lib_adj)
    harmonics = [int(x) for x in args.harmonics.split(",")]
    eps0 = 0.2
    ref = SNRef(eps0)
    log(f"Slab: 4b mixture with water gap [{X0}, {GAP_END}] cm, interface "
        f"x0 = {X0} vibrating along x (fuel negative side); S_N k = "
        f"{ref.k:.8f}; detector {DET} thermal; cancellation "
        f"{'1 cm' if args.cancellation else 'off'}")
    fwd02 = {}

    if "probe" in only or "8a" in only:
        for n in harmonics:
            field_ref, det_ref = ref.ratios(n)
            fw = run_forward(lib, n, eps0, ref.k, args, f"fwd_n{n}_e0.2")
            fwd02[n] = fw
            tag = f"{'probe' if 'probe' in only else '8a'} n={n}"
            check_forward(tag, fw, field_ref, det_ref, gate,
                          gated="8a" in only)
            SAVED[f"fwd_n{n}_e0.2_field"] = fw["field"]
            SAVED[f"fwd_n{n}_e0.2_dd"] = fw["dd"]
            SAVED[f"fwd_n{n}_e0.2_phiD"] = fw["phiD"]
            SAVED[f"fwd_n{n}_e0.2_times"] = fw["times"]
            SAVED[f"ref_n{n}_e0.2_field"] = field_ref
            SAVED[f"ref_n{n}_e0.2_det"] = det_ref
            if "8a" not in only:
                continue
            A, ta, ga = run_adjoint(lib_adj, n, eps0, ref.k, args,
                                    f"adj_n{n}_e0.2")
            P = contraction(fw, A, n, eps0, ref.k)
            cval = P.sum(axis=1)
            ok, info = paired_correlated(fw["dd"][:, None], cval[:, None])
            md, mc = fw["dd"].mean(), cval.mean()
            log(f"  8a n={n} second estimator: {'PASS' if ok else 'FAIL'} "
                f"paired t Re={info['t'][0]:+.2f} Im={info['t'][1]:+.2f}, "
                f"Hotelling p={info['hotelling'][2]:.3f}; forward dd "
                f"{fmt_c(md)} contraction {fmt_c(mc)} (P0 "
                f"{fmt_c(P[:, 0].mean())}, P1 {fmt_c(P[:, 1].mean())}, "
                f"fission {fmt_c(P[:, 2].mean())}); |fwd/contr|-1 = "
                f"{abs(md) / abs(mc) - 1:+.4f}; rel sigma fwd "
                f"{rel_sigma(fw['dd'][:, None])[0]:.4f} contr "
                f"{rel_sigma(cval[:, None])[0]:.4f}; adjoint T "
                f"{ta.sum():.0f}s, generations/batch {ga}")
            # contraction against S_N, through the same phi_0(D)
            m, sre, sim, _ = summarize((cval / fw["phiD"])[:, None])
            _, inf2 = std_criterion(m, sre, sim, np.array([det_ref]))
            log(f"    contraction / phi0(D) vs S_N (logged): "
                f"{fmt_info(inf2)}")
            gate.append((f"8a n={n} second estimator", ok))
            SAVED[f"adj_n{n}_e0.2_A"] = A
            SAVED[f"contr_n{n}_e0.2"] = P

    if "8b" in only:
        log("Test 8b: n = 1 response / eps")
        if 1 not in fwd02:
            field_ref, det_ref = ref.ratios(1)
            fwd02[1] = run_forward(lib, 1, eps0, ref.k, args, "fwd_n1_e0.2")
        base = fwd02[1]["det"] / eps0
        for eps in [float(x) for x in args.eps_8b.split(",")]:
            r = SNRef(eps)
            field_ref, det_ref = r.ratios(1)
            fw = run_forward(lib, 1, eps, r.k, args, f"fwd_n1_e{eps:g}")
            check_forward(f"8b eps={eps:g}", fw, field_ref, det_ref, gate,
                          gated=False)
            m, sre, sim, _ = summarize(fw["det"][:, None])
            ok1, i1 = std_criterion(m, sre, sim, np.array([det_ref]))
            ok2, i2 = two_sample_correlated((fw["det"] / eps)[:, None],
                                            base[:, None])
            log(f"  8b eps={eps:g}: detector vs S_N {'PASS' if ok1 else 'FAIL'}"
                f" ({fmt_info(i1)}); (dd/phi0)/eps vs eps=0.2: "
                f"{'PASS' if ok2 else 'FAIL'} Welch t Re={i2['z'][0]:+.2f} "
                f"Im={i2['z'][1]:+.2f}; ratio "
                f"{fmt_c((fw['det'] / eps).mean() / base.mean())}; S_N "
                f"ratio {fmt_c(det_ref / eps / (ref.ratios(1)[1] / eps0))}")
            gate.append((f"8b eps={eps:g} vs S_N", ok1))
            gate.append((f"8b eps={eps:g} linear", ok2))
            SAVED[f"fwd_n1_e{eps:g}_det"] = fw["det"]
            SAVED[f"fwd_n1_e{eps:g}_times"] = fw["times"]

    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 8 ({args.only}): {len(gate) - n_fail}/{len(gate)} gate "
        "checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_") + sfx
    with open(HERE / f"stage8_log_{tag}.md", "w") as fh:
        fh.write("# Stage 8 log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    np.savez(HERE / f"stage8_vals_{tag}.npz", **SAVED)
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
