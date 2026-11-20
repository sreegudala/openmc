"""Stage 7 gates: adjoint noise and reciprocity (design doc Stage 7,
amendment 2026-11-20; amendment 2026-11-20b for 7e).

Bilinear convention: the adjoint chi solves L(omega)^T chi = R at the same
+omega and is computed by the unchanged noise solver on the transposed
library of adjoint_noise_lib (N1 + N2). The detector signal is
delta d = R^T delta phi = chi^T delta S.

7a  Infinite medium, full operator, delayed neutrons on: the N3 material
    of run_stage7_prep (transposed multiplicity 2.4) and the 2-group CORTEX
    fuel, 0.1, 1, 10, 30, 100 Hz, against ref_inf's adjoint. STD for both
    estimators (N3b promoted to a gate).
7b  Reciprocity with an external source on the Stage 4b slab (two groups,
    100 cm, CORTEX fuel homogenised with water, k from S_N). Forward: unit
    isotropic fast source in S = [-10, -5] cm, signal = thermal flux in
    D = [20, 25] cm. Adjoint: transposed library, unit isotropic thermal
    source in D (R = 1 in D, the TRIGA convention S_i = V_i Sigma_R,i with
    Sigma_R = 1), response = fast flux in S. With equal volumes
    (5 x 10 x 10 cm^3) the two tallies are equal. 1, 30 (the point-kinetic
    break, beta_eff / (2 pi Lambda_eff) = 29 Hz) and 100 Hz. Gate:
    two-sample test (Welch on Re and Im, Hotelling) between the forward and
    adjoint replica sets.
7d  S_N cross-check of 7b: the forward and adjoint 5 cm mesh fields (both
    groups) against ref_sn1d forward and transposed solutions. STD-C.
lf  7b's low-frequency sanity check: a 40 cm subcritical slab of the same
    mixture with all neutrons prompt, noise keff = 1 (fixed-source
    multiplication). The adjoint noise solution at 0.1 Hz (eta 1000, so
    that copies are frequent) against a plain
    (non-noise) OpenMC fixed-source run on the same transposed library,
    two-sample STD-C; each also against ref_sn1d (STD-C).
7c  Reciprocity with the sampled source of Stage 6 on the 7b slab, in the
    1 cm cell P = [30, 31] cm: (abs) thermal absorber, dSigma_a,2 =
    0.01 Sigma_a,2 (total only); (scat) scattering density change
    dSigma_s = 0.01 Sigma_s in both groups with dSigma_t = dSigma_s. The
    forward run gives delta d (thermal flux in D) and the static moments
    of phi_0 in P (Y00 and Y1m, NoisePhaseFilter static bin); the 7b
    adjoint runs give the moments of chi in P. The contraction is
    pt_lib.pt_terms, copied from ce_mg_triga (P0 and P1 parts, factor 3
    and sign of the P1 term from the direction reversal), applied to the
    complex chi through adjoint_noise_lib.pt_terms_complex. Replica r of
    the forward set is paired with replica r of the adjoint set; the paired
    differences are i.i.d., gate: paired test. S_N (ref_sn1d) shows that
    with P1 the 1 cm single-bin contraction is exact to 3e-4 here; without
    P1 the scattering case would be 25 % off, so the P1 term is tested.
7f  7b with the 32-bin histogram scattering library of 5e (mu_bar 0.3);
    k from an OpenMC eigenvalue run. Two-sample forward against adjoint.
7e  Not run (amendment 2026-11-20b): the code rejects negative frequencies,
    and under the bilinear convention -omega is never needed.

Log: delta d, FOM forward against adjoint for delta d, generation counts.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for sub in ("common", "reference", "stage3", "stage4", "stage7"):
    sys.path.insert(0, str(ROOT / sub))

import harness  # noqa: E402
from adjoint_noise_lib import (effective_operator, infinite_k,  # noqa: E402
                               make_adjoint_noise_library, pt_terms_complex,
                               second_moment_radius)
from build import box_source, make_library, slab_model, x_mesh  # noqa: E402
from harness import (fmt_info_correlated, paired_correlated,  # noqa: E402
                     run_replicas, std_correlated, std_criterion, summarize,
                     two_sample_correlated)
from noise_common import Material  # noqa: E402
from ref_inf import cortex_fuel, solve as inf_solve, static_k  # noqa: E402
from ref_sn1d import Slab  # noqa: E402
import run_stage3 as s3  # noqa: E402
import run_stage4 as s4  # noqa: E402
from run_stage7_prep import n3_material  # noqa: E402

WORK = HERE / "work"
LOG = []
SAVED = {}

# 7b/7c/7d/7f slab: cell 1 is the perturbed cell P
W = 100.0
X0 = -W / 2
REGIONS = [(80.0, "m"), (1.0, "m"), (19.0, "m")]
SRC, DET, PERT = (-10.0, -5.0), (20.0, 25.0), (30.0, 31.0)
NB = 20                      # 5 cm mesh bins
I_SRC, I_DET = 8, 14         # mesh bins of SRC and DET
V_BIN, V_P = 500.0, 100.0    # cm^3 (transverse 10 x 10)
FAST, THERMAL = 0, 1         # library group order (fastest first)
E_FAST, E_THERMAL = 1, 0     # EnergyFilter bins (ascending energy)
# lf slab
W_LF, NB_LF, DET_LF = 40.0, 8, (5.0, 10.0)
F_LF = 0.1
# At 0.1 Hz and eta = 1 a copy happens once per ~1e6 collisions and the Im
# estimate rests on a few hundred events; eta does not change the expected
# value, so a large eta makes copies frequent (copy-loop gain ~1e-3).
ETA_LF = 1000.0
ADJ_SEED_OFFSET = 1000003    # adjoint replica seeds differ from forward


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def eta_slab(f):
    return 3.0 if f >= 100 else 1.0     # as 4a and 6b


def eta_inf(f):
    return 5.0 if f >= 100 else 1.0


def floats(s):
    return [float(x) for x in s.split(",")] if s else []


def fom_mod(vals, times):
    """1 / (relative sigma of the modulus^2 x total time) of a scalar."""
    m = vals.mean()
    s = np.hypot(vals.real.std(ddof=1), vals.imag.std(ddof=1)) \
        / np.sqrt(vals.size)
    return 1.0 / ((s / abs(m))**2 * times.sum()), s / abs(m)


def fmt_c(z):
    return f"{z.real:+.5e}{z.imag:+.5e}j"


# ---------------------------------------------------------------- models
def slab(lib, f, keff, particles, batches, source=None, sampled=None,
         pt_tally=False, noise=True, regions=REGIONS, x0=X0, nb=NB,
         args=None, eta=None):
    """Slab model with a mesh tally 'tl' [x][E] (+ noise phase if sampled)
    and optionally a tally 'pt' on cell P: [E][Y00, Y1-1, Y10, Y11]."""
    model, edges = slab_model(regions, x0=x0, library=lib)
    cells = list(model.geometry.root_universe.cells.values())
    s = model.settings
    s.particles = particles
    s.batches = batches
    if noise:
        s.noise = {"frequency": f, "eta": eta or eta_slab(f), "keff": keff}
        if sampled is not None:
            s.noise.update({"nskip": args.nskip,
                            "static_inactive": args.static_inactive,
                            "perturbations": [{"cells": [cells[1]],
                                               **sampled}]})
    if source is not None:
        lo, hi, g = source
        s.source = [box_source(lo, hi, g, 2)]
    else:
        s.source = [box_source(edges[0], edges[-1], FAST, 2)]
    phase = [openmc.NoisePhaseFilter()] if sampled is not None else []
    cp = [openmc.ComplexPartFilter()] if noise else []
    energy = openmc.EnergyFilter(np.arange(3.0))
    t = openmc.Tally(name="tl")
    t.filters = [openmc.MeshFilter(x_mesh(edges[0], edges[-1], nb)),
                 energy] + phase + cp
    t.scores = ["flux"]
    tallies = [t]
    if pt_tally:
        t2 = openmc.Tally(name="pt")
        t2.filters = [openmc.CellFilter([cells[1]]), energy,
                      openmc.SphericalHarmonicsFilter(1)] + phase + cp
        t2.scores = ["flux"]
        tallies.append(t2)
    model.tallies = openmc.Tallies(tallies)
    return model


def plain_mean(name):
    """Replica extra: batch mean of a tally without ComplexPartFilter."""
    return lambda sp, d: sp.get_tally(name=name).mean.ravel().copy()


def gens(out):
    g = [e["generations_per_batch"].mean() for e in out["extras"]
         if "generations_per_batch" in e]
    return f"generations/batch {np.mean(g):.0f}" if g else ""


def read_gens(sp, d):
    return s3.read_noise(harness.last_statepoint(d))


# ---------------------------------------------------------------- 7a
def test_7a(args, gate):
    log("Test 7a: infinite medium, adjoint on the N2 library vs ref_inf "
        "adjoint (delayed neutrons on)")
    s3.WORK = WORK
    s3.LOG = LOG
    parts = dict(kv.split("=") for kv in args.particles_7a.split(","))
    results = []
    for name, mat in (("n3", n3_material()), ("cortex", cortex_fuel())):
        a = argparse.Namespace(reps=args.reps, particles=int(parts[name]),
                               batches=args.batches, base_seed=args.base_seed)
        fwd = WORK / f"lib_7a_{name}_forward.h5"
        adj = WORK / f"lib_7a_{name}_adjoint.h5"
        make_library({"m": mat}, fwd)
        rep = make_adjoint_noise_library(fwd, adj)
        k, _ = static_k(mat)
        op = effective_operator(adj, "m")
        log(f"  {name}: k_inf {k:.10f} (transposed library "
            f"{infinite_k(op):.10f}); max multiplicity "
            f"{rep['max_multiplicity']:.3f}; second-moment radius "
            f"{second_moment_radius(op):.3f}")
        G = mat.G
        R = np.eye(G)[G - 1]
        for f in floats(args.freqs_7a):
            ref = inf_solve(mat, 2 * np.pi * f, R, k=k, adjoint=True)[::-1]
            e = s3.run_case(f"7a_{name}_{f:g}Hz", mat, adj, f, eta_inf(f), k,
                            ref, a, results, src_group=G - 1)
            gate.append((e["label"], e["ok"]))
            for est in ("tl", "col"):
                SAVED[f"{e['label']}_{est}"] = e[est]["vals"]
            SAVED[f"{e['label']}_ref"] = ref
    z = np.concatenate([np.concatenate([r["tl"]["z"], r["col"]["z"]])
                        for r in results])
    log(f"  7a pooled z: n={z.size} mean={z.mean():+.2f} "
        f"std={z.std(ddof=1):.2f} |z|>=3: {int(np.sum(np.abs(z) >= 3))}")


# ---------------------------------------------------------------- slab
class SlabCase:
    """Libraries and S_N references for the 7b/7c/7d/7f slab."""

    def __init__(self, args):
        self.args = args
        self.mix = s4.mix_4b()
        self.lib_f = WORK / "lib_7b_forward.h5"
        self.lib_a = WORK / "lib_7b_adjoint.h5"
        make_library({"m": self.mix}, self.lib_f)
        make_adjoint_noise_library(self.lib_f, self.lib_a)
        self.sl = Slab([(W, self.mix, 400)], x0=X0, n_angles=64)
        self.k, _ = self.sl.eigenvalue()
        self.edges = np.linspace(X0, -X0, NB + 1)
        self.adj = {}
        log(f"Slab (4b mixture, {W:.0f} cm, S_N 400 cells S64): k = "
            f"{self.k:.8f}; source {SRC} fast, detector {DET} thermal, "
            f"perturbed cell {PERT}")

    def binned(self, vals):
        return s4.bin_to_mesh(self.sl, vals, self.edges)[:, ::-1].reshape(-1)

    def sn(self, f, adjoint):
        w = 2 * np.pi * f
        if adjoint:
            q = self.sl.isotropic_source(DET, [0.0, 1.0])
        else:
            q = self.sl.isotropic_source(SRC, [1.0, 0.0])
        return self.binned(self.sl.solve(w, q, adjoint=adjoint))

    def adjoint_runs(self, f):
        """7b adjoint replicas at f (also used by 7c)."""
        if f not in self.adj:
            a = self.args
            self.adj[f] = run_replicas(
                lambda: slab(self.lib_a, f, self.k, a.particles, a.batches,
                             source=(DET[0], DET[1], THERMAL),
                             pt_tally=True),
                WORK / f"7b_adj_{f:g}", n_rep=a.reps,
                base_seed=a.base_seed + ADJ_SEED_OFFSET,
                tally_names=("tl", "pt"), extra=read_gens)
        return self.adj[f]


def reciprocity(label, f, out_f, out_a, ref_dd, reps):
    """Two-sample test of forward delta d against the adjoint response."""
    vf = out_f["values"]["tl"].reshape(reps, NB, 2)
    va = out_a["values"]["tl"].reshape(reps, NB, 2)
    dd_f = vf[:, I_DET, E_THERMAL]
    dd_a = va[:, I_SRC, E_FAST]
    ok, info = two_sample_correlated(dd_f[:, None], dd_a[:, None])
    mf, ma = dd_f.mean(), dd_a.mean()
    fom_f, rel_f = fom_mod(dd_f, out_f["times"])
    fom_a, rel_a = fom_mod(dd_a, out_a["times"])
    k_, T2, p = info["hotelling"]
    line = (f"  {label} {f:g} Hz: {'PASS' if ok else 'FAIL'} two-sample "
            f"Welch t Re={info['z'][0]:+.2f} Im={info['z'][1]:+.2f}, "
            f"Hotelling p={p:.3f}; forward {fmt_c(mf)} adjoint {fmt_c(ma)}; "
            f"|fwd/adj|-1 = {abs(mf) / abs(ma) - 1:+.4f}, phase diff "
            f"{np.degrees(np.angle(mf / ma)):+.3f} deg; sigma rel fwd "
            f"{rel_f:.4f} adj {rel_a:.4f}; FOM fwd {fom_f:.3g} adj "
            f"{fom_a:.3g} (adj/fwd {fom_a / fom_f:.2f})")
    if ref_dd is not None:
        for nm, v in (("fwd", dd_f), ("adj", dd_a)):
            m, sre, sim, _ = summarize(v[:, None])
            _, inf = std_criterion(m, sre, sim, np.array([ref_dd]))
            line += (f"; {nm} vs S_N z Re={inf['z'][0]:+.2f} "
                     f"Im={inf['z'][1]:+.2f}")
        line += f"; S_N {fmt_c(ref_dd)}"
    log(line)
    log(f"    T fwd {out_f['times'].sum():.0f}s adj "
        f"{out_a['times'].sum():.0f}s; fwd {gens(out_f)}; adj {gens(out_a)}")
    SAVED[f"{label}_{f:g}_fwd"] = vf
    SAVED[f"{label}_{f:g}_adj"] = va
    SAVED[f"{label}_{f:g}_times"] = np.stack([out_f["times"],
                                              out_a["times"]])
    return ok


def test_7b(args, sc, gate):
    log("Test 7b: reciprocity, external source (and 7d: S_N fields)")
    for f in floats(args.freqs_7b):
        out_f = run_replicas(
            lambda: slab(sc.lib_f, f, sc.k, args.particles, args.batches,
                         source=(SRC[0], SRC[1], FAST)),
            WORK / f"7b_fwd_{f:g}", n_rep=args.reps,
            base_seed=args.base_seed, extra=read_gens)
        out_a = sc.adjoint_runs(f)
        ref_f, ref_a = sc.sn(f, False), sc.sn(f, True)
        ref_dd = ref_f.reshape(NB, 2)[I_DET, E_THERMAL]
        log(f"  S_N reciprocity {f:g} Hz: forward {fmt_c(ref_dd)}, adjoint "
            f"{fmt_c(ref_a.reshape(NB, 2)[I_SRC, E_FAST])}")
        gate.append((f"7b {f:g} Hz",
                     reciprocity("7b", f, out_f, out_a, ref_dd, args.reps)))
        for nm, out, ref in (("forward", out_f, ref_f),
                             ("adjoint", out_a, ref_a)):
            vals = out["values"]["tl"].reshape(args.reps, -1)
            mask = np.abs(ref) > 0.1 * np.abs(ref).max()
            ok, info = std_correlated(vals, ref, mask)
            log(f"  7d {nm} {f:g} Hz: {'PASS' if ok else 'FAIL'} (STD-C) "
                f"{fmt_info_correlated(info)}")
            gate.append((f"7d {nm} {f:g} Hz", ok))
            SAVED[f"7d_{nm}_{f:g}_ref"] = ref


def test_7c(args, sc, gate):
    log("Test 7c: reciprocity, sampled source in P (contraction "
        "pt_lib.pt_terms, P0 + P1)")
    mix = sc.mix
    sig_s = mix.scatter.sum(axis=1)
    sig_a = mix.total - sig_s
    G = mix.G
    cases = {
        "abs": ({"total": [0.0, 0.01 * sig_a[1] / mix.total[1]]},
                np.array([0.0, 0.01 * sig_a[1]]), np.zeros((G, G))),
        "scat": ({"total": list(0.01 * sig_s / mix.total), "scatter": 0.01},
                 0.01 * sig_s, 0.01 * mix.scatter),
    }
    phi0_P = sc.sl.region_average(sc.sl.phi0, PERT) * (PERT[1] - PERT[0])
    for case, (pert, dSt, dN0) in cases.items():
        if case not in args.cases_7c.split(","):
            continue
        dN1 = np.zeros((G, G))      # isotropic library
        for f in floats(args.freqs_7c):
            w = 2 * np.pi * f
            out_f = run_replicas(
                lambda: slab(sc.lib_f, f, sc.k, args.particles_7c,
                             args.batches, sampled=pert, pt_tally=True,
                             args=args),
                WORK / f"7c_{case}_{f:g}", n_rep=args.reps,
                base_seed=args.base_seed + 2 * ADJ_SEED_OFFSET,
                tally_names=("tl", "pt"), extra=read_gens)
            out_a = sc.adjoint_runs(f)
            R = args.reps
            vt = out_f["values"]["tl"].reshape(R, NB, 2, 2)
            dd = vt[:, I_DET, E_THERMAL, 1]
            vp = out_f["values"]["pt"].reshape(R, 2, 4, 2)
            F = vp[:, ::-1, :, 0].real                 # (R, G, 4) library
            A = V_BIN * out_a["values"]["pt"].reshape(R, 2, 4)[:, ::-1, :]
            P = np.array([pt_terms_complex(F[r][None], A[r][None],
                                           [1.0 / V_P], dSt, dN0, dN1)
                          for r in range(R)])           # (R, 2)
            c = P.sum(axis=1)
            ok, info = paired_correlated(dd[:, None], c[:, None])
            md, mc = dd.mean(), c.mean()
            # S_N, normalised by the static thermal flux integral over P
            q = sc.sl.perturbation_source(w, PERT, d_total=dSt,
                                          d_scatter=dN0)
            sn = sc.sl.solve(w, q)
            sn_dd = sc.binned(sn).reshape(NB, 2)[I_DET, E_THERMAL] * 5.0
            sn_ratio = sn_dd / phi0_P[THERMAL]
            norm = F[:, THERMAL, 0]
            zs = []
            for v in (dd / norm, c / norm):
                m, sre, sim, _ = summarize(v[:, None])
                _, inf = std_criterion(m, sre, sim, np.array([sn_ratio]))
                zs.append(inf["z"])
            fom_f, rel_f = fom_mod(dd, out_f["times"])
            t_ = info["t"]
            log(f"  7c {case} {f:g} Hz: {'PASS' if ok else 'FAIL'} paired "
                f"t Re={t_[0]:+.2f} Im={t_[1]:+.2f}, Hotelling "
                f"p={info['hotelling'][2]:.3f}; forward delta d {fmt_c(md)} "
                f"contraction {fmt_c(mc)} (P0 {fmt_c(P[:, 0].mean())}, P1 "
                f"{fmt_c(P[:, 1].mean())}); |fwd/contr|-1 = "
                f"{abs(md) / abs(mc) - 1:+.4f}, phase diff "
                f"{np.degrees(np.angle(md / mc)):+.3f} deg; sigma rel fwd "
                f"{rel_f:.4f}")
            log(f"    vs S_N (delta d / static thermal flux in P, "
                f"{fmt_c(sn_ratio)}): forward z Re={zs[0][0]:+.2f} "
                f"Im={zs[0][1]:+.2f}, contraction z Re={zs[1][0]:+.2f} "
                f"Im={zs[1][1]:+.2f}; T fwd {out_f['times'].sum():.0f}s; "
                f"{gens(out_f)}")
            gate.append((f"7c {case} {f:g} Hz", ok))
            SAVED[f"7c_{case}_{f:g}_dd"] = dd
            SAVED[f"7c_{case}_{f:g}_F"] = F
            SAVED[f"7c_{case}_{f:g}_A"] = A
            SAVED[f"7c_{case}_{f:g}_P"] = P
            SAVED[f"7c_{case}_{f:g}_sn_ratio"] = sn_ratio


def test_lf(args, gate):
    mix = s4.mix_4b()
    mat = Material(mix.total, mix.scatter, mix.inv_v,
                   nu_fission_p=mix.nu_fission, chi_p=mix.chi_p)
    sl = Slab([(W_LF, mat, int(W_LF * 4))], x0=-W_LF / 2, n_angles=64)
    k, _ = sl.eigenvalue()
    log(f"Test 7b low frequency: {W_LF:.0f} cm slab, prompt only, S_N k = "
        f"{k:.6f}; noise keff = 1 (fixed-source multiplication), "
        f"{F_LF} Hz, eta {ETA_LF:g}; thermal adjoint source in {DET_LF}")
    lib_f = WORK / "lib_lf_forward.h5"
    lib_a = WORK / "lib_lf_adjoint.h5"
    make_library({"m": mat}, lib_f)
    make_adjoint_noise_library(lib_f, lib_a)
    edges = np.linspace(-W_LF / 2, W_LF / 2, NB_LF + 1)
    q = sl.isotropic_source(DET_LF, [0.0, 1.0])
    ref_n = s4.bin_to_mesh(sl, sl.solve(2 * np.pi * F_LF, q, k=1.0,
                                        adjoint=True), edges)
    ref_s = s4.bin_to_mesh(sl, sl.solve(0.0, q, k=1.0, adjoint=True), edges)
    ref_n, ref_s = (r[:, ::-1].reshape(-1) for r in (ref_n, ref_s))
    log(f"  S_N: max |Im/Re| of the 0.1 Hz adjoint "
        f"{np.max(np.abs(ref_n.imag / ref_n.real)):.2e}; max |noise/static "
        f"- 1| (Re) {np.max(np.abs(ref_n.real / ref_s.real - 1)):.2e}")
    regions = [(W_LF, "m")]
    src = (DET_LF[0], DET_LF[1], THERMAL)
    out_n = run_replicas(
        lambda: slab(lib_a, F_LF, 1.0, args.particles, args.batches,
                     source=src, regions=regions, x0=-W_LF / 2, nb=NB_LF,
                     eta=ETA_LF),
        WORK / "lf_noise", n_rep=args.reps,
        base_seed=args.base_seed + ADJ_SEED_OFFSET, extra=read_gens)
    out_s = run_replicas(
        lambda: slab(lib_a, 0.0, 1.0, args.particles, args.batches,
                     source=src, noise=False, regions=regions,
                     x0=-W_LF / 2, nb=NB_LF),
        WORK / "lf_static", n_rep=args.reps,
        base_seed=args.base_seed + 3 * ADJ_SEED_OFFSET, tally_names=(),
        extra=plain_mean("tl"))
    vn = out_n["values"]["tl"].reshape(args.reps, -1)
    vs = np.array(out_s["extras"]).astype(complex)
    mask = np.abs(ref_s) > 0.1 * np.abs(ref_s).max()
    # Re only: the noise solution has a small physical Im (S_N above),
    # which the static run lacks; Im is checked against S_N below.
    ok, info = two_sample_correlated(vn.real.astype(complex), vs, mask)
    k_, T2, p = info["hotelling"]
    rel = vn.mean(0).real[mask] / vs.mean(0).real[mask] - 1
    log(f"  lf noise vs static (Re): {'PASS' if ok else 'FAIL'} two-sample "
        f"STD-C n={info['z'].size} bad={info['n_bad']} max|t|="
        f"{np.abs(info['z']).max():.2f}; integral t Re="
        f"{info['z_int'][0]:+.2f} Im={info['z_int'][1]:+.2f}; Hotelling "
        f"k={k_} p={p:.3f}; Re noise/static - 1: max |.| "
        f"{np.abs(rel).max():.4f}")
    gate.append(("7b lf noise vs static", ok))
    for nm, vals, ref in (("noise", vn, ref_n), ("static", vs, ref_s)):
        ok, info = std_correlated(vals, ref, mask)
        log(f"  lf {nm} vs S_N: {'PASS' if ok else 'FAIL'} (STD-C) "
            f"{fmt_info_correlated(info)}")
        gate.append((f"7b lf {nm} vs S_N", ok))
    log(f"    T noise {out_n['times'].sum():.0f}s static "
        f"{out_s['times'].sum():.0f}s; noise {gens(out_n)}")
    SAVED.update({"lf_noise": vn, "lf_static": vs, "lf_ref_noise": ref_n,
                  "lf_ref_static": ref_s})


def test_7f(args, gate):
    mix = s4.mix_4b()
    lib_f = WORK / "lib_7f_forward.h5"
    lib_a = WORK / "lib_7f_adjoint.h5"
    make_library({"m": mix}, lib_f, histogram=(0.3, 32))
    make_adjoint_noise_library(lib_f, lib_a)
    s4.WORK = WORK
    k, sk = s4.openmc_keff(lib_f, REGIONS, X0, particles=args.keff_particles,
                           batches=args.keff_batches)
    log(f"Test 7f: 7b with the 32-bin histogram library (mu_bar 0.3); "
        f"OpenMC k = {k:.6f} +- {sk:.6f}")
    for f in floats(args.freqs_7f):
        out_f = run_replicas(
            lambda: slab(lib_f, f, k, args.particles, args.batches,
                         source=(SRC[0], SRC[1], FAST)),
            WORK / f"7f_fwd_{f:g}", n_rep=args.reps,
            base_seed=args.base_seed, extra=read_gens)
        out_a = run_replicas(
            lambda: slab(lib_a, f, k, args.particles, args.batches,
                         source=(DET[0], DET[1], THERMAL)),
            WORK / f"7f_adj_{f:g}", n_rep=args.reps,
            base_seed=args.base_seed + ADJ_SEED_OFFSET, extra=read_gens)
        gate.append((f"7f {f:g} Hz",
                     reciprocity("7f", f, out_f, out_a, None, args.reps)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--only", default="7a,7b,lf,7c,7f")
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=20000,
                   help="particles per batch, slab cases (7b, lf, 7f)")
    p.add_argument("--particles-7a", default="n3=80000,cortex=1000000",
                   help="particles per batch per 7a material")
    p.add_argument("--particles-7c", type=int, default=60000,
                   help="static particles per generation (sampled source)")
    p.add_argument("--batches", type=int, default=20)
    p.add_argument("--nskip", type=int, default=20)
    p.add_argument("--static-inactive", type=int, default=300)
    p.add_argument("--freqs-7a", default="0.1,1,10,30,100")
    p.add_argument("--freqs-7b", default="1,30,100")
    p.add_argument("--freqs-7c", default="1,30,100")
    p.add_argument("--cases-7c", default="abs,scat")
    p.add_argument("--freqs-7f", default="1,30")
    p.add_argument("--keff-particles", type=int, default=100000)
    p.add_argument("--keff-batches", type=int, default=300)
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

    if "7a" in only:
        test_7a(args, gate)
    if "7b" in only or "7c" in only:
        sc = SlabCase(args)
        if "7b" in only:
            test_7b(args, sc, gate)
        if "7c" in only:
            test_7c(args, sc, gate)
    if "lf" in only:
        test_lf(args, gate)
    if "7f" in only:
        test_7f(args, gate)

    n_fail = sum(not ok for _, ok in gate)
    log(f"\nStage 7 ({args.only}): {len(gate) - n_fail}/{len(gate)} gate "
        "checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAILED: {name}")
    tag = args.only.replace(",", "_") + sfx
    with open(HERE / f"stage7_log_{tag}.md", "w") as fh:
        fh.write("# Stage 7 log\n\n```\n" + "\n".join(LOG) + "\n```\n")
    np.savez(HERE / f"stage7_vals_{tag}.npz", **SAVED)
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
