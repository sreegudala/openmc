"""Stage 9b: TRIGA dry run (design doc, Stage 9; amendment 2026-11-26c).

ce_mg_triga's zoned TRIGA model with the N5 CASMO-70 library (vector form,
noise.keff 1.00137) and the N5 central absorber: absorption of the central
thimble water (cell 3789, r = 1.689 cm, z = -80 to 80 cm) raised by EPS.
Detectors in the pool from ce_mg_triga zones.json: point detectors det0-5
and arcs arc0-5 (r = 55.8 to 80.8 cm), 1/v response
R_g = v0 * (1/v)_g, v0 = 2.2e5 cm/s, from each detector's own water.

Modes:
  forward    R forward replicas at one frequency with cancellation mesh A
             (aligned: one half fuel element by one fuel segment per cell)
             or B (A with every interval halved)
             -> stage9b_fwd_{A|B}_{f}Hz.npz: detector flux per group (static,
             noise), total fission (static, noise), static moments F in the
             central cell (axial voxels x groups x P0 + P1), per-generation
             sum|w|, static k per generation; with --wwmap also the static
             flux on the window mesh
  windows    CADIS lower bounds C / phi0(r, G) for the adjoint runs from a
             importance map -> stage9b_ww.npz (1/phi and C per detector)
  adjoint    R adjoint replicas (N1 + N2 transposed library) at one
             frequency for one detector: source uniform in the detector,
             groups by R_g, total strength S = V sum_g R_g; windows cadis
             (stage9b_ww.npz) or flat (0.5-2 in one cell), with
             cancellation mesh A (--adj-cancel)
             -> stage9b_adj_{det}_{ww}_{f}Hz.npz: adjoint moments A in the
             central cell, per-generation sum|w|, and the transfer map
             G(omega) (adjoint flux on the window mesh, mean and sigma)
  linearity  one replica at EPS 0.01, 0.1 and 1 with one seed (local check)
  analyze    the gates and logged checks -> stage9b_log.md

Contraction (as 7c): delta d = <chi, dS> = pt_terms_complex(F, S A, 1/V_v,
dSt = EPS Sigma_a, 0, 0), P0 + P1, on axial voxels of the central cell, with
forward replica r paired with adjoint replica r (harness.paired_correlated).
"""

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import openmc
import openmc.mgxs

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "common"))
sys.path.insert(0, str(ROOT / "prep9b"))
sys.path.insert(0, str(ROOT / "reference"))
sys.path.insert(0, str(ROOT / "stage3"))

from harness import (last_statepoint, paired_correlated,  # noqa: E402
                     run_replicas)
from adjoint_noise_lib import pt_terms_complex  # noqa: E402
import ref_g0  # noqa: E402
import run_n5 as n5  # noqa: E402
from run_stage3 import gen_ratio  # noqa: E402

EPS = n5.EPS_ABS
KEFF = n5.K_FORWARD["vector"][0]
F7 = n5.F7
GATE_DETS = ("arc0", "arc3", "arc5")
V0 = 2.2e5                          # cm/s, 1/v response normalisation
CENTRAL_R = np.sqrt(2.85305881)     # surface 11318
CENTRAL_Z = (-80.0, 80.0)
NZ = 160                            # 1 cm axial voxels in the central cell
WW_STEP = 5.0                       # cm, window mesh
WW_RATIO = 5.0
ADJ_SEED_OFFSET = 100000
DRHO_N5 = -32.64e-5                 # prep9b/n5_log.md, IFP fit, EPS = 0.1
LOG = []


def log(msg=""):
    print(msg, flush=True)
    LOG.append(msg)


def eta(f):
    return 1.0 if f <= 30 else 3.0


# ------------------------------------------------------------------ data
def zones():
    z = json.loads((n5.triga_dir() / "model_zoned" / "zones.json")
                   .read_text())
    return z


def detectors():
    """Point detectors and arcs, in zones.json order."""
    z = zones()
    out = []
    for d in z["detectors"] + z["arcs"]:
        out.append({"name": d["name"], "cell_id": d["cell_id"],
                    "material": f"water_shell{d['shell']}_{d['material_id']}",
                    "volume": d.get("volume_cm3"), "r_cm": d["r_cm"],
                    "theta_deg": d.get("theta_deg")})
    return out


def vector_library():
    return (n5.lib_dir() / n5.LIBS["vector"]).resolve()


def library_data():
    """Group edges (ascending, eV), 1/v per detector material and Sigma_a,
    Sigma_t of the central water, library order (highest energy first)."""
    with h5py.File(vector_library(), "r") as f:
        edges = np.array(f.attrs["group structure"])

        def ds(name, key):
            g = f[name]
            return g[next(k for k in g if k != "kTs")][key][()]
        invv = {d["material"]: ds(d["material"], "inverse-velocity")
                for d in detectors()}
        sa = ds(n5.CENTRAL_MATERIAL, "absorption")
        st = ds(n5.CENTRAL_MATERIAL, "total")
    return edges, invv, sa, st


def response(det, invv):
    return V0 * invv[det["material"]]


# ---------------------------------------------------------------- meshes
def cancellation_mesh(kind):
    """A: x edges between the hexagonal lattice's element columns, y edges
    at the element centres (each cell holds one half element), z edges at
    the fuel-segment boundaries of both element types (zones.json). B: A
    with every interval halved."""
    segs = zones()["fuel_segments"]
    xy = np.array([s["xy_cm"] for s in segs])
    ux = np.unique(xy[:, 0].round(3))
    uy = np.unique(xy[:, 1].round(3))
    dx, dy = np.diff(ux).mean(), np.diff(uy).mean()
    x = np.concatenate([ux - dx / 2, [ux[-1] + dx / 2]])
    y = np.concatenate([[uy[0] - dy], uy, [uy[-1] + dy]])
    zb = np.sort(np.concatenate([[s["box_cm"][0][2], s["box_cm"][1][2]]
                                 for s in segs]))
    clusters = [[zb[0]]]
    for v in zb[1:]:
        if v - clusters[-1][-1] < 0.2:
            clusters[-1].append(v)
        else:
            clusters.append([v])
    z = np.array([np.mean(c) for c in clusters])
    if kind == "B":
        half = lambda e: np.sort(np.concatenate([e, (e[1:] + e[:-1]) / 2]))
        x, y, z = half(x), half(y), half(z)
    m = openmc.RectilinearMesh()
    m.x_grid, m.y_grid, m.z_grid = x, y, z
    return m


def central_zmesh():
    m = openmc.RegularMesh()
    m.lower_left = (-CENTRAL_R, -CENTRAL_R, CENTRAL_Z[0])
    m.upper_right = (CENTRAL_R, CENTRAL_R, CENTRAL_Z[1])
    m.dimension = (1, 1, NZ)
    return m


def voxel_volume():
    return np.pi * CENTRAL_R**2 * (CENTRAL_Z[1] - CENTRAL_Z[0]) / NZ


def window_mesh(model):
    bb = model.geometry.bounding_box
    lo, hi = np.array(bb.lower_left), np.array(bb.upper_right)
    n = np.maximum(1, np.ceil((hi - lo) / WW_STEP).astype(int))
    m = openmc.RegularMesh()
    m.lower_left, m.upper_right, m.dimension = lo, hi, n
    return m


def ww_energy():
    return np.array(openmc.mgxs.GROUP_STRUCTURES["CASMO-8"])


# ---------------------------------------------------------------- models
def base_model(library, particles, batches):
    model = n5.mg_model(library, particles, 0, batches, ww=None)
    s = model.settings
    s.run_mode = "fixed source"
    s.batches, s.inactive = batches, 0
    return model


def pt_tally(cells, edges, phase):
    t = openmc.Tally(name="pt")
    t.filters = ([openmc.CellFilter(cells),
                  openmc.MeshFilter(central_zmesh()),
                  openmc.EnergyFilter(edges),
                  openmc.SphericalHarmonicsFilter(1)]
                 + ([openmc.NoisePhaseFilter()] if phase else [])
                 + [openmc.ComplexPartFilter()])
    t.scores = ["flux"]
    return t


def forward_model(f, kind, particles, batches, nskip, static_inactive,
                  wwmap=False, eps=EPS):
    lib = vector_library()
    model = base_model(lib, particles, batches)
    cells = model.geometry.get_all_cells()
    cell = cells[n5.CENTRAL_CELL]
    if getattr(cell.fill, "_macroscopic", None) != n5.CENTRAL_MATERIAL:
        raise SystemExit(f"cell {n5.CENTRAL_CELL} is not "
                         f"{n5.CENTRAL_MATERIAL}")
    amp = n5.central_amplitudes(lib) * (eps / n5.EPS_ABS)
    s = model.settings
    s.noise = {"frequency": f, "eta": eta(f), "keff": KEFF, "nskip": nskip,
               "static_inactive": static_inactive,
               "cancellation_mesh": cancellation_mesh(kind),
               "perturbations": [{"cells": [cell],
                                  "total": [complex(a) for a in amp]}]}
    edges, _, _, _ = library_data()
    dets = [cells[d["cell_id"]] for d in detectors()]
    t_det = openmc.Tally(name="det")
    t_det.filters = [openmc.CellFilter(dets), openmc.EnergyFilter(edges),
                     openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
    t_det.scores = ["flux"]
    t_fis = openmc.Tally(name="fission")
    t_fis.filters = [openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
    t_fis.scores = ["fission"]
    tallies = [t_det, t_fis, pt_tally([cell], edges, True)]
    if wwmap:
        t = openmc.Tally(name="wwmap")
        t.filters = [openmc.MeshFilter(window_mesh(model)),
                     openmc.EnergyFilter(ww_energy()),
                     openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
        t.scores = ["flux"]
        tallies.append(t)
    model.tallies = openmc.Tallies(tallies)
    return model


def arc_source(det, edges, R):
    """Uniform in the annular sector of an arc (thickness 2.54 cm, height
    10 cm, zones.json theta range), groups with probability R_g."""
    r_c = det["r_cm"]
    th = np.radians(det["theta_deg"])
    space = openmc.stats.CylindricalIndependent(
        r=openmc.stats.PowerLaw(r_c - 1.27, r_c + 1.27, 1),
        phi=openmc.stats.Uniform(th[0], th[1]),
        z=openmc.stats.Uniform(-5.0, 5.0))
    mid = np.sqrt(np.maximum(edges[:-1], 1e-5) * edges[1:])[::-1]
    energy = openmc.stats.Discrete(mid, R / R.sum())
    return openmc.IndependentSource(space=space, energy=energy,
                                    angle=openmc.stats.Isotropic())


def adjoint_model(f, det, particles, batches, ww="cadis", cancel="A"):
    work = HERE / "work_stage9b_lib"
    work.mkdir(exist_ok=True)
    lib = n5.adjoint_library("vector", work)
    model = base_model(lib, particles, batches)
    s = model.settings
    s.noise = {"frequency": f, "eta": eta(f), "keff": KEFF}
    if cancel != "none":
        # the adjoint's transposed multiplicity makes its offspring heavy-
        # tailed; without cancellation the phase is lost (pilot log)
        s.noise["cancellation_mesh"] = cancellation_mesh(cancel)
    edges, invv, _, _ = library_data()
    s.source = [arc_source(det, edges, response(det, invv))]
    if ww == "none":
        wws = None
    elif ww == "flat":
        bb = model.geometry.bounding_box
        mesh = openmc.RegularMesh()
        mesh.lower_left, mesh.upper_right = bb.lower_left, bb.upper_right
        mesh.dimension = (1, 1, 1)
        wws = openmc.WeightWindows(mesh, lower_ww_bounds=[0.5],
                                   upper_ww_bounds=[2.0],
                                   particle_type="neutron")
    else:
        d = np.load(HERE / "stage9b_ww.npz")
        mesh = window_mesh(model)
        if tuple(mesh.dimension) != tuple(d["dimension"]):
            raise SystemExit("window mesh differs from stage9b_ww.npz")
        C = float(d["C"][list(d["names"]).index(det["name"])])
        lower = C * d["inv_phi"].astype(float)
        wws = openmc.WeightWindows(
            mesh, lower_ww_bounds=lower,
            upper_ww_bounds=lower * WW_RATIO,
            energy_bounds=ww_energy(), particle_type="neutron")
    if wws is not None:
        s.weight_windows = [wws]
        s.weight_windows_on = True
        s.max_history_splits = 10000
    cell = model.geometry.get_all_cells()[n5.CENTRAL_CELL]
    gmap = openmc.Tally(name="gmap")         # the transfer map G(omega)
    gmap.filters = [openmc.MeshFilter(window_mesh(model)),
                    openmc.EnergyFilter(ww_energy()),
                    openmc.ComplexPartFilter()]
    gmap.scores = ["flux"]
    model.tallies = openmc.Tallies([pt_tally([cell], edges, False), gmap])
    return model


# --------------------------------------------------------------- running
def read_noise(sp, d):
    """Replica extra: the statepoint's per-generation noise records."""
    with h5py.File(last_statepoint(d), "r") as f:
        g = f["noise"]
        keys = ("generations_per_batch", "gen_sum_abs_wgt", "gen_n_batches",
                "static_k")
        return {k: g[k][()] for k in keys if k in g}


def pad_stack(arrs):
    n = max(len(a) for a in arrs)
    out = np.zeros((len(arrs), n))
    for i, a in enumerate(arrs):
        out[i, :len(a)] = a
    return out


def run_forward(args):
    f, kind = args.freq, args.mesh
    work = HERE / f"work_stage9b_fwd_{kind}_{f:g}Hz"
    work.mkdir(exist_ok=True)
    names = ("det", "fission", "pt") + (("wwmap",) if args.wwmap else ())
    out = run_replicas(
        lambda: forward_model(f, kind, args.particles, args.batches,
                              args.nskip, args.static_inactive, args.wwmap),
        work / "reps", n_rep=args.reps, tally_names=names,
        base_seed=args.base_seed, extra=read_noise)
    R = args.reps
    v = out["values"]
    G = library_data()[0].size - 1
    nd = len(detectors())
    ex = out["extras"]
    res = {"det": v["det"].reshape(R, nd, G, 2),
           "fission": v["fission"].reshape(R, 2),
           "pt": v["pt"].reshape(R, NZ, G, 4, 2)[..., 0].real,
           "times": out["times"], "frequency": f,
           "gen_sum_abs_wgt": pad_stack([e["gen_sum_abs_wgt"] for e in ex]),
           "gen_n_batches": pad_stack([e["gen_n_batches"] for e in ex]),
           "gens_per_batch": pad_stack([e["generations_per_batch"]
                                        for e in ex]),
           "static_k": pad_stack([e["static_k"] for e in ex])}
    if args.wwmap:
        res["wwmap"] = v["wwmap"].reshape(R, -1, ww_energy().size - 1,
                                          2)[..., 0].real
    path = HERE / f"stage9b_fwd_{kind}_{f:g}Hz.npz"
    np.savez(path, **res)
    r = res["fission"][:, 1] / res["fission"][:, 0].real
    log(f"forward {kind} {f:g} Hz: dF/F0 {r.mean():.4e}, T "
        f"{out['times'].sum():.0f}s, generations/batch "
        f"{res['gens_per_batch'][res['gens_per_batch'] > 0].mean():.0f} "
        f"-> {path.name}")


def importance_map(args):
    """Forward flux phi(r, G) on the window mesh (cells x-fastest, CASMO-8
    ascending), the importance of adjoint particles.

    rr (default): ce_mg_triga's random-ray forward flux
    (results/step3/forward_flux.npz, step3a_adjoint.py forward: fission-free,
    from the core's CE source strengths, CASMO-8), read in place. It covers
    the pool, where a forward Monte Carlo map is barely scored.
    forward: the static flux of a stage9b forward run (--wwmap)."""
    if args.importance == "rr":
        src = n5.triga_dir() / "results" / "step3" / "forward_flux.npz"
        d = np.load(src)
        e8 = ww_energy()
        if not np.allclose(d["energy_edges"], e8):
            raise SystemExit(f"{src}: not CASMO-8")
        flux = d["flux"]                        # (nx, ny, nz, G)
        return src, (np.array(d["lower_left"]), np.array(d["upper_right"]),
                     flux.shape[:3]), \
            flux.transpose(2, 1, 0, 3).reshape(-1, flux.shape[3])
    src = HERE / args.wwmap_from
    return src, None, np.load(src)["wwmap"].mean(axis=0)


def run_windows(args):
    """CADIS lower bounds for the adjoint particles of detector det:
    lower = C / phi(r, G), upper = 5 lower, with phi from importance_map.
    Where the point-kinetic term dominates, the forward solution for the
    contraction's source has the shape of the static flux; a fission-free
    map differs from it inside the core, which costs efficiency, not bias.
    C puts an R-weighted source particle at its window centre (3 lower).
    Cells with phi = 0 get the group's largest lower bound, so every cell
    has a valid window: a cell without one would get neither windows nor
    the noise roulette. The infinite-medium tests behind this choice
    (energy-resolved importance and cancellation are both needed; CASMO-8
    is enough) are in stage9b_log_pilot.md."""
    src, geo, phi = importance_map(args)
    phi = np.where(phi > 0, phi, np.nan)
    fill = np.nanmin(phi, axis=0)
    n_fill = int(np.isnan(phi).sum())
    phi = np.where(np.isnan(phi), fill[None, :], phi)
    model = base_model(vector_library(), 10, 1)
    mesh = window_mesh(model)
    if geo is not None and not (
            np.allclose(geo[0], mesh.lower_left)
            and np.allclose(geo[1], mesh.upper_right)
            and tuple(geo[2]) == tuple(mesh.dimension)):
        raise SystemExit(f"{src} is not on the window mesh {mesh}")
    edges, invv, _, _ = library_data()
    e8 = ww_energy()
    det_list = {x["name"]: x for x in detectors()}
    names = args.detectors.split(",")
    inv = (1.0 / phi).reshape(*np.array(mesh.dimension)[::-1], -1) \
        .transpose(2, 1, 0, 3).astype(np.float32)
    Cs = {}
    for name in names:
        det = det_list[name]
        R = response(det, invv)[::-1]           # ascending energy
        mid = np.sqrt(np.maximum(edges[:-1], 1e-5) * edges[1:])
        G8 = np.clip(np.searchsorted(e8, mid) - 1, 0, e8.size - 2)
        R8 = np.bincount(G8, weights=R, minlength=e8.size - 1)
        # window cells inside the arc: sample points of the sector
        rng = np.random.default_rng(1)
        n = 20000
        rr = np.sqrt(rng.uniform((det["r_cm"] - 1.27)**2,
                                 (det["r_cm"] + 1.27)**2, n))
        tt = np.radians(rng.uniform(*det["theta_deg"], n))
        pts = np.stack([rr * np.cos(tt), rr * np.sin(tt),
                        rng.uniform(-5, 5, n)], axis=1)
        ll = np.array(mesh.lower_left)
        width = (np.array(mesh.upper_right) - ll) / np.array(mesh.dimension)
        ijk = np.floor((pts - ll) / width).astype(int)
        dim = np.array(mesh.dimension)
        flat = ijk[:, 0] + dim[0] * (ijk[:, 1] + dim[1] * ijk[:, 2])
        phi_src = (R8[None, :] * phi[flat]).sum(1).mean() / R8.sum()
        C = phi_src / 3.0
        Cs[name] = C
        log(f"windows {name}: from {src.name}, mesh {tuple(dim)} x "
            f"{e8.size - 1} groups, C {C:.3e}, {n_fill} unscored "
            f"(cell, group) entries filled; lower bounds "
            f"{C * inv.min():.2e} to {C * inv.max():.2e}")
    path = HERE / "stage9b_ww.npz"
    np.savez_compressed(path, inv_phi=inv, dimension=dim, source=str(src),
                        names=np.array(list(Cs)),
                        C=np.array(list(Cs.values())),
                        n_filled=n_fill)
    log(f"-> {path.name} (lower bound = C[detector] * inv_phi)")


def run_adjoint(args):
    f = args.freq
    det = {x["name"]: x for x in detectors()}[args.detector]
    tag = (f"{det['name']}_{args.ww}"
           + ("" if args.adj_cancel == "A" else f"_c{args.adj_cancel}")
           + f"_{f:g}Hz")
    work = HERE / f"work_stage9b_adj_{tag}"
    work.mkdir(exist_ok=True)
    out = run_replicas(
        lambda: adjoint_model(f, det, args.adj_particles, args.adj_batches,
                              args.ww, args.adj_cancel),
        work / "reps", n_rep=args.reps, tally_names=("pt", "gmap"),
        base_seed=args.base_seed + ADJ_SEED_OFFSET, extra=read_noise)
    R = args.reps
    G = library_data()[0].size - 1
    ex = out["extras"]
    gm = out["values"]["gmap"].reshape(R, -1, ww_energy().size - 1)
    path = HERE / f"stage9b_adj_{tag}.npz"
    np.savez(path, pt=out["values"]["pt"].reshape(R, NZ, G, 4),
             gmap_mean=gm.mean(0), gmap_sre=gm.real.std(0, ddof=1) / R**0.5,
             gmap_sim=gm.imag.std(0, ddof=1) / R**0.5,
             times=out["times"], frequency=f,
             gen_sum_abs_wgt=pad_stack([e["gen_sum_abs_wgt"] for e in ex]),
             gen_n_batches=pad_stack([e["gen_n_batches"] for e in ex]),
             gens_per_batch=pad_stack([e["generations_per_batch"]
                                       for e in ex]))
    log(f"adjoint {tag}: T {out['times'].sum():.0f}s -> {path.name}")


def run_linearity(args):
    """The first-order game is linear in the perturbation: the same seed at
    three strengths must give responses in the ratio of the strengths."""
    f = args.freq
    vals = {}
    for eps in (0.01, 0.1, 1.0):
        work = HERE / f"work_stage9b_lin_{eps:g}"
        work.mkdir(exist_ok=True)
        out = run_replicas(
            lambda: forward_model(f, "A", args.particles, args.batches,
                                  args.nskip, args.static_inactive, eps=eps),
            work / "reps", n_rep=1, tally_names=("det", "fission"),
            base_seed=args.base_seed)
        vals[eps] = (out["values"]["det"].reshape(-1, 2)[:, 1],
                     out["values"]["fission"].reshape(2)[1])
    for eps in (0.01, 1.0):
        rd = vals[eps][0] / vals[0.1][0] * (0.1 / eps)
        rd = rd[np.isfinite(rd)]
        rf = vals[eps][1] / vals[0.1][1] * (0.1 / eps)
        log(f"linearity eps {eps:g} against 0.1 ({f:g} Hz, one seed): "
            f"fission ratio - 1 = {abs(rf - 1):.1e}; detector bins max "
            f"|ratio - 1| = {np.abs(rd - 1).max():.1e} ({rd.size} bins)")
    (HERE / "stage9b_log_linearity.md").write_text(
        "# Stage 9b linearity\n\n```\n" + "\n".join(LOG) + "\n```\n")


# --------------------------------------------------------------- analyze
def contraction(F, A, S, dSt, nmerge=1):
    """P0 + P1 contraction per replica pair; voxels merged nmerge at a time
    (for the voxel-size check)."""
    R = F.shape[0]
    if nmerge > 1:
        F = F.reshape(R, NZ // nmerge, nmerge, *F.shape[2:]).sum(2)
        A = A.reshape(R, NZ // nmerge, nmerge, *A.shape[2:]).sum(2)
    w = np.full(F.shape[1], 1.0 / (voxel_volume() * nmerge))
    P = np.array([pt_terms_complex(F[r], S * A[r], w, dSt, np.zeros((dSt.size,
                  dSt.size)), np.zeros((dSt.size, dSt.size)))
                  for r in range(R)])
    return P                                    # (R, 2): P0, P1


def rel_sigma(x):
    """Relative standard error of a complex replica mean."""
    return (np.hypot(x.real.std(ddof=1), x.imag.std(ddof=1))
            / np.sqrt(x.size) / abs(x.mean()))


def sumw_check(d, label):
    S = d["gen_sum_abs_wgt"].sum(0)
    nb = d["gen_n_batches"].sum(0)
    ok_len = S.size > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        Sm = np.where(nb > 0, S / np.maximum(nb, 1), 0)
    ratio, npts = gen_ratio(S, nb) if ok_len else (np.nan, 0)
    peak = Sm.max() / Sm[0] if ok_len and Sm[0] > 0 else np.nan
    ok = bool(np.isfinite(ratio) and ratio < 1)
    log(f"    sum|w| {label}: generation ratio {ratio:.4f} ({npts} "
        f"generations fitted), peak/first {peak:.2f}, generations/batch "
        f"{d['gens_per_batch'][d['gens_per_batch'] > 0].mean():.0f} -> "
        f"{'ok' if ok else 'FAIL'}")
    return ok


def analyze(args):
    edges, invv, sa, st = library_data()
    dets = detectors()
    names = [d["name"] for d in dets]
    dSt = EPS * sa
    gate = []
    log("# Stage 9b: TRIGA dry run")
    log(f"EPS {EPS}, noise.keff {KEFF}, detectors {', '.join(names)}; "
        f"reciprocity detectors {', '.join(GATE_DETS)}")
    for f in F7:
        pa = HERE / f"stage9b_fwd_A_{f:g}Hz.npz"
        if not pa.exists():
            continue
        log(f"\n## {f:g} Hz (eta {eta(f):g})")
        fa = np.load(pa)
        R = fa["det"].shape[0]
        Rg = np.array([response(d, invv) for d in dets])[:, ::-1]  # ascend
        d0 = (fa["det"][..., 0].real * Rg).sum(-1)        # (R, nd)
        dd = (fa["det"][..., 1] * Rg).sum(-1)
        F0 = fa["fission"][:, 0].real
        dF = fa["fission"][:, 1] / F0
        # ---- reciprocity
        F = fa["pt"][:, :, ::-1, :]                       # library order
        fwd, con, ok_all = [], [], True
        for name in GATE_DETS:
            pj = HERE / f"stage9b_adj_{name}_{args.ww}_{f:g}Hz.npz"
            if not pj.exists():
                log(f"  {name}: no {pj.name}")
                ok_all = False
                continue
            aj = np.load(pj)
            det = dets[names.index(name)]
            S = det["volume"] * response(det, invv).sum()
            A = aj["pt"][:, :, ::-1, :]
            n = min(R, A.shape[0])
            P = contraction(F[:n], A[:n], S, dSt)
            c = P.sum(1)
            i = names.index(name)
            fwd.append(dd[:n, i])
            con.append(c)
            md, mc = dd[:n, i].mean(), c.mean()
            s_f, s_c = rel_sigma(dd[:n, i]), rel_sigma(c)
            T_f, T_a = fa["times"].sum(), aj["times"].sum()
            Pm = [contraction(F[:n], A[:n], S, dSt, m).sum(1).mean()
                  for m in (2, 4)]
            log(f"  {name}: forward {md:.4e} (sigma {s_f:.3f}),"
                f" contraction {mc:.4e} (sigma {s_c:.3f}; P1 share "
                f"{abs(P[:, 1].mean()) / abs(mc):.1e}; 2 cm and 4 cm voxels "
                f"change it by {abs(Pm[0] / mc - 1):.1e}, "
                f"{abs(Pm[1] / mc - 1):.1e}); |fwd/contr| - 1 "
                f"{abs(md) / abs(mc) - 1:+.4f}, phase diff "
                f"{np.degrees(np.angle(md / mc)):+.2f} deg; T fwd "
                f"{T_f:.0f}s adj {T_a:.0f}s; FOM adj/fwd "
                f"{(s_f**2 * T_f) / (s_c**2 * T_a):.2f}")
            gate.append((f"sum|w| adjoint {name} {f:g} Hz",
                         sumw_check(aj, f"adjoint {name}")))
        if fwd and ok_all:
            n = min(len(x) for x in fwd)
            ok, info = paired_correlated(np.stack([x[:n] for x in fwd], 1),
                                         np.stack([x[:n] for x in con], 1))
            t = info["t"]
            log(f"  reciprocity: {'PASS' if ok else 'FAIL'} paired t "
                + " ".join(f"{x:+.2f}" for x in t) + f", integral t "
                f"{info['t_int'][0]:+.2f} {info['t_int'][1]:+.2f}, "
                f"Hotelling p {info['hotelling'][2]:.3f}")
            gate.append((f"reciprocity {f:g} Hz", ok))
        else:
            gate.append((f"reciprocity {f:g} Hz", False))
        gate.append((f"sum|w| forward A {f:g} Hz", sumw_check(fa,
                                                                "forward A")))
        # ---- two meshes
        pb = HERE / f"stage9b_fwd_B_{f:g}Hz.npz"
        if pb.exists():
            # same seeds as A: the static phases are identical, so the
            # test is paired; delta d per static source particle (d0 at the
            # far detectors is too poorly sampled for delta d / d0)
            fb = np.load(pb)
            ddb = (fb["det"][..., 1] * Rg).sum(-1)
            dFb = fb["fission"][:, 1] / fb["fission"][:, 0].real
            n = min(R, ddb.shape[0])
            a = np.column_stack([dd, dF])[:n]
            b = np.column_stack([ddb, dFb])[:n]
            ok, info = paired_correlated(a, b)
            t = info["t"]
            log(f"  meshes A vs B (12 detectors' delta d and dF/F0): "
                f"{'PASS' if ok else 'FAIL'} max |paired t| "
                f"{np.abs(t).max():.2f}, integral {info['t_int'][0]:+.2f} "
                f"{info['t_int'][1]:+.2f}, Hotelling p "
                f"{info['hotelling'][2]:.3f}, pair correlation "
                f"{info['corr_pair']:.2f}; dF/F0 A {dF.mean():.4e} B "
                f"{dFb.mean():.4e}")
            gate.append((f"meshes {f:g} Hz", ok))
            gate.append((f"sum|w| forward B {f:g} Hz",
                         sumw_check(fb, "forward B")))
        else:
            log(f"  meshes: no {pb.name}")
            gate.append((f"meshes {f:g} Hz", False))
        # ---- logged
        g0 = ref_g0.g0(2 * np.pi * f, n5.BETA_EFF_K, n5.LAMBDA_K,
                       n5.LAMBDA_EFF) * DRHO_N5
        log(f"  logged: dF/F0 {dF.mean():.4e} against G0 drho_N5 {g0:.4e} "
            f"(ratio - 1 {abs(dF.mean() / g0) - 1:+.3f}, "
            f"{np.degrees(np.angle(dF.mean() / g0)):+.2f} deg)")
        # ratio of means: d0 at the far detectors is zero in some replicas
        for i, name in enumerate(names):
            m = dd[:, i].mean() / d0[:, i].mean()
            s0 = d0[:, i].std(ddof=1) / np.sqrt(R) / d0[:, i].mean()
            log(f"    {name}: dd/d0 {m:.4e} (sigma dd "
                f"{rel_sigma(dd[:, i]):.3f}"
                f", d0 {s0:.3f}); against G0 drho_N5 {abs(m / g0) - 1:+.3f}, "
                f"{np.degrees(np.angle(m / g0)):+.2f} deg")
        k = fa["static_k"][:, args.static_inactive:]
        km = np.array([x[x > 0].mean() for x in k])
        log(f"    static k (active static generations, replica means): "
            f"{km.mean():.5f} +- {km.std(ddof=1) / np.sqrt(km.size):.5f}; "
            f"against noise.keff {(km.mean() - KEFF) * 1e5:+.0f} pcm")
    n_ok = sum(ok for _, ok in gate)
    log(f"\nStage 9b: {n_ok}/{len(gate)} gate checks passed")
    for name, ok in gate:
        if not ok:
            log(f"  FAIL {name}")
    (HERE / "stage9b_log.md").write_text("\n".join(LOG) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["forward", "windows", "adjoint",
                                    "linearity", "analyze"])
    p.add_argument("--freq", type=float, default=1.0)
    p.add_argument("--mesh", choices=["A", "B"], default="A")
    p.add_argument("--wwmap", action="store_true")
    p.add_argument("--wwmap-from", default="stage9b_fwd_A_1Hz.npz")
    p.add_argument("--detectors", default=",".join(GATE_DETS))
    p.add_argument("--importance", choices=["rr", "forward"], default="rr")
    p.add_argument("--adj-cancel", choices=["A", "B", "none"], default="A")
    p.add_argument("--detector", default="arc0")
    p.add_argument("--ww", choices=["cadis", "flat", "none"],
                   default="cadis")
    p.add_argument("--reps", type=int, default=20)
    p.add_argument("--particles", type=int, default=10000)
    p.add_argument("--batches", type=int, default=25)
    p.add_argument("--nskip", type=int, default=5)
    p.add_argument("--static-inactive", type=int, default=100)
    p.add_argument("--adj-particles", type=int, default=10000)
    p.add_argument("--adj-batches", type=int, default=25)
    p.add_argument("--base-seed", type=int, default=1)
    args = p.parse_args()
    {"forward": run_forward, "windows": run_windows, "adjoint": run_adjoint,
     "linearity": run_linearity, "analyze": analyze}[args.mode](args)


if __name__ == "__main__":
    main()
