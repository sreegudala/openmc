"""C5G7 2-D noise problem of FEMFFUSION's examples/2D_C5G7_NN/NN1.2 as an
OpenMC multigroup model (Stage 9a).

All data are read from FEMFFUSION's own example files, so both codes use the
same numbers:
- input.xs.xml: 7-group C5G7 cross sections (9 mixes, 0-based ids). SigmaS
  is stored as S[to, from] (its column sums are the out-scatter).
- c5g7.prec: speeds, 8 delayed families (beta, lambda, delayed spectra) for
  the materials listed in PrecMat (1-based: mat k is mix k - 1).
- input.geom.xml: 51 x 51 pin cells of pitch 1.26 cm (3 x 3 assemblies:
  UO2, MOX / MOX, UO2, and water reflector), fuel radius 0.54 cm; pin types
  0 (water box) to 7 (the perturbed pin). Row i of the layout is y index i
  (y from the reflective side); FEMFFUSION's output grid has the same
  orientation (its relative noise peaks at row 15, column 18, the
  perturbed pin). Reflective at x = 0 and y = 0, vacuum at 64.26 cm.
- input.dxs.xml: the perturbation of mix 8, dSigma_t = dSigma_a =
  -0.05 i Sigma_a (5 % absorption, i.e. capture, since fission and
  scattering are unchanged) in one MOX 7.0 % pin at 1 Hz (c5g7.dyn.prm).

Quirk of the example, reproduced: mix 8 has MOX 7.0 % cross sections but
its kinetics (speeds, beta, delayed spectra; .prec mat 9) are UO2's.

Prompt spectrum as in FEMFFUSION (Materials::parse_prec):
chi_p = (chi - sum_k beta_k chi_d,k) / (1 - beta), so the total spectrum of
the static problem is chi. Mixes without precursors (the fission chamber)
have prompt fission only.

Group convention (common/build.py): data group g (0 = fastest) occupies
OpenMC energy (G-1-g, G-g] eV, so OpenMC energy bin b is data group 6 - b.
"""

import os
import re
from pathlib import Path

import numpy as np
import openmc

G = 7
K = 8
PITCH = 1.26
RADIUS = 0.54
N = 51
SIDE = N * PITCH          # 64.26
HALF_Z = 1.0
EXAMPLE = "examples/2D_C5G7_NN/NN1.2"


def femffusion_dir():
    d = os.environ.get("FEMFFUSION_DIR")
    if d:
        return Path(d)
    src = os.environ.get("OPENMC_DEV_SRC")
    if src and (Path(src).parent / "FEMFFUSION").is_dir():
        return Path(src).parent / "FEMFFUSION"
    return Path.home() / "dev" / "external" / "FEMFFUSION"


def example_dir():
    d = femffusion_dir() / EXAMPLE
    if not d.is_dir():
        raise SystemExit(f"FEMFFUSION example not found: {d} "
                         "(set FEMFFUSION_DIR)")
    return d


def _numbers(text):
    return np.array([float(x) for x in re.findall(
        r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", text)])


def _mixes(path, complex_values=False):
    s = Path(path).read_text()
    out = {}
    for m in re.finditer(r'<mix id="(\d+)".*?<name>(.*?)</name>(.*?)</mix>',
                         s, re.S):
        d = {"name": m.group(2).strip()}
        for t in re.finditer(r"<(\w+)>(.*?)</\1>", m.group(3), re.S):
            v = _numbers(t.group(2))
            if complex_values:
                v = v[0::2] + 1j * v[1::2]
            d[t.group(1)] = v
        out[int(m.group(1))] = d
    return out


def read_xs(d=None):
    """Cross sections per mix id: total, absorption, scatter[from, to],
    nu_fission, chi (data group order, 0 fastest)."""
    d = d or example_dir()
    xs = {}
    for i, m in _mixes(d / "input.xs.xml").items():
        S = m["SigmaS"].reshape(G, G)                 # [to, from]
        assert np.allclose(m["SigmaA"] + S.sum(0), m["SigmaT"], atol=1e-5)
        xs[i] = {"name": m["name"], "total": m["SigmaT"],
                 "absorption": m["SigmaA"], "scatter": S.T.copy(),
                 "nu_fission": m["NuSigF"], "chi": m["Chi"]}
    return xs


def read_prec(d=None):
    """Kinetics: speeds per mix (cm/s), and per precursor mix beta (K),
    lambda (K), delayed spectra chi_d (K, G)."""
    d = d or example_dir()
    lines = (d / "c5g7.prec").read_text().splitlines()
    rows = lambda i, n: [_numbers(l) for l in lines[i:i + n]
                         if not l.lstrip().startswith("#")]
    data = [l for l in lines]
    i = next(j for j, l in enumerate(data) if l.startswith("Velocities"))
    vel = [r for r in rows(i + 1, 12) if r.size > 1][:G]
    vel = np.array([r[1:] for r in vel])               # (G, n_mats)
    pm = next(l for l in data if l.startswith("PrecMat"))
    prec_mats = [int(x) - 1 for x in pm.split()[1:]]   # mix ids
    def table(key):
        j = next(j for j, l in enumerate(data) if l.startswith(key))
        t = [r for r in rows(j + 1, 12) if r.size > 1][:K]
        return np.array([r[1:] for r in t])             # (K, n_prec)
    beta = table("Delayed_Fractions")
    lam = table("Delayed_Decay_Constants")
    chi_d = {}
    for n, mix in enumerate(prec_mats):
        j = next(j for j, l in enumerate(data)
                 if l.strip() == f"Material {mix + 1}")
        t = [r for r in rows(j + 1, 10) if r.size > 1][:G]
        chi_d[mix] = np.array([r[1:] for r in t]).T     # (K, G)
    out = {"speed": {m: vel[:, m] for m in range(vel.shape[1])},
           "prec": {}}
    for n, mix in enumerate(prec_mats):
        out["prec"][mix] = {"beta": beta[:, n], "lam": lam[:, n],
                            "chi_d": chi_d[mix]}
    return out


def read_layout(d=None):
    """(51, 51) array of pin types, [row = y index, column = x index]."""
    d = d or example_dir()
    s = (d / "input.geom.xml").read_text()
    p = s[s.index('<plane id="0">'):]
    c = re.search(r"<components>(.*?)</components>", p, re.S).group(1)
    c = re.sub(r"<!--.*?-->", "", c, flags=re.S)
    rows = [r.split() for r in c.split(";") if r.split()]
    lay = np.array(rows, dtype=int)
    assert lay.shape == (N, N)
    return lay


def read_pins(d=None):
    """Pin type -> (inner mix, outer mix) (box pins: (mix, mix))."""
    d = d or example_dir()
    s = (d / "input.geom.xml").read_text()
    pins = {}
    for m in re.finditer(r'<pin id="(\d+)">(.*?)</pin>', s, re.S):
        mats = [int(x) for x in re.search(r"<materials>(.*?)</materials>",
                                          m.group(2)).group(1).split()]
        pins[int(m.group(1))] = (mats[0], mats[-1])
    return pins


def read_perturbation(d=None):
    """{mix: relative dSigma_t / Sigma_t per group (complex)}; checks that
    only Sigma_t = Sigma_a change."""
    d = d or example_dir()
    xs = read_xs(d)
    out = {}
    for i, m in _mixes(d / "input.dxs.xml", complex_values=True).items():
        dt = m["SigmaT"]
        if not np.any(dt != 0):
            continue
        assert np.allclose(m["SigmaA"], dt)
        for key in ("SigmaS", "NuSigF", "SigF", "Nu", "Chi"):
            assert not np.any(m[key] != 0), (i, key)
        out[i] = dt / xs[i]["total"]
    return out


def read_frequency(d=None):
    d = d or example_dir()
    s = (d / "c5g7.dyn.prm").read_text()
    return float(re.search(r"Frequency\s*=\s*([-\d.eE+]+)", s).group(1))


def library(path, d=None):
    """Write the 9-mix MGXS library (names 'mix0' ... 'mix8')."""
    d = d or example_dir()
    xs, kin = read_xs(d), read_prec(d)
    # XSdata group 0 is the highest energy, as data group 0
    groups = openmc.mgxs.EnergyGroups(np.arange(G + 1, dtype=float))
    xsdata = []
    for i, x in xs.items():
        fissile = x["nu_fission"].sum() > 0
        prec = kin["prec"].get(i)
        data = openmc.XSdata(f"mix{i}", groups,
                             num_delayed_groups=K if fissile else 0)
        data.order = 0
        data.set_total(x["total"])
        data.set_absorption(x["absorption"])
        data.set_scatter_matrix(x["scatter"][:, :, None])
        data.set_inverse_velocity(1.0 / kin["speed"][i])
        if fissile:
            if prec is not None:
                b, chi_d = prec["beta"], prec["chi_d"]
                chi_p = (x["chi"] - b @ chi_d) / (1 - b.sum())
                assert np.all(chi_p >= -1e-12), (i, chi_p)
                data.set_prompt_nu_fission((1 - b.sum()) * x["nu_fission"])
                data.set_delayed_nu_fission(b[:, None] * x["nu_fission"])
                data.set_chi_prompt(np.clip(chi_p, 0, None))
                data.set_chi_delayed(chi_d)
                data.set_decay_rate(prec["lam"])
            else:
                # no precursors (fission chamber): prompt only. OpenMC
                # needs the delayed arrays in a delayed-group library
                lam = next(iter(kin["prec"].values()))["lam"]
                data.set_prompt_nu_fission(x["nu_fission"])
                data.set_delayed_nu_fission(np.zeros((K, G)))
                data.set_chi_prompt(x["chi"])
                data.set_chi_delayed(np.tile(x["chi"], (K, 1)))
                data.set_decay_rate(lam)
            data.set_fission(x["nu_fission"] / 2.5)
        xsdata.append(data)
    lib = openmc.MGXSLibrary(groups, num_delayed_groups=K)
    lib.add_xsdatas(xsdata)
    lib.export_to_hdf5(str(path))


def model(lib, particles, batches, keff=1.0, run_mode="noise", inactive=100,
          eta=1.0, nskip=20, static_inactive=200, cancellation=None,
          d=None):
    """OpenMC model. run_mode 'eigenvalue' (k and the static flux) or
    'noise' (fixed source, the perturbation sampled as the noise source).
    cancellation: n for an n x n x 1 mesh over the core (51 = pin cells).
    The flux tally 'flux' is on the 51 x 51 pin-cell mesh, 7 groups."""
    d = d or example_dir()
    lay, pins = read_layout(d), read_pins(d)
    mats = {}
    for i in range(9):
        m = openmc.Material(name=f"mix{i}")
        m.set_density("macro", 1.0)
        m.add_macroscopic(f"mix{i}")
        mats[i] = m

    cyl = openmc.ZCylinder(r=RADIUS)
    universes, pert_cells = {}, []
    for t, (inner, outer) in pins.items():
        if inner == outer:
            cells = [openmc.Cell(name=f"pin{t} box", fill=mats[inner])]
        else:
            fuel = openmc.Cell(name=f"pin{t} inner", fill=mats[inner],
                               region=-cyl)
            cells = [fuel, openmc.Cell(name=f"pin{t} outer",
                                       fill=mats[outer], region=+cyl)]
        universes[t] = openmc.Universe(name=f"pin{t}", cells=cells)
    pert = read_perturbation(d)
    for t, (inner, _) in pins.items():
        if inner in pert:
            pert_cells.append(universes[t].cells[
                next(iter(universes[t].cells))])
    assert len(pert) == 1 and len(pert_cells) == 1
    assert (lay == next(t for t, (i, _) in pins.items()
                        if i in pert)).sum() == 1

    lat = openmc.RectLattice(name="core")
    lat.pitch = (PITCH, PITCH)
    lat.lower_left = (0.0, 0.0)
    lat.universes = [[universes[lay[N - 1 - r, c]] for c in range(N)]
                     for r in range(N)]          # OpenMC rows: top first
    x0 = openmc.XPlane(0.0, boundary_type="reflective")
    x1 = openmc.XPlane(SIDE, boundary_type="vacuum")
    y0 = openmc.YPlane(0.0, boundary_type="reflective")
    y1 = openmc.YPlane(SIDE, boundary_type="vacuum")
    z0 = openmc.ZPlane(-HALF_Z, boundary_type="reflective")
    z1 = openmc.ZPlane(HALF_Z, boundary_type="reflective")
    root = openmc.Cell(name="core", fill=lat,
                       region=+x0 & -x1 & +y0 & -y1 & +z0 & -z1)

    mdl = openmc.Model()
    mdl.geometry = openmc.Geometry([root])
    mdl.materials = openmc.Materials(list(mats.values()))
    mdl.materials.cross_sections = str(lib)
    s = mdl.settings
    s.energy_mode = "multi-group"
    s.particles = particles
    s.batches = batches
    fuel_box = openmc.stats.Box((0, 0, -HALF_Z), (2 * 17 * PITCH,
                                                  2 * 17 * PITCH, HALF_Z))
    s.source = [openmc.IndependentSource(
        space=fuel_box, angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([G - 0.5], [1.0]),
        constraints={"fissionable": True})]

    mesh = openmc.RegularMesh()
    mesh.lower_left = (0.0, 0.0, -HALF_Z)
    mesh.upper_right = (SIDE, SIDE, HALF_Z)
    mesh.dimension = (N, N, 1)
    efilter = openmc.EnergyFilter(np.arange(G + 1, dtype=float))
    t = openmc.Tally(name="flux")
    t.scores = ["flux"]
    if run_mode == "eigenvalue":
        s.run_mode = "eigenvalue"
        s.inactive = inactive
        t.filters = [openmc.MeshFilter(mesh), efilter]
        mdl.tallies = openmc.Tallies([t])
        return mdl

    s.run_mode = "fixed source"
    (mix, amp), = pert.items()
    noise = {"frequency": read_frequency(d), "eta": eta, "keff": keff,
             "nskip": nskip, "static_inactive": static_inactive,
             "perturbations": [{"cells": pert_cells,
                                "total": [complex(a) for a in amp]}]}
    if cancellation:
        cm = openmc.RegularMesh()
        cm.lower_left = (0.0, 0.0, -HALF_Z)
        cm.upper_right = (SIDE, SIDE, HALF_Z)
        cm.dimension = (cancellation, cancellation, 1)
        noise["cancellation_mesh"] = cm
    s.noise = noise
    t.filters = [openmc.MeshFilter(mesh), efilter,
                 openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
    mdl.tallies = openmc.Tallies([t])
    return mdl


def to_grid(v):
    """(..., 51 * 51 * 7) tally values (mesh bin ix + 51 iy, then energy
    bin) to (..., 7 data groups, 51 rows = y index, 51 columns = x index);
    energy bin b is data group 6 - b."""
    v = np.asarray(v)
    v = v.reshape(*v.shape[:-1], N, N, G)
    return np.moveaxis(v, -1, -3)[..., ::-1, :, :]
