"""CORTEX 2-D UOX assembly (design doc Section 3.4) as an OpenMC model.

Reflective 17x17 assembly, pitch 1.26 cm, square fuel pins of side 0.7314
cm, 25 water holes, 0.08 cm water blade, 2 groups, 1 precursor family. The
geometry and data follow Abeille's noise_oscillation.yaml. The perturbed
pin is the box x in [2.1543, 2.8857], y in [-2.8857, -2.1543] (lattice
column 10 from the left, row 6 from the bottom).

Group convention: the library from build.make_library puts data group g in
energy (G-1-g, G-g], so OpenMC energy bin 0 = (0, 1] eV is the thermal group
(data group 2) and bin 1 the fast group. Abeille's bins are ascending with
group index, so its bin 0 is the fast group.
"""

import sys
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "reference"))

from build import box_source, make_library  # noqa
from ref_inf import cortex_fuel, cortex_water  # noqa

KEFF = 0.99916          # Abeille noise input
PITCH = 1.26
HALF_PIN = 0.3657
HALF_LAT = 17 * PITCH / 2   # 10.71
HALF_BOX = 10.79
PERTURBED = (10, 6)     # (column from left, row from bottom)
EPS = {"total": 0.041, "fission": 0.021, "scatter": 0.034}

WATER_HOLES = [
    "11111111111111111",
    "11111111111111111",
    "11111211211211111",
    "11121111111112111",
    "11111111111111111",
    "11211211211211211",
    "11111111111111111",
    "11111111111111111",
    "11211211211211211",
    "11111111111111111",
    "11111111111111111",
    "11211211211211211",
    "11111111111111111",
    "11121111111112111",
    "11111211211211111",
    "11111111111111111",
    "11111111111111111",
]


def library(path):
    make_library({"fuel": cortex_fuel(), "water": cortex_water()}, path)


def model(lib, particles, batches, frequency=1.0, eta=1.0, nskip=20,
          static_inactive=200, amplitudes=None, cancellation=None,
          mesh_shape=(138, 138), run_mode="noise"):
    """Build the CORTEX model.

    amplitudes : dict of complex relative amplitudes dSigma_x / Sigma_x for
        'total', 'scatter', 'fission' (default: Abeille's square oscillation
        at its fundamental, -i pi epsilon_x).
    cancellation : int n for an n x n x 1 cancellation mesh, or None.
    run_mode : 'noise' or 'eigenvalue' (static check).
    """
    mats = {}
    for name in ("fuel", "water"):
        m = openmc.Material(name=name)
        m.set_density("macro", 1.0)
        m.add_macroscopic(name)
        mats[name] = m

    xl, xh = openmc.XPlane(-HALF_PIN), openmc.XPlane(HALF_PIN)
    yl, yh = openmc.YPlane(-HALF_PIN), openmc.YPlane(HALF_PIN)
    pin_box = +xl & -xh & +yl & -yh

    def pin(label):
        fuel = openmc.Cell(name=f"{label} fuel", fill=mats["fuel"],
                           region=pin_box)
        water = openmc.Cell(name=f"{label} water", fill=mats["water"],
                            region=~pin_box)
        return openmc.Universe(cells=[fuel, water]), fuel

    u_pin, _ = pin("pin")
    u_pert, perturbed_cell = pin("perturbed pin")
    u_hole = openmc.Universe(cells=[openmc.Cell(name="water hole",
                                                fill=mats["water"])])
    u_outer = openmc.Universe(cells=[openmc.Cell(name="blade",
                                                 fill=mats["water"])])
    lat = openmc.RectLattice(name="assembly")
    lat.pitch = (PITCH, PITCH)
    lat.lower_left = (-HALF_LAT, -HALF_LAT)
    rows = []
    for r, line in enumerate(WATER_HOLES):        # r = row from the top
        row = []
        for c, ch in enumerate(line):
            if ch == "2":
                row.append(u_hole)
            elif (c, 16 - r) == PERTURBED:
                row.append(u_pert)
            else:
                row.append(u_pin)
        rows.append(row)
    lat.universes = rows
    lat.outer = u_outer

    box = openmc.model.RectangularPrism(2 * HALF_BOX, 2 * HALF_BOX,
                                        boundary_type="reflective")
    zl = openmc.ZPlane(-2.0, boundary_type="reflective")
    zh = openmc.ZPlane(2.0, boundary_type="reflective")
    root = openmc.Cell(name="root", fill=lat, region=-box & +zl & -zh)

    mdl = openmc.Model()
    mdl.geometry = openmc.Geometry([root])
    mdl.materials = openmc.Materials(list(mats.values()))
    mdl.materials.cross_sections = str(lib)
    s = mdl.settings
    s.energy_mode = "multi-group"
    s.particles = particles
    s.batches = batches
    src = openmc.IndependentSource(
        space=openmc.stats.Box((-HALF_LAT, -HALF_LAT, -2), (HALF_LAT,
                                                            HALF_LAT, 2)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([1.5], [1.0]))   # fast group
    s.source = [src]

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-HALF_BOX, -HALF_BOX, -2)
    mesh.upper_right = (HALF_BOX, HALF_BOX, 2)
    mesh.dimension = (*mesh_shape, 1)
    efilter = openmc.EnergyFilter([0.0, 1.0, 2.0])
    if run_mode == "eigenvalue":
        s.run_mode = "eigenvalue"
        s.inactive = static_inactive
        t = openmc.Tally(name="flux")
        t.filters = [openmc.MeshFilter(mesh), efilter]
        t.scores = ["flux"]
        mdl.tallies = openmc.Tallies([t])
        return mdl

    s.run_mode = "fixed source"
    if amplitudes is None:
        amplitudes = {k: -1j * np.pi * v for k, v in EPS.items()}
    noise = {"frequency": frequency, "eta": eta, "keff": KEFF,
             "nskip": nskip, "static_inactive": static_inactive,
             "perturbations": [{"cells": [perturbed_cell], **amplitudes}]}
    if cancellation:
        cm = openmc.RegularMesh()
        cm.lower_left = (-HALF_BOX, -HALF_BOX, -2)
        cm.upper_right = (HALF_BOX, HALF_BOX, 2)
        cm.dimension = (cancellation, cancellation, 1)
        noise["cancellation_mesh"] = cm
    s.noise = noise
    t = openmc.Tally(name="flux")
    t.filters = [openmc.MeshFilter(mesh), efilter, openmc.NoisePhaseFilter(),
                 openmc.ComplexPartFilter()]
    t.scores = ["flux"]
    mdl.tallies = openmc.Tallies([t])
    return mdl
