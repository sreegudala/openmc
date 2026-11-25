"""OpenMC on the 1-D vibrating plate (run from noise_dev/; 10 replicas,
0.125 cm cancellation). Writes stage8/diag/work_plate/vals.npy."""
import sys
sys.path[:0] = ["stage8/diag", "stage8", "common", "reference", "stage4"]; sys.argv = sys.argv[:1]
import numpy as np, openmc
from pathlib import Path
import plate as P
import run_stage8 as S
from build import make_library, slab_model, box_source, x_mesh
from harness import run_replicas
W = Path("stage8/diag/work_plate"); W.mkdir(exist_ok=True)
lib = (W / "lib.h5").resolve(); make_library({"mix": P.mix, "water": P.wat}, lib)
names = ["mix", "water", "water", "mix", "mix", "mix", "water", "water", "mix"]
regions = [(w, nm) for (w, _, _), nm in zip(P.REG, names)]
sl, k, _ = P.sn(1)
def make():
    model, edges = slab_model(regions, x0=-50.0, library=lib)
    cells = list(model.geometry.root_universe.cells.values()); mats = {m.name: m for m in model.materials}
    s = model.settings; s.particles = 20000; s.batches = 60; s.source = [box_source(-50, 50, 0, 2)]
    s.noise = {"frequency": 1.0, "eta": 1.0, "keff": k, "nskip": 20, "static_inactive": 200,
               "cancellation_mesh": x_mesh(-50, 50, 800),
               "perturbations": [
                   {"cells": cells[2:4], "vibration": {"axis": "x", "position": P.F1, "amplitude": P.EPS, "harmonic": 1,
                                                        "negative_material": mats["water"], "positive_material": mats["mix"]}},
                   {"cells": cells[5:7], "vibration": {"axis": "x", "position": P.F2, "amplitude": P.EPS, "harmonic": 1,
                                                        "negative_material": mats["mix"], "positive_material": mats["water"]}}]}
    t = openmc.Tally(name="tl"); t.filters = [openmc.MeshFilter(x_mesh(-50, 50, 800)), openmc.EnergyFilter([0., 1., 2.]),
                                               openmc.NoisePhaseFilter(), openmc.ComplexPartFilter()]
    t.scores = ["flux"]; model.tallies = openmc.Tallies([t]); return model
out = run_replicas(make, W / "reps", n_rep=10, base_seed=11)
v = out["values"]["tl"].reshape(10, 800, 2, 2)
np.save(W / "vals.npy", v); print("T", out["times"].sum())
