"""Regenerate the sync fixture for common/adjoint_lib.py.

Writes a random 4-group library (forward.h5) with 8-bin histogram
scattering, upscatter, a non-unit multiplicity in one material and a
fissionable material with two delayed families, then transposes it with
make_adjoint_library imported from the ce_mg_triga repository itself
(adjoint_step3a.h5). test_adjoint_lib.py checks that the copy in
noise_dev/common gives an identical file.

The brief (noise_dev/briefs/noise_stage7_brief.md, N1) asks for the fixture
to come from the TRIGA repository; until it ships one, this script makes it
from that repository's code at the commit recorded in adjoint_lib.py. Rerun
it whenever adjoint_lib.py is refreshed:

    python make_fixture.py --triga ~/TRIGA_proposal_work/ce_mg_triga

"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import openmc
import openmc.mgxs

HERE = Path(__file__).resolve().parent
G = 4
N_BINS = 8


def random_library(path, seed=20261120):
    rng = np.random.default_rng(seed)
    groups = openmc.mgxs.EnergyGroups(np.array([0.0, 1.0, 1e3, 1e5, 2e7]))
    edges_mu = np.linspace(-1, 1, N_BINS + 1)
    xsdata = []
    for name, fissile, mult in (("fuel", True, False), ("mod", False, True)):
        x = openmc.XSdata(name, groups, num_delayed_groups=2 if fissile else 0)
        x.scatter_format = "histogram"
        x.order = N_BINS
        total = rng.uniform(0.3, 1.5, G)
        # full transfer matrix (downscatter and upscatter), random mu shape
        S = rng.uniform(0.0, 1.0, (G, G)) * np.where(
            np.arange(G)[None, :] >= np.arange(G)[:, None] - 1, 1.0, 0.1)
        frac = rng.uniform(0.6, 0.95, G)      # scatter fraction of total
        S *= (frac * total / S.sum(1))[:, None]
        shape = rng.uniform(0.2, 1.0, (G, G, N_BINS)) \
            * (1 + 0.5 * 0.5 * (edges_mu[:-1] + edges_mu[1:]))
        shape /= shape.sum(2, keepdims=True)
        M = np.ones((G, G))
        if mult:
            M[0, 0] = M[0, 1] = 1.3                # (n,2n)-like transfers
        x.set_total(total)
        x.set_absorption(total - S.sum(1))
        x.set_scatter_matrix(S[:, :, None] * M[:, :, None] * shape)
        if mult:
            x.set_multiplicity_matrix(M)
        x.set_inverse_velocity(1.0 / rng.uniform(2e5, 2e9, G))
        if fissile:
            nuf = rng.uniform(0.05, 0.4, G) * (total - S.sum(1))
            x.set_nu_fission(nuf)
            x.set_chi(np.array([0.7, 0.29, 0.01, 0.0]))
            x.set_beta(np.array([0.002, 0.004]))
            x.set_decay_rate(np.array([0.08, 1.2]))
            x.set_fission(nuf / 2.4)
        xsdata.append(x)
    lib = openmc.MGXSLibrary(groups, num_delayed_groups=2)
    lib.add_xsdatas(xsdata)
    lib.export_to_hdf5(str(path))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--triga", type=Path, required=True,
                   help="ce_mg_triga checkout (its step3a_adjoint.py is used)")
    a = p.parse_args()
    triga = a.triga.resolve()
    sys.path.insert(0, str(triga))
    import step3a_adjoint                     # noqa: E402 (TRIGA repo code)

    commit = subprocess.run(["git", "-C", str(triga), "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    random_library(HERE / "forward.h5")
    report = step3a_adjoint.make_adjoint_library(HERE / "forward.h5",
                                                 HERE / "adjoint_step3a.h5")
    (HERE / "SOURCE").write_text(f"ce_mg_triga {commit}\n{report}\n")
    print(f"fixture from ce_mg_triga {commit}: {report}")


if __name__ == "__main__":
    main()
