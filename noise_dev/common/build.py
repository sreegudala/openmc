"""Build OpenMC multigroup models from the reference-solver Material objects.

The same `noise_common.Material` instances drive both the reference solvers
(noise_dev/reference) and the OpenMC runs, so the data are identical.

Group convention: reference group index 0 is the *fastest* group, which is
also OpenMC's group 0 when the group edges are given in ascending order.
"""

import sys
from pathlib import Path

import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "reference"))

from noise_common import Material  # noqa: E402


def group_edges(G):
    """Ascending energy edges [eV]: group g occupies (G-1-g, G-g]."""
    return np.arange(G + 1, dtype=float)


def group_energy(g, G):
    """Representative energy of group g (for monoenergetic sources)."""
    return G - g - 0.5


def histogram_bins(mu_bar, n_bins):
    """Bin integrals of f(mu) = (1 + 3 mu_bar mu) / 2 over n equal mu bins."""
    mu = np.linspace(-1, 1, n_bins + 1)
    F = 0.5 * (mu + 1.5 * mu_bar * mu**2)
    return np.diff(F)


def make_library(materials, path, histogram=None):
    """Write an MGXS library with one macroscopic XSdata per material.

    Parameters
    ----------
    materials : dict of name -> Material
    path : str or Path
    histogram : tuple (mu_bar, n_bins), optional
        Write histogram angular scattering with a linearly anisotropic
        distribution of mean cosine mu_bar (same for every transfer).

    """
    G = next(iter(materials.values())).G
    groups = openmc.mgxs.EnergyGroups(group_edges(G))
    K = max(m.lam.size for m in materials.values())
    xsdata = []
    for name, m in materials.items():
        fissile = m.nu_fission.sum() > 0
        x = openmc.XSdata(name, groups,
                          num_delayed_groups=m.lam.size if fissile else 0)
        if histogram is None:
            x.order = 0
        else:
            x.scatter_format = "histogram"
            x.order = histogram[1]
        x.set_total(m.total)
        x.set_absorption(m.absorption)
        if histogram is None:
            x.set_scatter_matrix(m.scatter[:, :, None])
        else:
            mu_bar, nb = histogram
            x.set_scatter_matrix(m.scatter[:, :, None]
                                 * histogram_bins(mu_bar, nb)[None, None, :])
        x.set_inverse_velocity(m.inv_v)
        if fissile:
            if m.lam.size > 0:
                x.set_prompt_nu_fission(m.nu_fission_p)
                x.set_delayed_nu_fission(m.nu_fission_d)
                x.set_chi_prompt(m.chi_p)
                x.set_chi_delayed(m.chi_d)
                x.set_decay_rate(m.lam)
            else:
                x.set_nu_fission(m.nu_fission_p)
                x.set_chi(m.chi_p)
            # fission xs is optional; nu = 2.5 for scoring 'fission'
            x.set_fission(m.nu_fission / 2.5)
        xsdata.append(x)
    lib = openmc.MGXSLibrary(groups, num_delayed_groups=K)
    lib.add_xsdatas(xsdata)
    lib.export_to_hdf5(str(path))


def slab_model(regions, bc=("vacuum", "vacuum"), x0=0.0, transverse=10.0,
               library="mgxs.h5"):
    """1-D slab in x (reflective in y and z).

    Parameters
    ----------
    regions : list of (width, material_name)
    bc : boundary types at x-min and x-max

    Returns
    -------
    openmc.Model, list of x edges
    """
    mats = {}
    for _, name in regions:
        if name not in mats:
            m = openmc.Material(name=name)
            m.set_density("macro", 1.0)
            m.add_macroscopic(name)
            mats[name] = m
    edges = x0 + np.concatenate([[0.0], np.cumsum([w for w, _ in regions])])
    xs = [openmc.XPlane(e) for e in edges]
    xs[0].boundary_type = bc[0]
    xs[-1].boundary_type = bc[1]
    h = transverse / 2
    ylo = openmc.YPlane(-h, boundary_type="reflective")
    yhi = openmc.YPlane(h, boundary_type="reflective")
    zlo = openmc.ZPlane(-h, boundary_type="reflective")
    zhi = openmc.ZPlane(h, boundary_type="reflective")
    cells = []
    for i, (_, name) in enumerate(regions):
        c = openmc.Cell(fill=mats[name],
                        region=+xs[i] & -xs[i + 1] & +ylo & -yhi & +zlo & -zhi)
        cells.append(c)
    model = openmc.Model()
    model.geometry = openmc.Geometry(cells)
    model.materials = openmc.Materials(list(mats.values()))
    model.materials.cross_sections = str(library)
    model.settings.energy_mode = "multi-group"
    model.settings.run_mode = "fixed source"
    return model, edges


def box_source(xlo, xhi, group, G, transverse=10.0, strength=1.0, phase=0.0):
    h = transverse / 2
    src = openmc.IndependentSource(
        space=openmc.stats.Box((xlo, -h, -h), (xhi, h, h)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([group_energy(group, G)], [1.0]),
        strength=strength)
    src.phase = phase
    return src


def x_mesh(xlo, xhi, n, transverse=10.0):
    h = transverse / 2
    mesh = openmc.RegularMesh()
    mesh.lower_left = (xlo, -h, -h)
    mesh.upper_right = (xhi, h, h)
    mesh.dimension = (n, 1, 1)
    return mesh
