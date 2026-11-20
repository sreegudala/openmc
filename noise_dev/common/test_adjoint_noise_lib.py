"""Unit checks of adjoint_noise_lib (brief N2): no transport, no replicas.

    python -m pytest noise_dev/common/test_adjoint_noise_lib.py

"""

import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "reference"))

import adjoint_lib  # noqa: E402
from adjoint_noise_lib import (effective_operator, infinite_k,  # noqa: E402
                               make_adjoint_noise_library, production_matrices)
from build import make_library  # noqa: E402
sys.path.insert(0, str(HERE.parent / "stage7"))
from run_stage7_prep import n3_material, n3x_material  # noqa: E402

FIXTURE = HERE / "fixtures" / "adjoint_lib" / "forward.h5"


@pytest.fixture(scope="module", params=["fixture", "n3", "n3x"])
def forward(request, tmp_path_factory):
    """Forward libraries in two fission forms: chi + nu-fission + beta
    (fixture) and prompt/delayed vectors with chi-prompt/chi-delayed
    (noise_dev build.make_library)."""
    if request.param == "fixture":
        return FIXTURE
    path = tmp_path_factory.mktemp("lib") / f"{request.param}.h5"
    mat = n3_material() if request.param == "n3" else n3x_material()
    make_library({"m": mat}, path)
    return path


def _datasets(path):
    out = {}
    with h5py.File(path, "r") as f:
        G = int(f.attrs["energy_groups"])
        K = int(f.attrs["delayed_groups"])
        for name, group in f.items():
            if not isinstance(group, h5py.Group):
                continue
            nb, fmt = adjoint_lib._nbins(group)
            for key, data in group.items():
                if key == "kTs":
                    continue
                N, _ = adjoint_lib._dense(data["scatter_data"], G, nb)
                d = {"nu_scatter": N, "total": data["total"][()],
                     "inv_v": data["inverse-velocity"][()]}
                if group.attrs["fissionable"]:
                    d["Mp"], d["Md"] = production_matrices(data, G, K)
                    d["lam"] = data["decay-rate"][()]
                out[name, key] = d
    return out


def test_transpose_twice_is_identity(forward, tmp_path):
    adj = tmp_path / "adj.h5"
    twice = tmp_path / "twice.h5"
    make_adjoint_noise_library(forward, adj)
    make_adjoint_noise_library(adj, twice)
    a, b = _datasets(forward), _datasets(twice)
    assert a.keys() == b.keys()
    for key in a:
        assert a[key].keys() == b[key].keys(), key
        for name in a[key]:
            np.testing.assert_allclose(b[key][name], a[key][name],
                                       rtol=1e-14, atol=0, err_msg=str(key))


def test_k_equal(forward, tmp_path):
    adj = tmp_path / "adj.h5"
    report = make_adjoint_noise_library(forward, adj)
    assert report["fission_form"] == "matrix"
    for name in report["fissionable"]:
        kf = infinite_k(effective_operator(forward, name))
        ka = infinite_k(effective_operator(adj, name))
        assert ka == pytest.approx(kf, rel=1e-12)


def test_effective_operator_is_transpose(forward, tmp_path):
    """What OpenMC transports with the adjoint library is A^T exactly
    (scatter after N1's absorption rebalance and multiplicity, and
    production), including where the rebalance zeroes absorption."""
    adj = tmp_path / "adj.h5"
    report = make_adjoint_noise_library(forward, adj)
    assert report["max_multiplicity"] > 1.5
    with h5py.File(forward, "r") as f:
        names = [n for n, g in f.items() if isinstance(g, h5py.Group)]
    for name in names:
        of, oa = effective_operator(forward, name), effective_operator(adj,
                                                                       name)
        np.testing.assert_allclose(oa["total"], of["total"], rtol=1e-15)
        np.testing.assert_allclose(oa["scatter"], of["scatter"].T,
                                   rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(oa["Mp"], of["Mp"].T, rtol=1e-15)
        np.testing.assert_allclose(oa["Md"], of["Md"].transpose(0, 2, 1),
                                   rtol=1e-15)
        np.testing.assert_array_equal(oa["decay_rate"], of["decay_rate"])
        np.testing.assert_array_equal(oa["inv_v"], of["inv_v"])
