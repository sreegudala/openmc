"""Sync test: the copied adjoint_lib.py reproduces ce_mg_triga's step3a.

The fixture (fixtures/adjoint_lib/) holds a random 4-group library and its
transpose made by ce_mg_triga's own make_adjoint_library (commit in
fixtures/adjoint_lib/SOURCE). Run with

    python -m pytest noise_dev/common/test_adjoint_lib.py

"""

import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import adjoint_lib  # noqa: E402

FIXTURE = HERE / "fixtures" / "adjoint_lib"


def _items(f):
    """{path: (kind, value)} for every attribute and dataset in a file."""
    out = {}

    def visit(name, obj):
        for k, v in obj.attrs.items():
            out[f"{name}@{k}"] = v
        if isinstance(obj, h5py.Dataset):
            out[name] = obj[()]

    for k, v in f.attrs.items():
        out[f"/@{k}"] = v
    f.visititems(visit)
    return out


def assert_same_library(a, b):
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        ia, ib = _items(fa), _items(fb)
    assert sorted(ia) == sorted(ib)
    for k in ia:
        assert np.array_equal(np.asarray(ia[k]), np.asarray(ib[k])), k


def test_copy_matches_step3a(tmp_path):
    out = tmp_path / "adjoint.h5"
    report = adjoint_lib.make_adjoint_library(FIXTURE / "forward.h5", out)
    assert report["materials"] == 2
    assert_same_library(out, FIXTURE / "adjoint_step3a.h5")
