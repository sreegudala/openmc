"""Regenerate the sync fixture for common/pt_lib.py.

Takes pt_terms from ce_mg_triga's step5e_pool_pt.py at the commit recorded
in pt_lib.py (read with ``git show``, so the TRIGA working tree does not
matter), evaluates it on random inputs and writes inputs and outputs to
sync.npz, with the commit in SOURCE. test_pt_lib.py checks that the copy
in noise_dev/common gives identical outputs. Rerun whenever pt_lib.py is
refreshed:

    python make_fixture.py --triga ~/TRIGA_proposal_work/ce_mg_triga

"""

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
COPY = HERE.parent.parent / "pt_lib.py"


def triga_pt_terms(triga, commit):
    src = subprocess.run(["git", "-C", str(triga), "show",
                          f"{commit}:step5e_pool_pt.py"], check=True,
                         capture_output=True, text=True).stdout
    m = re.search(r"^def pt_terms\(.*?(?=^def )", src, re.S | re.M)
    ns = {"np": np}
    exec(m.group(0), ns)
    return ns["pt_terms"]


def cases(seed=20261120):
    """Random inputs covering both shapes of dSt and dN (shared and per
    voxel) and a voxel with weight 0."""
    rng = np.random.default_rng(seed)
    V, G = 5, 4
    out = []
    for per_voxel in (False, True):
        F = rng.random((V, G, 4)) - 0.3
        A = rng.random((V, G, 4)) - 0.3
        w = rng.random(V)
        w[2] = 0.0
        shape = (V, G) if per_voxel else (G,)
        dSt = rng.random(shape) - 0.5
        dN0 = rng.random(shape + (G,)) - 0.5
        dN1 = rng.random(shape + (G,)) - 0.5
        out.append({"F": F, "A": A, "w": w, "dSt": dSt, "dN0": dN0,
                    "dN1": dN1})
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--triga", required=True)
    args = p.parse_args()
    commit = re.search(r"Source commit: ([0-9a-f]{40})",
                       COPY.read_text()).group(1)
    fn = triga_pt_terms(args.triga, commit)
    data = {}
    for i, c in enumerate(cases()):
        for k, v in c.items():
            data[f"{i}_{k}"] = v
        data[f"{i}_out"] = np.array(fn(**c))
    np.savez(HERE / "sync.npz", **data)
    (HERE / "SOURCE").write_text(f"ce_mg_triga {commit} step5e_pool_pt.py "
                                 "pt_terms\n")
    print(f"wrote sync.npz from {commit}")


if __name__ == "__main__":
    main()
