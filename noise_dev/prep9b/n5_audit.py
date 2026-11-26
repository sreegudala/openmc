"""N5: audit of the regenerated ce_mg_triga CASMO-70 library (design doc,
Stage 9b preparation; brief noise_stage7_brief.md, N5).

Reads OpenMC MGXS library files in place (no copy) and writes a text
summary, so the library itself can stay where it is:

    python prep9b/n5_audit.py [-o OUT.md] LIB.h5 [LIB2.h5 ...]\n        (default OUT: prep9b/n5_audit_log.md)

Per library: group structure, delayed groups, and per XSdata (material) and
temperature:
- which datasets are present: total, absorption, scatter, multiplicity,
  nu-fission (vector or matrix form), prompt/delayed nu-fission, chi,
  chi-prompt, chi-delayed, beta, decay-rate, inverse-velocity;
- the delayed data N5 asks for: beta per family and its sum, lambda_k,
  chi_d normalisation, inverse velocity range;
- whether the multiplicity varies along a scatter-matrix row (OpenMC's MG
  sampling counts such a row's multiplicity twice; STATUS, open items), and
  the largest multiplicity;
- the balance sum_g' sigma_s(g -> g') + sigma_a against sigma_t.

Nothing here judges the k or kinetics checks; those need the TRIGA model.
"""

import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
LOG = []


def log(msg=""):
    print(msg, flush=True)
    LOG.append(msg)


def arr(g, name):
    return np.array(g[name][()]) if name in g else None


def fmt(v, n=8):
    v = np.ravel(v)
    s = ", ".join(f"{x:.4g}" for x in v[:n])
    return f"[{s}{', ...' if v.size > n else ''}] ({v.size})"


def audit_temperature(name, tname, t, G, K):
    present = sorted(t.keys())
    log(f"  {name} / {tname}: datasets {present}")
    tot = arr(t, "total")
    absn = arr(t, "absorption")
    # scatter: OpenMC stores scatter_data/scatter_matrix (flattened, with
    # g_min/g_max) and optionally multiplicity_matrix
    sd = t.get("scatter_data")
    if sd is not None:
        keys = sorted(sd.keys())
        log(f"    scatter_data: {keys}; attrs "
            f"{ {k: sd.attrs[k] for k in sd.attrs} }")
        if "g_min" in sd and "g_max" in sd and "scatter_matrix" in sd:
            gmin, gmax = sd["g_min"][()], sd["g_max"][()]
            flat = sd["scatter_matrix"][()]
            n_ang = flat.size // int(np.sum(gmax - gmin + 1))
            fmt_attr = t.parent.attrs.get("scatter_format", b"legendre")
            hist = (fmt_attr.decode() if isinstance(fmt_attr, bytes)
                    else str(fmt_attr)) == "histogram"
            # P0-equivalent (histogram: sum over bins; Legendre: order 0)
            S = np.zeros((G, G))
            mult = np.ones((G, G))
            mflat = sd["multiplicity_matrix"][()] \
                if "multiplicity_matrix" in sd else None
            i = j = 0
            for gin in range(G):
                for gout in range(gmin[gin] - 1, gmax[gin]):
                    block = flat[i:i + n_ang]
                    S[gin, gout] = block.sum() if hist else block[0]
                    if mflat is not None:
                        mult[gin, gout] = mflat[j]
                    i += n_ang
                    j += 1
            log(f"    scatter: {n_ang} angular entries per transfer")
            if mflat is not None:
                varies = []
                for gin in range(G):
                    row = mult[gin, gmin[gin] - 1:gmax[gin]]
                    w = S[gin, gmin[gin] - 1:gmax[gin]]
                    row = row[w > 0]
                    if row.size and np.ptp(row) > 1e-6:
                        varies.append(gin)
                log(f"    multiplicity: max {mult.max():.4f}; varies along "
                    f"rows {varies if varies else 'none'}")
            if tot is not None and absn is not None:
                out = S.sum(axis=1)    # nu-scatter row sums
                if mflat is not None:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        sc = np.where(mult > 0, S / mult, 0).sum(axis=1)
                else:
                    sc = out
                with np.errstate(divide="ignore", invalid="ignore"):
                    bal = np.abs(sc + absn - tot) / tot
                log(f"    balance |sum scatter + absorption - total| / total:"
                    f" max {np.nanmax(bal):.2e}")
                for g in np.where(~(bal <= 1e-3))[0][:10]:
                    log(f"      group {g}: total {tot[g]:.4e}, absorption "
                        f"{absn[g]:.4e}, scatter {sc[g]:.4e}, nu-scatter "
                        f"{out[g]:.4e}")
                if mflat is not None:
                    for gin in range(G):
                        row = mult[gin, gmin[gin] - 1:gmax[gin]]
                        w = S[gin, gmin[gin] - 1:gmax[gin]]
                        if w.sum() > 0 and np.ptp(row[w > 0]) > 1e-6:
                            mbar = (w * row).sum() / w.sum()
                            log(f"      row {gin}: multiplicity "
                                f"{row[w > 0].min():.3f}-{row[w > 0].max():.3f}"
                                f", nu-scatter-weighted mean {mbar:.4f}")
    for key in ("inverse-velocity", "decay-rate", "beta"):
        v = arr(t, key)
        if v is not None:
            log(f"    {key}: shape {v.shape}, {fmt(v)}")
            if key == "beta":
                b = v if v.ndim == 1 else v.sum(axis=-1) if \
                    v.shape[-1] == G else v
                log(f"      beta sum {np.sum(b) if b.ndim == 1 else b.sum(0)}")
    for key in ("nu-fission", "prompt-nu-fission", "delayed-nu-fission",
                "chi", "chi-prompt", "chi-delayed", "fission",
                "kappa-fission"):
        v = arr(t, key)
        if v is not None:
            form = "matrix" if v.ndim >= 2 and v.shape[-1] == G and \
                v.shape[-2] == G else "vector"
            extra = ""
            if key.startswith("chi"):
                extra = f", sums {fmt(v.sum(axis=-1))}"
            log(f"    {key}: shape {v.shape} ({form}){extra}")


def audit(path):
    log(f"\n## {path}")
    p = Path(path)
    log(f"size {p.stat().st_size / 1e6:.1f} MB")
    with h5py.File(path, "r") as f:
        attrs = {k: f.attrs[k] for k in f.attrs}
        G = int(attrs.get("energy_groups", 0))
        K = int(attrs.get("delayed_groups", 0))
        log(f"file attrs: energy_groups {G}, delayed_groups {K}, "
            f"version {attrs.get('version')}")
        if "group structure" in attrs:
            e = np.array(attrs["group structure"])
            log(f"group structure: {e.size} edges, {e.min():.4g} to "
                f"{e.max():.4g} eV")
        names = list(f.keys())
        log(f"{len(names)} XSdata: {names}")
        for name in names:
            g = f[name]
            a = {k: (v.decode() if isinstance(v, bytes) else v)
                 for k, v in g.attrs.items()}
            log(f"\n### {name}: attrs {a}")
            temps = [k for k in g.keys() if k.endswith("K")]
            log(f"  temperatures {temps}")
            for tname in temps[:1] if len(temps) > 1 else temps:
                audit_temperature(name, tname, g[tname], G, K)
            if len(temps) > 1:
                log(f"  (only {temps[0]} shown of {len(temps)})")


def main():
    args = sys.argv[1:]
    out = HERE / "n5_audit_log.md"
    if len(args) >= 2 and args[0] == "-o":
        out, args = Path(args[1]), args[2:]
    if not args:
        raise SystemExit(__doc__)
    log("# N5 library audit")
    for path in args:
        audit(path)
    out.write_text("\n".join(LOG) + "\n")
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
