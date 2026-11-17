"""Stage 1 gates: complex weights, complex tallies and covariance.

1a identity       noise off vs noise on at w=0 (same seed): Re bins equal the
                  baseline bit for bit in sum and sum_sq; Im bins exactly 0.
1b phase rotation fission off, source phase pi/3: Re = cos * baseline,
                  Im = sin * baseline to 1e-12 relative.
1c covariance     correlation of 1b equals 1 to 1e-10 and sigma_phi = 0; two
                  sources at phases 0 and pi/2: accumulated covariance equals
                  the per-batch-statepoint covariance to 1e-12.
1d plumbing       source-file round trip with wgt_im (Python and C++ writers,
                  old-format files still read), 1 vs 8 threads, event vs
                  history.

Problem P1: 2-group slab 40 cm, CORTEX fuel data, vacuum boundaries,
isotropic group-1 source in the central 2 cm (k_inf of the fuel is 0.62, so
fission stays on in 1a).
"""

import argparse
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import openmc
import openmc.noise

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent / "reference"))

from build import box_source, make_library, slab_model, x_mesh  # noqa
from ref_inf import cortex_fuel  # noqa

WORK = HERE / "work"
RESULTS = []


def record(name, ok, detail):
    RESULTS.append((name, ok, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def p1_model(fission=True, noise=None, complex_tally=False, phase=0.0,
             survival=False, particles=20000, batches=10, sources=None):
    model, _ = slab_model([(40.0, "fuel")], x0=-20.0)
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.seed = 12345
    model.settings.create_fission_neutrons = fission
    model.settings.survival_biasing = survival
    if survival:
        model.settings.cutoff = {"weight": 0.25, "weight_avg": 0.5}
    model.settings.source = sources or [box_source(-1, 1, 0, 2, phase=phase)]
    if noise is not None:
        model.settings.noise = noise
    mesh = x_mesh(-20, 20, 40)
    filters = [openmc.MeshFilter(mesh), openmc.EnergyFilter([0., 1., 2.])]
    if complex_tally:
        filters.append(openmc.ComplexPartFilter())
    scores = ["flux", "total", "absorption", "scatter", "fission",
              "nu-fission"]
    t1 = openmc.Tally(name="tl")
    t1.filters = filters
    t1.scores = scores
    t1.estimator = "tracklength"
    t2 = openmc.Tally(name="col")
    t2.filters = filters
    t2.scores = scores
    t2.estimator = "collision"
    model.tallies = openmc.Tallies([t1, t2])
    return model


def run(model, name, threads=1, event=False, **kw):
    d = WORK / name
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    model.materials.cross_sections = str(WORK / "mgxs.h5")
    if event:
        model.settings.event_based = True
    sp = model.run(cwd=d, threads=threads, output=False, **kw)
    return openmc.StatePoint(sp)


def split_parts(tally):
    """Return (sum_re, sum_im, sq_re, sq_im) reshaped to the baseline shape
    (complex filter is the last filter, stride 1)."""
    s = tally.sum
    q = tally.sum_sq
    return s[0::2], s[1::2], q[0::2], q[1::2]


def test_1a():
    print("Test 1a: identity at w = 0")
    for surv in (False, True):
        base = run(p1_model(survival=surv), f"1a_base_s{int(surv)}")
        noise = run(p1_model(survival=surv, noise={"frequency": 0.0},
                             complex_tally=True), f"1a_noise_s{int(surv)}")
        for tname in ("tl", "col"):
            b = base.get_tally(name=tname)
            n = noise.get_tally(name=tname)
            sre, sim, qre, qim = split_parts(n)
            same = (np.array_equal(sre, b.sum) and np.array_equal(qre, b.sum_sq))
            zero = (not np.any(sim)) and (not np.any(qim))
            nz = np.count_nonzero(b.sum)
            record(f"1a {tname} survival={surv} Re bit-identical", same,
                   f"{nz} non-zero bins; max |diff| "
                   f"{np.abs(sre - b.sum).max():.1e}")
            record(f"1a {tname} survival={surv} Im exactly zero", zero,
                   f"max |Im| {np.abs(sim).max():.1e}")


def test_1b_1c():
    print("Test 1b/1c: phase rotation and covariance")
    th = np.pi / 3
    for surv in (False, True):
        base = run(p1_model(fission=False, survival=surv), f"1b_base_s{int(surv)}")
        rot = run(p1_model(fission=False, survival=surv,
                           noise={"frequency": 0.0}, complex_tally=True,
                           phase=th), f"1b_rot_s{int(surv)}")
        for tname in ("tl", "col"):
            b = base.get_tally(name=tname)
            n = rot.get_tally(name=tname)
            sre, sim, _, _ = split_parts(n)
            m = b.sum != 0
            ere = np.abs(sre[m] / (np.cos(th) * b.sum[m]) - 1).max()
            eim = np.abs(sim[m] / (np.sin(th) * b.sum[m]) - 1).max()
            zero_ok = not np.any(sre[~m]) and not np.any(sim[~m])
            record(f"1b {tname} survival={surv} Re = cos*base", ere < 1e-12
                   and zero_ok, f"max rel err {ere:.1e}")
            record(f"1b {tname} survival={surv} Im = sin*base", eim < 1e-12,
                   f"max rel err {eim:.1e}")
            st = openmc.noise.complex_tally_stats(n)
            mm = st["amplitude"] > 0
            cerr = np.abs(st["corr"][mm] - 1).max()
            record(f"1c {tname} survival={surv} corr = 1", cerr < 1e-10,
                   f"max |corr - 1| {cerr:.1e}")
            # "Zero to rounding": a variance formed from running sums has a
            # floor of ~eps * mean^2, so sigma_phi cannot fall below about
            # sqrt(eps / (n - 1)) rad even for perfectly correlated parts.
            n_r = n.num_realizations
            floor = 10 * np.sqrt(np.finfo(float).eps / (n_r - 1))
            smax = st["phase_std"][mm].max()
            record(f"1c {tname} survival={surv} sigma_phi zero to rounding",
                   smax < floor, f"max sigma_phi {smax:.1e} rad (rounding "
                   f"floor 10*sqrt(eps/(n-1)) = {floor:.1e})")

    # Two sources, phases 0 and pi/2, in different regions
    srcs = [box_source(-10, -8, 0, 2, phase=0.0),
            box_source(6, 9, 1, 2, phase=np.pi / 2)]
    nb = 12
    model = p1_model(fission=False, noise={"frequency": 0.0},
                     complex_tally=True, sources=srcs, batches=nb)
    model.settings.statepoint = {"batches": list(range(1, nb + 1))}
    d = WORK / "1c_two"
    sp = run(model, "1c_two")
    for tname in ("tl", "col"):
        sums = []
        for k in range(1, nb + 1):
            with openmc.StatePoint(d / f"statepoint.{k:02d}.h5") as spk:
                sums.append(spk.get_tally(name=tname).sum.copy())
        vals = np.diff(np.array([np.zeros_like(sums[0])] + sums), axis=0)
        re, im = vals[:, 0::2], vals[:, 1::2]
        n = nb
        cov_direct = ((re * im).sum(0) / n - re.mean(0) * im.mean(0)) / (n - 1)
        st = openmc.noise.complex_tally_stats(sp.get_tally(name=tname))
        scale = np.sqrt(np.abs(st["re_std"] * st["im_std"])).max() ** 2
        err = np.abs(st["cov"] - cov_direct).max() / scale
        # also the raw cross sum
        cross = sp.get_tally(name=tname).sum_cross[0::2]
        cross_direct = (re * im).sum(0)
        cerr = np.abs(cross - cross_direct).max() / np.abs(cross_direct).max()
        record(f"1c two-source {tname} covariance vs per-batch", err < 1e-12
               and cerr < 1e-12, f"cov err/scale {err:.1e}, cross rel err "
               f"{cerr:.1e}; corr range [{np.nanmin(st['corr']):.2f}, "
               f"{np.nanmax(st['corr']):.2f}]")


def test_1d():
    print("Test 1d: plumbing")
    th = 0.7
    rng = np.random.default_rng(3)
    N = 5000
    pts = []
    for i in range(N):
        mu = 2 * rng.random() - 1
        ph = 2 * np.pi * rng.random()
        u = (mu, np.sqrt(1 - mu * mu) * np.cos(ph), np.sqrt(1 - mu * mu) * np.sin(ph))
        r = (2 * rng.random() - 1, 10 * rng.random() - 5, 10 * rng.random() - 5)
        pts.append((r, u))
    WORK.mkdir(exist_ok=True)
    # complex weights in the file
    pl_c = openmc.ParticleList([openmc.SourceParticle(
        r=r, u=u, E=1.5, wgt=np.cos(th), wgt_im=np.sin(th)) for r, u in pts])
    pl_c.export_to_hdf5(WORK / "src_complex.h5")
    back = openmc.ParticleList.from_hdf5(WORK / "src_complex.h5")
    ok = all(a.wgt == b.wgt and a.wgt_im == b.wgt_im and a.r == b.r
             for a, b in zip(pl_c, back))
    with h5py.File(WORK / "src_complex.h5") as f:
        names = f["source_bank"].dtype.names
    record("1d Python source file round trip with wgt_im", ok and
           "wgt_im" in names, f"fields {names}")
    # real file + source phase
    pl_r = openmc.ParticleList([openmc.SourceParticle(r=r, u=u, E=1.5)
                                for r, u in pts])
    pl_r.export_to_hdf5(WORK / "src_real.h5")
    with h5py.File(WORK / "src_real.h5") as f:
        names = f["source_bank"].dtype.names
    record("1d real-only file keeps legacy format", "wgt_im" not in names,
           f"fields {names}")

    def file_model(fname, phase):
        m = p1_model(fission=False, noise={"frequency": 0.0},
                     complex_tally=True, particles=N, batches=1)
        src = openmc.FileSource(str(WORK / fname))
        src.phase = phase
        m.settings.source = [src]
        return m

    a = run(file_model("src_complex.h5", 0.0), "1d_file_complex")
    b = run(file_model("src_real.h5", th), "1d_file_phase")
    ta, tb = a.get_tally(name="tl"), b.get_tally(name="tl")
    err = np.abs(ta.sum - tb.sum).max() / np.abs(tb.sum).max()
    record("1d C++ reads wgt_im from file (== real file with phase)",
           err < 1e-12, f"max rel diff {err:.1e}")
    # legacy file in non-noise mode still runs
    m = p1_model(fission=False, particles=N, batches=1)
    m.settings.source = [openmc.FileSource(str(WORK / "src_real.h5"))]
    c = run(m, "1d_file_legacy")
    tc = c.get_tally(name="tl")
    err = np.abs(tc.sum - tb.sum[0::2] / np.cos(th)).max() / np.abs(tc.sum).max()
    record("1d legacy file in non-noise mode", err < 1e-12,
           f"max rel diff vs rotated {err:.1e}")
    # C++ writer: surface source written in noise mode carries wgt_im, and
    # weights keep the source phase (wgt_im / wgt = tan theta)
    m = p1_model(fission=False, survival=True, noise={"frequency": 0.0},
                 complex_tally=True, phase=th, particles=4000, batches=1)
    surf = [s for s in m.geometry.get_all_surfaces().values()
            if isinstance(s, openmc.XPlane) and s.x0 == 20.0][0]
    m.settings.surf_source_write = {"surface_ids": [surf.id],
                                    "max_particles": 10000}
    run(m, "1d_write_surf")
    fname = WORK / "1d_write_surf" / "surface_source.h5"
    with h5py.File(fname) as f:
        names = f["source_bank"].dtype.names
    pl = openmc.ParticleList.from_hdf5(fname)
    w = np.array([[p.wgt, p.wgt_im] for p in pl])
    perr = np.abs(np.arctan2(w[:, 1], w[:, 0]) - th).max()
    record("1d C++ written source file carries wgt_im", "wgt_im" in names
           and len(pl) > 0 and perr < 1e-12,
           f"{len(pl)} sites, max |phase - theta| {perr:.1e}, fields {names}")

    # threads and event-based
    th3 = np.pi / 3

    def mk():
        return p1_model(fission=True, survival=True, noise={"frequency": 0.0},
                        complex_tally=True, phase=0.0)
    r1 = run(mk(), "1d_t1", threads=1).get_tally(name="tl").sum
    r8 = run(mk(), "1d_t8", threads=8).get_tally(name="tl").sum
    ev = run(mk(), "1d_ev", threads=8, event=True).get_tally(name="tl").sum
    e8 = np.abs(r8 - r1).max() / np.abs(r1).max()
    eev = np.abs(ev - r1).max() / np.abs(r1).max()
    record("1d 1 vs 8 threads (fission on)", e8 < 1e-12, f"max rel diff {e8:.1e}")
    record("1d event vs history", eev < 1e-12, f"max rel diff {eev:.1e}")
    mk2 = lambda: p1_model(fission=False, survival=True, noise={"frequency": 0.0},
                           complex_tally=True, phase=th3)
    r1 = run(mk2(), "1d_t1c", threads=1).get_tally(name="tl").sum
    r8 = run(mk2(), "1d_t8c", threads=8).get_tally(name="tl").sum
    ev = run(mk2(), "1d_evc", threads=8, event=True).get_tally(name="tl").sum
    e8 = np.abs(r8 - r1).max() / np.abs(r1).max()
    eev = np.abs(ev - r1).max() / np.abs(r1).max()
    record("1d 1 vs 8 threads (complex weights)", e8 < 1e-12,
           f"max rel diff {e8:.1e}")
    record("1d event vs history (complex weights)", eev < 1e-12,
           f"max rel diff {eev:.1e}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("tests", nargs="*", default=["1a", "1b", "1d"])
    args = p.parse_args()
    WORK.mkdir(exist_ok=True)
    make_library({"fuel": cortex_fuel()}, WORK / "mgxs.h5")
    if "1a" in args.tests:
        test_1a()
    if "1b" in args.tests:
        test_1b_1c()
    if "1d" in args.tests:
        test_1d()
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\nStage 1: {len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    with open(HERE / "stage1_log.md", "w") as f:
        f.write("# Stage 1 results\n\n| check | result | detail |\n|---|---|---|\n")
        for name, ok, detail in RESULTS:
            f.write(f"| {name} | {'PASS' if ok else 'FAIL'} | {detail} |\n")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
