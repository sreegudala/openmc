"""Unit tests for the frequency-domain noise Python API."""

import numpy as np
import pytest
import openmc
import openmc.noise
import lxml.etree as ET


def test_settings_noise_roundtrip(run_in_tmpdir):
    s = openmc.Settings()
    s.noise = {'frequency': 1.5, 'eta': 2.0, 'keff': 0.99916,
               'max_generations': 5000}
    s.export_to_xml()
    s2 = openmc.Settings.from_xml()
    assert s2.noise == s.noise
    assert isinstance(s2.noise['max_generations'], int)


def test_settings_noise_cancellation_mesh(run_in_tmpdir):
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-1., -1., -1.)
    mesh.upper_right = (1., 1., 1.)
    mesh.dimension = (4, 4, 1)
    s = openmc.Settings()
    s.noise = {'frequency': 10.0, 'cancellation_mesh': mesh,
               'cancellation_group': False, 'cancellation_roulette': False,
               'growth_limit': 20.0}
    s.export_to_xml()
    s2 = openmc.Settings.from_xml()
    assert s2.noise['cancellation_mesh'].id == mesh.id
    assert tuple(s2.noise['cancellation_mesh'].dimension) == (4, 4, 1)
    assert s2.noise['cancellation_group'] is False
    assert s2.noise['cancellation_roulette'] is False
    assert s2.noise['growth_limit'] == 20.0


@pytest.mark.parametrize('noise', [
    {'frequency': -1.0}, {'eta': 0.0}, {'keff': -1.0},
    {'max_generations': 0}, {'bogus': 1.0}, {'growth_limit': 1.0},
    {'cancellation_mesh': 3}, {'cancellation_group': 'yes'},
    {'cancellation_roulette': 1}])
def test_settings_noise_invalid(noise):
    s = openmc.Settings()
    with pytest.raises((ValueError, TypeError)):
        s.noise = noise


def test_source_phase_roundtrip():
    src = openmc.IndependentSource(strength=2.0)
    assert src.phase == 0.0
    elem = src.to_xml_element()
    assert elem.get('phase') is None
    src.phase = np.pi / 3
    elem = src.to_xml_element()
    new = openmc.SourceBase.from_xml_element(elem)
    assert new.phase == pytest.approx(np.pi / 3)
    with pytest.raises(TypeError):
        src.phase = 'x'


def test_file_source_phase(tmp_path):
    src = openmc.FileSource(tmp_path / 'source.h5')
    src.phase = 0.5
    new = openmc.SourceBase.from_xml_element(src.to_xml_element())
    assert isinstance(new, openmc.FileSource)
    assert new.phase == pytest.approx(0.5)


def test_complex_part_filter_xml():
    f = openmc.ComplexPartFilter(filter_id=7)
    assert list(f.bins) == ['real', 'imag']
    assert f.num_bins == 2
    elem = f.to_xml_element()
    assert elem.get('type') == 'complexpart'
    new = openmc.Filter.from_xml_element(elem)
    assert isinstance(new, openmc.ComplexPartFilter)
    assert new.id == 7
    with pytest.raises(ValueError):
        openmc.ComplexPartFilter(bins=['a', 'b'])


def test_source_particle_wgt_im_file(tmp_path):
    parts = [openmc.SourceParticle(E=1.0, wgt=0.5, wgt_im=-0.25),
             openmc.SourceParticle(E=2.0, wgt=1.0)]
    openmc.write_source_file(parts, tmp_path / 'c.h5')
    back = openmc.ParticleList.from_hdf5(tmp_path / 'c.h5')
    assert [p.wgt_im for p in back] == [-0.25, 0.0]
    assert [p.wgt for p in back] == [0.5, 1.0]
    # real-only lists keep the legacy layout
    openmc.write_source_file(parts[1:], tmp_path / 'r.h5')
    back = openmc.ParticleList.from_hdf5(tmp_path / 'r.h5')
    assert back[0].wgt_im == 0.0


class _FakeTally:
    """Minimal stand-in for openmc.Tally used by complex_tally_stats."""
    def __init__(self, filters, s, s2, sc, n):
        self.id = 1
        self.filters = filters
        self.sum = s
        self.sum_sq = s2
        self.sum_cross = sc
        self.num_realizations = n


def test_complex_tally_stats():
    rng = np.random.default_rng(0)
    n = 50
    # two mesh-like bins x complex part (complex filter last)
    re = rng.normal([1.0, -2.0], 0.1, size=(n, 2))
    im = 0.5 * re + rng.normal(0, 0.05, size=(n, 2))
    vals = np.empty((n, 4))
    vals[:, 0::2] = re
    vals[:, 1::2] = im
    s = vals.sum(0)[:, None, None]
    s2 = (vals**2).sum(0)[:, None, None]
    cross = np.empty(4)
    cross[0::2] = (re * im).sum(0)
    cross[1::2] = (re * im).sum(0)
    filters = [openmc.EnergyFilter([0., 1., 2.]), openmc.ComplexPartFilter()]
    t = _FakeTally(filters, s, s2, cross[:, None, None], n)
    st = openmc.noise.complex_tally_stats(t)
    assert st['re'].squeeze() == pytest.approx(re.mean(0))
    assert st['im'].squeeze() == pytest.approx(im.mean(0))
    cov = np.array([np.cov(re[:, i], im[:, i])[0, 1] / n for i in range(2)])
    assert st['cov'].squeeze() == pytest.approx(cov)
    assert st['re_std'].squeeze() == pytest.approx(re.std(0, ddof=1) / np.sqrt(n))
    # first-order propagation for the amplitude
    x, y = re.mean(0), im.mean(0)
    vx, vy = re.var(0, ddof=1) / n, im.var(0, ddof=1) / n
    A = np.hypot(x, y)
    sA = np.sqrt((x**2 * vx + y**2 * vy + 2 * x * y * cov)) / A
    assert st['amplitude_std'].squeeze() == pytest.approx(sA)
    assert st['phase'].squeeze() == pytest.approx(np.arctan2(y, x))


def test_complex_tally_stats_requires_filter():
    t = _FakeTally([openmc.EnergyFilter([0., 1.])], np.zeros((1, 1, 1)),
                   np.zeros((1, 1, 1)), None, 2)
    with pytest.raises(ValueError):
        openmc.noise.complex_tally_stats(t)


def test_settings_noise_perturbations(run_in_tmpdir):
    cell = openmc.Cell(cell_id=11)
    s = openmc.Settings()
    s.noise = {'frequency': 1.0, 'nskip': 3, 'static_inactive': 10,
               'perturbations': [
                   {'cells': [cell, 12], 'total': 1e-3 - 2e-4j},
                   {'cells': [13], 'scatter': [0.01, 0.02j],
                    'fission': [0.5, -0.25 + 1j]}]}
    s.export_to_xml()
    s2 = openmc.Settings.from_xml()
    n = s2.noise
    assert n['nskip'] == 3 and isinstance(n['nskip'], int)
    assert n['static_inactive'] == 10
    p0, p1 = n['perturbations']
    assert p0['cells'] == [11, 12]
    assert p0['total'] == pytest.approx([1e-3 - 2e-4j])
    assert p1['scatter'] == pytest.approx([0.01, 0.02j])
    assert p1['fission'] == pytest.approx([0.5, -0.25 + 1j])


@pytest.mark.parametrize('noise', [
    {'nskip': 0}, {'static_inactive': -1},
    {'perturbations': [{'total': 1e-3}]},
    {'perturbations': [{'cells': [1], 'bogus': 1.0}]},
    {'perturbations': [{'cells': ['a']}]}])
def test_settings_noise_perturbations_invalid(noise):
    s = openmc.Settings()
    with pytest.raises((ValueError, TypeError)):
        s.noise = noise


def test_noise_phase_filter_xml():
    f = openmc.NoisePhaseFilter(filter_id=8)
    assert list(f.bins) == ['static', 'noise']
    elem = f.to_xml_element()
    assert elem.get('type') == 'noisephase'
    new = openmc.Filter.from_xml_element(elem)
    assert isinstance(new, openmc.NoisePhaseFilter)
    assert new.id == 8
    with pytest.raises(ValueError):
        openmc.NoisePhaseFilter(bins=['a', 'b'])


def _ww_model(windows=True, phase=2.0, generator=False):
    """One-group reflective box (infinite medium), Sigma_t = 1, Sigma_a = 0.5,
    1/v = 1e-5 s/cm at omega = 5e4 rad/s, so omega/v = 0.5 and the volume
    integrated flux per source particle is exp(i phase) / (0.5 + 0.5 i)."""
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    xs = openmc.XSdata('m', groups)
    xs.order = 0
    xs.set_total([1.0])
    xs.set_absorption([0.5])
    xs.set_scatter_matrix([[[0.5]]])
    xs.set_inverse_velocity([1.0e-5])
    lib = openmc.MGXSLibrary(groups)
    lib.add_xsdata(xs)
    lib.export_to_hdf5('mgxs.h5')

    mat = openmc.Material()
    mat.add_macroscopic('m')
    mats = openmc.Materials([mat])
    mats.cross_sections = 'mgxs.h5'
    box = openmc.model.RectangularParallelepiped(
        -5, 5, -5, 5, -5, 5, boundary_type='reflective')
    geom = openmc.Geometry([openmc.Cell(fill=mat, region=-box)])

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-5, -5, -5)
    mesh.upper_right = (5, 5, 5)
    mesh.dimension = (1, 1, 1)
    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    settings.particles = 2000
    settings.batches = 10
    src = openmc.IndependentSource(
        space=openmc.stats.Point(), energy=openmc.stats.Discrete([1e6], [1]))
    src.phase = phase
    settings.source = src
    settings.noise = {'frequency': 5.0e4 / (2 * np.pi)}
    if windows:
        # Tight window: unit source weights are split (upper bound 0.1) and
        # survivors below 0.05 are rouletted at every collision
        settings.weight_windows = openmc.WeightWindows(
            mesh, lower_ww_bounds=[0.05], upper_bound_ratio=2.0,
            max_split=10)
    if generator:
        settings.weight_window_generators = openmc.WeightWindowGenerator(mesh)
    tally = openmc.Tally(name='flux')
    tally.filters = [openmc.ComplexPartFilter()]
    tally.scores = ['flux']
    tally.estimator = 'tracklength'
    return openmc.Model(geom, mats, settings, openmc.Tallies([tally]))


@pytest.mark.parametrize('windows', [False, True])
def test_noise_weight_windows_complex(run_in_tmpdir, windows):
    # The source phase makes Re w negative at birth: windows acting on Re w
    # would kill every particle (cutoff) and splits would lose Im w
    phase = 2.0
    sp = _ww_model(windows, phase).run()
    with openmc.StatePoint(sp) as s:
        st = openmc.noise.complex_tally_stats(s.get_tally(name='flux'))
    ref = np.exp(1j * phase) / (0.5 + 0.5j)
    v = st['value'].ravel()[0]
    assert abs(v.real - ref.real) < 4 * st['re_std'].ravel()[0]
    assert abs(v.imag - ref.imag) < 4 * st['im_std'].ravel()[0]
    assert abs(abs(v) / abs(ref) - 1) < 0.03


def test_noise_weight_windows_birth_roulette_threads(run_in_tmpdir):
    # Absolute windows above the unit source weight roulette every source
    # particle at birth, before its first step sets the random number
    # stream; with several threads a fresh Particle's stream was undefined
    # (segfault in prn, Stage 9b smoke test)
    phase = 2.0
    model = _ww_model(windows=True, phase=phase)
    ww = model.settings.weight_windows[0]
    ww.lower_ww_bounds = [2.0]
    ww.upper_ww_bounds = [4.0]
    sp = model.run(threads=4)
    with openmc.StatePoint(sp) as s:
        st = openmc.noise.complex_tally_stats(s.get_tally(name='flux'))
    ref = np.exp(1j * phase) / (0.5 + 0.5j)
    v = st['value'].ravel()[0]
    assert abs(v.real - ref.real) < 4 * st['re_std'].ravel()[0]
    assert abs(v.imag - ref.imag) < 4 * st['im_std'].ravel()[0]


def test_noise_refuses_weight_window_generators(run_in_tmpdir):
    # Generators tally a real flux; noise mode must refuse them
    with pytest.raises(RuntimeError, match='[Gg]enerators'):
        _ww_model(windows=False, generator=True).run()


def test_settings_noise_vibration_roundtrip(run_in_tmpdir):
    fuel = openmc.Material(material_id=21)
    s = openmc.Settings()
    vib = {'axis': 'y', 'position': -2.1543, 'amplitude': 0.2,
           'harmonic': 2, 'negative_material': fuel,
           'positive_material': 22}
    s.noise = {'frequency': 2.0,
               'perturbations': [{'cells': [5, 6], 'vibration': vib}]}
    s.export_to_xml()
    p = openmc.Settings.from_xml().noise['perturbations'][0]
    assert p['cells'] == [5, 6]
    assert p['vibration'] == {'axis': 'y', 'position': -2.1543,
                              'amplitude': 0.2, 'harmonic': 2,
                              'negative_material': 21,
                              'positive_material': 22}


@pytest.mark.parametrize('vib', [
    {'axis': 'w', 'position': 0.0, 'amplitude': 0.2, 'harmonic': 1,
     'negative_material': 1, 'positive_material': 2},
    {'axis': 'x', 'position': 0.0, 'amplitude': -0.2, 'harmonic': 1,
     'negative_material': 1, 'positive_material': 2},
    {'axis': 'x', 'position': 0.0, 'amplitude': 0.2, 'harmonic': 0,
     'negative_material': 1, 'positive_material': 2},
    {'axis': 'x', 'position': 0.0, 'amplitude': 0.2, 'harmonic': 1}])
def test_settings_noise_vibration_invalid(vib):
    s = openmc.Settings()
    with pytest.raises((ValueError, TypeError)):
        s.noise = {'perturbations': [{'cells': [1], 'vibration': vib}]}
    with pytest.raises(ValueError):
        s.noise = {'perturbations': [{'cells': [1], 'total': 0.1,
                                      'vibration': {}}]}


def _vibration_model(harmonic):
    """Homogeneous infinite medium of M (k_inf = 1, noise keff 1.5); the
    layer cells |x| < 0.2 carry a vibration between two other materials L
    (fissile) and R, so the sampled source has a closed-form expectation."""
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    lib = openmc.MGXSLibrary(groups)
    data = {'M': (1.0, 0.5, 0.5), 'L': (1.2, 0.4, 0.9), 'R': (0.7, 0.6, 0.0)}
    for name, (tot, scat, nuf) in data.items():
        xs = openmc.XSdata(name, groups)
        xs.order = 0
        xs.set_total([tot])
        xs.set_absorption([tot - scat])
        xs.set_scatter_matrix([[[scat]]])
        xs.set_inverse_velocity([1.0e-5])
        if nuf > 0:
            xs.set_nu_fission([nuf])
            xs.set_chi([1.0])
        lib.add_xsdata(xs)
    lib.export_to_hdf5('mgxs.h5')
    mats = {}
    for name in data:
        m = openmc.Material(name=name)
        m.add_macroscopic(name)
        mats[name] = m
    materials = openmc.Materials(mats.values())
    materials.cross_sections = 'mgxs.h5'
    planes = [openmc.XPlane(x) for x in (-5, -0.2, 0.0, 0.2, 5)]
    planes[0].boundary_type = planes[-1].boundary_type = 'reflective'
    yz = openmc.model.RectangularPrism(10, 10, axis='x',
                                       boundary_type='reflective')
    cells = [openmc.Cell(fill=mats['M'], region=+a & -b & -yz)
             for a, b in zip(planes[:-1], planes[1:])]
    # L and R must fill cells to have data: a separate reflective box that
    # no particle reaches
    far = [openmc.XPlane(x) for x in (20, 21, 22)]
    far[0].boundary_type = far[-1].boundary_type = 'reflective'
    side = [openmc.Cell(fill=mats[k], region=+a & -b & -yz)
            for k, a, b in (('L', far[0], far[1]), ('R', far[1], far[2]))]
    geom = openmc.Geometry(cells + side)
    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    settings.particles = 4000
    settings.batches = 10
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-5, -5, -5), (5, 5, 5)),
        energy=openmc.stats.Discrete([1e6], [1]))
    # Without delayed neutrons the |w| fission game at keff = k_inf is
    # divergent for any omega > 0; keff = 1.5 (subcritical noise problem)
    # and omega / v = 0.25 give a generation ratio of about 0.83
    settings.noise = {
        'frequency': 2.5e4 / (2 * np.pi), 'keff': 1.5, 'nskip': 3,
        'static_inactive': 10,
        'perturbations': [{'cells': cells[1:3], 'vibration': {
            'axis': 'x', 'position': 0.0, 'amplitude': 0.2,
            'harmonic': harmonic, 'negative_material': mats['L'],
            'positive_material': mats['R']}}]}
    tally = openmc.Tally(name='layer')
    tally.filters = [openmc.CellFilter(cells[1:3]), openmc.NoisePhaseFilter(),
                     openmc.ComplexPartFilter()]
    tally.scores = ['flux']
    return openmc.Model(geom, materials, settings, openmc.Tallies([tally]))


@pytest.mark.parametrize('harmonic', [1, 2])
def test_noise_vibration_source_integral(run_in_tmpdir, harmonic):
    # Sum of sampled source weights per static particle against
    # int c_n du / 2 x (-dSigma_t + dSigma_s + dnuSigma_f / k) x Phi_layer,
    # with int c_1 du = -i pi and int c_2 du = 0
    import h5py
    sp = _vibration_model(harmonic).run()
    with openmc.StatePoint(sp) as s:
        st = openmc.noise.complex_tally_stats(s.get_tally(name='layer'))
        n = s.n_particles
    phi = st['re'].reshape(2, 2)[:, 0].sum()          # static bin, both cells
    s_phi = np.hypot(*st['re_std'].reshape(2, 2)[:, 0])
    with h5py.File(sp, 'r') as f:
        ss = f['noise/source_sum'][()].reshape(-1, 2) / n
    src = ss[:, 0] + 1j * ss[:, 1]
    bracket = -(1.2 - 0.7) + (0.4 - 0.6) + 0.9 / 1.5
    ref = (-1j * np.pi / 2 if harmonic == 1 else 0.0) * bracket * phi
    s_src = np.array([src.real.std(ddof=1), src.imag.std(ddof=1)]) \
        / np.sqrt(src.size)
    s_ref = abs(bracket * np.pi / 2) * s_phi if harmonic == 1 else 0.0
    m = src.mean()
    # c_2 is real, so Im is exactly 0 on both sides for n = 2. The net source
    # is a small difference of large terms (about 9 % sigma here; a run with
    # 20000 x 40 gives 1.016 +- 0.020 of the expectation).
    assert abs(m.real - ref.real) <= 4 * np.hypot(s_src[0], s_ref)
    assert abs(m.imag - ref.imag) <= 4 * np.hypot(s_src[1], s_ref)
