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


def test_noise_refuses_weight_window_generators(run_in_tmpdir):
    # Generators tally a real flux; noise mode must refuse them
    with pytest.raises(RuntimeError, match='[Gg]enerators'):
        _ww_model(windows=False, generator=True).run()
