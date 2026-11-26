"""Tests of the STD tolerance rules in harness.std_criterion (amendment
2026-11-26b: the tolerance must be exceeded at k sigma for an entry to fail).

    python -m pytest noise_dev/common/test_harness.py

"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from harness import std_criterion  # noqa: E402


def _case(rel_amp, sigma_rel, phase_deg=0.0, sigma_phase_deg=0.1):
    """One entry with reference 1 and a given relative amplitude error,
    relative amplitude sigma and phase (sigma split so that the amplitude
    and phase sigmas come out as given, for a mean on the real axis)."""
    mean = np.array([(1 + rel_amp) * np.exp(1j * np.radians(phase_deg))])
    A = abs(mean[0])
    s_re = np.array([sigma_rel])                     # along the amplitude
    s_im = np.array([np.radians(sigma_phase_deg) * A])
    return mean, s_re, s_im, np.array([1.0 + 0j])


def test_default_rule_unchanged():
    # 5.5 % off with sigma 1.5 %: evaluated (1.5 % < 5/3 %), fails today
    m, sr, si, ref = _case(-0.055, 0.015)
    ok, info = std_criterion(m, sr, si, ref, 0.05, 2.0)
    assert info["amp_eval"] == 1
    assert not info["ok_amp"]
    assert info["n_amp_fail"] == 1


def test_three_sigma_rule_tolerates_fluctuation():
    # the 9a misses: 5.3-5.7 % off at sigma 1.4-1.5 %
    m, sr, si, ref = _case(-0.057, 0.015)
    _, info = std_criterion(m, sr, si, ref, 0.05, 2.0, significance=3)
    assert info["ok_amp"] and info["n_amp_fail"] == 0


def test_three_sigma_rule_catches_real_excess():
    # 6 % off at sigma 0.2 %: beyond 5 % + 3 sigma = 5.6 %
    m, sr, si, ref = _case(0.06, 0.002)
    _, info = std_criterion(m, sr, si, ref, 0.05, 2.0, significance=3)
    assert not info["ok_amp"] and info["n_amp_fail"] == 1
    # phase: 2.5 deg off at sigma 0.1 deg fails; at sigma 0.2 deg passes
    m, sr, si, ref = _case(0.0, 0.002, phase_deg=2.5, sigma_phase_deg=0.1)
    _, info = std_criterion(m, sr, si, ref, 0.05, 2.0, significance=3)
    assert not info["ok_ph"]
    m, sr, si, ref = _case(0.0, 0.002, phase_deg=2.5, sigma_phase_deg=0.2)
    _, info = std_criterion(m, sr, si, ref, 0.05, 2.0, significance=3)
    assert info["ok_ph"]


def test_evaluation_threshold_unchanged():
    # sigma 2 % > 5/3 %: not evaluated under either rule
    m, sr, si, ref = _case(-0.2, 0.02)
    for k in (None, 3):
        _, info = std_criterion(m, sr, si, ref, 0.05, 2.0, significance=k)
        assert info["amp_eval"] == 0 and info["ok_amp"]
