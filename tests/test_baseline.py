"""Smoke + validation tests for the Nigeria baseline."""
import os
import numpy as np
import pytest
import sciris as sc
import run_sim as rs

# data/nigeria_asr_cancer_incidence.csv. Note: this is the repo's calibration
# target; a separate Globocan figure (~26.2) is higher — the model is calibrated
# to case counts + genotype distributions, which yield ASR ~18.
TARGET_ASR_2020 = 18.4


def _has_v3_pars(path='results/nigeria_pars.obj'):
    """v2 pars files have top-level 'genotype_pars' / 'sev_dist' / 'hiv_pars';
    v3 uses nested 'network' / 'cross_immunity' / genotype keys."""
    if not os.path.exists(path):
        return False
    try:
        pars = sc.loadobj(path)
    except Exception:
        return False
    v2_keys = {'genotype_pars', 'sev_dist', 'hiv_pars', 'f_partners', 'm_partners'}
    return isinstance(pars, dict) and not (set(pars) & v2_keys)


HAS_V3_PARS = _has_v3_pars()


def test_make_sim_debug_runs():
    """Fast smoke test (debug=1) suitable for CI."""
    sim = rs.make_sim(debug=1, stop=2000)
    sim.run()
    assert sim.timevec[-1] >= 2000


@pytest.mark.skipif(not HAS_V3_PARS, reason='awaiting v3 recalibration artifacts')
def test_pars_load():
    pars = sc.loadobj('results/nigeria_pars.obj')
    assert isinstance(pars, dict) and len(pars) > 0


@pytest.mark.skipif(not HAS_V3_PARS, reason='awaiting v3 recalibration artifacts')
def test_baseline_asr_in_range():
    """Full-resolution baseline ASR near the calibration target (18.4/100k)."""
    pars = sc.loadobj('results/nigeria_pars.obj')
    sim = rs.make_sim(pars=pars, stop=2020, debug=0)
    sim.run()
    asr = float(np.asarray(sim.results['all_hpv']['asr_cancer_incidence'])[-1])
    assert 0.75 * TARGET_ASR_2020 <= asr <= 1.25 * TARGET_ASR_2020, \
        f'ASR={asr:.1f} vs target {TARGET_ASR_2020}'
