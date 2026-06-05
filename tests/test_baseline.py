"""Smoke + validation tests for the Nigeria baseline."""
import numpy as np
import sciris as sc
import run_sim as rs

# data/nigeria_asr_cancer_incidence.csv. Note: this is the repo's calibration
# target; a separate Globocan figure (~26.2) is higher — the model is calibrated
# to case counts + genotype distributions, which yield ASR ~18.
TARGET_ASR_2020 = 18.4


def test_pars_load():
    pars = sc.loadobj('results/nigeria_pars.obj')
    assert isinstance(pars, dict) and 'beta' in pars


def test_make_sim_debug_runs():
    """Fast smoke test (debug=1) suitable for CI."""
    sim = rs.make_sim(debug=1, calib_pars=None, end=2000)
    sim.run()
    assert sim.results['year'][-1] >= 2000


def test_baseline_asr_in_range():
    """Full-resolution baseline ASR near the calibration target (18.4/100k)."""
    pars = sc.loadobj('results/nigeria_pars.obj')
    pars.pop('hiv_pars', None)
    sim = rs.make_sim(calib_pars=pars, end=2020, debug=0)
    sim.run()
    asr = float(np.array(sim.results['asr_cancer_incidence'])[-1])
    assert 0.75 * TARGET_ASR_2020 <= asr <= 1.25 * TARGET_ASR_2020, \
        f'ASR={asr:.1f} vs target {TARGET_ASR_2020}'
