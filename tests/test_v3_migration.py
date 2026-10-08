"""Smoke test: the v3.2 port runs end-to-end without a calibration artifact.

Runs on CI with no v3 pars file; the baseline/scenarios tests cover the
richer paths once ``results/nigeria_pars.obj`` is regenerated under v3.
"""
import hpvsim as hpv
import run_sim as rs
import run_scenarios as rsc


def test_hpvsim_version():
    assert hpv.__version__.startswith('3.2'), f'need v3.2, got {hpv.__version__}'


def test_baseline_sim_runs_debug():
    sim = rs.make_sim(debug=1, stop=1990)
    sim.run()
    assert 'all_hpv' in sim.results
    assert 'asr_cancer_incidence' in sim.results['all_hpv']


def test_scenario_sim_runs_debug():
    """A WHO-arm sim with no calibration overrides should still build + run.

    Scenario interventions are defined over YEARS=[2020, 2100]; starsim
    requires those years to fit within sim [start, stop], so the sim has
    to run to 2100.
    """
    sim = rs.make_sim(interventions=rsc.who_interventions(), debug=1, stop=2100)
    sim.run()
    assert sim.timevec[-1] >= 2100
