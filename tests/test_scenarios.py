"""Smoke test: both scenario arms build and run (debug mode, full 2100 horizon)."""
import sciris as sc
import run_sim as rs
import run_scenarios as rsc


def test_scenarios_build_and_run():
    cp = sc.loadobj('results/nigeria_pars.obj')
    cp.pop('hiv_pars', None)
    for make in (rsc.baseline_interventions, rsc.who_interventions):
        sim = rs.make_sim(calib_pars=cp, interventions=make(), end=2100, debug=1)
        sim.run()
        assert 'asr_cancer_incidence' in sim.results
