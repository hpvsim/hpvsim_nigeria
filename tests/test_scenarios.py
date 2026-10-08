"""Smoke test: both scenario arms build and run (debug mode, full 2100 horizon)."""
import os
import pytest
import sciris as sc
import run_sim as rs
import run_scenarios as rsc


def _has_v3_pars(path='results/nigeria_pars.obj'):
    if not os.path.exists(path):
        return False
    try:
        pars = sc.loadobj(path)
    except Exception:
        return False
    v2_keys = {'genotype_pars', 'sev_dist', 'hiv_pars', 'f_partners', 'm_partners'}
    return isinstance(pars, dict) and not (set(pars) & v2_keys)


HAS_V3_PARS = _has_v3_pars()


def test_scenario_builders():
    """Intervention lists construct without a sim or calib pars."""
    for make in (rsc.baseline_interventions, rsc.who_interventions):
        intvs = make()
        labels = [iv.label for iv in intvs]
        assert 'routine_vx' in labels and 'screening' in labels and 'triage' in labels


@pytest.mark.skipif(not HAS_V3_PARS, reason='awaiting v3 recalibration artifacts')
def test_scenarios_build_and_run():
    cp = sc.loadobj('results/nigeria_pars.obj')
    for make in (rsc.baseline_interventions, rsc.who_interventions):
        sim = rs.make_sim(pars=cp, interventions=make(), stop=2100, debug=1)
        sim.run()
        assert 'asr_cancer_incidence' in sim.results['all_hpv']
