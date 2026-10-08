"""Smoke test: scenario arms build and run."""
import os
import pytest
import sciris as sc
import run_sim as rs
import run_scenarios as rsc


def _has_v3_pars(path='results/nigeria_pars.obj'):
    """True if a nested-v3 pars file exists at ``path``."""
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
    """All three scenario builders construct a plausible intervention list."""
    assert rsc.no_interventions() == []
    for make in (rsc.status_quo_interventions, rsc.who_90_70_90_interventions):
        intvs = make()
        labels = [iv.label for iv in intvs]
        assert 'screening' in labels and 'triage' in labels
        assert any('adol_vx' in n for n in [iv.name for iv in intvs])


@pytest.mark.skipif(not HAS_V3_PARS, reason='awaiting v3 recalibration artifacts')
def test_scenarios_build_and_run():
    cp = sc.loadobj('results/nigeria_pars.obj')
    for make in rsc.SCENARIO_BUILDERS.values():
        sim = rs.make_sim(pars=cp, interventions=make(), stop=2100, debug=1)
        sim.run()
        assert 'asr_cancer_incidence' in sim.results['all_hpv']
