"""
Nigeria Baseline vs WHO 90-70-90 scale-up scenarios (hpvsim v3.2).

Encodes the partner's parameter table:
  - Routine HPV vaccination of girls 9-10 (begins 2025).
  - Screening of women 30-50 (VIA- or HPV-test-based).
  - Triage (colposcopy linkage) of screen-positives.
  - Ablation / excision of precancer; radiation as cancer treatment.

"Gradual increase ... by 2030" rows ramp linearly 2025->2030 then hold flat;
colposcopy linkage uses a fast (2024->2025) ramp to approximate a step change.
Treatment coverages (ablation/excision/cancer-tx) are held at their scenario
level (treat_num takes a scalar probability). These are designed for a 2100
horizon; adjust coverages/timings to the question at hand.
"""
import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv
import run_sim as rs

HIST = 2020
RAMP0, RAMP1 = 2025, 2030
END_SIM = 2100
YEARS = np.arange(HIST, END_SIM + 1)
SCREEN_GAP_YEARS = 5


def _not_recently_screened(sim, intv_name='screening', gap_years=SCREEN_GAP_YEARS):
    """UIDs of alive people NOT screened in the last ``gap_years`` years.

    In v3 the screening intervention stores ``ti_screened`` + ``screened``
    (a BoolState); eligibility must subtract recently-screened UIDs from
    the alive population. Age filter (30-50) is applied by the screening
    intervention itself via ``age_range``.
    """
    alive = sim.people.alive.uids
    intv = sim.interventions.get(intv_name, None)
    if intv is None:
        return alive
    screened_uids = ss.uids(intv.screened.uids)
    if not len(screened_uids):
        return alive
    dt_year = sim.t.dt_year
    ti_thresh = sim.ti - int(gap_years / dt_year)
    ti_arr = np.asarray(intv.ti_screened[screened_uids])
    recent = np.asarray(screened_uids)[ti_arr > ti_thresh]
    return alive.remove(ss.uids(recent))


def _ramp(v0, v1, y0=RAMP0, y1=RAMP1):
    """Coverage flat at v0 to y0, linear to v1 by y1, flat after (aligned to YEARS)."""
    return np.interp(YEARS, [HIST, y0, y1, END_SIM], [v0, v0, v1, v1])


def _annual_screen_prob(cov):
    """Convert target screening coverage of women 30-50 to a per-year screening
    probability (denominator = half the 20-year age span). The 5-year re-screen
    interval is enforced separately by the screening eligibility function."""
    return 1 - (1 - np.asarray(cov)) ** (1 / ((50 - 30) / 2))


def make_vaccination(cov0, cov1, start=2025):
    """Routine HPV vaccination of girls 9-10, beginning `start` (0 before)."""
    cov = np.interp(YEARS, [HIST, start - 1, start, RAMP1, END_SIM], [0, 0, cov0, cov1, cov1])
    return hpv.routine_vx(prob=cov, years=YEARS, product='nonavalent',
                          age_range=(9, 10), name='routine_vx', label='routine_vx')


def make_st(screen_cov0, screen_cov1, product, linkage0, linkage1, linkage_fast,
            treat_coverage, cancer_tx_coverage):
    """Screen -> triage -> ablation/excision/radiation cascade with ramped coverage.

    In v3 ``sim.interventions`` is an ndict keyed on ``name=``, so each
    intervention needs a unique name even when the class is the same
    (three ``treat_num`` instances here).
    """
    scov = _ramp(screen_cov0, screen_cov1)
    screening = hpv.routine_screening(
        prob=_annual_screen_prob(scov), years=YEARS, eligibility=_not_recently_screened,
        product=product, age_range=[30, 50], name='screening', label='screening')

    pos = lambda sim: ss.uids(sim.interventions['screening'].outcomes['positive'])
    lcov = _ramp(linkage0, linkage1, y0=2024, y1=2025) if linkage_fast else _ramp(linkage0, linkage1)
    triage = hpv.routine_triage(prob=lcov, years=YEARS, annual_prob=False,
                                product='tx_assigner', eligibility=pos,
                                name='triage', label='triage')

    # Note: intervention ``name=`` must differ from its product's module name
    # (v3 guardrail; see interventions._check_name_collision). Use *_program.
    # Eligibility callbacks must return ``ss.uids`` (not raw np arrays or lists).
    abl_elig = lambda sim: ss.uids(sim.interventions['triage'].outcomes['ablation'])
    ablation = hpv.treat_num(prob=treat_coverage, product='ablation',
                             eligibility=abl_elig, name='ablation_program', label='ablation')
    exc_elig = lambda sim: ss.uids(np.unique(np.concatenate([
        np.asarray(sim.interventions['triage'].outcomes['excision']),
        np.asarray(sim.interventions['ablation_program'].outcomes['unsuccessful']),
    ])).astype(int))
    excision = hpv.treat_num(prob=treat_coverage, product='excision',
                             eligibility=exc_elig, name='excision_program', label='excision')
    rad_elig = lambda sim: ss.uids(sim.interventions['triage'].outcomes['radiation'])
    radiation = hpv.treat_num(prob=cancer_tx_coverage, product=hpv.radiation(),
                              eligibility=rad_elig, name='radiation_program', label='radiation')

    return [screening, triage, ablation, excision, radiation]


def baseline_interventions():
    return [make_vaccination(0.84, 0.84)] + make_st(
        screen_cov0=0.09, screen_cov1=0.09, product='via',
        linkage0=0.20, linkage1=0.20, linkage_fast=False,
        treat_coverage=0.50, cancer_tx_coverage=0.44)


def who_interventions():
    return [make_vaccination(0.84, 0.90)] + make_st(
        screen_cov0=0.09, screen_cov1=0.70, product='hpv',
        linkage0=0.20, linkage1=0.70, linkage_fast=True,
        treat_coverage=0.90, cancer_tx_coverage=0.90)


def run_scenarios(calib_pars=None, end=END_SIM, n_seeds=3, do_save=True):
    """Run each scenario across `n_seeds` seeds; return {name: MultiSim}.

    Interventions are rebuilt per seed (they bind to a sim on init).
    If ``calib_pars`` is None, loads from ``results/nigeria_pars.obj``.
    """
    if calib_pars is None:
        calib_pars = sc.loadobj('results/nigeria_pars.obj')
    builders = {'Baseline': baseline_interventions, 'WHO': who_interventions}
    msims = sc.objdict()
    for name, build in builders.items():
        sims = [rs.make_sim(pars=calib_pars, interventions=build(), stop=end, seed=s)
                for s in range(n_seeds)]
        msim = ss.MultiSim(sims=sims)
        msim.run()
        msims[name] = msim
    if do_save:
        sc.saveobj('results/nigeria_scenarios.obj', msims)
    return msims


if __name__ == '__main__':
    T = sc.timer()
    msims = run_scenarios(end=END_SIM, n_seeds=3)
    T.toc('Done')
