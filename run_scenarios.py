"""
Nigeria Baseline vs WHO scale-up scenarios (hpvsim v2.3.0).

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
import hpvsim as hpv
import run_sim as rs

HIST = 2020
RAMP0, RAMP1 = 2025, 2030
END_SIM = 2100
YEARS = np.arange(HIST, END_SIM + 1)


def _ramp(v0, v1, y0=RAMP0, y1=RAMP1):
    """Coverage flat at v0 to y0, linear to v1 by y1, flat after (aligned to YEARS)."""
    return np.interp(YEARS, [HIST, y0, y1, END_SIM], [v0, v0, v1, v1])


def _annual_screen_prob(cov):
    """Convert target screening coverage of women 30-50 to a per-year screening
    probability (hpvsim convention from the source repos: denominator = half the
    20-year age span). The 5-year re-screen interval is enforced separately by the
    screening eligibility function."""
    return 1 - (1 - np.asarray(cov)) ** (1 / ((50 - 30) / 2))


def make_vaccination(cov0, cov1, start=2025):
    """Routine HPV vaccination of girls 9-10, beginning `start` (0 before)."""
    cov = np.interp(YEARS, [HIST, start - 1, start, RAMP1, END_SIM], [0, 0, cov0, cov1, cov1])
    return hpv.routine_vx(prob=cov, years=YEARS, product='nonavalent',
                          age_range=(9, 10), label='routine_vx')


def make_st(screen_cov0, screen_cov1, product, linkage0, linkage1, linkage_fast,
            treat_coverage, cancer_tx_coverage):
    """Screen -> triage -> ablation/excision/radiation cascade with ramped coverage."""
    scov = _ramp(screen_cov0, screen_cov1)
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | \
        (sim.t > (sim.people.date_screened + 5 / sim['dt']))
    screening = hpv.routine_screening(
        prob=_annual_screen_prob(scov), years=YEARS, eligibility=screen_eligible,
        product=product, age_range=[30, 50], label='screening')

    pos = lambda sim: sim.get_intervention('screening').outcomes['positive']
    lcov = _ramp(linkage0, linkage1, y0=2024, y1=2025) if linkage_fast else _ramp(linkage0, linkage1)
    triage = hpv.routine_triage(prob=lcov, years=YEARS, annual_prob=False,
                                product='tx_assigner', eligibility=pos, label='triage')

    abl_elig = lambda sim: sim.get_intervention('triage').outcomes['ablation']
    ablation = hpv.treat_num(prob=treat_coverage, product='ablation',
                             eligibility=abl_elig, label='ablation')
    exc_elig = lambda sim: list(set(sim.get_intervention('triage').outcomes['excision'].tolist() +
                                    sim.get_intervention('ablation').outcomes['unsuccessful'].tolist()))
    excision = hpv.treat_num(prob=treat_coverage, product='excision',
                             eligibility=exc_elig, label='excision')
    rad_elig = lambda sim: sim.get_intervention('triage').outcomes['radiation']
    radiation = hpv.treat_num(prob=cancer_tx_coverage, product=hpv.radiation(),
                              eligibility=rad_elig, label='radiation')

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


def run_scenarios(end=END_SIM, n_seeds=3, do_save=True):
    """Run each scenario across `n_seeds` seeds; return {name: MultiSim}.

    Interventions are rebuilt per seed (they bind to a sim on init).
    """
    calib_pars = sc.loadobj('results/nigeria_pars.obj')
    calib_pars.pop('hiv_pars', None)
    builders = {'Baseline': baseline_interventions, 'WHO': who_interventions}
    msims = sc.objdict()
    for name, build in builders.items():
        sims = [rs.make_sim(calib_pars=calib_pars, interventions=build(), end=end, seed=s)
                for s in range(n_seeds)]
        msim = hpv.MultiSim(sims)
        msim.run()
        msims[name] = msim
    if do_save:
        sc.saveobj('results/nigeria_scenarios.obj', msims)
    return msims


if __name__ == '__main__':
    T = sc.timer()
    msims = run_scenarios(end=END_SIM, n_seeds=3)
    T.toc('Done')
