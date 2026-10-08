"""Nigeria HPV scenarios: no-interventions, status-quo, WHO 90-70-90 by 2030.

Screening uses a 10-year rescreen gap over ages 30-50. Vaccine product is
quadrivalent (Nigeria's Gavi-supplied programme).
"""
import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv

import run_sim as rs


# -- Timeline ----------------------------------------------------------------
HIST = 2020
SQ_END = 2025   # last year of historical adol rollout
WHO_END = 2030  # year by which WHO 90-70-90 targets are reached
END_SIM = 2100
YEARS = np.arange(HIST, END_SIM + 1)

# -- Vaccination -------------------------------------------------------------
VX_PRODUCT = 'quadrivalent'
VX_SP = 0.95
VX_AGE_RANGE = (9, 10)

# Nigeria historical adol vax rollout 2023-2025 (aggregate coverage).
BASE_ADOL_RAMP = {2023: 0.27, 2024: 0.60, 2025: 0.60}
SQ_PLATEAU = 0.60
WHO_VX_TARGET = 0.90

# -- Screening ---------------------------------------------------------------
SCREEN_AGE_LO, SCREEN_AGE_HI = 30, 50
SCREEN_WINDOW_YEARS = SCREEN_AGE_HI - SCREEN_AGE_LO
SCREEN_GAP_YEARS = 10
SQ_SCREEN_LIFETIME = 0.15
WHO_SCREEN_LIFETIME = 0.70

# -- Treatment ---------------------------------------------------------------
SQ_TREAT_COV = 0.50
WHO_TREAT_COV = 0.90


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _annual_from_lifetime(lifetime_cov, n_years=SCREEN_WINDOW_YEARS):
    """Per-year prob equivalent to lifetime coverage ``C`` over ``N`` years:
    inverts ``1 - (1 - p)^N == C``."""
    c = float(min(max(lifetime_cov, 0.0), 1.0))
    if c <= 0: return 0.0
    if c >= 1.0: return 1.0
    return 1.0 - (1.0 - c) ** (1.0 / n_years)


def _linear_ramp(v0, v1, y0, y1, years=YEARS):
    """Flat ``v0`` through ``y0``, linear to ``v1`` by ``y1``, flat thereafter."""
    return np.interp(years, [years[0], y0, y1, years[-1]], [v0, v0, v1, v1])


def _not_recently_screened(sim, intv_name='screening', gap_years=SCREEN_GAP_YEARS):
    """Eligibility callback: alive UIDs not screened in the last ``gap_years``.

    v3 stores screening history on ``intv.screened`` + ``intv.ti_screened``.
    Must return ``ss.uids`` (v3.2 guardrail).
    """
    alive = sim.people.alive.uids
    intv = sim.interventions.get(intv_name, None)
    if intv is None:
        return alive
    screened_uids = ss.uids(intv.screened.uids)
    if not len(screened_uids):
        return alive
    ti_thresh = sim.ti - int(gap_years / sim.t.dt_year)
    ti_arr = np.asarray(intv.ti_screened[screened_uids])
    recent = np.asarray(screened_uids)[ti_arr > ti_thresh]
    return alive.remove(ss.uids(recent))


def _make_tx_assigner():
    """Override shipped ``tx_assigner`` routing.

    Shipped product sends only 5% of pre-cancer screen-positives to ablation
    and 0% to excision (the rest default to 'none' / untreated), which swamps
    any screening effect. Override: 80% ablation / 20% excision for pre-cancer,
    100% radiation for cancer. Susceptible left alone (false positive).
    """
    prod = hpv.products.dx(name='tx_assigner')
    df = prod.df

    def _set(state, ablation=0.0, excision=0.0, radiation=0.0):
        mask = df['state'] == state
        df.loc[mask & (df['result'] == 'ablation'),  'probability'] = ablation
        df.loc[mask & (df['result'] == 'excision'),  'probability'] = excision
        df.loc[mask & (df['result'] == 'radiation'), 'probability'] = radiation
        df.loc[mask & (df['result'] == 'none'),      'probability'] = \
            1.0 - ablation - excision - radiation

    for pre in ('latent', 'precin', 'cin'):
        _set(pre, ablation=0.80, excision=0.20)
    _set('cancerous', radiation=1.00)
    return prod


# ---------------------------------------------------------------------------
# Intervention builders
# ---------------------------------------------------------------------------

def _historical_adol_vx(name_prefix='hist'):
    """Nigeria 2023-2025 adol routine vx rollout."""
    years = list(BASE_ADOL_RAMP.keys())
    probs = list(BASE_ADOL_RAMP.values())
    prod = hpv.vx(name=VX_PRODUCT, sterilizing_p=VX_SP)
    # Product name must differ from intervention name (v3.2 guardrail).
    prod.name = f'vx_{name_prefix}_hist'
    return hpv.routine_vx(
        name=f'{name_prefix}_adol_vx_hist',
        product=prod,
        age_range=VX_AGE_RANGE,
        years=years, prob=probs,
    )


def _post_2025_adol_vx(name_prefix, start_cov, end_cov, ramp_end_year):
    """Post-2025 adol vx: scalar if start==end, else linear ramp then flat."""
    years = list(range(SQ_END + 1, END_SIM + 1))
    if abs(end_cov - start_cov) < 1e-9:
        probs = [start_cov] * len(years)
    else:
        n_ramp = ramp_end_year - SQ_END
        probs = []
        for y in years:
            if y <= ramp_end_year:
                probs.append(start_cov + (end_cov - start_cov) * (y - SQ_END) / n_ramp)
            else:
                probs.append(end_cov)
    prod = hpv.vx(name=VX_PRODUCT, sterilizing_p=VX_SP)
    prod.name = f'vx_{name_prefix}_post'
    return hpv.routine_vx(
        name=f'{name_prefix}_post2025_vx',
        product=prod,
        age_range=VX_AGE_RANGE,
        years=years, prob=probs,
    )


def _screen_treat_cascade(baseline_lifetime, target_lifetime, ramp_end_year,
                          baseline_treat, target_treat):
    """Screen-triage-treat cascade with linear ramps 2025 -> ramp_end_year.

    Triage prob is the treatment-cascade coverage dial (captures LTFU).
    Treatment products routed via the overridden tx_assigner.
    """
    scov = _linear_ramp(baseline_lifetime, target_lifetime, SQ_END, ramp_end_year)
    scov_per_year = np.array([_annual_from_lifetime(c) for c in scov])

    tcov = _linear_ramp(baseline_treat, target_treat, SQ_END, ramp_end_year)

    screening = hpv.routine_screening(
        name='screening', label='screening',
        prob=scov_per_year, years=YEARS,
        product='hpv',
        age_range=[SCREEN_AGE_LO, SCREEN_AGE_HI],
        eligibility=_not_recently_screened,
    )

    triage = hpv.routine_triage(
        name='triage', label='triage',
        prob=tcov, years=YEARS, annual_prob=False,
        product=_make_tx_assigner(),
        eligibility=lambda sim: ss.uids(sim.interventions['screening'].outcomes['positive']),
    )

    # Intervention name must differ from product module name (v3.2 guardrail).
    abl_prod = hpv.products.tx(name='ablation'); abl_prod.name = 'prod_ablation'
    exc_prod = hpv.products.tx(name='excision'); exc_prod.name = 'prod_excision'
    rad_prod = hpv.radiation();                  rad_prod.name = 'prod_radiation'

    ablation = hpv.treat_num(
        name='ablation_program', label='ablation', product=abl_prod, prob=1.0,
        eligibility=lambda sim: ss.uids(sim.interventions['triage'].outcomes['ablation']),
    )
    excision = hpv.treat_num(
        name='excision_program', label='excision', product=exc_prod, prob=1.0,
        eligibility=lambda sim: ss.uids(np.unique(np.concatenate([
            np.asarray(sim.interventions['triage'].outcomes['excision']),
            np.asarray(sim.interventions['ablation_program'].outcomes['unsuccessful']),
        ])).astype(int)),
    )
    radiation = hpv.treat_num(
        name='radiation_program', label='radiation', product=rad_prod, prob=1.0,
        eligibility=lambda sim: ss.uids(sim.interventions['triage'].outcomes['radiation']),
    )
    return [screening, triage, ablation, excision, radiation]


# ---------------------------------------------------------------------------
# Scenario definitions
# ---------------------------------------------------------------------------

def no_interventions():
    """Pure natural history."""
    return []


def status_quo_interventions():
    """Historical adol rollout continuing at 60%, 15% lifetime screening, 50% treat cov."""
    vx = [_historical_adol_vx('sq'),
          _post_2025_adol_vx('sq', SQ_PLATEAU, SQ_PLATEAU, SQ_END)]
    st = _screen_treat_cascade(
        baseline_lifetime=SQ_SCREEN_LIFETIME,
        target_lifetime=SQ_SCREEN_LIFETIME,
        ramp_end_year=SQ_END,
        baseline_treat=SQ_TREAT_COV,
        target_treat=SQ_TREAT_COV,
    )
    return vx + st


def who_90_70_90_interventions():
    """Linear scale-up 2025 -> 2030 to 90/70/90, held through 2100."""
    vx = [_historical_adol_vx('who'),
          _post_2025_adol_vx('who', SQ_PLATEAU, WHO_VX_TARGET, WHO_END)]
    st = _screen_treat_cascade(
        baseline_lifetime=SQ_SCREEN_LIFETIME,
        target_lifetime=WHO_SCREEN_LIFETIME,
        ramp_end_year=WHO_END,
        baseline_treat=SQ_TREAT_COV,
        target_treat=WHO_TREAT_COV,
    )
    return vx + st


SCENARIO_BUILDERS = {
    'No interventions': no_interventions,
    'Status quo': status_quo_interventions,
    'WHO 90-70-90 by 2030': who_90_70_90_interventions,
}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def run_scenarios(calib_pars=None, end=END_SIM, n_seeds=3, do_save=True):
    """Run all scenarios across ``n_seeds`` seeds; interventions rebuilt per seed."""
    if calib_pars is None:
        calib_pars = sc.loadobj('results/nigeria_pars.obj')
    msims = sc.objdict()
    for name, build in SCENARIO_BUILDERS.items():
        sims = [rs.make_sim(pars=calib_pars, interventions=build(),
                             stop=end, seed=s)
                for s in range(n_seeds)]
        msim = ss.MultiSim(sims=sims)
        msim.run()
        msims[name] = msim
    if do_save:
        sc.saveobj('results/nigeria_scenarios.obj', msims)
    return msims


# Aliases consumed by tests.
baseline_interventions = status_quo_interventions
who_interventions = who_90_70_90_interventions


if __name__ == '__main__':
    T = sc.timer()
    msims = run_scenarios(end=END_SIM, n_seeds=3)
    T.toc('Done')
