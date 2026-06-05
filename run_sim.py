"""
Define the HPVsim simulation for Nigeria (hpvsim v2.3.0).

Nigeria-specific parameters are inlined here (single-country repo convention).
Uses hpvsim default mixing and initial conditions. Parameters match the
calibration in results/nigeria_pars.obj (ported from hpvsim_pxv_younger).
"""
import os
os.environ.update(
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1',
)

import numpy as np
import sciris as sc
import hpvsim as hpv

LOCATION = 'nigeria'


# %% v2.3 compatibility helpers (layer_probs / cross-layer probs are annual in v2.3+)
def _to_annual_prob(p, dt):
    p = np.clip(p, 0, 1 - 1e-10)
    return 1 - (1 - p) ** (1 / dt)


def _layer_probs_to_annual(layer_probs, dt):
    out = {}
    for lkey, lp in layer_probs.items():
        lp_new = np.asarray(lp).copy().astype(float)
        for row in [1, 2]:
            lp_new[row, :] = _to_annual_prob(lp_new[row, :], dt)
        out[lkey] = lp_new
    return out


def _convert_calib_pars_to_annual(calib_pars, dt):
    if calib_pars is None:
        return calib_pars
    out = dict(calib_pars)
    for key in ('m_cross_layer', 'f_cross_layer'):
        if key in out and out[key] is not None:
            out[key] = float(_to_annual_prob(out[key], dt))
    if 'layer_probs' in out and out['layer_probs'] is not None:
        out['layer_probs'] = _layer_probs_to_annual(out['layer_probs'], dt)
    return out


# Marital ('m') and casual ('c') partnership probabilities by age, fitted to 2018
# Nigeria DHS; ported from hpvsim_pxv_younger. Rows: [age bins], [female], [male].
def _nigeria_layer_probs():
    return dict(
        m=np.array([
            [0, 5, 10, 15,  20,  25,   30,   35,  40,  45,  50,  55,  60,  65,    70,   75],
            [0, 0, 0, 0.1, 0.1, 0.15, 0.15, 0.15, 0.2, 0.3, 0.4, 0.4, 0.2, 0.07, 0.035, 0.007],
            [0, 0, 0, 0.1, 0.1, 0.15, 0.15, 0.2,  0.2, 0.4, 0.4, 0.4, 0.2, 0.1,  0.05,  0.01]]),
        c=np.array([
            [0, 5, 10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,   65,   70,   75],
            [0, 0, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4, 0.7, 0.7, 0.6, 0.2, 0.10, 0.02, 0.02, 0.02],
            [0, 0, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.6, 0.5, 0.2, 0.02, 0.02, 0.02, 0.02]]),
    )


# Sexual debut, fitted to 2018 Nigeria DHS; ported from hpvsim_pxv_younger.
_DEBUT = dict(
    f=dict(dist='lognormal', par1=16., par2=4),
    m=dict(dist='lognormal', par1=18., par2=4),
)
_M_PARTNERS = dict(m=dict(dist='poisson1', par1=0.01), c=dict(dist='poisson1', par1=0.2))
_F_PARTNERS = dict(m=dict(dist='poisson1', par1=0.01), c=dict(dist='poisson1', par1=0.2))


# %% Simulation factory
def make_sim(calib=False, calib_pars=None, debug=0, interventions=None, analyzers=None,
             seed=1, end=None, datafile=None):
    if end is None:
        end = 2100
    if calib:
        end = 2020

    dt = [0.25, 1.0][debug]
    layer_probs = _layer_probs_to_annual(_nigeria_layer_probs(), dt)

    pars = dict(
        n_agents=[20e3, 1e3][debug],
        dt=dt,
        start=[1960, 1980][debug],
        end=end,
        genotypes=[16, 18, 'hi5', 'ohr'],
        location=LOCATION,
        debut=_DEBUT,
        layer_probs=layer_probs,
        f_partners=_F_PARTNERS,
        m_partners=_M_PARTNERS,
        ms_agent_ratio=100,
        verbose=0.0,
    )

    if calib_pars is not None:
        calib_pars = _convert_calib_pars_to_annual(calib_pars, dt)
        pars = sc.mergedicts(pars, calib_pars)

    if analyzers is None:
        analyzers = []
    return hpv.Sim(pars=pars, interventions=interventions, analyzers=analyzers,
                   datafile=datafile, rand_seed=seed)


def run_sim(calib_pars=None, interventions=None, analyzers=None, debug=0, seed=1,
            verbose=0.2, do_save=True, end=2100, do_shrink=True):
    if calib_pars is None:
        calib_pars = sc.loadobj('results/nigeria_pars.obj')
        calib_pars.pop('hiv_pars', None)
    sim = make_sim(debug=debug, end=end, interventions=interventions,
                   analyzers=analyzers, calib_pars=calib_pars)
    sim['rand_seed'] = seed
    sim.label = f'nigeria--{seed}'
    sim['verbose'] = verbose
    sim.run()
    if do_shrink:
        sim.shrink()
    if do_save:
        sim.save('results/nigeria.sim')
    return sim


if __name__ == '__main__':
    T = sc.timer()
    sim = run_sim(end=2020, debug=0, do_save=False)
    sim.plot()
    T.toc('Done')
