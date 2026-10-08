"""
Define the HPVsim simulation for Nigeria (hpvsim v3.2).

Nigeria-specific parameters are inlined here (single-country repo convention).
Parameter medians were fit to the 2018 Nigeria DHS; the calibration priors
live in run_calibration.py and the committed best-par set in
results/nigeria_pars.obj (regenerated under v3).
"""
import os
os.environ.update(
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1',
)

import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv


LOCATION = 'nigeria'


# Marital ('m') and casual ('c') partnership probabilities by age, fitted to
# 2018 Nigeria DHS. Rows: [age bins], [female], [male]. Values adopted from
# hpvsim_pxv_younger/model.py (which notes the casual probs are bumped toward
# Kaz-style near-saturation for young ages to give HPV-transmission headroom
# during calibration).
def _nigeria_layer_probs():
    return dict(
        marital=np.array([
            [0, 5, 10, 15,  20,  25,  30,  35,   40,  45,  50,  55,  60,  65,    70,    75],
            [0, 0, 0,  0.1, 0.5, 0.5, 0.4, 0.15, 0.2, 0.3, 0.4, 0.4, 0.2, 0.07, 0.035, 0.007],
            [0, 0, 0,  0.1, 0.5, 0.5, 0.4, 0.2,  0.2, 0.4, 0.4, 0.4, 0.2, 0.1,  0.05,  0.01]]),
        casual=np.array([
            [0, 5, 10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,   70,   75],
            [0, 0, 0.1, 0.8, 0.8, 0.6, 0.5, 0.4, 0.4, 0.3, 0.3, 0.2, 0.1, 0.02, 0.02, 0.02],
            [0, 0, 0.0, 0.5, 0.6, 0.6, 0.7, 0.6, 0.5, 0.5, 0.4, 0.3, 0.1, 0.02, 0.02, 0.02]]),
    )


def network_pars():
    """Nigeria-specific network pars to layer over ``hpv.NetworkPars`` defaults.

    Adopted from hpvsim_pxv_younger/model.py (DHS-fit paper setup), with the
    hpv.NetworkPars default partners_marital kept (not overridden here).
    """
    lp = _nigeria_layer_probs()
    pars = dict(
        debut_f=ss.lognorm_ex(loc=16.0, scale=2.0),
        debut_m=ss.lognorm_ex(loc=19.0, scale=3.0),
        m_partners_casual=0.2,
        f_partners_casual=0.2,
        layer_probs_marital=lp['marital'],
        layer_probs_casual=lp['casual'],
    )
    return pars


def make_sim(pars=None, debug=0, n_agents=None, dt=None, start=None, stop=2100,
             genotypes=None, ms_agent_ratio=100,
             interventions=None, analyzers=None, seed=1, calib=False):
    """Build the Nigeria sim."""
    if n_agents is None:
        n_agents = [20_000, 1_000][debug]
    if dt is None:
        dt = [0.25, 1.0][debug]
    if start is None:
        start = [1960, 1980][debug]
    if genotypes is None:
        genotypes = [16, 18, 'hi5', 'ohr']
    if calib:
        stop = 2020
    pars = sc.mergedicts(network_pars(), pars)

    sim = hpv.Sim(
        location=LOCATION,
        n_agents=n_agents,
        dt=dt,
        start=start,
        stop=stop,
        genotypes=genotypes,
        ms_agent_ratio=ms_agent_ratio,
        interventions=interventions,
        analyzers=analyzers,
        rand_seed=seed,
        pars=pars,
    )
    return sim


def run_sim(pars=None, seed=1, do_save=False, do_shrink=True, **kwargs):
    sim = make_sim(pars=pars, seed=seed, **kwargs)
    sim.label = f'nigeria--{seed}'
    sim.run()
    if do_shrink:
        sim.shrink()
    if do_save:
        sim.save('results/nigeria.sim')
    return sim


if __name__ == '__main__':
    T = sc.timer()
    sim = run_sim(stop=2020)
    sim.plot()
    T.toc('Done')
