"""Calibrate HPVsim Nigeria.

Calibration targets are passed as long-format CSVs via ``data=``;
``hpv.Calibration`` derives age bins and years from the data and attaches
a ``by_age`` analyzer. ``calib_pars`` uses the nested v3 form: scopes nest
by module, leaves are ``[best, low, high]``.
"""
import os
os.environ.update(
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1',
)

import sciris as sc
import hpvsim as hpv

import run_sim as rs
import utils as ut


to_run = [
    # 'run_calibration',   # VM only
    'plot_calibration',
]
debug = False
do_save = True
n_trials = [2000, 2][debug]
n_workers = [64, 2][debug]

# Top-N trials kept in the shrunken (committable) calib object.
N_KEEP = 100


def make_calib_pars():
    """Nested [best, low, high] specs for each calibration parameter."""
    pars = dict(
        beta=[0.34, 0.15, 0.50],
        m_cross_layer=[0.60, 0.10, 0.90],
        f_cross_layer=[0.70, 0.20, 0.95],
        network=dict(
            m_partners_casual=[0.37, 0.10, 0.70],
            f_partners_casual=[0.58, 0.10, 0.90],
        ),
        cross_immunity=dict(rel_sev=dict(loc=[0.80, 0.30, 2.00])),
    )
    for g in ['hi5', 'ohr']:
        pars[g] = dict(
            cancer_fn=dict(transform_prob=[1.0e-3, 0.3e-3, 3.0e-3]),
            cin_fn=dict(k=[0.20, 0.10, 0.35]),
            dur_cin=dict(mean=[4.5, 3.0, 7.0], std=[20, 10, 30]),
        )
    return pars


def run_calib(n_trials=None, n_workers=None, do_plot=False, do_save=True, filestem=''):
    calib = hpv.Calibration(
        rs.make_sim(calib=True), calib_pars=make_calib_pars(),
        data=ut.make_datafiles(), total_trials=n_trials, n_workers=n_workers,
    )
    try:
        calib.calibrate()
    except Exception as e:
        print(f'calibrate() raised: {e}; saving partial results anyway')
    if do_save:
        # Full object -> raw_results/ (gitignored); shrunken -> results/ (committed).
        sc.saveobj(f'raw_results/nigeria_calib{filestem}.obj', calib)
        shrunk = calib.shrink(n_results=N_KEEP)
        sc.saveobj(f'results/nigeria_calib{filestem}.obj', shrunk)
        sc.saveobj(f'results/nigeria_pars{filestem}.obj', calib.best_pars)
    if do_plot:
        fig = hpv.plot_calibration(calib)
        fig.savefig(f'figures/nigeria_calib{filestem}.png')
    if getattr(calib, 'best_pars', None) is not None:
        print(f'Best pars: {calib.best_pars}')
    return calib


def load_calib(do_plot=True, filestem=''):
    calib = sc.load(f'results/nigeria_calib{filestem}.obj')
    if do_plot:
        fig = hpv.plot_calibration(calib)
        fig.savefig(f'figures/nigeria_calib{filestem}.png')
    return calib


if __name__ == '__main__':
    T = sc.timer()
    if 'run_calibration' in to_run:
        run_calib(n_trials=n_trials, n_workers=n_workers, do_save=do_save)
    if 'plot_calibration' in to_run:
        load_calib(do_plot=True)
    T.toc('Done')
