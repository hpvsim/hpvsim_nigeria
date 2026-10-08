"""
Calibrate HPVsim Nigeria (hpvsim v3.2).

Heavy calibration (`run_calib`) is fast ONLY on multi-core VMs — never local.
Plotting/extraction (`load_calib`) runs locally.

Calibration targets are provided as long-format CSVs via ``data=``;
``hpv.Calibration`` derives age bins + years from the data, attaches a
``by_age`` analyzer named ``'calib_by_age'``, and extends ``sim.stop`` past
the latest data year as needed.

``calib_pars`` uses the nested v3 form (see hpv.Calibration docstring):
scopes nest by module, leaves are ``[best, low, high]`` lists collapsed
by ``sc.flattendict`` before Optuna sees them.
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


# Set by user before running
to_run = [
    # 'run_calibration',   # uncomment to RUN (VM only)
    'plot_calibration',     # uncomment to PLOT/extract (local)
]
debug = False
do_save = True
n_trials = [2000, 2][debug]  # 2k converges — best came in by trial ~1k on the 5k run
n_workers = [64, 2][debug]

# Top-N trials to keep in the shrunken (committable) calib object.
N_KEEP = 100


def make_calib_pars():
    """Nested [best, low, high] specs for each calibration parameter.

    Priors are centred on hpvsim_pxv_younger's v3.1 Nigeria best-fit pars and
    widened to accommodate hpvsim v3.2's two cancer-affecting changes: higher
    default post-clearance immunity (uniform(0.5, 0.95) replacing Beta mean
    0.35) and smooth age_risk ramp (replacing the step at 30). Both reduce
    cancer burden, so beta/f_cross_layer/rel_sev ceilings go up to let the
    optimizer compensate.
    """
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
        # Full object (Optuna study + all trials) -> raw_results/ (gitignored).
        # Shrunken (top-N trials only) -> results/ (committable).
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
    """Load the shrunken calib from results/ (committed)."""
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
