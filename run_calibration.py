"""
Calibrate HPVsim Nigeria.

Heavy calibration (`run_calibration`) is fast ONLY on multi-core VMs — never local.
Plotting/extraction (`plot_calibration`) runs locally.
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
n_trials = [4000, 2][debug]
n_workers = [40, 1][debug]
storage = None


def make_priors():
    default = dict(
        cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3, 2e-4]),
        cin_fn=dict(k=[.15, .1, .25, 0.01]),
        dur_cin=dict(par1=[4.5, 3.5, 5.5, 0.5], par2=[20, 16, 24, 0.5]),
    )
    return dict(hi5=sc.dcp(default), ohr=sc.dcp(default))


def run_calib(n_trials=None, n_workers=None, do_plot=False, do_save=True, filestem=''):
    sim = rs.make_sim(calib=True)
    datafiles = ut.make_datafiles()
    calib_pars = dict(
        beta=[0.2, 0.1, 0.34, 0.02],
        m_cross_layer=[0.3, 0.1, 0.7, 0.05],
        m_partners=dict(c=dict(par1=[0.2, 0.1, 0.6, 0.02])),
        f_cross_layer=[0.1, 0.05, 0.5, 0.05],
        f_partners=dict(c=dict(par1=[0.2, 0.1, 0.6, 0.02])),
        sev_dist=dict(par1=[1, 0.5, 1.5, 0.01]),
    )
    calib = hpv.Calibration(
        sim, calib_pars=calib_pars, genotype_pars=make_priors(),
        name='nigeria_calib', datafiles=datafiles,
        extra_sim_result_keys=['cancers', 'cancer_incidence', 'asr_cancer_incidence'],
        total_trials=n_trials, n_workers=n_workers, storage=storage,
    )
    calib.calibrate()
    if do_plot:
        calib.plot(do_save=True, fig_path='figures/nigeria_calib.png')
    if do_save:
        sc.saveobj(f'raw_results/nigeria_calib{filestem}.obj', calib)
    print(f'Best pars: {calib.best_pars}')
    return sim, calib


def load_calib(do_plot=True, which_pars=0, save_pars=True, filestem=''):
    calib = sc.load(f'raw_results/nigeria_calib{filestem}.obj')
    if do_plot:
        ut.set_font()
        fig = calib.plot(res_to_plot=200, plot_type='sns.boxplot', do_save=False)
        fig.suptitle('Calibration results, Nigeria')
        fig.tight_layout()
        fig.savefig(f'figures/nigeria_calib{filestem}.png')
    if save_pars:
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=which_pars)
        trial_pars = sc.autolist()
        for i in range(100):
            trial_pars += calib.trial_pars_to_sim_pars(which_pars=i)
        sc.save(f'results/nigeria_pars{filestem}.obj', calib_pars)
        sc.save(f'results/nigeria_pars{filestem}_all.obj', trial_pars)
    return calib


if __name__ == '__main__':
    T = sc.timer()
    if 'run_calibration' in to_run:
        run_calib(n_trials=n_trials, n_workers=n_workers, do_save=do_save)
    if 'plot_calibration' in to_run:
        load_calib(do_plot=True, save_pars=True)
    T.toc('Done')
