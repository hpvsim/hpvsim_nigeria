"""Plot cervical cancer incidence (ASR) for Baseline vs WHO scenarios, Nigeria.

Runs each scenario across several seeds and plots the median with a 10-90%
uncertainty band, reproducing the partner's figure format with the
v2.3.0-calibrated model.
"""
import numpy as np
import matplotlib.pyplot as plt
import run_scenarios as rsc
import utils as ut

N_SEEDS = 3
COLORS = {'Baseline': '#c0392b', 'WHO': '#2980b9'}


def _band(msim):
    """Return (year, median, low, high) ASR across the MultiSim's seeds."""
    yr = np.array(msim.sims[0].results['year'])
    arrs = np.array([np.array(s.results['asr_cancer_incidence']) for s in msim.sims])
    return yr, np.median(arrs, 0), np.percentile(arrs, 10, 0), np.percentile(arrs, 90, 0)


def main(n_seeds=N_SEEDS):
    msims = rsc.run_scenarios(end=2100, n_seeds=n_seeds, do_save=True)
    ut.set_font(13)
    fig, ax = plt.subplots(figsize=(9, 5))
    for name, msim in msims.items():
        yr, med, lo, hi = _band(msim)
        m = yr >= 2006
        ax.fill_between(yr[m], lo[m], hi[m], color=COLORS[name], alpha=0.2)
        ax.plot(yr[m], med[m], lw=2.5, color=COLORS[name], label=name)
        print(f'{name}: ASR 2020={med[np.argmin(abs(yr-2020))]:.1f} '
              f'2050={med[np.argmin(abs(yr-2050))]:.1f} 2100={med[np.argmin(abs(yr-2100))]:.1f}')
    ax.set_title('HPVsim v2.3.0: Cervical cancer incidence per 100k women (ASR) — Nigeria',
                 fontsize=13, fontweight='bold')
    ax.set_xlabel('Year')
    ax.set_ylabel('ASR cancer incidence (per 100,000)')
    ax.set_xlim(2006, 2100)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, title=f'(median of {n_seeds} seeds, 10-90% band)')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('figures/nigeria_scenarios_asr.png', dpi=200)
    print('saved figures/nigeria_scenarios_asr.png')


if __name__ == '__main__':
    main()
