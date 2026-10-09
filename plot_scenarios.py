"""Plot cervical cancer incidence (ASR) across Nigeria scenarios.

Reads a thin CSV (``results/nigeria_scenarios.csv``) if present; otherwise
runs the full scenarios + writes the CSV + plots. The CSV holds just
scenario x year x {median, low, high} for ASR — small enough to commit.
The full MultiSim .obj (~250 MB) stays gitignored and is regeneratable.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import run_scenarios as rsc
import utils as ut

N_SEEDS = 3
CSV_PATH = 'results/nigeria_scenarios.csv'
COLORS = {
    'No interventions':      '#7f8c8d',
    'Status quo':            '#c0392b',
    'WHO 90-70-90 by 2030':  '#2980b9',
}


def _band(msim):
    """Year, median, low (10%), high (90%) ASR across the MultiSim's seeds."""
    annual_sims = [s.results['all_hpv']['asr_cancer_incidence'].annualize() for s in msim.sims]
    yr = np.asarray(annual_sims[0].timevec.years, dtype=float)
    arrs = np.array([np.asarray(a.values) for a in annual_sims])
    return yr, np.median(arrs, 0), np.percentile(arrs, 10, 0), np.percentile(arrs, 90, 0)


def extract_plot_data(msims, out_path=CSV_PATH):
    """Collapse a {scenario: MultiSim} dict into a long-format CSV."""
    rows = []
    for name, msim in msims.items():
        yr, med, lo, hi = _band(msim)
        for i, y in enumerate(yr):
            rows.append(dict(scenario=name, year=float(y),
                             median=float(med[i]), low=float(lo[i]), high=float(hi[i])))
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    df.to_csv(out_path, index=False)
    return df


def _load_plot_data(csv_path=CSV_PATH):
    df = pd.read_csv(csv_path)
    out = {}
    for name, g in df.groupby('scenario', sort=False):
        g = g.sort_values('year')
        out[name] = (g.year.to_numpy(), g['median'].to_numpy(),
                     g.low.to_numpy(), g.high.to_numpy())
    return out


def main(n_seeds=N_SEEDS, msims=None, use_cached=True):
    """Plot ASR trajectories. Prefers the committed CSV; falls back to
    running the full scenarios if no CSV + no msims are supplied."""
    if msims is None and use_cached and os.path.exists(CSV_PATH):
        series = _load_plot_data()
    else:
        if msims is None:
            msims = rsc.run_scenarios(end=2100, n_seeds=n_seeds, do_save=True)
        extract_plot_data(msims)
        series = {name: _band(msim) for name, msim in msims.items()}

    ut.set_font(13)
    fig, ax = plt.subplots(figsize=(9, 5))
    for name, (yr, med, lo, hi) in series.items():
        m = yr >= 2006
        ax.fill_between(yr[m], lo[m], hi[m], color=COLORS.get(name, '#333'), alpha=0.2)
        ax.plot(yr[m], med[m], lw=2.5, color=COLORS.get(name, '#333'), label=name)
        print(f'{name}: ASR 2020={med[np.argmin(abs(yr-2020))]:.1f} '
              f'2050={med[np.argmin(abs(yr-2050))]:.1f} 2100={med[np.argmin(abs(yr-2100))]:.1f}')
    ax.set_title('Cervical cancer incidence per 100k women (ASR) — Nigeria',
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
