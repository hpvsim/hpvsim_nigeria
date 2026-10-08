# hpvsim_nigeria

An [HPVsim](https://hpvsim.org) model of cervical cancer for Nigeria. Built on
**hpvsim v3.2**.

## Install

```bash
pip install -r requirements.txt
```

## What's here

| File | Purpose |
|------|---------|
| `run_sim.py` | Defines the Nigeria simulation (`make_sim`, `run_sim`), with network pars fit to the 2018 Nigeria DHS. |
| `run_calibration.py` | Calibrates the model to Nigeria data (`hpv.Calibration`, 14-param nested prior). |
| `run_scenarios.py` | Three scenarios: no interventions, status quo, WHO 90-70-90 by 2030. |
| `plot_scenarios.py` | Renders the ASR median + 10-90% band across scenarios. |
| `utils.py` | Fonts and calibration datafiles. |
| `data/` | Calibration targets (cancer cases, CIN/cancer genotype distributions, ASR, HPV prevalence). |
| `results/nigeria_pars.obj` | Best-fit parameter set (v3.2). |
| `results/nigeria_calib.obj` | Shrunken calibration object (top-100 trials). |
| `results/nigeria_pars_all.obj` | Top-100 flat parsets for scenario uncertainty propagation. |
| `tests/` | Smoke + validation tests. |

## Data provenance

- `nigeria_cancer_cases.csv`, `nigeria_asr_cancer_incidence.csv` — cervical cancer
  incidence (IARC/Globocan).
- `nigeria_cancer_types.csv`, `nigeria_cin_types.csv` — HPV genotype distribution in
  cancers / CIN (ICO/IARC HPV Information Centre).
- `nigeria_hpv_prevalence.csv` — HPV/precancer prevalence (validation check).

## How to run

```bash
python run_sim.py               # single baseline run + plot (local)

# Calibration — RUN only on a multi-core VM (edit `to_run` in the file):
python run_calibration.py       # 'plot_calibration' extracts/plots locally;
                                # 'run_calibration' fits (VM only, ~15 min)

python plot_scenarios.py        # 3 scenarios × 3 seeds + ASR figure (~8 min local)
```

> **Calibration compute:** the calibration runs 2000 Optuna trials × 64 workers.
> Fast only on a multi-core VM. Use `plot_calibration` locally to extract the
> best parameters into `results/nigeria_pars.obj`.

## Calibration status

The parameter set in `results/nigeria_pars.obj` reproduces the repo's ASR target
under hpvsim v3.2 (model ≈ 19.2 vs target 18.4 per 100,000, 2020) — see
`tests/test_baseline.py`.

> **Note on absolute incidence:** this model is calibrated to cervical-cancer **case
> counts and genotype distributions**, which yield an ASR of ~18/100,000. A separate
> Globocan headline figure for Nigeria is higher (~26/100,000). The difference is a
> choice of calibration target, not a version effect.

## Scenarios

| Scenario | Vax | Screening | Treatment |
|----------|-----|-----------|-----------|
| No interventions | — | — | — |
| Status quo | Historical 2023-25 ramp to 60%, held | 15% lifetime | 50% cascade coverage |
| WHO 90-70-90 by 2030 | Linear ramp 60→90% by 2030 | Linear ramp 15→70% lifetime | Linear ramp 50→90% |

Linear scale-ups all run 2025 → 2030, then hold through 2100.

## Testing

```bash
pytest tests/                   # full suite incl. baseline ASR validation
```
