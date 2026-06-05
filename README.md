# hpvsim_nigeria

An [HPVsim](https://hpvsim.org) model of cervical cancer for Nigeria. Built on
**hpvsim v2.3.0**. Ported from the `hpvsim_pxv_younger` analysis repo.

## Install

```bash
pip install -r requirements.txt
```

Requires `hpvsim==2.3.0`.

## What's here

| File | Purpose |
|------|---------|
| `run_sim.py` | Defines the Nigeria simulation (`make_sim`, `run_sim`). Nigeria-specific parameters (sexual behaviour, layer probabilities) are inlined here; mixing and initial conditions use hpvsim defaults. |
| `run_calibration.py` | Calibrates the model to Nigeria data (`hpv.Calibration`). |
| `run_scenarios.py` | Baseline vs WHO scale-up scenarios. |
| `utils.py` | Fonts and calibration datafiles. |
| `data/` | Calibration targets (cancer cases, CIN/cancer genotype distributions, ASR, HPV prevalence). |
| `results/nigeria_pars.obj` | Calibrated parameter set (validated under v2.3.0). |
| `tests/` | Smoke + validation tests. |

## Data provenance

- `nigeria_cancer_cases.csv`, `nigeria_asr_cancer_incidence.csv` — cervical cancer
  incidence (IARC/Globocan).
- `nigeria_cancer_types.csv`, `nigeria_cin_types.csv` — HPV genotype distribution in
  cancers / CIN (ICO/IARC HPV Information Centre).
- `nigeria_hpv_prevalence.csv` — HPV/precancer prevalence (validation check).

Data and the calibration framework were ported from `hpvsim_pxv_younger`.

## How to run

```bash
python run_sim.py                 # single baseline run + plot (local)

# Calibration — RUN only on a multi-core VM (edit `to_run` in the file):
python run_calibration.py         # 'plot_calibration' extracts/plots locally;
                                  # 'run_calibration' fits (VM only)

python run_scenarios.py           # Baseline vs WHO comparison (2100 horizon)
```

> **Calibration compute:** the calibration is only fast on multi-core machines. Run the
> `run_calibration` step on a VM; use `plot_calibration` locally to extract the best
> parameters into `results/nigeria_pars.obj`.

## Calibration status

The parameter set in `results/nigeria_pars.obj` reproduces the repo's ASR target
(model ≈ 17.9 vs target 18.4 per 100,000, 2020) under hpvsim v2.3.0 — see
`tests/test_baseline.py`.

> **Note on absolute incidence:** this model is calibrated to cervical-cancer **case
> counts and genotype distributions**, which yield an ASR of ~18/100,000. A separate
> Globocan headline figure for Nigeria is higher (~26/100,000). The difference is a
> choice of calibration target, not a version effect — running the same calibrated
> parameters on hpvsim v2.2.6 and v2.3.0 gives identical results.

## Testing

```bash
pytest tests/                     # full suite, incl. ~1-2 min ASR validation
```
