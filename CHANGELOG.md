# Changelog

## 2026-10-08 — port to hpvsim v3.2
- Ported `run_sim.py`, `run_calibration.py`, `run_scenarios.py`, `plot_scenarios.py`
  and tests to the v3.2 API (Kazakhstan-style `network_pars()`, `ss.MultiSim`,
  `sim.interventions[...]`, `ss.uids` eligibility returns, `intv.screened` /
  `ti_screened` for screening history).
- Recalibrated under v3.2: 2000-trial Optuna run, nested calib pars with widened
  bounds; baseline ASR 19.2 vs Globocan target 18.4 (ratio 1.04). Two-tier
  artifacts: full calib in `raw_results/` (gitignored), shrunken top-100 in
  `results/` (committed).
- Redesigned scenarios as three arms: no interventions, status quo, WHO 90-70-90
  by 2030 (linear 2025 → 2030 scale-up to 90% vax / 70% lifetime screening /
  90% treatment cascade).

## 2026-06-05 — v2.3.0 modernization
- Stood up a clean Nigeria repo on hpvsim v2.3.0; Nigeria parameters folded into
  `make_sim` (single-country convention), using hpvsim default mixing/initialisation.
- Ported data and calibration framework from `hpvsim_pxv_younger`.
- Validated the existing calibration under v2.3.0 (baseline ASR ≈ 17.9 vs target 18.4)
  — no recalibration required.
- Added Baseline/WHO scale-up scenarios, smoke + validation tests, CI, README,
  CHANGELOG, and packaging.
