# Changelog

## 2026-06-05 — v2.3.0 modernization
- Stood up a clean Nigeria repo on hpvsim v2.3.0; Nigeria parameters folded into
  `make_sim` (single-country convention), using hpvsim default mixing/initialisation.
- Ported data and calibration framework from `hpvsim_pxv_younger`.
- Validated the existing calibration under v2.3.0 (baseline ASR ≈ 17.9 vs target 18.4)
  — no recalibration required.
- Added Baseline/WHO scale-up scenarios, smoke + validation tests, CI, README,
  CHANGELOG, and packaging.
