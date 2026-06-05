# Project Engineering Score

- **Project**: `/Users/robynstuart/gf/hpvsim_nigeria`
- **Tier**: 2 (Small-scale project used by multiple people or projects)
- **Overall Score**: 76/100
- **Status**: PASS
- **Date**: 2026-06-05
- **Version**: idm-eng-plugin:eng-quality-checker v1.3_2026.04.13
- **Time spent**: 140s

## Summary

| Category | Score | Weight |
| -- | -- | -- |
| Quality | 75/100 | 40% |
| Usability | 75/100 | 40% |
| Safety | 80/100 | 20% |
| **Total** | **76/100** | 100% |

| Metric | Score | Notes |
| -- | -- | -- |
| correct | 7/10 | 4 tests (pars load, debug smoke, full-res ASR validation vs 18.4, both scenarios); CI runs full suite; no serious bugs. Few edge-case tests. |
| clear | 8/10 | Logical structure, docstrings throughout, README file table, cited constants. |
| concise | 10/10 | No duplication; scenarios cleanly parameterized; libraries used appropriately. |
| simple | 8/10 | Three intuitive entry points with defaults; `to_run` comment pattern is minor friction; no arg validation. |
| powerful | 7/10 | Good kwargs; some assumptions (n_trials, coverage levels) inlined rather than args. |
| performant | 7/10 | Thread-pinning, debug mode, n_workers scaling; no benchmarks. |
| documented | 7/10 | Thorough README + docstrings; no tutorial. |
| accessible | 9/10 | Public repo, MIT license, CHANGELOG, CI, 1-command install; no CONTRIBUTING. |
| compliant | 10/10 | MIT, no secrets, all deps permissive, public data. |
| reproducible | 5/10 | Git + `rand_seed` determinism, version-bounded deps; no version tag or lock file. |

A clean Tier 2 handover repo: modular, documented, MIT-licensed, tested, CI-wired, with the calibration validated under v2.3.0 (ASR 17.9 vs target 18.4). Strongest on compliance and conciseness; biggest gap is reproducibility (no version tag/lock) and the absence of a tutorial. No correctness or compliance failures.

## Recommendations

1. **[reproducible] — Tag a release** *(quick; automated: yes)*
   Add a semantic-version git tag (e.g. `v1.0.0`) once merged.
2. **[correct] — Add edge-case tests** *(medium; automated: yes)*
   Test `_to_annual_prob`/`_layer_probs_to_annual` at boundary values, and a zero-coverage intervention (no effect).
3. **[documented / simple] — Add a short tutorial** *(medium; automated: no)*
   A `tutorial.ipynb` walking through run → scenarios → plot, to aid the handover.
4. **[simple] — CLI flags for calibration modes** *(medium; automated: yes)*
   Replace the `to_run` comment/uncomment pattern in `run_calibration.py` with `argparse`.
5. **[reproducible] — Lock file** *(quick; automated: yes)*
   Add `requirements-lock.txt` (pip freeze of the validated environment).

## Full Results

```yaml
project: /Users/robynstuart/gf/hpvsim_nigeria
tier: 2
overall_score: 76
failed: false
quality:
  correct:    {score: 7, weight: 7, reason: "4 tests covering main workflows incl. full-res ASR validation; CI runs full suite; no serious bugs; few edge-case tests."}
  clear:      {score: 8, weight: 2, reason: "Logical structure, docstrings throughout, README file table, inline-cited constants."}
  concise:    {score: 10, weight: 1, reason: "No duplication; scenarios cleanly parameterized; libraries used appropriately."}
usability:
  simple:     {score: 8, weight: 3, reason: "Three intuitive entry points with defaults; to_run comment pattern minor friction; no arg validation."}
  powerful:   {score: 7, weight: 2, reason: "Good kwargs; some assumptions inlined rather than args."}
  performant: {score: 7, weight: 2, reason: "Thread-pinning, debug mode, n_workers scaling; no benchmarks."}
  documented: {score: 7, weight: 2, reason: "Thorough README + docstrings; no tutorial."}
  accessible: {score: 9, weight: 1, reason: "Public repo, MIT, CHANGELOG, CI, 1-command install; no CONTRIBUTING."}
safety:
  compliant:    {score: 10, weight: 6, reason: "MIT, no secrets, permissive deps, public data."}
  reproducible: {score: 5, weight: 4, reason: "Git + rand_seed determinism, version-bounded deps; no version tag/lock file."}
```
