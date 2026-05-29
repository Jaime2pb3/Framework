# pfv_robust — Robust statistics & honest reporting for Isomeric Polarization

Companion layer for `pfv_molecule_v11.x` (TwoQuarks / PfV). It does **not**
change how the model probes an LLM; it fixes how the model *measures* and
*reports* what it found. Pure Python standard library — no `numpy`/`scipy`,
so it runs anywhere (including CI without the scientific stack).

## Why

A re-audit of the v11.4b Anthropic runs (`C5 reasoning_drift`, 15 and 5 views)
found five estimation/calibration problems. All are reproduced by
`diagnose_runs.py` against the bundled `data/*.json`:

| # | Problem | Evidence from the runs |
|---|---------|------------------------|
| 1 | **ΔL3 is noise-limited.** `ΔL3 = L3_exp − L3_base` subtracts two noisy estimators. | Baseline carries **55%** (15v) and **69%** (5v) of the ΔL3 variance budget. |
| 2 | **Headline ρ is not robust.** | ρ = **−0.672** (15v) vs **+0.378** (5v) — *sign flip* across runs. |
| 3 | **p-values ignore autocorrelation.** Naive p uses n = 17. | lag-1 autocorr 0.375 → **effective n ≈ 7.7**. |
| 4 | **Γ_CP is a phantom number.** Reported even with zero recovery events. | `gamma_cp` = drift_velocity / 0.02 (the floor); carries no recovery info. |
| 5 | **Wrong-sign "drift".** A large negative ρ on a *drift* case = exp/base *converging*. | 15v flagged **ANTI_DRIFT** (stabilising), not drift. |

Plus: the corrector **can never fire** — `WARN = 0.39` but observed
`max_stress` was 0.34 / 0.30, so `n_corrections = 0` and the whole
"correction effectiveness" half of the model goes unmeasured.

## What it provides (`ip_robust_stats.py`)

- `rho_with_ci(series)` — ρ with **percentile-bootstrap CI** + **permutation
  p-value** (valid under autocorrelation) + **effective-n** context. The
  significance gate is permutation-based, so a genuine monotone trend is not
  penalised for being autocorrelated.
- `baseline_noise_share(l3_exp, l3_base)` — fraction of ΔL3 variance that is
  just baseline; flags `noise_limited` when ≥ 40%.
- `drift_sign_label(...)` — `DRIFT` / `ANTI_DRIFT` / `INCONCLUSIVE`, CI-aware.
- `honest_gamma_cp(drift_velocity, recovery_signals)` — returns
  `gamma_cp = None` and `UNOBSERVED_RECOVERY` when no real recovery event
  exists. Never a misleading number.
- `warn_threshold_percentile(stress_sample, q=0.90)` — data-driven WARN gate
  (P90 of observed stress ≈ **0.25** here) instead of an unreachable 0.39.
- `standardize_metric(value, reference)` — view-count-invariant z-score of a
  flavor against a per-N null distribution. The principled fix for the
  v11.4b "re-weight `BOTTOM_WEIGHTS`" symptom: standardize instead of
  re-weight, so Down/Strange/Up are comparable across view counts.
- `assess_series(...)` — one-call bundle returning the full robust verdict.

## Run it

```bash
# self-test
python3 pfv_robust/ip_robust_stats.py

# re-audit the bundled runs (or pass your own result JSONs)
python3 pfv_robust/diagnose_runs.py
python3 pfv_robust/diagnose_runs.py path/to/result.json --md REPORT.md
```

`REPORT.md` in this folder is a checked-in sample of that output.

## How to wire it into the core model

1. **Reporting:** after a probe completes, call `assess_series(...)` with the
   `delta_l3_series`, `l3_exp`/`l3_base`, `drift_velocity`, the real
   `antibottom_samples` R-values, and `stress_series`. Store the returned
   verdict next to `spearman_rho` in the result dict; treat `significant` and
   `drift_label` as the headline, not the raw ρ/p.
2. **Γ_CP:** replace the floored `recovery_velocity` computation in
   `BottomAntibottomEngine.summarize` with `honest_gamma_cp`.
3. **WARN gate:** seed `CORRECTOR_THRESHOLD` from
   `warn_threshold_percentile` over a held-out stress sample per
   (provider, case-family), so the corrector can actually be exercised.
4. **Flavor comparability:** standardize each flavor with `standardize_metric`
   against a per-view-count null before feeding `BottomOrchestrator`, and
   retire the per-N `BOTTOM_WEIGHTS` re-calibration.

The underlying probe behaviour (Principle 1: detector is passive) is
untouched.
