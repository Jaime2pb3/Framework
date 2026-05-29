# Framework — Isomeric Polarization (PfV)

Modular stability control for LLMs: measure semantic drift across multi-view
responses (decomposed into 6 quark-analog "flavors"), and **correct it before
it fails** — with a stateless, signal-driven controller (it acts on signal, it
does **not** learn; no parameter or memory is updated at runtime).

## Files

- **`Molecule.py`** — the model, **v12.2**, packaged as a single self-contained
  file (the robust-stats layer is inlined; no external import needed). Requires
  `numpy`, `scipy`, `scikit-learn` for the model itself.
- **`pfv_robust/`** — the standalone robust-stats layer + a re-audit tool:
  - `ip_robust_stats.py` (pure stdlib): bootstrap-CI + permutation-p Spearman,
    effective-n, baseline-noise share, sign-aware drift label, honest Γ_CP,
    data-driven WARN percentile, view-count-invariant standardization.
  - `diagnose_runs.py`, `data/`, `REPORT.md`: before/after re-audit of the
    v11.4b runs. `Molecule.py` embeds a mirror of `ip_robust_stats.py`.

## What v12 adds over v11.4b (probe behaviour unchanged — detector stays passive)

| Area | Change |
|------|--------|
| **WARN gate** | `--warn_threshold` makes the corrector trigger calibrable. The fixed 0.39 was unreachable (observed max stress ~0.30–0.34) so it never fired. Each run reports the P90 data-driven suggestion. |
| **Honest Γ_CP** | `None` + `UNOBSERVED_RECOVERY` when no real recovery event exists, instead of `drift_velocity/floor` (zero recovery information). |
| **Robust reporting** | ρ with bootstrap CI + permutation p + effective-n, baseline-noise share, and a `DRIFT/ANTI_DRIFT/INCONCLUSIVE` label, added to every result. |
| **Closed-loop corrector** | Bounded escalation (`MAX_CORRECTION_ATTEMPTS`): act → regenerate → re-measure effective recovery; if not `CORRECTED`, escalate to a harder reset; stop on success. Stateless. |
| **Precise-moment trigger** | Fires at drift **onset** (rising + accelerating + surprise) as well as on `≥WARN`; suppressed when self-stabilizing (`dPt<0`). `--no_early_warning` to disable. |
| **CRN (v12.2)** | `--paired_seeds`: BASE shares EXP's seed stream to cancel baseline sampling noise in ΔL3. Only useful with **OpenAI** (honors `seed`); no-op on Anthropic, degenerate on mock (a runtime warning is printed). |

## Run

```bash
# flagship case, corrector actually activable, predictive trigger on:
python3 Molecule.py --provider anthropic --cases C5_reasoning_drift \
        --views 15 --seed_base 1054 --warn_threshold 0.25 --plots

# OpenAI with CRN variance reduction on the baseline arm:
python3 Molecule.py --provider openai --cases C5_reasoning_drift \
        --views 15 --warn_threshold 0.25 --paired_seeds

# re-audit existing result JSONs with robust stats:
python3 pfv_robust/diagnose_runs.py path/to/result.json --md REPORT.md
```

Key flags: `--warn_threshold FLOAT`, `--no_early_warning`, `--paired_seeds`,
`--fixed_temp`, `--views`, `--repeats`, `--list_cases`.
