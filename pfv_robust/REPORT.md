# PfV robust re-audit

### C5_reasoning_drift  (n_views=15)

**Reported by model:**
- spearman_rho = -0.6716  (p = 0.003156)
- gamma_cp = 1.186  → regime = UNOBSERVED_RECOVERY
- max_bottom_state = WATCH  | corrections = 0  | max_stress = 0.341

**Robust verdict:**
- ρ = -0.6716  | boot 95% CI = (-0.8835, -0.2921)  | permutation p = 0.005  → **significant = True**
- n = 17 but effective n ≈ 7.73 (lag-1 autocorr = 0.3749)  → naive p assumed too many independent points
- drift label = **ANTI_DRIFT** (positive ρ = accumulating drift; negative = stabilising)
- baseline noise share of ΔL3 variance = **55%**  (noise-limited = True)
- honest Γ_CP = None  (recovery_observed = False, events = 0)  → UNOBSERVED_RECOVERY
- data-driven WARN (P90 of observed stress) = **0.253**  (model used fixed 0.39, which max_stress never reached)

### C5_reasoning_drift_12t  (n_views=5)

**Reported by model:**
- spearman_rho = 0.3775  (p = 0.135273)
- gamma_cp = 1.1273  → regime = UNOBSERVED_RECOVERY
- max_bottom_state = WATCH  | corrections = 0  | max_stress = 0.3004

**Robust verdict:**
- ρ = 0.3775  | boot 95% CI = (-0.1599, 0.7869)  | permutation p = 0.1348  → **significant = False**
- n = 17 but effective n ≈ 17.0 (lag-1 autocorr = -0.0006)  → naive p assumed too many independent points
- drift label = **INCONCLUSIVE** (positive ρ = accumulating drift; negative = stabilising)
- baseline noise share of ΔL3 variance = **69%**  (noise-limited = True)
- honest Γ_CP = None  (recovery_observed = False, events = 0)  → UNOBSERVED_RECOVERY
- data-driven WARN (P90 of observed stress) = **0.2244**  (model used fixed 0.39, which max_stress never reached)

### Cross-run robustness

- n_views=15: ρ=-0.6716  significant=True
- n_views=5: ρ=0.3775  significant=False

> ⚠️ ρ **changes sign** across runs — the headline statistic is not robust to view count / case variation.
