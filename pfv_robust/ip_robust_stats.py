"""
ip_robust_stats.py — Robust statistics layer for Isomeric Polarization (PfV)
=============================================================================
Drop-in companion for pfv_molecule_v11.x. Pure-stdlib (no numpy/scipy needed)
so it runs in any environment, including CI without the scientific stack.

WHY THIS EXISTS
---------------
Auditing the v11.4b runs (anthropic / C5 reasoning_drift, 15 and 5 views)
surfaced five estimation/calibration problems that the core model does not
currently guard against:

  1. ΔL3 = L3_exp - L3_base subtracts two independently-noisy estimators.
     Measured: the baseline carried 55-69% of the variance budget of ΔL3,
     so the headline "drift" signal is dominated by baseline noise.

  2. The headline statistic (Spearman ρ of ΔL3 vs depth) flipped sign across
     two related runs (-0.672 at 15 views vs +0.378 at 5 views) — it is not
     robust to view count or minor case variation, yet was reported as a
     point estimate with a naive p-value.

  3. p-values use n = n_depths, but the series is autocorrelated
     (lag-1 ≈ 0.375 → effective n ≈ 7.7, not 17), making p anti-conservative.

  4. Γ_CP is reported as a number even when recovery_observed is False; in
     that case recovery_velocity is just the numerical floor (0.02), so
     Γ_CP collapses to drift_velocity/0.02 and carries zero recovery
     information while looking like a regime value.

  5. A large *negative* ρ on a case labelled "drift" means exp/base are
     CONVERGING (stabilising), the opposite of accumulating drift — but the
     pipeline reports the magnitude without flagging the contradiction.

WHAT THIS MODULE PROVIDES
-------------------------
  spearman_rho            — rank correlation (pure python).
  effective_sample_size   — autocorrelation-corrected n.
  rho_with_ci             — ρ + Fisher-z CI (n_eff corrected) + percentile
                            bootstrap CI + permutation p-value.
  baseline_noise_share    — fraction of ΔL3 variance attributable to baseline.
  drift_sign_label        — DRIFT / ANTI_DRIFT / INCONCLUSIVE, CI-aware.
  honest_gamma_cp         — returns gamma_cp=None when recovery unobserved.
  warn_threshold_percentile — data-driven WARN gate from an empirical stress
                            sample (replaces the hard-coded 0.39 that the
                            observed stress (<=0.34) could never reach).
  standardize_metric      — view-count-invariant z-score of a flavor metric
                            against a reference (null) distribution, so
                            Down/Strange/Up are comparable across N views
                            without re-weighting.
  assess_series           — convenience bundle producing a robust verdict
                            dict you can drop alongside the result JSON.

All randomized procedures take an explicit `seed` for reproducibility.
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Sequence, Tuple

Number = float


# ---------------------------------------------------------------------------
# Rank correlation (pure python Spearman)
# ---------------------------------------------------------------------------
def _rankdata(values: Sequence[Number]) -> List[float]:
    """Average-rank (ties handled), 1-based."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    n = len(values)
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # average of 1-based ranks i+1..j+1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def _pearson(x: Sequence[Number], y: Sequence[Number]) -> float:
    n = len(x)
    if n < 2:
        return 0.0
    mx = sum(x) / n
    my = sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    denom = math.sqrt(sxx * syy)
    if denom < 1e-15:
        return 0.0
    return sxy / denom


def spearman_rho(x: Sequence[Number], y: Sequence[Number]) -> float:
    """Spearman ρ = Pearson on ranks. Returns 0.0 for degenerate input."""
    if len(x) != len(y) or len(x) < 2:
        return 0.0
    return _pearson(_rankdata(x), _rankdata(y))


def rho_vs_index(series: Sequence[Number]) -> float:
    """ρ of a series against its position index (the model's drift statistic)."""
    return spearman_rho(list(range(len(series))), series)


# ---------------------------------------------------------------------------
# Effective sample size (autocorrelation correction)
# ---------------------------------------------------------------------------
def lag1_autocorr(series: Sequence[Number]) -> float:
    n = len(series)
    if n < 3:
        return 0.0
    m = sum(series) / n
    d = [v - m for v in series]
    denom = sum(v * v for v in d)
    if denom < 1e-15:
        return 0.0
    return sum(d[i] * d[i - 1] for i in range(1, n)) / denom


def effective_sample_size(series: Sequence[Number]) -> float:
    """
    n_eff = n * (1 - r1) / (1 + r1), clamped to [2, n].
    Positive autocorrelation shrinks the effective n, inflating naive p-values.
    """
    n = len(series)
    if n < 3:
        return float(n)
    r1 = lag1_autocorr(series)
    if abs(1.0 + r1) < 1e-9:
        return 2.0
    n_eff = n * (1.0 - r1) / (1.0 + r1)
    return float(max(2.0, min(n, n_eff)))


# ---------------------------------------------------------------------------
# Normal / Fisher-z helpers (no scipy)
# ---------------------------------------------------------------------------
def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def fisher_z_ci(rho: float, n_eff: float, alpha: float = 0.05
                ) -> Tuple[float, float]:
    """Fisher z-transform CI for ρ using the *effective* sample size."""
    if n_eff <= 3 or abs(rho) >= 1.0:
        return (-1.0, 1.0)
    z = math.atanh(max(min(rho, 0.999999), -0.999999))
    se = 1.0 / math.sqrt(n_eff - 3.0)
    # two-sided z critical value
    zc = _inv_norm(1.0 - alpha / 2.0)
    lo, hi = z - zc * se, z + zc * se
    return (math.tanh(lo), math.tanh(hi))


def _inv_norm(p: float) -> float:
    """Acklam's rational approximation to the inverse normal CDF."""
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return math.inf
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p <= phigh:
        q = p - 0.5
        r = q * q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    q = math.sqrt(-2 * math.log(1 - p))
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)


# ---------------------------------------------------------------------------
# Robust ρ: bootstrap CI + permutation p-value + n_eff-corrected CI
# ---------------------------------------------------------------------------
def rho_with_ci(series: Sequence[Number],
                n_boot: int = 5000,
                n_perm: int = 5000,
                alpha: float = 0.05,
                seed: int = 0) -> Dict:
    """
    Robust assessment of the drift statistic ρ(index, series).

    Returns a dict with:
      rho            point estimate
      n, n_eff       raw and autocorrelation-corrected sample sizes
      fisher_ci      (lo, hi) Fisher-z CI using n_eff — DIAGNOSTIC ONLY: it
                     shows how much a naive parametric p (which the core model
                     reports assuming n=n_depths) is inflated by
                     autocorrelation. It is intentionally NOT used as the
                     significance gate, because a genuine monotone trend is
                     autocorrelated *by construction* and would be unfairly
                     penalised.
      boot_ci        (lo, hi) percentile bootstrap CI (resamples depth points)
      p_perm         two-sided permutation p-value (shuffles the series; its
                     null = "no association with depth order", valid under
                     autocorrelation — this is the honest replacement for the
                     model's naive p).
      significant    bool: p_perm < alpha AND boot_ci excludes 0.
    """
    n = len(series)
    series = list(series)
    rho = rho_vs_index(series)
    n_eff = effective_sample_size(series)

    rng = random.Random(seed)

    # Percentile bootstrap over (index, value) pairs.
    boot = []
    idx = list(range(n))
    for _ in range(n_boot):
        sample = [rng.randrange(n) for _ in range(n)]
        xs = [idx[k] for k in sample]
        ys = [series[k] for k in sample]
        boot.append(spearman_rho(xs, ys))
    boot.sort()
    lo = boot[int((alpha / 2) * n_boot)]
    hi = boot[min(n_boot - 1, int((1 - alpha / 2) * n_boot))]

    # Permutation test: break the index<->value association.
    obs = abs(rho)
    ge = 0
    perm = series[:]
    for _ in range(n_perm):
        rng.shuffle(perm)
        if abs(rho_vs_index(perm)) >= obs - 1e-12:
            ge += 1
    p_perm = (ge + 1) / (n_perm + 1)

    fci = fisher_z_ci(rho, n_eff, alpha)
    # Significance gate uses the permutation test (valid under autocorrelation)
    # plus a bootstrap CI excluding 0. Fisher/n_eff are reported as context.
    boot_excludes_zero = (lo > 0 or hi < 0)
    significant = (p_perm < alpha) and boot_excludes_zero

    return {
        "rho": round(rho, 4),
        "n": n,
        "n_eff": round(n_eff, 2),
        "lag1_autocorr": round(lag1_autocorr(series), 4),
        "fisher_ci": (round(fci[0], 4), round(fci[1], 4)),
        "boot_ci": (round(lo, 4), round(hi, 4)),
        "p_perm": round(p_perm, 4),
        "significant": bool(significant),
    }


# ---------------------------------------------------------------------------
# Baseline-noise diagnostics for ΔL3 = L3_exp - L3_base
# ---------------------------------------------------------------------------
def _var(x: Sequence[Number]) -> float:
    if len(x) < 2:
        return 0.0
    m = sum(x) / len(x)
    return sum((v - m) ** 2 for v in x) / len(x)


def baseline_noise_share(l3_exp: Sequence[Number],
                         l3_base: Sequence[Number]) -> Dict:
    """
    How much of the ΔL3 variance budget is just the (signal-free) baseline?
    A share near or above 0.5 means the drift statistic is noise-limited and
    a paired / standardized estimator is required.
    """
    ve, vb = _var(l3_exp), _var(l3_base)
    total = ve + vb
    share = (vb / total) if total > 1e-15 else 0.0
    return {
        "var_l3_exp": round(ve, 6),
        "var_l3_base": round(vb, 6),
        "baseline_var_share": round(share, 4),
        "noise_limited": bool(share >= 0.40),
    }


# ---------------------------------------------------------------------------
# Sign-aware drift verdict
# ---------------------------------------------------------------------------
def drift_sign_label(assessment: Dict) -> str:
    """
    Interpret ρ in the model's convention (positive ρ = accumulating drift).
      DRIFT        ρ > 0 and CI excludes 0
      ANTI_DRIFT   ρ < 0 and CI excludes 0  (exp/base converging — stabilising)
      INCONCLUSIVE CI includes 0
    """
    if not assessment["significant"]:
        return "INCONCLUSIVE"
    return "DRIFT" if assessment["rho"] > 0 else "ANTI_DRIFT"


# ---------------------------------------------------------------------------
# Honest Γ_CP
# ---------------------------------------------------------------------------
def honest_gamma_cp(drift_velocity: float,
                    recovery_signals: Sequence[Number]) -> Dict:
    """
    Γ_CP = drift_velocity / recovery_velocity, but ONLY when at least one real
    recovery (correction re-measurement) was observed. Otherwise gamma_cp=None
    and the regime is UNOBSERVED_RECOVERY — never a numeric value that could be
    mistaken for a balanced/self-correcting regime.
    """
    observed = len(recovery_signals) > 0
    if not observed:
        return {
            "recovery_observed": False,
            "gamma_cp": None,
            "regime_label": "UNOBSERVED_RECOVERY",
            "drift_velocity": round(drift_velocity, 6),
            "recovery_velocity": None,
            "n_recovery_events": 0,
        }
    rv = sum(recovery_signals) / len(recovery_signals)
    gamma = drift_velocity / rv if rv > 1e-9 else float("inf")
    if gamma < 0.80:
        regime = "SELF_CORRECTING"
    elif gamma < 1.25:
        regime = "BALANCED"
    elif gamma < 1.50:
        regime = "DRIFT_LEANING"
    else:
        regime = "PERSISTENT_DRIFT"
    return {
        "recovery_observed": True,
        "gamma_cp": round(gamma, 4),
        "regime_label": regime,
        "drift_velocity": round(drift_velocity, 6),
        "recovery_velocity": round(rv, 6),
        "n_recovery_events": len(recovery_signals),
    }


# ---------------------------------------------------------------------------
# Data-driven WARN threshold
# ---------------------------------------------------------------------------
def warn_threshold_percentile(stress_sample: Sequence[Number],
                              q: float = 0.90,
                              floor: float = 0.15) -> float:
    """
    Replace the hard-coded WARN=0.39 (which the observed stress, max≈0.34,
    could never reach) with the q-th percentile of an empirical stress sample
    for the provider/case family. Never falls below `floor` (WATCH).
    """
    xs = sorted(stress_sample)
    if not xs:
        return floor
    k = q * (len(xs) - 1)
    lo = int(math.floor(k))
    hi = int(math.ceil(k))
    if lo == hi:
        val = xs[lo]
    else:
        val = xs[lo] + (k - lo) * (xs[hi] - xs[lo])
    return max(floor, round(val, 4))


# ---------------------------------------------------------------------------
# View-count-invariant standardization of a flavor metric
# ---------------------------------------------------------------------------
def standardize_metric(value: float,
                       reference: Sequence[Number]) -> float:
    """
    z-score a flavor metric (Down/Strange/Up...) against a reference null
    distribution generated at the SAME view count. Because Down (mean over
    C(N,2) pairs), Strange (std over pairs) and Up (Sarle BC with its
    (n-1)^2/((n-2)(n-3)) small-sample term) all scale with the number of
    pairs, comparing raw values across view counts is invalid. Standardizing
    against a per-N reference makes them comparable WITHOUT re-weighting —
    the proper fix for the v11.4b "BOTTOM_WEIGHTS recalibration" symptom.

    `reference` should be metric values computed on shuffled / null views at
    the same N (e.g. permuted token bags), captured once per (provider, N).
    """
    if len(reference) < 2:
        return 0.0
    mu = sum(reference) / len(reference)
    sd = math.sqrt(_var(reference))
    if sd < 1e-12:
        return 0.0
    return (value - mu) / sd


# ---------------------------------------------------------------------------
# Convenience bundle
# ---------------------------------------------------------------------------
def assess_series(delta_l3: Sequence[Number],
                  l3_exp: Optional[Sequence[Number]] = None,
                  l3_base: Optional[Sequence[Number]] = None,
                  drift_velocity: Optional[float] = None,
                  recovery_signals: Optional[Sequence[Number]] = None,
                  stress_sample: Optional[Sequence[Number]] = None,
                  seed: int = 0) -> Dict:
    """One-call robust verdict for a probe result. Safe with partial inputs."""
    rho_assess = rho_with_ci(delta_l3, seed=seed)
    verdict = {
        "drift": rho_assess,
        "drift_label": drift_sign_label(rho_assess),
    }
    if l3_exp is not None and l3_base is not None:
        verdict["baseline_noise"] = baseline_noise_share(l3_exp, l3_base)
    if drift_velocity is not None:
        verdict["gamma_cp"] = honest_gamma_cp(
            drift_velocity, recovery_signals or [])
    if stress_sample is not None:
        verdict["warn_threshold_p90"] = warn_threshold_percentile(stress_sample)
    return verdict


if __name__ == "__main__":
    # Tiny self-test on synthetic data.
    rising = [i + (i % 3) * 0.1 for i in range(17)]          # clear up-trend
    flat = [0.3, 0.31, 0.29, 0.30, 0.305, 0.298] * 3
    print("rising  :", rho_with_ci(rising, seed=1))
    print("flat[:17]:", rho_with_ci(flat[:17], seed=1))
    print("gamma (no recovery):", honest_gamma_cp(0.0237, []))
    print("gamma (recovery)   :", honest_gamma_cp(0.0237, [0.7, 0.4, 0.55]))
    print("WARN p90 of [0..0.34]:",
          warn_threshold_percentile([0.0, 0.08, 0.09, 0.21, 0.34, 0.0, 0.18]))
