"""
diagnose_runs.py — Re-audit pfv_molecule_v11.4b result JSONs with robust stats.

Loads one or more PfV result JSON files and, for each probe result, prints
(and optionally writes to markdown) a before/after comparison:

  • the headline statistic the model reported (spearman_rho / p / gamma_cp /
    regime), vs
  • the robust verdict from ip_robust_stats (CI-aware ρ, permutation p,
    effective n, baseline-noise share, sign-aware drift label, honest Γ_CP,
    data-driven WARN threshold).

Usage
-----
  python diagnose_runs.py                      # uses bundled data/*.json
  python diagnose_runs.py path/to/result.json ...
  python diagnose_runs.py --md report.md       # also write markdown report
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ip_robust_stats as rs  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_FILES = [
    os.path.join(HERE, "data", "run_15views.json"),
    os.path.join(HERE, "data", "run_5views.json"),
]


def extract_recovery_signals(result: Dict) -> List[float]:
    """Real recovery events only (corrector actually re-measured)."""
    bab = result.get("bab_log", {})
    samples = bab.get("antibottom_samples", []) or []
    return [s.get("R", 0.0) for s in samples]


def diagnose_result(top: Dict, result: Dict, seed: int = 0) -> Dict:
    n_views = top.get("n_views")
    states = result.get("molecule_states", [])
    l3_exp = [s.get("l3_exp", s.get("l3", 0.0)) for s in states]
    l3_base = [s.get("l3_base", 0.0) for s in states]
    delta = result.get("delta_l3_series", [])
    stress = result.get("stress_series", [])

    assess = rs.assess_series(
        delta_l3=delta,
        l3_exp=l3_exp,
        l3_base=l3_base,
        drift_velocity=result.get("drift_velocity", 0.0),
        recovery_signals=extract_recovery_signals(result),
        stress_sample=stress,
        seed=seed,
    )
    return {
        "case_id": result.get("case_id"),
        "n_views": n_views,
        "reported": {
            "spearman_rho": result.get("spearman_rho"),
            "spearman_p": result.get("spearman_p"),
            "gamma_cp": result.get("gamma_cp"),
            "regime_label": result.get("regime_label"),
            "max_bottom_state": result.get("max_bottom_state"),
            "n_corrections": result.get("n_corrections"),
            "max_stress": round(max(stress), 4) if stress else None,
        },
        "robust": assess,
    }


def fmt_block(d: Dict) -> List[str]:
    r = d["reported"]
    rob = d["robust"]
    drift = rob["drift"]
    bn = rob.get("baseline_noise", {})
    g = rob.get("gamma_cp", {})
    lines = []
    lines.append(f"### {d['case_id']}  (n_views={d['n_views']})")
    lines.append("")
    lines.append("**Reported by model:**")
    lines.append(f"- spearman_rho = {r['spearman_rho']}  (p = {r['spearman_p']})")
    lines.append(f"- gamma_cp = {r['gamma_cp']}  → regime = {r['regime_label']}")
    lines.append(f"- max_bottom_state = {r['max_bottom_state']}  | "
                 f"corrections = {r['n_corrections']}  | "
                 f"max_stress = {r['max_stress']}")
    lines.append("")
    lines.append("**Robust verdict:**")
    lines.append(f"- ρ = {drift['rho']}  | boot 95% CI = {drift['boot_ci']}  | "
                 f"permutation p = {drift['p_perm']}  → "
                 f"**significant = {drift['significant']}**")
    lines.append(f"- n = {drift['n']} but effective n ≈ {drift['n_eff']} "
                 f"(lag-1 autocorr = {drift['lag1_autocorr']})  → "
                 f"naive p assumed too many independent points")
    lines.append(f"- drift label = **{rob['drift_label']}** "
                 f"(positive ρ = accumulating drift; negative = stabilising)")
    if bn:
        lines.append(f"- baseline noise share of ΔL3 variance = "
                     f"**{bn['baseline_var_share']:.0%}**  "
                     f"(noise-limited = {bn['noise_limited']})")
    if g:
        lines.append(f"- honest Γ_CP = {g['gamma_cp']}  "
                     f"(recovery_observed = {g['recovery_observed']}, "
                     f"events = {g['n_recovery_events']})  → {g['regime_label']}")
    if "warn_threshold_p90" in rob:
        lines.append(f"- data-driven WARN (P90 of observed stress) = "
                     f"**{rob['warn_threshold_p90']}**  "
                     f"(model used fixed 0.39, which max_stress never reached)")
    lines.append("")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description="Robust re-audit of PfV result JSONs")
    ap.add_argument("files", nargs="*", default=None,
                    help="result JSON files (default: bundled data/)")
    ap.add_argument("--md", default=None, help="also write a markdown report")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    files = args.files if args.files else DEFAULT_FILES
    out: List[str] = ["# PfV robust re-audit", ""]
    diagnosed = []
    for fn in files:
        with open(fn) as fh:
            top = json.load(fh)
        for result in top.get("results", []):
            d = diagnose_result(top, result, seed=args.seed)
            diagnosed.append(d)
            out.extend(fmt_block(d))

    # Cross-run note if two reasoning_drift runs are present.
    rhos = [(d["n_views"], d["robust"]["drift"]["rho"],
             d["robust"]["drift"]["significant"]) for d in diagnosed]
    if len(rhos) >= 2:
        signs = {(r > 0) for _, r, _ in rhos}
        out.append("### Cross-run robustness")
        out.append("")
        for nv, r, sig in rhos:
            out.append(f"- n_views={nv}: ρ={r}  significant={sig}")
        out.append("")
        if len(signs) > 1:
            out.append("> ⚠️ ρ **changes sign** across runs — the headline "
                       "statistic is not robust to view count / case variation.")
        out.append("")

    text = "\n".join(out)
    print(text)
    if args.md:
        with open(args.md, "w") as fh:
            fh.write(text)
        print(f"\n[written] {args.md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
