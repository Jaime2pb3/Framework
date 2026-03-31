# TwoQuarks — Quark Flavors

**A Modular Framework for Adaptive Stability Control in Sequence Models Under Regime Uncertainty**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

→ **[twoquarks.com](https://twoquarks.com)** | **[Preprint](https://twoquarks.com/preprint.pdf)**

---

## Overview

Six modular control mechanisms — quark flavors — that monitor internal pre-instability signals and apply targeted interventions at inference time, without modifying model parameters, policies, or training objectives.

Each flavor operates independently and addresses a distinct class of instability.

---

## Flavors

### ↓ DOWN — Early Tension
Monitors TD-variability among functionally similar states via the **Neighborhood TD-Variability Monitor (NTVM)**. When states that should agree start producing inconsistent learning signals, something is already wrong — before any external metric shows it.

### ~ STRANGE — Hidden Shifts
Tracks swarm disagreement across policies operating on the same input via the **Swarm Disagreement & Regime Boundary Detector (SDRBD)**. When the swarm splits, a regime boundary is nearby — even if the environment gives no explicit signal.

### ↑ UP — Variance Awareness
Penalizes internal disagreement across ensemble value estimates via the **Variance-Penalized Ensemble Controller (VPEC)**. Confidence without consensus is a warning, not a signal to act on.

### ◎ CHARM — Coherence Under Drift
Preserves internal coherence under gradual structural erosion via the **Potential-Based Coherence Stabilizer (PBCS)**. Charm doesn't react to collapse — it makes collapse harder to reach.

### ⊤ TOP — Transient Coupling
Under extreme stress, antagonistically coupled policies briefly synchronize — a state called **toponium**. Top activates at that moment, then disengages. It exists only in the window between pressure and collapse.

### ⊥ BOTTOM — Stress Aggregation
Defines collapse events, records instability signal activations, and aggregates contributions from all other mechanisms into a single stress indicator via the **Stability Instrumentation Layer (SIL)**.

---

## Experimental Results

| Agent | Mean Return | Failure Rate | P95 Return |
|-------|------------|-------------|------------|
| HF Levo | 0.538 | 0.0817 | 2.00 |
| LevoParadoxIsomer | 0.546 | 0.0833 | 2.35 |
| LevoParadoxPPOHybrid | 0.256 | 0.1400 | 1.10 |

Evaluated over 1,200 episodes in the Paradox environment (deceptive reward signals + delayed feedback).

---

## Structure

```
quarks/
├── down/        # NTVM — TD-variability, paradox environments
├── strange/     # SDRBD — swarm disagreement, regime detection
├── up/          # VPEC — variance-penalized ensemble, bellman shift
├── charm/       # PBCS — coherence stabilizer, enchanted valley
├── top/         # BSPC/ARCC — toponium, gauge field environments
└── bottom/      # SIL — stress aggregation, baseline
```

---

## Citation

```bibtex
@techreport{ledesma2026twoquarks,
  author      = {Ledesma, Luis Jaime},
  title       = {A Modular Framework for Adaptive Stability Control in Sequence Models Under Regime Uncertainty},
  institution = {TwoQuarks Research},
  year        = {2026},
  url         = {https://twoquarks.com/preprint.pdf}
}
```

---

**TwoQuarks Research** · [twoquarks.com](https://twoquarks.com) · research@twoquarks.com
