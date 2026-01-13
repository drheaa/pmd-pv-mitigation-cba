# Methodology

This document explains what this project is doing and why, in plain language. The goal is to keep the workflow clear and reproducible, even if you are still learning power systems.

---

## 1) What problem are we studying?

Low-voltage (LV) distribution networks were originally designed for power to flow one way:
from the transformer to houses.

With lots of rooftop PV, that assumption breaks.

Typical issues in LV networks with high PV include:

- **Voltage rise**: PV injects power locally, pushing voltages above allowed limits.
- **Phase unbalance (VUF)**: PV and loads are often single-phase and unevenly spread, creating unbalanced voltages.
- **Local congestion**: the “bad” behaviour usually starts near the end of feeders and spreads inward.

This project models those effects using realistic feeders and realistic customer load allocations.

---

## 2) What datasets are we using?

We use public datasets:

1. **Six UK-based Low Voltage Distribution Network Models**
   - provides feeder network models (topology, lines, transformers, phases)

2. **D-Suite Modelling: Customer Load Profile Allocations**
   - provides customer load profile allocations so loads are time-varying and realistic

We treat these as read-only inputs and do not modify the raw data.

---

## 3) What does OPF mean here?

**OPF = Optimal Power Flow**

In this project, OPF is a tool to answer:

> Given a network, loads, and PV (and possibly mitigation devices), can the network operate within limits?

and if yes,

> What voltages happen, and how much device effort (reactive power, etc.) is needed?

Important: we are not doing real-time control design.
We are doing planning-style analysis: feasibility + constraint pressure + comparison across scenarios.

We run **unbalanced AC OPF** because LV networks are unbalanced (phases are not identical).

---

## 4) What are the “limits” we care about?

### 4.1 Voltage magnitude limits
Each bus and phase has a voltage magnitude (in per-unit).

We check whether voltage is within bounds, for example:

- lower bound: 0.94 pu
- upper bound: 1.10 pu

(Exact values can be configured.)

### 4.2 Voltage Unbalance Factor (VUF)
VUF measures how unbalanced the voltage is.

A common definition is:

- VUF = |negative sequence| / |positive sequence|

Higher VUF = more unbalance.

We typically compare against a limit like 2% (0.02).

Why it matters:
- High VUF can cause heating and problems for motors and equipment.
- VUF issues can appear even when voltage magnitudes look fine.

---

## 5) What is a “scenario”?

A scenario is one specific case we run, defined by:

- which feeder we are using
- load level or time slice (if doing time series)
- PV penetration or PV size rule
- PV placement rule (which buses get PV)
- whether mitigation devices are included and how they are sized

Each scenario has a `scenario_id` so results are comparable and traceable.

---

## 6) Workflow stages

### Stage A — Baseline (loads only)
Purpose:
- verify the feeder parses and solves
- establish the baseline voltage + VUF behaviour

Outputs:
- baseline bus metrics (voltages, margins, VUF)
- baseline plots

---

### Stage B — PV integration studies
Purpose:
- add PV in a realistic way
- increase PV penetration / size until limits start to get tight or violated

PV modelling approach (planning-level):
- PV is modelled as generation at customer buses.
- It can be single-phase or three-phase depending on assumptions.
- PV size can follow a rule (fixed size, penetration fraction, or dataset-informed).

Outputs:
- PV sweep curves (how vmin/vmax/VUF change with PV level)
- a set of “stress scenarios” to carry forward

---

### Stage C — Congestion identification (PV-driven)
Purpose:
- locate where the network first becomes problematic
- quantify severity
- define a “congestion region” along the feeder

What we compute for each scenario:
- min/max voltage (per bus and phase)
- number of buses/phases violating voltage bounds
- VUF per bus and number violating VUF limit
- explicit margins (how close to limits even if not violating)

Congestion region logic (simple):
- find buses with violations (or very small margins)
- map them by distance from the source
- describe the contiguous area near the end (or any cluster) as the region

Outputs:
- `scenario_summary.csv`
- `scenario_bus_metrics__<scenario_id>.csv`
- congestion plots (voltage, margin, VUF vs distance)

---

### Stage D — STATCOM mitigation
Purpose:
- test whether reactive power support can reduce PV-driven voltage rise and/or reduce unbalance pressure

What a STATCOM is (in this project):
- a device that injects or absorbs reactive power (Q)
- it does not supply real energy
- it can be modelled as a controllable generator with P ≈ 0 and bounded Q

In the OPF:
- the solver chooses Q injection (often per phase)
- we check whether voltage and/or VUF improves
- we log whether STATCOM hits its Q limits (binding)

Outputs:
- “before vs after” plots and metrics
- STATCOM dispatch logs (q per phase)
- summary table for CBA

---

### Stage E — SOP mitigation (later stage)
Purpose:
- test stronger flexibility: shifting power between locations (and often also providing reactive support)

An SOP acts like a controllable link between two points, helping:
- move power away from constrained areas
- reduce voltage rise/unbalance indirectly

We only do SOP after STATCOM because SOP is generally more complex and costly.

---

### Stage F — Battery mitigation (later stage)
Purpose:
- absorb PV during peaks and reduce export
- provide time-shifting value (especially if time series is included)

In a snapshot study, batteries can still be modelled as:
- charging (absorbing power) at stressed times
- discharging at other times if running multiple time steps

---

### Stage G — Cost–Benefit Analysis (CBA)
Purpose:
- compare mitigation options fairly

Inputs:
- technical metrics from each mitigation run (violations reduced, margins increased, device effort)
- cost assumptions from `configs/costs.toml`

Outputs:
- cost vs benefit tables
- recommendation: which option provides the best value for the chosen stress scenarios

---

## 7) Reproducibility rules

- Raw datasets are stored in `data/raw/` and never edited.
- All assumptions live in `configs/`.
- Scripts write outputs into `results/<stage>/tables` and `results/<stage>/figures`.
- Downstream scripts read upstream CSVs instead of re-deriving logic manually.

This keeps the pipeline clean and prevents “messy glue” later.

---
