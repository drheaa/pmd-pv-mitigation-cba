# LV PV Integration and Mitigation using PowerModelsDistribution

This repository contains research code for analysing the impact of high rooftop PV penetration on low-voltage (LV) distribution networks, and for testing mitigation options such as STATCOMs, SOPs, and batteries.

The work is exploratory and research-focused, and is intended to replicate and extend a D-Suite-style (https://www.spenergynetworks.co.uk/userfiles/file/D-Suite-D1.3-Report-Cost-Benefit-Analysis.pdf) from SP Energy Networks, planning workflow using unbalanced AC optimal power flow (OPF).

This project is being developed as part of a student research project and is not production-grade software.

---

## Motivation

High penetrations of rooftop PV can create operational challenges in LV networks, including:

- voltage rise above statutory limits  
- phase unbalance and high Voltage Unbalance Factor (VUF)  
- spatially localised congestion near the end of feeders  

Distribution planners typically address these issues using power-electronic devices (e.g. STATCOMs, SOPs) or energy storage.

This project aims to:
1. quantify PV-driven congestion in realistic LV feeders  
2. identify where and why violations occur  
3. evaluate the technical effectiveness of mitigation devices  
4. compare mitigation options using a simple cost–benefit framework  

---

## Scope of the study

The project follows this workflow:

1. **Baseline modelling**  
   - Unbalanced LV feeder with customer loads only  
   - AC OPF used to compute voltages and unbalance  

2. **PV integration studies**  
   - Rooftop PV added at customer buses  
   - PV penetration increased to create stress scenarios  
   - Voltage and VUF congestion identified  

3. **Congestion identification**  
   - Voltage magnitude violations  
   - VUF violations  
   - Spatial congestion regions along feeders  

4. **Mitigation modelling**  
   - STATCOM placement and sizing  
   - SOP modelling  
   - Battery-based mitigation  

5. **Cost–benefit analysis (CBA)**  
   - Comparison of mitigation options under identical stress scenarios  

---

## Datasets used

This project uses **publicly available datasets**:

- **Six UK-based Low Voltage Distribution Network Models**  
  https://doi.org/10.25405/data.ncl.27175317  

- **D-Suite Modelling: Customer Load Profile Allocations**  
  https://doi.org/10.25405/data.ncl.27175359  

These datasets are used without modification and remain the property of their original authors.

All results and interpretations in this repository are the responsibility of the author and do not represent the views of the dataset creators.

---

## Tools and methods

- Julia  
- PowerModelsDistribution.jl  
- JuMP + Ipopt  
- Unbalanced AC OPF (IVR / explicit neutral formulations)  

The focus is on **planning-level analysis**, not real-time control.

---

## Repository structure (high level)

pmd_pv_mitigation_cba/
├─ Project.toml
├─ Manifest.toml                 # optional to commit; if you want strict reproducibility, commit it
├─ README.md
├─ .gitignore
├─ LICENSE
│
├─ data/
│  ├─ raw/                       # read-only: the two Newcastle datasets unpacked
│  │  ├─ six_uk_lv_models/
│  │  └─ dsuite_load_allocations/
│  ├─ external/                  # anything else (weather, PV shapes) later
│  └─ processed/                 # cached cleaned versions for faster runs
│
├─ configs/
│  ├─ config.toml                # global defaults: voltage limits, VUF limits, solver, etc.
│  ├─ feeders.toml               # which feeders to run, paths to Master.dss, sbase, etc.
│  ├─ pv_profiles.toml           # PV shape assumptions (or mapping to dataset if available)
│  ├─ costs.toml                 # STATCOM/SOP/Battery cost parameters for CBA
│  └─ scenarios.csv              # scenario registry: feeder_id, scenario_id, pv_level, etc.
│
├─ src/                          # reusable Julia code (your "library")
│  ├─ IO.jl                      # parse_file wrappers, path helpers, saving/loading tables
│  ├─ Network.jl                 # transform_data_model, voltage limits, kron flags, etc.
│  ├─ Metrics.jl                 # voltage metrics, margin metrics, VUF, violations, region detection
│  ├─ Plotting.jl                # all plotting functions (voltage profile, histograms, congestion plots)
│  ├─ PV.jl                      # PV insertion functions (single-phase, 3-phase, placement rules)
│  ├─ STATCOM.jl                 # STATCOM insertion, q-limits, extracting q dispatch
│  ├─ SOP.jl                     # SOP model insertion (later)
│  ├─ Battery.jl                 # battery model insertion (later)
│  ├─ OPF.jl                     # solve wrappers, consistent solver settings, status parsing
│  └─ Pipeline.jl                # stage runner utilities, scenario loop, logging
│
├─ scripts/                      # runnable entry points (no notebook needed)
│  ├─ 00_setup_verify.jl         # checks paths, parses one feeder, prints sanity checks
│  ├─ 01_baseline_opf.jl         # baseline solve + baseline plots
│  ├─ 02_pv_sensitivity.jl       # fixed PV size, move PV location rules (source/mid/end/weakest)
│  ├─ 03_pv_severity_sweep.jl    # sweep PV penetration/size at fixed rule
│  ├─ 04_congestion_identify.jl  # produces scenario_summary + per-scenario bus metrics + plots
│  ├─ 05_statcom_mitigation.jl   # reads congestion outputs, runs STATCOM candidate tests
│  ├─ 06_sop_mitigation.jl       # reads congestion outputs, runs SOP (later)
│  ├─ 07_battery_mitigation.jl   # reads congestion outputs, runs battery (later)
│  ├─ 08_cba.jl                  # reads all mitigation summaries and computes CBA tables/plots
│  └─ 99_make_report_figs.jl     # optional: curated figures for report
│
├─ results/
│  ├─ baseline/
│  │  ├─ figures/
│  │  └─ tables/
│  ├─ pv_sensitivity/
│  │  ├─ figures/
│  │  └─ tables/
│  ├─ pv_severity/
│  │  ├─ figures/
│  │  └─ tables/
│  ├─ congestion/
│  │  ├─ figures/
│  │  └─ tables/
│  ├─ statcom/
│  │  ├─ figures/
│  │  └─ tables/
│  ├─ sop/
│  │  ├─ figures/
│  │  └─ tables/
│  ├─ battery/
│  │  ├─ figures/
│  │  └─ tables/
│  └─ cba/
│     ├─ figures/
│     └─ tables/
│
├─ docs/
│  ├─ methodology.md             # explain assumptions, what each stage does, in student language
│  ├─ dataset_notes.md           # how raw dataset maps to feeders & load allocations
│  └─ figures_index.md           # where each figure is produced
│
└─ notebooks/                    # optional explorations (not the main pipeline)
   ├─ scratch_voltage_debug.ipynb
   └─ scratch_pv_profiles.ipynb
 

Each analysis stage writes structured CSV outputs so downstream stages can run without manual intervention.

---

## Disclaimer

This repository is part of an academic research project.

- The code is experimental and under active development  
- Results should not be used for operational or regulatory decisions  
- Simplifying assumptions are made for tractability  

Any errors are the responsibility of the author.

---

## Credits and acknowledgements

- Network and load datasets provided by Newcastle University and collaborators  
- Modelling framework based on PowerModelsDistribution.jl  

This project is informed by discussions and guidance from researchers in the power systems community.

---

## Author

Rhea  
(Student researcher, learning power systems through modelling)
