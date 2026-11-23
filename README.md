# Fossil-concretion observability benchmark

This repository contains the reduced one-dimensional reactive-transport and poromechanics benchmark, its observation-operator and identifiability analyses, the outputs cited in the accompanying manuscript, and the figure-generation files.

## Scope

The code is a synthetic numerical benchmark. It is not a specimen-specific reconstruction, a calibrated reactive-transport model, or a substitute for PHREEQC speciation. Model inputs and kinetics are assumptions. The manuscript reports no new experimental data.

## Contents

- `coupled_1d_model.py`, `verified_1d_reference_model.py`, and `run_refinement_check.py`: model and verification calculations.
- `observation_operator.py`, `identifiability_screen.py`, `postprocess_identifiability.py`, `extended_identifiability.py`, `recompute_rank_corrected.py`, and `save_svd_results.py`: observation, sensitivity, and rank analyses.
- `tolerances.yaml` and `threshold_justification_preresults.md`: frozen numerical and distinguishability criteria.
- `results/`: CSV and JSON outputs used in the manuscript.
- `figures/`: figure-generation scripts and image outputs.
- `data/Model_Input_Register_and_Source_Trace.xlsx`: parameter provenance register.
- `USER_GUIDE.md`: installation and reproduction instructions.

## Quick start

Use Python 3.10 or later. Install dependencies with `python -m pip install -r requirements.txt`. Follow `USER_GUIDE.md` for the full reproduction sequence.

## Data and interpretation

Published observations are catalogued in `results/published_observation_ledger.csv`; derived values and conflicts are labelled in the ledger. The Conc4-A/B total-sulphur conflict is retained in `results/source_conflict_ledger.csv` and excluded from sulphur mass balance. The LHS design is a sensitivity coverage design, not a probabilistic sample. No parameters are fitted.

## Licence

The intended licence is MIT, subject to confirmation by the code copyright holder(s). Complete the copyright-holder line in `LICENSE` and obtain any needed co-author approval before making the licence effective.
