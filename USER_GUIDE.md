# User guide

## Requirements

- Python 3.10 or later
- NumPy, SciPy, and PyYAML, installed from `requirements.txt`
- A terminal or command prompt

## Install

From the repository root, create and activate an environment. On Windows PowerShell:

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Reproduce the model and analyses

Run commands from the repository root in this order:

```bash
python verified_1d_reference_model.py --out results_reference
python coupled_1d_model.py
python run_refinement_check.py
python observation_operator.py
python identifiability_screen.py
python postprocess_identifiability.py
python extended_identifiability.py
```

The scripts use the frozen settings in `tolerances.yaml`. Model and analysis outputs are written to `results/` unless a script reports a different output directory. Existing `results/` files are the outputs used in the manuscript and should be copied elsewhere before rerunning if they must be preserved.

Figure scripts are in `figures/`. Run a figure script from the repository root; inspect its header for any output-path options before execution. PNG files included beside the scripts are the plotted outputs.

## Main outputs

- `results/verification_summary.json` and `results/reaction_balance.csv`: limiting-case, reaction-balance, and conservation summaries.
- `results/coupled_timeseries.csv`, `results/coupled_node_profiles.csv`, and `results/stage1_stage2_contrasts.csv`: model trajectories and coupling comparisons.
- `results/observation_operator_catalog.csv` and `results/observation_fit_eligibility.csv`: observation mapping and inverse-fit eligibility.
- `results/corrected_lhs_rank_summary.json`, `results/identifiability_svd_refined.json`, and related CSV files: rank and sensitivity summaries.
- `results/corrected_carbon_closure_summary.json` and `results/carbon_closure_lhs_sensitivity_corrected.csv`: corrected isotope closure outputs.
- `results/barite_model_form_screen.csv` and `results/bone_Ba_source_rate_sweep.csv`: barite model-form screening.

## Interpretation limits

The calculations test numerical consistency and sensitivity of a synthetic benchmark. They do not demonstrate that the assumed rates, chemistry, permeability, pressure, stiffness, or reaction pathways applied to the fossil. No fit or calibration was performed.
