#!/usr/bin/env python3
"""Compare coupling signal with the predefined space/time refinement envelope."""
import csv
import json
from pathlib import Path
from coupled_1d_model import TOL, run_case

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
METRICS = ["mean_porosity", "total_calcite_mol_m2", "total_barite_mol_m2", "mean_pH", "mean_strain"]
rows = []
for scenario in ("S1", "S2", "S3"):
    for law in ("kozeny_carman", "power_law"):
        cases = {}
        for label, nx, dt in (("coarse", 20, 0.2), ("fine", 40, 0.05)):
            for feedback in (False, True):
                result = run_case(scenario=scenario, permeability_law=law, feedback=feedback,
                                  nx=nx, dt=dt, final_time=60.0)
                cases[(label, feedback)] = result["summary"]
        for metric in METRICS:
            a = cases[("fine", False)][metric]
            b = cases[("fine", True)][metric]
            scale = max(abs(a), 1e-12)
            feedback_relative = abs(b - a) / scale
            refinement_relative = max(
                abs(cases[("fine", mode)][metric] - cases[("coarse", mode)][metric]) / scale
                for mode in (False, True)
            )
            threshold = TOL["coupling_decision"]["stage1_vs_stage2_relative_output_difference"]
            rows.append({
                "scenario": scenario, "permeability_law": law, "metric": metric,
                "stage1_fine": a, "stage2_fine": b,
                "stage1_stage2_relative_difference": feedback_relative,
                "coarse_fine_refinement_envelope": refinement_relative,
                "threshold": threshold,
                "above_threshold": feedback_relative > threshold,
                "above_refinement_envelope": feedback_relative > refinement_relative,
                "retain_two_way": feedback_relative > threshold and feedback_relative > refinement_relative,
            })
with (OUT / "coupling_refinement_check.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader(); w.writerows(rows)
count = sum(r["retain_two_way"] for r in rows)
summary = {
    "comparisons": len(rows),
    "feedback_effects_above_1pct": sum(r["above_threshold"] for r in rows),
    "feedback_effects_above_refinement_envelope": sum(r["above_refinement_envelope"] for r in rows),
    "effects_meeting_both_rules": count,
    "maximum_feedback_relative_difference": max(r["stage1_stage2_relative_difference"] for r in rows),
    "maximum_refinement_envelope": max(r["coarse_fine_refinement_envelope"] for r in rows),
    "decision": "retain Stage 1 one-way coupling for this synthetic benchmark" if count == 0 else "retain Stage 2 for qualifying outputs",
    "comparison_design": "20 cells, dt=0.2 s versus 40 cells, dt=0.05 s; final time 60 s; same prescribed inputs",
}
(OUT / "coupling_refinement_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

# Keep the main verification summary synchronized with the completed refinement
# decision rather than leaving its provisional status after a rerun.
verification_path = OUT / "verification_summary.json"
if verification_path.exists():
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    verification["coupling_decision"] = summary["decision"]
    verification_path.write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")

# Replace the provisional note in the initial coarse-grid contrast table with
# the final rule-based decision for each scenario and permeability law.
path = OUT / "stage1_stage2_contrasts.csv"
with path.open(newline="", encoding="utf-8") as f:
    contrasts = list(csv.DictReader(f))
for r in contrasts:
    matched = [x for x in rows if x["scenario"] == r["scenario"]
               and x["permeability_law"] == r["permeability_law"]]
    r["two_way_decision"] = "retain one-way for this benchmark" if not any(x["retain_two_way"] for x in matched) else "retain two-way for qualifying output"
with path.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=contrasts[0].keys())
    w.writeheader(); w.writerows(contrasts)
print(summary)
