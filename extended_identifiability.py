#!/usr/bin/env python3
"""Model-form and LHS rank screens; no fitting or fossil-specific reconstruction."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from coupled_1d_model import KINETIC_BASE, run_case
from identifiability_screen import PARAMS, H, OBS_SCALE, transform
from observation_operator import carbon_isotope_operator, mineral_pattern_operator

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
LAWS = ("kozeny_carman", "power_law")
NX, DT, FINAL = 20, 0.2, 60.0
P_H_BONE = 0.80
P_H_MATRIX = 0.20
ACID_PH = 6.5
ACID_WIDTH = 0.25
BA_SOURCE_RATE = 7.5e-4  # synthetic concentration units s-1 for 60 s
N_LHS = 32
LHS_SEED = 20261111


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def model_form_screen() -> list[dict]:
    cases = [
        ("original_benchmark_bone_Ba_enrichment", False, 0.0, 10.0, 1.0),
        # First four isolate each proposed switch under the benchmark's
        # background Ba concentration. The last two test a source-only
        # endmember in which matrix/inlet Ba is very low.
        ("control_uniform_Ba", False, 0.0, 1.0, 1.0),
        ("a_acidity_switch", True, 0.0, 1.0, 1.0),
        ("b_local_source_at_benchmark_background", False, BA_SOURCE_RATE, 1.0, 1.0),
        ("c_acidity_plus_source_at_benchmark_background", True, BA_SOURCE_RATE, 1.0, 1.0),
        ("b_local_source_low_background_Ba", False, 1.0e-4, 1.0, 1.0e-3),
        ("c_acidity_plus_source_low_background_Ba", True, 1.0e-4, 1.0, 1.0e-3),
    ]
    rows = []
    for law in LAWS:
        for name, acid, ba_rate, initial_bone_factor, ba_scale in cases:
            r = run_case(scenario="S3", permeability_law=law, feedback=False,
                         nx=NX, dt=DT, final_time=FINAL,
                         bone_ba_initial_factor=initial_bone_factor,
                         bone_ba_source_rate=ba_rate,
                         concentration_scales={"Ba": ba_scale},
                         acid_calcite_switch=acid,
                         acid_pH_threshold=ACID_PH,
                         acid_transition_width=ACID_WIDTH)
            op = mineral_pattern_operator(r["node_profile"])
            bone_cells = [p for p in r["node_profile"] if float(p["x_m"]) < 0.25e-3]
            matrix_cells = [p for p in r["node_profile"] if float(p["x_m"]) >= 0.25e-3]
            calcite_bone = sum(float(p["calcite_mol_m3_bulk"]) * float(p["dx_m"]) for p in bone_cells)
            calcite_matrix = sum(float(p["calcite_mol_m3_bulk"]) * float(p["dx_m"]) for p in matrix_cells)
            barite_bone = sum(float(p["barite_mol_m3_bulk"]) * float(p["dx_m"]) for p in bone_cells)
            barite_matrix = sum(float(p["barite_mol_m3_bulk"]) * float(p["dx_m"]) for p in matrix_cells)
            qualifies = (op["p_barite_detect_bone_proxy"] >= P_H_BONE
                         and op["p_barite_detect_matrix_proxy"] <= P_H_MATRIX)
            rows.append({
                "case": name, "permeability_law": law, "acid_switch": acid,
                "acid_pH_threshold": ACID_PH if acid else "off",
                "acid_transition_width_pH": ACID_WIDTH if acid else "off",
                "initial_bone_Ba_factor": initial_bone_factor,
                "global_Ba_concentration_scale": ba_scale,
                "bone_local_Ba_source_rate_synthetic_units_s": ba_rate,
                "p_barite_bone_proxy": op["p_barite_detect_bone_proxy"],
                "p_barite_matrix_proxy": op["p_barite_detect_matrix_proxy"],
                "p_calcite_matrix_proxy": op["p_calcite_detect_matrix_proxy"],
                "barite_localisation_criterion": f"p_bone >= {P_H_BONE:.2f} and p_matrix <= {P_H_MATRIX:.2f}",
                "criterion_met": qualifies,
                "barite_bone_inventory_mol_m2": barite_bone,
                "barite_matrix_inventory_mol_m2": barite_matrix,
                "calcite_bone_inventory_mol_m2": calcite_bone,
                "calcite_matrix_inventory_mol_m2": calcite_matrix,
                "interpretation": "synthetic model-form scenario only; not a fossil-specific result",
            })
    write_csv(OUT / "barite_model_form_screen.csv", rows)
    return rows


def source_rate_sweep() -> list[dict]:
    """Explore a source-only, low-background-Ba endmember without fitting."""
    rates = (1.0e-5, 3.0e-5, 1.0e-4, 2.0e-4, 4.0e-4, 6.0e-4)
    rows = []
    for law in LAWS:
        for acid in (False, True):
            for rate in rates:
                r = run_case(scenario="S3", permeability_law=law, feedback=False,
                             nx=NX, dt=DT, final_time=FINAL,
                             bone_ba_initial_factor=1.0,
                             bone_ba_source_rate=rate,
                             concentration_scales={"Ba": 1.0e-3},
                             acid_calcite_switch=acid,
                             acid_pH_threshold=ACID_PH,
                             acid_transition_width=ACID_WIDTH,
                             compute_mechanics=False)
                op = mineral_pattern_operator(r["node_profile"])
                rows.append({"permeability_law": law, "acid_calcite_switch": acid,
                             "global_Ba_background_scale": 1.0e-3,
                             "bone_Ba_source_rate_synthetic_units_s": rate,
                             "p_barite_bone_proxy": op["p_barite_detect_bone_proxy"],
                             "p_barite_matrix_proxy": op["p_barite_detect_matrix_proxy"],
                             "barite_localisation_criterion": f"p_bone >= {P_H_BONE:.2f} and p_matrix <= {P_H_MATRIX:.2f}",
                             "criterion_met": op["p_barite_detect_bone_proxy"] >= P_H_BONE and op["p_barite_detect_matrix_proxy"] <= P_H_MATRIX,
                             "total_barite_mol_m2": r["summary"]["total_barite_mol_m2"],
                             "interpretation": "sensitivity sweep in an assumed low-background-Ba endmember; not calibrated"})
    write_csv(OUT / "bone_Ba_source_rate_sweep.csv", rows)
    return rows


def chemistry_predict(z: np.ndarray, law: str, dt: float) -> np.ndarray:
    p = transform(z)
    rates = {r["key"]: p[r["name"]] for r in PARAMS if r["kind"] == "kinetic"}
    scales = {r["key"]: p[r["name"]] for r in PARAMS if r["kind"] == "concentration"}
    model = run_case(scenario="S3", permeability_law=law, feedback=False,
                     nx=NX, dt=dt, final_time=FINAL, kinetics=rates,
                     concentration_scales=scales, diffusivity=p["ln_D"],
                     permeability_reference=p["ln_k0"], compute_mechanics=False)
    op = mineral_pattern_operator(model["node_profile"])
    carbon = carbon_isotope_operator(
        f_org=p["f_org"], delta_org=p["delta_org_permil"],
        delta_background_calcite=p["delta_background_permil"],
        epsilon_srb=p["epsilon_SRB_permil"],
        epsilon_calcite_hco3=p["epsilon_calcite_HCO3_permil"])
    return np.array([op["p_barite_detect_bone_proxy"],
                     op["p_barite_detect_matrix_proxy"],
                     op["p_calcite_detect_matrix_proxy"], carbon], dtype=float)


def jacobian_at(z: np.ndarray, law: str) -> np.ndarray:
    s = np.zeros((4, len(PARAMS)), dtype=float)
    # Use a single timestep for every column of this sample's Jacobian, selected
    # from the largest D reached by its central-difference perturbations.
    center = transform(z)
    max_d = center["ln_D"] * math.exp(H * math.log(3.0))
    dt = 0.1 if max_d > 2.0e-9 else DT
    for j, p in enumerate(PARAMS):
        zp, zm = z.copy(), z.copy()
        zp[j] += H
        zm[j] -= H
        if p["kind"] == "carbon":
            plus, minus = transform(zp), transform(zm)
            op = lambda a: carbon_isotope_operator(
                f_org=a["f_org"], delta_org=a["delta_org_permil"],
                delta_background_calcite=a["delta_background_permil"],
                epsilon_srb=a["epsilon_SRB_permil"],
                epsilon_calcite_hco3=a["epsilon_calcite_HCO3_permil"])
            s[3, j] = (op(plus) - op(minus)) / (2 * H * OBS_SCALE[3])
        else:
            yp = chemistry_predict(zp, law, dt)
            ym = chemistry_predict(zm, law, dt)
            s[:, j] = (yp - ym) / (2 * H * OBS_SCALE)
    return s


def lhs_rank_screen() -> tuple[list[dict], list[dict], list[dict]]:
    # The design distribution spans +/- one previously declared parameter
    # screening unit. It is a Latin-hypercube coverage design, not a posterior.
    u = qmc.LatinHypercube(d=len(PARAMS), seed=LHS_SEED).random(N_LHS)
    design = 2.0 * u - 1.0
    all_rows, rank_rows, delta_rows = [], [], []
    for law in LAWS:
        for k, z in enumerate(design, start=1):
            s = jacobian_at(z, law)
            singular = np.linalg.svd(s, compute_uv=False)
            numerical_rank = int(np.sum(singular / singular[0] >= 1e-8)) if singular[0] else 0
            practical_rank = int(np.sum(singular >= 1.0))
            p = transform(z)
            center = transform(z)
            max_d = center["ln_D"] * math.exp(H * math.log(3.0))
            dt_used = 0.1 if max_d > 2.0e-9 else DT
            pred = chemistry_predict(z, law, dt_used)[3]
            all_rows.append({
                "sample": k, "seed": LHS_SEED, "permeability_law": law,
                "numerical_rank_rel_1e-8": numerical_rank,
                "practical_rank_sigma_ge_1": practical_rank,
                "dt_s": dt_used,
                "singular_values": ";".join(f"{v:.9g}" for v in singular),
                "f_org": p["f_org"], "epsilon_SRB_permil": p["epsilon_SRB_permil"],
                "epsilon_calcite_HCO3_permil": p["epsilon_calcite_HCO3_permil"],
                "delta_background_permil": p["delta_background_permil"],
                "delta_org_permil": p["delta_org_permil"],
                "predicted_delta13C_calcite_permil": pred,
                "fit_performed": False,
            })
            if law == LAWS[0]:
                delta_rows.append({"sample": k, "seed": LHS_SEED,
                                   "predicted_delta13C_calcite_permil": pred,
                                   "f_org": p["f_org"],
                                   "epsilon_SRB_permil": p["epsilon_SRB_permil"],
                                   "epsilon_calcite_HCO3_permil": p["epsilon_calcite_HCO3_permil"],
                                   "delta_background_permil": p["delta_background_permil"],
                                   "delta_org_permil": p["delta_org_permil"],
                                   "within_published_mean_plusminus_2p1_spread": abs(pred + 13.7) <= 2.1})
            print(f"LHS {law}: {k}/{N_LHS}", flush=True)
    write_csv(OUT / "prior_envelope_lhs_sensitivity.csv", all_rows)
    for law in LAWS:
        subset = [r for r in all_rows if r["permeability_law"] == law]
        for key in ("numerical_rank_rel_1e-8", "practical_rank_sigma_ge_1"):
            counts = {}
            for r in subset:
                v = int(r[key]); counts[str(v)] = counts.get(str(v), 0) + 1
            rank_rows.append({"permeability_law": law, "rank_type": key,
                              "rank_counts": json.dumps(counts, sort_keys=True),
                              "sample_count": len(subset), "LHS_seed": LHS_SEED})
    write_csv(OUT / "prior_envelope_rank_distribution.csv", rank_rows)
    write_csv(OUT / "carbon_closure_lhs_sensitivity.csv", delta_rows)
    return all_rows, rank_rows, delta_rows


def main() -> None:
    mf = model_form_screen()
    print("model-form cases complete; criterion met:",
          [(r["permeability_law"], r["case"], r["criterion_met"]) for r in mf], flush=True)
    print("source-rate sweep:", [(r["permeability_law"], r["acid_calcite_switch"], r["bone_Ba_source_rate_synthetic_units_s"], r["criterion_met"]) for r in source_rate_sweep()], flush=True)
    lhs_rows, rank_rows, _ = lhs_rank_screen()
    # Summarise carbon propagation over the LHS design and identify rank stability.
    d = list(csv.DictReader((OUT / "carbon_closure_lhs_sensitivity.csv").open()))
    values = np.array([float(r["predicted_delta13C_calcite_permil"]) for r in d])
    carbon_modes = []
    for r in d:
        f = float(r["f_org"])
        eo = float(r["delta_org_permil"])
        eb = float(r["delta_background_permil"])
        es = float(r["epsilon_SRB_permil"])
        ec = float(r["epsilon_calcite_HCO3_permil"])
        derivatives = np.array([
            0.25 * (eo + es + ec - eb),
            2.0 * f,
            0.2 * f,
            2.0 * (1.0 - f),
            0.5 * f,
        ]) / 2.1
        carbon_modes.append(float(np.linalg.norm(derivatives)))
    within = sum(r["within_published_mean_plusminus_2p1_spread"].lower() == "true" for r in d)
    s = {"design": "32-point LHS, seed 20261111; coordinates uniformly span +/- one frozen screening unit",
         "fit_performed": False, "delta13C_target_permil": -13.7,
         "scale_used_permil": 2.1,
         "predicted_delta13C_min_permil": float(values.min()),
         "predicted_delta13C_max_permil": float(values.max()),
         "predicted_delta13C_median_permil": float(np.median(values)),
         "samples_within_target_plusminus_scale": int(within),
         "fractionation_envelope_epsilon_SRB_permil": [-2.0, 2.0],
         "fractionation_envelope_epsilon_calcite_HCO3_permil": [0.8, 1.2],
         "carbon_closure_local_rank": 1,
         "carbon_direction_singular_value_min": float(min(carbon_modes)),
         "carbon_direction_singular_value_max": float(max(carbon_modes)),
         "carbon_direction_practically_distinguishable_all_samples": bool(min(carbon_modes) >= 1.0),
         "joint_practical_rank_counts": {
             f"{law}:rank_{rank}": sum(1 for r in lhs_rows
                                        if r["permeability_law"] == law
                                        and int(r["practical_rank_sigma_ge_1"]) == rank)
             for law in LAWS for rank in sorted({int(r["practical_rank_sigma_ge_1"])
                                                  for r in lhs_rows if r["permeability_law"] == law})
         },
         "sample_count": len(values),
         "interpretation": "screening envelope propagation; not a fitted or probabilistic posterior interval"}
    (OUT / "carbon_closure_lhs_summary.json").write_text(json.dumps(s, indent=2) + "\n")
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
