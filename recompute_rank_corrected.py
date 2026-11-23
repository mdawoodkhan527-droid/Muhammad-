#!/usr/bin/env python3
"""Recompute the paired LHS rank screen after correcting carbon closure."""
import csv, json, os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from extended_identifiability import jacobian_at, chemistry_predict, transform, PARAMS, H, OBS_SCALE, LAWS, N_LHS, LHS_SEED

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"

def work(task):
    law, k, z = task
    s = jacobian_at(np.asarray(z), law)
    singular = np.linalg.svd(s, compute_uv=False)
    nr = int(np.sum(singular / singular[0] >= 1e-8)) if singular[0] else 0
    pr = int(np.sum(singular >= 1.0))
    p = transform(np.asarray(z))
    max_d = p["ln_D"] * np.exp(H * np.log(3.0))
    dt_used = 0.1 if max_d > 2.0e-9 else 0.2
    pred = chemistry_predict(np.asarray(z), law, dt_used)[3]
    return {"sample": k, "seed": LHS_SEED, "permeability_law": law,
            "numerical_rank_rel_1e-8": nr, "practical_rank_sigma_ge_1": pr,
            "dt_s": dt_used,
            "singular_values": ";".join(f"{v:.9g}" for v in singular),
            "predicted_delta13C_calcite_permil": pred,
            "f_org": p["f_org"], "epsilon_SRB_permil": p["epsilon_SRB_permil"],
            "epsilon_calcite_HCO3_permil": p["epsilon_calcite_HCO3_permil"],
            "delta_background_permil": p["delta_background_permil"],
            "delta_org_permil": p["delta_org_permil"], "fit_performed": False}

def main():
    u=qmc.LatinHypercube(d=len(PARAMS), seed=LHS_SEED).random(N_LHS)
    design=2.0*u-1.0
    tasks=[(law,k+1,z.tolist()) for law in LAWS for k,z in enumerate(design)]
    workers=min(8, os.cpu_count() or 1)
    print(f"Recomputing {len(tasks)} paired-design Jacobians with {workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows=list(pool.map(work,tasks,chunksize=1))
    with (OUT/"corrected_lhs_rank_screen.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary={}
    for law in LAWS:
        subset=[r for r in rows if r["permeability_law"]==law]
        summary[law]={}
        for key in ("numerical_rank_rel_1e-8","practical_rank_sigma_ge_1"):
            counts={str(int(v)):sum(int(r[key])==int(v) for r in subset) for v in sorted(set(int(r[key]) for r in subset))}
            summary[law][key]=counts
    (OUT/"corrected_lhs_rank_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    # Carbon output is identical under both permeability laws; save one paired design.
    c=[r for r in rows if r["permeability_law"]==LAWS[0]]
    vals=np.array([float(r["predicted_delta13C_calcite_permil"]) for r in c])
    stats={"n":len(vals),"minimum_permil":float(vals.min()),"median_permil":float(np.median(vals)),"maximum_permil":float(vals.max()),
           "within_published_mean_plusminus_2p1_n":int(np.sum(np.abs(vals+13.7)<=2.1)),"published_mean_permil":-13.7,"published_spread_permil":2.1,
           "operator":"delta_mix = f_org*(delta_org + epsilon_SRB) + (1-f_org)*delta_background; delta_calcite = delta_mix + epsilon_calcite_HCO3"}
    (OUT/"corrected_carbon_closure_summary.json").write_text(json.dumps(stats,indent=2)+"\n")
    print(json.dumps({"rank":summary,"carbon":stats},indent=2),flush=True)

if __name__=="__main__": main()
