#!/usr/bin/env python3
"""Prefit local SVD screen for the frozen mineral/13C observation operators."""
from __future__ import annotations
import csv, json, math
from pathlib import Path
import numpy as np
from coupled_1d_model import run_case, KINETIC_BASE
from observation_operator import mineral_pattern_operator, carbon_isotope_operator

ROOT=Path(__file__).resolve().parent; OUT=ROOT/"results"; OUT.mkdir(exist_ok=True)
# Coordinates are normalized so one unit is the declared screening scale in protocol.
PARAMS=[
 {"name":"ln_k_SRB","kind":"kinetic","key":"k_srb","base":KINETIC_BASE["k_srb"],"scale_ln":math.log(3)},
 {"name":"ln_k_SOX","kind":"kinetic","key":"k_sox","base":KINETIC_BASE["k_sox"],"scale_ln":math.log(3)},
 {"name":"ln_k_SDB","kind":"kinetic","key":"k_sdb","base":KINETIC_BASE["k_sdb"],"scale_ln":math.log(3)},
 {"name":"ln_k_calcite","kind":"kinetic","key":"k_calcite","base":KINETIC_BASE["k_calcite"],"scale_ln":math.log(3)},
 {"name":"ln_k_barite","kind":"kinetic","key":"k_barite","base":KINETIC_BASE["k_barite"],"scale_ln":math.log(3)},
 {"name":"ln_D","kind":"diffusivity","base":1.0e-9,"scale_ln":math.log(3)},
 {"name":"ln_k0","kind":"permeability","base":1.0e-18,"scale_ln":math.log(10)},
 {"name":"ln_Ba0","kind":"concentration","key":"Ba","base":1.0,"scale_ln":math.log(2)},
 {"name":"ln_SO40","kind":"concentration","key":"SO4","base":1.0,"scale_ln":math.log(2)},
 {"name":"ln_OM0","kind":"concentration","key":"OM","base":1.0,"scale_ln":math.log(2)},
 {"name":"f_org","kind":"carbon","base":0.5,"scale_abs":0.25},
 {"name":"epsilon_SRB_permil","kind":"carbon","base":0.0,"scale_abs":2.0},
 {"name":"epsilon_calcite_HCO3_permil","kind":"carbon","base":1.0,"scale_abs":0.2},
 {"name":"delta_background_permil","kind":"carbon","base":-3.0,"scale_abs":2.0},
 {"name":"delta_org_permil","kind":"carbon","base":-31.7,"scale_abs":0.5},
]
H=0.05; OBS_SCALE=np.array([1.,1.,1.,2.1])
LAW_LABEL={"kozeny_carman":"Kozeny-Carman","power_law":"Power law"}
MEASUREMENTS={
 "ln_k_SRB":"Time-resolved sulphate and dissolved inorganic carbon profiles under constrained OM supply",
 "ln_k_SOX":"Co-located sulfate, sulfide and S0 inventories through time, with an electron-acceptor budget",
 "ln_k_SDB":"Paired time-resolved S0, sulfate and sulfide abundance and isotope measurements",
 "ln_k_calcite":"Quantitative calcite volume with co-located pH, alkalinity, Ca and DIC",
 "ln_k_barite":"Quantitative barite volume plus porewater Ba and sulfate at registered bone/matrix locations",
 "ln_D":"Spatial porewater tracer profiles with a known length scale and measured diffusion coefficient",
 "ln_k0":"Hydraulic conductivity or permeability measurement on the relevant matrix",
 "ln_Ba0":"Porewater Ba concentration or a source-resolved Ba inventory",
 "ln_SO40":"Porewater sulfate concentration or a source-resolved sulfate inventory",
 "ln_OM0":"Labile organic-carbon inventory and independently measured degradation rate",
 "f_org":"Paired δ13C of DIC, carbonate, organic substrate and source endmembers from the same microenvironment",
 "epsilon_SRB_permil":"Culture or natural-system δ13C of consumed substrate and released DIC for the relevant SRB pathway",
 "epsilon_calcite_HCO3_permil":"Paired δ13C measurements of bicarbonate and precipitated calcite at constrained temperature",
 "delta_background_permil":"Co-located background DIC δ13C measurement; shale carbonate alone is not a fluid measurement",
 "delta_org_permil":"Direct δ13C measurement of the reactive organic substrate, not only bulk rib OM",
}

def transform(z):
    p={}
    for i,r in enumerate(PARAMS):
        if "scale_ln" in r: p[r["name"]]=r["base"]*math.exp(float(z[i])*r["scale_ln"])
        else: p[r["name"]]=r["base"]+float(z[i])*r["scale_abs"]
    if not 0 < p["f_org"] < 1: raise ValueError("f_org perturbation outside (0,1)")
    return p

def predict(z, law, nx, dt):
    p=transform(z)
    rates={r["key"]:p[r["name"]] for r in PARAMS if r["kind"]=="kinetic"}
    scales={r["key"]:p[r["name"]] for r in PARAMS if r["kind"]=="concentration"}
    model=run_case(scenario="S3",permeability_law=law,feedback=False,nx=nx,dt=dt,final_time=60.0,
                   kinetics=rates,concentration_scales=scales,diffusivity=p["ln_D"],
                   permeability_reference=p["ln_k0"],compute_mechanics=False)
    op=mineral_pattern_operator(model["node_profile"])
    carbon=carbon_isotope_operator(f_org=p["f_org"],delta_org=p["delta_org_permil"],
                    delta_background_calcite=p["delta_background_permil"],
                    epsilon_srb=p["epsilon_SRB_permil"],epsilon_calcite_hco3=p["epsilon_calcite_HCO3_permil"])
    return np.array([op["p_barite_detect_bone_proxy"],op["p_barite_detect_matrix_proxy"],
                     op["p_calcite_detect_matrix_proxy"],carbon],dtype=float)

def finite_difference_matrix(law,nx,dt):
    n=len(PARAMS); S=np.zeros((4,n)); zero=np.zeros(n)
    baseline=transform(zero)
    base_model=run_case(scenario="S3",permeability_law=law,feedback=False,nx=nx,dt=dt,final_time=60.0,
              kinetics={r["key"]:baseline[r["name"]] for r in PARAMS if r["kind"]=="kinetic"},
              concentration_scales={r["key"]:baseline[r["name"]] for r in PARAMS if r["kind"]=="concentration"},
              diffusivity=baseline["ln_D"],permeability_reference=baseline["ln_k0"])
    base_mineral=mineral_pattern_operator(base_model["node_profile"])
    for j in range(n):
        if PARAMS[j]["kind"]=="carbon":
            zp=zero.copy(); zm=zero.copy(); zp[j]=H; zm[j]=-H
            pp=transform(zp); pm=transform(zm)
            op=np.array([0.,0.,0.,carbon_isotope_operator(f_org=pp["f_org"],delta_org=pp["delta_org_permil"],delta_background_calcite=pp["delta_background_permil"],epsilon_srb=pp["epsilon_SRB_permil"],epsilon_calcite_hco3=pp["epsilon_calcite_HCO3_permil"])])
            om=np.array([0.,0.,0.,carbon_isotope_operator(f_org=pm["f_org"],delta_org=pm["delta_org_permil"],delta_background_calcite=pm["delta_background_permil"],epsilon_srb=pm["epsilon_SRB_permil"],epsilon_calcite_hco3=pm["epsilon_calcite_HCO3_permil"])])
            S[:,j]=(op-om)/(2*H*OBS_SCALE)
            continue
        zp=zero.copy(); zm=zero.copy(); zp[j]=H; zm[j]=-H
        pplus=transform(zp); pminus=transform(zm)
        def chemistry_profile(p):
            model=run_case(scenario="S3",permeability_law=law,feedback=False,nx=nx,dt=dt,final_time=60.0,
                   kinetics={r["key"]:p[r["name"]] for r in PARAMS if r["kind"]=="kinetic"},
                   concentration_scales={r["key"]:p[r["name"]] for r in PARAMS if r["kind"]=="concentration"},
                   diffusivity=p["ln_D"],permeability_reference=p["ln_k0"],compute_mechanics=False)
            op=mineral_pattern_operator(model["node_profile"])
            return np.array([op["p_barite_detect_bone_proxy"],op["p_barite_detect_matrix_proxy"],op["p_calcite_detect_matrix_proxy"],0.])
        plus=chemistry_profile(pplus); minus=chemistry_profile(pminus)
        S[:,j]=(plus-minus)/(2*H*OBS_SCALE)
    return S

def top_terms(v, cutoff=0.18):
    order=np.argsort(np.abs(v))[::-1]
    return [{"parameter":PARAMS[i]["name"],"coefficient":float(v[i])} for i in order if abs(v[i])>=cutoff]

def analyze(law,nx,dt):
    S=finite_difference_matrix(law,nx,dt)
    U,s,Vt=np.linalg.svd(S,full_matrices=True)
    rank=int(np.sum((s/s[0])>=1e-8)) if len(s) and s[0]>0 else 0
    practical=[float(v) for v in s if v>=1.0]
    directions=[]
    for i,val in enumerate(s):
        directions.append({"mode":i+1,"singular_value":float(val),"relative_singular_value":float(val/s[0]) if s[0] else 0,
                           "distinguishable":bool(val>=1.0),"dominant_parameter_combination":top_terms(Vt[i])})
    null=[]
    for i in range(rank,len(PARAMS)):
        v=Vt[i]
        null.append({"null_mode":i-rank+1,"singular_value":float(s[i]) if i<len(s) else 0.0,
                     "dominant_parameter_combination":top_terms(v),
                     "measurement_to_break_direction":"; ".join(sorted({MEASUREMENTS[PARAMS[j]["name"]] for j in np.argsort(np.abs(v))[::-1][:4]}))})
    return {"permeability_law":LAW_LABEL[law],"nx":nx,"dt_s":dt,"observation_count":4,
            "parameter_count":len(PARAMS),"matrix_rank":rank,"structural_rank_ceiling":4,
            "null_space_dimension":len(PARAMS)-rank,"distinguishable_direction_count":len(practical),
            "singular_values":s.tolist(),"directions":directions,"null_directions":null,
            "S_whitened":S.tolist()}

def main():
    # Store mask-boundary envelope at each target resolution before interpreting SVD.
    envelope=[]
    for nx,dt in [(20,0.2),(40,0.05)]:
        for law in ("kozeny_carman","power_law"):
            base=run_case(scenario="S3",permeability_law=law,feedback=False,nx=nx,dt=dt,final_time=60.0)
            vals={"p_barite_detect_bone_proxy":[],"p_barite_detect_matrix_proxy":[],"p_calcite_detect_matrix_proxy":[]}
            for b in np.linspace(.20,.30,21):
                q=mineral_pattern_operator(base["node_profile"],boundary_fraction=(float(b),float(b)),boundary_nodes=1)
                for k in vals: vals[k].append(q[k])
            for k,v in vals.items():
                envelope.append({"permeability_law":LAW_LABEL[law],"nx":nx,"dt_s":dt,"operator":k,
                                "boundary_low":.20,"boundary_high":.30,"minimum":min(v),"maximum":max(v),
                                "mean_over_boundary":float(np.mean(v))})
    with (OUT/"mask_boundary_envelope.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=envelope[0].keys());w.writeheader();w.writerows(envelope)
    results=[analyze(law,nx,dt) for law in ("kozeny_carman","power_law") for nx,dt in [(20,.2),(40,.05)]]
    (OUT/"identifiability_svd.json").write_text(json.dumps({"protocol":"identifiability_protocol.md","fit_performed":False,"parameter_names":[x["name"] for x in PARAMS],"results":results},indent=2)+"\n",encoding="utf-8")
    flat=[]
    for r in results:
        flat.append({"permeability_law":r["permeability_law"],"nx":r["nx"],"matrix_rank":r["matrix_rank"],
                     "null_space_dimension":r["null_space_dimension"],"distinguishable_direction_count":r["distinguishable_direction_count"],
                     "singular_values":";".join(f"{v:.6g}" for v in r["singular_values"])})
    with (OUT/"identifiability_rank_summary.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=flat[0].keys());w.writeheader();w.writerows(flat)
    print(json.dumps(flat,indent=2))

if __name__=="__main__": main()
