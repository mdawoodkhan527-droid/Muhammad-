#!/usr/bin/env python3
"""Observation-to-state map for published Jian et al. (2026) constraints.

The catalogue deliberately separates direct measurements, derived quantities,
source conflicts and proposed mechanisms. Only mineral inventory integration
is currently implemented as a forward operator, and it requires an explicit
sample mask. A presence/absence result additionally requires a user-supplied
detection threshold. No threshold or fossil-zone mapping is inferred here.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Sequence

OBSERVATION_CATALOG = [
    {"id":"PETRO_BARITE_CAVITY", "observable":"Barite observed in bone marrow cavities", "status":"direct qualitative petrographic observation", "source":"Jian et al. (2026), mineralogical figures and text", "model_state":"barite_mol_m3_bulk", "operator":"prior-averaged logistic presence probability over an uncertain 1D mask and threshold fraction", "availability":"implemented for synthetic bone-proxy mask only; not registered fossil anatomy", "uniqueness":"non-unique: saturation, Ba supply, sulfate supply, nucleation and preservation can yield the same presence pattern", "inverse_use":"not admitted to likelihood until fossil-to-grid registration and detection limit are supported"},
    {"id":"PETRO_CALCITE_MATRIX", "observable":"Calcite-dominant concretion matrix", "status":"direct qualitative petrographic/mineralogical observation", "source":"Jian et al. (2026), mineralogical figures and text", "model_state":"calcite_mol_m3_bulk", "operator":"prior-averaged logistic presence probability over an uncertain 1D mask and threshold fraction", "availability":"implemented for synthetic matrix-proxy mask only; not registered fossil anatomy", "uniqueness":"non-unique: pH, DIC, Ca, kinetics, nucleation and subsequent alteration can produce similar calcite abundance", "inverse_use":"not admitted to likelihood until fossil-to-grid registration and detection limit are supported"},
    {"id":"PETRO_FAP_REPLACEMENT", "observable":"Fluorapatite replacement of collagen", "status":"direct qualitative mineralogical observation", "source":"Jian et al. (2026), mineralogical figures and text", "model_state":"none; FAP absent from coupled solver", "operator":"none implemented", "availability":"unavailable", "uniqueness":"non-unique: precursor composition, F supply, pH, dissolution-reprecipitation and kinetic inhibition are unresolved", "inverse_use":"exclude from current inverse problem"},
    {"id":"TIC", "observable":"Total inorganic carbon", "status":"direct measurement", "source":"Jian et al. (2026), Table 2; source-trace workbook cell references", "model_state":"calcite carbon inventory", "operator":"requires carbonate phase inventory, sample support, dry bulk mass and unit conversion", "availability":"not identifiable from current output; sample support and bulk-density conversion absent", "uniqueness":"non-unique with respect to carbonate phase mix and sampled volume", "inverse_use":"comparison only after sample support and conversion are defined"},
    {"id":"D13C_CARB", "observable":"Carbonate δ13C", "status":"direct measurement", "source":"Jian et al. (2026), Table 2", "model_state":"none; isotope ratios absent", "operator":"post-processing two-source isotope mixture with explicit endmembers and separate SRB and calcite-HCO3 offsets", "availability":"implemented as a conditional closure; it is not coupled to modeled carbon-source provenance", "uniqueness":"non-unique: organic fraction, source endmembers, SRB fractionation and carbonate-DIC fractionation trade off", "inverse_use":"secondary closure only; no fitting in the current screen"},
    {"id":"D18O_CARB", "observable":"Carbonate δ18O", "status":"direct measurement", "source":"Jian et al. (2026), Table 2", "model_state":"none; oxygen isotopes absent", "operator":"requires water δ18O, temperature, mineral-water fractionation and alteration model", "availability":"unavailable", "uniqueness":"non-unique: temperature and water δ18O are covariant; alteration adds further ambiguity", "inverse_use":"exclude from current inverse problem"},
    {"id":"D13C_ORG", "observable":"Organic matter δ13C", "status":"direct measurement", "source":"Jian et al. (2026), Table 2", "model_state":"OM concentration only; no isotope ratio", "operator":"requires source-specific organic carbon isotope pools and fractionation", "availability":"unavailable", "uniqueness":"non-unique: substrate mixture and microbial fractionation are confounded", "inverse_use":"exclude from current inverse problem"},
    {"id":"TOTAL_S", "observable":"Total sulphur abundance", "status":"direct measurement with two unresolved source conflicts", "source":"Jian et al. (2026), Table 2; Mendeley Table 2 and Rock-Eval workbooks", "model_state":"aqueous SO4, H2S/HS and barite only; pyrite/TRIS and other S solids absent", "operator":"requires all sulphur pools and the exact sample support", "availability":"not identifiable for whole-rock total S", "uniqueness":"non-unique among unmodelled pyrite, sulphate minerals, aqueous sulphur and sampled mass", "inverse_use":"exclude Conc4-A/B from total-S calibration; other rows comparison only after phase closure"},
    {"id":"TRIS_ABUNDANCE", "observable":"TRIS abundance", "status":"direct measurement", "source":"Jian et al. (2026), Table 2", "model_state":"none; pyrite/TRIS phase absent", "operator":"none implemented", "availability":"unavailable", "uniqueness":"non-unique without iron, pyrite precipitation and sulfur pool closure", "inverse_use":"exclude from current inverse problem"},
    {"id":"SULPHATE_ABUNDANCE", "observable":"Sulphate abundance", "status":"derived by total-S minus TRIS difference", "source":"Jian et al. (2026), Table 2 footnotes and mass balance", "model_state":"aqueous SO4 and barite", "operator":"model dissolved plus precipitated sulphate-S over matching sample support", "availability":"partial state overlap only; calculated rock abundance cannot be equated to dissolved SO4", "uniqueness":"non-unique due phase allocation and the Conc4 total-S conflicts", "inverse_use":"comparison-only derived value; never count as independent datum alongside total S and TRIS"},
    {"id":"D34S_TOTAL", "observable":"Total sulphur δ34S", "status":"direct measurement", "source":"Jian et al. (2026), Table 2", "model_state":"none; no isotope-resolved S pools", "operator":"requires isotope-resolved inventories for every S-bearing phase", "availability":"unavailable", "uniqueness":"non-unique with source mixing and phase allocation", "inverse_use":"exclude from current inverse problem"},
    {"id":"D34S_TRIS", "observable":"TRIS δ34S", "status":"direct measurement", "source":"Jian et al. (2026), Table 2", "model_state":"none; TRIS mineral phase and isotope ratio absent", "operator":"requires isotope-resolved TRIS production and fractionation", "availability":"unavailable", "uniqueness":"non-unique between fractionation, substrate supply, oxidation and mixing", "inverse_use":"exclude from current inverse problem"},
    {"id":"D34S_SULPHATE", "observable":"Sulphate δ34S", "status":"approximate calculated value, not direct measurement", "source":"Jian et al. (2026), Table 2 isotope mass balance", "model_state":"none; isotope ratios absent", "operator":"requires isotope-resolved sulfate and complete paired sulfur balance", "availability":"unavailable", "uniqueness":"non-unique; derived value inherits total-S/TRIS uncertainty and source conflicts", "inverse_use":"exclude from likelihood; retain only as secondary comparison with provenance flag"},
    {"id":"ROCKEVAL_S1_S2_TOC_TMAX", "observable":"Rock-Eval S1, S2, TOC and Tmax", "status":"direct measurement", "source":"Jian et al. (2026), Table 1 and published workbook", "model_state":"OM concentration only", "operator":"requires thermal maturation and Rock-Eval response model", "availability":"unavailable", "uniqueness":"non-unique: present-day OM concentration does not determine pyrolysis response or Tmax", "inverse_use":"do not fit current model"},
    {"id":"HI_OI", "observable":"Hydrogen index and oxygen index", "status":"calculated Rock-Eval indices", "source":"Jian et al. (2026), Table 1; derived from Rock-Eval measurements", "model_state":"none", "operator":"calculated from measured Rock-Eval yields and TOC, not from current state", "availability":"not a model observation operator", "uniqueness":"not independent of S1/S2/TOC inputs", "inverse_use":"never count as additional independent observations"},
    {"id":"N_AND_D15N", "observable":"Total nitrogen and δ15N", "status":"direct measurement", "source":"Jian et al. (2026), Table 2", "model_state":"none; nitrogen absent", "operator":"none implemented", "availability":"unavailable", "uniqueness":"not applicable to current model state", "inverse_use":"exclude from current inverse problem"},
    {"id":"RIB_D34S_OFFSET", "observable":"Rib sulphate-TRIS δ34S offset", "status":"derived from measured TRIS δ34S and approximate calculated sulphate δ34S", "source":"Jian et al. (2026), Table 2; pairing must remain sample-specific", "model_state":"none; isotope-resolved sulfate/TRIS absent", "operator":"paired isotope mass balance, requiring covariance and source-specific fractions", "availability":"unavailable", "uniqueness":"non-unique; do not substitute quadrature of unpaired standard deviations for paired uncertainty", "inverse_use":"secondary consistency check only; not independent of its component values"},
    {"id":"D13C_MIXING_FRACTION", "observable":"Organic-derived DIC fraction", "status":"derived two-endmember mixing estimate", "source":"Derived from carbonate δ13C and selected endmembers; not a Jian direct measurement", "model_state":"OM and DIC concentrations, but no isotope pools", "operator":"same post-processing equation as carbonate δ13C closure", "availability":"conditionally computable only if endmembers and offsets are specified", "uniqueness":"non-unique with endmember choice and unmodelled microbial/carbonate fractionation", "inverse_use":"do not use as independent calibration target"},
    {"id":"D18O_WATER_RANGE", "observable":"Inferred water δ18O range", "status":"derived inverse range", "source":"Derived from carbonate δ18O over assumed temperature range", "model_state":"none; oxygen isotopes and temperature history absent", "operator":"carbonate-water fractionation relation over an explicit temperature prior", "availability":"unavailable in current solver", "uniqueness":"non-unique temperature-water covariance", "inverse_use":"do not use as independent calibration target"},
    {"id":"VFA_PHPATHWAY", "observable":"VFA-mediated pH suppression and anaerobic sulphide oxidation pathway", "status":"proposed mechanism, not measured observation", "source":"Jian et al. (2026), discussion and reaction scheme", "model_state":"synthetic OM/H+/S species and scenario switches", "operator":"scenario-conditioned model state only; no measured VFA or pH target", "availability":"mechanism switch only", "uniqueness":"highly non-unique: microbial rates, electron acceptors, transport and substrate supply are unconstrained", "inverse_use":"hypothesis scenario; never treat as observed or calibrated"},
    {"id":"MECHANICAL_PROPERTIES", "observable":"Burial stress, elastic moduli, strain and cementation timing", "status":"not measured in the source study", "source":"No specimen-specific mechanical data identified", "model_state":"synthetic stress, stiffness, strain and porosity feedback", "operator":"none available without independent mechanical observations or literature priors", "availability":"unavailable", "uniqueness":"structurally non-identifiable from mineral/isotope observations alone", "inverse_use":"do not infer fossil mechanics from current data"},
]


def phase_inventory_operator(profile: Sequence[Mapping[str, float]], *, cells: Iterable[int],
                             phase: str, detection_threshold_mol_m2: float | None = None) -> dict[str, object]:
    """Integrate phase inventory over a user-specified sample mask.

    `calcite_mol_m3_bulk` or `barite_mol_m3_bulk` is mol per bulk cubic metre;
    integration over cell thickness returns mol m-2. The sample mask and
    optional detection threshold are external measurement-design inputs.
    """
    field = {"calcite":"calcite_mol_m3_bulk", "barite":"barite_mol_m3_bulk"}.get(phase)
    if field is None:
        raise ValueError("phase must be 'calcite' or 'barite'")
    selected = set(int(i) for i in cells)
    chosen = [r for r in profile if int(r["cell"]) in selected]
    if not chosen:
        raise ValueError("cells must select at least one profile cell")
    if len(chosen) != len(selected):
        raise ValueError("sample mask contains a cell index absent from the profile")
    if any(float(r[field]) < -1e-12 for r in chosen):
        raise ValueError("negative mineral inventory in supplied model profile")
    amount = sum(float(r[field]) * float(r["dx_m"]) for r in chosen)
    result = {"phase":phase, "selected_cells":sorted(selected), "inventory_mol_m2":amount,
              "detection_threshold_mol_m2":detection_threshold_mol_m2,
              "presence_prediction":None if detection_threshold_mol_m2 is None else amount >= detection_threshold_mol_m2,
              "interpretation":"inventory only; presence is conditional on an externally specified detection threshold and sample mask"}
    if detection_threshold_mol_m2 is not None and detection_threshold_mol_m2 < 0:
        raise ValueError("detection threshold must be non-negative")
    return result


def mineral_pattern_operator(profile: Sequence[Mapping[str, float]], *,
                             boundary_fraction: tuple[float,float]=(0.20,0.30),
                             threshold_fraction: tuple[float,float]=(0.01,0.50),
                             boundary_nodes: int=41, threshold_nodes: int=41,
                             logistic_width_ln: float=0.5) -> dict[str,float]:
    """Prior-average soft mineral detection across uncertain 1D masks/thresholds."""
    import numpy as np
    rows=sorted(profile,key=lambda r:int(r["cell"]))
    x=np.asarray([float(r["x_m"]) for r in rows]); dx=np.asarray([float(r["dx_m"]) for r in rows])
    length=float(np.sum(dx))
    phase_fields={"barite":"barite_mol_m3_bulk","calcite":"calcite_mol_m3_bulk"}
    fields={p:np.asarray([float(r[col]) for r in rows]) for p,col in phase_fields.items()}
    if any(np.any(v < -1e-12) for v in fields.values()):
        raise ValueError("negative phase inventory in model profile")
    boundaries=np.linspace(*boundary_fraction,boundary_nodes)
    qvals=np.geomspace(*threshold_fraction,threshold_nodes)
    scores={"p_barite_detect_bone_proxy":[],"p_barite_detect_matrix_proxy":[],"p_calcite_detect_matrix_proxy":[]}
    for b in boundaries:
        bone=x < b*length; matrix=~bone
        for phase,mask,key in [("barite",bone,"p_barite_detect_bone_proxy"),
                               ("barite",matrix,"p_barite_detect_matrix_proxy"),
                               ("calcite",matrix,"p_calcite_detect_matrix_proxy")]:
            total=float(np.sum(fields[phase]*dx)); amount=float(np.sum(fields[phase][mask]*dx[mask]))
            if total <= 0 or amount <= 0:
                scores[key].append(0.0); continue
            fraction=min(max(amount/total,1e-300),1.0)
            # q is a fraction of this phase's total model inventory. The logistic
            # transition smooths a threshold with no known analytical detection limit.
            probs=1.0/(1.0+np.exp(-np.clip((np.log(fraction)-np.log(qvals))/logistic_width_ln,-700,700)))
            scores[key].append(float(np.mean(probs)))
    return {k:float(np.mean(v)) for k,v in scores.items()}


def carbon_isotope_operator(*, f_org: float, delta_org: float=-31.7,
                             delta_background_calcite: float=-3.0,
                             epsilon_srb: float=0.0,
                             epsilon_calcite_hco3: float=1.0) -> float:
    """Post-process a two-source δ13C mix; no isotope pools are in the solver."""
    if not 0.0 <= f_org <= 1.0:
        raise ValueError("f_org must lie in [0, 1]")
    delta_org_derived_dIC=delta_org+epsilon_srb
    # Carbonate-DIC fractionation acts on the mixed DIC pool, not only on
    # the organic-derived fraction.
    delta_mixed_dIC=(f_org*delta_org_derived_dIC
                     +(1.0-f_org)*delta_background_calcite)
    return float(delta_mixed_dIC+epsilon_calcite_hco3)


def likelihood_eligibility() -> list[dict[str, object]]:
    """Return current fit eligibility without silently promoting derived data."""
    return [{"id":r["id"], "eligible_now":False,
             "reason":"no site-registered sample support and calibrated forward map in the synthetic benchmark"}
            for r in OBSERVATION_CATALOG]


def write_catalog(path: str | None = None) -> None:
    import csv
    from pathlib import Path
    target = Path(path) if path else Path(__file__).resolve().parent / "results" / "observation_operator_catalog.csv"
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(OBSERVATION_CATALOG[0].keys()))
        writer.writeheader(); writer.writerows(OBSERVATION_CATALOG)


if __name__ == "__main__":
    import csv
    import json
    write_catalog()
    root = __import__("pathlib").Path(__file__).resolve().parent
    out = root / "results"
    eligible = likelihood_eligibility()
    with (out / "observation_fit_eligibility.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(eligible[0].keys()))
        w.writeheader(); w.writerows(eligible)
    status_counts = {}
    for row in OBSERVATION_CATALOG:
        status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
    summary = {"operator_mappings":len(OBSERVATION_CATALOG), "current_inverse_fit_eligible":sum(bool(x["eligible_now"]) for x in eligible),
               "evidence_status_counts":status_counts,
               "derived_non_unique_flagged":sum("derived" in str(x["status"]).lower() or "non-unique" in str(x["uniqueness"]).lower() for x in OBSERVATION_CATALOG),
               "source_conflicts_preserved":True,
               "interpretation":"No published quantity is currently used in a likelihood because no fossil-registered sample support and calibrated forward map exist in the synthetic benchmark."}
    (out / "observation_operator_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(OBSERVATION_CATALOG)} mappings; current inverse-fit eligible count = 0.")
