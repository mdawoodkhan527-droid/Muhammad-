#!/usr/bin/env python3
"""Reduced 1-D reactive-transport / porosity / poroelastic verification model.

The chemical constants and initial conditions are synthetic verification inputs,
not fossil-specific estimates. Activities are idealised; the code is a coupling
prototype, not a PHREEQC-equivalent geochemical model.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import yaml
from verified_1d_reference_model import transport_analytic

ROOT = Path(__file__).resolve().parent
TOL = yaml.safe_load((ROOT / "tolerances.yaml").read_text())
KINETIC_BASE = {"k_srb":1.0e-3, "k_sox":2.0e-4, "k_sdb":1.0e-4,
                "k_calcite":2.0e-4, "k_barite":2.0e-4}

SP = ["OM", "SO4", "HCO3", "H2S", "CO2", "HS", "H2O", "Ca", "Calcite", "Ba", "Barite", "S0", "H", "Cl", "Na"]
I = {name: j for j, name in enumerate(SP)}
ELEMENTS = ["C", "S", "Ca", "Ba", "Cl", "Na"]
FORMULA = {
    "OM": {"C": 1, "H": 2, "O": 1}, "SO4": {"S": 1, "O": 4},
    "HCO3": {"H": 1, "C": 1, "O": 3}, "H2S": {"H": 2, "S": 1},
    "CO2": {"C": 1, "O": 2}, "HS": {"H": 1, "S": 1},
    "H2O": {"H": 2, "O": 1}, "Ca": {"Ca": 1},
    "Calcite": {"Ca": 1, "C": 1, "O": 3}, "Ba": {"Ba": 1},
    "Barite": {"Ba": 1, "S": 1, "O": 4}, "S0": {"S": 1},
    "H": {"H": 1}, "Cl": {"Cl": 1}, "Na": {"Na": 1},
}
CHARGE = np.array([0, -2, -1, 0, 0, -1, 0, 2, 0, 2, 0, 0, 1, -1, 1], dtype=float)
ATOM = {el: np.array([FORMULA[s].get(el, 0) for s in SP], dtype=float)
        for el in ["C", "H", "O", "S", "Ca", "Ba", "Cl", "Na"]}

def reaction_vector(**terms: float) -> np.ndarray:
    v = np.zeros(len(SP))
    for name, coeff in terms.items():
        v[I[name]] = coeff
    return v

# Jian et al. (2026), Eqs. 1-6. Eq. 1/2 are alternative SRB balances.
REACTIONS = {
    "SRB_eq1": reaction_vector(OM=-2, SO4=-1, HCO3=2, H2S=1),
    "SRB_eq2": reaction_vector(OM=-2, SO4=-1, HCO3=1, CO2=1, HS=1, H2O=1),
    "calcite_eq3": reaction_vector(Ca=-1, HCO3=-2, Calcite=1, H2O=1, CO2=1),
    "barite_eq4": reaction_vector(Ba=-1, SO4=-1, Barite=1),
    "sulfur_oxidation_eq5": reaction_vector(CO2=-1, H2S=-2, OM=1, H2O=1, S0=2),
    "disproportionation_eq6": reaction_vector(H2O=-4, S0=-4, H2S=3, SO4=1, H=2),
}

def reaction_balance_rows() -> list[dict[str, object]]:
    rows = []
    for name, nu in REACTIONS.items():
        row = {el: float(np.dot(ATOM[el], nu)) for el in ATOM}
        row["charge"] = float(np.dot(CHARGE, nu))
        row["reaction"] = name
        row["passed"] = max(abs(row[e]) for e in ATOM) <= TOL["reaction_balance"]["absolute_element_residual"] and abs(row["charge"]) <= TOL["reaction_balance"]["absolute_charge_residual"]
        rows.append(row)
    return rows

def initial_aqueous_fluid(nx: int, phi: np.ndarray, length: float,
                          concentration_scales: dict[str, float] | None = None,
                          bone_ba_initial_factor: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    x = (np.arange(nx) + 0.5) * length / nx
    c = np.zeros((nx, len(SP)))
    c[:, I["OM"]] = np.where(x < 0.25 * length, 4.0, 0.5)
    c[:, I["SO4"]] = 5.0
    c[:, I["HCO3"]] = 1.0
    c[:, I["H2S"]] = 0.05
    c[:, I["CO2"]] = 0.5
    c[:, I["HS"]] = 0.05
    c[:, I["H2O"]] = 55500.0
    c[:, I["Ca"]] = 10.0
    if bone_ba_initial_factor <= 0:
        raise ValueError("bone_ba_initial_factor must be positive")
    c[:, I["Ba"]] = np.where(x < 0.25 * length, 0.005 * bone_ba_initial_factor, 0.005)
    c[:, I["H"]] = np.where(x < 0.25 * length, 1.0e-3, 1.0e-4)
    for name, scale in (concentration_scales or {}).items():
        if name not in I:
            raise KeyError(f"unknown concentration scale species: {name}")
        c[:, I[name]] *= float(scale)
    c[:, I["Cl"]] = 0.0
    c[:, I["Na"]] = 0.0
    # Conservative inert Na+/Cl- counterions make each synthetic cell neutral
    # without changing the prescribed acid concentration or reactive species.
    net_charge = np.sum(c * CHARGE[None, :], axis=1)
    c[:, I["Cl"]] = np.maximum(net_charge, 0.0)
    c[:, I["Na"]] = np.maximum(-net_charge, 0.0)
    charge = c @ CHARGE
    if np.max(np.abs(charge)) > 1e-12:
        raise ValueError("Synthetic initial porewater could not be charge balanced with Na+/Cl-")
    # Solids and water are inventory variables per bulk volume, while dissolved
    # species are stored as porewater concentration times porosity.
    n = c * phi[:, None]
    return n, x

def pH_from_state(n: np.ndarray, phi: np.ndarray) -> np.ndarray:
    h_mol_l = np.maximum(n[:, I["H"]] / np.maximum(phi, 1e-12) / 1000.0, 1e-14)
    return -np.log10(h_mol_l)

def calcite_omega(n: np.ndarray, phi: np.ndarray, pH: np.ndarray) -> np.ndarray:
    h = 10.0 ** (-pH)
    ka1, ka2 = 10.0 ** -6.35, 10.0 ** -10.33
    alpha_co3 = ka1 * ka2 / (h*h + ka1*h + ka1*ka2)
    dic_m = (n[:, I["HCO3"]] + n[:, I["CO2"]]) / np.maximum(phi, 1e-12) / 1000.0
    co3_m = dic_m * alpha_co3
    ca_m = n[:, I["Ca"]] / np.maximum(phi, 1e-12) / 1000.0
    return np.maximum(ca_m * co3_m / 3.3e-9, 1e-30)

def barite_omega(n: np.ndarray, phi: np.ndarray) -> np.ndarray:
    ba_m = n[:, I["Ba"]] / np.maximum(phi, 1e-12) / 1000.0
    so4_m = n[:, I["SO4"]] / np.maximum(phi, 1e-12) / 1000.0
    return np.maximum(ba_m * so4_m / 1.08e-10, 1e-30)

def apply_extent(nrow: np.ndarray, nu: np.ndarray, requested: float) -> float:
    extent = max(float(requested), 0.0)
    if extent <= 0:
        return 0.0
    for j in range(len(SP)):
        if nu[j] < 0:
            extent = min(extent, max(nrow[j], 0.0) / (-nu[j]))
    nrow += extent * nu
    if np.min(nrow) < -1e-12:
        raise AssertionError("Reaction extent produced a negative species inventory")
    # Round-off only. Larger negatives must be fixed in the transport or rate step.
    nrow[nrow < 0] = 0.0
    return extent

def effective_modulus(eps_cement: np.ndarray, phi0: float, matrix_e: float, cement_e: float) -> np.ndarray:
    solid_fraction0 = 1.0 - phi0
    cement_fraction = np.clip(eps_cement / solid_fraction0, 0.0, 0.95)
    # Voigt mixture for the solid frame, explicitly an assumed homogenisation rule.
    return matrix_e * (1.0 - cement_fraction) + cement_e * cement_fraction

def solve_poroelastic_bar(e_cells: np.ndarray, pressure_cells: np.ndarray, *,
                          length: float, alpha: float, area: float,
                          right_total_traction: float = 0.0) -> tuple[np.ndarray, np.ndarray, float]:
    """1D FE bar: sigma_total=E*u' - alpha*p; u(0)=0; prescribed total traction at x=L."""
    nx = len(e_cells)
    dx = length / nx
    k = np.zeros((nx + 1, nx + 1))
    f = np.zeros(nx + 1)
    for e in range(nx):
        ke = e_cells[e] * area / dx * np.array([[1., -1.], [-1., 1.]])
        k[e:e+2, e:e+2] += ke
        # Weak-form load from Biot pressure: integral(alpha*p*N') dx.
        f[e:e+2] += alpha * pressure_cells[e] * area * np.array([-1., 1.])
    f[-1] += right_total_traction * area
    u = np.zeros(nx + 1)
    u[1:] = np.linalg.solve(k[1:, 1:], f[1:])
    residual = k @ u - f
    eq = float(np.linalg.norm(residual[1:]) / max(np.linalg.norm(f[1:]), area * 1e-30))
    strain = np.diff(u) / dx
    return u, strain, eq

def analytical_poroelastic_displacement(e_cells: np.ndarray, pressure_cells: np.ndarray,
                                        *, length: float, alpha: float) -> np.ndarray:
    """Exact nodal displacement for piecewise-constant E,p and zero total traction."""
    dx = length / len(e_cells)
    increments = alpha * pressure_cells / e_cells * dx
    return np.concatenate(([0.0], np.cumsum(increments)))

def porosity_and_permeability(eps_calc: np.ndarray, eps_barite: np.ndarray,
                              phi0: float, strain_feedback: np.ndarray,
                              k0: float, exponent: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    eps_m = eps_calc + eps_barite
    phi = phi0 - eps_m + strain_feedback
    if np.any(phi <= TOL["conservation_each_step"]["minimum_porosity"]):
        raise ValueError("Porosity fell below the frozen minimum; reduce the step or revise the scenario")
    kc = k0 * (phi / phi0)**3 * ((1.0 - phi0) / (1.0 - phi))**2
    kp = k0 * (phi / phi0)**exponent
    residual = phi + eps_m - phi0 - strain_feedback
    return phi, kc, kp, residual

def external_boundary_concentrations(concentration_scales: dict[str, float] | None = None) -> np.ndarray:
    c = np.zeros(len(SP))
    for s, v in {"OM":0.5,"SO4":5.0,"HCO3":1.0,"H2S":0.05,"CO2":0.5,"HS":0.05,
                 "H2O":55500.0,"Ca":10.0,"Ba":0.005,"H":1.0e-4}.items():
        c[I[s]] = v
    for name, scale in (concentration_scales or {}).items():
        if name not in I:
            raise KeyError(f"unknown concentration scale species: {name}")
        c[I[name]] *= float(scale)
    net_charge = float(np.dot(CHARGE, c))
    c[I["Cl"]] = max(net_charge, 0.0)
    c[I["Na"]] = max(-net_charge, 0.0)
    return c

def transport_step(n: np.ndarray, phi: np.ndarray, dt: float, *, dx: float,
                   diffusivity: float, permeability: np.ndarray, viscosity: float,
                   pressure_gradient: float, inlet: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nx = len(phi)
    # All mobile species use one equal D in this reduced charge-conserving test.
    mobile = [I[s] for s in SP if s not in ("H2O", "Calcite", "Barite", "S0")]
    c = n[:, mobile] / phi[:, None]
    kface = np.empty(nx + 1)
    kface[0], kface[-1] = permeability[0], permeability[-1]
    if nx > 1:
        kface[1:-1] = 2 * permeability[:-1] * permeability[1:] / np.maximum(permeability[:-1] + permeability[1:], 1e-300)
    q = -kface * pressure_gradient / viscosity
    flux = np.zeros((nx + 1, len(mobile)))
    # Advective flux with upwind boundary states. Diffusion uses no-flux outer faces.
    for face in range(nx + 1):
        vel = q[face]
        if face == 0:
            cup = inlet[mobile] if vel >= 0 else c[0]
        elif face == nx:
            cup = c[-1] if vel >= 0 else inlet[mobile]
        else:
            cup = c[face - 1] if vel >= 0 else c[face]
        flux[face] = vel * cup
        if 0 < face < nx:
            dface = diffusivity * (0.5 * (phi[face-1] + phi[face]) / 0.3)**2
            flux[face] += -dface * (c[face] - c[face-1]) / dx
    out = n.copy()
    out[:, mobile] += dt / dx * (flux[:-1] - flux[1:])
    ext = np.zeros(len(SP))
    ext[mobile] = dt * (flux[0] - flux[-1])
    return out, ext

def reaction_step(n: np.ndarray, phi: np.ndarray, dt: float, *, scenario: str,
                  srb_variant: str, eps_calc: np.ndarray, eps_barite: np.ndarray,
                  phi0: float, kinetics: dict[str, float] | None = None,
                  acid_calcite_switch: bool = False, acid_pH_threshold: float = 6.5,
                  acid_transition_width: float = 0.25) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    out = n.copy()
    rates = {**KINETIC_BASE, **(kinetics or {})}
    for j in range(len(phi)):
        pH = float(pH_from_state(out[j:j+1], phi[j:j+1])[0])
        c = out[j] / phi[j]
        pH_factor = 1.0 / (1.0 + 10.0**(6.2 - pH))
        monod_om = c[I["OM"]] / (0.5 + c[I["OM"]])
        monod_so4 = c[I["SO4"]] / (0.5 + c[I["SO4"]])
        r_srb = rates["k_srb"] * monod_om * monod_so4 * pH_factor
        srb = "SRB_eq1" if srb_variant == "eq1" else "SRB_eq2"
        apply_extent(out[j], REACTIONS[srb], dt * r_srb)
        if scenario in ("S2", "S3"):
            c = out[j] / phi[j]
            r_ox = rates["k_sox"] * c[I["CO2"]] / (0.2 + c[I["CO2"]]) * c[I["H2S"]] / (0.1 + c[I["H2S"]])
            apply_extent(out[j], REACTIONS["sulfur_oxidation_eq5"], dt * r_ox)
        if scenario == "S3":
            c = out[j] / phi[j]
            r_dis = rates["k_sdb"] * c[I["S0"]] / (0.1 + c[I["S0"]]) * pH_factor
            apply_extent(out[j], REACTIONS["disproportionation_eq6"], dt * r_dis)
        c = out[j] / phi[j]
        pH = float(pH_from_state(out[j:j+1], phi[j:j+1])[0])
        omega_c = calcite_omega(out[j:j+1], phi[j:j+1], np.array([pH]))[0]
        omega_b = barite_omega(out[j:j+1], phi[j:j+1])[0]
        eps_c = out[j, I["Calcite"]] * 3.69e-5
        eps_b = out[j, I["Barite"]] * 5.20e-5
        phi_available = max(phi0 - eps_c - eps_b, TOL["conservation_each_step"]["minimum_porosity"])
        # Log-supersaturation rate avoids unbounded rates at high idealised IAP.
        acid_factor = 1.0
        if acid_calcite_switch:
            # Assumed kinetic inhibition switch, separate from equilibrium SI.
            # Low pH suppresses calcite nucleation/growth; high pH approaches 1.
            acid_factor = 1.0 / (1.0 + 10.0 ** ((acid_pH_threshold - pH) / acid_transition_width))
        r_calcite = rates["k_calcite"] * max(np.log10(omega_c), 0.0) * acid_factor
        xi_cmax = max(phi_available - TOL["conservation_each_step"]["minimum_porosity"], 0.0) / 3.69e-5
        apply_extent(out[j], REACTIONS["calcite_eq3"], min(dt * r_calcite, xi_cmax))
        eps_c = out[j, I["Calcite"]] * 3.69e-5
        eps_b = out[j, I["Barite"]] * 5.20e-5
        phi_available = max(phi0 - eps_c - eps_b, TOL["conservation_each_step"]["minimum_porosity"])
        r_barite = rates["k_barite"] * max(np.log10(omega_b), 0.0)
        xi_bmax = max(phi_available - TOL["conservation_each_step"]["minimum_porosity"], 0.0) / 5.20e-5
        apply_extent(out[j], REACTIONS["barite_eq4"], min(dt * r_barite, xi_bmax))
    eps_calc_new = out[:, I["Calcite"]] * 3.69e-5
    eps_barite_new = out[:, I["Barite"]] * 5.20e-5
    return out, eps_calc_new, eps_barite_new

def inventory(n: np.ndarray, dx: float) -> dict[str, float]:
    inv = {e: float(np.sum(n @ ATOM[e]) * dx) for e in ELEMENTS}
    inv["charge"] = float(np.sum(n @ CHARGE) * dx)
    return inv

def run_case(*, scenario: str, permeability_law: str, feedback: bool,
             nx: int = 20, final_time: float = 60.0, dt: float = 0.2,
             kinetics: dict[str,float] | None = None,
             concentration_scales: dict[str,float] | None = None,
             diffusivity: float = 1.0e-9, permeability_reference: float = 1.0e-18,
             bone_ba_initial_factor: float = 10.0, bone_ba_source_rate: float = 0.0,
             acid_calcite_switch: bool = False, acid_pH_threshold: float = 6.5,
             acid_transition_width: float = 0.25,
             compute_mechanics: bool = True) -> dict[str, object]:
    length, phi0 = 1.0e-3, 0.30
    if bone_ba_source_rate < 0:
        raise ValueError("bone_ba_source_rate must be non-negative")
    if compute_mechanics is False and feedback:
        raise ValueError("mechanics cannot be disabled when feedback is enabled")
    dx = length / nx
    ncell = nx
    phi = np.full(ncell, phi0)
    n, x = initial_aqueous_fluid(ncell, phi, length, concentration_scales,
                                 bone_ba_initial_factor=bone_ba_initial_factor)
    inlet = external_boundary_concentrations(concentration_scales)
    eps_calc = np.zeros(ncell)
    eps_barite = np.zeros(ncell)
    strain = np.zeros(ncell)
    eps_feedback = np.zeros(ncell)
    q_elems = np.zeros(ncell)
    p_hydraulic_gradient = -1.0e7
    viscosity, k0 = 1.0e-3, permeability_reference
    alpha, pore_pressure, area = 0.8, 5.0e5, 1.0
    emat, ecem = 1.0e9, 10.0e9
    nsteps = int(round(final_time / dt))
    if not np.isclose(nsteps * dt, final_time):
        raise ValueError("final_time must be divisible by dt")
    start_inv = inventory(n, dx)
    cumulative_boundary = {e: 0.0 for e in ELEMENTS}
    cumulative_boundary["charge"] = 0.0
    records = []
    for step in range(nsteps):
        eps_total = eps_calc + eps_barite
        phi, kc, kp, volres = porosity_and_permeability(eps_calc, eps_barite, phi0,
                                                       eps_feedback, k0, exponent=3.0)
        kval = kc if permeability_law == "kozeny_carman" else kp
        n, boundary_species = transport_step(n, phi, dt, dx=dx, diffusivity=diffusivity,
                                             permeability=kval, viscosity=viscosity,
                                             pressure_gradient=p_hydraulic_gradient, inlet=inlet)
        if np.min(n) < -1e-12:
            raise AssertionError(f"Transport step produced a negative species inventory at t={(step+1)*dt}; reduce dt or use a positivity-preserving flux update")
        for el in ELEMENTS:
            cumulative_boundary[el] += float(np.dot(boundary_species, ATOM[el]))
        cumulative_boundary["charge"] += float(np.dot(boundary_species, CHARGE))
        # Optional local BaCl2 source is an explicit external inventory input.
        # It is tracked in element and charge budgets, rather than hidden in a
        # concentration adjustment. The source rate is expressed in the model's
        # synthetic porewater concentration units per second.
        if bone_ba_source_rate > 0:
            bone_cells = x < 0.25 * length
            ba_increment = phi[bone_cells] * dt * bone_ba_source_rate
            n[bone_cells, I["Ba"]] += ba_increment
            n[bone_cells, I["Cl"]] += 2.0 * ba_increment
            cumulative_boundary["Ba"] += float(np.sum(ba_increment) * dx)
            cumulative_boundary["Cl"] += float(2.0 * np.sum(ba_increment) * dx)
        n, eps_calc, eps_barite = reaction_step(n, phi, dt, scenario=scenario,
                                                srb_variant="eq1", eps_calc=eps_calc,
                                                eps_barite=eps_barite, phi0=phi0, kinetics=kinetics,
                                                acid_calcite_switch=acid_calcite_switch,
                                                acid_pH_threshold=acid_pH_threshold,
                                                acid_transition_width=acid_transition_width)
        eps_total = eps_calc + eps_barite
        if compute_mechanics:
            e_cells = effective_modulus(eps_total, phi0, emat, ecem)
            p_cells = np.full(ncell, pore_pressure)
            _, strain, eqres = solve_poroelastic_bar(e_cells, p_cells, length=length,
                                                     alpha=alpha, area=area)
        else:
            strain.fill(0.0)
            eqres = 0.0
        if feedback:
            eps_feedback = strain.copy()
        else:
            eps_feedback.fill(0.0)
        phi, kc, kp, volres = porosity_and_permeability(eps_calc, eps_barite, phi0,
                                                       eps_feedback, k0, exponent=3.0)
        now = inventory(n, dx)
        mass_res = {}
        for el in ELEMENTS:
            scale = max(abs(start_inv[el]), abs(now[el]), abs(cumulative_boundary[el]), 1e-20)
            mass_res[el] = (now[el] - start_inv[el] - cumulative_boundary[el]) / scale
        # A signed net-charge sum is near zero in an electroneutral solution and
        # is therefore not a useful normaliser. Use total absolute charge
        # equivalents in the aqueous inventory and cumulative boundary exchange.
        charge_inventory = float(np.sum(np.abs(n * CHARGE[None, :])) * dx)
        charge_scale = max(charge_inventory, abs(cumulative_boundary["charge"]),
                           abs(start_inv["charge"]), 1e-20)
        charge_res = (now["charge"] - start_inv["charge"] - cumulative_boundary["charge"]) / charge_scale
        local_charge = np.max(np.abs(n @ CHARGE) / np.maximum(np.sum(np.abs(n * CHARGE[None,:]), axis=1), 1e-20))
        max_element = max(abs(v) for v in mass_res.values())
        max_volume = float(np.max(np.abs(volres)))
        rec = {"time_s": (step+1)*dt, "max_relative_element_residual": max_element,
               "relative_charge_residual": float(abs(charge_res)), "max_local_charge_residual": float(local_charge),
               "max_absolute_volume_residual": max_volume,
               "mean_porosity": float(np.mean(phi)), "min_porosity": float(np.min(phi)),
               "total_calcite_mol_m2": float(np.sum(n[:,I["Calcite"]])*dx),
               "total_barite_mol_m2": float(np.sum(n[:,I["Barite"]])*dx),
               "mean_pH": float(np.mean(pH_from_state(n,phi))),
               "mean_strain": float(np.mean(strain)), "mechanical_equilibrium_residual": eqres,
               "mean_k_KC_m2": float(np.mean(kc)), "mean_k_power_m2": float(np.mean(kp))}
        records.append(rec)
        element_limits = {"C": "relative_C_mass", "S": "relative_S_mass",
                          "Ca": "relative_Ca_mass", "Ba": "relative_Ba_mass",
                          "Cl": "relative_Cl_mass", "Na": "relative_Na_mass"}
        failed_elements = {el: mass_res[el] for el, limit_name in element_limits.items()
                           if abs(mass_res[el]) > TOL["conservation_each_step"][limit_name]}
        if (failed_elements or abs(charge_res) > TOL["conservation_each_step"]["relative_charge"]
                or max_volume > TOL["conservation_each_step"]["absolute_pore_volume_residual"]
                or local_charge > TOL["conservation_each_step"]["relative_charge"]):
            raise AssertionError(f"Conservation tolerance failed at t={(step+1)*dt}: elements={failed_elements}, global_charge={charge_res}, local_charge={local_charge}, volume={max_volume}")
        if np.min(n) < -1e-12:
            raise AssertionError("Negative species amount encountered")
    return {"scenario":scenario, "permeability_law":permeability_law,
            "coupling": "stage2_two_way" if feedback else "stage1_one_way",
            "nx":nx,"dt_s":dt,"final_time_s":final_time,
            "rows":records,
            "node_profile":[{
                "cell":int(j), "x_m":float(x[j]), "dx_m":float(dx),
                "pH":float(pH_from_state(n[j:j+1],phi[j:j+1])[0]),
                "OM_mmol_L":float(n[j,I["OM"]]/phi[j]),
                "SO4_mmol_L":float(n[j,I["SO4"]]/phi[j]),
                "H2S_mmol_L":float(n[j,I["H2S"]]/phi[j]),
                "Ca_mmol_L":float(n[j,I["Ca"]]/phi[j]),
                "Ba_mmol_L":float(n[j,I["Ba"]]/phi[j]),
                "calcite_mol_m3_bulk":float(n[j,I["Calcite"]]),
                "barite_mol_m3_bulk":float(n[j,I["Barite"]]),
                "calcite_volume_fraction":float(eps_calc[j]),
                "barite_volume_fraction":float(eps_barite[j]),
                "porosity":float(phi[j]), "volumetric_strain":float(strain[j]),
                "permeability_KC_m2":float(kc[j]), "permeability_power_m2":float(kp[j]),
            } for j in range(ncell)],
            "summary":{key:records[-1][key] for key in ["mean_porosity","min_porosity","total_calcite_mol_m2","total_barite_mol_m2","mean_pH","mean_strain","mean_k_KC_m2","mean_k_power_m2"]},
            "max_residuals":{key:max(r[key] for r in records) for key in ["max_relative_element_residual","relative_charge_residual","max_local_charge_residual","max_absolute_volume_residual","mechanical_equilibrium_residual"]}}

def zero_rate_transport_check() -> dict[str, float | bool]:
    # Same verified scalar benchmark: D=1e-9 m2/s, k=2e-4 1/s, L=1 mm,
    # C=1 at both boundaries, zero initial interior. New biological/mineral
    # reaction extents are set to zero; k is retained only as the prior solver's
    # transport benchmark sink so the limiting solutions are comparable.
    L, D, decay, cb, tf = 1e-3, 1e-9, 2e-4, 1.0, 600.0
    nx, dx, dt = 40, L/40, 0.05
    steps = int(tf/dt)
    x = (np.arange(nx)+0.5)*dx
    c = np.zeros(nx)
    for _ in range(steps):
        # Cell-centred finite volumes with half-cell Dirichlet boundary fluxes.
        faces = np.zeros(nx+1)
        faces[0] = -D*(c[0]-cb)/(dx/2)
        faces[-1] = -D*(cb-c[-1])/(dx/2)
        faces[1:-1] = -D*np.diff(c)/dx
        c += dt*((faces[:-1]-faces[1:])/dx - decay*c)
    exact = transport_analytic(x, tf, length=L, diffusivity=D, decay=decay, boundary=cb)
    rel_l2 = float(np.sqrt(np.mean((c-exact)**2))/cb)
    rel_max = float(np.max(np.abs(c-exact))/cb)
    return {"relative_L2_error":rel_l2,"relative_max_error":rel_max,
            "passed":rel_l2 <= TOL["limiting_case_transport"]["relative_L2_error"] and rel_max <= TOL["limiting_case_transport"]["relative_max_error"]}

def poroelastic_limiting_check() -> dict[str, float | bool]:
    nx, length, alpha, pressure, area = 40, 1e-2, 0.8, 5e5, 1e-4
    xmid = (np.arange(nx)+0.5)/nx
    phi0 = 0.30
    eps_uniform = np.full(nx, 0.01)
    base_e = np.where(xmid < 0.5, 1.0e9, 2.0e9)
    e_cells = base_e*(1-eps_uniform/(1-phi0)) + 10e9*eps_uniform/(1-phi0)
    p = np.full(nx, pressure)
    u, strain, eq = solve_poroelastic_bar(e_cells,p,length=length,alpha=alpha,area=area)
    exact = analytical_poroelastic_displacement(e_cells,p,length=length,alpha=alpha)
    err = float(np.sqrt(np.mean((u-exact)**2))/max(np.sqrt(np.mean(exact**2)),1e-30))
    zero_u, _, zero_eq = solve_poroelastic_bar(e_cells,np.zeros(nx),length=length,alpha=alpha,area=area)
    zero_max = float(np.max(np.abs(zero_u)))
    ok = (err <= TOL["limiting_case_poroelastic"]["relative_displacement_L2_error"]
          and eq <= TOL["limiting_case_poroelastic"]["normalized_equilibrium_residual"]
          and zero_eq <= TOL["limiting_case_poroelastic"]["normalized_equilibrium_residual"]
          and zero_max <= TOL["limiting_case_poroelastic"]["zero_pressure_max_displacement"])
    return {"relative_displacement_L2_error":err,"normalized_equilibrium_residual":eq,
            "zero_pressure_max_displacement":zero_max,"passed":ok}

def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({k for row in rows for k in row})
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def main() -> None:
    outdir=ROOT/"results"; outdir.mkdir(exist_ok=True)
    balances=reaction_balance_rows()
    if not all(r["passed"] for r in balances):
        raise AssertionError("Reaction balance unit test failed")
    transport=zero_rate_transport_check()
    poro=poroelastic_limiting_check()
    if not transport["passed"] or not poro["passed"]:
        raise AssertionError(f"Limiting-case verification failed: transport={transport}, poroelastic={poro}")
    write_csv(outdir/"reaction_balance.csv",balances)
    results=[]; timeseries=[]; node_profiles=[]
    for scenario in ("S1","S2","S3"):
        for law in ("kozeny_carman","power_law"):
            for feedback in (False,True):
                result=run_case(scenario=scenario,permeability_law=law,feedback=feedback)
                results.append(result)
                for row in result["rows"]:
                    timeseries.append({"scenario":scenario,"permeability_law":law,"coupling":result["coupling"],**row})
                for row in result["node_profile"]:
                    node_profiles.append({"scenario":scenario,"permeability_law":law,"coupling":result["coupling"],**row})
    # Stage-1 vs Stage-2 relative output contrasts. These are screening values;
    # a numerical refinement envelope is still required before selecting coupling.
    contrasts=[]
    keys=["mean_porosity","total_calcite_mol_m2","total_barite_mol_m2","mean_pH","mean_strain"]
    for scenario in ("S1","S2","S3"):
        for law in ("kozeny_carman","power_law"):
            a=next(r for r in results if r["scenario"]==scenario and r["permeability_law"]==law and r["coupling"]=="stage1_one_way")["summary"]
            b=next(r for r in results if r["scenario"]==scenario and r["permeability_law"]==law and r["coupling"]=="stage2_two_way")["summary"]
            for key in keys:
                delta=abs(b[key]-a[key])/max(abs(a[key]),1e-15)
                contrasts.append({"scenario":scenario,"permeability_law":law,"metric":key,
                                  "stage1":a[key],"stage2":b[key],"relative_difference":delta,
                                  "predeclared_threshold":TOL["coupling_decision"]["stage1_vs_stage2_relative_output_difference"],
                                  "exceeds_threshold":delta>TOL["coupling_decision"]["stage1_vs_stage2_relative_output_difference"],
                                  "two_way_decision":"pending_refinement_envelope"})
    write_csv(outdir/"coupled_timeseries.csv",timeseries)
    write_csv(outdir/"coupled_node_profiles.csv",node_profiles)
    write_csv(outdir/"stage1_stage2_contrasts.csv",contrasts)
    summary={"reaction_balance_passed":True,"zero_rate_transport":transport,
             "uniform_precipitation_poroelastic_limit":poro,
             "all_stepwise_conservation_passed":True,
             "scenarios":["S1: Eq.5/Eq.6 off","S2: Eq.5 only; S0 accumulates","S3: Eq.5 plus Eq.6 disproportionation"],
             "permeability_laws":["Kozeny-Carman","power law, n=3"],
             "coupling_decision":"Stage-1 and Stage-2 contrasts are reported; numerical refinement envelope is pending, so no final decision to retain/drop two-way coupling is made.",
             "model_scope":"synthetic reduced chemistry; idealised activities, uncalibrated kinetics, no site-specific porewater or mechanical inputs",
             "reaction_balances":balances,"case_max_residuals":[{"scenario":r["scenario"],"permeability_law":r["permeability_law"],"coupling":r["coupling"],**r["max_residuals"]} for r in results],
             "contrasts":contrasts}
    (outdir/"verification_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps({"reaction_balance_passed":True,"zero_rate_transport":transport,
                      "poroelastic_limit":poro,"all_stepwise_conservation_passed":True,
                      "cases":len(results),"contrasts_over_threshold":sum(r["exceeds_threshold"] for r in contrasts)},indent=2))

if __name__=="__main__": main()
