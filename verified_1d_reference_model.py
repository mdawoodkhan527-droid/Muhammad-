#!/usr/bin/env python3
"""Reduced 1-D transport and mechanics reference model with analytical checks.

This is a numerical verification prototype, not a fossil-specific simulation.
The transport test solves diffusion plus first-order removal in a slab with
fixed equal boundary concentrations. The mechanics test solves a heterogeneous
linear-elastic bar under uniform body force and end traction.

Run: python3 verified_1d_reference_model.py --out results
Dependencies: numpy, scipy
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


def transport_analytic(x: np.ndarray, t: float, *, length: float,
                       diffusivity: float, decay: float, boundary: float,
                       terms: int = 160) -> np.ndarray:
    """Series solution for D*C_xx - k*C, zero initial interior, C(0)=C(L)=Cb."""
    lam = np.sqrt(decay / diffusivity)
    steady = boundary * np.cosh(lam * (x - length / 2.0)) / np.cosh(lam * length / 2.0)
    result = steady.copy()
    for n in range(1, terms + 1):
        wave = n * np.pi / length
        coefficient_integral, _ = quad(
            lambda z: -boundary * np.cosh(lam * (z - length / 2.0))
            / np.cosh(lam * length / 2.0) * np.sin(wave * z),
            0.0, length, epsabs=1e-11, epsrel=1e-11, limit=300)
        coefficient = (2.0 / length) * coefficient_integral
        result += coefficient * np.sin(wave * x) * np.exp(
            -(diffusivity * wave**2 + decay) * t)
    return result


def transport_cn(intervals: int, dt: float, final_time: float, *, length: float,
                 diffusivity: float, decay: float, boundary: float) -> tuple[np.ndarray, np.ndarray]:
    """Crank-Nicolson solution on a uniform grid with fixed Dirichlet endpoints."""
    steps_float = final_time / dt
    steps = int(round(steps_float))
    if not np.isclose(steps * dt, final_time, rtol=0.0, atol=1e-12 * max(final_time, 1.0)):
        raise ValueError("final_time must be an integer multiple of dt")
    dx = length / intervals
    x = np.linspace(0.0, length, intervals + 1)
    n = intervals - 1
    alpha = diffusivity / dx**2
    main = np.full(n, -2.0 * alpha - decay)
    off = np.full(n - 1, alpha)
    operator = diags([off, main, off], [-1, 0, 1], format="csc")
    identity = diags(np.ones(n), 0, format="csc")
    lhs = identity - 0.5 * dt * operator
    rhs_matrix = identity + 0.5 * dt * operator
    boundary_source = np.zeros(n)
    boundary_source[0] += alpha * boundary
    boundary_source[-1] += alpha * boundary
    c = np.zeros(n)
    for _ in range(steps):
        c = spsolve(lhs, rhs_matrix @ c + dt * boundary_source)
    full = np.empty(intervals + 1)
    full[0], full[-1] = boundary, boundary
    full[1:-1] = c
    return x, full


def transport_verification() -> list[dict[str, float]]:
    # Dimensionless groups: diffusion time L^2/D = 1000 s; k*L^2/D = 0.2.
    length, diffusivity, decay, boundary, final_time = 1.0e-3, 1.0e-9, 2.0e-4, 1.0, 600.0
    rows = []
    # Spatial refinement at fixed time step.
    for intervals in [20, 40, 80, 160]:
        dt = 0.5
        x, numerical = transport_cn(intervals, dt, final_time, length=length,
                                    diffusivity=diffusivity, decay=decay, boundary=boundary)
        exact = transport_analytic(x, final_time, length=length, diffusivity=diffusivity,
                                   decay=decay, boundary=boundary)
        rel_l2 = float(np.sqrt(np.mean((numerical - exact)**2)) / boundary)
        max_abs = float(np.max(np.abs(numerical - exact)) / boundary)
        rows.append({"component": "transport", "study": "spatial_refinement", "intervals": intervals, "dt_s": dt,
                     "relative_L2_error": rel_l2, "relative_max_error": max_abs,
                     "criterion": "relative L2 < 0.02 and relative max < 0.03",
                     "passed": rel_l2 < 0.02 and max_abs < 0.03})
    # Time-step refinement at fixed mesh.
    for dt in [8.0, 4.0, 2.0, 1.0, 0.5]:
        intervals = 80
        x, numerical = transport_cn(intervals, dt, final_time, length=length,
                                    diffusivity=diffusivity, decay=decay, boundary=boundary)
        exact = transport_analytic(x, final_time, length=length, diffusivity=diffusivity,
                                   decay=decay, boundary=boundary)
        rel_l2 = float(np.sqrt(np.mean((numerical - exact)**2)) / boundary)
        max_abs = float(np.max(np.abs(numerical - exact)) / boundary)
        rows.append({"component": "transport", "study": "time_refinement", "intervals": intervals, "dt_s": dt,
                     "relative_L2_error": rel_l2, "relative_max_error": max_abs,
                     "criterion": "relative L2 < 0.02 and relative max < 0.03",
                     "passed": rel_l2 < 0.02 and max_abs < 0.03})
    return rows


def mechanical_fem(intervals: int, *, length: float, traction: float,
                   body_force_density: float, area: float,
                   interface_fraction: float, young_left: float,
                   young_right: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Linear FEM bar: d(sigma)/dx + b = 0, sigma(L)=traction, u(0)=0."""
    x = np.linspace(0.0, length, intervals + 1)
    dx = length / intervals
    interface = interface_fraction * length
    stiffness = np.zeros((intervals + 1, intervals + 1))
    load = np.zeros(intervals + 1)
    for e in range(intervals):
        left, right = x[e], x[e + 1]
        left_overlap = max(0.0, min(right, interface) - left)
        right_overlap = dx - left_overlap
        # Integrate E exactly within each element, including elements crossed
        # by the material interface: K_e = A/dx^2 * integral(E dx).
        integral_e_dx = young_left * left_overlap + young_right * right_overlap
        ke = integral_e_dx * area / dx**2 * np.array([[1.0, -1.0], [-1.0, 1.0]])
        stiffness[e:e+2, e:e+2] += ke
        # Consistent element load for constant body force per unit volume.
        load[e:e+2] += body_force_density * area * dx / 2.0
    load[-1] += traction * area
    displacement = np.zeros(intervals + 1)
    displacement[1:] = np.linalg.solve(stiffness[1:, 1:], load[1:])
    residual = stiffness @ displacement - load
    # The first node is fixed, so its residual is the support reaction.
    free_residual = float(np.linalg.norm(residual[1:]) / max(np.linalg.norm(load[1:]), 1e-30))
    return x, displacement, np.array([free_residual])


def mechanical_analytic(x: np.ndarray, *, length: float, traction: float,
                        body_force_density: float, interface_fraction: float,
                        young_left: float, young_right: float) -> np.ndarray:
    """Exact displacement from integrating strain=(traction+b*(L-x))/E(x)."""
    interface = interface_fraction * length

    def strain(s: float) -> float:
        young = young_left if s < interface else young_right
        return (traction + body_force_density * (length - s)) / young

    return np.array([quad(strain, 0.0, float(xi), points=[interface],
                          epsabs=1e-15, epsrel=1e-12, limit=100)[0]
                     for xi in x])


def mechanics_verification() -> list[dict[str, float]]:
    length, area = 1.0e-2, 1.0e-4
    traction, body = 5.0e5, 2.0e7
    # The discontinuity is not aligned with most grid nodes, so elements crossed
    # by the interface test exact coefficient integration.
    interface_fraction, e_left, e_right = 0.437, 1.2e9, 2.4e9
    rows = []
    solutions = {}
    for intervals in [20, 40, 80, 160]:
        x, numerical, residual = mechanical_fem(
            intervals, length=length, traction=traction,
            body_force_density=body, area=area,
            interface_fraction=interface_fraction,
            young_left=e_left, young_right=e_right)
        exact = mechanical_analytic(
            x, length=length, traction=traction,
            body_force_density=body, interface_fraction=interface_fraction,
            young_left=e_left, young_right=e_right)
        solutions[intervals] = numerical
        scale = max(float(np.max(np.abs(exact))), 1e-30)
        rel_l2 = float(np.sqrt(np.mean((numerical - exact)**2)) / scale)
        max_abs = float(np.max(np.abs(numerical - exact)) / scale)
        eq_resid = float(residual[0])
        rows.append({"component": "mechanics", "study": "heterogeneous_mesh_refinement", "intervals": intervals, "dt_s": 0.0,
                     "relative_L2_error": rel_l2, "relative_max_error": max_abs,
                     "relative_equilibrium_residual": eq_resid,
                     "criterion": "relative L2 < 0.005 and equilibrium residual < 1e-10",
                     "passed": rel_l2 < 0.005 and eq_resid < 1e-10})
    # Compare the 80-element result with the 160-element result at shared nodes.
    mesh_change = float(np.linalg.norm(solutions[160][::2] - solutions[80])
                        / np.linalg.norm(solutions[160][::2]))
    rows.append({"component": "mechanics", "study": "fine_mesh_change", "intervals": 160, "dt_s": 0.0,
                 "relative_mesh_change": mesh_change,
                 "criterion": "relative change between 80 and 160 elements < 0.01",
                 "passed": mesh_change < 0.01})
    return rows


def write_outputs(out_dir: Path) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = transport_verification() + mechanics_verification()
    csv_path = out_dir / "verification_results.csv"
    columns = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    transport_rows = [r for r in rows if r["component"] == "transport"]
    mech_fine = [r for r in rows if r["component"] == "mechanics"
                 and r["study"] == "heterogeneous_mesh_refinement" and r["intervals"] in (80, 160)]
    mesh_check = [r for r in rows if r.get("study") == "fine_mesh_change"]
    suite_passed = (all(bool(r["passed"]) for r in transport_rows)
                    and all(bool(r["passed"]) for r in mech_fine)
                    and len(mesh_check) == 1 and bool(mesh_check[0]["passed"]))
    summary = {
        "model_scope": "reduced 1-D numerical verification prototype; not fossil-specific",
        "transport_equation": "dC/dt = D*d2C/dx2 - k*C, C(0,t)=C(L,t)=Cb, C(x,0)=0 in the interior",
        "mechanics_equation": "d(sigma)/dx + b = 0, sigma=E(x)*du/dx, u(0)=0, sigma(L)=traction",
        "transport_criterion": "relative L2 error < 0.02 and relative maximum error < 0.03",
        "mechanics_criterion": "relative L2 error < 0.005 and normalized free-node equilibrium residual < 1e-10",
        "all_passed": suite_passed,
        "coarse_mechanics_levels_within_tolerance": all(bool(r["passed"]) for r in rows
            if r["component"] == "mechanics" and r["study"] == "heterogeneous_mesh_refinement" and r["intervals"] < 80),
        "rows": rows,
        "published_data_used": False,
        "fossil_parameter_identifiability_assessed": False
    }
    (out_dir / "verification_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    lines = [
        "# Reduced 1D model verification report", "",
        "This is a numerical qualification prototype. It is not calibrated to or presented as a reconstruction of the ichthyosaur fossil.", "",
        "## Transport benchmark", "",
        "The Crank–Nicolson solver was compared with the analytical eigenfunction solution for one-dimensional diffusion with first-order removal, zero initial interior concentration, and equal fixed boundary concentrations. The model uses a uniform grid and reports relative L2 and maximum errors. Benchmark values were L = 1 mm, D = 1 × 10⁻⁹ m² s⁻¹, k = 2 × 10⁻⁴ s⁻¹, Cb = 1 mol m⁻³, and t = 600 s; these are verification values, not fossil estimates.", "",
        "## Mechanics benchmark", "",
        "A linear finite-element bar with piecewise-constant Young's modulus was compared with the exact displacement obtained by integrating strain under uniform body-force density and end traction. The material interface falls inside elements on the tested meshes, and element stiffness integrates the modulus exactly across the interface. Parameters were L = 10 mm, A = 1 × 10⁻⁴ m², E = 1.2 and 2.4 GPa, body-force density = 2 × 10⁷ N m⁻³, and end traction = 0.5 MPa; these are verification values, not fossil estimates. The report includes relative displacement errors, normalized equilibrium residuals, and change between the two finest meshes.", "",
        "## Results", "",
        "| Component | Study | Elements/intervals | Time step (s) | Relative L2 error | Relative max error | Equilibrium residual | Mesh change | Meets criterion |",
        "|---|---|---:|---:|---:|---:|---:|---:|:---:|"
    ]
    for row in rows:
        def cell(key: str, precision: str = ".6g") -> str:
            value = row.get(key)
            return "n/a" if value is None else format(value, precision)
        lines.append(f"| {row['component']} | {row['study']} | {row['intervals']} | {cell('dt_s', '.4g')} | {cell('relative_L2_error')} | {cell('relative_max_error')} | {cell('relative_equilibrium_residual', '.3g')} | {cell('relative_mesh_change', '.3g')} | {'PASS' if row['passed'] else 'FAIL'} |")
    lines += ["", f"**Overall benchmark suite: {'PASS' if suite_passed else 'FAIL'}.** The mechanics suite passes only if the two finest levels meet the error criterion and their shared-node change is below 1%; coarser levels remain visible to show refinement behaviour.", ""]
    lines += ["", "## Interpretation and limits", "",
              "Passing these tests verifies only the implemented 1D transport and linear-elastic mechanics components for the stated benchmark equations. It does not validate the sulfur-cycle chemistry, mineral precipitation kinetics, porosity–permeability feedback, poroelastic coupling, or fossil-specific properties. No published fossil measurements were used in these tests. The separate IDENTIFIABILITY_SCREEN.md is a structural data audit, not a numerical identifiability result.", "",
              "The transport benchmark is a method-of-lines finite-difference solver with Crank–Nicolson time integration. The mechanics benchmark is a 1D linear finite-element bar. Both use SI units internally.", ""]
    (out_dir / "VERIFICATION_REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results", help="directory for CSV and report outputs")
    args = parser.parse_args()
    summary = write_outputs(Path(args.out))
    print(json.dumps({"all_passed": summary["all_passed"],
                      "transport": [r for r in summary["rows"] if r["component"] == "transport"],
                      "mechanics": [r for r in summary["rows"] if r["component"] == "mechanics"]}, indent=2))
    if not summary["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
