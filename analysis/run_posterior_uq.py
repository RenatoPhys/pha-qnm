#!/usr/bin/env python3
"""Controlled posterior propagation for critical and spinodal observables.

The 1600 HDF5 draws are already posterior samples, so their base weights are
equal.  We compress the 1589 successful-CEP draws to deterministic medoids in
the joint parameter/CEP space and assign each medoid the population of its
Voronoi cell.  The MAP draw is forced into the controlled ensemble.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import json
import math
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

import charged_hydrodynamics as hydro
from pha_qnm.thermodynamics import (
    BackgroundNumerics,
    PHAParameters,
    PHAModel,
    locate_cusp_critical_point,
    solve_background,
)


ROOT = Path(__file__).resolve().parents[1]
HDF5 = ROOT / "data" / "raw" / "Bayesian_polyhyper_muses.hdf5"
RESULTS = ROOT / "results" / "python"
PARAMETERS = ("kappa2", "gamma", "b2", "b4", "b6",
              "c1", "c2", "c3", "d1", "d2")
CRITICAL_OFFSETS = (0.50, 0.35, 0.25, 0.18, 0.12, 0.08, 0.05, 0.03,
                    0.02, 0.012, 0.008, 0.005)


def load_successful_draws() -> list[dict[str, Any]]:
    rows = []
    with h5py.File(HDF5, "r") as handle:
        for name, sample in handle["posterior_samples"].items():
            critical = sample["critical_point"].attrs
            if critical.get("status", "failed") != "success" or "Tc_in_MeV" not in critical:
                continue
            parameters = sample["parameters"].attrs
            row = {
                "sample": name,
                "sample_index": int(name.removeprefix("sample")),
                "log_likelihood": float(sample.attrs["log_likelihood"]),
                "metadata_T_c_MeV": float(critical["Tc_in_MeV"]),
                "metadata_mu_c_MeV": float(critical["muc_in_MeV"]),
                "Lambda_MeV": float(parameters["Lambda"]),
            }
            row.update({key: float(parameters[key]) for key in PARAMETERS})
            rows.append(row)
    return sorted(rows, key=lambda row: row["sample_index"])


def select_medoids(draws: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
    if not 2 <= count <= len(draws):
        raise ValueError("posterior sample count is outside the available range")
    features = np.asarray([
        [row[key] for key in ("metadata_T_c_MeV", "metadata_mu_c_MeV",
                              "Lambda_MeV", *PARAMETERS)]
        for row in draws
    ])
    center = np.median(features, axis=0)
    scale = np.subtract(*np.percentile(features, [75.0, 25.0], axis=0))
    scale = np.where(scale > 1.0e-14, scale, np.std(features, axis=0))
    standardized = (features - center) / scale
    _u, singular, vh = np.linalg.svd(standardized, full_matrices=False)
    cumulative = np.cumsum(singular**2) / np.sum(singular**2)
    dimensions = max(2, int(np.searchsorted(cumulative, 0.95) + 1))
    reduced = standardized @ vh[:dimensions].T
    centroids, labels = kmeans2(reduced, count, minit="++", iter=100, seed=271828)
    selected_indices = []
    for cluster in range(count):
        members = np.flatnonzero(labels == cluster)
        if members.size == 0:
            raise RuntimeError("empty posterior cluster")
        distances = np.linalg.norm(reduced[members] - centroids[cluster], axis=1)
        selected_indices.append(int(members[np.argmin(distances)]))

    map_index = next(i for i, row in enumerate(draws) if row["sample_index"] == 74)
    if map_index not in selected_indices:
        map_cluster = int(labels[map_index])
        selected_indices[map_cluster] = map_index
    medoid_features = reduced[selected_indices]
    assignment = np.argmin(
        np.linalg.norm(reduced[:, None, :] - medoid_features[None, :, :], axis=2),
        axis=1,
    )
    population = np.bincount(assignment, minlength=count)
    selected = []
    for medoid, size in zip(selected_indices, population):
        row = dict(draws[medoid])
        row["posterior_weight"] = float(size / len(draws))
        row["represented_draws"] = int(size)
        row["selection_pca_dimensions"] = dimensions
        selected.append(row)
    return sorted(selected, key=lambda row: row["sample_index"])


def parameters_from_row(row: dict[str, Any]) -> PHAParameters:
    return PHAParameters(
        Lambda_MeV=row["Lambda_MeV"],
        **{key: row[key] for key in PARAMETERS},
    )


def _charge_at_mu(phi0: float, target_mu: float, model: PHAModel,
                  numerics: BackgroundNumerics,
                  preferred_charge: float | None = None) -> tuple[float, dict[str, Any]]:
    cache: dict[float, dict[str, Any]] = {}

    def state(charge: float) -> dict[str, Any]:
        key = float(charge)
        if key not in cache:
            cache[key] = solve_background(phi0, key, model=model, numerics=numerics)
        return cache[key]

    grid = np.linspace(2.0e-4, 0.78, 33)
    values = []
    for charge in grid:
        try:
            values.append(state(float(charge))["mu_MeV"] - target_mu)
        except (RuntimeError, ValueError, OverflowError):
            values.append(float("nan"))
    brackets = []
    for left, right, fleft, fright in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if np.isfinite(fleft) and np.isfinite(fright) and fleft * fright <= 0.0:
            brackets.append((float(left), float(right)))
    if not brackets:
        raise RuntimeError(f"cannot bracket mu={target_mu:.3f} MeV at phi0={phi0:.5f}")
    if preferred_charge is not None:
        brackets.sort(key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - preferred_charge))
    charge = brentq(lambda value: state(value)["mu_MeV"] - target_mu,
                    *brackets[0], xtol=2.0e-10, rtol=2.0e-10)
    return float(charge), state(float(charge))


def _best_cusp(model: PHAModel) -> dict[str, Any]:
    numerics = BackgroundNumerics(phiA_tolerance=1.0e-5, ricci_tolerance=1.0e-5)
    candidates = []
    accepted = False
    for phi_step, charge_step in ((1.0e-2, 1.0e-3),
                                  (1.5e-2, 1.5e-3),
                                  (2.0e-2, 2.0e-3),
                                  (3.0e-2, 3.0e-3)):
        guesses = [(3.53, 0.346), (3.0, 0.42), (4.1, 0.27)]
        if candidates:
            best_point = min(candidates, key=lambda item: item[0])[1]
            guesses.insert(0, (best_point["phi0"], best_point["charge_fraction"]))
        for guess in guesses:
            try:
                result = locate_cusp_critical_point(
                    initial_guess=guess, model=model, background_numerics=numerics,
                    phi_step=phi_step, charge_step=charge_step,
                )
                residual = float(np.linalg.norm(result["cusp_equation_residual"]))
                candidates.append((residual, result))
                if residual < 2.0e-4:
                    accepted = True
                    break
            except (RuntimeError, ValueError, OverflowError, np.linalg.LinAlgError):
                continue
        if accepted:
            break
    if not candidates:
        raise RuntimeError("all independent cusp searches failed")
    residual, result = min(candidates, key=lambda item: item[0])
    result = dict(result)
    result["posterior_acceptance_residual"] = residual
    result["status"] = "success" if residual < 1.0e-3 else "failed"
    return result


def _critical_path(cusp: dict[str, Any], model: PHAModel) -> list[dict[str, float]]:
    numerics = hydro.TIGHT_BACKGROUND_NUMERICS
    phi_c = float(cusp["phi0"])
    mu_c = float(cusp["mu_B_c_MeV"])
    T_c = float(cusp["T_c_MeV"])
    preferred = float(cusp["charge_fraction"])
    points = []
    for offset in CRITICAL_OFFSETS:
        phi0 = phi_c - offset
        charge, equilibrium = _charge_at_mu(
            phi0, mu_c, model, numerics, preferred_charge=preferred
        )
        state = hydro.charged_hydro_state(
            phi0, charge, phi_step=1.0e-3, charge_step=1.0e-4,
            backend="reference_tight", model=model,
        )
        distance = math.hypot(
            equilibrium["T_MeV"] / T_c - 1.0,
            equilibrium["mu_MeV"] / mu_c - 1.0,
        )
        points.append({
            "phi_offset": offset,
            "phi0": phi0,
            "charge_fraction": charge,
            "T_MeV": equilibrium["T_MeV"],
            "mu_MeV": equilibrium["mu_MeV"],
            "reduced_distance": distance,
            "D_B_hat": hydro.baryon_diffusion_hat(state),
            "chi_B_over_T2": state.chi_B_over_T2,
            "thermo_condition": state.thermo_jacobian_condition,
            "maxwell_error": state.maxwell_relation_relative_error,
        })
        preferred = charge
    usable = [point for point in points
              if point["D_B_hat"] > 0.0 and point["chi_B_over_T2"] > 0.0
              and point["reduced_distance"] > 0.0]
    if len(usable) < 8:
        raise RuntimeError("insufficient positive critical-path points")
    windows = []
    for count in range(4, 9):
        fit_points = usable[-count:]
        log_r = np.log([point["reduced_distance"] for point in fit_points])
        log_D = np.log([point["D_B_hat"] for point in fit_points])
        log_chi = np.log([point["chi_B_over_T2"] for point in fit_points])
        D_r_slope, D_r_intercept = np.polyfit(log_r, log_D, 1)
        chi_r_slope, chi_r_intercept = np.polyfit(log_r, log_chi, 1)
        D_chi_slope, D_chi_intercept = np.polyfit(log_chi, log_D, 1)
        windows.append({
            "fit_point_count": count,
            "D_vs_distance_exponent": float(D_r_slope),
            "chi_vs_distance_exponent": float(chi_r_slope),
            "D_vs_chi_exponent": float(D_chi_slope),
            "z_eta0": float(2.0 - 2.0 * D_chi_slope),
            "D_vs_distance_intercept": float(D_r_intercept),
            "chi_vs_distance_intercept": float(chi_r_intercept),
            "D_vs_chi_intercept": float(D_chi_intercept),
        })
    central = {key: float(np.median([window[key] for window in windows]))
               for key in ("D_vs_distance_exponent", "chi_vs_distance_exponent",
                           "D_vs_chi_exponent", "z_eta0")}
    return points, {
        **central,
        "fit_windows": windows,
        "z_window_half_range": float(
            0.5 * np.ptp([window["z_eta0"] for window in windows])
        ),
    }


def refine_critical_sample(row: dict[str, Any]) -> dict[str, Any]:
    if row.get("status") != "success":
        return row
    try:
        model = PHAModel(parameters_from_row(row))
        cusp = {
            "phi0": row["cusp"]["phi0"],
            "charge_fraction": row["cusp"]["charge_fraction"],
            "T_c_MeV": row["cusp"]["T_c_MeV"],
            "mu_B_c_MeV": row["cusp"]["mu_B_c_MeV"],
        }
        path, fits = _critical_path(cusp, model)
        updated = dict(row)
        updated["critical_path"] = path
        updated["critical_fit"] = fits
        return updated
    except Exception as error:
        updated = dict(row)
        updated["status"] = "failed"
        updated["error"] = f"critical refinement: {type(error).__name__}: {error}"
        return updated


def _spinodal_line(cusp: dict[str, Any], model: PHAModel) -> dict[str, Any]:
    numerics = hydro.TIGHT_BACKGROUND_NUMERICS
    mu_target = float(cusp["mu_B_c_MeV"] + 50.0)
    phi_c = float(cusp["phi0"])
    preferred = float(cusp["charge_fraction"])
    rows = []
    for phi0 in np.linspace(phi_c - 0.75, phi_c + 0.95, 35):
        try:
            charge, state = _charge_at_mu(
                float(phi0), mu_target, model, numerics, preferred_charge=preferred
            )
        except RuntimeError:
            continue
        rows.append({"phi0": float(phi0), "charge_fraction": charge,
                     "T_MeV": float(state["T_MeV"]),
                     "mu_MeV": float(state["mu_MeV"])})
        preferred = charge
    if len(rows) < 15:
        raise RuntimeError("insufficient constant-mu points for spinodal search")
    phi = np.asarray([row["phi0"] for row in rows])
    temperature = np.asarray([row["T_MeV"] for row in rows])
    derivative = np.gradient(temperature, phi, edge_order=2)
    changes = np.flatnonzero(derivative[:-1] * derivative[1:] < 0.0)
    if changes.size < 2:
        raise RuntimeError("two spinodal folds were not resolved")
    derivative_spline = CubicSpline(phi, derivative)
    T_spline = CubicSpline(phi, temperature)
    q_spline = CubicSpline(phi, [row["charge_fraction"] for row in rows])
    folds = [brentq(derivative_spline, phi[index], phi[index + 1])
             for index in changes[:2]]
    midpoint = 0.5 * sum(folds)
    midpoint_charge, _ = _charge_at_mu(
        midpoint, mu_target, model, numerics,
        preferred_charge=float(q_spline(midpoint)),
    )
    midpoint_hydro = hydro.charged_hydro_state(
        midpoint, midpoint_charge, phi_step=5.0e-3, charge_step=5.0e-4,
        backend="reference_tight", model=model,
    )
    return {
        "mu_MeV": mu_target,
        "phi0_fold_hot": float(folds[0]),
        "phi0_fold_cold": float(folds[1]),
        "T_fold_low_MeV": float(min(T_spline(folds[0]), T_spline(folds[1]))),
        "T_fold_high_MeV": float(max(T_spline(folds[0]), T_spline(folds[1]))),
        "spinodal_temperature_width_MeV": float(abs(T_spline(folds[1]) - T_spline(folds[0]))),
        "midpoint_phi0": float(midpoint),
        "midpoint_charge_fraction": float(midpoint_charge),
        "midpoint_D_B_hat": hydro.baryon_diffusion_hat(midpoint_hydro),
        "grid_point_count": len(rows),
    }


def process_sample(row: dict[str, Any]) -> dict[str, Any]:
    result = {key: row[key] for key in row}
    try:
        model = PHAModel(parameters_from_row(row))
        cusp = _best_cusp(model)
        if cusp["status"] != "success":
            raise RuntimeError(
                f"cusp residual {cusp['posterior_acceptance_residual']:.3e}"
            )
        path, fits = _critical_path(cusp, model)
        spinodal = _spinodal_line(cusp, model)
        result.update({
            "status": "success",
            "cusp": {
                "phi0": cusp["phi0"],
                "charge_fraction": cusp["charge_fraction"],
                "T_c_MeV": cusp["T_c_MeV"],
                "mu_B_c_MeV": cusp["mu_B_c_MeV"],
                "equation_residual": cusp["posterior_acceptance_residual"],
                "metadata_delta_T_MeV": (
                    cusp["T_c_MeV"] - row["metadata_T_c_MeV"]
                ),
                "metadata_delta_mu_MeV": (
                    cusp["mu_B_c_MeV"] - row["metadata_mu_c_MeV"]
                ),
            },
            "critical_path": path,
            "critical_fit": fits,
            "spinodal": spinodal,
        })
    except Exception as error:  # preserve failed samples in the audit artifact
        result.update({"status": "failed", "error": f"{type(error).__name__}: {error}"})
    return result


def weighted_quantile(values: np.ndarray, weights: np.ndarray,
                      quantiles: tuple[float, ...]) -> list[float]:
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative /= np.sum(weights)
    return [float(np.interp(quantile, cumulative, values)) for quantile in quantiles]


def summarize(results: list[dict[str, Any]]) -> dict[str, Any]:
    successful = [row for row in results if row["status"] == "success"]
    failed = [row for row in results if row["status"] != "success"]
    observables = {
        "T_c_MeV": lambda row: row["cusp"]["T_c_MeV"],
        "mu_c_MeV": lambda row: row["cusp"]["mu_B_c_MeV"],
        "z_eta0": lambda row: row["critical_fit"]["z_eta0"],
        "D_vs_distance_exponent": lambda row: row["critical_fit"]["D_vs_distance_exponent"],
        "spinodal_temperature_width_MeV": lambda row: row["spinodal"]["spinodal_temperature_width_MeV"],
        "spinodal_midpoint_D_B_hat": lambda row: row["spinodal"]["midpoint_D_B_hat"],
    }
    intervals = {}
    for name, getter in observables.items():
        values = np.asarray([getter(row) for row in successful], dtype=float)
        weights = np.asarray([row["posterior_weight"] for row in successful], dtype=float)
        intervals[name] = dict(zip(
            ("q2p5", "median", "q97p5"),
            weighted_quantile(values, weights, (0.025, 0.5, 0.975)),
        ))
    return {
        "available_successful_hdf5_draws": 1589,
        "selected_medoids": len(results),
        "successful_medoids": len(successful),
        "represented_posterior_weight": float(sum(
            row["posterior_weight"] for row in successful
        )),
        "failed_samples": [row["sample"] for row in failed],
        "credible_intervals": intervals,
        "spinodal_diffusion_unstable_fraction": float(sum(
            row["posterior_weight"] for row in successful
            if row["spinodal"]["midpoint_D_B_hat"] < 0.0
        ) / max(sum(row["posterior_weight"] for row in successful), 1.0e-300)),
        "selection": (
            "deterministic weighted medoids of the successful HDF5 posterior "
            "in 95%-variance PCA parameter/CEP space, with sample74 forced in"
        ),
    }


def write_selection(rows: list[dict[str, Any]]) -> None:
    path = RESULTS / "posterior_selection.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None,
                        help="debug only: process the first N selected medoids")
    parser.add_argument("--select-only", action="store_true")
    parser.add_argument("--refine-existing", action="store_true",
                        help="recompute only the asymptotic critical paths in the saved ensemble")
    parser.add_argument("--retry-failed", action="store_true",
                        help="rerun only failed medoids with the adaptive cusp search")
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    if args.refine_existing:
        path = RESULTS / "posterior_uq_samples.json"
        existing = json.loads(path.read_text(encoding="utf-8"))
        if args.workers == 1:
            completed = [refine_critical_sample(row) for row in existing]
        else:
            completed = []
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                futures = {executor.submit(refine_critical_sample, row): row
                           for row in existing}
                for index, future in enumerate(as_completed(futures), 1):
                    result = future.result(); completed.append(result)
                    print(f"refine {index}/{len(existing)} {result['sample']} "
                          f"{result['status']}", flush=True)
        completed.sort(key=lambda row: row["sample_index"])
        path.write_text(json.dumps(completed, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
        summary = summarize(completed)
        (RESULTS / "posterior_uq_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return
    if args.retry_failed:
        path = RESULTS / "posterior_uq_samples.json"
        existing = json.loads(path.read_text(encoding="utf-8"))
        failures = [row for row in existing if row.get("status") != "success"]
        retained = [row for row in existing if row.get("status") == "success"]
        retried = []
        if args.workers == 1:
            retried = [process_sample(row) for row in failures]
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                futures = {executor.submit(process_sample, row): row for row in failures}
                for index, future in enumerate(as_completed(futures), 1):
                    result = future.result(); retried.append(result)
                    print(f"retry {index}/{len(failures)} {result['sample']} "
                          f"{result['status']}", flush=True)
        completed = sorted(retained + retried, key=lambda row: row["sample_index"])
        path.write_text(json.dumps(completed, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
        summary = summarize(completed)
        (RESULTS / "posterior_uq_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return
    selected = select_medoids(load_successful_draws(), args.samples)
    write_selection(selected)
    if args.select_only:
        print(json.dumps({"selected": len(selected),
                          "posterior_weight": sum(row["posterior_weight"] for row in selected)},
                         indent=2))
        return
    if args.limit is not None:
        selected = selected[:args.limit]
    completed = []
    if args.workers == 1:
        for index, row in enumerate(selected, 1):
            completed.append(process_sample(row))
            print(f"posterior {index}/{len(selected)} {row['sample']} "
                  f"{completed[-1]['status']}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(process_sample, row): row for row in selected}
            for index, future in enumerate(as_completed(futures), 1):
                result = future.result()
                completed.append(result)
                print(f"posterior {index}/{len(selected)} {result['sample']} "
                      f"{result['status']}", flush=True)
    completed.sort(key=lambda row: row["sample_index"])
    (RESULTS / "posterior_uq_samples.json").write_text(
        json.dumps(completed, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = summarize(completed)
    (RESULTS / "posterior_uq_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["successful_medoids"] < max(20, math.ceil(0.8 * len(completed))):
        raise SystemExit("fewer than 20 controlled posterior medoids passed")


if __name__ == "__main__":
    main()
