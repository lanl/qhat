#!/usr/bin/env python3
"""LiH single-case Trotter step/time sweep.

Place this file in ``analysis/`` next to ``benchmark_L_sweep_trotter.py``.
It deliberately reuses the benchmark module's QHAT/OpenFermion loading,
Jordan-Wigner ordering, Pauli matrices, and first/second-order product formulas.

The default case is the small LiH configuration

    Li-H_1.32_sto-6g_as-002-002.tensors.npz

produced from

    Li-H_1.32_sto-6g_as-002-002_JW.config

by ``hamiltonian_generator/build_L_sweep_tensors.py``.

The script isolates the Trotter-step effect by using the raw JW insertion order
and varying only total evolution time T, step count r, and therefore dt=T/r.
It writes one CSV and separate PNGs for operator error, state-vector error,
phase-aligned state error, the overlap error 1-|<psi_exact|psi_Trot>|,
and explicit repeated-step runtime.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import math
import statistics
import time
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

import benchmark_L_sweep_trotter as baseline


DEFAULT_TENSOR = Path(
    "hamiltonian_generator/library/Li-H/1.32/sto-6g/"
    "Li-H_1.32_sto-6g_as-002-002.tensors.npz"
)
DEFAULT_STEPS = [1, 2, 4, 8, 16, 32, 64, 128, 256]
DEFAULT_TIMES = [1.0]
DEFAULT_OUTPUT_DIR = Path("analysis/lih_trotter_step_case_study_overlap")

CSV_FIELDS = [
    "case_id",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "number_of_pauli_terms",
    "formula_order",
    "evolution_time",
    "trotter_steps",
    "trotter_dt",
    "nominal_exponential_count",
    "exact_evolution_time_seconds",
    "explicit_state_trotter_runtime_seconds",
    "unitary_matrix_power_runtime_seconds",
    "operator_norm_error",
    "state_vector_2norm_error",
    "state_overlap_abs",
    "one_minus_abs_overlap",
    "phase_aligned_state_2norm_error",
    "operator_local_convergence_order",
    "state_vector_local_convergence_order",
    "overlap_error_local_convergence_order",
    "phase_aligned_local_convergence_order",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Single LiH case study: sweep Trotter steps r and total time T, "
            "then plot cost and error."
        )
    )
    parser.add_argument(
        "--tensor",
        type=Path,
        default=DEFAULT_TENSOR,
        help="Input *.tensors.npz file.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=DEFAULT_STEPS,
        help="Trotter step counts r.",
    )
    parser.add_argument(
        "--times",
        type=float,
        nargs="+",
        default=DEFAULT_TIMES,
        help="Total evolution times T.",
    )
    parser.add_argument(
        "--formula-order",
        type=int,
        choices=[1, 2],
        default=1,
        help="Product-formula order. Use 1 for the main step-count case study.",
    )
    parser.add_argument(
        "--runtime-repeats",
        type=int,
        default=7,
        help=(
            "Repeat explicit state stepping this many times and report the "
            "median classical runtime."
        ),
    )
    parser.add_argument(
        "--fit-points",
        type=int,
        default=4,
        help="Number of largest-r points used for the printed asymptotic fit.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for CSV and PNG outputs.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=baseline.DEFAULT_TOLERANCE,
        help="Coefficient/numerical tolerance.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.tensor.is_file():
        raise FileNotFoundError(
            f"Tensor file not found: {args.tensor}\n"
            "Generate the LiH tensor first with build_L_sweep_tensors.py."
        )
    if any(step <= 0 for step in args.steps):
        raise ValueError("Every Trotter step count must be positive.")
    if any(total_time <= 0.0 for total_time in args.times):
        raise ValueError("Every total evolution time must be positive.")
    if args.runtime_repeats <= 0:
        raise ValueError("--runtime-repeats must be positive.")
    if args.fit_points < 2:
        raise ValueError("--fit-points must be at least 2.")


def prepare_case(
    tensor_path: Path,
    tolerance: float,
) -> tuple[
    baseline.CaseMetadata,
    np.ndarray,
    np.ndarray,
    list[tuple[float, np.ndarray]],
    int,
]:
    """Build the same raw-JW dense objects used by the repository benchmark."""
    interaction, n_qubits = baseline.load_interaction_operator(tensor_path)
    metadata = baseline.parse_case_metadata(tensor_path, n_qubits)

    fermion_hamiltonian = baseline.clean_fermion_operator(
        baseline.get_fermion_operator(interaction),
        tolerance,
    )
    full_jw_hamiltonian = baseline.jordan_wigner(fermion_hamiltonian)
    full_jw_hamiltonian.compress(abs_tol=tolerance)

    # ``build_orderings`` gained an ``n_qubits`` argument on newer L-sweep
    # revisions.  Detect the local repository signature so this case-study
    # driver works with both the older and newer benchmark implementations.
    build_orderings_kwargs = {
        "fermion_hamiltonian": fermion_hamiltonian,
        "full_jw_hamiltonian": full_jw_hamiltonian,
        "requested_orderings": ["jw_raw"],
        "tolerance": tolerance,
    }
    if "n_qubits" in inspect.signature(baseline.build_orderings).parameters:
        build_orderings_kwargs["n_qubits"] = n_qubits

    orderings, _, final_coefficients = baseline.build_orderings(
        **build_orderings_kwargs
    )

    raw_pauli_keys = orderings["jw_raw"].pauli_keys
    matrix_cache = {
        key: baseline.pauli_matrix_from_key(key, n_qubits)
        for key in raw_pauli_keys
    }

    dimension = 2**n_qubits
    identity = np.eye(dimension, dtype=complex)
    exact_hamiltonian = np.zeros((dimension, dimension), dtype=complex)

    ordered_terms: list[tuple[float, np.ndarray]] = []
    for pauli_key in raw_pauli_keys:
        coefficient = baseline.real_coefficient(
            final_coefficients[pauli_key],
            tolerance,
        )
        pauli_matrix = matrix_cache[pauli_key]
        exact_hamiltonian += coefficient * pauli_matrix
        ordered_terms.append((coefficient, pauli_matrix))

    initial_state = baseline.build_hartree_fock_state(
        n_qubits,
        metadata.active_occupied,
    )

    return (
        metadata,
        exact_hamiltonian,
        initial_state,
        ordered_terms,
        len(raw_pauli_keys),
    )


def build_one_slice(
    formula_order: int,
    ordered_terms: Sequence[tuple[float, np.ndarray]],
    dt: float,
    identity: np.ndarray,
) -> tuple[np.ndarray, int]:
    if formula_order == 1:
        one_slice = baseline.first_order_slice(ordered_terms, dt, identity)
        exponentials_per_step = len(ordered_terms)
    elif formula_order == 2:
        one_slice = baseline.second_order_slice(ordered_terms, dt, identity)
        exponentials_per_step = max(0, 2 * len(ordered_terms) - 1)
    else:
        raise ValueError(f"Unsupported formula order: {formula_order}")
    return one_slice, exponentials_per_step


def explicit_state_runtime(
    one_slice: np.ndarray,
    initial_state: np.ndarray,
    trotter_steps: int,
    repeats: int,
) -> tuple[float, np.ndarray]:
    """Time actual repeated step application instead of matrix_power."""
    samples: list[float] = []
    final_state = initial_state.copy()

    for _ in range(repeats):
        state = initial_state.copy()
        start = time.perf_counter()
        for _step in range(trotter_steps):
            state = one_slice @ state
        samples.append(time.perf_counter() - start)
        final_state = state

    return statistics.median(samples), final_state


def local_order(
    previous_r: int,
    previous_error: float,
    current_r: int,
    current_error: float,
) -> float:
    if previous_error <= 0.0 or current_error <= 0.0:
        return math.nan
    if current_r == previous_r:
        return math.nan
    return math.log(previous_error / current_error) / math.log(
        current_r / previous_r
    )


def run_sweep(args: argparse.Namespace) -> list[dict[str, float | int | str]]:
    (
        metadata,
        exact_hamiltonian,
        initial_state,
        ordered_terms,
        number_of_pauli_terms,
    ) = prepare_case(args.tensor, args.tolerance)

    dimension = exact_hamiltonian.shape[0]
    identity = np.eye(dimension, dtype=complex)
    rows: list[dict[str, float | int | str]] = []

    steps = sorted(set(args.steps))
    times = sorted(set(args.times))

    print(f"Case:               {metadata.case_id}")
    print(f"Qubits:             {metadata.n_qubits}")
    print(f"Raw JW Pauli terms: {number_of_pauli_terms}")
    print(f"Formula order:      {args.formula_order}")
    print(f"Total times T:      {times}")
    print(f"Trotter steps r:    {steps}")
    print()

    for total_time in times:
        exact_start = time.perf_counter()
        exact_unitary = baseline.expm(-1j * total_time * exact_hamiltonian)
        exact_evolution_runtime = time.perf_counter() - exact_start
        exact_state = exact_unitary @ initial_state

        previous: dict[str, tuple[int, float]] = {}

        for trotter_steps in steps:
            dt = total_time / trotter_steps
            one_slice, exponentials_per_step = build_one_slice(
                args.formula_order,
                ordered_terms,
                dt,
                identity,
            )

            state_runtime, approximate_state = explicit_state_runtime(
                one_slice,
                initial_state,
                trotter_steps,
                args.runtime_repeats,
            )

            matrix_power_start = time.perf_counter()
            trotter_unitary = np.linalg.matrix_power(
                one_slice,
                trotter_steps,
            )
            matrix_power_runtime = time.perf_counter() - matrix_power_start

            operator_error = baseline.spectral_norm(
                trotter_unitary - exact_unitary
            )
            state_vector_error = float(
                np.linalg.norm(exact_state - approximate_state, ord=2)
            )

            # Phase-insensitive state comparison.  Normalize before taking the
            # overlap so tiny floating-point norm drift cannot make |overlap| > 1.
            exact_norm = float(np.linalg.norm(exact_state, ord=2))
            approximate_norm = float(np.linalg.norm(approximate_state, ord=2))
            if exact_norm == 0.0 or approximate_norm == 0.0:
                raise ValueError("Cannot compute overlap for a zero-norm state.")

            overlap_abs = float(
                abs(np.vdot(exact_state, approximate_state))
                / (exact_norm * approximate_norm)
            )
            overlap_abs = min(1.0, max(0.0, overlap_abs))

            # Requested metric: 1 - |<psi_exact|psi_Trot>|.
            one_minus_abs_overlap = max(0.0, 1.0 - overlap_abs)

            # For normalized states this is the minimum 2-norm distance after
            # removing an irrelevant global phase:
            #   min_phi ||psi_exact - exp(i phi) psi_Trot||_2
            #       = sqrt(2 - 2 |<psi_exact|psi_Trot>|).
            phase_aligned_state_error = math.sqrt(
                max(0.0, 2.0 - 2.0 * overlap_abs)
            )

            row: dict[str, float | int | str] = {
                "case_id": metadata.case_id,
                "molecule": metadata.molecule,
                "bond_length": metadata.bond_length,
                "basis": metadata.basis,
                "active_occupied": metadata.active_occupied,
                "active_vacant": metadata.active_vacant,
                "n_qubits": metadata.n_qubits,
                "number_of_pauli_terms": number_of_pauli_terms,
                "formula_order": args.formula_order,
                "evolution_time": total_time,
                "trotter_steps": trotter_steps,
                "trotter_dt": dt,
                "nominal_exponential_count": (
                    trotter_steps * exponentials_per_step
                ),
                "exact_evolution_time_seconds": exact_evolution_runtime,
                "explicit_state_trotter_runtime_seconds": state_runtime,
                "unitary_matrix_power_runtime_seconds": matrix_power_runtime,
                "operator_norm_error": operator_error,
                "state_vector_2norm_error": state_vector_error,
                "state_overlap_abs": overlap_abs,
                "one_minus_abs_overlap": one_minus_abs_overlap,
                "phase_aligned_state_2norm_error": phase_aligned_state_error,
                "operator_local_convergence_order": math.nan,
                "state_vector_local_convergence_order": math.nan,
                "overlap_error_local_convergence_order": math.nan,
                "phase_aligned_local_convergence_order": math.nan,
            }

            metric_map = {
                "operator": (
                    "operator_norm_error",
                    "operator_local_convergence_order",
                ),
                "state_vector": (
                    "state_vector_2norm_error",
                    "state_vector_local_convergence_order",
                ),
                "overlap_error": (
                    "one_minus_abs_overlap",
                    "overlap_error_local_convergence_order",
                ),
                "phase_aligned": (
                    "phase_aligned_state_2norm_error",
                    "phase_aligned_local_convergence_order",
                ),
            }

            for key, (error_column, order_column) in metric_map.items():
                current_error = float(row[error_column])
                if key in previous:
                    previous_r, previous_error = previous[key]
                    row[order_column] = local_order(
                        previous_r,
                        previous_error,
                        trotter_steps,
                        current_error,
                    )
                previous[key] = (trotter_steps, current_error)

            rows.append(row)

            print(
                f"T={total_time:g}  r={trotter_steps:4d}  "
                f"dt={dt:.6g}  exp={int(row['nominal_exponential_count']):6d}  "
                f"op={operator_error:.6e}  "
                f"state2={state_vector_error:.6e}  "
                f"1-|ov|={one_minus_abs_overlap:.6e}  "
                f"phase2={phase_aligned_state_error:.6e}  "
                f"runtime={state_runtime:.3e}s"
            )

    return rows


def write_csv(
    rows: Sequence[dict[str, float | int | str]],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def rows_for_time(
    rows: Sequence[dict[str, float | int | str]],
    total_time: float,
) -> list[dict[str, float | int | str]]:
    return sorted(
        [row for row in rows if float(row["evolution_time"]) == total_time],
        key=lambda row: int(row["trotter_steps"]),
    )


def save_loglog_plot(
    rows: Sequence[dict[str, float | int | str]],
    metric: str,
    ylabel: str,
    output_path: Path,
) -> None:
    plt.figure(figsize=(7.2, 5.0))
    for total_time in sorted({float(row["evolution_time"]) for row in rows}):
        subset = rows_for_time(rows, total_time)
        x = [int(row["trotter_steps"]) for row in subset]
        y = [max(float(row[metric]), np.finfo(float).tiny) for row in subset]
        plt.plot(x, y, marker="o", label=f"T={total_time:g}")

    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.xlabel("Trotter steps r")
    plt.ylabel(ylabel)
    plt.title("LiH Trotter error versus step count")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def save_runtime_plot(
    rows: Sequence[dict[str, float | int | str]],
    output_path: Path,
) -> None:
    plt.figure(figsize=(7.2, 5.0))
    for total_time in sorted({float(row["evolution_time"]) for row in rows}):
        subset = rows_for_time(rows, total_time)
        x = [int(row["trotter_steps"]) for row in subset]
        y = [
            max(
                float(row["explicit_state_trotter_runtime_seconds"]),
                np.finfo(float).tiny,
            )
            for row in subset
        ]
        plt.plot(x, y, marker="o", label=f"T={total_time:g}")

    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.xlabel("Trotter steps r")
    plt.ylabel("Median explicit state-stepping runtime [s]")
    plt.title("LiH classical repeated-step cost versus r")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def save_exponential_count_plot(
    rows: Sequence[dict[str, float | int | str]],
    output_path: Path,
) -> None:
    first_time = min(float(row["evolution_time"]) for row in rows)
    subset = rows_for_time(rows, first_time)
    x = [int(row["trotter_steps"]) for row in subset]
    y = [int(row["nominal_exponential_count"]) for row in subset]

    plt.figure(figsize=(7.2, 5.0))
    plt.plot(x, y, marker="o")
    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.xlabel("Trotter steps r")
    plt.ylabel("Nominal Pauli exponentials")
    plt.title("LiH Trotter cost proxy versus step count")
    plt.grid(True, which="both", alpha=0.25)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def fit_scaling(
    subset: Sequence[dict[str, float | int | str]],
    metric: str,
    fit_points: int,
) -> float:
    usable = [row for row in subset if float(row[metric]) > 0.0]
    usable = usable[-min(fit_points, len(usable)):]
    if len(usable) < 2:
        return math.nan

    log_r = np.log([float(row["trotter_steps"]) for row in usable])
    log_error = np.log([float(row[metric]) for row in usable])
    slope, _intercept = np.polyfit(log_r, log_error, 1)
    return float(-slope)


def print_scaling_summary(
    rows: Sequence[dict[str, float | int | str]],
    formula_order: int,
    fit_points: int,
) -> None:
    expected_state_order = float(formula_order)
    expected_overlap_error_order = float(2 * formula_order)

    print()
    print("Asymptotic fits from the largest-r points")
    print("error ~ r^(-p)")
    print(
        f"Expected generic p: operator/state ~= {expected_state_order:g}; "
        f"1-|overlap| ~= {expected_overlap_error_order:g}"
    )

    for total_time in sorted({float(row["evolution_time"]) for row in rows}):
        subset = rows_for_time(rows, total_time)
        operator_p = fit_scaling(
            subset,
            "operator_norm_error",
            fit_points,
        )
        state_p = fit_scaling(
            subset,
            "state_vector_2norm_error",
            fit_points,
        )
        overlap_error_p = fit_scaling(
            subset,
            "one_minus_abs_overlap",
            fit_points,
        )
        phase_aligned_p = fit_scaling(
            subset,
            "phase_aligned_state_2norm_error",
            fit_points,
        )
        print(
            f"T={total_time:g}: "
            f"p_operator={operator_p:.3f}, "
            f"p_state2={state_p:.3f}, "
            f"p_phase_aligned={phase_aligned_p:.3f}, "
            f"p_1-|overlap|={overlap_error_p:.3f}"
        )


def main() -> None:
    args = parse_args()
    validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = run_sweep(args)

    csv_path = args.output_dir / "lih_trotter_step_overlap_sweep.csv"
    write_csv(rows, csv_path)

    save_loglog_plot(
        rows,
        "operator_norm_error",
        r"Operator error $\|U_{\rm Trot}-U_{\rm exact}\|_2$",
        args.output_dir / "operator_error_vs_steps.png",
    )
    save_loglog_plot(
        rows,
        "state_vector_2norm_error",
        r"State-vector error $\|\psi_{\rm Trot}-\psi_{\rm exact}\|_2$",
        args.output_dir / "state_vector_error_vs_steps.png",
    )
    save_loglog_plot(
        rows,
        "one_minus_abs_overlap",
        r"State overlap error $1-|\langle\psi_{\rm exact}|\psi_{\rm Trot}\rangle|$",
        args.output_dir / "one_minus_abs_overlap_vs_steps.png",
    )
    save_loglog_plot(
        rows,
        "phase_aligned_state_2norm_error",
        r"Phase-aligned state error $\min_\phi\|\psi_{\rm exact}-e^{i\phi}\psi_{\rm Trot}\|_2$",
        args.output_dir / "phase_aligned_state_error_vs_steps.png",
    )
    save_runtime_plot(
        rows,
        args.output_dir / "runtime_vs_steps.png",
    )
    save_exponential_count_plot(
        rows,
        args.output_dir / "exponential_count_vs_steps.png",
    )

    print_scaling_summary(rows, args.formula_order, args.fit_points)

    print()
    print(f"CSV:     {csv_path}")
    print(f"Figures: {args.output_dir}")


if __name__ == "__main__":
    main()

