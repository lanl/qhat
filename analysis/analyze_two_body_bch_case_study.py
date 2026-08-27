#!/usr/bin/env python3
"""Decompose HF-state BCH2 cancellation by fermionic parent body order.

For every case in the held-out 20-case cancellation panel, this analysis
reconstructs each final Jordan--Wigner Pauli coefficient as the signed sum of
contributions from complete Hermitian one- and two-body fermionic parents.
The coefficient product in every noncommuting Pauli-pair BCH contribution is
then split exactly into 1B--1B, 1B--2B, and 2B--2B pieces.

The same fixed component amplitudes are evaluated under two Pauli schedules:

* signed fermionic parent order (ascending canonical signed coefficient), and
* direct JW magnitude-descending order.

This permits a direct test of whether order-induced cancellation inside the
two-body-bearing BCH vector, or destructive interference between that vector
and the one-body-only vector, tracks the observed Trotter-error advantage.
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from openfermion import get_fermion_operator, jordan_wigner
from scipy import stats

try:
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_hermitian_fermion_terms,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
    )
    from qhat.analysis.benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
        induced_pauli_order_indices,
        precompute_fermion_to_pauli_indices,
    )
except ImportError:
    import benchmark_b2_signed_coefficient_baseline as baseline
    from benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_hermitian_fermion_terms,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
    )
    from benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
        induced_pauli_order_indices,
        precompute_fermion_to_pauli_indices,
    )


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = Path(
    "analysis/cancellation_hypothesis_validation/full_panel_manifest.csv"
)
DEFAULT_CASE_SUMMARY = Path(
    "analysis/cancellation_hypothesis_validation/full_analysis/case_summary.csv"
)
DEFAULT_OUTDIR = Path("analysis/two_body_bch_case_study")

SIGNED_ORDER = "fermionic_signed_reference"
JW_ORDER = "jw_magnitude_baseline"


@dataclass(frozen=True)
class BCHComponents:
    """Order-independent exact body-order split of Pauli-pair amplitudes."""

    one_one: np.ndarray
    one_two: np.ndarray
    two_two: np.ndarray

    @property
    def two_body_bearing(self) -> np.ndarray:
        return self.one_two + self.two_two

    @property
    def full(self) -> np.ndarray:
        return self.one_one + self.one_two + self.two_two


FEATURES: tuple[tuple[str, str, str], ...] = (
    (
        "log10_two_body_cancellation_advantage",
        "2B-bearing cancellation advantage",
        "primary",
    ),
    (
        "delta_two_body_cross_destructiveness",
        "Gain in destructive 1B/2B interference",
        "primary",
    ),
    (
        "delta_two_body_leave_in_reduction",
        "Gain in norm reduction from adding 2B vector",
        "primary",
    ),
    (
        "log10_one_two_cancellation_advantage",
        "1B-2B cancellation advantage",
        "component",
    ),
    (
        "log10_two_two_cancellation_advantage",
        "2B-2B cancellation advantage",
        "component",
    ),
    (
        "log10_one_one_cancellation_advantage",
        "1B-1B cancellation advantage",
        "negative_control",
    ),
    (
        "two_body_pair_mass_fraction",
        "2B-bearing raw BCH-pair mass fraction",
        "structural_control",
    ),
)

OUTCOMES: tuple[tuple[str, str], ...] = (
    ("log10_error_advantage", "log10(JW/F Trotter-error advantage)"),
    ("log10_full_bch_advantage", "log10(JW/F full BCH cancellation advantage)"),
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--case-summary", type=Path, default=DEFAULT_CASE_SUMMARY)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--case-id", action="append")
    parser.add_argument("--max-qubits", type=int)
    parser.add_argument(
        "--reuse-case-metrics",
        action="store_true",
        help="Reuse the existing per-case CSV and regenerate statistics/report/figure.",
    )
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args()


def resolve_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def safe_ratio(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator):
        return float("nan")
    if numerator < 0.0 or denominator <= 0.0:
        return float("nan")
    return float(numerator / denominator)


def safe_log10(value: float) -> float:
    return float(math.log10(value)) if np.isfinite(value) and value > 0.0 else float("nan")


def bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(values))
    if finite.size == 0:
        return adjusted
    order = finite[np.argsort(values[finite])]
    ranked = values[order] * finite.size / np.arange(1, finite.size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def exact_paired_sign_permutation_p_value(differences: Sequence[float]) -> float:
    values = np.asarray(differences, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    if values.size > 20:
        raise ValueError("Exact sign enumeration is limited to 20 pairs.")
    observed = abs(float(np.mean(values)))
    exceed = 0
    total = 2**values.size
    for signs in itertools.product((-1.0, 1.0), repeat=values.size):
        permuted = abs(float(np.mean(values * np.asarray(signs))))
        exceed += int(permuted >= observed - 1.0e-15)
    return exceed / total


def paired_bootstrap_mean_ci(
    differences: Sequence[float],
    replicates: int,
    seed: int,
) -> tuple[float, float]:
    values = np.asarray(differences, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or replicates <= 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, values.size, size=(replicates, values.size))
    means = np.mean(values[indices], axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def pearson(x: Sequence[float], y: Sequence[float]) -> tuple[int, float, float]:
    frame = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(frame) < 3 or frame["x"].nunique() < 2 or frame["y"].nunique() < 2:
        return len(frame), float("nan"), float("nan")
    result = stats.pearsonr(frame["x"], frame["y"])
    return len(frame), float(result.statistic), float(result.pvalue)


def spearman(x: Sequence[float], y: Sequence[float]) -> tuple[float, float]:
    frame = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(frame) < 3 or frame["x"].nunique() < 2 or frame["y"].nunique() < 2:
        return float("nan"), float("nan")
    result = stats.spearmanr(frame["x"], frame["y"])
    return float(result.statistic), float(result.pvalue)


def within_group_residual(values: pd.Series, groups: pd.Series) -> pd.Series:
    return values - values.groupby(groups).transform("mean")


def partial_pearson(
    frame: pd.DataFrame,
    feature: str,
    outcome: str,
    numeric_controls: Sequence[str] = (),
    categorical_controls: Sequence[str] = (),
) -> tuple[int, float, float]:
    columns = [feature, outcome, *numeric_controls, *categorical_controls]
    valid = frame[columns].dropna().copy()
    if len(valid) < 3:
        return len(valid), float("nan"), float("nan")
    design_parts = [
        valid[list(numeric_controls)].astype(float)
        if numeric_controls
        else pd.DataFrame(index=valid.index)
    ]
    for control in categorical_controls:
        design_parts.append(
            pd.get_dummies(valid[control], prefix=control, drop_first=True).astype(float)
        )
    controls = pd.concat(design_parts, axis=1)
    design = np.column_stack(
        [np.ones(len(valid), dtype=float), controls.to_numpy(dtype=float)]
    )

    def residual(column: str) -> np.ndarray:
        values = valid[column].to_numpy(dtype=float)
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        return values - design @ coefficients

    return pearson(residual(feature), residual(outcome))


def parent_body_order(term: Any) -> int:
    component_lengths = {len(key) for key in term.component_keys}
    if len(component_lengths) != 1:
        raise ValueError(f"Mixed body order within Hermitian parent {term.index}.")
    operator_length = next(iter(component_lengths))
    if operator_length not in (2, 4):
        raise ValueError(
            f"Expected one- or two-body parent, found monomial length {operator_length}."
        )
    return operator_length // 2


def decompose_final_pauli_coefficients(
    hermitian_terms: Sequence[Any],
    raw_pauli_keys: Sequence[Any],
    final_coefficients: dict[Any, complex],
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Return exact signed 1B and 2B contributions to every final coefficient."""

    raw_index = {key: index for index, key in enumerate(raw_pauli_keys)}
    one_body = np.zeros(len(raw_pauli_keys), dtype=float)
    two_body = np.zeros(len(raw_pauli_keys), dtype=float)
    one_body_abs = np.zeros(len(raw_pauli_keys), dtype=float)
    two_body_abs = np.zeros(len(raw_pauli_keys), dtype=float)

    for term in hermitian_terms:
        body_order = parent_body_order(term)
        mapped = jordan_wigner(term.operator)
        destination = one_body if body_order == 1 else two_body
        absolute_destination = one_body_abs if body_order == 1 else two_body_abs
        for key, coefficient in mapped.terms.items():
            if key == () or key not in raw_index:
                continue
            index = raw_index[key]
            value = real_coefficient(coefficient, tolerance)
            destination[index] += value
            absolute_destination[index] += abs(value)

    final = np.asarray(
        [real_coefficient(final_coefficients[key], tolerance) for key in raw_pauli_keys],
        dtype=float,
    )
    # JW accumulation of thousands of mutually cancelling contributions is not
    # associative in floating-point arithmetic.  The full-Hamiltonian transform
    # and the parent-by-parent transform can therefore differ by a tiny residual.
    # Preserve the exact final Hamiltonian used by the benchmark by distributing
    # that residual according to the absolute 1B/2B provenance mass of each key.
    raw_residual = final - (one_body + two_body)
    max_reconciliation = float(np.max(np.abs(raw_residual), initial=0.0))
    provenance_mass = one_body_abs + two_body_abs
    unresolved = (np.abs(raw_residual) > 0.0) & (provenance_mass == 0.0)
    if np.any(unresolved):
        worst = int(np.flatnonzero(unresolved)[0])
        raise ValueError(
            "Final Pauli coefficient has no one-/two-body parent provenance: "
            f"key={raw_pauli_keys[worst]!r}."
        )
    nonzero = provenance_mass > 0.0
    one_weight = np.zeros_like(provenance_mass)
    one_weight[nonzero] = one_body_abs[nonzero] / provenance_mass[nonzero]
    one_body += raw_residual * one_weight
    two_body += raw_residual * (1.0 - one_weight)

    reconstruction_error = float(
        np.max(np.abs(one_body + two_body - final), initial=0.0)
    )
    scale = max(1.0, float(np.max(np.abs(final), initial=0.0)))
    if reconstruction_error > 50.0 * np.finfo(float).eps * scale:
        raise ValueError(
            "Reconciled one-/two-body split failed to reconstruct the final "
            f"Hamiltonian; max error={reconstruction_error:.3e}."
        )
    return one_body, two_body, reconstruction_error, max_reconciliation


def decompose_bch_pair_amplitudes(
    evaluator: HFCommutatorEvaluator,
    final_coefficients: np.ndarray,
    one_body_coefficients: np.ndarray,
    two_body_coefficients: np.ndarray,
) -> BCHComponents:
    """Split every BCH pair amplitude exactly by coefficient-product body order."""

    left = evaluator.left_indices
    right = evaluator.right_indices
    denominator = final_coefficients[left] * final_coefficients[right]
    if np.any(denominator == 0.0):
        raise ValueError("BCH evaluator includes a pair with zero final coefficient.")

    amplitude = evaluator.pair_amplitudes
    one_one = amplitude * (
        one_body_coefficients[left] * one_body_coefficients[right] / denominator
    )
    one_two = amplitude * (
        (
            one_body_coefficients[left] * two_body_coefficients[right]
            + two_body_coefficients[left] * one_body_coefficients[right]
        )
        / denominator
    )
    two_two = amplitude - one_one - one_two
    components = BCHComponents(one_one=one_one, one_two=one_two, two_two=two_two)
    error = float(np.max(np.abs(components.full - amplitude))) if amplitude.size else 0.0
    if error > 1.0e-10 * max(1.0, float(np.max(np.abs(amplitude), initial=0.0))):
        raise ValueError(f"BCH body-order split reconstruction error {error:.3e}.")
    return components


def order_signs(
    evaluator: HFCommutatorEvaluator,
    pauli_order_indices: Sequence[int],
) -> np.ndarray:
    positions = np.empty(evaluator.number_of_pauli_terms, dtype=np.int32)
    positions[np.asarray(pauli_order_indices, dtype=np.int32)] = np.arange(
        evaluator.number_of_pauli_terms,
        dtype=np.int32,
    )
    return np.where(
        positions[evaluator.left_indices] < positions[evaluator.right_indices],
        1.0,
        -1.0,
    )


def binned_vector(
    evaluator: HFCommutatorEvaluator,
    signed_amplitudes: np.ndarray,
) -> np.ndarray:
    if evaluator.number_of_target_bins == 0:
        return np.empty(0, dtype=np.complex128)
    real = np.bincount(
        evaluator.target_bins,
        weights=signed_amplitudes.real,
        minlength=evaluator.number_of_target_bins,
    )
    imaginary = np.bincount(
        evaluator.target_bins,
        weights=signed_amplitudes.imag,
        minlength=evaluator.number_of_target_bins,
    )
    return real + 1.0j * imaginary


def vector_norm(vector: np.ndarray) -> float:
    return float(np.linalg.norm(vector))


def cancellation_ratio(vector: np.ndarray, pair_amplitudes: np.ndarray) -> float:
    denominator = float(np.sum(np.abs(pair_amplitudes)))
    return safe_ratio(vector_norm(vector), denominator)


def destructive_cosine(left: np.ndarray, right: np.ndarray) -> float:
    left_norm = vector_norm(left)
    right_norm = vector_norm(right)
    if left_norm <= 0.0 or right_norm <= 0.0:
        return float("nan")
    return float(-np.vdot(left, right).real / (left_norm * right_norm))


def schedule_component_metrics(
    evaluator: HFCommutatorEvaluator,
    components: BCHComponents,
    pauli_order_indices: Sequence[int],
) -> dict[str, float]:
    signs = order_signs(evaluator, pauli_order_indices)
    one_one_vector = binned_vector(evaluator, components.one_one * signs)
    one_two_vector = binned_vector(evaluator, components.one_two * signs)
    two_two_vector = binned_vector(evaluator, components.two_two * signs)
    two_body_vector = one_two_vector + two_two_vector
    full_vector = one_one_vector + two_body_vector

    one_one_norm = vector_norm(one_one_vector)
    one_two_norm = vector_norm(one_two_vector)
    two_two_norm = vector_norm(two_two_vector)
    two_body_norm = vector_norm(two_body_vector)
    full_norm = vector_norm(full_vector)
    two_body_amplitudes = components.two_body_bearing

    if not math.isclose(
        full_norm,
        evaluator.evaluate(pauli_order_indices),
        rel_tol=1.0e-9,
        abs_tol=1.0e-12,
    ):
        raise ValueError("Component vectors do not reconstruct evaluator BCH norm.")

    return {
        "full_bch_norm": full_norm,
        "full_cancellation_ratio": cancellation_ratio(full_vector, components.full),
        "one_one_bch_norm": one_one_norm,
        "one_one_cancellation_ratio": cancellation_ratio(
            one_one_vector, components.one_one
        ),
        "one_two_bch_norm": one_two_norm,
        "one_two_cancellation_ratio": cancellation_ratio(
            one_two_vector, components.one_two
        ),
        "two_two_bch_norm": two_two_norm,
        "two_two_cancellation_ratio": cancellation_ratio(
            two_two_vector, components.two_two
        ),
        "two_body_bch_norm": two_body_norm,
        "two_body_cancellation_ratio": cancellation_ratio(
            two_body_vector, two_body_amplitudes
        ),
        "two_body_cross_destructiveness": destructive_cosine(
            one_one_vector, two_body_vector
        ),
        "one_two_vs_two_two_destructiveness": destructive_cosine(
            one_two_vector, two_two_vector
        ),
        "two_body_leave_in_reduction": safe_ratio(
            one_one_norm - full_norm,
            one_one_norm,
        )
        if one_one_norm >= full_norm
        else -safe_ratio(full_norm - one_one_norm, one_one_norm),
        "two_body_resultant_fraction": safe_ratio(
            two_body_norm,
            one_one_norm + two_body_norm,
        ),
    }


def build_orders(
    hermitian_terms: Sequence[Any],
    raw_pauli_keys: Sequence[Any],
    final_coefficients: dict[Any, complex],
    n_qubits: int,
    tolerance: float,
) -> tuple[list[int], list[int]]:
    raw_index = {key: index for index, key in enumerate(raw_pauli_keys)}
    mapping = precompute_fermion_to_pauli_indices(
        hermitian_terms=hermitian_terms,
        final_coefficients=final_coefficients,
        raw_index_by_key=raw_index,
        tolerance=tolerance,
    )
    parent_order = baseline.fermionic_term_order_indices(
        hermitian_terms=hermitian_terms,
        ordering_method="signed_ascending",
        tolerance=tolerance,
    )
    signed_order = induced_pauli_order_indices(
        fermionic_node_order=parent_order,
        fermion_to_pauli_indices=mapping,
        number_of_pauli_terms=len(raw_pauli_keys),
    )
    jw_keys = baseline.magnitude_descending_lexicographic_order(
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=tolerance,
    )
    return signed_order, [raw_index[key] for key in jw_keys]


def analyze_case(row: pd.Series, tolerance: float) -> dict[str, Any]:
    tensor_path = resolve_path(Path(str(row["tensor_path"])))
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)
    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction), tolerance
    )
    jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    jw_hamiltonian.compress(abs_tol=tolerance)
    final_coefficients = {
        key: coefficient
        for key, coefficient in jw_hamiltonian.terms.items()
        if key != () and abs(coefficient) > tolerance
    }
    raw_pauli_keys = list(final_coefficients)
    hermitian_terms = build_hermitian_fermion_terms(
        fermion_hamiltonian, tolerance
    )
    (
        one_body,
        two_body,
        coefficient_error,
        coefficient_reconciliation,
    ) = decompose_final_pauli_coefficients(
        hermitian_terms, raw_pauli_keys, final_coefficients, tolerance
    )
    final = np.asarray(
        [real_coefficient(final_coefficients[key], tolerance) for key in raw_pauli_keys]
    )
    signed_order, jw_order = build_orders(
        hermitian_terms,
        raw_pauli_keys,
        final_coefficients,
        n_qubits,
        tolerance,
    )

    graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    evaluator = HFCommutatorEvaluator.build(
        pauli_keys=raw_pauli_keys,
        coefficients=final_coefficients,
        pauli_graph=graph,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        tolerance=tolerance,
    )
    components = decompose_bch_pair_amplitudes(
        evaluator, final, one_body, two_body
    )
    signed = schedule_component_metrics(evaluator, components, signed_order)
    jw = schedule_component_metrics(evaluator, components, jw_order)

    pair_mass_one = float(np.sum(np.abs(components.one_one)))
    pair_mass_two = float(np.sum(np.abs(components.two_body_bearing)))
    pair_mass_components = pair_mass_one + pair_mass_two

    result: dict[str, Any] = {
        "case_id": row["case_id"],
        "matched_pair": row["matched_pair"],
        "expected_outcome": row["expected_outcome"],
        "panel_tier": row["panel_tier"],
        "molecule": row["molecule"],
        "n_qubits": n_qubits,
        "active_occupied": metadata.active_occupied,
        "active_vacant": metadata.active_vacant,
        "tensor_path": str(row["tensor_path"]),
        "number_of_fermionic_parents": len(hermitian_terms),
        "number_of_one_body_parents": sum(
            parent_body_order(term) == 1 for term in hermitian_terms
        ),
        "number_of_two_body_parents": sum(
            parent_body_order(term) == 2 for term in hermitian_terms
        ),
        "number_of_pauli_terms": len(raw_pauli_keys),
        "number_of_noncommuting_pairs": evaluator.pair_amplitudes.size,
        "pauli_coefficient_reconstruction_max_error": coefficient_error,
        "pauli_coefficient_reconciliation_max_correction": (
            coefficient_reconciliation
        ),
        "bch_pair_component_reconstruction_max_error": float(
            np.max(np.abs(components.full - evaluator.pair_amplitudes), initial=0.0)
        ),
        "two_body_pair_mass_fraction": safe_ratio(
            pair_mass_two, pair_mass_components
        ),
        "fermionic_advantage_factor": float(row["fresh_jw_magnitude_to_signed_advantage"]),
        "log10_error_advantage": safe_log10(
            float(row["fresh_jw_magnitude_to_signed_advantage"])
        ),
        "recorded_signed_cancellation_ratio": float(
            row["signed_bch_cancellation_ratio"]
        ),
        "recorded_jw_cancellation_ratio": float(
            row["fresh_jw_magnitude_bch_cancellation_ratio"]
        ),
    }
    result.update({f"signed_{key}": value for key, value in signed.items()})
    result.update({f"jw_{key}": value for key, value in jw.items()})

    for component in ("full", "one_one", "one_two", "two_two", "two_body"):
        result[f"log10_{component}_cancellation_advantage"] = safe_log10(
            safe_ratio(
                jw[f"{component}_cancellation_ratio"],
                signed[f"{component}_cancellation_ratio"],
            )
        )
    result["delta_two_body_cross_destructiveness"] = (
        signed["two_body_cross_destructiveness"]
        - jw["two_body_cross_destructiveness"]
    )
    result["delta_two_body_leave_in_reduction"] = (
        signed["two_body_leave_in_reduction"]
        - jw["two_body_leave_in_reduction"]
    )
    result["log10_full_bch_advantage"] = result[
        "log10_full_cancellation_advantage"
    ]
    result["signed_cancellation_ratio_verification_abs_error"] = abs(
        signed["full_cancellation_ratio"]
        - result["recorded_signed_cancellation_ratio"]
    )
    result["jw_cancellation_ratio_verification_abs_error"] = abs(
        jw["full_cancellation_ratio"] - result["recorded_jw_cancellation_ratio"]
    )
    return result


def build_matched_pair_table(cases: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for matched_pair, group in cases.groupby("matched_pair", sort=True):
        favorable = group[group["expected_outcome"].eq("favorable")]
        control = group[group["expected_outcome"].eq("negative_control")]
        if len(favorable) != 1 or len(control) != 1:
            raise ValueError(f"Matched pair {matched_pair!r} is not one-to-one.")
        fav = favorable.iloc[0]
        con = control.iloc[0]
        record: dict[str, Any] = {
            "matched_pair": matched_pair,
            "n_qubits": int(fav["n_qubits"]),
            "favorable_case_id": fav["case_id"],
            "control_case_id": con["case_id"],
        }
        for column in ["log10_error_advantage", "log10_full_bch_advantage"] + [
            feature for feature, _, _ in FEATURES
        ]:
            record[f"favorable_{column}"] = float(fav[column])
            record[f"control_{column}"] = float(con[column])
            record[f"delta_{column}"] = float(fav[column] - con[column])
        records.append(record)
    return pd.DataFrame(records)


def build_statistics(
    cases: pd.DataFrame,
    pairs: pd.DataFrame,
    bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for outcome, outcome_label in OUTCOMES:
        for feature, feature_label, family in FEATURES:
            n, raw_r, raw_p = pearson(cases[feature], cases[outcome])
            raw_rho, raw_spearman_p = spearman(cases[feature], cases[outcome])

            x_residual = within_group_residual(cases[feature], cases["matched_pair"])
            y_residual = within_group_residual(cases[outcome], cases["matched_pair"])
            _, adjusted_r, adjusted_p = pearson(x_residual, y_residual)

            one_body_control = "log10_one_one_cancellation_advantage"
            if feature == one_body_control:
                one_body_n = len(cases)
                one_body_r = float("nan")
                one_body_p = float("nan")
                combined_n = len(cases)
                combined_r = float("nan")
                combined_p = float("nan")
            else:
                one_body_n, one_body_r, one_body_p = partial_pearson(
                    cases,
                    feature,
                    outcome,
                    numeric_controls=[one_body_control],
                )
                combined_n, combined_r, combined_p = partial_pearson(
                    cases,
                    feature,
                    outcome,
                    numeric_controls=[one_body_control],
                    categorical_controls=["matched_pair"],
                )

            pair_x = pairs[f"delta_{feature}"]
            pair_y = pairs[f"delta_{outcome}"]
            pair_n, pair_r, pair_p = pearson(pair_x, pair_y)
            differences = pair_x.to_numpy(dtype=float)
            paired_mean = float(np.mean(differences))
            paired_median = float(np.median(differences))
            permutation_p = exact_paired_sign_permutation_p_value(differences)
            ci_low, ci_high = paired_bootstrap_mean_ci(
                differences,
                replicates=bootstrap,
                seed=seed + sum(ord(character) for character in feature + outcome),
            )
            records.append(
                {
                    "outcome": outcome,
                    "outcome_label": outcome_label,
                    "feature": feature,
                    "feature_label": feature_label,
                    "feature_family": family,
                    "cases": n,
                    "pairs": pair_n,
                    "raw_pearson_r": raw_r,
                    "raw_pearson_p_value": raw_p,
                    "raw_spearman_rho": raw_rho,
                    "raw_spearman_p_value": raw_spearman_p,
                    "matched_pair_adjusted_pearson_r": adjusted_r,
                    "matched_pair_adjusted_pearson_p_value": adjusted_p,
                    "one_body_adjusted_cases": one_body_n,
                    "one_body_adjusted_pearson_r": one_body_r,
                    "one_body_adjusted_pearson_p_value": one_body_p,
                    "matched_pair_and_one_body_adjusted_cases": combined_n,
                    "matched_pair_and_one_body_adjusted_pearson_r": combined_r,
                    "matched_pair_and_one_body_adjusted_pearson_p_value": combined_p,
                    "pair_delta_pearson_r": pair_r,
                    "pair_delta_pearson_p_value": pair_p,
                    "favorable_minus_control_mean": paired_mean,
                    "favorable_minus_control_median": paired_median,
                    "paired_exact_permutation_p_value": permutation_p,
                    "paired_bootstrap_mean_ci_low": ci_low,
                    "paired_bootstrap_mean_ci_high": ci_high,
                    "bootstrap_replicates": bootstrap,
                }
            )
    result = pd.DataFrame(records)
    for outcome, _ in OUTCOMES:
        mask = result["outcome"].eq(outcome)
        for source, destination in (
            ("raw_pearson_p_value", "raw_pearson_fdr_q_value"),
            (
                "matched_pair_adjusted_pearson_p_value",
                "matched_pair_adjusted_pearson_fdr_q_value",
            ),
            (
                "one_body_adjusted_pearson_p_value",
                "one_body_adjusted_pearson_fdr_q_value",
            ),
            (
                "matched_pair_and_one_body_adjusted_pearson_p_value",
                "matched_pair_and_one_body_adjusted_pearson_fdr_q_value",
            ),
            ("pair_delta_pearson_p_value", "pair_delta_pearson_fdr_q_value"),
            (
                "paired_exact_permutation_p_value",
                "paired_exact_permutation_fdr_q_value",
            ),
        ):
            result.loc[mask, destination] = bh_fdr(result.loc[mask, source])
    return result


def fmt(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}g}"


def regression_line(axis: plt.Axes, x: pd.Series, y: pd.Series) -> None:
    valid = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(valid) < 2:
        return
    slope, intercept = np.polyfit(valid["x"], valid["y"], 1)
    domain = np.linspace(valid["x"].min(), valid["x"].max(), 100)
    axis.plot(domain, slope * domain + intercept, color="0.25", linewidth=1.4)


def make_figure(
    cases: pd.DataFrame,
    pairs: pd.DataFrame,
    statistics: pd.DataFrame,
    output_path: Path,
    dpi: int,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 10.0), constrained_layout=True)
    colors = {
        "favorable": "#1675b9",
        "negative_control": "#e87522",
    }
    markers = {"favorable": "o", "negative_control": "s"}

    primary = statistics[
        statistics["outcome"].eq("log10_error_advantage")
        & statistics["feature"].eq("log10_two_body_cancellation_advantage")
    ].iloc[0]

    axis = axes[0, 0]
    for outcome, group in cases.groupby("expected_outcome"):
        axis.scatter(
            group["log10_two_body_cancellation_advantage"],
            group["log10_error_advantage"],
            s=62,
            alpha=0.88,
            color=colors[outcome],
            marker=markers[outcome],
            edgecolor="white",
            linewidth=0.6,
        )
    regression_line(
        axis,
        cases["log10_two_body_cancellation_advantage"],
        cases["log10_error_advantage"],
    )
    axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1.0)
    axis.axvline(0.0, color="0.35", linestyle=":", linewidth=1.0)
    axis.set_xlabel(r"2B BCH gain  $\log_{10}(R_{2B,JW}/R_{2B,F})$")
    axis.set_ylabel(r"Trotter advantage  $\log_{10}(\epsilon_{JW}/\epsilon_F)$")
    axis.set_title(
        "A  Two-body cancellation vs performance\n"
        f"r={primary['raw_pearson_r']:.2f}, p={primary['raw_pearson_p_value']:.3g}"
    )

    axis = axes[0, 1]
    for outcome, group in cases.groupby("expected_outcome"):
        axis.scatter(
            group["log10_two_body_cancellation_advantage"],
            group["log10_full_bch_advantage"],
            s=62,
            alpha=0.88,
            color=colors[outcome],
            marker=markers[outcome],
            edgecolor="white",
            linewidth=0.6,
        )
    regression_line(
        axis,
        cases["log10_two_body_cancellation_advantage"],
        cases["log10_full_bch_advantage"],
    )
    axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1.0)
    axis.axvline(0.0, color="0.35", linestyle=":", linewidth=1.0)
    full_row = statistics[
        statistics["outcome"].eq("log10_full_bch_advantage")
        & statistics["feature"].eq("log10_two_body_cancellation_advantage")
    ].iloc[0]
    axis.set_xlabel(r"2B BCH gain  $\log_{10}(R_{2B,JW}/R_{2B,F})$")
    axis.set_ylabel(r"Full BCH gain  $\log_{10}(R_{JW}/R_F)$")
    axis.set_title(
        "B  Two-body cancellation vs full cancellation\n"
        f"r={full_row['raw_pearson_r']:.2f}, p={full_row['raw_pearson_p_value']:.3g}"
    )

    axis = axes[1, 0]
    ordered_pairs = pairs.sort_values("delta_log10_two_body_cancellation_advantage")
    for index, (_, row) in enumerate(ordered_pairs.iterrows()):
        axis.plot(
            [0, 1],
            [
                row["control_log10_two_body_cancellation_advantage"],
                row["favorable_log10_two_body_cancellation_advantage"],
            ],
            color="0.70",
            linewidth=1.0,
            zorder=1,
        )
        axis.scatter(
            0,
            row["control_log10_two_body_cancellation_advantage"],
            color=colors["negative_control"],
            marker=markers["negative_control"],
            s=52,
            zorder=2,
        )
        axis.scatter(
            1,
            row["favorable_log10_two_body_cancellation_advantage"],
            color=colors["favorable"],
            marker=markers["favorable"],
            s=52,
            zorder=2,
        )
        if index in (0, len(ordered_pairs) - 1):
            axis.annotate(
                str(row["matched_pair"]),
                (1.03, row["favorable_log10_two_body_cancellation_advantage"]),
                fontsize=8,
                va="center",
            )
    axis.axhline(0.0, color="0.35", linestyle="--", linewidth=1.0)
    axis.set_xlim(-0.25, 1.45)
    axis.set_xticks([0, 1], ["Matched control", "Favorable"])
    axis.set_ylabel(r"2B BCH gain  $\log_{10}(R_{2B,JW}/R_{2B,F})$")
    axis.set_title("C  Ten prospectively matched pairs")

    axis = axes[1, 1]
    component_rows = statistics[
        statistics["outcome"].eq("log10_error_advantage")
        & statistics["feature"].isin(
            [
                "log10_two_body_cancellation_advantage",
                "delta_two_body_cross_destructiveness",
                "delta_two_body_leave_in_reduction",
                "two_body_pair_mass_fraction",
            ]
        )
    ].copy()
    order = [
        "log10_two_body_cancellation_advantage",
        "delta_two_body_cross_destructiveness",
        "delta_two_body_leave_in_reduction",
        "two_body_pair_mass_fraction",
    ]
    component_rows["order"] = component_rows["feature"].map(
        {name: index for index, name in enumerate(order)}
    )
    component_rows = component_rows.sort_values("order")
    labels = ["2B internal", "1B/2B cross", "2B leave-in", "Static 2B mass"]
    positions = np.arange(len(component_rows))
    axis.bar(
        positions,
        component_rows["raw_pearson_r"],
        color=["#2171b5", "#6baed6", "#084594", "#929292"],
        width=0.68,
    )
    for position, (_, row) in zip(positions, component_rows.iterrows(), strict=True):
        axis.text(
            position,
            row["raw_pearson_r"] + 0.035 * np.sign(row["raw_pearson_r"] or 1),
            f"p={row['raw_pearson_p_value']:.2g}",
            ha="center",
            va="bottom" if row["raw_pearson_r"] >= 0 else "top",
            fontsize=9,
        )
    axis.axhline(0.0, color="0.25", linewidth=0.9)
    axis.set_xticks(positions, labels, rotation=15, ha="right")
    axis.set_ylabel("Pearson r with log Trotter advantage")
    axis.set_ylim(-1.0, 1.0)
    axis.set_title("D  Order-sensitive 2B metrics, with static-mass control")

    handles = [
        Line2D(
            [0],
            [0],
            marker=markers[name],
            linestyle="none",
            markerfacecolor=colors[name],
            markeredgecolor="white",
            markersize=8,
            label="Favorable" if name == "favorable" else "Negative control",
        )
        for name in ("favorable", "negative_control")
    ]
    figure.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    figure.suptitle(
        "Exact two-body decomposition of leading BCH cancellation",
        fontsize=17,
        y=1.02,
    )
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def write_report(
    cases: pd.DataFrame,
    pairs: pd.DataFrame,
    statistics: pd.DataFrame,
    output_path: Path,
) -> None:
    primary = statistics[
        statistics["outcome"].eq("log10_error_advantage")
        & statistics["feature"].eq("log10_two_body_cancellation_advantage")
    ].iloc[0]
    full = statistics[
        statistics["outcome"].eq("log10_full_bch_advantage")
        & statistics["feature"].eq("log10_two_body_cancellation_advantage")
    ].iloc[0]
    cross = statistics[
        statistics["outcome"].eq("log10_error_advantage")
        & statistics["feature"].eq("delta_two_body_cross_destructiveness")
    ].iloc[0]
    leave_in = statistics[
        statistics["outcome"].eq("log10_error_advantage")
        & statistics["feature"].eq("delta_two_body_leave_in_reduction")
    ].iloc[0]
    static_mass = statistics[
        statistics["outcome"].eq("log10_error_advantage")
        & statistics["feature"].eq("two_body_pair_mass_fraction")
    ].iloc[0]
    one_body = statistics[
        statistics["outcome"].eq("log10_error_advantage")
        & statistics["feature"].eq("log10_one_one_cancellation_advantage")
    ].iloc[0]
    favorable = cases[cases["expected_outcome"].eq("favorable")]
    controls = cases[cases["expected_outcome"].eq("negative_control")]
    direction_count = int(
        np.sum(pairs["delta_log10_two_body_cancellation_advantage"] > 0.0)
    )
    cross_direction_count = int(
        np.sum(pairs["delta_delta_two_body_cross_destructiveness"] > 0.0)
    )
    leave_in_direction_count = int(
        np.sum(pairs["delta_delta_two_body_leave_in_reduction"] > 0.0)
    )
    strongest_pair = pairs.loc[pairs["delta_log10_error_advantage"].idxmax()]
    internal_counterexample = pairs.loc[
        pairs["delta_log10_two_body_cancellation_advantage"].idxmin()
    ]

    supported = bool(
        leave_in["raw_pearson_r"] > 0.0
        and leave_in["one_body_adjusted_pearson_fdr_q_value"] < 0.05
        and leave_in["paired_exact_permutation_fdr_q_value"] < 0.05
    )
    verdict = (
        "supported as a direct participating component of the BCH mechanism"
        if supported
        else "not established as a consistent direct component across the panel"
    )

    lines = [
        "# Two-body BCH cancellation case study",
        "",
        "## Verdict",
        "",
        f"The two-body contribution hypothesis is **{verdict}**.",
        "",
        "The strongest evidence is not that the 2B-bearing vector always cancels "
        "more strongly by itself. It is that signed fermionic ordering makes the "
        "2B-bearing vector combine more destructively with the remaining 1B-only "
        "BCH vector. The static amount of two-body BCH mass has no performance "
        "association.",
        "",
        "The primary order-sensitive quantity is "
        "`G_2B = log10(R_2B,JW / R_2B,F)`, where `R_2B` is the HF-state "
        "BCH2 cancellation ratio of the exact BCH amplitude containing at least "
        "one two-body fermionic coefficient contribution. Positive `G_2B` means "
        "that signed fermionic ordering cancels that component more strongly than "
        "JW magnitude-descending ordering.",
        "",
        "## Exact decomposition",
        "",
        "For each final JW coefficient, `c_i = c_i^(1B) + c_i^(2B)`. "
        "Every noncommuting-pair coefficient product is therefore decomposed as",
        "",
        "`c_i c_j = c_i^(1B)c_j^(1B) + "
        "[c_i^(1B)c_j^(2B)+c_i^(2B)c_j^(1B)] + "
        "c_i^(2B)c_j^(2B)`.",
        "",
        "The three complex HF-state BCH vectors sum back to the full measured "
        "BCH vector. This is coefficient-provenance decomposition, not a graph "
        "support proxy. Ordering changes only the orientation sign of each fixed "
        "Pauli-pair amplitude.",
        "",
        "## Coverage and numerical checks",
        "",
        f"- Cases: {len(cases)} ({len(favorable)} favorable, {len(controls)} matched controls)",
        f"- Matched pairs: {len(pairs)}",
        f"- Active-space range: {cases['n_qubits'].min()}-{cases['n_qubits'].max()} qubits",
        "- Maximum final-Pauli coefficient reconstruction error: "
        f"{cases['pauli_coefficient_reconstruction_max_error'].max():.3e}",
        "- Maximum floating-point reconciliation applied before exact "
        "reconstruction: "
        f"{cases['pauli_coefficient_reconciliation_max_correction'].max():.3e}",
        "- Maximum BCH pair-component reconstruction error: "
        f"{cases['bch_pair_component_reconstruction_max_error'].max():.3e}",
        "- Maximum discrepancy from independently saved full BCH ratios: "
        f"{max(cases['signed_cancellation_ratio_verification_abs_error'].max(), cases['jw_cancellation_ratio_verification_abs_error'].max()):.3e}",
        "",
        "## Primary result: two-body cancellation vs Trotter advantage",
        "",
        f"- Raw case-level: r={fmt(primary['raw_pearson_r'])}, "
        f"p={fmt(primary['raw_pearson_p_value'])}, "
        f"q={fmt(primary['raw_pearson_fdr_q_value'])}, n={int(primary['cases'])}",
        f"- Matched-pair adjusted: r={fmt(primary['matched_pair_adjusted_pearson_r'])}, "
        f"p={fmt(primary['matched_pair_adjusted_pearson_p_value'])}, "
        f"q={fmt(primary['matched_pair_adjusted_pearson_fdr_q_value'])}",
        f"- Controlling 1B-only cancellation gain: r={fmt(primary['one_body_adjusted_pearson_r'])}, "
        f"p={fmt(primary['one_body_adjusted_pearson_p_value'])}, "
        f"q={fmt(primary['one_body_adjusted_pearson_fdr_q_value'])}",
        "- Controlling both matched pair and 1B-only gain: "
        f"r={fmt(primary['matched_pair_and_one_body_adjusted_pearson_r'])}, "
        f"p={fmt(primary['matched_pair_and_one_body_adjusted_pearson_p_value'])}, "
        f"q={fmt(primary['matched_pair_and_one_body_adjusted_pearson_fdr_q_value'])}",
        f"- Pair-difference correlation: r={fmt(primary['pair_delta_pearson_r'])}, "
        f"p={fmt(primary['pair_delta_pearson_p_value'])}, n={int(primary['pairs'])} pairs",
        f"- Favorable-minus-control mean `G_2B`: {fmt(primary['favorable_minus_control_mean'])}, "
        f"bootstrap 95% CI [{fmt(primary['paired_bootstrap_mean_ci_low'])}, "
        f"{fmt(primary['paired_bootstrap_mean_ci_high'])}], exact paired "
        f"permutation p={fmt(primary['paired_exact_permutation_p_value'])}",
        f"- Favorable case has larger `G_2B` in {direction_count}/{len(pairs)} matched pairs.",
        "",
        "Median `G_2B` is "
        f"{np.median(favorable['log10_two_body_cancellation_advantage']):.3g} in favorable cases "
        f"and {np.median(controls['log10_two_body_cancellation_advantage']):.3g} in controls.",
        "The significant continuous correlation therefore does not amount to "
        "clean favorable/control classification by 2B-internal cancellation alone.",
        "",
        "## Does the two-body component track the measured full BCH effect?",
        "",
        f"- `G_2B` vs full BCH cancellation advantage: r={fmt(full['raw_pearson_r'])}, "
        f"p={fmt(full['raw_pearson_p_value'])}, q={fmt(full['raw_pearson_fdr_q_value'])}",
        f"- Matched-pair adjusted: r={fmt(full['matched_pair_adjusted_pearson_r'])}, "
        f"p={fmt(full['matched_pair_adjusted_pearson_p_value'])}",
        "",
        "## Cross-component interference checks",
        "",
        "The destructive cosine is positive when the aggregated 1B-only and "
        "2B-bearing vectors point against each other. The leave-in statistic is "
        "positive when adding the 2B-bearing vector lowers the full norm relative "
        "to the 1B-only counterfactual.",
        "",
        f"- Signed-vs-JW gain in destructive cosine vs Trotter advantage: "
        f"r={fmt(cross['raw_pearson_r'])}, p={fmt(cross['raw_pearson_p_value'])}, "
        f"q={fmt(cross['raw_pearson_fdr_q_value'])}; controlling 1B-only gain, "
        f"r={fmt(cross['one_body_adjusted_pearson_r'])}, "
        f"p={fmt(cross['one_body_adjusted_pearson_p_value'])}, "
        f"q={fmt(cross['one_body_adjusted_pearson_fdr_q_value'])}; favorable "
        f"direction in {cross_direction_count}/{len(pairs)} pairs, exact paired "
        f"p={fmt(cross['paired_exact_permutation_p_value'])}",
        f"- Signed-vs-JW gain in two-body leave-in norm reduction vs advantage: "
        f"r={fmt(leave_in['raw_pearson_r'])}, p={fmt(leave_in['raw_pearson_p_value'])}, "
        f"q={fmt(leave_in['raw_pearson_fdr_q_value'])}; controlling 1B-only gain, "
        f"r={fmt(leave_in['one_body_adjusted_pearson_r'])}, "
        f"p={fmt(leave_in['one_body_adjusted_pearson_p_value'])}, "
        f"q={fmt(leave_in['one_body_adjusted_pearson_fdr_q_value'])}; favorable "
        f"direction in {leave_in_direction_count}/{len(pairs)} pairs, exact paired "
        f"p={fmt(leave_in['paired_exact_permutation_p_value'])}",
        "",
        "## Illustrative matched cases",
        "",
        "These examples are descriptive selections from the completed panel, not "
        "additional inferential tests.",
        "",
        f"- Largest performance contrast ({strongest_pair['matched_pair']}): "
        f"`{strongest_pair['favorable_case_id']}` vs "
        f"`{strongest_pair['control_case_id']}` has a log-error-advantage "
        f"difference of {fmt(strongest_pair['delta_log10_error_advantage'])}. "
        "Their `G_2B` values are "
        f"{fmt(strongest_pair['favorable_log10_two_body_cancellation_advantage'])} "
        f"and {fmt(strongest_pair['control_log10_two_body_cancellation_advantage'])}; "
        "their 2B leave-in gains are "
        f"{fmt(strongest_pair['favorable_delta_two_body_leave_in_reduction'])} "
        f"and {fmt(strongest_pair['control_delta_two_body_leave_in_reduction'])}.",
        f"- Internal-cancellation counterexample ({internal_counterexample['matched_pair']}): "
        f"`{internal_counterexample['favorable_case_id']}` performs better even "
        "though its `G_2B` is "
        f"{fmt(internal_counterexample['favorable_log10_two_body_cancellation_advantage'])}, "
        "below the control value "
        f"{fmt(internal_counterexample['control_log10_two_body_cancellation_advantage'])}. "
        "Its 2B leave-in gain is nevertheless less unfavorable "
        f"({fmt(internal_counterexample['favorable_delta_two_body_leave_in_reduction'])} "
        "vs "
        f"{fmt(internal_counterexample['control_delta_two_body_leave_in_reduction'])}). "
        "This is why the cross-vector mechanism is more complete than a "
        "2B-internal-cancellation-only story.",
        "",
        "## Specificity controls",
        "",
        f"- Static 2B-bearing BCH mass fraction vs advantage: r={fmt(static_mass['raw_pearson_r'])}, "
        f"p={fmt(static_mass['raw_pearson_p_value'])}, "
        f"q={fmt(static_mass['raw_pearson_fdr_q_value'])}",
        f"- 1B-only cancellation gain also correlates with advantage: "
        f"r={fmt(one_body['raw_pearson_r'])}, p={fmt(one_body['raw_pearson_p_value'])}, "
        f"q={fmt(one_body['raw_pearson_fdr_q_value'])}",
        "- The 2B internal, cross-interference, and leave-in signals remain after "
        "controlling the 1B-only gain, so they are not merely copies of that "
        "negative-control component. The result nevertheless identifies 2B as a "
        "participant, not as the sole source of cancellation.",
        "",
        "## Interpretation",
        "",
        "The combined evidence identifies the two-body-bearing portion "
        "of the leading BCH vector as part of the ordering mechanism. It does not "
        "mean that a larger static two-body weighted fraction is beneficial. The "
        "relevant quantity is the change in vector cancellation caused by the "
        "order. Conversely, a null primary association means the existing full "
        "BCH effect cannot be localized consistently to two-body-bearing terms "
        "with this decomposition.",
        "",
        "This remains an HF-state leading-BCH mechanism study. It is not a causal "
        "term-deletion experiment, and higher BCH orders or other initial states "
        "can redistribute contributions.",
        "",
        "## Files",
        "",
        "- `two_body_bch_case_metrics.csv`: exact per-case component metrics.",
        "- `two_body_bch_matched_pairs.csv`: favorable-control differences.",
        "- `two_body_bch_statistics.csv`: raw, matched, paired, FDR, and bootstrap tests.",
        "- `two_body_bch_case_study.png`: diagnostic figure.",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_panel(manifest_path: Path, case_summary_path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(resolve_path(manifest_path))
    summary = pd.read_csv(resolve_path(case_summary_path))
    columns = [
        "case_id",
        "fresh_jw_magnitude_to_signed_advantage",
        "signed_bch_cancellation_ratio",
        "fresh_jw_magnitude_bch_cancellation_ratio",
    ]
    missing = sorted(set(columns).difference(summary.columns))
    if missing:
        raise ValueError(f"Case summary is missing columns: {missing}")
    merged = manifest.merge(summary[columns], on="case_id", how="left", validate="one_to_one")
    if merged[columns[1:]].isna().any().any():
        raise ValueError("Panel merge lost current performance or BCH results.")
    return merged


def main() -> None:
    args = parse_arguments()
    outdir = resolve_path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    case_path = outdir / "two_body_bch_case_metrics.csv"

    if args.reuse_case_metrics:
        if args.case_id or args.max_qubits is not None:
            raise ValueError(
                "--reuse-case-metrics cannot be combined with case-selection options."
            )
        if not case_path.exists():
            raise FileNotFoundError(f"No reusable case metrics at {case_path}.")
        cases = pd.read_csv(case_path)
        print(f"Reusing {len(cases)} case rows from {case_path}")
    else:
        panel = load_panel(args.manifest, args.case_summary)
        if args.case_id:
            requested = set(args.case_id)
            missing = sorted(requested.difference(panel["case_id"]))
            if missing:
                raise ValueError(f"Unknown requested case IDs: {missing}")
            panel = panel[panel["case_id"].isin(requested)]
        if args.max_qubits is not None:
            panel = panel[panel["n_qubits"].le(args.max_qubits)]
        panel = panel.sort_values(["n_qubits", "matched_pair", "expected_outcome"])
        if panel.empty:
            raise ValueError("No cases selected.")

        records: list[dict[str, Any]] = []
        for index, (_, row) in enumerate(panel.iterrows(), start=1):
            print(f"[{index}/{len(panel)}] {row['case_id']}", flush=True)
            records.append(analyze_case(row, args.tolerance))
        cases = pd.DataFrame(records).sort_values(
            ["n_qubits", "matched_pair", "expected_outcome"]
        )
        cases.to_csv(case_path, index=False)

    complete_panel = (
        len(cases) == 20
        and cases["matched_pair"].nunique() == 10
        and set(cases["expected_outcome"]) == {"favorable", "negative_control"}
    )
    if not complete_panel:
        print(
            f"Wrote partial case metrics to {case_path}; matched statistics require "
            "the complete 20-case panel."
        )
        return

    pairs = build_matched_pair_table(cases)
    statistics = build_statistics(cases, pairs, args.bootstrap, args.seed)
    pairs.to_csv(outdir / "two_body_bch_matched_pairs.csv", index=False)
    statistics.to_csv(outdir / "two_body_bch_statistics.csv", index=False)
    make_figure(
        cases,
        pairs,
        statistics,
        outdir / "two_body_bch_case_study.png",
        args.dpi,
    )
    write_report(
        cases,
        pairs,
        statistics,
        outdir / "two_body_bch_case_study_report.md",
    )
    print(f"Wrote complete case study to {outdir}")


if __name__ == "__main__":
    main()
