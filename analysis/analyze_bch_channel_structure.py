#!/usr/bin/env python3
"""Case-study test of BCH-channel structure behind fermionic ordering gains.

This analysis reuses QHAT's existing HFCommutatorEvaluator.  For every
noncommuting final-JW Pauli pair, that evaluator already stores

    * the pair's complex [H_i, H_j]|HF> amplitude, and
    * a target-bin label for the computational-basis error channel reached.

The script compares two schedules of the same final JW Pauli Hamiltonian:

    fermionic_signed_coefficient_lexicographic
    jw_magnitude_descending_lexicographic

It separates three ideas:

1. Collision opportunity (Hamiltonian-level, order invariant):
   How much BCH pair weight lands in channels that receive >= 2 contributions?

2. Realized within-channel cancellation (order dependent):
   For channel q, compare sum_e |a_e| with |sum_e s_e a_e|, where s_e is the
   pair orientation induced by the Trotter schedule.

3. Leading BCH norm and the already-measured Trotter advantage.

Primary hypothesis:

    signed_minus_jwmag_l1_channel_cancellation > 0

should accompany

    actual_log10_jwmag_over_fermionic_advantage > 0.

Run from the QHAT repository root.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
from openfermion import get_fermion_operator, jordan_wigner
from scipy import stats

try:
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        load_interaction_operator,
        parse_case_metadata,
    )
    from qhat.analysis.benchmark_b2_coloring_robustness import HFCommutatorEvaluator
except ImportError:
    import benchmark_b2_signed_coefficient_baseline as baseline
    from benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        load_interaction_operator,
        parse_case_metadata,
    )
    from benchmark_b2_coloring_robustness import HFCommutatorEvaluator


DEFAULT_INPUT = Path(
    "analysis/fermionic_aware_performance/fermionic_aware_case_performance.csv"
)
DEFAULT_OUTDIR = Path("analysis/bch_channel_structure")

# Small matched-size panel for the first mechanism check.  All are HGBS-5,
# 16-qubit cases.  It intentionally contains strong wins, near-neutral cases,
# and losses against JW descending-magnitude ordering.
DEFAULT_CASES = (
    "F-F_1.28_hgbs-5_as-014-002",
    "H2O_s-1.00_hgbs-5_as-006-010",
    "NH3_s-1.00_hgbs-5_as-006-010",
    "O-O_1.26_hgbs-5_as-006-010",
    "C-C_1.50_hgbs-5_as-006-010",
    "B-B_1.70_hgbs-5_as-006-010",
    "CH4_s-1.00_hgbs-5_as-006-010",
)

FERMIONIC_ORDER = "fermionic_signed_coefficient_lexicographic"
JW_MAGNITUDE_ORDER = "jw_magnitude_descending_lexicographic"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--case",
        action="append",
        dest="cases",
        help=(
            "Exact case_id from fermionic_aware_case_performance.csv. "
            "Repeat for multiple cases. If omitted, use the built-in matched "
            "16-qubit case-study panel."
        ),
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Coefficient/commutator tolerance. Default: repository tolerance.",
    )
    parser.add_argument(
        "--top-channels",
        type=int,
        default=20,
        help=(
            "Also write a compact CSV containing this many highest-weight "
            "BCH channels per case. Default: 20."
        ),
    )
    return parser.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator <= 0.0 or not math.isfinite(denominator):
        return math.nan
    return float(numerator / denominator)


def pair_signs(
    evaluator: HFCommutatorEvaluator,
    pauli_order_indices: Sequence[int],
) -> np.ndarray:
    """Return +/-1 orientation for each precomputed unordered Pauli pair."""
    order = np.asarray(pauli_order_indices, dtype=np.int32)
    if order.size != evaluator.number_of_pauli_terms:
        raise ValueError(
            "Pauli-order length does not match evaluator: "
            f"{order.size} != {evaluator.number_of_pauli_terms}."
        )
    positions = np.empty(evaluator.number_of_pauli_terms, dtype=np.int32)
    positions[order] = np.arange(order.size, dtype=np.int32)
    return np.where(
        positions[evaluator.left_indices] < positions[evaluator.right_indices],
        1.0,
        -1.0,
    )


def channel_arrays(
    evaluator: HFCommutatorEvaluator,
    pauli_order_indices: Sequence[int],
) -> dict[str, np.ndarray | float]:
    """Aggregate oriented BCH pair amplitudes inside HF output channels."""
    n_bins = evaluator.number_of_target_bins
    if evaluator.pair_amplitudes.size == 0:
        zeros = np.zeros(n_bins, dtype=float)
        return {
            "signs": np.empty(0, dtype=float),
            "counts": zeros.astype(np.int64),
            "raw_abs_by_channel": zeros.copy(),
            "net_complex_by_channel": np.zeros(n_bins, dtype=np.complex128),
            "net_abs_by_channel": zeros.copy(),
            "cancellation_by_channel": zeros.copy(),
            "pair_abs_sum": 0.0,
            "bch_norm": 0.0,
            "l1_net_sum": 0.0,
            "l1_channel_cancellation": 0.0,
        }

    signs = pair_signs(evaluator, pauli_order_indices)
    oriented = evaluator.pair_amplitudes * signs
    absolute_weights = np.abs(evaluator.pair_amplitudes)
    bins = evaluator.target_bins

    counts = np.bincount(bins, minlength=n_bins).astype(np.int64)
    raw_abs_by_channel = np.bincount(
        bins,
        weights=absolute_weights,
        minlength=n_bins,
    )
    net_real = np.bincount(
        bins,
        weights=oriented.real,
        minlength=n_bins,
    )
    net_imag = np.bincount(
        bins,
        weights=oriented.imag,
        minlength=n_bins,
    )
    net_complex = net_real + 1.0j * net_imag
    net_abs = np.abs(net_complex)

    cancellation = np.zeros(n_bins, dtype=float)
    active = raw_abs_by_channel > 0.0
    cancellation[active] = 1.0 - (
        net_abs[active] / raw_abs_by_channel[active]
    )
    # Roundoff can create tiny excursions outside [0, 1].
    cancellation = np.clip(cancellation, 0.0, 1.0)

    pair_abs_sum = float(np.sum(absolute_weights))
    l1_net_sum = float(np.sum(net_abs))
    bch_norm = float(np.linalg.norm(net_complex))
    l1_cancellation = (
        1.0 - l1_net_sum / pair_abs_sum if pair_abs_sum > 0.0 else 0.0
    )

    # Independent QA against the repository evaluator's own implementation.
    evaluator_norm = float(evaluator.evaluate(pauli_order_indices))
    if not math.isclose(
        bch_norm,
        evaluator_norm,
        rel_tol=1.0e-11,
        abs_tol=1.0e-13,
    ):
        raise RuntimeError(
            "Channel reconstruction does not reproduce HFCommutatorEvaluator: "
            f"{bch_norm:.16e} != {evaluator_norm:.16e}."
        )

    return {
        "signs": signs,
        "counts": counts,
        "raw_abs_by_channel": raw_abs_by_channel,
        "net_complex_by_channel": net_complex,
        "net_abs_by_channel": net_abs,
        "cancellation_by_channel": cancellation,
        "pair_abs_sum": pair_abs_sum,
        "bch_norm": bch_norm,
        "l1_net_sum": l1_net_sum,
        "l1_channel_cancellation": float(l1_cancellation),
    }


def invariant_channel_metrics(
    evaluator: HFCommutatorEvaluator,
    raw_abs_by_channel: np.ndarray,
    counts: np.ndarray,
) -> dict[str, float | int]:
    """Order-invariant descriptors of BCH-channel collision opportunity."""
    number_pairs = int(evaluator.pair_amplitudes.size)
    active = counts > 0
    colliding = counts >= 2
    number_channels = int(np.sum(active))
    number_collision_channels = int(np.sum(colliding))
    pair_abs_sum = float(np.sum(np.abs(evaluator.pair_amplitudes)))

    if pair_abs_sum > 0.0:
        p = raw_abs_by_channel[active] / pair_abs_sum
        effective_channels = float(1.0 / np.sum(p * p))
        max_channel_weight_fraction = float(np.max(p))
        collision_weight_fraction = float(
            np.sum(raw_abs_by_channel[colliding]) / pair_abs_sum
        )
        weighted_mean_multiplicity = float(
            np.sum(raw_abs_by_channel * counts) / pair_abs_sum
        )
    else:
        effective_channels = 0.0
        max_channel_weight_fraction = 0.0
        collision_weight_fraction = 0.0
        weighted_mean_multiplicity = 0.0

    return {
        "number_noncommuting_bch_pairs": number_pairs,
        "number_bch_channels": number_channels,
        "number_collision_channels": number_collision_channels,
        "collision_channel_fraction": safe_ratio(
            number_collision_channels,
            number_channels,
        ),
        "collision_pair_fraction": safe_ratio(
            int(np.sum(counts[colliding])),
            number_pairs,
        ),
        "collision_weight_fraction": collision_weight_fraction,
        "mean_pairs_per_channel": safe_ratio(number_pairs, number_channels),
        "max_pairs_in_channel": int(np.max(counts)) if number_channels else 0,
        "weighted_mean_channel_multiplicity": weighted_mean_multiplicity,
        "effective_weighted_channel_count": effective_channels,
        "effective_channel_fraction": safe_ratio(
            effective_channels,
            number_channels,
        ),
        "max_channel_weight_fraction": max_channel_weight_fraction,
    }


def ordering_metrics(
    prefix: str,
    arrays: dict[str, np.ndarray | float],
) -> dict[str, float]:
    raw_abs = np.asarray(arrays["raw_abs_by_channel"], dtype=float)
    cancellation = np.asarray(arrays["cancellation_by_channel"], dtype=float)
    pair_abs_sum = float(arrays["pair_abs_sum"])
    bch_norm = float(arrays["bch_norm"])
    colliding = raw_abs > 0.0

    # Restrict these descriptive summaries to channels that actually carry BCH
    # weight; weighted quantities remain the primary statistics.
    active_cancel = cancellation[colliding]
    weighted_cancel = (
        float(np.sum(raw_abs * cancellation) / pair_abs_sum)
        if pair_abs_sum > 0.0
        else 0.0
    )
    strong_50 = (
        float(np.sum(raw_abs[cancellation >= 0.50]) / pair_abs_sum)
        if pair_abs_sum > 0.0
        else 0.0
    )
    strong_90 = (
        float(np.sum(raw_abs[cancellation >= 0.90]) / pair_abs_sum)
        if pair_abs_sum > 0.0
        else 0.0
    )

    ratio = safe_ratio(bch_norm, pair_abs_sum)
    strength = -math.log10(max(ratio, 1.0e-300)) if math.isfinite(ratio) else math.nan

    return {
        f"{prefix}_pair_abs_sum": pair_abs_sum,
        f"{prefix}_bch2_hf_state_norm": bch_norm,
        f"{prefix}_bch_norm_over_pair_abs_sum": ratio,
        f"{prefix}_minus_log10_bch_norm_over_pair_abs_sum": strength,
        f"{prefix}_l1_channel_cancellation": float(
            arrays["l1_channel_cancellation"]
        ),
        f"{prefix}_weighted_channel_cancellation": weighted_cancel,
        f"{prefix}_median_channel_cancellation": (
            float(np.median(active_cancel)) if active_cancel.size else 0.0
        ),
        f"{prefix}_strong50_cancel_weight_fraction": strong_50,
        f"{prefix}_strong90_cancel_weight_fraction": strong_90,
    }


def build_case_rows(
    performance_row: pd.Series,
    tolerance: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    tensor_path = Path(str(performance_row["tensor_path"]))
    if not tensor_path.exists():
        raise FileNotFoundError(
            f"Tensor path from performance CSV does not exist: {tensor_path}"
        )

    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)
    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        tolerance,
    )
    jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    jw_hamiltonian.compress(abs_tol=tolerance)
    final_coefficients = {
        key: coefficient
        for key, coefficient in jw_hamiltonian.terms.items()
        if key != () and abs(coefficient) > tolerance
    }
    raw_pauli_keys = list(final_coefficients)
    if not raw_pauli_keys:
        raise ValueError(f"No identity-free JW Pauli terms for {metadata.case_id}.")

    orderings = baseline.build_deterministic_orderings(
        fermion_hamiltonian=fermion_hamiltonian,
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=tolerance,
    )
    raw_index_by_key = {key: index for index, key in enumerate(raw_pauli_keys)}
    signed_indices = [
        raw_index_by_key[key]
        for key in orderings[FERMIONIC_ORDER]
    ]
    jwmag_indices = [
        raw_index_by_key[key]
        for key in orderings[JW_MAGNITUDE_ORDER]
    ]

    pauli_graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    evaluator = HFCommutatorEvaluator.build(
        pauli_keys=raw_pauli_keys,
        coefficients=final_coefficients,
        pauli_graph=pauli_graph,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        tolerance=tolerance,
    )

    signed = channel_arrays(evaluator, signed_indices)
    jwmag = channel_arrays(evaluator, jwmag_indices)

    signed_raw_abs = np.asarray(signed["raw_abs_by_channel"], dtype=float)
    jwmag_raw_abs = np.asarray(jwmag["raw_abs_by_channel"], dtype=float)
    if not np.allclose(signed_raw_abs, jwmag_raw_abs, rtol=0.0, atol=1.0e-14):
        raise RuntimeError("Raw channel weights unexpectedly depend on ordering.")
    counts = np.asarray(signed["counts"], dtype=np.int64)

    signed_signs = np.asarray(signed["signs"], dtype=float)
    jwmag_signs = np.asarray(jwmag["signs"], dtype=float)
    pair_weights = np.abs(evaluator.pair_amplitudes)
    pair_abs_sum = float(np.sum(pair_weights))
    flipped = signed_signs != jwmag_signs
    weighted_flip_fraction = (
        float(np.sum(pair_weights[flipped]) / pair_abs_sum)
        if pair_abs_sum > 0.0
        else 0.0
    )

    signed_net = np.asarray(signed["net_abs_by_channel"], dtype=float)
    jwmag_net = np.asarray(jwmag["net_abs_by_channel"], dtype=float)
    channel_scale = np.maximum(signed_raw_abs, 1.0)
    compare_tol = 100.0 * np.finfo(float).eps * channel_scale
    signed_better = signed_net < (jwmag_net - compare_tol)
    jwmag_better = jwmag_net < (signed_net - compare_tol)

    total_channel_weight = float(np.sum(signed_raw_abs))
    signed_better_weight_fraction = (
        float(np.sum(signed_raw_abs[signed_better]) / total_channel_weight)
        if total_channel_weight > 0.0
        else 0.0
    )
    jwmag_better_weight_fraction = (
        float(np.sum(signed_raw_abs[jwmag_better]) / total_channel_weight)
        if total_channel_weight > 0.0
        else 0.0
    )

    actual_jw = float(performance_row["jw_magnitude_one_minus_overlap"])
    actual_ferm = float(performance_row["fermionic_aware_one_minus_overlap"])
    actual_advantage = safe_ratio(actual_jw, actual_ferm)
    actual_log_advantage = (
        math.log10(actual_advantage)
        if actual_advantage > 0.0 and math.isfinite(actual_advantage)
        else math.nan
    )

    signed_l1_cancel = float(signed["l1_channel_cancellation"])
    jwmag_l1_cancel = float(jwmag["l1_channel_cancellation"])
    signed_bch_norm = float(signed["bch_norm"])
    jwmag_bch_norm = float(jwmag["bch_norm"])
    bch_advantage = safe_ratio(jwmag_bch_norm, signed_bch_norm)
    bch_log_advantage = (
        math.log10(bch_advantage)
        if bch_advantage > 0.0 and math.isfinite(bch_advantage)
        else math.nan
    )

    summary: dict[str, Any] = {
        "case_id": metadata.case_id,
        "tensor_path": str(tensor_path),
        "molecule": metadata.molecule,
        "bond_length": metadata.bond_length,
        "basis": metadata.basis,
        "active_occupied": metadata.active_occupied,
        "active_vacant": metadata.active_vacant,
        "n_qubits": metadata.n_qubits,
        "number_pauli_terms": len(raw_pauli_keys),
        "pauli_noncommutation_edges": pauli_graph.number_of_edges(),
        "actual_jw_magnitude_one_minus_overlap": actual_jw,
        "actual_fermionic_one_minus_overlap": actual_ferm,
        "actual_jwmag_over_fermionic_advantage": actual_advantage,
        "actual_log10_jwmag_over_fermionic_advantage": actual_log_advantage,
        "actual_outcome": performance_row.get("outcome", ""),
        "weighted_pair_orientation_flip_fraction_signed_vs_jwmag": (
            weighted_flip_fraction
        ),
        "signed_minus_jwmag_l1_channel_cancellation": (
            signed_l1_cancel - jwmag_l1_cancel
        ),
        "signed_better_channel_weight_fraction": signed_better_weight_fraction,
        "jwmag_better_channel_weight_fraction": jwmag_better_weight_fraction,
        "signed_over_jwmag_bch_norm": safe_ratio(
            signed_bch_norm,
            jwmag_bch_norm,
        ),
        "jwmag_over_signed_bch_norm_advantage": bch_advantage,
        "log10_jwmag_over_signed_bch_norm_advantage": bch_log_advantage,
        "coefficient_tolerance": tolerance,
    }
    summary.update(invariant_channel_metrics(evaluator, signed_raw_abs, counts))
    summary.update(ordering_metrics("signed", signed))
    summary.update(ordering_metrics("jwmag", jwmag))

    signed_cancel = np.asarray(signed["cancellation_by_channel"], dtype=float)
    jwmag_cancel = np.asarray(jwmag["cancellation_by_channel"], dtype=float)
    active_channels = np.flatnonzero(counts > 0)
    ranked = sorted(
        active_channels.tolist(),
        key=lambda channel: (-signed_raw_abs[channel], channel),
    )
    rank_by_channel = {channel: rank + 1 for rank, channel in enumerate(ranked)}

    channel_rows: list[dict[str, Any]] = []
    for channel in active_channels:
        channel_rows.append(
            {
                "case_id": metadata.case_id,
                "channel_bin": int(channel),
                "channel_rank_by_raw_abs_weight": rank_by_channel[int(channel)],
                "pair_count": int(counts[channel]),
                "raw_abs_weight": float(signed_raw_abs[channel]),
                "raw_abs_weight_fraction": safe_ratio(
                    float(signed_raw_abs[channel]),
                    total_channel_weight,
                ),
                "is_collision_channel": bool(counts[channel] >= 2),
                "signed_net_abs": float(signed_net[channel]),
                "jwmag_net_abs": float(jwmag_net[channel]),
                "signed_channel_cancellation": float(signed_cancel[channel]),
                "jwmag_channel_cancellation": float(jwmag_cancel[channel]),
                "signed_minus_jwmag_channel_cancellation": float(
                    signed_cancel[channel] - jwmag_cancel[channel]
                ),
                "signed_better_than_jwmag": bool(signed_better[channel]),
                "jwmag_better_than_signed": bool(jwmag_better[channel]),
            }
        )

    return summary, channel_rows


def correlation_rows(summary: pd.DataFrame) -> pd.DataFrame:
    tests = (
        (
            "collision_weight_fraction",
            "signed_minus_jwmag_l1_channel_cancellation",
            "Does collision opportunity enable signed-order cancellation leverage?",
        ),
        (
            "signed_minus_jwmag_l1_channel_cancellation",
            "actual_log10_jwmag_over_fermionic_advantage",
            "Primary structural-mechanism test",
        ),
        (
            "log10_jwmag_over_signed_bch_norm_advantage",
            "actual_log10_jwmag_over_fermionic_advantage",
            "Leading BCH mechanism check",
        ),
        (
            "signed_better_channel_weight_fraction",
            "actual_log10_jwmag_over_fermionic_advantage",
            "Do wins improve more high-weight BCH channels?",
        ),
        (
            "weighted_pair_orientation_flip_fraction_signed_vs_jwmag",
            "signed_minus_jwmag_l1_channel_cancellation",
            "Are orientation changes sufficient, or must they be favorable?",
        ),
    )
    rows: list[dict[str, Any]] = []
    for x_name, y_name, purpose in tests:
        x = pd.to_numeric(summary[x_name], errors="coerce")
        y = pd.to_numeric(summary[y_name], errors="coerce")
        mask = np.isfinite(x.to_numpy()) & np.isfinite(y.to_numpy())
        xv = x.to_numpy(dtype=float)[mask]
        yv = y.to_numpy(dtype=float)[mask]
        row: dict[str, Any] = {
            "x": x_name,
            "y": y_name,
            "purpose": purpose,
            "n": int(xv.size),
            "pearson_r": math.nan,
            "pearson_p": math.nan,
            "spearman_rho": math.nan,
            "spearman_p": math.nan,
        }
        if xv.size >= 3 and np.std(xv) > 0.0 and np.std(yv) > 0.0:
            pearson = stats.pearsonr(xv, yv)
            spearman = stats.spearmanr(xv, yv)
            row.update(
                {
                    "pearson_r": float(pearson.statistic),
                    "pearson_p": float(pearson.pvalue),
                    "spearman_rho": float(spearman.statistic),
                    "spearman_p": float(spearman.pvalue),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_arguments()
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    if args.top_channels < 0:
        raise ValueError("--top-channels cannot be negative.")
    if not args.input.exists():
        raise FileNotFoundError(args.input)

    performance = pd.read_csv(args.input)
    required = {
        "case_id",
        "tensor_path",
        "valid_comparison",
        "jw_magnitude_one_minus_overlap",
        "fermionic_aware_one_minus_overlap",
    }
    missing = required.difference(performance.columns)
    if missing:
        raise ValueError(f"Performance CSV missing columns: {sorted(missing)}")

    performance = performance[as_bool(performance["valid_comparison"])].copy()
    requested_cases = list(args.cases) if args.cases else list(DEFAULT_CASES)

    selected_rows: list[pd.Series] = []
    for case_id in requested_cases:
        matches = performance[performance["case_id"] == case_id]
        if len(matches) != 1:
            raise ValueError(
                f"Expected exactly one valid row for {case_id!r}; found {len(matches)}."
            )
        selected_rows.append(matches.iloc[0])

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, Any]] = []
    all_channel_rows: list[dict[str, Any]] = []

    print("BCH channel-structure case study")
    print("================================")
    print(f"Cases: {len(selected_rows)}")
    print(f"Output: {args.outdir}")

    for index, performance_row in enumerate(selected_rows, start=1):
        case_id = str(performance_row["case_id"])
        print(f"[{index}/{len(selected_rows)}] {case_id}")
        summary, channels = build_case_rows(performance_row, args.tolerance)
        summary_rows.append(summary)
        all_channel_rows.extend(channels)
        print(
            "    actual JWmag/fermionic = "
            f"{summary['actual_jwmag_over_fermionic_advantage']:.6g}; "
            "BCH JWmag/signed = "
            f"{summary['jwmag_over_signed_bch_norm_advantage']:.6g}; "
            "Delta L1 cancellation = "
            f"{summary['signed_minus_jwmag_l1_channel_cancellation']:+.6f}; "
            "collision weight = "
            f"{summary['collision_weight_fraction']:.6f}"
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_path = args.outdir / "bch_channel_case_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    channel_df = pd.DataFrame(all_channel_rows)
    channel_path = args.outdir / "bch_channel_details.csv"
    channel_df.to_csv(channel_path, index=False)

    if args.top_channels > 0 and not channel_df.empty:
        top_df = (
            channel_df.sort_values(
                ["case_id", "channel_rank_by_raw_abs_weight"],
                kind="stable",
            )
            .groupby("case_id", sort=False)
            .head(args.top_channels)
            .reset_index(drop=True)
        )
    else:
        top_df = channel_df.iloc[0:0].copy()
    top_path = args.outdir / "bch_channel_top_channels.csv"
    top_df.to_csv(top_path, index=False)

    correlations = correlation_rows(summary_df)
    correlation_path = args.outdir / "bch_channel_correlations.csv"
    correlations.to_csv(correlation_path, index=False)

    display_columns = [
        "case_id",
        "actual_jwmag_over_fermionic_advantage",
        "collision_weight_fraction",
        "signed_l1_channel_cancellation",
        "jwmag_l1_channel_cancellation",
        "signed_minus_jwmag_l1_channel_cancellation",
        "jwmag_over_signed_bch_norm_advantage",
        "signed_better_channel_weight_fraction",
    ]
    print()
    print("Case summary")
    print("------------")
    print(summary_df[display_columns].to_string(index=False))

    print()
    print("Correlations")
    print("------------")
    print(
        correlations[
            ["purpose", "n", "pearson_r", "pearson_p", "spearman_rho", "spearman_p"]
        ].to_string(index=False)
    )

    print()
    print("Wrote:")
    print(f"  {summary_path}")
    print(f"  {channel_path}")
    print(f"  {top_path}")
    print(f"  {correlation_path}")


if __name__ == "__main__":
    main()
