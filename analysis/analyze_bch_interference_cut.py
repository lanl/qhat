#!/usr/bin/env python3
"""Smoking-gun BCH interference-cut case study for fermionic ordering.

Goal
----
For each HF BCH output channel q, compare the JW descending-magnitude schedule
against the fermionic signed-coefficient schedule by splitting the JW-oriented
pair contributions into two groups:

    U_q = contributions whose pair orientation is unchanged
    F_q = contributions whose pair orientation flips

Then

    A_q(JW)     = U_q + F_q
    A_q(signed) = U_q - F_q

and therefore the channel-by-channel change in squared BCH state-error norm is
exactly

    |A_q(JW)|^2 - |A_q(signed)|^2 = 4 Re(conj(U_q) F_q).

Positive values mean the fermionic sign/orientation changes reduce the leading
BCH error in that channel. Negative values mean they destroy cancellation that
JW magnitude already had.

The script verifies this identity numerically, sums it over all orthogonal HF
output channels, and checks that the sum reproduces

    ||B_JW |HF>||^2 - ||B_signed |HF>||^2.

Run from the QHAT repository root, for example:

    python analysis/analyze_bch_interference_cut.py

The built-in panel is the same matched 16-qubit HGBS-5 panel used in the BCH
channel-structure case study: strong wins, near-neutral controls, and losses.
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
DEFAULT_OUTDIR = Path("analysis/bch_interference_cut")

# Fixed-size smoking-gun panel. Keep this panel unchanged for the first run.
DEFAULT_CASES = (
    "F-F_1.28_hgbs-5_as-014-002",       # very strong fermionic win
    "H2O_s-1.00_hgbs-5_as-006-010",    # strong fermionic win
    "NH3_s-1.00_hgbs-5_as-006-010",    # moderate fermionic win
    "O-O_1.26_hgbs-5_as-006-010",      # near neutral
    "C-C_1.50_hgbs-5_as-006-010",      # near neutral
    "B-B_1.70_hgbs-5_as-006-010",      # fermionic loss
    "CH4_s-1.00_hgbs-5_as-006-010",    # strong fermionic loss
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
            "Exact case_id from fermionic_aware_case_performance.csv. Repeat "
            "for multiple cases. If omitted, use the fixed seven-case panel."
        ),
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Coefficient/commutator tolerance.",
    )
    parser.add_argument(
        "--top-channels",
        type=int,
        default=25,
        help=(
            "Write this many channels per case ranked by absolute interference-"
            "cut contribution. Default: 25."
        ),
    )
    return parser.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0.0 or not math.isfinite(denominator):
        return math.nan
    return float(numerator / denominator)


def safe_log10_ratio(numerator: float, denominator: float) -> float:
    ratio = safe_ratio(numerator, denominator)
    if ratio <= 0.0 or not math.isfinite(ratio):
        return math.nan
    return float(math.log10(ratio))


def pair_signs(
    evaluator: HFCommutatorEvaluator,
    pauli_order_indices: Sequence[int],
) -> np.ndarray:
    """Return +/-1 pair orientation for every precomputed unordered Pauli pair."""
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


def bincount_complex(
    bins: np.ndarray,
    values: np.ndarray,
    minlength: int,
) -> np.ndarray:
    real = np.bincount(bins, weights=values.real, minlength=minlength)
    imag = np.bincount(bins, weights=values.imag, minlength=minlength)
    return real + 1.0j * imag


def build_case_objects(
    performance_row: pd.Series,
    tolerance: float,
) -> tuple[Any, HFCommutatorEvaluator, list[int], list[int], int, int]:
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

    return (
        metadata,
        evaluator,
        signed_indices,
        jwmag_indices,
        len(raw_pauli_keys),
        pauli_graph.number_of_edges(),
    )


def decompose_interference_cut(
    evaluator: HFCommutatorEvaluator,
    signed_indices: Sequence[int],
    jwmag_indices: Sequence[int],
) -> dict[str, np.ndarray | float | int]:
    """Compute U/F decomposition relative to the JW-magnitude pair orientation."""
    n_bins = evaluator.number_of_target_bins
    pair_amplitudes = np.asarray(evaluator.pair_amplitudes, dtype=np.complex128)
    bins = np.asarray(evaluator.target_bins, dtype=np.int64)

    signed_sign = pair_signs(evaluator, signed_indices)
    jwmag_sign = pair_signs(evaluator, jwmag_indices)
    flipped = signed_sign != jwmag_sign
    unchanged = ~flipped

    # Orient every unordered pair according to JW magnitude.  Signed ordering
    # keeps unchanged terms and negates exactly the flipped terms.
    jw_oriented_pairs = pair_amplitudes * jwmag_sign

    u_pairs = np.where(unchanged, jw_oriented_pairs, 0.0 + 0.0j)
    f_pairs = np.where(flipped, jw_oriented_pairs, 0.0 + 0.0j)

    U = bincount_complex(bins, u_pairs, n_bins)
    F = bincount_complex(bins, f_pairs, n_bins)
    jw_channel = U + F
    signed_channel = U - F

    jw_sq = np.abs(jw_channel) ** 2
    signed_sq = np.abs(signed_channel) ** 2
    delta_direct = jw_sq - signed_sq
    delta_cut = 4.0 * np.real(np.conj(U) * F)

    # Channel occupancy and pair-weight descriptors.
    pair_abs = np.abs(pair_amplitudes)
    counts = np.bincount(bins, minlength=n_bins).astype(np.int64)
    flipped_counts = np.bincount(
        bins,
        weights=flipped.astype(np.int64),
        minlength=n_bins,
    ).astype(np.int64)
    unchanged_counts = counts - flipped_counts
    raw_abs = np.bincount(bins, weights=pair_abs, minlength=n_bins)
    flipped_abs = np.bincount(
        bins,
        weights=pair_abs * flipped.astype(float),
        minlength=n_bins,
    )
    unchanged_abs = raw_abs - flipped_abs

    max_channel_identity_residual = float(
        np.max(np.abs(delta_direct - delta_cut)) if n_bins else 0.0
    )

    jw_norm = float(np.linalg.norm(jw_channel))
    signed_norm = float(np.linalg.norm(signed_channel))
    evaluator_jw_norm = float(evaluator.evaluate(jwmag_indices))
    evaluator_signed_norm = float(evaluator.evaluate(signed_indices))

    for label, reconstructed, evaluator_value in (
        ("JW magnitude", jw_norm, evaluator_jw_norm),
        ("fermionic signed", signed_norm, evaluator_signed_norm),
    ):
        if not math.isclose(
            reconstructed,
            evaluator_value,
            rel_tol=1.0e-11,
            abs_tol=1.0e-13,
        ):
            raise RuntimeError(
                f"{label} channel reconstruction failed: "
                f"{reconstructed:.16e} != {evaluator_value:.16e}."
            )

    total_delta_direct = float(np.sum(delta_direct))
    total_delta_cut = float(np.sum(delta_cut))
    total_norm_sq_difference = float(jw_norm * jw_norm - signed_norm * signed_norm)

    identity_scale = max(
        abs(total_delta_direct),
        abs(total_delta_cut),
        abs(total_norm_sq_difference),
        1.0,
    )
    if max_channel_identity_residual > 5.0e-11 * identity_scale:
        raise RuntimeError(
            "Per-channel interference identity failed; max residual = "
            f"{max_channel_identity_residual:.6e}."
        )
    if not math.isclose(
        total_delta_cut,
        total_norm_sq_difference,
        rel_tol=5.0e-11,
        abs_tol=5.0e-13,
    ):
        raise RuntimeError(
            "Summed cut interference does not reproduce BCH norm-square "
            "difference: "
            f"{total_delta_cut:.16e} != {total_norm_sq_difference:.16e}."
        )

    pair_abs_sum = float(np.sum(pair_abs))
    weighted_flip_fraction = (
        float(np.sum(pair_abs[flipped]) / pair_abs_sum)
        if pair_abs_sum > 0.0
        else 0.0
    )

    return {
        "U": U,
        "F": F,
        "jw_channel": jw_channel,
        "signed_channel": signed_channel,
        "jw_sq": jw_sq,
        "signed_sq": signed_sq,
        "delta_direct": delta_direct,
        "delta_cut": delta_cut,
        "counts": counts,
        "flipped_counts": flipped_counts,
        "unchanged_counts": unchanged_counts,
        "raw_abs": raw_abs,
        "flipped_abs": flipped_abs,
        "unchanged_abs": unchanged_abs,
        "flipped": flipped,
        "weighted_flip_fraction": weighted_flip_fraction,
        "jw_norm": jw_norm,
        "signed_norm": signed_norm,
        "total_delta_cut": total_delta_cut,
        "total_norm_sq_difference": total_norm_sq_difference,
        "max_channel_identity_residual": max_channel_identity_residual,
        "number_pairs": int(pair_amplitudes.size),
    }


def build_case_rows(
    performance_row: pd.Series,
    tolerance: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    (
        metadata,
        evaluator,
        signed_indices,
        jwmag_indices,
        number_pauli_terms,
        number_noncommutation_edges,
    ) = build_case_objects(performance_row, tolerance)

    cut = decompose_interference_cut(evaluator, signed_indices, jwmag_indices)

    U = np.asarray(cut["U"], dtype=np.complex128)
    F = np.asarray(cut["F"], dtype=np.complex128)
    jw_channel = np.asarray(cut["jw_channel"], dtype=np.complex128)
    signed_channel = np.asarray(cut["signed_channel"], dtype=np.complex128)
    jw_sq = np.asarray(cut["jw_sq"], dtype=float)
    signed_sq = np.asarray(cut["signed_sq"], dtype=float)
    delta_direct = np.asarray(cut["delta_direct"], dtype=float)
    delta_cut = np.asarray(cut["delta_cut"], dtype=float)
    counts = np.asarray(cut["counts"], dtype=np.int64)
    flipped_counts = np.asarray(cut["flipped_counts"], dtype=np.int64)
    unchanged_counts = np.asarray(cut["unchanged_counts"], dtype=np.int64)
    raw_abs = np.asarray(cut["raw_abs"], dtype=float)
    flipped_abs = np.asarray(cut["flipped_abs"], dtype=float)
    unchanged_abs = np.asarray(cut["unchanged_abs"], dtype=float)

    active = counts > 0
    active_channels = np.flatnonzero(active)

    positive = delta_cut > 0.0
    negative = delta_cut < 0.0
    positive_cut_sum = float(np.sum(delta_cut[positive]))
    negative_cut_abs_sum = float(-np.sum(delta_cut[negative]))
    absolute_cut_sum = positive_cut_sum + negative_cut_abs_sum
    total_delta_cut = float(cut["total_delta_cut"])

    positive_cut_share = (
        positive_cut_sum / absolute_cut_sum if absolute_cut_sum > 0.0 else 0.5
    )
    net_cut_balance = (
        total_delta_cut / absolute_cut_sum if absolute_cut_sum > 0.0 else 0.0
    )

    abs_cut = np.abs(delta_cut)
    ranked = sorted(
        active_channels.tolist(),
        key=lambda q: (-abs_cut[q], q),
    )
    rank_by_channel = {q: rank + 1 for rank, q in enumerate(ranked)}
    top5_abs_cut_fraction = (
        float(np.sum(abs_cut[ranked[:5]]) / np.sum(abs_cut[active]))
        if active_channels.size and np.sum(abs_cut[active]) > 0.0
        else 0.0
    )
    top10_abs_cut_fraction = (
        float(np.sum(abs_cut[ranked[:10]]) / np.sum(abs_cut[active]))
        if active_channels.size and np.sum(abs_cut[active]) > 0.0
        else 0.0
    )

    actual_jw = float(performance_row["jw_magnitude_one_minus_overlap"])
    actual_ferm = float(performance_row["fermionic_aware_one_minus_overlap"])
    actual_advantage = safe_ratio(actual_jw, actual_ferm)
    actual_log_advantage = safe_log10_ratio(actual_jw, actual_ferm)

    jw_norm = float(cut["jw_norm"])
    signed_norm = float(cut["signed_norm"])
    bch_advantage = safe_ratio(jw_norm, signed_norm)
    bch_log_advantage = safe_log10_ratio(jw_norm, signed_norm)

    jw_sq_total = jw_norm * jw_norm
    signed_sq_total = signed_norm * signed_norm
    normalized_cut_score = safe_ratio(
        jw_sq_total - signed_sq_total,
        jw_sq_total + signed_sq_total,
    )
    jw_error_removed_fraction = safe_ratio(
        jw_sq_total - signed_sq_total,
        jw_sq_total,
    )

    summary: dict[str, Any] = {
        "case_id": metadata.case_id,
        "tensor_path": str(performance_row["tensor_path"]),
        "molecule": metadata.molecule,
        "bond_length": metadata.bond_length,
        "basis": metadata.basis,
        "active_occupied": metadata.active_occupied,
        "active_vacant": metadata.active_vacant,
        "n_qubits": metadata.n_qubits,
        "number_pauli_terms": number_pauli_terms,
        "number_noncommutation_edges": number_noncommutation_edges,
        "number_bch_pairs": int(cut["number_pairs"]),
        "number_active_bch_channels": int(active_channels.size),
        "actual_jw_magnitude_one_minus_overlap": actual_jw,
        "actual_fermionic_one_minus_overlap": actual_ferm,
        "actual_jwmag_over_fermionic_advantage": actual_advantage,
        "actual_log10_jwmag_over_fermionic_advantage": actual_log_advantage,
        "actual_outcome": performance_row.get("outcome", ""),
        "jwmag_bch_norm": jw_norm,
        "signed_bch_norm": signed_norm,
        "jwmag_over_signed_bch_norm_advantage": bch_advantage,
        "log10_jwmag_over_signed_bch_norm_advantage": bch_log_advantage,
        "jwmag_bch_norm_squared": jw_sq_total,
        "signed_bch_norm_squared": signed_sq_total,
        "bch_norm_squared_difference_jw_minus_signed": jw_sq_total - signed_sq_total,
        "summed_interference_cut": total_delta_cut,
        "normalized_interference_cut_score": normalized_cut_score,
        "jw_bch_error_removed_fraction": jw_error_removed_fraction,
        "positive_interference_cut_sum": positive_cut_sum,
        "negative_interference_cut_abs_sum": negative_cut_abs_sum,
        "positive_interference_cut_share": positive_cut_share,
        "net_interference_cut_balance": net_cut_balance,
        "number_beneficial_channels": int(np.sum(positive & active)),
        "number_harmful_channels": int(np.sum(negative & active)),
        "number_zero_cut_channels": int(np.sum((delta_cut == 0.0) & active)),
        "top5_abs_interference_cut_fraction": top5_abs_cut_fraction,
        "top10_abs_interference_cut_fraction": top10_abs_cut_fraction,
        "weighted_pair_orientation_flip_fraction_signed_vs_jwmag": float(
            cut["weighted_flip_fraction"]
        ),
        "max_channel_identity_residual": float(
            cut["max_channel_identity_residual"]
        ),
        "total_identity_residual": float(
            total_delta_cut - float(cut["total_norm_sq_difference"])
        ),
        "coefficient_tolerance": tolerance,
    }

    channel_rows: list[dict[str, Any]] = []
    total_raw_abs = float(np.sum(raw_abs))
    total_abs_cut = float(np.sum(abs_cut[active]))

    for q in active_channels:
        raw_weight = float(raw_abs[q])
        cut_value = float(delta_cut[q])
        channel_rows.append(
            {
                "case_id": metadata.case_id,
                "channel_bin": int(q),
                "rank_by_abs_interference_cut": rank_by_channel[int(q)],
                "pair_count": int(counts[q]),
                "unchanged_pair_count": int(unchanged_counts[q]),
                "flipped_pair_count": int(flipped_counts[q]),
                "raw_abs_weight": raw_weight,
                "raw_abs_weight_fraction": safe_ratio(raw_weight, total_raw_abs),
                "unchanged_pair_abs_weight": float(unchanged_abs[q]),
                "flipped_pair_abs_weight": float(flipped_abs[q]),
                "channel_flipped_weight_fraction": safe_ratio(
                    float(flipped_abs[q]),
                    raw_weight,
                ),
                "U_real": float(U[q].real),
                "U_imag": float(U[q].imag),
                "U_abs": float(abs(U[q])),
                "F_real": float(F[q].real),
                "F_imag": float(F[q].imag),
                "F_abs": float(abs(F[q])),
                "jw_channel_real": float(jw_channel[q].real),
                "jw_channel_imag": float(jw_channel[q].imag),
                "jw_channel_abs": float(abs(jw_channel[q])),
                "signed_channel_real": float(signed_channel[q].real),
                "signed_channel_imag": float(signed_channel[q].imag),
                "signed_channel_abs": float(abs(signed_channel[q])),
                "jw_channel_norm_squared": float(jw_sq[q]),
                "signed_channel_norm_squared": float(signed_sq[q]),
                "direct_norm_squared_difference": float(delta_direct[q]),
                "interference_cut_4Re_UstarF": cut_value,
                "abs_interference_cut": float(abs_cut[q]),
                "abs_interference_cut_fraction": safe_ratio(
                    float(abs_cut[q]),
                    total_abs_cut,
                ),
                "fermionic_beneficial_channel": bool(cut_value > 0.0),
                "fermionic_harmful_channel": bool(cut_value < 0.0),
                "channel_identity_residual": float(delta_direct[q] - delta_cut[q]),
            }
        )

    return summary, channel_rows


def correlation_rows(summary: pd.DataFrame) -> pd.DataFrame:
    tests = (
        (
            "log10_jwmag_over_signed_bch_norm_advantage",
            "actual_log10_jwmag_over_fermionic_advantage",
            "BCH-to-finite-step bridge",
        ),
        (
            "normalized_interference_cut_score",
            "actual_log10_jwmag_over_fermionic_advantage",
            "Does net interference-cut direction track Trotter advantage?",
        ),
        (
            "positive_interference_cut_share",
            "actual_log10_jwmag_over_fermionic_advantage",
            "Are winning cases dominated by beneficial cut interference?",
        ),
        (
            "top5_abs_interference_cut_fraction",
            "actual_log10_jwmag_over_fermionic_advantage",
            "Is the ordering effect concentrated in a few channels?",
        ),
    )

    rows: list[dict[str, Any]] = []
    for x_name, y_name, purpose in tests:
        x = pd.to_numeric(summary[x_name], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(summary[y_name], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        xv = x[mask]
        yv = y[mask]
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

    print("BCH interference-cut smoking-gun case study")
    print("===========================================")
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
            "cut score = "
            f"{summary['normalized_interference_cut_score']:+.6f}; "
            "identity residual = "
            f"{summary['total_identity_residual']:.3e}"
        )

    summary_df = pd.DataFrame(summary_rows)
    channels_df = pd.DataFrame(all_channel_rows)
    correlations_df = correlation_rows(summary_df)

    summary_path = args.outdir / "bch_interference_cut_case_summary.csv"
    channels_path = args.outdir / "bch_interference_cut_channels.csv"
    top_path = args.outdir / "bch_interference_cut_top_channels.csv"
    correlation_path = args.outdir / "bch_interference_cut_correlations.csv"

    summary_df.to_csv(summary_path, index=False)
    channels_df.to_csv(channels_path, index=False)

    if args.top_channels > 0 and not channels_df.empty:
        top_df = (
            channels_df.sort_values(
                ["case_id", "rank_by_abs_interference_cut"],
                kind="stable",
            )
            .groupby("case_id", sort=False)
            .head(args.top_channels)
            .reset_index(drop=True)
        )
    else:
        top_df = channels_df.iloc[0:0].copy()
    top_df.to_csv(top_path, index=False)
    correlations_df.to_csv(correlation_path, index=False)

    display_columns = [
        "case_id",
        "actual_jwmag_over_fermionic_advantage",
        "jwmag_over_signed_bch_norm_advantage",
        "normalized_interference_cut_score",
        "positive_interference_cut_share",
        "top5_abs_interference_cut_fraction",
        "total_identity_residual",
    ]
    print()
    print("Case summary")
    print("------------")
    print(summary_df[display_columns].to_string(index=False))

    print()
    print("Correlations")
    print("------------")
    print(
        correlations_df[
            ["purpose", "n", "pearson_r", "pearson_p", "spearman_rho", "spearman_p"]
        ].to_string(index=False)
    )

    print()
    print("Wrote:")
    print(f"  {summary_path}")
    print(f"  {channels_path}")
    print(f"  {top_path}")
    print(f"  {correlation_path}")


if __name__ == "__main__":
    main()
