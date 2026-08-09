#!/usr/bin/env python3
"""Mechanism study for fermionic-informed first-order Trotter ordering.

This script is intentionally separate from the existing L-sweep benchmarks.
It keeps the final identity-free Jordan-Wigner Pauli Hamiltonian fixed and
changes only the order of those final Pauli exponentials.

The reference schedule is the current best deterministic rule:

    fermionic_signed_reference
        Complete Hermitian fermionic terms are ordered by increasing signed
        canonical coefficient, and the final combined JW Pauli order is
        induced by first occurrence.

The ablations isolate different pieces of that structure:

    signed_parent_blocks_randomized
        Build the reference parent-owned Pauli blocks once, then shuffle only
        the order of those whole blocks. Descendants stay contiguous and their
        internal order is unchanged. This isolates global/inter-parent order.

    signed_parent_within_randomized
        Keep the signed parent-block order fixed and shuffle only the Pauli
        descendants inside each reference block. This isolates internal
        descendant order while preserving parent contiguity.

    signed_parent_descendants_round_robin
        Keep the signed parent priority and the internal descendant order, but
        interleave one Pauli from each parent block in round-robin fashion.
        This deliberately destroys parent-block contiguity.

Two deterministic controls are also included:

    fermionic_magnitude_reference
        Complete Hermitian fermionic terms ordered by decreasing absolute
        canonical coefficient magnitude.

    jw_signed_baseline
        Final JW Pauli terms ordered directly by increasing signed coefficient.

Every schedule is evaluated on the same final Pauli set and coefficients.
Besides Trotter state error, the script records state-dependent BCH diagnostics,
a cancellation ratio, Pauli-pair orientation changes relative to the signed
fermionic reference, and parent-fragmentation metrics.

Run from the QHAT repository root after placing this file in analysis/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from openfermion import get_fermion_operator, jordan_wigner

try:
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_hermitian_fermion_terms,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        load_interaction_operator,
        parse_case_metadata,
        validate_pauli_order,
    )
    from qhat.analysis.benchmark_b2_active_spaces_matrix_free import (
        build_hartree_fock_state,
        compare_states,
        compile_ordered_terms,
        exact_reference_state,
        evolve_trotter_state,
        number_sector_basis_indices,
        warm_up_numba,
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
        validate_pauli_order,
    )
    from benchmark_b2_active_spaces_matrix_free import (
        build_hartree_fock_state,
        compare_states,
        compile_ordered_terms,
        exact_reference_state,
        evolve_trotter_state,
        number_sector_basis_indices,
        warm_up_numba,
    )
    from benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
        induced_pauli_order_indices,
        precompute_fermion_to_pauli_indices,
    )


PauliKey = tuple[tuple[int, str], ...]

REFERENCE = "fermionic_signed_reference"
DETERMINISTIC_SCHEDULES = (
    REFERENCE,
    "fermionic_magnitude_reference",
    "jw_signed_baseline",
    "signed_parent_descendants_round_robin",
)
RANDOM_SCHEDULES = (
    "signed_parent_blocks_randomized",
    "signed_parent_within_randomized",
)
ALL_SCHEDULES = DETERMINISTIC_SCHEDULES + RANDOM_SCHEDULES

FIELDNAMES = [
    "status",
    "error_message",
    "case_id",
    "tensor_path",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "number_of_fermionic_terms",
    "number_of_pauli_terms",
    "number_of_shared_pauli_terms",
    "shared_pauli_fraction",
    "number_of_unowned_pauli_terms",
    "schedule",
    "sample_index",
    "random_seed",
    "ordering_definition",
    "parent_order_hash",
    "pauli_order_hash",
    "reference_parent_bucket_sizes",
    "mean_parent_runs",
    "max_parent_runs",
    "fraction_parents_contiguous",
    "mean_parent_span_ratio",
    "weighted_parent_span_ratio",
    "orientation_flip_fraction_vs_signed_reference",
    "weighted_orientation_flip_fraction_vs_signed_reference",
    "cross_owner_weighted_flip_fraction_vs_signed_reference",
    "bch2_hf_state_norm",
    "bch_norm_ratio_to_signed_reference",
    "bch_pair_abs_sum",
    "bch_cancellation_ratio",
    "bch_cancellation_ratio_to_signed_reference",
    "bch_same_owner_norm",
    "bch_cross_owner_norm",
    "same_owner_pair_fraction",
    "same_owner_abs_weight_fraction",
    "unowned_pair_abs_weight_fraction",
    "leading_state_error_norm_estimate",
    "trotter_steps",
    "trotter_dt",
    "evolution_time",
    "nominal_exponential_count",
    "exact_sector_dimension",
    "exact_build_time_seconds",
    "exact_evolution_time_seconds",
    "trotter_runtime_seconds",
    "state_overlap_abs",
    "one_minus_overlap",
    "one_minus_overlap_ratio_to_signed_reference",
    "state_infidelity",
    "state_infidelity_ratio_to_signed_reference",
    "state_vector_2norm_error",
    "phase_aligned_state_2norm_error",
    "particle_number_leakage",
    "spin_sector_leakage",
    "coefficient_tolerance",
]

ResumeKey = tuple[str, str, int, int, int, float]


@dataclass(frozen=True)
class ScheduleSpec:
    name: str
    sample_index: int
    random_seed: int
    pauli_order_indices: tuple[int, ...]
    parent_order: tuple[int, ...]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Ablate fermionic parent ordering and Pauli-descendant structure "
            "while keeping the final JW Hamiltonian fixed."
        )
    )
    parser.add_argument(
        "--tensor",
        type=Path,
        action="append",
        required=True,
        help=(
            "Exact *.tensors.npz case path. Repeat --tensor to run multiple "
            "cases in one command."
        ),
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=[100],
        help="First-order Trotter step counts. Default: 100.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        dest="evolution_time",
        help="Total evolution time. Default: 1.0.",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=20,
        help="Samples for each randomized ablation family. Default: 20.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260807,
        help="Base random seed. Default: 20260807.",
    )
    parser.add_argument(
        "--schedules",
        nargs="+",
        choices=ALL_SCHEDULES,
        default=list(ALL_SCHEDULES),
        help="Subset of schedule families to evaluate.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/fermionic_structure_ablation_results.csv"),
        help="Append-only result CSV.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Coefficient and commutator tolerance.",
    )
    parser.add_argument(
        "--parallel-threshold",
        type=int,
        default=2**16,
        help="Use the parallel Pauli kernel at or above this vector size.",
    )
    parser.add_argument(
        "--no-spin-sector",
        action="store_true",
        help="Use only particle-number restriction for the exact reference.",
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.tensor:
        raise ValueError("At least one --tensor path is required.")
    for path in args.tensor:
        if not path.exists():
            raise FileNotFoundError(path)
    if not args.steps or any(value <= 0 for value in args.steps):
        raise ValueError("--steps must contain positive integers.")
    if args.evolution_time <= 0.0:
        raise ValueError("--time must be positive.")
    if args.samples < 0:
        raise ValueError("--samples cannot be negative.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    if args.parallel_threshold <= 0:
        raise ValueError("--parallel-threshold must be positive.")
    if len(args.schedules) != len(set(args.schedules)):
        raise ValueError("--schedules must not contain duplicates.")
    if REFERENCE not in args.schedules:
        raise ValueError(
            f"{REFERENCE!r} must be included because all ratios use it."
        )


def hash_integer_order(values: Sequence[int]) -> str:
    payload = ",".join(str(int(value)) for value in values)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def blank_row() -> dict[str, Any]:
    return {field: "" for field in FIELDNAMES}


def safe_ratio(numerator: float, denominator: float) -> str | float:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return ""
    if denominator <= 0.0:
        return ""
    return numerator / denominator


def resume_key(
    case_id: str,
    schedule: str,
    sample_index: int,
    random_seed: int,
    steps: int,
    evolution_time: float,
) -> ResumeKey:
    return (
        case_id,
        schedule,
        sample_index,
        random_seed,
        steps,
        evolution_time,
    )


def load_completed(output_path: Path) -> set[ResumeKey]:
    completed: set[ResumeKey] = set()
    if not output_path.exists() or output_path.stat().st_size == 0:
        return completed
    with output_path.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != FIELDNAMES:
            raise ValueError(
                "Existing output schema does not match this script: "
                f"{output_path}"
            )
        for row in reader:
            if row.get("status") != "success":
                continue
            try:
                completed.add(
                    resume_key(
                        row["case_id"],
                        row["schedule"],
                        int(row["sample_index"]),
                        int(row["random_seed"]),
                        int(row["trotter_steps"]),
                        float(row["evolution_time"]),
                    )
                )
            except (KeyError, TypeError, ValueError):
                continue
    return completed


def schedule_seed(
    base_seed: int,
    case_id: str,
    schedule: str,
    sample_index: int,
) -> int:
    schedule_code = {
        "signed_parent_blocks_randomized": 1,
        "signed_parent_within_randomized": 2,
    }[schedule]
    case_code = zlib.crc32(case_id.encode("utf-8")) & 0xFFFFFFFF
    sequence = np.random.SeedSequence(
        [base_seed, case_code, schedule_code, sample_index]
    )
    return int(sequence.generate_state(1, dtype=np.uint32)[0])


def build_reference_parent_buckets(
    parent_order: Sequence[int],
    fermion_to_pauli_indices: Sequence[Sequence[int]],
    number_of_pauli_terms: int,
) -> tuple[list[tuple[int, list[int]]], np.ndarray, list[int]]:
    """Assign each final Pauli string to its first parent in reference order.

    This mirrors the repository's first-occurrence induced-order convention.
    If a final Pauli string is reached by multiple complete fermionic terms, its
    owner is the first of those parents in the signed reference order.
    """
    owner_by_pauli = np.full(number_of_pauli_terms, -1, dtype=np.int32)
    buckets: list[tuple[int, list[int]]] = []

    for parent in parent_order:
        bucket: list[int] = []
        for pauli_index in fermion_to_pauli_indices[parent]:
            index = int(pauli_index)
            if owner_by_pauli[index] < 0:
                owner_by_pauli[index] = int(parent)
                bucket.append(index)
        if bucket:
            buckets.append((int(parent), bucket))

    fallback = [
        index
        for index in range(number_of_pauli_terms)
        if owner_by_pauli[index] < 0
    ]
    return buckets, owner_by_pauli, fallback


def flatten_reference_buckets(
    buckets: Sequence[tuple[int, Sequence[int]]],
    fallback: Sequence[int],
) -> list[int]:
    return [
        int(pauli)
        for _, bucket in buckets
        for pauli in bucket
    ] + [int(pauli) for pauli in fallback]


def randomized_block_order(
    buckets: Sequence[tuple[int, Sequence[int]]],
    fallback: Sequence[int],
    rng: np.random.Generator,
) -> tuple[list[int], list[int]]:
    order = list(range(len(buckets)))
    rng.shuffle(order)
    pauli_order = [
        int(pauli)
        for block_index in order
        for pauli in buckets[block_index][1]
    ]
    pauli_order.extend(int(pauli) for pauli in fallback)
    parent_order = [int(buckets[index][0]) for index in order]
    return pauli_order, parent_order


def randomized_within_blocks(
    buckets: Sequence[tuple[int, Sequence[int]]],
    fallback: Sequence[int],
    rng: np.random.Generator,
) -> tuple[list[int], list[int]]:
    pauli_order: list[int] = []
    parent_order: list[int] = []
    for parent, bucket in buckets:
        shuffled = list(map(int, bucket))
        rng.shuffle(shuffled)
        pauli_order.extend(shuffled)
        parent_order.append(int(parent))
    pauli_order.extend(int(pauli) for pauli in fallback)
    return pauli_order, parent_order


def round_robin_descendants(
    buckets: Sequence[tuple[int, Sequence[int]]],
    fallback: Sequence[int],
) -> tuple[list[int], list[int]]:
    """Interleave reference blocks while preserving order inside every block."""
    pauli_order: list[int] = []
    max_size = max((len(bucket) for _, bucket in buckets), default=0)
    for offset in range(max_size):
        for _, bucket in buckets:
            if offset < len(bucket):
                pauli_order.append(int(bucket[offset]))
    pauli_order.extend(int(pauli) for pauli in fallback)
    return pauli_order, [int(parent) for parent, _ in buckets]


def build_schedule_specs(
    case_id: str,
    args: argparse.Namespace,
    raw_pauli_keys: Sequence[PauliKey],
    final_coefficients: dict[PauliKey, complex],
    n_qubits: int,
    hermitian_terms: Sequence[Any],
    fermion_to_pauli_indices: Sequence[Sequence[int]],
    signed_parent_order: Sequence[int],
    magnitude_parent_order: Sequence[int],
    reference_buckets: Sequence[tuple[int, Sequence[int]]],
    fallback: Sequence[int],
) -> list[ScheduleSpec]:
    number_of_pauli_terms = len(raw_pauli_keys)
    raw_index_by_key = {
        key: index for index, key in enumerate(raw_pauli_keys)
    }

    reference_indices = induced_pauli_order_indices(
        fermionic_node_order=signed_parent_order,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        number_of_pauli_terms=number_of_pauli_terms,
    )
    reconstructed_reference = flatten_reference_buckets(
        reference_buckets,
        fallback,
    )
    if reference_indices != reconstructed_reference:
        raise RuntimeError(
            "Reference owner buckets do not reconstruct the repository's "
            "first-occurrence fermionic signed order."
        )

    magnitude_indices = induced_pauli_order_indices(
        fermionic_node_order=magnitude_parent_order,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        number_of_pauli_terms=number_of_pauli_terms,
    )

    jw_signed_keys = baseline.signed_coefficient_lexicographic_order(
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=args.tolerance,
    )
    jw_signed_indices = [raw_index_by_key[key] for key in jw_signed_keys]

    round_robin_indices, round_robin_parents = round_robin_descendants(
        reference_buckets,
        fallback,
    )

    deterministic = {
        REFERENCE: ScheduleSpec(
            name=REFERENCE,
            sample_index=-1,
            random_seed=-1,
            pauli_order_indices=tuple(reference_indices),
            parent_order=tuple(int(value) for value in signed_parent_order),
        ),
        "fermionic_magnitude_reference": ScheduleSpec(
            name="fermionic_magnitude_reference",
            sample_index=-1,
            random_seed=-1,
            pauli_order_indices=tuple(magnitude_indices),
            parent_order=tuple(int(value) for value in magnitude_parent_order),
        ),
        "jw_signed_baseline": ScheduleSpec(
            name="jw_signed_baseline",
            sample_index=-1,
            random_seed=-1,
            pauli_order_indices=tuple(jw_signed_indices),
            parent_order=tuple(),
        ),
        "signed_parent_descendants_round_robin": ScheduleSpec(
            name="signed_parent_descendants_round_robin",
            sample_index=-1,
            random_seed=-1,
            pauli_order_indices=tuple(round_robin_indices),
            parent_order=tuple(round_robin_parents),
        ),
    }

    specs: list[ScheduleSpec] = []
    for name in DETERMINISTIC_SCHEDULES:
        if name in args.schedules:
            specs.append(deterministic[name])

    for name in RANDOM_SCHEDULES:
        if name not in args.schedules:
            continue
        for sample_index in range(args.samples):
            seed = schedule_seed(
                args.seed,
                case_id,
                name,
                sample_index,
            )
            rng = np.random.default_rng(seed)
            if name == "signed_parent_blocks_randomized":
                indices, parent_order = randomized_block_order(
                    reference_buckets,
                    fallback,
                    rng,
                )
            elif name == "signed_parent_within_randomized":
                indices, parent_order = randomized_within_blocks(
                    reference_buckets,
                    fallback,
                    rng,
                )
            else:
                raise RuntimeError(f"Unhandled randomized schedule {name!r}.")
            specs.append(
                ScheduleSpec(
                    name=name,
                    sample_index=sample_index,
                    random_seed=seed,
                    pauli_order_indices=tuple(indices),
                    parent_order=tuple(parent_order),
                )
            )

    for spec in specs:
        validate_pauli_order(
            spec.name,
            [raw_pauli_keys[index] for index in spec.pauli_order_indices],
            raw_pauli_keys,
        )

    return specs


def ordering_definition(name: str) -> str:
    definitions = {
        REFERENCE: (
            "complete Hermitian fermionic terms ordered by increasing signed "
            "canonical coefficient; final combined JW Pauli order induced by "
            "first occurrence"
        ),
        "fermionic_magnitude_reference": (
            "complete Hermitian fermionic terms ordered by decreasing absolute "
            "canonical coefficient magnitude; final combined JW Pauli order "
            "induced by first occurrence"
        ),
        "jw_signed_baseline": (
            "final JW Pauli terms ordered directly by increasing signed "
            "coefficient with dense Pauli lexicographic tie-breaking"
        ),
        "signed_parent_blocks_randomized": (
            "reference signed-fermionic parent-owned Pauli blocks preserved "
            "internally and kept contiguous, but whole parent blocks randomly "
            "permuted"
        ),
        "signed_parent_within_randomized": (
            "reference signed-fermionic parent-block order kept fixed and "
            "contiguous, but Pauli descendants randomly permuted inside each "
            "parent-owned block"
        ),
        "signed_parent_descendants_round_robin": (
            "reference signed-fermionic parent priority and within-parent "
            "descendant order retained, but descendants interleaved round-robin "
            "across parent blocks to destroy contiguity"
        ),
    }
    return definitions[name]


def pair_signs(
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


def vector_norm_for_pair_mask(
    evaluator: HFCommutatorEvaluator,
    signed_amplitudes: np.ndarray,
    mask: np.ndarray,
) -> float:
    if signed_amplitudes.size == 0 or not np.any(mask):
        return 0.0
    bins = evaluator.target_bins[mask]
    values = signed_amplitudes[mask]
    real_parts = np.bincount(
        bins,
        weights=values.real,
        minlength=evaluator.number_of_target_bins,
    )
    imaginary_parts = np.bincount(
        bins,
        weights=values.imag,
        minlength=evaluator.number_of_target_bins,
    )
    return float(
        math.sqrt(
            float(
                np.dot(real_parts, real_parts)
                + np.dot(imaginary_parts, imaginary_parts)
            )
        )
    )


def commutator_diagnostics(
    evaluator: HFCommutatorEvaluator,
    pauli_order_indices: Sequence[int],
    reference_signs: np.ndarray,
    owner_by_pauli: np.ndarray,
) -> dict[str, float]:
    if evaluator.pair_amplitudes.size == 0:
        return {
            "bch2_hf_state_norm": 0.0,
            "bch_pair_abs_sum": 0.0,
            "bch_cancellation_ratio": 0.0,
            "orientation_flip_fraction_vs_signed_reference": 0.0,
            "weighted_orientation_flip_fraction_vs_signed_reference": 0.0,
            "cross_owner_weighted_flip_fraction_vs_signed_reference": 0.0,
            "bch_same_owner_norm": 0.0,
            "bch_cross_owner_norm": 0.0,
            "same_owner_pair_fraction": 0.0,
            "same_owner_abs_weight_fraction": 0.0,
            "unowned_pair_abs_weight_fraction": 0.0,
        }

    signs = pair_signs(evaluator, pauli_order_indices)
    signed_amplitudes = evaluator.pair_amplitudes * signs
    absolute_weights = np.abs(evaluator.pair_amplitudes)
    pair_abs_sum = float(np.sum(absolute_weights))
    bch_norm = evaluator.evaluate(pauli_order_indices)

    flipped = signs != reference_signs
    flip_fraction = float(np.mean(flipped))
    weighted_flip = (
        float(np.sum(absolute_weights[flipped])) / pair_abs_sum
        if pair_abs_sum > 0.0
        else 0.0
    )

    left_owner = owner_by_pauli[evaluator.left_indices]
    right_owner = owner_by_pauli[evaluator.right_indices]
    owned = (left_owner >= 0) & (right_owner >= 0)
    same_owner = owned & (left_owner == right_owner)
    cross_owner = owned & (left_owner != right_owner)
    unowned = ~owned

    same_owner_norm = vector_norm_for_pair_mask(
        evaluator,
        signed_amplitudes,
        same_owner,
    )
    cross_owner_norm = vector_norm_for_pair_mask(
        evaluator,
        signed_amplitudes,
        cross_owner,
    )

    same_owner_pair_fraction = float(np.mean(same_owner))
    same_owner_abs_weight_fraction = (
        float(np.sum(absolute_weights[same_owner])) / pair_abs_sum
        if pair_abs_sum > 0.0
        else 0.0
    )
    unowned_abs_weight_fraction = (
        float(np.sum(absolute_weights[unowned])) / pair_abs_sum
        if pair_abs_sum > 0.0
        else 0.0
    )
    cross_abs_sum = float(np.sum(absolute_weights[cross_owner]))
    cross_flipped = cross_owner & flipped
    cross_weighted_flip = (
        float(np.sum(absolute_weights[cross_flipped])) / cross_abs_sum
        if cross_abs_sum > 0.0
        else 0.0
    )

    return {
        "bch2_hf_state_norm": bch_norm,
        "bch_pair_abs_sum": pair_abs_sum,
        "bch_cancellation_ratio": (
            bch_norm / pair_abs_sum if pair_abs_sum > 0.0 else 0.0
        ),
        "orientation_flip_fraction_vs_signed_reference": flip_fraction,
        "weighted_orientation_flip_fraction_vs_signed_reference": weighted_flip,
        "cross_owner_weighted_flip_fraction_vs_signed_reference": (
            cross_weighted_flip
        ),
        "bch_same_owner_norm": same_owner_norm,
        "bch_cross_owner_norm": cross_owner_norm,
        "same_owner_pair_fraction": same_owner_pair_fraction,
        "same_owner_abs_weight_fraction": same_owner_abs_weight_fraction,
        "unowned_pair_abs_weight_fraction": unowned_abs_weight_fraction,
    }


def parent_fragmentation_metrics(
    pauli_order_indices: Sequence[int],
    owner_by_pauli: np.ndarray,
) -> dict[str, float]:
    positions = np.empty(len(pauli_order_indices), dtype=np.int32)
    positions[np.asarray(pauli_order_indices, dtype=np.int32)] = np.arange(
        len(pauli_order_indices),
        dtype=np.int32,
    )

    parents = sorted(
        int(value)
        for value in np.unique(owner_by_pauli)
        if int(value) >= 0
    )
    if not parents:
        return {
            "mean_parent_runs": 0.0,
            "max_parent_runs": 0.0,
            "fraction_parents_contiguous": 0.0,
            "mean_parent_span_ratio": 0.0,
            "weighted_parent_span_ratio": 0.0,
        }

    runs_values: list[int] = []
    span_ratios: list[float] = []
    sizes: list[int] = []
    for parent in parents:
        members = np.flatnonzero(owner_by_pauli == parent)
        if members.size == 0:
            continue
        member_positions = np.sort(positions[members])
        gaps = np.diff(member_positions)
        runs = 1 + int(np.count_nonzero(gaps > 1))
        span = int(member_positions[-1] - member_positions[0] + 1)
        size = int(member_positions.size)
        runs_values.append(runs)
        span_ratios.append(span / size)
        sizes.append(size)

    weights = np.asarray(sizes, dtype=float)
    ratios = np.asarray(span_ratios, dtype=float)
    return {
        "mean_parent_runs": float(np.mean(runs_values)),
        "max_parent_runs": float(np.max(runs_values)),
        "fraction_parents_contiguous": float(
            np.mean(np.asarray(runs_values) == 1)
        ),
        "mean_parent_span_ratio": float(np.mean(ratios)),
        "weighted_parent_span_ratio": float(
            np.average(ratios, weights=weights)
        ),
    }


def parent_multiplicity_statistics(
    fermion_to_pauli_indices: Sequence[Sequence[int]],
    number_of_pauli_terms: int,
) -> tuple[int, float]:
    multiplicity = np.zeros(number_of_pauli_terms, dtype=np.int32)
    for mapped in fermion_to_pauli_indices:
        if not mapped:
            continue
        unique = np.unique(np.asarray(mapped, dtype=np.int32))
        multiplicity[unique] += 1
    shared = int(np.count_nonzero(multiplicity > 1))
    return shared, shared / number_of_pauli_terms if number_of_pauli_terms else 0.0


def benchmark_case(
    tensor_path: Path,
    args: argparse.Namespace,
    completed: set[ResumeKey],
    writer: csv.DictWriter,
    output_stream: Any,
) -> None:
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)

    print()
    print("============================================================")
    print(f"Case: {metadata.case_id}")
    print(f"Tensor: {tensor_path}")
    print("============================================================")

    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        args.tolerance,
    )
    full_jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    full_jw_hamiltonian.compress(abs_tol=args.tolerance)
    final_coefficients = {
        key: coefficient
        for key, coefficient in full_jw_hamiltonian.terms.items()
        if key != () and abs(coefficient) > args.tolerance
    }
    raw_pauli_keys = list(final_coefficients)
    if not raw_pauli_keys:
        raise ValueError("The identity-free JW Hamiltonian has no Pauli terms.")
    raw_index_by_key = {
        key: index for index, key in enumerate(raw_pauli_keys)
    }

    print("  Building complete Hermitian fermionic terms...")
    hermitian_terms = build_hermitian_fermion_terms(
        fermion_hamiltonian,
        args.tolerance,
    )
    fermion_to_pauli_indices = precompute_fermion_to_pauli_indices(
        hermitian_terms=hermitian_terms,
        final_coefficients=final_coefficients,
        raw_index_by_key=raw_index_by_key,
        tolerance=args.tolerance,
    )

    signed_parent_order = baseline.fermionic_term_order_indices(
        hermitian_terms=hermitian_terms,
        ordering_method="signed_ascending",
        tolerance=args.tolerance,
    )
    magnitude_parent_order = baseline.fermionic_term_order_indices(
        hermitian_terms=hermitian_terms,
        ordering_method="magnitude_descending",
        tolerance=args.tolerance,
    )

    reference_buckets, owner_by_pauli, fallback = build_reference_parent_buckets(
        signed_parent_order,
        fermion_to_pauli_indices,
        len(raw_pauli_keys),
    )
    reference_bucket_sizes = [len(bucket) for _, bucket in reference_buckets]
    number_unowned = len(fallback)
    number_shared, shared_fraction = parent_multiplicity_statistics(
        fermion_to_pauli_indices,
        len(raw_pauli_keys),
    )

    print(
        "  Parent ownership: "
        f"{len(reference_buckets)} nonempty signed-reference blocks, "
        f"{number_shared} shared Pauli strings, "
        f"{number_unowned} unowned final Pauli strings."
    )

    specs = build_schedule_specs(
        case_id=metadata.case_id,
        args=args,
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        hermitian_terms=hermitian_terms,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        signed_parent_order=signed_parent_order,
        magnitude_parent_order=magnitude_parent_order,
        reference_buckets=reference_buckets,
        fallback=fallback,
    )

    reference_spec = next(spec for spec in specs if spec.name == REFERENCE)

    print("  Building Pauli noncommutation graph and BCH evaluator...")
    jw_graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    bch_evaluator = HFCommutatorEvaluator.build(
        pauli_keys=raw_pauli_keys,
        coefficients=final_coefficients,
        pauli_graph=jw_graph,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        tolerance=args.tolerance,
    )
    reference_signs = pair_signs(
        bch_evaluator,
        reference_spec.pauli_order_indices,
    )
    reference_commutator = commutator_diagnostics(
        bch_evaluator,
        reference_spec.pauli_order_indices,
        reference_signs,
        owner_by_pauli,
    )
    reference_bch_norm = float(reference_commutator["bch2_hf_state_norm"])
    reference_cancellation_ratio = float(
        reference_commutator["bch_cancellation_ratio"]
    )

    print("  Building exact symmetry-sector reference...")
    (
        exact_sector_state,
        exact_basis_indices,
        exact_build_time,
        exact_evolution_time,
    ) = exact_reference_state(
        fermion_hamiltonian=fermion_hamiltonian,
        n_qubits=n_qubits,
        n_electrons=metadata.active_occupied,
        evolution_time=args.evolution_time,
        tolerance=args.tolerance,
        spin_preserving=not args.no_spin_sector,
    )
    number_basis_indices = number_sector_basis_indices(
        n_qubits,
        metadata.active_occupied,
    )
    initial_state = build_hartree_fock_state(
        n_qubits,
        metadata.active_occupied,
    )
    compiled_raw_terms = compile_ordered_terms(
        raw_pauli_keys,
        final_coefficients,
        n_qubits,
        args.tolerance,
    )

    # Process the signed reference first for every step count so error ratios are
    # available immediately for the subsequent ablations.
    ordered_specs = [reference_spec] + [
        spec for spec in specs if spec.name != REFERENCE
    ]

    for steps in args.steps:
        reference_one_minus_overlap: float | None = None
        reference_infidelity: float | None = None

        # If a signed reference row already exists, recover its metrics so later
        # rows can still get ratio columns without rerunning it.
        if args.output.exists() and args.output.stat().st_size > 0:
            with args.output.open("r", newline="", encoding="utf-8") as stream:
                for old_row in csv.DictReader(stream):
                    if old_row.get("status") != "success":
                        continue
                    if old_row.get("case_id") != metadata.case_id:
                        continue
                    if old_row.get("schedule") != REFERENCE:
                        continue
                    try:
                        if int(old_row["trotter_steps"]) != steps:
                            continue
                        if not math.isclose(
                            float(old_row["evolution_time"]),
                            args.evolution_time,
                            rel_tol=0.0,
                            abs_tol=1.0e-14,
                        ):
                            continue
                        reference_one_minus_overlap = float(
                            old_row["one_minus_overlap"]
                        )
                        reference_infidelity = float(old_row["state_infidelity"])
                        break
                    except (KeyError, TypeError, ValueError):
                        continue

        for spec in ordered_specs:
            key = resume_key(
                metadata.case_id,
                spec.name,
                spec.sample_index,
                spec.random_seed,
                steps,
                args.evolution_time,
            )
            if key in completed:
                print(
                    f"  SKIP {spec.name} sample={spec.sample_index} "
                    f"steps={steps}"
                )
                continue

            print(
                f"  RUN  {spec.name} sample={spec.sample_index} steps={steps}"
            )
            row = blank_row()
            row.update(
                {
                    "status": "success",
                    "case_id": metadata.case_id,
                    "tensor_path": str(tensor_path),
                    "molecule": metadata.molecule,
                    "bond_length": metadata.bond_length,
                    "basis": metadata.basis,
                    "active_occupied": metadata.active_occupied,
                    "active_vacant": metadata.active_vacant,
                    "n_qubits": n_qubits,
                    "number_of_fermionic_terms": len(hermitian_terms),
                    "number_of_pauli_terms": len(raw_pauli_keys),
                    "number_of_shared_pauli_terms": number_shared,
                    "shared_pauli_fraction": shared_fraction,
                    "number_of_unowned_pauli_terms": number_unowned,
                    "schedule": spec.name,
                    "sample_index": spec.sample_index,
                    "random_seed": spec.random_seed,
                    "ordering_definition": ordering_definition(spec.name),
                    "parent_order_hash": (
                        hash_integer_order(spec.parent_order)
                        if spec.parent_order
                        else ""
                    ),
                    "pauli_order_hash": hash_integer_order(
                        spec.pauli_order_indices
                    ),
                    "reference_parent_bucket_sizes": json.dumps(
                        reference_bucket_sizes,
                        separators=(",", ":"),
                    ),
                    "trotter_steps": steps,
                    "trotter_dt": args.evolution_time / steps,
                    "evolution_time": args.evolution_time,
                    "exact_sector_dimension": exact_sector_state.size,
                    "exact_build_time_seconds": exact_build_time,
                    "exact_evolution_time_seconds": exact_evolution_time,
                    "coefficient_tolerance": args.tolerance,
                }
            )

            try:
                fragmentation = parent_fragmentation_metrics(
                    spec.pauli_order_indices,
                    owner_by_pauli,
                )
                row.update(fragmentation)

                commutator = commutator_diagnostics(
                    bch_evaluator,
                    spec.pauli_order_indices,
                    reference_signs,
                    owner_by_pauli,
                )
                row.update(commutator)
                bch_norm = float(commutator["bch2_hf_state_norm"])
                cancellation_ratio = float(
                    commutator["bch_cancellation_ratio"]
                )
                row["bch_norm_ratio_to_signed_reference"] = safe_ratio(
                    bch_norm,
                    reference_bch_norm,
                )
                row[
                    "bch_cancellation_ratio_to_signed_reference"
                ] = safe_ratio(
                    cancellation_ratio,
                    reference_cancellation_ratio,
                )
                row["leading_state_error_norm_estimate"] = (
                    0.5
                    * args.evolution_time**2
                    * bch_norm
                    / steps
                )

                ordered_terms = [
                    compiled_raw_terms[index]
                    for index in spec.pauli_order_indices
                ]
                trotter_start = time.perf_counter()
                approximate_state, nominal_exponential_count = (
                    evolve_trotter_state(
                        initial_state=initial_state,
                        terms=ordered_terms,
                        formula_order=1,
                        trotter_steps=steps,
                        evolution_time=args.evolution_time,
                        parallel_threshold=args.parallel_threshold,
                    )
                )
                row["trotter_runtime_seconds"] = (
                    time.perf_counter() - trotter_start
                )
                row["nominal_exponential_count"] = nominal_exponential_count

                metrics = compare_states(
                    exact_sector_state=exact_sector_state,
                    exact_basis_indices=exact_basis_indices,
                    approximate_full_state=approximate_state,
                    number_basis_indices=number_basis_indices,
                )
                for field in (
                    "state_overlap_abs",
                    "state_infidelity",
                    "state_vector_2norm_error",
                    "phase_aligned_state_2norm_error",
                    "particle_number_leakage",
                    "spin_sector_leakage",
                ):
                    row[field] = metrics.get(field, "")

                overlap = float(metrics["state_overlap_abs"])
                one_minus_overlap = max(0.0, 1.0 - overlap)
                infidelity = float(metrics["state_infidelity"])
                row["one_minus_overlap"] = one_minus_overlap

                if spec.name == REFERENCE:
                    reference_one_minus_overlap = one_minus_overlap
                    reference_infidelity = infidelity
                    row["one_minus_overlap_ratio_to_signed_reference"] = 1.0
                    row["state_infidelity_ratio_to_signed_reference"] = 1.0
                else:
                    if reference_one_minus_overlap is not None:
                        row[
                            "one_minus_overlap_ratio_to_signed_reference"
                        ] = safe_ratio(
                            one_minus_overlap,
                            reference_one_minus_overlap,
                        )
                    if reference_infidelity is not None:
                        row[
                            "state_infidelity_ratio_to_signed_reference"
                        ] = safe_ratio(
                            infidelity,
                            reference_infidelity,
                        )

                writer.writerow(row)
                output_stream.flush()
                completed.add(key)

                print(
                    "       "
                    f"1-overlap={one_minus_overlap:.6e}  "
                    f"BCH={bch_norm:.6e}  "
                    f"cancel={cancellation_ratio:.6e}  "
                    f"flip={float(row['weighted_orientation_flip_fraction_vs_signed_reference']):.3f}"
                )
            except Exception as error:
                row["status"] = "failed"
                row["error_message"] = f"{type(error).__name__}: {error}"
                writer.writerow(row)
                output_stream.flush()
                print(f"       FAILED: {type(error).__name__}: {error}")


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    warm_up_numba()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    completed = load_completed(args.output)
    write_header = not args.output.exists() or args.output.stat().st_size == 0

    with args.output.open("a", newline="", encoding="utf-8") as output_stream:
        writer = csv.DictWriter(output_stream, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
            output_stream.flush()

        for tensor_path in args.tensor:
            benchmark_case(
                tensor_path=tensor_path,
                args=args,
                completed=completed,
                writer=writer,
                output_stream=output_stream,
            )

    print()
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
