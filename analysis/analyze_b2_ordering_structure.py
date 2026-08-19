#!/usr/bin/env python3
"""Diagnose why B2 fermionic coloring succeeds or fails by active space.

Run this file from the QHAT repository root after placing it in ``analysis/``.
It reuses the exact graph and ordering builders from the ``L-sweep`` branch and
reports, for selected B2/STO-6G active spaces:

* fermionic-term and final JW Pauli-term counts;
* graph vertices, edges, density, colors, and exact color-group sizes;
* the complete Pauli ordering produced by raw JW, JW coloring, and fermionic
  coloring;
* an order-independent coefficient-weighted noncommutation sum;
* two order-dependent leading-BCH commutator measures:
  - Pauli-coefficient l2 norm of sum_{j<k} [H_j, H_k];
  - norm of the same operator acting on the Hartree-Fock state.

The order-dependent measures are the important diagnostics for first-order
Trotter error. A simple edge count or sum of absolute pair weights cannot
capture cancellation caused by changing the term order.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence

import networkx as nx
from openfermion import get_fermion_operator, jordan_wigner

try:
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_fermionic_noncommutation_graph,
        build_hermitian_fermion_terms,
        build_orderings,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
        pauli_keys_commute,
        real_coefficient,
    )
except ImportError:
    from benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_fermionic_noncommutation_graph,
        build_hermitian_fermion_terms,
        build_orderings,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        discover_tensor_paths,
        load_interaction_operator,
        parse_case_metadata,
        pauli_keys_commute,
        real_coefficient,
    )


PauliKey = tuple[tuple[int, str], ...]
ORDERING_NAMES = ("jw_raw", "jw_coloring", "fermionic_coloring")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare graph structure, exact color-group sizes, generated "
            "orders, and leading commutator measures for B2 active spaces."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/b2_active_space_library"),
        help="Root containing B2 *.tensors.npz files.",
    )
    parser.add_argument(
        "--qubits",
        type=int,
        nargs="+",
        default=[12, 16, 18],
        help="Active-space sizes to compare. Default: 12 16 18.",
    )
    parser.add_argument(
        "--bond-length",
        default="1.70",
        help="Select one B2 bond-length directory. Default: 1.70.",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("analysis/b2_ordering_structure"),
        help="Directory for summary, group, and ordering files.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Coefficient and commutator compression tolerance.",
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.library.exists():
        raise FileNotFoundError(args.library)
    if not args.qubits or any(value <= 0 for value in args.qubits):
        raise ValueError("--qubits must contain positive integers.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")


def select_cases(
    library: Path,
    requested_qubits: Sequence[int],
    bond_length: str,
) -> list[Path]:
    requested = set(requested_qubits)
    selected: dict[int, list[Path]] = defaultdict(list)

    for path in discover_tensor_paths(library, case_pattern=None, limit=None):
        try:
            interaction, n_qubits = load_interaction_operator(path)
            del interaction
            metadata = parse_case_metadata(path, n_qubits)
        except Exception:
            continue

        if metadata.molecule != "B-B":
            continue
        if metadata.basis.lower() != "sto-6g":
            continue
        if str(metadata.bond_length) != str(bond_length):
            continue
        if n_qubits in requested:
            selected[n_qubits].append(path)

    missing = sorted(requested.difference(selected))
    if missing:
        raise FileNotFoundError(
            "No matching B2/STO-6G tensor files for active-space qubits "
            f"{missing} at bond length {bond_length}."
        )

    ambiguous = {
        n_qubits: paths
        for n_qubits, paths in selected.items()
        if len(paths) != 1
    }
    if ambiguous:
        details = "\n".join(
            f"  {n_qubits}: {[str(path) for path in paths]}"
            for n_qubits, paths in sorted(ambiguous.items())
        )
        raise ValueError(
            "Expected exactly one tensor file per active-space size:\n"
            f"{details}"
        )

    return [selected[n_qubits][0] for n_qubits in sorted(requested)]


def deterministic_color_groups(graph: nx.Graph) -> list[list[int]]:
    """Reproduce benchmark_L_sweep_trotter.deterministic_coloring_order."""
    coloring = nx.greedy_color(graph, strategy="largest_first")
    if not coloring:
        return []

    return [
        [node for node in graph.nodes if coloring[node] == color]
        for color in sorted(set(coloring.values()))
    ]


def graph_density(graph: nx.Graph) -> float:
    vertices = graph.number_of_nodes()
    if vertices < 2:
        return 0.0
    return 2.0 * graph.number_of_edges() / (vertices * (vertices - 1))


def group_statistics(groups: Sequence[Sequence[int]]) -> dict[str, Any]:
    sizes = [len(group) for group in groups]
    if not sizes:
        return {
            "group_sizes": [],
            "minimum_group_size": 0,
            "maximum_group_size": 0,
            "mean_group_size": 0.0,
            "group_size_std": 0.0,
        }

    mean = sum(sizes) / len(sizes)
    variance = sum((size - mean) ** 2 for size in sizes) / len(sizes)
    return {
        "group_sizes": sizes,
        "minimum_group_size": min(sizes),
        "maximum_group_size": max(sizes),
        "mean_group_size": mean,
        "group_size_std": math.sqrt(variance),
    }


def format_pauli_key(key: PauliKey) -> str:
    if not key:
        return "I"
    return " ".join(f"{label}{qubit}" for qubit, label in key)


_SINGLE_QUBIT_PRODUCT: dict[tuple[str, str], tuple[complex, str | None]] = {
    ("X", "X"): (1.0 + 0.0j, None),
    ("Y", "Y"): (1.0 + 0.0j, None),
    ("Z", "Z"): (1.0 + 0.0j, None),
    ("X", "Y"): (1.0j, "Z"),
    ("Y", "X"): (-1.0j, "Z"),
    ("Y", "Z"): (1.0j, "X"),
    ("Z", "Y"): (-1.0j, "X"),
    ("Z", "X"): (1.0j, "Y"),
    ("X", "Z"): (-1.0j, "Y"),
}


def multiply_pauli_keys(left: PauliKey, right: PauliKey) -> tuple[complex, PauliKey]:
    left_map = dict(left)
    right_map = dict(right)
    phase = 1.0 + 0.0j
    product: list[tuple[int, str]] = []

    for qubit in sorted(left_map.keys() | right_map.keys()):
        left_label = left_map.get(qubit)
        right_label = right_map.get(qubit)

        if left_label is None:
            product.append((qubit, right_label))
        elif right_label is None:
            product.append((qubit, left_label))
        else:
            local_phase, local_label = _SINGLE_QUBIT_PRODUCT[
                (left_label, right_label)
            ]
            phase *= local_phase
            if local_label is not None:
                product.append((qubit, local_label))

    return phase, tuple(product)


def apply_pauli_to_basis_index(
    pauli_key: PauliKey,
    basis_index: int,
    n_qubits: int,
) -> tuple[int, complex]:
    target = basis_index
    phase = 1.0 + 0.0j

    for qubit, label in pauli_key:
        bit_position = n_qubits - 1 - qubit
        bit_mask = 1 << bit_position
        occupied = bool(target & bit_mask)

        if label == "X":
            target ^= bit_mask
        elif label == "Y":
            phase *= -1.0j if occupied else 1.0j
            target ^= bit_mask
        elif label == "Z":
            if occupied:
                phase = -phase
        else:
            raise ValueError(f"Unsupported Pauli label {label!r}.")

    return target, phase


def hartree_fock_basis_index(n_qubits: int, n_electrons: int) -> int:
    return sum(
        1 << (n_qubits - 1 - qubit)
        for qubit in range(n_electrons)
    )


def leading_commutator_metrics(
    order: Sequence[PauliKey],
    coefficients: dict[PauliKey, complex],
    n_qubits: int,
    n_electrons: int,
    tolerance: float,
) -> dict[str, float | int]:
    """Compute matrix-free first-order BCH diagnostics for one term order.

    The accumulated operator is C2 = sum_{j<k} [H_j, H_k]. The Pauli l2 norm
    is sqrt(sum_P |alpha_P|^2) after combining equal Pauli products. It equals
    the Frobenius norm divided by sqrt(2**n_qubits).
    """
    real_coefficients = {
        key: real_coefficient(value, tolerance)
        for key, value in coefficients.items()
    }
    accumulated: dict[PauliKey, complex] = defaultdict(complex)
    pairwise_absolute_weight = 0.0
    noncommuting_pairs = 0

    for left_index, left_key in enumerate(order):
        left_coefficient = real_coefficients[left_key]
        for right_key in order[left_index + 1 :]:
            if pauli_keys_commute(left_key, right_key):
                continue

            right_coefficient = real_coefficients[right_key]
            product_phase, product_key = multiply_pauli_keys(
                left_key,
                right_key,
            )
            commutator_coefficient = (
                2.0
                * left_coefficient
                * right_coefficient
                * product_phase
            )
            accumulated[product_key] += commutator_coefficient
            pairwise_absolute_weight += abs(commutator_coefficient)
            noncommuting_pairs += 1

    accumulated = {
        key: value
        for key, value in accumulated.items()
        if abs(value) > tolerance
    }
    pauli_l2 = math.sqrt(
        sum(abs(value) ** 2 for value in accumulated.values())
    )

    hf_index = hartree_fock_basis_index(n_qubits, n_electrons)
    hf_amplitudes: dict[int, complex] = defaultdict(complex)
    for pauli_key, coefficient in accumulated.items():
        target, pauli_phase = apply_pauli_to_basis_index(
            pauli_key,
            hf_index,
            n_qubits,
        )
        hf_amplitudes[target] += coefficient * pauli_phase

    hf_state_norm = math.sqrt(
        sum(abs(value) ** 2 for value in hf_amplitudes.values())
    )

    return {
        "noncommuting_pauli_pairs": noncommuting_pairs,
        "pairwise_absolute_commutator_weight": pairwise_absolute_weight,
        "combined_commutator_pauli_terms": len(accumulated),
        "bch2_pauli_l2": pauli_l2,
        "bch2_hf_state_norm": hf_state_norm,
        "leading_hf_error_coefficient": 0.5 * hf_state_norm,
    }


def write_ordering_csv(
    output_path: Path,
    orderings: dict[str, Any],
    coefficients: dict[PauliKey, complex],
    tolerance: float,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = max(len(orderings[name].pauli_keys) for name in ORDERING_NAMES)
    fieldnames = ["position"]
    for name in ORDERING_NAMES:
        fieldnames.extend([f"{name}_pauli", f"{name}_coefficient"])

    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for position in range(rows):
            row: dict[str, Any] = {"position": position}
            for name in ORDERING_NAMES:
                keys = orderings[name].pauli_keys
                if position >= len(keys):
                    row[f"{name}_pauli"] = ""
                    row[f"{name}_coefficient"] = ""
                    continue
                key = keys[position]
                row[f"{name}_pauli"] = format_pauli_key(key)
                row[f"{name}_coefficient"] = real_coefficient(
                    coefficients[key],
                    tolerance,
                )
            writer.writerow(row)


def json_ready_fermionic_group(
    color: int,
    nodes: Sequence[int],
    hermitian_terms: Sequence[Any],
) -> dict[str, Any]:
    return {
        "color": color,
        "size": len(nodes),
        "nodes": list(nodes),
        "terms": [
            {
                "node": node,
                "first_raw_index": hermitian_terms[node].first_raw_index,
                "operator": str(hermitian_terms[node].operator),
            }
            for node in nodes
        ],
    }


def json_ready_jw_group(
    color: int,
    nodes: Sequence[int],
    raw_pauli_keys: Sequence[PauliKey],
    coefficients: dict[PauliKey, complex],
    tolerance: float,
) -> dict[str, Any]:
    return {
        "color": color,
        "size": len(nodes),
        "nodes": list(nodes),
        "terms": [
            {
                "node": node,
                "pauli": format_pauli_key(raw_pauli_keys[node]),
                "coefficient": real_coefficient(
                    coefficients[raw_pauli_keys[node]],
                    tolerance,
                ),
            }
            for node in nodes
        ],
    }


def analyze_case(
    tensor_path: Path,
    output_directory: Path,
    tolerance: float,
) -> list[dict[str, Any]]:
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)

    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        tolerance,
    )
    full_jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    full_jw_hamiltonian.compress(abs_tol=tolerance)

    orderings, hermitian_terms, final_coefficients = build_orderings(
        fermion_hamiltonian=fermion_hamiltonian,
        full_jw_hamiltonian=full_jw_hamiltonian,
        requested_orderings=ORDERING_NAMES,
        tolerance=tolerance,
        n_qubits=n_qubits,
    )
    if not hermitian_terms:
        hermitian_terms = build_hermitian_fermion_terms(
            fermion_hamiltonian,
            tolerance,
        )

    raw_pauli_keys = list(final_coefficients)
    jw_graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    fermionic_graph, _ = build_fermionic_noncommutation_graph(
        hermitian_terms,
        tolerance,
    )
    jw_groups = deterministic_color_groups(jw_graph)
    fermionic_groups = deterministic_color_groups(fermionic_graph)

    expected_jw_order = [node for group in jw_groups for node in group]
    actual_jw_order = [
        raw_pauli_keys.index(key)
        for key in orderings["jw_coloring"].pauli_keys
    ]
    if expected_jw_order != actual_jw_order:
        raise RuntimeError("Reconstructed JW color order does not match benchmark order.")

    case_directory = output_directory / metadata.case_id
    case_directory.mkdir(parents=True, exist_ok=True)

    group_payload = {
        "case_id": metadata.case_id,
        "tensor_path": str(tensor_path),
        "n_qubits": n_qubits,
        "active_occupied": metadata.active_occupied,
        "active_vacant": metadata.active_vacant,
        "jw_coloring": [
            json_ready_jw_group(
                color,
                nodes,
                raw_pauli_keys,
                final_coefficients,
                tolerance,
            )
            for color, nodes in enumerate(jw_groups)
        ],
        "fermionic_coloring": [
            json_ready_fermionic_group(
                color,
                nodes,
                hermitian_terms,
            )
            for color, nodes in enumerate(fermionic_groups)
        ],
    }
    with (case_directory / "color_groups.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(group_payload, stream, indent=2)

    write_ordering_csv(
        case_directory / "pauli_orderings.csv",
        orderings,
        final_coefficients,
        tolerance,
    )

    graph_by_ordering = {
        "jw_raw": None,
        "jw_coloring": jw_graph,
        "fermionic_coloring": fermionic_graph,
    }
    groups_by_ordering = {
        "jw_raw": [],
        "jw_coloring": jw_groups,
        "fermionic_coloring": fermionic_groups,
    }

    rows: list[dict[str, Any]] = []
    for ordering_name in ORDERING_NAMES:
        graph = graph_by_ordering[ordering_name]
        groups = groups_by_ordering[ordering_name]
        group_stats = group_statistics(groups)
        commutator_metrics = leading_commutator_metrics(
            order=orderings[ordering_name].pauli_keys,
            coefficients=final_coefficients,
            n_qubits=n_qubits,
            n_electrons=metadata.active_occupied,
            tolerance=tolerance,
        )

        rows.append(
            {
                "case_id": metadata.case_id,
                "tensor_path": str(tensor_path),
                "n_qubits": n_qubits,
                "active_occupied": metadata.active_occupied,
                "active_vacant": metadata.active_vacant,
                "number_of_fermionic_terms": len(hermitian_terms),
                "number_of_pauli_terms": len(final_coefficients),
                "ordering": ordering_name,
                "graph_vertices": (
                    graph.number_of_nodes() if graph is not None else ""
                ),
                "graph_edges": (
                    graph.number_of_edges() if graph is not None else ""
                ),
                "graph_density": (
                    graph_density(graph) if graph is not None else ""
                ),
                "number_of_colors": len(groups) if graph is not None else "",
                "group_sizes": json.dumps(group_stats["group_sizes"]),
                "minimum_group_size": (
                    group_stats["minimum_group_size"] if graph is not None else ""
                ),
                "maximum_group_size": (
                    group_stats["maximum_group_size"] if graph is not None else ""
                ),
                "mean_group_size": (
                    group_stats["mean_group_size"] if graph is not None else ""
                ),
                "group_size_std": (
                    group_stats["group_size_std"] if graph is not None else ""
                ),
                **commutator_metrics,
            }
        )

    with (case_directory / "summary.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(rows, stream, indent=2)

    return rows


def write_summary_csv(output_path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    args.output_directory.mkdir(parents=True, exist_ok=True)

    tensor_paths = select_cases(
        args.library,
        args.qubits,
        args.bond_length,
    )
    all_rows: list[dict[str, Any]] = []

    for index, tensor_path in enumerate(tensor_paths, start=1):
        print("=" * 88)
        print(f"[{index}/{len(tensor_paths)}] {tensor_path}")
        case_rows = analyze_case(
            tensor_path,
            args.output_directory,
            args.tolerance,
        )
        all_rows.extend(case_rows)
        for row in case_rows:
            print(
                f"  {row['ordering']:20s} "
                f"BCH2_Pauli_L2={row['bch2_pauli_l2']:.6e} "
                f"BCH2_HF_norm={row['bch2_hf_state_norm']:.6e}"
            )

    summary_path = args.output_directory / "b2_ordering_structure_summary.csv"
    write_summary_csv(summary_path, all_rows)
    print("=" * 88)
    print(f"Summary:   {summary_path}")
    print(f"Case files: {args.output_directory}")


if __name__ == "__main__":
    main()
