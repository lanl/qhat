#!/usr/bin/env python3
"""
Benchmark Trotter orderings for polyatomic geometry-sweep Hamiltonians.

This is the polyatomic companion to the L-sweep benchmark. It compares three
orderings of the same final, combined Jordan-Wigner Hamiltonian:

    jw_raw
        Raw OpenFermion Jordan-Wigner insertion order.

    jw_coloring
        Greedy coloring of the Pauli-string noncommutation graph.

    fermionic_coloring
        Greedy coloring of Hermitian fermionic terms, followed by the induced
        ordering of the final combined JW Pauli strings.

Expected tensor-library layout:

    polyatomic_library/<molecule>/s-<scale>/<basis>/*.tensors.npz

Each tensor file must contain ``constant``, ``one_body``, and ``two_body``.
The identity term is omitted from both the exact and product-formula
Hamiltonians because it contributes only a global phase. First-, second-, and
fourth-order Suzuki formulas are supported.

In addition to the original error columns, this version records atom count,
geometry scale, graph density, Pauli/fermionic term ratio, and direct error
ratios relative to JW coloring. A ratio below 1 means the current ordering is
more accurate than JW coloring for the same case and Trotter parameters.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import networkx as nx
import numpy as np
from scipy.linalg import expm

from openfermion import (
    FermionOperator,
    InteractionOperator,
    get_fermion_operator,
    jordan_wigner,
    normal_ordered,
)
from openfermion.utils import commutator, hermitian_conjugated


DEFAULT_TOLERANCE = 1.0e-12
SUPPORTED_ORDERINGS = (
    "jw_raw",
    "jw_coloring",
    "fermionic_coloring",
)

CASE_PATTERN = re.compile(
    r"^(?P<prefix>.+)_as-(?P<active_occupied>\d{3})-"
    r"(?P<active_vacant>\d{3})(?:_(?P<source_mapping>JW|BK))?"
    r"\.tensors\.npz$",
    re.IGNORECASE,
)
SCALE_PATTERN = re.compile(r"^s-(?P<scale>\d+(?:\.\d+)?)$", re.IGNORECASE)
FORMULA_TOKEN_PATTERN = re.compile(r"([A-Z][a-z]?)(\d*)")


@dataclass(frozen=True)
class CaseMetadata:
    case_id: str
    tensor_path: Path
    molecule: str
    atom_count: int
    geometry_tag: str
    geometry_scale: float
    basis: str
    source_mapping: str
    active_occupied: int
    active_vacant: int
    n_qubits: int


@dataclass
class HermitianFermionTerm:
    index: int
    first_raw_index: int
    operator: FermionOperator
    component_keys: tuple[tuple[tuple[int, int], ...], ...]


@dataclass
class OrderingResult:
    name: str
    pauli_keys: list[tuple[tuple[int, str], ...]]
    graph_vertices: int | None = None
    graph_edges: int | None = None
    number_of_colors: int | None = None
    graph_time_seconds: float | None = None
    coloring_time_seconds: float | None = None
    build_time_seconds: float | None = None


# -----------------------------------------------------------------------------
# Input discovery and metadata
# -----------------------------------------------------------------------------


def chemical_formula_atom_count(formula: str) -> int:
    """Return the total atom count encoded by a simple molecular formula."""
    tokens = FORMULA_TOKEN_PATTERN.findall(formula)
    reconstructed = "".join(element + count for element, count in tokens)
    if not tokens or reconstructed != formula:
        raise ValueError(
            f"Cannot infer atom count from molecule directory name {formula!r}."
        )
    return sum(int(count) if count else 1 for _, count in tokens)


def parse_case_metadata(tensor_path: Path, n_qubits: int) -> CaseMetadata:
    match = CASE_PATTERN.match(tensor_path.name)
    if match is None:
        raise ValueError(
            "Tensor filename does not match the expected "
            f"'*_as-XXX-XXX[ _JW].tensors.npz' format: {tensor_path.name}"
        )

    active_occupied = int(match.group("active_occupied"))
    active_vacant = int(match.group("active_vacant"))

    if active_occupied + active_vacant != n_qubits:
        raise ValueError(
            "Active-space counts do not match the tensor dimension: "
            f"{active_occupied} + {active_vacant} != {n_qubits} "
            f"for {tensor_path}"
        )

    try:
        molecule = tensor_path.parents[2].name
        geometry_tag = tensor_path.parents[1].name
        basis = tensor_path.parent.name
    except IndexError as error:
        raise ValueError(
            "Expected tensor path layout "
            "polyatomic_library/<molecule>/s-<scale>/<basis>/<file>."
        ) from error

    scale_match = SCALE_PATTERN.fullmatch(geometry_tag)
    if scale_match is None:
        raise ValueError(
            f"Expected geometry directory 's-<scale>', found {geometry_tag!r}."
        )

    return CaseMetadata(
        case_id=tensor_path.name.removesuffix(".tensors.npz"),
        tensor_path=tensor_path,
        molecule=molecule,
        atom_count=chemical_formula_atom_count(molecule),
        geometry_tag=geometry_tag,
        geometry_scale=float(scale_match.group("scale")),
        basis=basis,
        source_mapping=(match.group("source_mapping") or "JW").upper(),
        active_occupied=active_occupied,
        active_vacant=active_vacant,
        n_qubits=n_qubits,
    )


def tensor_sort_key(path: Path) -> tuple[int, str]:
    """Sort small active spaces first so ``--limit 1`` is a cheap test."""
    match = CASE_PATTERN.match(path.name)
    if match is None:
        return (10**9, str(path))

    n_qubits = int(match.group("active_occupied")) + int(
        match.group("active_vacant")
    )
    return (n_qubits, str(path))


def discover_tensor_paths(
    library: Path,
    case_pattern: str | None,
    molecules: Sequence[str] | None,
    bases: Sequence[str] | None,
    scales: Sequence[float] | None,
    max_qubits: int | None,
    limit: int | None,
) -> list[Path]:
    paths = sorted(library.rglob("*.tensors.npz"), key=tensor_sort_key)

    molecule_filter = {item.lower() for item in molecules or []}
    basis_filter = {item.lower() for item in bases or []}
    scale_filter = list(scales or [])

    selected: list[Path] = []
    for path in paths:
        match = CASE_PATTERN.match(path.name)
        if match is None:
            continue

        n_qubits = int(match.group("active_occupied")) + int(
            match.group("active_vacant")
        )
        try:
            molecule = path.parents[2].name
            geometry_tag = path.parents[1].name
            basis = path.parent.name
        except IndexError:
            continue

        scale_match = SCALE_PATTERN.fullmatch(geometry_tag)
        if scale_match is None:
            continue
        scale = float(scale_match.group("scale"))

        if case_pattern and case_pattern not in str(path):
            continue
        if molecule_filter and molecule.lower() not in molecule_filter:
            continue
        if basis_filter and basis.lower() not in basis_filter:
            continue
        if scale_filter and not any(
            math.isclose(scale, requested, rel_tol=0.0, abs_tol=1.0e-9)
            for requested in scale_filter
        ):
            continue
        if max_qubits is not None and n_qubits > max_qubits:
            continue

        selected.append(path)

    if limit is not None:
        selected = selected[:limit]

    return selected


def load_interaction_operator(
    tensor_path: Path,
) -> tuple[InteractionOperator, int]:
    with np.load(tensor_path, allow_pickle=False) as data:
        required = {"constant", "one_body", "two_body"}
        missing = required.difference(data.files)
        if missing:
            raise ValueError(
                f"Missing arrays {sorted(missing)} in {tensor_path}."
            )

        constant = complex(np.asarray(data["constant"]).item())
        one_body = np.asarray(data["one_body"], dtype=complex)
        two_body = np.asarray(data["two_body"], dtype=complex)

    if one_body.ndim != 2 or one_body.shape[0] != one_body.shape[1]:
        raise ValueError(f"Invalid one-body shape {one_body.shape}.")

    n_qubits = int(one_body.shape[0])
    expected_two_body_shape = (n_qubits,) * 4
    if two_body.shape != expected_two_body_shape:
        raise ValueError(
            f"Invalid two-body shape {two_body.shape}; expected "
            f"{expected_two_body_shape}."
        )

    interaction = InteractionOperator(constant, one_body, two_body)
    return interaction, n_qubits


# -----------------------------------------------------------------------------
# Fermionic terms and graphs
# -----------------------------------------------------------------------------


def clean_fermion_operator(
    operator: FermionOperator,
    tolerance: float,
) -> FermionOperator:
    cleaned = normal_ordered(operator)
    cleaned.compress(abs_tol=tolerance)
    return cleaned


def adjoint_key_and_phase(
    term_key: tuple[tuple[int, int], ...],
    tolerance: float,
) -> tuple[tuple[tuple[int, int], ...], complex]:
    unit_term = FermionOperator(term_key, 1.0)
    adjoint = clean_fermion_operator(
        hermitian_conjugated(unit_term),
        tolerance,
    )

    if len(adjoint.terms) != 1:
        raise ValueError(
            "Expected the adjoint of one normal-ordered monomial to "
            f"contain one monomial, but found {len(adjoint.terms)} "
            f"for {term_key}."
        )

    return next(iter(adjoint.terms.items()))


def build_hermitian_fermion_terms(
    fermion_hamiltonian: FermionOperator,
    tolerance: float,
) -> list[HermitianFermionTerm]:
    """
    Pair each normal-ordered monomial with its Hermitian adjoint.

    The returned list contains one object for each complete Hermitian term,
    not one row per monomial. The order is determined by first appearance in
    the OpenFermion ``FermionOperator.terms`` insertion order.
    """
    raw_items = list(fermion_hamiltonian.terms.items())
    raw_index = {key: index for index, (key, _) in enumerate(raw_items)}

    visited: set[tuple[tuple[int, int], ...]] = set()
    hermitian_terms: list[HermitianFermionTerm] = []
    reconstruction = FermionOperator()

    for term_key, coefficient in raw_items:
        if term_key in visited or abs(coefficient) <= tolerance:
            continue

        partner_key, adjoint_phase = adjoint_key_and_phase(
            term_key,
            tolerance,
        )

        if partner_key not in fermion_hamiltonian.terms:
            raise ValueError(
                "Missing Hermitian-adjoint partner for fermionic key "
                f"{term_key}; expected {partner_key}."
            )

        partner_coefficient = fermion_hamiltonian.terms[partner_key]
        expected_partner = complex(coefficient).conjugate() * adjoint_phase

        if abs(partner_coefficient - expected_partner) > tolerance:
            raise ValueError(
                "Hermitian partner coefficient mismatch:\n"
                f"  key: {term_key}\n"
                f"  partner: {partner_key}\n"
                f"  actual partner coefficient: {partner_coefficient}\n"
                f"  expected partner coefficient: {expected_partner}"
            )

        term_operator = FermionOperator(term_key, coefficient)
        component_keys = [term_key]

        if partner_key != term_key:
            term_operator += FermionOperator(
                partner_key,
                partner_coefficient,
            )
            component_keys.append(partner_key)

        hermiticity_error = clean_fermion_operator(
            term_operator - hermitian_conjugated(term_operator),
            tolerance,
        )
        if hermiticity_error.terms:
            raise ValueError(
                "Constructed fermionic term is not Hermitian: "
                f"{term_operator}"
            )

        first_raw_index = min(raw_index[key] for key in component_keys)
        hermitian_terms.append(
            HermitianFermionTerm(
                index=len(hermitian_terms),
                first_raw_index=first_raw_index,
                operator=term_operator,
                component_keys=tuple(component_keys),
            )
        )

        reconstruction += term_operator
        visited.update(component_keys)

    reconstruction_error = clean_fermion_operator(
        fermion_hamiltonian - reconstruction,
        tolerance,
    )
    if reconstruction_error.terms:
        raise ValueError(
            "Hermitian-term grouping did not reconstruct the full "
            f"fermionic Hamiltonian:\n{reconstruction_error}"
        )

    # The identity contributes only a global phase and has no graph edges.
    return [
        term
        for term in hermitian_terms
        if term.component_keys != ((),)
    ]


def fermionic_terms_noncommute(
    left: FermionOperator,
    right: FermionOperator,
    tolerance: float,
) -> bool:
    result = clean_fermion_operator(
        commutator(left, right),
        tolerance,
    )
    return bool(result.terms)


def build_fermionic_noncommutation_graph(
    terms: Sequence[HermitianFermionTerm],
    tolerance: float,
) -> tuple[nx.Graph, float]:
    start = time.perf_counter()

    graph = nx.Graph()
    graph.add_nodes_from(range(len(terms)))

    for left_index in range(len(terms)):
        for right_index in range(left_index + 1, len(terms)):
            if fermionic_terms_noncommute(
                terms[left_index].operator,
                terms[right_index].operator,
                tolerance,
            ):
                graph.add_edge(left_index, right_index)

    return graph, time.perf_counter() - start


# -----------------------------------------------------------------------------
# Pauli keys, graphs, and deterministic coloring orders
# -----------------------------------------------------------------------------


def pauli_keys_commute(
    left: tuple[tuple[int, str], ...],
    right: tuple[tuple[int, str], ...],
) -> bool:
    left_map = dict(left)
    right_map = dict(right)

    anti_commuting_positions = sum(
        left_map[qubit] != right_map[qubit]
        for qubit in left_map.keys() & right_map.keys()
    )
    return anti_commuting_positions % 2 == 0


def build_pauli_noncommutation_graph(
    pauli_keys: Sequence[tuple[tuple[int, str], ...]],
) -> tuple[nx.Graph, float]:
    start = time.perf_counter()

    graph = nx.Graph()
    graph.add_nodes_from(range(len(pauli_keys)))

    for left_index in range(len(pauli_keys)):
        for right_index in range(left_index + 1, len(pauli_keys)):
            if not pauli_keys_commute(
                pauli_keys[left_index],
                pauli_keys[right_index],
            ):
                graph.add_edge(left_index, right_index)

    return graph, time.perf_counter() - start


def deterministic_coloring_order(
    graph: nx.Graph,
) -> tuple[list[int], int, float]:
    """
    Use ``largest_first`` coloring, then order color groups by color label.

    Within each color group, preserve the original node index. This same rule
    is used for both JW and fermionic coloring so that the graph level is the
    only structural difference between the two methods.
    """
    start = time.perf_counter()
    coloring = nx.greedy_color(graph, strategy="largest_first")
    coloring_time = time.perf_counter() - start

    if not coloring:
        return [], 0, coloring_time

    colors = sorted(set(coloring.values()))
    ordered_nodes = [
        node
        for color in colors
        for node in graph.nodes
        if coloring[node] == color
    ]
    return ordered_nodes, len(colors), coloring_time


def validate_pauli_order(
    name: str,
    order: Sequence[tuple[tuple[int, str], ...]],
    raw_pauli_keys: Sequence[tuple[tuple[int, str], ...]],
) -> None:
    if len(order) != len(set(order)):
        raise ValueError(f"Ordering {name!r} contains duplicate Pauli keys.")

    if set(order) != set(raw_pauli_keys):
        missing = set(raw_pauli_keys).difference(order)
        extra = set(order).difference(raw_pauli_keys)
        raise ValueError(
            f"Ordering {name!r} is not a permutation of the final JW terms. "
            f"Missing={len(missing)}, extra={len(extra)}."
        )


def build_orderings(
    fermion_hamiltonian: FermionOperator,
    full_jw_hamiltonian: Any,
    requested_orderings: Sequence[str],
    tolerance: float,
) -> tuple[
    dict[str, OrderingResult],
    list[HermitianFermionTerm],
    dict[tuple[tuple[int, str], ...], complex],
]:
    final_coefficients = {
        key: coefficient
        for key, coefficient in full_jw_hamiltonian.terms.items()
        if key != () and abs(coefficient) > tolerance
    }
    raw_pauli_keys = list(final_coefficients)

    results: dict[str, OrderingResult] = {}
    results["jw_raw"] = OrderingResult(
        name="jw_raw",
        pauli_keys=list(raw_pauli_keys),
        build_time_seconds=0.0,
    )

    hermitian_terms: list[HermitianFermionTerm] = []

    if "jw_coloring" in requested_orderings:
        build_start = time.perf_counter()
        jw_graph, graph_time = build_pauli_noncommutation_graph(
            raw_pauli_keys
        )
        node_order, number_of_colors, coloring_time = (
            deterministic_coloring_order(jw_graph)
        )
        pauli_order = [raw_pauli_keys[node] for node in node_order]
        results["jw_coloring"] = OrderingResult(
            name="jw_coloring",
            pauli_keys=pauli_order,
            graph_vertices=jw_graph.number_of_nodes(),
            graph_edges=jw_graph.number_of_edges(),
            number_of_colors=number_of_colors,
            graph_time_seconds=graph_time,
            coloring_time_seconds=coloring_time,
            build_time_seconds=time.perf_counter() - build_start,
        )

    if "fermionic_coloring" in requested_orderings:
        build_start = time.perf_counter()
        hermitian_terms = build_hermitian_fermion_terms(
            fermion_hamiltonian,
            tolerance,
        )
        fermionic_graph, graph_time = (
            build_fermionic_noncommutation_graph(
                hermitian_terms,
                tolerance,
            )
        )
        node_order, number_of_colors, coloring_time = (
            deterministic_coloring_order(fermionic_graph)
        )

        ordered_fermionic_terms = [
            hermitian_terms[node] for node in node_order
        ]

        induced_pauli_order: list[tuple[tuple[int, str], ...]] = []
        seen: set[tuple[tuple[int, str], ...]] = set()

        for fermionic_term in ordered_fermionic_terms:
            mapped_term = jordan_wigner(fermionic_term.operator)
            mapped_term.compress(abs_tol=tolerance)

            for pauli_key, coefficient in mapped_term.terms.items():
                if pauli_key == () or abs(coefficient) <= tolerance:
                    continue
                if pauli_key not in final_coefficients:
                    # This Pauli key cancels in the fully combined Hamiltonian.
                    continue
                if pauli_key not in seen:
                    seen.add(pauli_key)
                    induced_pauli_order.append(pauli_key)

        # This should normally add nothing, but it makes the behavior explicit
        # and guarantees a full permutation if a term was removed by tolerance.
        for pauli_key in raw_pauli_keys:
            if pauli_key not in seen:
                seen.add(pauli_key)
                induced_pauli_order.append(pauli_key)

        results["fermionic_coloring"] = OrderingResult(
            name="fermionic_coloring",
            pauli_keys=induced_pauli_order,
            graph_vertices=fermionic_graph.number_of_nodes(),
            graph_edges=fermionic_graph.number_of_edges(),
            number_of_colors=number_of_colors,
            graph_time_seconds=graph_time,
            coloring_time_seconds=coloring_time,
            build_time_seconds=time.perf_counter() - build_start,
        )

    for ordering_name in requested_orderings:
        if ordering_name not in results:
            raise ValueError(f"Failed to build ordering {ordering_name!r}.")
        validate_pauli_order(
            ordering_name,
            results[ordering_name].pauli_keys,
            raw_pauli_keys,
        )

    return results, hermitian_terms, final_coefficients


# -----------------------------------------------------------------------------
# Dense matrices and Trotter formulas
# -----------------------------------------------------------------------------


def pauli_matrix_from_key(
    pauli_key: tuple[tuple[int, str], ...],
    n_qubits: int,
) -> np.ndarray:
    single_qubit = {
        "I": np.eye(2, dtype=complex),
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.array([[1, 0], [0, -1]], dtype=complex),
    }

    explicit = dict(pauli_key)
    result = np.array([[1.0 + 0.0j]])

    # Qubit 0 is the most-significant tensor factor, matching QHAT's
    # computational-basis initial-state convention.
    for qubit in range(n_qubits):
        result = np.kron(
            result,
            single_qubit[explicit.get(qubit, "I")],
        )

    return result


def real_coefficient(
    coefficient: complex,
    tolerance: float,
) -> float:
    coefficient = complex(coefficient)
    if abs(coefficient.imag) > tolerance:
        raise ValueError(
            "Expected a real coefficient for a Hermitian Pauli term, "
            f"but found {coefficient}."
        )
    return float(coefficient.real)


def pauli_exponential(
    coefficient: float,
    pauli_matrix: np.ndarray,
    duration: float,
    identity: np.ndarray,
) -> np.ndarray:
    angle = coefficient * duration
    return (
        math.cos(angle) * identity
        - 1j * math.sin(angle) * pauli_matrix
    )


def first_order_slice(
    ordered_terms: Sequence[tuple[float, np.ndarray]],
    dt: float,
    identity: np.ndarray,
) -> np.ndarray:
    unitary = identity.copy()
    for coefficient, pauli_matrix in ordered_terms:
        # The list order is the chronological application order.
        unitary = pauli_exponential(
            coefficient,
            pauli_matrix,
            dt,
            identity,
        ) @ unitary
    return unitary


def second_order_slice(
    ordered_terms: Sequence[tuple[float, np.ndarray]],
    dt: float,
    identity: np.ndarray,
) -> np.ndarray:
    """Symmetric second-order formula using 2m-1 exponentials per step."""
    if not ordered_terms:
        return identity.copy()

    unitary = identity.copy()

    for coefficient, pauli_matrix in ordered_terms[:-1]:
        unitary = pauli_exponential(
            coefficient,
            pauli_matrix,
            dt / 2.0,
            identity,
        ) @ unitary

    center_coefficient, center_matrix = ordered_terms[-1]
    unitary = pauli_exponential(
        center_coefficient,
        center_matrix,
        dt,
        identity,
    ) @ unitary

    for coefficient, pauli_matrix in reversed(ordered_terms[:-1]):
        unitary = pauli_exponential(
            coefficient,
            pauli_matrix,
            dt / 2.0,
            identity,
        ) @ unitary

    return unitary


def fourth_order_slice(
    ordered_terms: Sequence[tuple[float, np.ndarray]],
    dt: float,
    identity: np.ndarray,
) -> np.ndarray:
    """
    Symmetric fourth-order Suzuki formula.

    This uses the five-stage recursive composition

        S4(dt) = S2(p dt) S2(p dt) S2((1 - 4p) dt)
                 S2(p dt) S2(p dt),

    where p = 1 / (4 - 4**(1/3)). The middle stage has a negative
    duration, which is expected for this fourth-order construction.

    The returned matrix represents one complete fourth-order Trotter
    slice. Adjacent exponentials are not algebraically fused when
    reporting the nominal exponential count.
    """
    if not ordered_terms:
        return identity.copy()

    p = 1.0 / (4.0 - 4.0 ** (1.0 / 3.0))
    stage_scales = (p, p, 1.0 - 4.0 * p, p, p)

    unitary = identity.copy()
    for scale in stage_scales:
        stage = second_order_slice(
            ordered_terms,
            scale * dt,
            identity,
        )
        unitary = stage @ unitary

    return unitary


def build_hartree_fock_state(
    n_qubits: int,
    active_occupied: int,
) -> np.ndarray:
    if not 0 <= active_occupied <= n_qubits:
        raise ValueError(
            f"Invalid occupied count {active_occupied} for {n_qubits} qubits."
        )

    basis_index = sum(
        1 << (n_qubits - 1 - qubit)
        for qubit in range(active_occupied)
    )
    state = np.zeros(2**n_qubits, dtype=complex)
    state[basis_index] = 1.0
    return state


def spectral_norm(matrix: np.ndarray) -> float:
    return float(np.linalg.norm(matrix, ord=2))


def state_infidelity(
    exact_state: np.ndarray,
    approximate_state: np.ndarray,
) -> float:
    fidelity = abs(np.vdot(exact_state, approximate_state)) ** 2
    fidelity = min(1.0, max(0.0, float(fidelity.real)))
    return 1.0 - fidelity


# -----------------------------------------------------------------------------
# CSV and resume helpers
# -----------------------------------------------------------------------------


FIELDNAMES = [
    "status",
    "error_message",
    "case_id",
    "tensor_path",
    "molecule",
    "atom_count",
    "geometry_tag",
    "geometry_scale",
    "basis",
    "source_mapping",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "identity_coefficient_real",
    "identity_coefficient_imag",
    "number_of_fermionic_terms",
    "number_of_pauli_terms",
    "pauli_to_fermionic_term_ratio",
    "ordering",
    "formula_order",
    "trotter_steps",
    "trotter_dt",
    "evolution_time",
    "nominal_exponential_count",
    "ordering_build_time_seconds",
    "graph_vertices",
    "graph_edges",
    "graph_density",
    "number_of_colors",
    "mean_vertices_per_color",
    "graph_time_seconds",
    "coloring_time_seconds",
    "exact_evolution_time_seconds",
    "trotter_runtime_seconds",
    "operator_norm_error",
    "state_infidelity",
    "state_vector_2norm_error",
    "operator_error_ratio_to_jw_raw",
    "state_infidelity_ratio_to_jw_raw",
    "operator_error_ratio_to_jw_coloring",
    "state_infidelity_ratio_to_jw_coloring",
    "coefficient_tolerance",
]


ErrorCache = dict[
    tuple[str, str, int, int, float],
    tuple[float, float, float],
]


def load_resume_data(
    output_path: Path,
) -> tuple[set[tuple[str, str, int, int, float]], ErrorCache]:
    completed: set[tuple[str, str, int, int, float]] = set()
    error_cache: ErrorCache = {}

    if not output_path.exists():
        return completed, error_cache

    with output_path.open("r", newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        for row in reader:
            if row.get("status") != "success":
                continue
            try:
                key = (
                    row["case_id"],
                    row["ordering"],
                    int(row["formula_order"]),
                    int(row["trotter_steps"]),
                    float(row["evolution_time"]),
                )
                operator_text = row.get("operator_norm_error", "").strip()
                operator_error = float(operator_text) if operator_text else math.nan
                infidelity = float(row["state_infidelity"])
                state_vector_error = float(row["state_vector_2norm_error"])
            except (KeyError, TypeError, ValueError):
                continue

            completed.add(key)
            error_cache[key] = (
                operator_error,
                infidelity,
                state_vector_error,
            )

    return completed, error_cache


def blank_row() -> dict[str, Any]:
    return {field: "" for field in FIELDNAMES}


def add_case_metadata(
    row: dict[str, Any],
    metadata: CaseMetadata,
) -> None:
    row.update(
        {
            "case_id": metadata.case_id,
            "tensor_path": str(metadata.tensor_path),
            "molecule": metadata.molecule,
            "atom_count": metadata.atom_count,
            "geometry_tag": metadata.geometry_tag,
            "geometry_scale": metadata.geometry_scale,
            "basis": metadata.basis,
            "source_mapping": metadata.source_mapping,
            "active_occupied": metadata.active_occupied,
            "active_vacant": metadata.active_vacant,
            "n_qubits": metadata.n_qubits,
        }
    )


# -----------------------------------------------------------------------------
# One-case benchmark
# -----------------------------------------------------------------------------


def benchmark_case(
    tensor_path: Path,
    requested_orderings: Sequence[str],
    formula_orders: Sequence[int],
    steps_list: Sequence[int],
    evolution_time: float,
    tolerance: float,
    compute_operator_norm: bool,
    completed: set[tuple[str, str, int, int, float]],
    error_cache: ErrorCache,
    writer: csv.DictWriter,
    output_file: Any,
) -> None:
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)

    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        tolerance,
    )

    full_jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    full_jw_hamiltonian.compress(abs_tol=tolerance)

    identity_coefficient = complex(
        full_jw_hamiltonian.terms.get((), 0.0)
    )

    orderings, hermitian_terms, final_coefficients = build_orderings(
        fermion_hamiltonian=fermion_hamiltonian,
        full_jw_hamiltonian=full_jw_hamiltonian,
        requested_orderings=requested_orderings,
        tolerance=tolerance,
    )

    raw_pauli_keys = orderings["jw_raw"].pauli_keys
    matrix_cache = {
        key: pauli_matrix_from_key(key, n_qubits)
        for key in raw_pauli_keys
    }

    dimension = 2**n_qubits
    identity = np.eye(dimension, dtype=complex)

    exact_hamiltonian = np.zeros(
        (dimension, dimension),
        dtype=complex,
    )
    for pauli_key, coefficient in final_coefficients.items():
        exact_hamiltonian += real_coefficient(
            coefficient,
            tolerance,
        ) * matrix_cache[pauli_key]

    hermiticity_error = spectral_norm(
        exact_hamiltonian - exact_hamiltonian.conjugate().T
    )
    if hermiticity_error > 100.0 * tolerance:
        raise ValueError(
            "JW Hamiltonian matrix is not Hermitian; error="
            f"{hermiticity_error}."
        )

    exact_start = time.perf_counter()
    exact_unitary = expm(-1j * evolution_time * exact_hamiltonian)
    exact_evolution_time = time.perf_counter() - exact_start

    initial_state = build_hartree_fock_state(
        n_qubits,
        metadata.active_occupied,
    )
    exact_state = exact_unitary @ initial_state

    # Process both references before fermionic coloring so all ratios are
    # available immediately. Any additional future ordering follows afterward.
    preferred_order = ["jw_raw", "jw_coloring", "fermionic_coloring"]
    processing_order = [
        name for name in preferred_order if name in requested_orderings
    ] + [
        name for name in requested_orderings if name not in preferred_order
    ]

    number_of_fermionic_terms = (
        len(hermitian_terms)
        if hermitian_terms
        else len(build_hermitian_fermion_terms(
            fermion_hamiltonian,
            tolerance,
        ))
    )
    number_of_pauli_terms = len(raw_pauli_keys)
    pauli_to_fermionic_term_ratio = (
        number_of_pauli_terms / number_of_fermionic_terms
        if number_of_fermionic_terms
        else math.nan
    )

    for formula_order in formula_orders:
        for trotter_steps in steps_list:
            dt = evolution_time / trotter_steps

            for ordering_name in processing_order:
                result_key = (
                    metadata.case_id,
                    ordering_name,
                    formula_order,
                    trotter_steps,
                    evolution_time,
                )
                if result_key in completed:
                    print(
                        "  SKIP completed: "
                        f"{ordering_name}, order={formula_order}, "
                        f"steps={trotter_steps}"
                    )
                    continue

                ordering = orderings[ordering_name]
                ordered_terms = [
                    (
                        real_coefficient(
                            final_coefficients[pauli_key],
                            tolerance,
                        ),
                        matrix_cache[pauli_key],
                    )
                    for pauli_key in ordering.pauli_keys
                ]

                trotter_start = time.perf_counter()

                if formula_order == 1:
                    one_slice = first_order_slice(
                        ordered_terms,
                        dt,
                        identity,
                    )
                    nominal_exponential_count = (
                        trotter_steps * number_of_pauli_terms
                    )
                elif formula_order == 2:
                    one_slice = second_order_slice(
                        ordered_terms,
                        dt,
                        identity,
                    )
                    nominal_exponential_count = trotter_steps * (
                        2 * number_of_pauli_terms - 1
                    )
                elif formula_order == 4:
                    one_slice = fourth_order_slice(
                        ordered_terms,
                        dt,
                        identity,
                    )
                    nominal_exponential_count = trotter_steps * (
                        5 * (2 * number_of_pauli_terms - 1)
                    )
                else:
                    raise ValueError(
                        f"Unsupported formula order {formula_order}."
                    )

                trotter_unitary = np.linalg.matrix_power(
                    one_slice,
                    trotter_steps,
                )
                trotter_runtime = time.perf_counter() - trotter_start

                approximate_state = trotter_unitary @ initial_state
                infidelity = state_infidelity(
                    exact_state,
                    approximate_state,
                )
                state_vector_error = float(
                    np.linalg.norm(
                        exact_state - approximate_state,
                        ord=2,
                    )
                )

                if compute_operator_norm:
                    operator_error: float | str = spectral_norm(
                        trotter_unitary - exact_unitary
                    )
                else:
                    operator_error = ""

                current_key = (
                    metadata.case_id,
                    ordering_name,
                    formula_order,
                    trotter_steps,
                    evolution_time,
                )
                error_cache[current_key] = (
                    float(operator_error) if operator_error != "" else math.nan,
                    infidelity,
                    state_vector_error,
                )

                def ratios_to(reference_name: str) -> tuple[float | str, float | str]:
                    reference_key = (
                        metadata.case_id,
                        reference_name,
                        formula_order,
                        trotter_steps,
                        evolution_time,
                    )
                    reference = error_cache.get(reference_key)
                    if reference is None:
                        return "", ""

                    reference_operator, reference_infidelity, _ = reference

                    if ordering_name == reference_name:
                        operator_ratio: float | str = 1.0
                        infidelity_ratio: float | str = 1.0
                    else:
                        if (
                            operator_error == ""
                            or not math.isfinite(reference_operator)
                            or reference_operator <= tolerance
                        ):
                            operator_ratio = ""
                        else:
                            operator_ratio = float(operator_error) / reference_operator

                        if reference_infidelity <= tolerance:
                            infidelity_ratio = ""
                        else:
                            infidelity_ratio = infidelity / reference_infidelity

                    return operator_ratio, infidelity_ratio

                operator_ratio_raw, infidelity_ratio_raw = ratios_to("jw_raw")
                operator_ratio_jw_color, infidelity_ratio_jw_color = ratios_to(
                    "jw_coloring"
                )

                row = blank_row()
                add_case_metadata(row, metadata)
                row.update(
                    {
                        "status": "success",
                        "identity_coefficient_real": (
                            identity_coefficient.real
                        ),
                        "identity_coefficient_imag": (
                            identity_coefficient.imag
                        ),
                        "number_of_fermionic_terms": (
                            number_of_fermionic_terms
                        ),
                        "number_of_pauli_terms": number_of_pauli_terms,
                        "pauli_to_fermionic_term_ratio": (
                            pauli_to_fermionic_term_ratio
                        ),
                        "ordering": ordering_name,
                        "formula_order": formula_order,
                        "trotter_steps": trotter_steps,
                        "trotter_dt": dt,
                        "evolution_time": evolution_time,
                        "nominal_exponential_count": (
                            nominal_exponential_count
                        ),
                        "ordering_build_time_seconds": (
                            ordering.build_time_seconds
                        ),
                        "graph_vertices": ordering.graph_vertices,
                        "graph_edges": ordering.graph_edges,
                        "graph_density": (
                            nx.density(nx.empty_graph(0))
                            if ordering.graph_vertices == 0
                            else (
                                2.0 * ordering.graph_edges
                                / (ordering.graph_vertices * (ordering.graph_vertices - 1))
                                if ordering.graph_vertices is not None
                                and ordering.graph_edges is not None
                                and ordering.graph_vertices > 1
                                else ""
                            )
                        ),
                        "number_of_colors": ordering.number_of_colors,
                        "mean_vertices_per_color": (
                            ordering.graph_vertices / ordering.number_of_colors
                            if ordering.graph_vertices is not None
                            and ordering.number_of_colors
                            else ""
                        ),
                        "graph_time_seconds": (
                            ordering.graph_time_seconds
                        ),
                        "coloring_time_seconds": (
                            ordering.coloring_time_seconds
                        ),
                        "exact_evolution_time_seconds": (
                            exact_evolution_time
                        ),
                        "trotter_runtime_seconds": trotter_runtime,
                        "operator_norm_error": operator_error,
                        "state_infidelity": infidelity,
                        "state_vector_2norm_error": (
                            state_vector_error
                        ),
                        "operator_error_ratio_to_jw_raw": (
                            operator_ratio_raw
                        ),
                        "state_infidelity_ratio_to_jw_raw": (
                            infidelity_ratio_raw
                        ),
                        "operator_error_ratio_to_jw_coloring": (
                            operator_ratio_jw_color
                        ),
                        "state_infidelity_ratio_to_jw_coloring": (
                            infidelity_ratio_jw_color
                        ),
                        "coefficient_tolerance": tolerance,
                    }
                )

                writer.writerow(row)
                output_file.flush()
                completed.add(result_key)

                print(
                    "  DONE "
                    f"{ordering_name:20s} "
                    f"order={formula_order} "
                    f"steps={trotter_steps:4d} "
                    f"op_error={operator_error} "
                    f"infidelity={infidelity:.6e}"
                )


# -----------------------------------------------------------------------------
# Command line
# -----------------------------------------------------------------------------


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare raw JW, JW coloring, and fermionic-color-induced JW "
            "Trotter orderings for polyatomic geometry-sweep tensor files."
        )
    )

    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/polyatomic_library"),
        help=(
            "Root directory containing molecule/s-<scale>/<basis>/"
            "*.tensors.npz files."
        ),
    )
    parser.add_argument(
        "--molecules",
        nargs="+",
        default=None,
        help="Optional molecule directories to include, e.g. BeH2 H2O NH3 CH4.",
    )
    parser.add_argument(
        "--bases",
        nargs="+",
        default=None,
        help="Optional basis filters, e.g. sto-6g hgbs-5.",
    )
    parser.add_argument(
        "--scales",
        type=float,
        nargs="+",
        default=None,
        help="Optional geometry-scale filters, e.g. 0.8 1.0 1.5.",
    )
    parser.add_argument(
        "--max-qubits",
        type=int,
        default=8,
        help=(
            "Skip active spaces larger than this value (default: 8). "
            "Dense exact evolution grows exponentially."
        ),
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=[10, 100, 1000],
        help="Fixed Trotter step counts.",
    )
    parser.add_argument(
        "--formula-orders",
        type=int,
        nargs="+",
        choices=[1, 2, 4],
        default=[1, 2, 4],
        help="Product-formula orders to benchmark: 1, 2, and/or 4.",
    )
    parser.add_argument(
        "--orderings",
        nargs="+",
        choices=SUPPORTED_ORDERINGS,
        default=list(SUPPORTED_ORDERINGS),
        help="Ordering methods to benchmark.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        dest="evolution_time",
        help="Total evolution time.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N selected tensor files after sorting.",
    )
    parser.add_argument(
        "--case-pattern",
        type=str,
        default=None,
        help="Keep only tensor paths containing this substring.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/polyatomic_trotter_results.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Coefficient and commutator compression tolerance.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Append to an existing CSV and skip completed rows.",
    )
    parser.add_argument(
        "--skip-operator-norm",
        action="store_true",
        help=(
            "Skip the spectral operator norm. State errors are much cheaper "
            "for a first full-pipeline test."
        ),
    )

    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.library.exists():
        raise FileNotFoundError(args.library)
    if any(step <= 0 for step in args.steps):
        raise ValueError("All Trotter step counts must be positive.")
    if args.evolution_time <= 0.0:
        raise ValueError("Evolution time must be positive.")
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be positive.")
    if args.max_qubits <= 0:
        raise ValueError("--max-qubits must be positive.")
    if args.scales is not None and any(scale <= 0.0 for scale in args.scales):
        raise ValueError("All --scales values must be positive.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)

    requested_orderings = list(dict.fromkeys(args.orderings))
    if "jw_raw" not in requested_orderings:
        requested_orderings.insert(0, "jw_raw")
        print(
            "Added jw_raw automatically because it is required for baseline "
            "error ratios."
        )
    if (
        "fermionic_coloring" in requested_orderings
        and "jw_coloring" not in requested_orderings
    ):
        requested_orderings.insert(1, "jw_coloring")
        print(
            "Added jw_coloring automatically because it is the direct "
            "reference for fermionic coloring."
        )

    tensor_paths = discover_tensor_paths(
        library=args.library,
        case_pattern=args.case_pattern,
        molecules=args.molecules,
        bases=args.bases,
        scales=args.scales,
        max_qubits=args.max_qubits,
        limit=args.limit,
    )

    if not tensor_paths:
        raise FileNotFoundError(
            f"No matching *.tensors.npz files found under {args.library}. "
            "Confirm that the tensor-building stage has already been run."
        )

    print(f"Tensor cases selected: {len(tensor_paths)}")
    print(f"Molecule filter: {args.molecules or 'all'}")
    print(f"Basis filter: {args.bases or 'all'}")
    print(f"Scale filter: {args.scales or 'all'}")
    print(f"Maximum qubits: {args.max_qubits}")
    print(f"Orderings: {requested_orderings}")
    print(f"Formula orders: {args.formula_orders}")
    print(f"Trotter steps: {args.steps}")
    print(f"Evolution time: {args.evolution_time}")
    print(f"Output: {args.output}")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.resume:
        completed, error_cache = load_resume_data(args.output)
        file_mode = "a"
        write_header = not args.output.exists() or args.output.stat().st_size == 0
        print(f"Resume rows already completed: {len(completed)}")
    else:
        completed = set()
        error_cache = {}
        file_mode = "w"
        write_header = True

    successful_cases = 0
    failed_cases = 0

    with args.output.open(
        file_mode,
        newline="",
        encoding="utf-8",
    ) as output_file:
        writer = csv.DictWriter(output_file, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
            output_file.flush()

        for case_index, tensor_path in enumerate(tensor_paths, start=1):
            print()
            print("=" * 80)
            print(f"[{case_index}/{len(tensor_paths)}] {tensor_path}")

            try:
                benchmark_case(
                    tensor_path=tensor_path,
                    requested_orderings=requested_orderings,
                    formula_orders=args.formula_orders,
                    steps_list=args.steps,
                    evolution_time=args.evolution_time,
                    tolerance=args.tolerance,
                    compute_operator_norm=not args.skip_operator_norm,
                    completed=completed,
                    error_cache=error_cache,
                    writer=writer,
                    output_file=output_file,
                )
            except Exception as error:
                failed_cases += 1
                print(
                    f"FAILED {tensor_path}: "
                    f"{type(error).__name__}: {error}"
                )

                failure_row = blank_row()
                failure_row.update(
                    {
                        "status": "failed",
                        "error_message": f"{type(error).__name__}: {error}",
                        "case_id": tensor_path.name.removesuffix(
                            ".tensors.npz"
                        ),
                        "tensor_path": str(tensor_path),
                        "coefficient_tolerance": args.tolerance,
                    }
                )
                writer.writerow(failure_row)
                output_file.flush()
            else:
                successful_cases += 1

    print()
    print("=" * 80)
    print(f"Cases selected:   {len(tensor_paths)}")
    print(f"Cases successful: {successful_cases}")
    print(f"Cases failed:     {failed_cases}")
    print(f"Results:          {args.output}")
    print("=" * 80)


if __name__ == "__main__":
    main()
