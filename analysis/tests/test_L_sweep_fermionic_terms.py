"""Focused tests for L-sweep Pauli and fermionic-term schedules."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from openfermion import FermionOperator, hermitian_conjugated, jordan_wigner
from scipy.linalg import expm

matplotlib.use("Agg")

from analysis.benchmark_L_sweep_trotter import (  # noqa: E402
    FIELDNAMES,
    LEGACY_FIELDNAMES,
    SCHEMA_VERSION,
    SUPPORTED_ORDERINGS,
    HamiltonianFactor,
    build_hermitian_fermion_terms,
    build_orderings,
    clean_fermion_operator,
    deterministic_coloring_order,
    first_order_slice,
    fourth_order_slice,
    load_resume_data,
    nominal_exponential_count,
    normalize_requested_orderings,
    parse_arguments,
    prepare_resume_output,
    second_order_slice,
)
from analysis.visualize_L_sweep_trotter import (  # noqa: E402
    build_case_ranking,
    load_results,
    plot_case_grid,
    plot_complexity_vs_error,
    plot_win_rates,
    summary_target,
)


TOLERANCE = 1.0e-12


def small_fermion_hamiltonian() -> FermionOperator:
    """Two-mode Hermitian fixture with identity, diagonal, and hopping terms."""
    hamiltonian = FermionOperator((), 0.7)
    hamiltonian += FermionOperator("0^ 0", 0.5)
    hamiltonian += FermionOperator("1^ 1", -0.2)
    hamiltonian += FermionOperator("0^ 1", 0.3)
    hamiltonian += FermionOperator("1^ 0", 0.3)
    hamiltonian += FermionOperator("0^ 1^ 1 0", 0.4)
    return clean_fermion_operator(hamiltonian, TOLERANCE)


def all_schedules():
    fermion_hamiltonian = small_fermion_hamiltonian()
    jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    jw_hamiltonian.compress(abs_tol=TOLERANCE)
    schedules, terms, coefficients = build_orderings(
        fermion_hamiltonian,
        jw_hamiltonian,
        SUPPORTED_ORDERINGS,
        TOLERANCE,
        2,
    )
    return fermion_hamiltonian, jw_hamiltonian, schedules, terms, coefficients


def test_hermitian_monomial_pairing_and_nonidentity_reconstruction() -> None:
    hamiltonian = small_fermion_hamiltonian()
    terms = build_hermitian_fermion_terms(hamiltonian, TOLERANCE)

    hopping_keys = {
        ((0, 1), (1, 0)),
        ((1, 1), (0, 0)),
    }
    hopping = [term for term in terms if set(term.component_keys) == hopping_keys]
    assert len(hopping) == 1
    assert not (
        hopping[0].operator - hermitian_conjugated(hopping[0].operator)
    ).terms
    assert all(term.component_keys != ((),) for term in terms)

    reconstructed = FermionOperator()
    for term in terms:
        reconstructed += term.operator
    expected = FermionOperator()
    for key, coefficient in hamiltonian.terms.items():
        if key:
            expected += FermionOperator(key, coefficient)
    difference = clean_fermion_operator(
        reconstructed - expected,
        TOLERANCE,
    )
    assert not difference.terms


def test_all_factorizations_reconstruct_same_identity_free_jw_matrix() -> None:
    fermion_hamiltonian, jw_hamiltonian, schedules, terms, _ = all_schedules()

    assert len(terms) > 0
    assert fermion_hamiltonian.terms[()] != jw_hamiltonian.terms[()]
    for schedule in schedules.values():
        assert schedule.mapping == "jordan_wigner"
        assert schedule.number_of_factors > 0
        assert schedule.factor_reconstruction_error < 100 * TOLERANCE
        reconstructed = sum(
            (factor.hamiltonian_matrix for factor in schedule.ordered_factors),
            np.zeros((4, 4), dtype=complex),
        )
        assert np.linalg.norm(
            reconstructed - reconstructed.conjugate().T,
            ord=2,
        ) < 100 * TOLERANCE
    assert all(
        abs(np.trace(factor.matrix)) < 100 * TOLERANCE
        for factor in schedules[
            "fermionic_term_coloring"
        ].ordered_factors
    )


def test_legacy_is_pauli_permutation_and_new_method_is_not_flattened() -> None:
    _, _, schedules, terms, coefficients = all_schedules()
    raw_keys = list(coefficients)
    legacy = schedules["fermionic_coloring"]
    genuine = schedules["fermionic_term_coloring"]

    assert legacy.factorization_level == "pauli_string"
    assert len(legacy.pauli_keys) == len(set(legacy.pauli_keys))
    assert set(legacy.pauli_keys) == set(raw_keys)
    assert all(
        factor.exponential_type == "analytic_pauli"
        for factor in legacy.ordered_factors
    )

    assert genuine.factorization_level == "fermionic_term"
    assert genuine.pauli_keys == []
    assert genuine.number_of_factors == len(terms)
    assert all(
        factor.exponential_type == "dense_hermitian_expm"
        for factor in genuine.ordered_factors
    )
    # The hopping term maps to a multi-Pauli matrix but remains one factor.
    assert any(
        len(jordan_wigner(term.operator).terms) > 1
        for term in genuine.fermionic_terms
    )


def test_coloring_order_is_deterministic_and_preserves_node_order_in_color() -> None:
    graph = nx.Graph()
    graph.add_nodes_from(range(6))
    graph.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 4), (3, 5)])
    first = deterministic_coloring_order(graph)
    second = deterministic_coloring_order(graph)
    assert first[:2] == second[:2]
    assert sorted(first[0]) == list(range(6))


@pytest.mark.parametrize(
    "slice_builder",
    [first_order_slice, second_order_slice, fourth_order_slice],
)
def test_all_formula_orders_accept_general_hermitian_factors(
    slice_builder,
) -> None:
    _, _, schedules, _, _ = all_schedules()
    identity = np.eye(4, dtype=complex)
    for schedule in schedules.values():
        unitary = slice_builder(schedule.ordered_factors, 0.05, identity)
        assert np.linalg.norm(unitary.conjugate().T @ unitary - identity) < 1e-12


def test_product_multiplication_direction_detects_reversal() -> None:
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    z = np.array([[1, 0], [0, -1]], dtype=complex)
    identity = np.eye(2, dtype=complex)
    factors = [
        HamiltonianFactor(x, "dense_hermitian_expm"),
        HamiltonianFactor(z, "dense_hermitian_expm"),
    ]
    actual = first_order_slice(factors, 0.3, identity)
    chronological = expm(-0.3j * z) @ expm(-0.3j * x)
    reversed_product = expm(-0.3j * x) @ expm(-0.3j * z)
    assert np.allclose(actual, chronological)
    assert not np.allclose(actual, reversed_product)


@pytest.mark.parametrize(
    ("factor_count", "order", "steps", "expected"),
    [(4, 1, 3, 12), (4, 2, 3, 21), (4, 4, 3, 105), (3, 2, 2, 10)],
)
def test_nominal_counts_use_selected_factor_count(
    factor_count: int,
    order: int,
    steps: int,
    expected: int,
) -> None:
    assert nominal_exponential_count(factor_count, order, steps) == expected


def test_cli_default_and_explicit_subset(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["benchmark"])
    assert parse_arguments().orderings == list(SUPPORTED_ORDERINGS)

    monkeypatch.setattr(
        sys,
        "argv",
        ["benchmark", "--orderings", "fermionic_term_coloring"],
    )
    assert parse_arguments().orderings == ["fermionic_term_coloring"]
    assert normalize_requested_orderings(
        ["fermionic_term_coloring"]
    ) == ["jw_raw", "fermionic_term_coloring"]


def csv_success_row(ordering: str) -> dict[str, object]:
    row: dict[str, object] = {field: "" for field in FIELDNAMES}
    row.update(
        {
            "status": "success",
            "case_id": "case",
            "ordering": ordering,
            "formula_order": 1,
            "trotter_steps": 2,
            "evolution_time": 1.0,
            "operator_norm_error": 0.1,
            "state_infidelity": 0.01,
            "number_of_pauli_terms": 4,
            "schema_version": SCHEMA_VERSION,
        }
    )
    return row


def test_old_csv_resume_compatibility_and_migration(tmp_path: Path) -> None:
    output = tmp_path / "legacy.csv"
    old_row = csv_success_row("fermionic_coloring")
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=LEGACY_FIELDNAMES)
        writer.writeheader()
        writer.writerow({field: old_row[field] for field in LEGACY_FIELDNAMES})

    completed, errors = load_resume_data(output)
    key = ("case", "fermionic_coloring", 1, 2, 1.0)
    assert key in completed
    assert errors[key] == (0.1, 0.01)
    assert prepare_resume_output(output) == ("a", False)

    with output.open("r", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
        assert stream.seekable()
    assert list(rows[0]) == FIELDNAMES
    assert rows[0]["factorization_level"] == "pauli_string"
    assert rows[0]["ordering_method"] == "fermionic_induced_greedy_coloring"
    assert rows[0]["schema_version"] == "1"


def test_new_schema_resume_and_reordered_header_rejection(tmp_path: Path) -> None:
    current = tmp_path / "current.csv"
    with current.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerow(csv_success_row("fermionic_term_coloring"))
    assert prepare_resume_output(current) == ("a", False)
    completed, _ = load_resume_data(current)
    assert ("case", "fermionic_term_coloring", 1, 2, 1.0) in completed

    reordered = tmp_path / "reordered.csv"
    reordered_fields = FIELDNAMES[1:] + FIELDNAMES[:1]
    with reordered.open("w", newline="", encoding="utf-8") as stream:
        csv.DictWriter(stream, fieldnames=reordered_fields).writeheader()
    with pytest.raises(ValueError, match="reordered"):
        prepare_resume_output(reordered)


def plotting_rows(methods: tuple[str, ...]) -> pd.DataFrame:
    graph_values = {
        "jw_raw": (np.nan, np.nan, np.nan),
        "jw_coloring": (4, 3, 2),
        "fermionic_coloring": (3, 2, 2),
        "fermionic_term_coloring": (3, 2, 2),
    }
    rows = []
    for formula_order in (1, 2):
        for steps in (1, 2):
            for method_index, method in enumerate(methods):
                vertices, edges, colors = graph_values[method]
                rows.append(
                    {
                        "status": "success",
                        "case_id": "H-H_0.70_sto-6g_as-001-001",
                        "molecule": "H-H",
                        "bond_length": 0.70,
                        "basis": "sto-6g",
                        "active_occupied": 1,
                        "active_vacant": 1,
                        "n_qubits": 2,
                        "ordering": method,
                        "formula_order": formula_order,
                        "trotter_steps": steps,
                        "evolution_time": 1.0,
                        "state_infidelity": (method_index + 1) * 1e-4 / steps,
                        "state_vector_2norm_error": (
                            (method_index + 1) * 1e-3 / steps
                        ),
                        "graph_vertices": vertices,
                        "graph_edges": edges,
                        "number_of_colors": colors,
                        "number_of_fermionic_terms": 3,
                        "number_of_pauli_terms": 4,
                        "number_of_factors": 3 if "term" in method else 4,
                    }
                )
    return pd.DataFrame(rows)


@pytest.mark.parametrize(
    "methods",
    [
        ("jw_raw", "jw_coloring", "fermionic_coloring"),
        SUPPORTED_ORDERINGS,
    ],
)
def test_plotting_supports_old_and_new_method_sets(
    methods: tuple[str, ...],
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "results.csv"
    plotting_rows(methods).to_csv(csv_path, index=False)
    data = load_results([csv_path])
    target = summary_target(data)
    expected_target = (
        "fermionic_term_coloring"
        if "fermionic_term_coloring" in methods
        else "fermionic_coloring"
    )
    assert target == expected_target
    ranking = build_case_ranking(data, "state_infidelity", 1e-12)
    assert not ranking.empty

    figure_dir = tmp_path / "figures"
    figure_dir.mkdir()
    assert plot_case_grid(
        data,
        "H-H_0.70_sto-6g_as-001-001",
        "state_infidelity",
        1e-16,
        figure_dir,
        50,
    ).exists()
    assert plot_win_rates(
        data,
        "state_infidelity",
        1e-12,
        figure_dir,
        50,
        target,
    ).exists()
    assert plot_complexity_vs_error(
        data,
        "state_infidelity",
        1e-12,
        figure_dir,
        50,
        target,
    ).exists()


def test_plotting_tolerates_a_partial_method_set(tmp_path: Path) -> None:
    csv_path = tmp_path / "partial.csv"
    plotting_rows(("jw_raw",)).to_csv(csv_path, index=False)
    data = load_results([csv_path])
    assert summary_target(data) is None
    ranking = build_case_ranking(data, "state_infidelity", 1e-12)
    assert not ranking.empty
    figure_dir = tmp_path / "partial-figures"
    figure_dir.mkdir()
    assert plot_win_rates(
        data,
        "state_infidelity",
        1e-12,
        figure_dir,
        50,
    ).exists()
    assert plot_complexity_vs_error(
        data,
        "state_infidelity",
        1e-12,
        figure_dir,
        50,
    ).exists()


def test_first_order_error_converges_for_noncommuting_factors() -> None:
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    z = np.array([[1, 0], [0, -1]], dtype=complex)
    identity = np.eye(2, dtype=complex)
    factors = [
        HamiltonianFactor(0.4 * x, "dense_hermitian_expm"),
        HamiltonianFactor(-0.3 * z, "dense_hermitian_expm"),
    ]
    exact = expm(-1j * (0.4 * x - 0.3 * z))
    errors = []
    for steps in (1, 8, 32):
        one_slice = first_order_slice(factors, 1.0 / steps, identity)
        approximate = np.linalg.matrix_power(one_slice, steps)
        errors.append(np.linalg.norm(approximate - exact, ord=2))
    assert errors[2] < errors[1] < errors[0]
