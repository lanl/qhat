from __future__ import annotations

import math

import numpy as np
from openfermion import FermionOperator

try:
    from qhat.analysis.analyze_fermionic_graph_centralization import (
        orbital_graph_metrics,
        project_terms_to_orbital_weights,
    )
    from qhat.analysis.benchmark_L_sweep_trotter import HermitianFermionTerm
except ImportError:
    from analysis.analyze_fermionic_graph_centralization import (
        orbital_graph_metrics,
        project_terms_to_orbital_weights,
    )
    from analysis.benchmark_L_sweep_trotter import HermitianFermionTerm


def test_bipartite_freeman_uniform_partitions_are_zero() -> None:
    weights = np.array(
        [
            [0.0, 0.0, 2.0, 2.0, 2.0],
            [0.0, 0.0, 2.0, 2.0, 2.0],
            [2.0, 2.0, 0.0, 0.0, 0.0],
            [2.0, 2.0, 0.0, 0.0, 0.0],
            [2.0, 2.0, 0.0, 0.0, 0.0],
        ]
    )
    metrics = orbital_graph_metrics(weights, active_occupied=2)
    assert metrics["bipartite_weighted_freeman_centralization"] == 0.0
    assert metrics["occupied_weighted_freeman_centralization"] == 0.0
    assert metrics["virtual_weighted_freeman_centralization"] == 0.0


def test_bipartite_freeman_detects_a_hub_on_either_side() -> None:
    weights = np.zeros((6, 6), dtype=float)
    # Three occupied nodes all connect only to virtual node 3. Occupied
    # strengths are uniform, while the virtual partition is a perfect star.
    for occupied, weight in enumerate((1.0, 2.0, 4.0)):
        weights[occupied, 3] = weight
        weights[3, occupied] = weight
    metrics = orbital_graph_metrics(weights, active_occupied=3)
    assert math.isclose(
        metrics["bipartite_weighted_freeman_centralization"],
        1.0,
    )
    assert math.isclose(metrics["virtual_weighted_freeman_centralization"], 1.0)


def test_term_projection_builds_occupied_virtual_edges() -> None:
    hopping = HermitianFermionTerm(
        index=0,
        first_raw_index=0,
        operator=(
            FermionOperator(((0, 1), (1, 0)), 2.0)
            + FermionOperator(((1, 1), (0, 0)), 2.0)
        ),
        component_keys=(((0, 1), (1, 0)), ((1, 1), (0, 0))),
    )
    four_orbital = HermitianFermionTerm(
        index=1,
        first_raw_index=2,
        operator=(
            FermionOperator(((0, 1), (1, 1), (2, 0), (3, 0)), 3.0)
            + FermionOperator(((3, 1), (2, 1), (1, 0), (0, 0)), 3.0)
        ),
        component_keys=(
            ((0, 1), (1, 1), (2, 0), (3, 0)),
            ((3, 1), (2, 1), (1, 0), (0, 0)),
        ),
    )
    onsite = HermitianFermionTerm(
        index=2,
        first_raw_index=4,
        operator=FermionOperator(((2, 1), (2, 0)), 5.0),
        component_keys=(((2, 1), (2, 0)),),
    )

    weights, metadata = project_terms_to_orbital_weights(
        [hopping, four_orbital, onsite],
        n_qubits=4,
        active_occupied=2,
        tolerance=1.0e-12,
    )

    # The hopping parent is occupied-only and therefore excluded. The
    # four-orbital parent contributes its canonical |c|=3 to all four O-V
    # support pairs. The onsite parent is also noncrossing.
    assert math.isclose(weights[0, 2], 3.0)
    assert math.isclose(weights[1, 3], 3.0)
    assert math.isclose(np.triu(weights, k=1).sum(), 12.0)
    assert math.isclose(metadata["total_parent_coefficient_mass"], 10.0)
    assert math.isclose(metadata["mixed_parent_coefficient_mass"], 3.0)
    assert metadata["number_of_mixed_occupied_virtual_parents"] == 1
    assert metadata["number_of_noncrossing_parents"] == 2
