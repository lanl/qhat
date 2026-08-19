#!/usr/bin/env python3
"""Pre-BCH structural diagnostics for fermionic-aware Trotter ordering.

Purpose
-------
Identify *pre-BCH* Hamiltonian/graph features that distinguish the favorable
H2O/NH3 signed-fermionic parent ordering from the unfavorable O2 8+10 HGBS-5
case, with BeH2 used as an intermediate/check case.

The predictor features in this script use only information available before a
BCH error is evaluated:

* complete Hermitian fermionic-parent coefficients and coefficient signs,
* one- vs two-body structure,
* excitation rank relative to the HF occupied/virtual split,
* orbital-index overlap between fermionic parents,
* coefficient-weighted parent connectivity,
* occupied/virtual orbital-incidence concentration,
* low-rank/concentration structure of occupied-to-virtual parent coupling,
* graph degree/strength, components, and signed-parent-order locality.

No commutator norm, BCH vector, BCH cancellation metric, or Trotter error is
used to construct any structural feature.  Trotter ``one_minus_overlap`` is
used only as an *external outcome label* when a deterministic-results CSV is
available, so that structural features can be correlated with whether
``fermionic_signed_reference`` beats the fair JW baseline

    min(JW signed-coefficient, JW magnitude-descending).

Run from the QHAT repository root, for example:

    python analysis/analyze_pre_bch_structure.pyh

or with explicit result files:

    python analysis/analyze_pre_bch_structure.pyh \
      --results-csv analysis/polyatomic_reference_deterministic_ordering.csv \
      --results-csv analysis/deterministic_orderings_hgbs5_results.csv

Default outputs are written to ``analysis/pre_bch_structure/``:

    pre_bch_case_features.csv
    pre_bch_feature_ranking.csv
    pre_bch_parent_type_summary.csv
    pre_bch_report.md
    pre_bch_top_features.png          (when matplotlib is available)

Implementation note
-------------------
The input ``*.tensors.npz`` files contain the number-conserving one- and
Two-body InteractionOperator tensors.  To keep this analysis lightweight and
strictly pre-BCH, the script reconstructs the normal-ordered rank-1/rank-2
fermionic monomials directly from those tensors, pairs each monomial with its
Hermitian adjoint, and defines the same canonical signed parent coefficient
used by the deterministic ordering code: the coefficient of the
lexicographically smallest component of each complete Hermitian parent.

For the current four cases this direct reconstruction is intentionally checked
against the expected parent counts when known.  It does not require building a
JW Hamiltonian or a BCH evaluator.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np


DEFAULT_TOLERANCE = 1.0e-12
DEFAULT_STEPS = 100
DEFAULT_TIME = 1.0
DEFAULT_OUTPUT_DIR = Path("analysis/pre_bch_structure")

# The exact four HGBS-5 cases requested for this mechanism check.
DEFAULT_CASES = (
    (
        "H2O",
        Path(
            "hamiltonian_generator/polyatomic_library/H2O/s-1.00/hgbs-5/"
            "H2O_s-1.00_hgbs-5_as-008-010.tensors.npz"
        ),
        8,
        10,
        1367,
        "favorable_anchor",
    ),
    (
        "NH3",
        Path(
            "hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/"
            "NH3_s-1.00_hgbs-5_as-008-010.tensors.npz"
        ),
        8,
        10,
        3899,
        "favorable_anchor",
    ),
    (
        "BeH2",
        Path(
            "hamiltonian_generator/polyatomic_library/BeH2/s-1.00/hgbs-5/"
            "BeH2_s-1.00_hgbs-5_as-006-012.tensors.npz"
        ),
        6,
        12,
        1131,
        "intermediate_check",
    ),
    (
        "O2",
        Path(
            "hamiltonian_generator/o2_active_space_library/O-O/1.26/hgbs-5/"
            "O-O_1.26_hgbs-5_as-008-010.tensors.npz"
        ),
        8,
        10,
        1293,
        "unfavorable_anchor",
    ),
)

FERMIONIC_ALIASES = {
    "fermionic_signed_reference",
    "fermionic_signed_coefficient_lexicographic",
}
JW_SIGNED_ALIASES = {
    "jw_signed_baseline",
    "signed_coefficient_lexicographic",
}
JW_MAGNITUDE_ALIASES = {
    "jw_magnitude_descending_lexicographic",
}

# Human-readable descriptions used in the ranking/report.  Features not listed
# here still appear in the CSV; they simply get a generated description.
FEATURE_DESCRIPTIONS: dict[str, str] = {
    "coefficient_effective_fraction": (
        "Participation-ratio effective number of |parent coefficients| divided by parent count"
    ),
    "coefficient_gini": "Gini concentration of absolute fermionic-parent coefficients",
    "coefficient_entropy": "Normalized Shannon entropy of absolute fermionic-parent coefficients",
    "coefficient_top1_mass_fraction": "Fraction of total |coefficient| mass in the largest parent",
    "coefficient_top10_mass_fraction": "Fraction of total |coefficient| mass in the 10 largest parents",
    "coefficient_top1pct_mass_fraction": "Fraction of total |coefficient| mass in the largest 1% of parents",
    "positive_parent_fraction": "Fraction of canonical parent coefficients that are positive",
    "coefficient_sign_mass_imbalance": (
        "Absolute imbalance between positive and negative |coefficient| mass"
    ),
    "excitation_rank2_coefficient_mass_fraction": (
        "Fraction of |coefficient| mass in parents transferring two particles across the occupied/virtual split"
    ),
    "excitation_rank1_coefficient_mass_fraction": (
        "Fraction of |coefficient| mass in parents transferring one particle across the occupied/virtual split"
    ),
    "diagonal_coefficient_mass_fraction": (
        "Fraction of |coefficient| mass in occupation-basis diagonal parents"
    ),
    "mixed_occ_virtual_coefficient_mass_fraction": (
        "Fraction of |coefficient| mass in parents touching both occupied and virtual orbitals"
    ),
    "overlap_edge_density": "Density of the parent graph whose edges share at least one orbital index",
    "degree_cv": "Coefficient of variation of parent overlap-graph degree",
    "degree_gini": "Gini concentration of parent overlap-graph degree",
    "max_degree_fraction": "Largest parent degree divided by N-1",
    "shared_index_strength_cv": "Coefficient of variation of total shared-index count per parent",
    "shared_index_strength_gini": "Gini concentration of total shared-index count per parent",
    "shared2_edge_fraction": "Fraction of overlap edges whose parents share at least two orbital indices",
    "same_excitation_rank_edge_fraction": "Fraction of overlap edges connecting the same excitation-rank class",
    "same_sign_overlap_edge_fraction": "Unweighted fraction of overlap edges joining same-sign parent coefficients",
    "same_sign_coeffproduct_edge_fraction": (
        "|c_i c_j|-weighted fraction of overlap edges joining same-sign parent coefficients"
    ),
    "largest_component_fraction": "Fraction of parents in the largest orbital-overlap connected component",
    "component_count_fraction": "Connected-component count divided by parent count",
    "signed_order_mean_edge_span_fraction": (
        "Mean parent-order separation of overlap edges, normalized by N-1"
    ),
    "signed_order_weighted_edge_span_fraction": (
        "|c_i c_j|-weighted mean separation of overlap edges in signed parent order"
    ),
    "signed_order_edge_within_1pct_fraction": (
        "Fraction of overlap edges whose endpoints lie within 1% of the signed parent sequence"
    ),
    "signed_order_edge_within_5pct_fraction": (
        "Fraction of overlap edges whose endpoints lie within 5% of the signed parent sequence"
    ),
    "signed_order_edge_within_10pct_fraction": (
        "Fraction of overlap edges whose endpoints lie within 10% of the signed parent sequence"
    ),
    "virtual_orbital_coefficient_mass_gini": (
        "Gini concentration of parent |coefficient| incidence across virtual orbitals"
    ),
    "virtual_orbital_coefficient_mass_entropy": (
        "Normalized entropy of parent |coefficient| incidence across virtual orbitals"
    ),
    "virtual_orbital_coefficient_mass_top1_fraction": (
        "Largest virtual-orbital share of coefficient-weighted parent incidence"
    ),
    "occupied_orbital_coefficient_mass_gini": (
        "Gini concentration of parent |coefficient| incidence across occupied orbitals"
    ),
    "ov_coupling_top_singular_fraction": (
        "Largest singular-value share of the coefficient-weighted occupied-to-virtual incidence network"
    ),
    "ov_coupling_effective_rank_fraction": (
        "Participation-ratio effective rank of occupied-to-virtual incidence coupling, normalized by matrix rank limit"
    ),
    "parent_orbital_incidence_top_singular_fraction": (
        "Largest singular-value share of the |coefficient|-weighted parent-orbital incidence matrix"
    ),
    "parent_orbital_incidence_effective_rank_fraction": (
        "Effective rank fraction of the |coefficient|-weighted parent-orbital incidence matrix"
    ),
}

# Features that make the most physically interpretable pre-BCH report.  The
# complete ranking CSV still contains every finite scalar feature.
REPORT_FEATURE_WHITELIST = {
    "coefficient_effective_fraction",
    "coefficient_gini",
    "coefficient_entropy",
    "coefficient_top1_mass_fraction",
    "coefficient_top10_mass_fraction",
    "coefficient_sign_mass_imbalance",
    "excitation_rank1_coefficient_mass_fraction",
    "excitation_rank2_coefficient_mass_fraction",
    "diagonal_coefficient_mass_fraction",
    "mixed_occ_virtual_coefficient_mass_fraction",
    "overlap_edge_density",
    "degree_cv",
    "degree_gini",
    "shared_index_strength_cv",
    "shared2_edge_fraction",
    "same_excitation_rank_edge_fraction",
    "same_sign_coeffproduct_edge_fraction",
    "largest_component_fraction",
    "signed_order_mean_edge_span_fraction",
    "signed_order_weighted_edge_span_fraction",
    "signed_order_edge_within_1pct_fraction",
    "signed_order_edge_within_5pct_fraction",
    "virtual_orbital_coefficient_mass_gini",
    "virtual_orbital_coefficient_mass_entropy",
    "virtual_orbital_coefficient_mass_top1_fraction",
    "occupied_orbital_coefficient_mass_gini",
    "ov_coupling_top_singular_fraction",
    "ov_coupling_effective_rank_fraction",
    "parent_orbital_incidence_top_singular_fraction",
    "parent_orbital_incidence_effective_rank_fraction",
}


FermionKey = tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class CaseSpec:
    molecule: str
    tensor_path: Path
    active_occupied: int
    active_vacant: int
    expected_parent_count: int | None
    contrast_role: str

    @property
    def n_qubits(self) -> int:
        return self.active_occupied + self.active_vacant

    @property
    def case_id(self) -> str:
        name = self.tensor_path.name
        if name.endswith(".tensors.npz"):
            return name[: -len(".tensors.npz")]
        return self.tensor_path.stem


@dataclass(frozen=True)
class FermionicParent:
    components: tuple[FermionKey, ...]
    canonical_coefficient: float
    support: tuple[int, ...]
    occupied_support: tuple[int, ...]
    virtual_support: tuple[int, ...]
    body_rank: int
    excitation_rank: int
    diagonal: bool

    @property
    def canonical_key(self) -> FermionKey:
        return min(self.components)


@dataclass
class Outcome:
    case_id: str
    fermionic_error: float | None = None
    jw_signed_error: float | None = None
    jw_magnitude_error: float | None = None
    best_jw_error: float | None = None
    best_jw_ordering: str = ""
    advantage: float | None = None
    log10_advantage: float | None = None
    classification: str = "unavailable"
    source: str = ""
    baseline_complete: bool = False


@dataclass(frozen=True)
class FeatureRankingRow:
    feature: str
    description: str
    category: str
    h2o: float
    nh3: float
    beh2: float
    o2: float
    favorable_mean: float
    favorable_std: float
    o2_minus_favorable: float
    symmetric_relative_gap: float
    range_normalized_gap: float
    favorable_consistency: float
    anchor_separation_score: float
    spearman_vs_log10_advantage: float
    combined_score: float
    direction: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute pre-BCH fermionic-parent structural features for the HGBS-5 "
            "H2O/NH3/BeH2/O2 mechanism cases."
        )
    )
    parser.add_argument(
        "--tensor-root",
        type=Path,
        default=Path("."),
        help=(
            "Repository root under which default tensor paths are resolved "
            "(default: current directory)."
        ),
    )
    parser.add_argument(
        "--results-csv",
        type=Path,
        action="append",
        default=[],
        help=(
            "Deterministic-ordering/performance CSV used only for external outcome labels. "
            "May be repeated. If omitted, likely analysis CSVs are auto-discovered."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=DEFAULT_STEPS,
        help="Trotter-step setting used only when selecting outcome-label rows (default: 100).",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=DEFAULT_TIME,
        dest="evolution_time",
        help="Evolution time used only when selecting outcome-label rows (default: 1.0).",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Tensor/coefficient tolerance (default: 1e-12).",
    )
    parser.add_argument(
        "--floor",
        type=float,
        default=1.0e-12,
        help="Outcome numerical floor below which win/loss is not classified (default: 1e-12).",
    )
    parser.add_argument(
        "--tie-rtol",
        type=float,
        default=1.0e-3,
        help="Relative tolerance around advantage=1 used to call an outcome a tie (default: 1e-3).",
    )
    parser.add_argument(
        "--allow-incomplete-jw-baseline",
        action="store_true",
        help=(
            "Allow outcome labeling when only one of JW-signed/JW-magnitude is present. "
            "By default the fair best-JW baseline requires both."
        ),
    )
    parser.add_argument(
        "--top-features",
        type=int,
        default=12,
        help="Number of top interpretable features printed/written in the report (default: 12).",
    )
    parser.add_argument(
        "--manual-advantage",
        action="append",
        default=[],
        metavar="CASE_ID=RATIO",
        help=(
            "Optional externally computed best-JW/fermionic advantage override. "
            "Useful when the latest outcome table is not in a long-form CSV. May be repeated."
        ),
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run lightweight algebra/metric self-tests and exit.",
    )
    return parser.parse_args()


def gini(values: Sequence[float] | np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = np.abs(x[np.isfinite(x)])
    if x.size == 0:
        return float("nan")
    total = float(x.sum())
    if total <= 0.0:
        return 0.0
    y = np.sort(x)
    n = y.size
    indices = np.arange(1, n + 1, dtype=float)
    return float(2.0 * np.dot(indices, y) / (n * total) - (n + 1.0) / n)


def normalized_entropy(values: Sequence[float] | np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x) & (x > 0.0)]
    if x.size <= 1:
        return 0.0
    p = x / x.sum()
    return float(-np.dot(p, np.log(p)) / math.log(x.size))


def effective_count(values: Sequence[float] | np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x) & (x >= 0.0)]
    total = float(x.sum())
    square_sum = float(np.dot(x, x))
    if square_sum <= 0.0:
        return 0.0
    return total * total / square_sum


def inversion_parity_for_descending(values: Sequence[int]) -> int:
    """Return +1/-1 parity for sorting values into descending order."""
    swaps = 0
    values = list(values)
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            if values[i] < values[j]:
                swaps += 1
    return -1 if swaps % 2 else 1


def canonicalize_number_conserving_key(key: FermionKey) -> tuple[FermionKey | None, int]:
    """Canonicalize a rank-1/rank-2 number-conserving monomial.

    InteractionOperator terms already place all creation operators before all
    annihilation operators. Normal ordering therefore only needs the canonical
    within-creation and within-annihilation index order. Repeated equal-mode
    creation or annihilation operators vanish by fermionic antisymmetry.
    """
    creators = [index for index, action in key if action == 1]
    annihilators = [index for index, action in key if action == 0]
    if len(creators) != len(annihilators):
        raise ValueError(f"Expected number-conserving term, got {key!r}")
    if len(set(creators)) != len(creators) or len(set(annihilators)) != len(annihilators):
        return None, 0

    sign = inversion_parity_for_descending(creators)
    sign *= inversion_parity_for_descending(annihilators)
    creators = sorted(creators, reverse=True)
    annihilators = sorted(annihilators, reverse=True)
    canonical = tuple((index, 1) for index in creators) + tuple(
        (index, 0) for index in annihilators
    )
    return canonical, sign


def adjoint_canonical_key(key: FermionKey) -> tuple[FermionKey | None, int]:
    raw_adjoint = tuple((index, 1 - action) for index, action in reversed(key))
    return canonicalize_number_conserving_key(raw_adjoint)


def real_coefficient(value: complex, tolerance: float, context: str) -> float:
    value = complex(value)
    scale = max(1.0, abs(value.real))
    if abs(value.imag) > max(100.0 * tolerance, 1.0e-10) * scale:
        raise ValueError(f"Non-negligible imaginary coefficient in {context}: {value!r}")
    return float(value.real)


def reconstruct_normal_ordered_terms(
    tensor_path: Path,
    tolerance: float,
) -> tuple[dict[FermionKey, complex], int]:
    with np.load(tensor_path) as arrays:
        one_body = np.asarray(arrays["one_body"])
        two_body = np.asarray(arrays["two_body"])

    if one_body.ndim != 2 or one_body.shape[0] != one_body.shape[1]:
        raise ValueError(f"Unexpected one_body shape in {tensor_path}: {one_body.shape}")
    n_qubits = int(one_body.shape[0])
    if two_body.shape != (n_qubits, n_qubits, n_qubits, n_qubits):
        raise ValueError(f"Unexpected two_body shape in {tensor_path}: {two_body.shape}")

    terms: dict[FermionKey, complex] = {}

    for p, q in np.argwhere(np.abs(one_body) > tolerance):
        key: FermionKey = ((int(p), 1), (int(q), 0))
        terms[key] = terms.get(key, 0.0j) + complex(one_body[p, q])

    for p, q, r, s in np.argwhere(np.abs(two_body) > tolerance):
        raw_key: FermionKey = (
            (int(p), 1),
            (int(q), 1),
            (int(r), 0),
            (int(s), 0),
        )
        key, sign = canonicalize_number_conserving_key(raw_key)
        if key is None:
            continue
        terms[key] = terms.get(key, 0.0j) + sign * complex(two_body[p, q, r, s])

    return ({key: value for key, value in terms.items() if abs(value) > tolerance}, n_qubits)


def reconstruct_hermitian_parents(
    tensor_path: Path,
    active_occupied: int,
    tolerance: float,
) -> tuple[list[FermionicParent], int, int]:
    terms, n_qubits = reconstruct_normal_ordered_terms(tensor_path, tolerance)
    visited: set[FermionKey] = set()
    parents: list[FermionicParent] = []

    for key, coefficient in terms.items():
        if key in visited:
            continue
        adjoint_key, adjoint_phase = adjoint_canonical_key(key)
        if adjoint_key is None:
            raise ValueError(f"Unexpected vanishing adjoint for nonzero key {key!r}")
        if adjoint_key not in terms:
            raise ValueError(
                f"Hermitian partner missing in {tensor_path.name}: {key!r} -> {adjoint_key!r}"
            )

        expected_partner = np.conjugate(coefficient) * adjoint_phase
        actual_partner = terms[adjoint_key]
        hermitian_error = abs(actual_partner - expected_partner)
        hermitian_scale = max(1.0, abs(actual_partner), abs(expected_partner))
        if hermitian_error > max(100.0 * tolerance, 1.0e-10) * hermitian_scale:
            raise ValueError(
                f"Hermitian-pair mismatch in {tensor_path.name}: {key!r}; "
                f"expected {expected_partner!r}, got {actual_partner!r}"
            )

        components = (key,) if adjoint_key == key else tuple(sorted((key, adjoint_key)))
        canonical_key = min(components)
        canonical_coefficient = real_coefficient(
            terms[canonical_key],
            tolerance,
            f"canonical parent {canonical_key!r}",
        )

        creators = [index for index, action in canonical_key if action == 1]
        annihilators = [index for index, action in canonical_key if action == 0]
        support = tuple(sorted(set(creators + annihilators)))
        occupied_support = tuple(index for index in support if index < active_occupied)
        virtual_support = tuple(index for index in support if index >= active_occupied)

        # Excitation rank is the absolute net transfer across the HF occupied /
        # virtual partition.  For number-conserving rank-1/rank-2 parents this
        # is 0, 1, or 2 and is invariant under taking the Hermitian adjoint.
        virtual_creators = sum(index >= active_occupied for index in creators)
        virtual_annihilators = sum(index >= active_occupied for index in annihilators)
        excitation_rank = abs(virtual_creators - virtual_annihilators)
        diagonal = sorted(creators) == sorted(annihilators)

        parents.append(
            FermionicParent(
                components=components,
                canonical_coefficient=canonical_coefficient,
                support=support,
                occupied_support=occupied_support,
                virtual_support=virtual_support,
                body_rank=len(creators),
                excitation_rank=excitation_rank,
                diagonal=diagonal,
            )
        )
        visited.update(components)

    # This is the actual signed-parent sequence used by the deterministic
    # fermionic-signed ordering: ascending canonical signed coefficient with
    # a lexicographic fermionic-term tie break.
    parents.sort(key=lambda parent: (parent.canonical_coefficient, parent.components))
    return parents, n_qubits, len(terms)


class DisjointSet:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        a = self.find(left)
        b = self.find(right)
        if a == b:
            return
        if self.rank[a] < self.rank[b]:
            a, b = b, a
        self.parent[b] = a
        if self.rank[a] == self.rank[b]:
            self.rank[a] += 1


def connected_component_sizes_from_incidence(incidence: np.ndarray) -> list[int]:
    n_parents, n_orbitals = incidence.shape
    dsu = DisjointSet(n_parents)
    for orbital in range(n_orbitals):
        members = np.flatnonzero(incidence[:, orbital])
        if members.size <= 1:
            continue
        anchor = int(members[0])
        for member in members[1:]:
            dsu.union(anchor, int(member))
    counts: dict[int, int] = {}
    for parent_index in range(n_parents):
        root = dsu.find(parent_index)
        counts[root] = counts.get(root, 0) + 1
    return sorted(counts.values(), reverse=True)


def safe_fraction(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator > 0.0 else 0.0


def coefficient_mass_fraction(mask: np.ndarray, absolute_coefficients: np.ndarray) -> float:
    return safe_fraction(float(absolute_coefficients[mask].sum()), float(absolute_coefficients.sum()))


def compute_case_features(
    case: CaseSpec,
    tensor_path: Path,
    tolerance: float,
) -> tuple[dict[str, float | int | str], list[dict[str, float | int | str]]]:
    parents, n_qubits, normal_ordered_count = reconstruct_hermitian_parents(
        tensor_path,
        case.active_occupied,
        tolerance,
    )
    if n_qubits != case.n_qubits:
        raise ValueError(
            f"{case.case_id}: tensor has {n_qubits} orbitals/qubits but path metadata implies {case.n_qubits}"
        )

    n_parents = len(parents)
    coefficients = np.asarray([parent.canonical_coefficient for parent in parents], dtype=float)
    absolute_coefficients = np.abs(coefficients)
    signs = np.sign(coefficients)
    coefficient_total = float(absolute_coefficients.sum())

    incidence = np.zeros((n_parents, n_qubits), dtype=np.uint8)
    for parent_index, parent in enumerate(parents):
        incidence[parent_index, list(parent.support)] = 1

    # Parent graph: edge iff two parents share at least one orbital index.
    # This is deliberately a structural proxy, not a commutator graph.
    shared_indices = incidence @ incidence.T
    np.fill_diagonal(shared_indices, 0)
    adjacency = shared_indices > 0
    edge_count = float(adjacency.sum()) / 2.0
    degree = adjacency.sum(axis=1).astype(float)
    shared_strength = shared_indices.sum(axis=1).astype(float)

    features: dict[str, float | int | str] = {
        "case_id": case.case_id,
        "molecule": case.molecule,
        "tensor_path": str(case.tensor_path),
        "basis": "hgbs-5",
        "active_occupied": case.active_occupied,
        "active_vacant": case.active_vacant,
        "n_qubits": n_qubits,
        "contrast_role": case.contrast_role,
        "number_of_fermionic_parents": n_parents,
        "number_of_normal_ordered_monomials": normal_ordered_count,
        "expected_parent_count": case.expected_parent_count if case.expected_parent_count is not None else "",
        "parent_count_delta": (
            n_parents - case.expected_parent_count if case.expected_parent_count is not None else ""
        ),
    }

    # Coefficient distribution.
    features["coefficient_effective_fraction"] = effective_count(absolute_coefficients) / max(1, n_parents)
    features["coefficient_gini"] = gini(absolute_coefficients)
    features["coefficient_entropy"] = normalized_entropy(absolute_coefficients)
    sorted_coefficients = np.sort(absolute_coefficients)
    features["coefficient_top1_mass_fraction"] = safe_fraction(
        float(sorted_coefficients[-1:].sum()), coefficient_total
    )
    features["coefficient_top10_mass_fraction"] = safe_fraction(
        float(sorted_coefficients[-min(10, n_parents) :].sum()), coefficient_total
    )
    top_1pct_count = max(1, math.ceil(0.01 * n_parents))
    features["coefficient_top1pct_mass_fraction"] = safe_fraction(
        float(sorted_coefficients[-top_1pct_count:].sum()), coefficient_total
    )
    positive_mass = float(absolute_coefficients[coefficients > 0.0].sum())
    negative_mass = float(absolute_coefficients[coefficients < 0.0].sum())
    features["positive_parent_fraction"] = float(np.mean(coefficients > 0.0))
    features["coefficient_sign_mass_imbalance"] = safe_fraction(
        abs(positive_mass - negative_mass), coefficient_total
    )

    # Parent type/excitation composition.
    body_ranks = np.asarray([parent.body_rank for parent in parents], dtype=int)
    excitation_ranks = np.asarray([parent.excitation_rank for parent in parents], dtype=int)
    diagonal_mask = np.asarray([parent.diagonal for parent in parents], dtype=bool)
    mixed_occ_virtual = np.asarray(
        [bool(parent.occupied_support and parent.virtual_support) for parent in parents],
        dtype=bool,
    )
    support_sizes = incidence.sum(axis=1).astype(int)

    for rank in (1, 2):
        mask = body_ranks == rank
        features[f"body_rank{rank}_parent_fraction"] = float(mask.mean())
        features[f"body_rank{rank}_coefficient_mass_fraction"] = coefficient_mass_fraction(
            mask, absolute_coefficients
        )
    for rank in (0, 1, 2):
        mask = excitation_ranks == rank
        features[f"excitation_rank{rank}_parent_fraction"] = float(mask.mean())
        features[f"excitation_rank{rank}_coefficient_mass_fraction"] = coefficient_mass_fraction(
            mask, absolute_coefficients
        )
    features["diagonal_parent_fraction"] = float(diagonal_mask.mean())
    features["diagonal_coefficient_mass_fraction"] = coefficient_mass_fraction(
        diagonal_mask, absolute_coefficients
    )
    features["mixed_occ_virtual_parent_fraction"] = float(mixed_occ_virtual.mean())
    features["mixed_occ_virtual_coefficient_mass_fraction"] = coefficient_mass_fraction(
        mixed_occ_virtual, absolute_coefficients
    )
    features["mean_parent_support_size"] = float(support_sizes.mean())
    for size in (1, 2, 3, 4):
        mask = support_sizes == size
        features[f"support_size{size}_parent_fraction"] = float(mask.mean())
        features[f"support_size{size}_coefficient_mass_fraction"] = coefficient_mass_fraction(
            mask, absolute_coefficients
        )

    # Parent orbital-overlap graph/network structure.
    possible_edges = n_parents * (n_parents - 1) / 2.0
    features["overlap_edge_density"] = safe_fraction(edge_count, possible_edges)
    features["mean_degree_fraction"] = safe_fraction(float(degree.mean()), max(1, n_parents - 1))
    features["degree_cv"] = safe_fraction(float(degree.std()), float(degree.mean()))
    features["degree_gini"] = gini(degree)
    features["max_degree_fraction"] = safe_fraction(float(degree.max()), max(1, n_parents - 1))
    features["shared_index_strength_cv"] = safe_fraction(
        float(shared_strength.std()), float(shared_strength.mean())
    )
    features["shared_index_strength_gini"] = gini(shared_strength)
    shared2_edges = float(np.count_nonzero(shared_indices >= 2)) / 2.0
    features["shared2_edge_fraction"] = safe_fraction(shared2_edges, edge_count)

    same_rank_edges = 0.0
    for rank in (0, 1, 2):
        indicator = (excitation_ranks == rank).astype(float)
        same_rank_edges += float(indicator @ adjacency @ indicator) / 2.0
    features["same_excitation_rank_edge_fraction"] = safe_fraction(same_rank_edges, edge_count)

    sign_quadratic = float(signs @ adjacency @ signs) / 2.0
    same_sign_edges = (edge_count + sign_quadratic) / 2.0
    features["same_sign_overlap_edge_fraction"] = safe_fraction(same_sign_edges, edge_count)

    total_coefficient_edge_weight = float(absolute_coefficients @ adjacency @ absolute_coefficients) / 2.0
    signed_coefficients = absolute_coefficients * signs
    signed_edge_weight = float(signed_coefficients @ adjacency @ signed_coefficients) / 2.0
    same_sign_coefficient_edge_weight = (total_coefficient_edge_weight + signed_edge_weight) / 2.0
    features["same_sign_coeffproduct_edge_fraction"] = safe_fraction(
        same_sign_coefficient_edge_weight, total_coefficient_edge_weight
    )
    features["overlap_coeffproduct_weight_per_parent"] = safe_fraction(
        total_coefficient_edge_weight, n_parents
    )

    component_sizes = connected_component_sizes_from_incidence(incidence)
    features["number_of_components"] = len(component_sizes)
    features["component_count_fraction"] = safe_fraction(len(component_sizes), n_parents)
    features["largest_component_fraction"] = safe_fraction(component_sizes[0], n_parents)

    # Locality in the signed coefficient parent order.  Parent list is already
    # sorted by the same signed-parent key used by the deterministic ordering.
    edge_span_sum = 0.0
    weighted_edge_span_sum = 0.0
    edge_weight_sum = 0.0
    local_edge_counts: dict[int, float] = {}
    local_weight_counts: dict[int, float] = {}
    cutoffs = {
        max(1, int(math.ceil(0.01 * n_parents))),
        max(1, int(math.ceil(0.05 * n_parents))),
        max(1, int(math.ceil(0.10 * n_parents))),
    }
    for cutoff in cutoffs:
        local_edge_counts[cutoff] = 0.0
        local_weight_counts[cutoff] = 0.0

    for separation in range(1, n_parents):
        edge_mask = np.diagonal(adjacency, offset=separation)
        count = float(edge_mask.sum())
        if count <= 0.0:
            continue
        edge_span_sum += separation * count
        pair_weights = absolute_coefficients[:-separation] * absolute_coefficients[separation:]
        selected_weight = float(pair_weights[edge_mask].sum())
        weighted_edge_span_sum += separation * selected_weight
        edge_weight_sum += selected_weight
        for cutoff in cutoffs:
            if separation <= cutoff:
                local_edge_counts[cutoff] += count
                local_weight_counts[cutoff] += selected_weight

    features["signed_order_mean_edge_span_fraction"] = safe_fraction(
        edge_span_sum, edge_count * max(1, n_parents - 1)
    )
    features["signed_order_weighted_edge_span_fraction"] = safe_fraction(
        weighted_edge_span_sum, edge_weight_sum * max(1, n_parents - 1)
    )
    for percent in (1, 5, 10):
        cutoff = max(1, int(math.ceil((percent / 100.0) * n_parents)))
        features[f"signed_order_edge_within_{percent}pct_fraction"] = safe_fraction(
            local_edge_counts[cutoff], edge_count
        )
        features[f"signed_order_weighted_edge_within_{percent}pct_fraction"] = safe_fraction(
            local_weight_counts[cutoff], edge_weight_sum
        )

    # Orbital concentration.  A parent contributes |c_i| to every orbital it
    # touches; this detects coefficient hubs without evaluating a commutator.
    orbital_incidence = incidence.sum(axis=0).astype(float)
    orbital_coefficient_mass = absolute_coefficients @ incidence
    features["orbital_incidence_gini"] = gini(orbital_incidence)
    features["orbital_coefficient_mass_gini"] = gini(orbital_coefficient_mass)
    features["orbital_coefficient_mass_entropy"] = normalized_entropy(orbital_coefficient_mass)
    features["orbital_coefficient_mass_top1_fraction"] = safe_fraction(
        float(orbital_coefficient_mass.max()), float(orbital_coefficient_mass.sum())
    )

    orbital_blocks = {
        "occupied": slice(0, case.active_occupied),
        "virtual": slice(case.active_occupied, n_qubits),
    }
    for label, block in orbital_blocks.items():
        block_incidence = orbital_incidence[block]
        block_coefficient_mass = orbital_coefficient_mass[block]
        features[f"{label}_orbital_incidence_gini"] = gini(block_incidence)
        features[f"{label}_orbital_coefficient_mass_gini"] = gini(block_coefficient_mass)
        features[f"{label}_orbital_coefficient_mass_entropy"] = normalized_entropy(
            block_coefficient_mass
        )
        features[f"{label}_orbital_coefficient_mass_top1_fraction"] = safe_fraction(
            float(block_coefficient_mass.max()), float(block_coefficient_mass.sum())
        )

    # Occupied-to-virtual structural coupling matrix.  Each parent contributes
    # |c_i| to every occupied/virtual support pair it jointly touches.  SVD
    # concentration asks whether that structural coupling is spread over many
    # channels (high effective rank) or dominated by a few hubs (low rank).
    occupied_incidence = incidence[:, : case.active_occupied].astype(float)
    virtual_incidence = incidence[:, case.active_occupied :].astype(float)
    ov_coupling = occupied_incidence.T @ (absolute_coefficients[:, None] * virtual_incidence)
    ov_singular_values = np.linalg.svd(ov_coupling, compute_uv=False)
    ov_singular_sum = float(ov_singular_values.sum())
    ov_rank_limit = max(1, min(ov_coupling.shape))
    features["ov_coupling_frobenius_norm"] = float(np.linalg.norm(ov_coupling))
    features["ov_coupling_top_singular_fraction"] = safe_fraction(
        float(ov_singular_values[0]) if ov_singular_values.size else 0.0,
        ov_singular_sum,
    )
    features["ov_coupling_singular_entropy"] = normalized_entropy(ov_singular_values)
    features["ov_coupling_effective_rank"] = effective_count(ov_singular_values)
    features["ov_coupling_effective_rank_fraction"] = (
        effective_count(ov_singular_values) / ov_rank_limit
    )

    # A second low-rank diagnostic on the full parent-orbital incidence matrix.
    weighted_incidence = np.sqrt(absolute_coefficients)[:, None] * incidence.astype(float)
    incidence_singular_values = np.linalg.svd(weighted_incidence, compute_uv=False)
    incidence_singular_sum = float(incidence_singular_values.sum())
    incidence_rank_limit = max(1, min(weighted_incidence.shape))
    features["parent_orbital_incidence_top_singular_fraction"] = safe_fraction(
        float(incidence_singular_values[0]), incidence_singular_sum
    )
    features["parent_orbital_incidence_effective_rank_fraction"] = (
        effective_count(incidence_singular_values) / incidence_rank_limit
    )
    features["parent_orbital_incidence_singular_entropy"] = normalized_entropy(
        incidence_singular_values
    )

    parent_type_rows: list[dict[str, float | int | str]] = []
    for body_rank in (1, 2):
        for excitation_rank in (0, 1, 2):
            mask = (body_ranks == body_rank) & (excitation_ranks == excitation_rank)
            if not np.any(mask):
                continue
            parent_type_rows.append(
                {
                    "case_id": case.case_id,
                    "molecule": case.molecule,
                    "body_rank": body_rank,
                    "excitation_rank": excitation_rank,
                    "parent_count": int(mask.sum()),
                    "parent_fraction": float(mask.mean()),
                    "absolute_coefficient_mass": float(absolute_coefficients[mask].sum()),
                    "coefficient_mass_fraction": coefficient_mass_fraction(mask, absolute_coefficients),
                    "positive_parent_fraction": float(np.mean(coefficients[mask] > 0.0)),
                    "mean_support_size": float(support_sizes[mask].mean()),
                }
            )

    return features, parent_type_rows


def normalize_case_id(value: str) -> str:
    text = Path(str(value).strip()).name
    if text.endswith(".tensors.npz"):
        text = text[: -len(".tensors.npz")]
    if text.endswith("_JW"):
        text = text[:-3]
    return text


def float_or_none(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        result = float(text)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def int_or_none(value: object) -> int | None:
    parsed = float_or_none(value)
    if parsed is None:
        return None
    rounded = int(round(parsed))
    if abs(parsed - rounded) > 1.0e-8:
        return None
    return rounded


def row_matches_settings(row: Mapping[str, str], steps: int, evolution_time: float) -> bool:
    status = str(row.get("status", "success")).strip().lower()
    if status and status not in {"success", "ok", "complete", "completed"}:
        return False

    step_fields = ("trotter_steps", "steps", "r")
    for field in step_fields:
        if field in row and str(row[field]).strip():
            row_steps = int_or_none(row[field])
            if row_steps is not None and row_steps != steps:
                return False
            break

    time_fields = ("evolution_time", "time", "total_time")
    for field in time_fields:
        if field in row and str(row[field]).strip():
            row_time = float_or_none(row[field])
            if row_time is not None and not math.isclose(
                row_time, evolution_time, rel_tol=1.0e-10, abs_tol=1.0e-12
            ):
                return False
            break
    return True


def outcome_file_priority(path: Path, explicit: bool) -> int:
    if explicit:
        return 100
    text = str(path).lower()
    if "fermionic_aware_performance" in text:
        return 90
    if "polyatomic_reference_deterministic" in text:
        return 85
    if "deterministic_orderings_hgbs5" in text:
        return 80
    if "deterministic" in text:
        return 70
    if "ablation" in text:
        return 20
    return 10


def discover_results_csvs(repo_root: Path, output_dir: Path) -> list[Path]:
    candidates: list[Path] = []
    preferred = (
        repo_root / "analysis/fermionic_aware_performance",
        repo_root / "analysis/polyatomic_reference_deterministic_ordering.csv",
        repo_root / "analysis/deterministic_orderings_hgbs5_results.csv",
        repo_root / "analysis/fermionic_structure_ablation_extension_hgbs5.csv",
    )
    for path in preferred:
        if path.is_file():
            candidates.append(path)
        elif path.is_dir():
            candidates.extend(sorted(path.glob("*.csv")))

    analysis_dir = repo_root / "analysis"
    if analysis_dir.is_dir():
        for path in analysis_dir.rglob("*.csv"):
            try:
                if output_dir.resolve() in path.resolve().parents:
                    continue
            except OSError:
                pass
            lowered = path.name.lower()
            if any(token in lowered for token in ("deterministic", "ordering", "performance", "ablation")):
                candidates.append(path)

    unique: dict[Path, None] = {}
    for path in candidates:
        if path.is_file():
            unique[path.resolve()] = None
    return list(unique)


def parse_manual_advantages(values: Sequence[str]) -> dict[str, float]:
    parsed: dict[str, float] = {}
    for item in values:
        if "=" not in item:
            raise ValueError(f"--manual-advantage must be CASE_ID=RATIO, got {item!r}")
        case_text, ratio_text = item.split("=", 1)
        case_id = normalize_case_id(case_text)
        ratio = float(ratio_text)
        if not math.isfinite(ratio) or ratio <= 0.0:
            raise ValueError(f"Invalid advantage ratio in {item!r}")
        parsed[case_id] = ratio
    return parsed


def load_outcomes(
    cases: Sequence[CaseSpec],
    repo_root: Path,
    explicit_csvs: Sequence[Path],
    steps: int,
    evolution_time: float,
    floor: float,
    tie_rtol: float,
    allow_incomplete_jw_baseline: bool,
    manual_advantages: Mapping[str, float],
    output_dir: Path,
) -> tuple[dict[str, Outcome], list[str]]:
    outcomes = {case.case_id: Outcome(case_id=case.case_id) for case in cases}
    case_ids = set(outcomes)
    warnings: list[str] = []

    candidate_files: list[tuple[Path, int]] = []
    if explicit_csvs:
        for path in explicit_csvs:
            resolved = path if path.is_absolute() else repo_root / path
            if not resolved.is_file():
                warnings.append(f"Outcome CSV not found: {resolved}")
                continue
            candidate_files.append((resolved, outcome_file_priority(resolved, explicit=True)))
    else:
        for path in discover_results_csvs(repo_root, output_dir):
            candidate_files.append((path, outcome_file_priority(path, explicit=False)))

    # Keep one value per (case, ordering-kind), favoring explicit/current
    # deterministic sources over lower-priority legacy/ablation files.
    observed: dict[tuple[str, str], tuple[int, float, str]] = {}
    direct_advantage: dict[str, tuple[int, float, str]] = {}

    for path, priority in candidate_files:
        try:
            with path.open("r", newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                if not reader.fieldnames:
                    continue
                fields = {field.strip() for field in reader.fieldnames if field}
                for row in reader:
                    if not row_matches_settings(row, steps, evolution_time):
                        continue
                    case_value = row.get("case_id") or row.get("tensor_path") or row.get("case")
                    if not case_value:
                        continue
                    case_id = normalize_case_id(case_value)
                    if case_id not in case_ids:
                        continue

                    # Support already-aggregated performance tables.
                    advantage_candidates = [
                        "fermionic_advantage",
                        "fermionic_advantage_best_jw_over_fermionic",
                        "advantage_vs_best_jw",
                        "best_jw_over_fermionic",
                        "advantage",
                    ]
                    advantage_field = next(
                        (field for field in advantage_candidates if field in fields),
                        None,
                    )
                    if advantage_field is None:
                        advantage_field = next(
                            (
                                field
                                for field in fields
                                if "advantage" in field.lower()
                                and "best" in field.lower()
                                and "jw" in field.lower()
                                and "median" not in field.lower()
                                and "geometric" not in field.lower()
                            ),
                            None,
                        )
                    if advantage_field:
                        advantage = float_or_none(row.get(advantage_field))
                        if advantage is not None and advantage > 0.0:
                            old = direct_advantage.get(case_id)
                            if old is None or priority > old[0]:
                                direct_advantage[case_id] = (priority, advantage, str(path))

                    # Support aggregated tables with explicit error columns.
                    aggregate_error_fields = {
                        "fermionic": (
                            "fermionic_error",
                            "fermionic_one_minus_overlap",
                            "fermionic_signed_error",
                        ),
                        "jw_signed": ("jw_signed_error", "jw_signed_one_minus_overlap"),
                        "jw_magnitude": (
                            "jw_magnitude_error",
                            "jw_magnitude_one_minus_overlap",
                            "jw_magnitude_descending_error",
                        ),
                        "best_jw": ("best_jw_error", "best_jw_one_minus_overlap"),
                    }
                    for kind, aliases in aggregate_error_fields.items():
                        for alias in aliases:
                            if alias in row:
                                value = float_or_none(row.get(alias))
                                if value is not None and value >= 0.0:
                                    old = observed.get((case_id, kind))
                                    if old is None or priority > old[0]:
                                        observed[(case_id, kind)] = (priority, value, str(path))
                                break

                    ordering = str(
                        row.get("ordering")
                        or row.get("schedule")
                        or row.get("ordering_name")
                        or ""
                    ).strip()
                    error = float_or_none(
                        row.get("one_minus_overlap")
                        or row.get("overlap_error")
                        or row.get("one_minus_abs_overlap")
                    )
                    if not ordering or error is None or error < 0.0:
                        continue
                    if ordering in FERMIONIC_ALIASES:
                        kind = "fermionic"
                    elif ordering in JW_SIGNED_ALIASES:
                        kind = "jw_signed"
                    elif ordering in JW_MAGNITUDE_ALIASES:
                        kind = "jw_magnitude"
                    else:
                        continue
                    old = observed.get((case_id, kind))
                    if old is None or priority > old[0]:
                        observed[(case_id, kind)] = (priority, error, str(path))
        except (OSError, UnicodeError, csv.Error) as exc:
            warnings.append(f"Could not read outcome CSV {path}: {exc}")

    for case in cases:
        outcome = outcomes[case.case_id]
        entries = {
            kind: observed.get((case.case_id, kind))
            for kind in ("fermionic", "jw_signed", "jw_magnitude", "best_jw")
        }
        if entries["fermionic"]:
            outcome.fermionic_error = entries["fermionic"][1]
        if entries["jw_signed"]:
            outcome.jw_signed_error = entries["jw_signed"][1]
        if entries["jw_magnitude"]:
            outcome.jw_magnitude_error = entries["jw_magnitude"][1]

        if entries["best_jw"]:
            outcome.best_jw_error = entries["best_jw"][1]
            outcome.best_jw_ordering = "precomputed_best_jw"
            outcome.baseline_complete = True
        elif outcome.jw_signed_error is not None and outcome.jw_magnitude_error is not None:
            if outcome.jw_signed_error <= outcome.jw_magnitude_error:
                outcome.best_jw_error = outcome.jw_signed_error
                outcome.best_jw_ordering = "jw_signed"
            else:
                outcome.best_jw_error = outcome.jw_magnitude_error
                outcome.best_jw_ordering = "jw_magnitude_descending"
            outcome.baseline_complete = True
        elif allow_incomplete_jw_baseline:
            available = [
                ("jw_signed", outcome.jw_signed_error),
                ("jw_magnitude_descending", outcome.jw_magnitude_error),
            ]
            available = [(name, value) for name, value in available if value is not None]
            if available:
                outcome.best_jw_ordering, outcome.best_jw_error = min(
                    available, key=lambda pair: pair[1]
                )
                outcome.baseline_complete = False

        sources = sorted(
            {
                entry[2]
                for entry in entries.values()
                if entry is not None
            }
        )
        outcome.source = ";".join(sources)

        direct = direct_advantage.get(case.case_id)
        if direct is not None:
            outcome.advantage = direct[1]
            outcome.log10_advantage = math.log10(direct[1])
            outcome.source = direct[2]
            # A case-level best-JW/fermionic advantage is already the fair
            # baseline result; its source is responsible for that aggregation.
            outcome.baseline_complete = True
        elif outcome.fermionic_error is not None and outcome.best_jw_error is not None:
            if outcome.fermionic_error <= floor or outcome.best_jw_error <= floor:
                outcome.classification = "numerical_floor"
            elif outcome.fermionic_error > 0.0:
                outcome.advantage = outcome.best_jw_error / outcome.fermionic_error
                outcome.log10_advantage = math.log10(outcome.advantage)

        if outcome.advantage is not None and outcome.classification != "numerical_floor":
            if math.isclose(outcome.advantage, 1.0, rel_tol=tie_rtol, abs_tol=0.0):
                outcome.classification = "tie"
            elif outcome.advantage > 1.0:
                outcome.classification = "win"
            else:
                outcome.classification = "loss"

    for case_id, ratio in manual_advantages.items():
        if case_id not in outcomes:
            warnings.append(f"Manual advantage ignored for unknown target case: {case_id}")
            continue
        outcome = outcomes[case_id]
        outcome.advantage = ratio
        outcome.log10_advantage = math.log10(ratio)
        outcome.classification = (
            "tie"
            if math.isclose(ratio, 1.0, rel_tol=tie_rtol, abs_tol=0.0)
            else ("win" if ratio > 1.0 else "loss")
        )
        outcome.source = "manual_override"
        outcome.baseline_complete = True

    return outcomes, warnings


def rankdata(values: Sequence[float]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(x.size, dtype=float)
    start = 0
    while start < x.size:
        end = start + 1
        while end < x.size and x[order[end]] == x[order[start]]:
            end += 1
        average_rank = 0.5 * ((start + 1) + end)
        ranks[order[start:end]] = average_rank
        start = end
    return ranks


def spearman_correlation(left: Sequence[float], right: Sequence[float]) -> float:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    if float(rx.std()) == 0.0 or float(ry.std()) == 0.0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def feature_category(feature: str) -> str:
    if feature.startswith("coefficient_") or feature.startswith("positive_"):
        return "coefficient"
    if feature.startswith("excitation_") or feature.startswith("body_") or feature.startswith("diagonal_") or feature.startswith("mixed_") or feature.startswith("support_"):
        return "parent_type"
    if feature.startswith("virtual_") or feature.startswith("occupied_") or feature.startswith("orbital_"):
        return "orbital_concentration"
    if feature.startswith("ov_coupling") or feature.startswith("parent_orbital"):
        return "network_rank"
    if feature.startswith("signed_order"):
        return "ordering_locality"
    return "parent_graph"


def numeric_feature_names(case_features: Sequence[Mapping[str, object]]) -> list[str]:
    excluded = {
        "case_id",
        "molecule",
        "tensor_path",
        "basis",
        "contrast_role",
        "active_occupied",
        "active_vacant",
        "n_qubits",
        "number_of_fermionic_parents",
        "number_of_normal_ordered_monomials",
        "expected_parent_count",
        "parent_count_delta",
        "number_of_components",
    }
    names = set(case_features[0])
    for row in case_features[1:]:
        names &= set(row)
    result: list[str] = []
    for name in sorted(names - excluded):
        values = []
        valid = True
        for row in case_features:
            try:
                value = float(row[name])
            except (TypeError, ValueError):
                valid = False
                break
            if not math.isfinite(value):
                valid = False
                break
            values.append(value)
        if valid and max(values) - min(values) > 1.0e-15:
            result.append(name)
    return result


def rank_features(
    case_features: Sequence[Mapping[str, object]],
    outcomes: Mapping[str, Outcome],
) -> list[FeatureRankingRow]:
    by_molecule = {str(row["molecule"]): row for row in case_features}
    required = {"H2O", "NH3", "BeH2", "O2"}
    if set(by_molecule) != required:
        raise ValueError(f"Expected exactly {sorted(required)}, got {sorted(by_molecule)}")

    molecule_order = ["H2O", "NH3", "BeH2", "O2"]
    advantages = []
    for molecule in molecule_order:
        case_id = str(by_molecule[molecule]["case_id"])
        outcome = outcomes.get(case_id)
        advantages.append(
            outcome.log10_advantage if outcome is not None and outcome.log10_advantage is not None else float("nan")
        )

    rows: list[FeatureRankingRow] = []
    for feature in numeric_feature_names(case_features):
        values = np.asarray([float(by_molecule[molecule][feature]) for molecule in molecule_order])
        h2o, nh3, beh2, o2 = map(float, values)
        favorable_values = np.asarray([h2o, nh3], dtype=float)
        favorable_mean = float(favorable_values.mean())
        favorable_std = float(favorable_values.std())
        full_range = float(values.max() - values.min())
        gap = o2 - favorable_mean
        range_normalized_gap = abs(gap) / full_range if full_range > 0.0 else 0.0
        favorable_consistency = (
            1.0 - abs(h2o - nh3) / full_range if full_range > 0.0 else 1.0
        )
        favorable_consistency = min(1.0, max(0.0, favorable_consistency))
        symmetric_denominator = 0.5 * (abs(o2) + abs(favorable_mean)) + 1.0e-15
        symmetric_relative_gap = abs(gap) / symmetric_denominator

        # Anchor score favors: (a) H2O/NH3 consistency, (b) large separation
        # from O2 relative to the four-case range, and (c) nontrivial relative
        # effect size. The final factor prevents a tiny 1-3% numerical shift
        # from outranking a physically much larger concentration difference.
        relative_effect_factor = math.sqrt(min(1.0, symmetric_relative_gap / 0.25))
        anchor_separation_score = (
            range_normalized_gap * favorable_consistency * relative_effect_factor
        )
        spearman = spearman_correlation(values, advantages)
        if math.isfinite(spearman):
            combined_score = 0.75 * anchor_separation_score + 0.25 * abs(spearman)
        else:
            combined_score = anchor_separation_score
        direction = "O2 higher" if gap > 0.0 else "H2O/NH3 higher"

        rows.append(
            FeatureRankingRow(
                feature=feature,
                description=FEATURE_DESCRIPTIONS.get(feature, feature.replace("_", " ")),
                category=feature_category(feature),
                h2o=h2o,
                nh3=nh3,
                beh2=beh2,
                o2=o2,
                favorable_mean=favorable_mean,
                favorable_std=favorable_std,
                o2_minus_favorable=gap,
                symmetric_relative_gap=symmetric_relative_gap,
                range_normalized_gap=range_normalized_gap,
                favorable_consistency=favorable_consistency,
                anchor_separation_score=anchor_separation_score,
                spearman_vs_log10_advantage=spearman,
                combined_score=combined_score,
                direction=direction,
            )
        )

    rows.sort(key=lambda row: (row.combined_score, row.anchor_separation_score), reverse=True)
    return rows


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ranking_to_dict(row: FeatureRankingRow) -> dict[str, object]:
    return {
        "feature": row.feature,
        "description": row.description,
        "category": row.category,
        "H2O": row.h2o,
        "NH3": row.nh3,
        "BeH2": row.beh2,
        "O2": row.o2,
        "H2O_NH3_mean": row.favorable_mean,
        "H2O_NH3_std": row.favorable_std,
        "O2_minus_H2O_NH3_mean": row.o2_minus_favorable,
        "symmetric_relative_gap": row.symmetric_relative_gap,
        "range_normalized_gap": row.range_normalized_gap,
        "H2O_NH3_consistency": row.favorable_consistency,
        "anchor_separation_score": row.anchor_separation_score,
        "spearman_vs_log10_advantage": row.spearman_vs_log10_advantage,
        "combined_score": row.combined_score,
        "direction": row.direction,
    }


def outcome_columns(outcome: Outcome) -> dict[str, object]:
    return {
        "fermionic_one_minus_overlap": outcome.fermionic_error if outcome.fermionic_error is not None else "",
        "jw_signed_one_minus_overlap": outcome.jw_signed_error if outcome.jw_signed_error is not None else "",
        "jw_magnitude_one_minus_overlap": outcome.jw_magnitude_error if outcome.jw_magnitude_error is not None else "",
        "best_jw_one_minus_overlap": outcome.best_jw_error if outcome.best_jw_error is not None else "",
        "best_jw_ordering": outcome.best_jw_ordering,
        "fermionic_advantage_best_jw_over_fermionic": outcome.advantage if outcome.advantage is not None else "",
        "log10_fermionic_advantage": outcome.log10_advantage if outcome.log10_advantage is not None else "",
        "outcome_classification": outcome.classification,
        "best_jw_baseline_complete": int(outcome.baseline_complete),
        "outcome_source": outcome.source,
    }


def choose_interpretable_top_features(
    ranking: Sequence[FeatureRankingRow],
    count: int,
) -> list[FeatureRankingRow]:
    preferred = [row for row in ranking if row.feature in REPORT_FEATURE_WHITELIST]
    # Avoid filling the report with several near-duplicate variants of the same
    # signal. Keep at most two features from one category until every important
    # category has a chance to appear.
    selected: list[FeatureRankingRow] = []
    per_category: dict[str, int] = {}
    for row in preferred:
        if per_category.get(row.category, 0) >= 3:
            continue
        selected.append(row)
        per_category[row.category] = per_category.get(row.category, 0) + 1
        if len(selected) >= count:
            break
    return selected


def format_value(value: float) -> str:
    magnitude = abs(value)
    if magnitude != 0.0 and (magnitude < 1.0e-3 or magnitude >= 1.0e3):
        return f"{value:.3e}"
    return f"{value:.4f}"


def write_report(
    path: Path,
    case_features: Sequence[Mapping[str, object]],
    outcomes: Mapping[str, Outcome],
    ranking: Sequence[FeatureRankingRow],
    warnings: Sequence[str],
    top_count: int,
) -> None:
    by_molecule = {str(row["molecule"]): row for row in case_features}
    top = choose_interpretable_top_features(ranking, top_count)
    outcome_available = sum(
        1
        for row in case_features
        if outcomes[str(row["case_id"])].advantage is not None
    )

    lines: list[str] = []
    lines.append("# Pre-BCH Structural Features for Fermionic-Aware Ordering")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(
        "All predictor features below are computed **before** BCH evaluation and use only "
        "fermionic parent coefficients, orbital indices, occupied/virtual excitation structure, "
        "and parent/orbital graph structure. No BCH norm, commutator-vector norm, or BCH "
        "cancellation quantity is used as a feature."
    )
    lines.append("")
    lines.append(
        "When an outcome CSV is available, the external target is `one_minus_overlap` at "
        "T=1 and r=100, with the fair JW baseline `min(JW signed, JW magnitude-descending)`."
    )
    lines.append("")
    lines.append("## Four-case outcome check")
    lines.append("")
    lines.append("| Molecule | active occ+virt | parents | best-JW / fermionic | label |")
    lines.append("|---|---:|---:|---:|---|")
    for molecule in ("H2O", "NH3", "BeH2", "O2"):
        row = by_molecule[molecule]
        outcome = outcomes[str(row["case_id"])]
        advantage = "n/a" if outcome.advantage is None else f"{outcome.advantage:.4g}"
        lines.append(
            f"| {molecule} | {row['active_occupied']}+{row['active_vacant']} | "
            f"{row['number_of_fermionic_parents']} | {advantage} | {outcome.classification} |"
        )
    lines.append("")
    if outcome_available < 4:
        lines.append(
            f"> Outcome correlation is partial: fair best-JW/fermionic advantages were found for "
            f"{outcome_available}/4 cases. The H2O/NH3-vs-O2 structural separation ranking still "
            "runs independently of BCH and independently of missing outcome rows."
        )
        lines.append("")

    lines.append("## Best pre-BCH separators")
    lines.append("")
    lines.append(
        "The ranking favors features for which H2O and NH3 are mutually consistent but O2 is "
        "well separated; when all four external advantages are available, a small Spearman "
        "correlation term is added. With only four cases these are descriptive diagnostics, not "
        "a fitted predictive model."
    )
    lines.append("")
    lines.append("| feature | H2O | NH3 | BeH2 | O2 | direction | score |")
    lines.append("|---|---:|---:|---:|---:|---|---:|")
    for row in top:
        lines.append(
            f"| `{row.feature}` | {format_value(row.h2o)} | {format_value(row.nh3)} | "
            f"{format_value(row.beh2)} | {format_value(row.o2)} | {row.direction} | "
            f"{row.combined_score:.3f} |"
        )
    lines.append("")

    # Pull the mechanistically clearest concentration/rank diagnostics even if
    # one falls just outside the numerical top-N ranking.
    ranking_by_name = {row.feature: row for row in ranking}
    key_names = [
        "virtual_orbital_coefficient_mass_entropy",
        "virtual_orbital_coefficient_mass_gini",
        "virtual_orbital_coefficient_mass_top1_fraction",
        "ov_coupling_top_singular_fraction",
        "ov_coupling_effective_rank_fraction",
        "excitation_rank2_coefficient_mass_fraction",
    ]
    keys = [ranking_by_name[name] for name in key_names if name in ranking_by_name]

    lines.append("## Structural interpretation")
    lines.append("")
    if keys:
        virt_entropy = ranking_by_name.get("virtual_orbital_coefficient_mass_entropy")
        virt_gini = ranking_by_name.get("virtual_orbital_coefficient_mass_gini")
        virt_top = ranking_by_name.get("virtual_orbital_coefficient_mass_top1_fraction")
        ov_top = ranking_by_name.get("ov_coupling_top_singular_fraction")
        ov_eff = ranking_by_name.get("ov_coupling_effective_rank_fraction")
        exc2 = ranking_by_name.get("excitation_rank2_coefficient_mass_fraction")

        if virt_entropy and virt_gini and virt_top:
            lines.append(
                "1. **O2 is much more concentrated on a small set of virtual-orbital hubs.** "
                f"Its virtual coefficient-mass entropy is {virt_entropy.o2:.3f}, versus "
                f"{virt_entropy.h2o:.3f}/{virt_entropy.nh3:.3f} for H2O/NH3; its virtual "
                f"Gini is {virt_gini.o2:.3f}, versus {virt_gini.h2o:.3f}/{virt_gini.nh3:.3f}; "
                f"and one virtual orbital carries {100.0*virt_top.o2:.1f}% of the weighted "
                f"incidence, versus {100.0*virt_top.h2o:.1f}%/{100.0*virt_top.nh3:.1f}%."
            )
            lines.append("")
        if ov_top and ov_eff:
            lines.append(
                "2. **O2's occupied-to-virtual parent network is lower-rank/more hub dominated.** "
                f"The top singular channel carries {100.0*ov_top.o2:.1f}% for O2 versus "
                f"{100.0*ov_top.h2o:.1f}%/{100.0*ov_top.nh3:.1f}% for H2O/NH3. The normalized "
                f"effective rank is {ov_eff.o2:.3f}, versus {ov_eff.h2o:.3f}/{ov_eff.nh3:.3f}."
            )
            lines.append("")
        if exc2:
            lines.append(
                "3. **O2 also puts more coefficient mass into rank-2 occupied/virtual transfer "
                "parents.** The fraction is "
                f"{100.0*exc2.o2:.3f}% for O2, compared with "
                f"{100.0*exc2.h2o:.3f}%/{100.0*exc2.nh3:.3f}% for H2O/NH3."
            )
            lines.append("")

    lines.append(
        "**Working pre-BCH hypothesis:** signed fermionic parent ordering is more likely to help "
        "when significant parent coefficient weight is distributed across many virtual orbitals "
        "and many occupied-to-virtual structural channels, rather than being concentrated in a "
        "small, low-rank set of hubs. A distributed network gives the parent sequence more "
        "structural degrees of freedom to change downstream interference; an O2-like hub-dominated "
        "network provides fewer such degrees of freedom."
    )
    lines.append("")
    lines.append(
        "BeH2 is useful as a check because it is mixed: its virtual-orbital coefficient mass is "
        "fairly distributed like H2O/NH3, while its occupied-to-virtual coupling is more "
        "concentrated. That is exactly the kind of intermediate structural signature expected for "
        "a case whose fermionic advantage is weaker or less robust."
    )
    lines.append("")
    lines.append("## Caveats / next test")
    lines.append("")
    lines.append(
        "This four-case analysis is hypothesis generation, not a proof. The next validation should "
        "apply these same pre-BCH features to the larger HGBS-5 sweep (especially Be2, Li2, F2, "
        "H2O, NH3, CH4, N2, and O2) and test whether virtual-hub concentration and O-V effective "
        "rank predict the externally measured best-JW/fermionic advantage out of sample."
    )

    if warnings:
        lines.append("")
        lines.append("## Warnings")
        lines.append("")
        for warning in warnings:
            lines.append(f"- {warning}")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_plot(
    path: Path,
    ranking: Sequence[FeatureRankingRow],
    top_count: int,
) -> str | None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return "matplotlib not installed; skipped pre_bch_top_features.png"

    top = choose_interpretable_top_features(ranking, min(top_count, 10))
    if not top:
        return "no ranked features available; skipped pre_bch_top_features.png"

    molecules = ["H2O", "NH3", "BeH2", "O2"]
    matrix = np.asarray([[row.h2o, row.nh3, row.beh2, row.o2] for row in top], dtype=float)
    # Normalize each feature independently so heterogeneous structural metrics
    # can share a diagnostic heatmap. This normalization is only for the plot.
    centers = matrix.mean(axis=1, keepdims=True)
    scales = matrix.std(axis=1, keepdims=True)
    scales[scales == 0.0] = 1.0
    z = (matrix - centers) / scales

    fig_height = max(4.8, 0.52 * len(top) + 1.8)
    fig, ax = plt.subplots(figsize=(8.8, fig_height))
    image = ax.imshow(z, aspect="auto")
    ax.set_xticks(np.arange(len(molecules)))
    ax.set_xticklabels(molecules)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels([row.feature for row in top], fontsize=8)
    ax.set_title("Pre-BCH structural separators (row-wise z score)")
    ax.set_xlabel("HGBS-5 case")
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("within-feature z score")
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return None


def run_self_test() -> None:
    # Canonicalization/parity checks.
    key, sign = canonicalize_number_conserving_key(((0, 1), (1, 1), (2, 0), (3, 0)))
    assert key == ((1, 1), (0, 1), (3, 0), (2, 0))
    assert sign == 1
    key2, sign2 = canonicalize_number_conserving_key(((1, 1), (0, 1), (2, 0), (3, 0)))
    assert key2 == ((1, 1), (0, 1), (3, 0), (2, 0))
    assert sign2 == -1
    vanished, vanished_sign = canonicalize_number_conserving_key(
        ((1, 1), (1, 1), (2, 0), (3, 0))
    )
    assert vanished is None and vanished_sign == 0

    # Distribution metrics.
    assert math.isclose(gini([1, 1, 1, 1]), 0.0, abs_tol=1.0e-15)
    assert math.isclose(normalized_entropy([1, 1, 1, 1]), 1.0, abs_tol=1.0e-15)
    assert math.isclose(effective_count([1, 1, 1, 1]), 4.0, abs_tol=1.0e-15)
    assert math.isclose(spearman_correlation([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
    assert math.isclose(spearman_correlation([1, 2, 3, 4], [40, 30, 20, 10]), -1.0)
    print("SELF-TEST PASS: pre-BCH canonicalization and structural metrics verified.")


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.steps <= 0:
        raise ValueError("--steps must be positive")
    if args.evolution_time <= 0.0:
        raise ValueError("--time must be positive")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive")

    repo_root = args.tensor_root.resolve()
    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        CaseSpec(
            molecule=molecule,
            tensor_path=relative_path,
            active_occupied=active_occupied,
            active_vacant=active_vacant,
            expected_parent_count=expected_parent_count,
            contrast_role=contrast_role,
        )
        for (
            molecule,
            relative_path,
            active_occupied,
            active_vacant,
            expected_parent_count,
            contrast_role,
        ) in DEFAULT_CASES
    ]

    missing = [case for case in cases if not (repo_root / case.tensor_path).is_file()]
    if missing:
        details = "\n".join(f"  - {repo_root / case.tensor_path}" for case in missing)
        raise FileNotFoundError(
            "Missing required HGBS-5 tensor cases. Run from the QHAT repository root "
            f"or pass --tensor-root. Missing:\n{details}"
        )

    case_features: list[dict[str, float | int | str]] = []
    parent_type_rows: list[dict[str, float | int | str]] = []
    warnings: list[str] = []

    print("Pre-BCH structural analysis")
    print("=" * 88)
    print("Predictors: coefficients + orbital/excitation/connectivity/network structure only")
    print("BCH/commutator features: DISABLED")
    print()

    for case in cases:
        tensor_path = repo_root / case.tensor_path
        print(f"[{case.molecule}] {case.case_id}", flush=True)
        features, type_rows = compute_case_features(case, tensor_path, args.tolerance)
        case_features.append(features)
        parent_type_rows.extend(type_rows)
        parent_count = int(features["number_of_fermionic_parents"])
        delta = features["parent_count_delta"]
        print(
            f"  parents={parent_count}  monomials={features['number_of_normal_ordered_monomials']}  "
            f"virtual-entropy={float(features['virtual_orbital_coefficient_mass_entropy']):.6f}  "
            f"OV-effective-rank-frac={float(features['ov_coupling_effective_rank_fraction']):.6f}"
        )
        if isinstance(delta, int) and delta != 0:
            warning = (
                f"{case.case_id}: reconstructed parent count {parent_count} differs from the "
                f"reference count {case.expected_parent_count} by {delta:+d}. This can indicate "
                "a tensor revision; structural features use the supplied tensor exactly."
            )
            warnings.append(warning)
            print(f"  WARNING: {warning}")

    manual_advantages = parse_manual_advantages(args.manual_advantage)
    outcomes, outcome_warnings = load_outcomes(
        cases=cases,
        repo_root=repo_root,
        explicit_csvs=args.results_csv,
        steps=args.steps,
        evolution_time=args.evolution_time,
        floor=args.floor,
        tie_rtol=args.tie_rtol,
        allow_incomplete_jw_baseline=args.allow_incomplete_jw_baseline,
        manual_advantages=manual_advantages,
        output_dir=output_dir,
    )
    warnings.extend(outcome_warnings)

    # Attach external outcomes after all pre-BCH features are finalized.
    combined_case_rows: list[dict[str, object]] = []
    for feature_row in case_features:
        row = dict(feature_row)
        row.update(outcome_columns(outcomes[str(row["case_id"])]))
        combined_case_rows.append(row)

    ranking = rank_features(case_features, outcomes)
    ranking_rows = [ranking_to_dict(row) for row in ranking]

    case_csv = output_dir / "pre_bch_case_features.csv"
    ranking_csv = output_dir / "pre_bch_feature_ranking.csv"
    type_csv = output_dir / "pre_bch_parent_type_summary.csv"
    report_md = output_dir / "pre_bch_report.md"
    plot_png = output_dir / "pre_bch_top_features.png"

    case_fieldnames: list[str] = []
    for row in combined_case_rows:
        for key in row:
            if key not in case_fieldnames:
                case_fieldnames.append(key)
    write_csv(case_csv, combined_case_rows, case_fieldnames)

    ranking_fieldnames = list(ranking_rows[0]) if ranking_rows else []
    if ranking_rows:
        write_csv(ranking_csv, ranking_rows, ranking_fieldnames)

    type_fieldnames = list(parent_type_rows[0]) if parent_type_rows else []
    if parent_type_rows:
        write_csv(type_csv, parent_type_rows, type_fieldnames)

    write_report(
        report_md,
        case_features,
        outcomes,
        ranking,
        warnings,
        args.top_features,
    )
    plot_warning = make_plot(plot_png, ranking, args.top_features)
    if plot_warning:
        warnings.append(plot_warning)

    print()
    print("External outcome labels (best JW / fermionic):")
    for case in cases:
        outcome = outcomes[case.case_id]
        advantage = "n/a" if outcome.advantage is None else f"{outcome.advantage:.6g}"
        completeness = "complete-best-JW" if outcome.baseline_complete else "incomplete/unavailable"
        print(
            f"  {case.molecule:5s} advantage={advantage:>10s}  "
            f"label={outcome.classification:>15s}  baseline={completeness}"
        )

    print()
    print("Top interpretable pre-BCH separators (H2O/NH3 vs O2):")
    for rank, row in enumerate(choose_interpretable_top_features(ranking, args.top_features), start=1):
        correlation = (
            "n/a"
            if not math.isfinite(row.spearman_vs_log10_advantage)
            else f"{row.spearman_vs_log10_advantage:+.3f}"
        )
        print(
            f"  {rank:2d}. {row.feature:52s} score={row.combined_score:.3f}  "
            f"{row.direction:14s}  rho={correlation}"
        )

    print()
    print("Wrote:")
    for path in (case_csv, ranking_csv, type_csv, report_md):
        print(f"  {path}")
    if plot_png.is_file():
        print(f"  {plot_png}")
    if warnings:
        print()
        print("Warnings:")
        for warning in warnings:
            print(f"  - {warning}")


if __name__ == "__main__":
    main()
