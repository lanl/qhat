#!/usr/bin/env python3
"""Test whether weighted fermionic-graph centralization predicts ordering gain.

The analysis uses active spin orbitals as graph vertices, partitioned by the
Hartree--Fock occupied/virtual split.  Each complete Hermitian fermionic parent
contributes its absolute canonical coefficient to every occupied--virtual pair
it jointly touches.  This is the graph version of the occupied-to-virtual
incidence network that motivated the qualitative hub-dominance hypothesis.

The primary graph statistic applies strength-weighted Freeman degree
centralization separately to the two sides of this bipartite graph and uses the
larger value.  Partition-wise normalization prevents active-space imbalance
from being mistaken for a hub.  The response is the requested
JW-magnitude/fermionic one-minus-overlap ratio; values above one favor the
fermionic-induced ordering.
"""

# ruff: noqa: E402  # Select a writable Matplotlib cache before importing it.

from __future__ import annotations

import argparse
import gc
import math
import os
import pickle
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

_MPL_CACHE = Path(tempfile.gettempdir()) / "qhat-matplotlib-cache"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openfermion import InteractionOperator, get_fermion_operator
from scipy import stats

try:
    from qhat.analysis.benchmark_L_sweep_trotter import (
        HermitianFermionTerm,
        build_hermitian_fermion_terms,
        clean_fermion_operator,
        load_interaction_operator,
    )
except ImportError:
    from benchmark_L_sweep_trotter import (
        HermitianFermionTerm,
        build_hermitian_fermion_terms,
        clean_fermion_operator,
        load_interaction_operator,
    )


DEFAULT_INPUT = Path(
    "analysis/fermionic_aware_performance/fermionic_aware_case_performance.csv"
)
DEFAULT_OUTDIR = Path("analysis/fermionic_graph_centralization")
DEFAULT_TOLERANCE = 1.0e-12
PRIMARY_METRIC = "bipartite_weighted_freeman_centralization"

MATCHED_MOLECULES = ("B-H", "Be-H", "Li-H")

REQUIRED_COLUMNS = {
    "case_id",
    "tensor_path",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "numerical_error_floor",
    "jw_magnitude_one_minus_overlap",
    "fermionic_aware_one_minus_overlap",
    "valid_comparison",
}


@dataclass
class MolecularDataCache:
    """Keep at most one large full-space molecular pickle resident."""

    path: Path | None = None
    molecule: Any | None = None

    def load(self, path: Path) -> Any:
        if self.path != path:
            self.molecule = None
            gc.collect()
            with path.open("rb") as stream:
                self.molecule = pickle.load(stream)  # noqa: S301 - trusted repo data
            self.path = path
        return self.molecule


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--basis", default="hgbs-5")
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "pdf", "svg"),
        default=["png", "pdf"],
    )
    return parser.parse_args()


def hermitian_term_mass(
    term: HermitianFermionTerm,
    tolerance: float,
) -> float:
    """Return the absolute canonical coefficient of a Hermitian parent.

    The deterministic fermionic ordering uses the lexicographically smallest
    component as the parent representative.  A non-self-adjoint parent and its
    adjoint have equal coefficient magnitude, so this avoids doubling only
    those parents that happen to require two normal-ordered components.
    """
    canonical_key = min(term.component_keys)
    coefficient = complex(term.operator.terms[canonical_key])
    return float(abs(coefficient)) if abs(coefficient) > tolerance else 0.0


def hermitian_term_support(term: HermitianFermionTerm) -> tuple[int, ...]:
    """Return the sorted spin-orbital support of a Hermitian term."""
    return tuple(
        sorted(
            {
                orbital
                for component in term.component_keys
                for orbital, _ in component
            }
        )
    )


def project_terms_to_orbital_weights(
    terms: Sequence[HermitianFermionTerm],
    n_qubits: int,
    active_occupied: int,
    tolerance: float,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Project fermionic parents onto occupied-to-virtual orbital edges.

    A parent of absolute canonical coefficient ``a`` contributes ``a`` to
    every occupied/virtual support pair.  Parents confined to one partition do
    not create an edge in this graph.  This matches the coefficient-weighted
    occupied-to-virtual incidence construction used by the motivating pre-BCH
    structural analysis.
    """
    if active_occupied <= 0 or active_occupied >= n_qubits:
        raise ValueError(
            "active_occupied must leave nonempty occupied and virtual partitions."
        )
    weights = np.zeros((n_qubits, n_qubits), dtype=float)
    parent_coefficient_mass = 0.0
    mixed_parent_mass = 0.0
    noncrossing_parent_mass = 0.0
    noncrossing_terms = 0
    projected_terms = 0
    maximum_support = 0

    for term in terms:
        mass = hermitian_term_mass(term, tolerance)
        if mass <= tolerance:
            continue
        support = hermitian_term_support(term)
        maximum_support = max(maximum_support, len(support))
        if any(index < 0 or index >= n_qubits for index in support):
            raise ValueError(
                f"Fermionic support {support} exceeds {n_qubits} orbitals."
            )
        parent_coefficient_mass += mass
        occupied_support = [index for index in support if index < active_occupied]
        virtual_support = [index for index in support if index >= active_occupied]
        if not occupied_support or not virtual_support:
            noncrossing_parent_mass += mass
            noncrossing_terms += 1
            continue

        for left in occupied_support:
            for right in virtual_support:
                weights[left, right] += mass
                weights[right, left] += mass
        mixed_parent_mass += mass
        projected_terms += 1

    metadata: dict[str, float | int] = {
        "number_of_hermitian_fermionic_terms": len(terms),
        "number_of_mixed_occupied_virtual_parents": projected_terms,
        "number_of_noncrossing_parents": noncrossing_terms,
        "total_parent_coefficient_mass": parent_coefficient_mass,
        "mixed_parent_coefficient_mass": mixed_parent_mass,
        "noncrossing_parent_coefficient_mass": noncrossing_parent_mass,
        "maximum_fermionic_support_size": maximum_support,
    }
    return weights, metadata


def strength_gini(strengths: np.ndarray) -> float:
    total = float(strengths.sum())
    n_nodes = len(strengths)
    if n_nodes == 0 or total <= 0.0:
        return math.nan
    differences = np.abs(strengths[:, None] - strengths[None, :]).sum()
    return float(differences / (2.0 * n_nodes * total))


def partition_weighted_freeman_centralization(
    strengths: np.ndarray,
    total_weight: float,
) -> float:
    """Freeman strength centralization within one bipartite partition.

    The maximum at fixed incident weight is a partition star: one node carries
    all weight and every other node carries none.  Its Freeman numerator is
    ``(m - 1) W``.
    """
    n_nodes = len(strengths)
    if n_nodes <= 1 or total_weight <= 0.0:
        return math.nan
    numerator = float(np.sum(float(strengths.max()) - strengths))
    return float(np.clip(numerator / ((n_nodes - 1) * total_weight), 0.0, 1.0))


def partition_unweighted_freeman_centralization(
    degrees: np.ndarray,
    opposite_nodes: int,
) -> float:
    n_nodes = len(degrees)
    if n_nodes <= 1 or opposite_nodes <= 0:
        return math.nan
    numerator = float(np.sum(float(degrees.max()) - degrees))
    return float(
        np.clip(numerator / ((n_nodes - 1) * opposite_nodes), 0.0, 1.0)
    )


def orbital_graph_metrics(
    weights: np.ndarray,
    active_occupied: int,
) -> dict[str, float | int]:
    """Compute bipartite weighted Freeman centralization metrics."""
    if weights.ndim != 2 or weights.shape[0] != weights.shape[1]:
        raise ValueError("weights must be a square matrix.")
    if not np.allclose(weights, weights.T, rtol=0.0, atol=1.0e-14):
        raise ValueError("weights must be symmetric.")
    if np.any(weights < -1.0e-14):
        raise ValueError("weights must be nonnegative.")
    if not np.allclose(np.diag(weights), 0.0, rtol=0.0, atol=1.0e-14):
        raise ValueError("weights cannot contain self-loops.")

    n_nodes = weights.shape[0]
    if active_occupied <= 0 or active_occupied >= n_nodes:
        raise ValueError("active_occupied must define two nonempty partitions.")
    occupied_count = active_occupied
    virtual_count = n_nodes - active_occupied
    within_occupied = weights[:active_occupied, :active_occupied]
    within_virtual = weights[active_occupied:, active_occupied:]
    if np.any(within_occupied > 1.0e-14) or np.any(within_virtual > 1.0e-14):
        raise ValueError("occupied-to-virtual graph cannot contain within-part edges.")

    upper = weights[np.triu_indices(n_nodes, k=1)]
    total_weight = float(upper.sum())
    strengths = weights.sum(axis=1)
    max_strength = float(strengths.max(initial=0.0))
    strength_numerator = float(np.sum(max_strength - strengths))

    whole_graph_centralization = math.nan
    if n_nodes > 2 and total_weight > 0.0:
        whole_graph_centralization = strength_numerator / (
            (n_nodes - 2) * total_weight
        )
        whole_graph_centralization = float(
            np.clip(whole_graph_centralization, 0.0, 1.0)
        )

    occupied_strengths = strengths[:active_occupied]
    virtual_strengths = strengths[active_occupied:]
    occupied_centralization = partition_weighted_freeman_centralization(
        occupied_strengths,
        total_weight,
    )
    virtual_centralization = partition_weighted_freeman_centralization(
        virtual_strengths,
        total_weight,
    )
    finite_partition_values = [
        value
        for value in (occupied_centralization, virtual_centralization)
        if math.isfinite(value)
    ]
    bipartite_centralization = (
        max(finite_partition_values) if finite_partition_values else math.nan
    )
    mean_partition_centralization = (
        float(np.mean(finite_partition_values))
        if finite_partition_values
        else math.nan
    )

    adjacency = weights > 0.0
    degrees = adjacency.sum(axis=1).astype(float)
    edge_count = int(np.triu(adjacency, k=1).sum())
    occupied_degree_centralization = partition_unweighted_freeman_centralization(
        degrees[:active_occupied],
        virtual_count,
    )
    virtual_degree_centralization = partition_unweighted_freeman_centralization(
        degrees[active_occupied:],
        occupied_count,
    )

    total_strength = float(strengths.sum())
    effective_orbitals = math.nan
    squared_strength = float(np.dot(strengths, strengths))
    if squared_strength > 0.0:
        effective_orbitals = total_strength**2 / squared_strength

    density = (
        2.0 * edge_count / (n_nodes * (n_nodes - 1))
        if n_nodes > 1
        else math.nan
    )
    occupied_hub_fraction = (
        float(occupied_strengths.max(initial=0.0)) / total_weight
        if total_weight > 0.0
        else math.nan
    )
    virtual_hub_fraction = (
        float(virtual_strengths.max(initial=0.0)) / total_weight
        if total_weight > 0.0
        else math.nan
    )

    return {
        "graph_nodes": n_nodes,
        "occupied_graph_nodes": occupied_count,
        "virtual_graph_nodes": virtual_count,
        "graph_edges": edge_count,
        "graph_density": density,
        "graph_total_edge_weight": total_weight,
        "maximum_orbital_strength": max_strength,
        "occupied_hub_incident_weight_fraction": occupied_hub_fraction,
        "virtual_hub_incident_weight_fraction": virtual_hub_fraction,
        "effective_strength_orbitals": effective_orbitals,
        "strength_gini": strength_gini(strengths),
        "bipartite_weighted_freeman_centralization": bipartite_centralization,
        "mean_partition_weighted_freeman_centralization": (
            mean_partition_centralization
        ),
        "occupied_weighted_freeman_centralization": occupied_centralization,
        "virtual_weighted_freeman_centralization": virtual_centralization,
        "whole_graph_same_weight_star_centralization": whole_graph_centralization,
        "occupied_unweighted_degree_centralization": (
            occupied_degree_centralization
        ),
        "virtual_unweighted_degree_centralization": virtual_degree_centralization,
        "hub_orbital_index": int(np.argmax(strengths)) if n_nodes else -1,
    }


def build_orbital_graph_metrics(
    interaction: InteractionOperator,
    n_qubits: int,
    active_occupied: int,
    tolerance: float,
) -> dict[str, float | int]:
    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        tolerance,
    )
    terms = build_hermitian_fermion_terms(fermion_hamiltonian, tolerance)
    weights, projection = project_terms_to_orbital_weights(
        terms,
        n_qubits,
        active_occupied,
        tolerance,
    )
    return {
        **projection,
        **orbital_graph_metrics(weights, active_occupied),
    }


def parent_molecular_pickle(tensor_path: Path) -> Path:
    prefix = tensor_path.name.split("_as-", maxsplit=1)[0]
    return tensor_path.parent / f"{prefix}.pickle"


def active_interaction_from_molecular_data(
    molecular_data: Any,
    active_occupied: int,
    active_vacant: int,
) -> InteractionOperator:
    """Reconstruct QHAT's active-space InteractionOperator in memory."""
    frozen_occupied_spin = int(molecular_data.n_electrons) - active_occupied
    n_active_spin = active_occupied + active_vacant
    if frozen_occupied_spin < 0 or frozen_occupied_spin % 2:
        raise ValueError("Invalid active occupied count for molecular data.")
    if n_active_spin <= 0 or n_active_spin % 2:
        raise ValueError("Active spin-orbital count must be positive and even.")

    first_active_spatial = frozen_occupied_spin // 2
    n_active_spatial = n_active_spin // 2
    interaction = molecular_data.get_molecular_hamiltonian(
        occupied_indices=range(first_active_spatial),
        active_indices=range(
            first_active_spatial,
            first_active_spatial + n_active_spatial,
        ),
    )
    if interaction.one_body_tensor.shape != (n_active_spin, n_active_spin):
        raise ValueError(
            "Reconstructed interaction has the wrong active-space shape: "
            f"{interaction.one_body_tensor.shape}, expected "
            f"({n_active_spin}, {n_active_spin})."
        )
    return interaction


def load_case_interaction(
    row: pd.Series,
    cache: MolecularDataCache,
) -> tuple[InteractionOperator | None, str, str]:
    tensor_path = Path(str(row["tensor_path"]))
    if tensor_path.exists():
        interaction, n_qubits = load_interaction_operator(tensor_path)
        if n_qubits != int(row["n_qubits"]):
            raise ValueError(
                f"Tensor size mismatch for {row['case_id']}: "
                f"{n_qubits} != {row['n_qubits']}."
            )
        return interaction, "active_tensor", ""

    molecular_pickle = parent_molecular_pickle(tensor_path)
    if molecular_pickle.exists():
        molecular_data = cache.load(molecular_pickle)
        interaction = active_interaction_from_molecular_data(
            molecular_data,
            int(row["active_occupied"]),
            int(row["active_vacant"]),
        )
        return interaction, "reconstructed_from_full_molecular_pickle", ""

    reason = (
        f"missing active tensor ({tensor_path}) and full molecular pickle "
        f"({molecular_pickle})"
    )
    return None, "unavailable", reason


def read_cases(path: Path, basis: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Input is missing columns: {', '.join(missing)}")
    selected = frame[frame["basis"].str.lower() == basis.lower()].copy()
    if selected.empty:
        raise ValueError(f"No {basis} rows found in {path}.")
    return selected.sort_values(
        ["molecule", "n_qubits", "active_occupied", "active_vacant"]
    ).reset_index(drop=True)


def analyze_cases(
    cases: pd.DataFrame,
    tolerance: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    cache = MolecularDataCache()
    for _, source in cases.iterrows():
        base = {
            "case_id": source["case_id"],
            "tensor_path": source["tensor_path"],
            "molecule": source["molecule"],
            "bond_length": source["bond_length"],
            "basis": source["basis"],
            "active_occupied": int(source["active_occupied"]),
            "active_vacant": int(source["active_vacant"]),
            "n_qubits": int(source["n_qubits"]),
            "numerical_error_floor": float(source["numerical_error_floor"]),
            "jw_magnitude_one_minus_overlap": source[
                "jw_magnitude_one_minus_overlap"
            ],
            "fermionic_aware_one_minus_overlap": source[
                "fermionic_aware_one_minus_overlap"
            ],
            "input_valid_comparison": bool(source["valid_comparison"]),
        }
        print(f"{source['case_id']}: building coefficient-weighted graph")
        try:
            interaction, graph_source, unavailable_reason = load_case_interaction(
                source,
                cache,
            )
            if interaction is None:
                rows.append(
                    {
                        **base,
                        "graph_status": "unavailable",
                        "graph_source": graph_source,
                        "graph_unavailable_reason": unavailable_reason,
                        "ratio_valid": False,
                    }
                )
                print("  unavailable")
                continue

            metrics = build_orbital_graph_metrics(
                interaction,
                int(source["n_qubits"]),
                int(source["active_occupied"]),
                tolerance,
            )
            jw_error = float(source["jw_magnitude_one_minus_overlap"])
            fermionic_error = float(source["fermionic_aware_one_minus_overlap"])
            error_floor = float(source["numerical_error_floor"])
            ratio_valid = bool(source["valid_comparison"]) and (
                math.isfinite(jw_error)
                and math.isfinite(fermionic_error)
                and jw_error > error_floor
                and fermionic_error > error_floor
            )
            ratio = jw_error / fermionic_error if ratio_valid else math.nan
            rows.append(
                {
                    **base,
                    "graph_status": "success",
                    "graph_source": graph_source,
                    "graph_unavailable_reason": "",
                    **metrics,
                    "ratio_valid": ratio_valid,
                    "jw_magnitude_to_fermionic_ratio": ratio,
                    "log10_jw_magnitude_to_fermionic_ratio": (
                        math.log10(ratio) if ratio_valid else math.nan
                    ),
                }
            )
            print(
                "  "
                f"source={graph_source}, Cw={metrics[PRIMARY_METRIC]:.6f}, "
                f"ratio={ratio if ratio_valid else 'excluded'}"
            )
        except Exception as error:  # Preserve full-sweep auditability.
            rows.append(
                {
                    **base,
                    "graph_status": "error",
                    "graph_source": "error",
                    "graph_unavailable_reason": f"{type(error).__name__}: {error}",
                    "ratio_valid": False,
                }
            )
            print(f"  ERROR: {type(error).__name__}: {error}")

    return pd.DataFrame(rows).sort_values(
        ["molecule", "n_qubits", "active_occupied", "active_vacant"]
    ).reset_index(drop=True)


def fixed_effect_residuals(
    frame: pd.DataFrame,
    values: np.ndarray,
    columns: Sequence[str],
) -> tuple[np.ndarray, int]:
    design_parts = [np.ones((len(frame), 1), dtype=float)]
    for column in columns:
        indicators = pd.get_dummies(frame[column], dtype=float)
        if indicators.shape[1] > 1:
            design_parts.append(indicators.iloc[:, 1:].to_numpy(dtype=float))
    design = np.column_stack(design_parts)
    fitted = design @ np.linalg.lstsq(design, values, rcond=None)[0]
    n_controls = int(np.linalg.matrix_rank(design) - 1)
    return values - fitted, n_controls


def correlation_statistics(
    frame: pd.DataFrame,
    metric: str,
    scope: str,
    fixed_effects: Sequence[str] = (),
) -> dict[str, Any]:
    usable = frame[
        frame["ratio_valid"]
        & np.isfinite(frame[metric])
        & np.isfinite(frame["log10_jw_magnitude_to_fermionic_ratio"])
    ].copy()
    x = usable[metric].to_numpy(dtype=float)
    y = usable["log10_jw_magnitude_to_fermionic_ratio"].to_numpy(dtype=float)
    n_controls = 0
    if fixed_effects:
        x, n_controls_x = fixed_effect_residuals(usable, x, fixed_effects)
        y, n_controls_y = fixed_effect_residuals(usable, y, fixed_effects)
        if n_controls_x != n_controls_y:
            raise RuntimeError("Fixed-effect design ranks unexpectedly differ.")
        n_controls = n_controls_x

    pearson_r = math.nan
    pearson_p = math.nan
    spearman_rho = math.nan
    spearman_p = math.nan
    slope = math.nan
    intercept = math.nan
    if len(usable) >= 3 and np.std(x) > 0.0 and np.std(y) > 0.0:
        pearson_r = float(np.corrcoef(x, y)[0, 1])
        degrees_freedom = len(usable) - n_controls - 2
        if degrees_freedom > 0 and abs(pearson_r) < 1.0:
            statistic = pearson_r * math.sqrt(
                degrees_freedom / (1.0 - pearson_r**2)
            )
            pearson_p = float(2.0 * stats.t.sf(abs(statistic), degrees_freedom))
        elif degrees_freedom > 0:
            pearson_p = 0.0
        spearman = stats.spearmanr(x, y)
        spearman_rho = float(spearman.statistic)
        if not fixed_effects:
            spearman_p = float(spearman.pvalue)
        elif degrees_freedom > 0 and abs(spearman_rho) < 1.0:
            rank_statistic = spearman_rho * math.sqrt(
                degrees_freedom / (1.0 - spearman_rho**2)
            )
            spearman_p = float(
                2.0 * stats.t.sf(abs(rank_statistic), degrees_freedom)
            )
        elif degrees_freedom > 0:
            spearman_p = 0.0
        slope, intercept = (float(value) for value in np.polyfit(x, y, 1))

    return {
        "scope": scope,
        "centralization_metric": metric,
        "fixed_effects": "+".join(fixed_effects) if fixed_effects else "none",
        "cases": len(usable),
        "molecules": int(usable["molecule"].nunique()),
        "active_space_sizes": int(usable["n_qubits"].nunique()),
        "n_controls": n_controls,
        "pearson_r_vs_log10_advantage": pearson_r,
        "pearson_p_value": pearson_p,
        "spearman_rho_vs_advantage": spearman_rho,
        "spearman_p_value": spearman_p,
        "ols_log10_advantage_slope": slope,
        "ols_intercept": intercept,
        "predicted_negative_direction": bool(pearson_r < 0.0),
        "negative_and_p_below_0_05": bool(pearson_r < 0.0 and pearson_p < 0.05),
    }


def build_correlation_table(cases: pd.DataFrame) -> pd.DataFrame:
    matched = cases[cases["molecule"].isin(MATCHED_MOLECULES)]
    direct = cases[cases["graph_source"] == "active_tensor"]
    rows = [
        correlation_statistics(matched, PRIMARY_METRIC, "matched_heteronuclear"),
        correlation_statistics(
            matched,
            PRIMARY_METRIC,
            "matched_heteronuclear",
            fixed_effects=("n_qubits",),
        ),
        correlation_statistics(
            matched,
            PRIMARY_METRIC,
            "matched_heteronuclear",
            fixed_effects=("molecule", "n_qubits"),
        ),
        correlation_statistics(cases, PRIMARY_METRIC, "full_hgbs5_sweep"),
        correlation_statistics(
            cases,
            PRIMARY_METRIC,
            "full_hgbs5_sweep",
            fixed_effects=("n_qubits",),
        ),
        correlation_statistics(
            cases,
            PRIMARY_METRIC,
            "full_hgbs5_sweep",
            fixed_effects=("molecule", "n_qubits"),
        ),
        correlation_statistics(
            direct,
            PRIMARY_METRIC,
            "direct_active_tensors_only",
        ),
        correlation_statistics(
            direct,
            PRIMARY_METRIC,
            "direct_active_tensors_only",
            fixed_effects=("n_qubits",),
        ),
    ]

    # Transparent alternatives check which side of the bipartite graph, if
    # either, carries an association and whether partition normalization
    # changes the result.
    for metric in (
        "occupied_weighted_freeman_centralization",
        "virtual_weighted_freeman_centralization",
        "mean_partition_weighted_freeman_centralization",
        "whole_graph_same_weight_star_centralization",
    ):
        rows.append(correlation_statistics(cases, metric, "full_hgbs5_sweep"))
    return pd.DataFrame(rows)


def build_coverage_table(cases: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for molecule, group in cases.groupby("molecule", sort=True):
        rows.append(
            {
                "molecule": molecule,
                "input_cases": len(group),
                "graph_cases": int((group["graph_status"] == "success").sum()),
                "direct_tensor_cases": int(
                    (group["graph_source"] == "active_tensor").sum()
                ),
                "reconstructed_cases": int(
                    (
                        group["graph_source"]
                        == "reconstructed_from_full_molecular_pickle"
                    ).sum()
                ),
                "unavailable_cases": int(
                    (group["graph_status"] == "unavailable").sum()
                ),
                "error_cases": int((group["graph_status"] == "error").sum()),
                "valid_ratio_cases": int(group["ratio_valid"].sum()),
            }
        )
    return pd.DataFrame(rows)


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "grid.linewidth": 0.7,
        }
    )


def correlation_lookup(
    correlations: pd.DataFrame,
    scope: str,
    fixed_effects: str,
) -> pd.Series:
    selected = correlations[
        (correlations["scope"] == scope)
        & (correlations["centralization_metric"] == PRIMARY_METRIC)
        & (correlations["fixed_effects"] == fixed_effects)
    ]
    if len(selected) != 1:
        raise RuntimeError(
            f"Expected one correlation row for {scope}/{fixed_effects}."
        )
    return selected.iloc[0]


def make_figure(
    cases: pd.DataFrame,
    correlations: pd.DataFrame,
    output_stem: Path,
    formats: Sequence[str],
    dpi: int,
) -> None:
    configure_plot_style()
    valid = cases[cases["ratio_valid"]].copy()
    matched = valid[valid["molecule"].isin(MATCHED_MOLECULES)]
    full_stats = correlation_lookup(correlations, "full_hgbs5_sweep", "none")
    full_size_stats = correlation_lookup(
        correlations,
        "full_hgbs5_sweep",
        "n_qubits",
    )
    matched_stats = correlation_lookup(
        correlations,
        "matched_heteronuclear",
        "none",
    )
    matched_size_stats = correlation_lookup(
        correlations,
        "matched_heteronuclear",
        "n_qubits",
    )

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.25), constrained_layout=True)
    colors = {"B-H": "#264653", "Be-H": "#e76f51", "Li-H": "#2a9d8f"}
    for molecule in MATCHED_MOLECULES:
        subset = matched[matched["molecule"] == molecule].sort_values("n_qubits")
        axes[0].plot(
            subset[PRIMARY_METRIC],
            subset["jw_magnitude_to_fermionic_ratio"],
            marker="o",
            linewidth=1.2,
            markersize=4.5,
            color=colors[molecule],
            label=molecule.replace("-", ""),
        )
    axes[0].axhline(1.0, color="#555555", linestyle="--", linewidth=1.0)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Bipartite weighted Freeman centralization")
    axes[0].set_ylabel(r"$E_{\mathrm{JW,mag}}/E_{\mathrm{fermionic}}$")
    axes[0].set_title(
        "Matched heteronuclear cases\n"
        f"raw r = {matched_stats['pearson_r_vs_log10_advantage']:.2f}; "
        "size-FE r = "
        f"{matched_size_stats['pearson_r_vs_log10_advantage']:.2f}"
    )
    axes[0].legend(frameon=False, title="Molecule")

    scatter = axes[1].scatter(
        valid[PRIMARY_METRIC],
        valid["jw_magnitude_to_fermionic_ratio"],
        c=valid["n_qubits"],
        cmap="viridis",
        s=30,
        alpha=0.82,
        edgecolors="white",
        linewidths=0.35,
    )
    if len(valid) >= 2:
        x_values = valid[PRIMARY_METRIC].to_numpy(dtype=float)
        y_values = valid[
            "log10_jw_magnitude_to_fermionic_ratio"
        ].to_numpy(dtype=float)
        slope, intercept = np.polyfit(x_values, y_values, 1)
        grid = np.linspace(x_values.min(), x_values.max(), 200)
        axes[1].plot(grid, 10.0 ** (intercept + slope * grid), color="#9b2226")
    axes[1].axhline(1.0, color="#555555", linestyle="--", linewidth=1.0)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Bipartite weighted Freeman centralization")
    axes[1].set_ylabel(r"$E_{\mathrm{JW,mag}}/E_{\mathrm{fermionic}}$")
    axes[1].set_title(
        "Available HGBS-5 sweep\n"
        f"raw r = {full_stats['pearson_r_vs_log10_advantage']:.2f}; "
        f"size-FE r = {full_size_stats['pearson_r_vs_log10_advantage']:.2f}"
    )
    colorbar = fig.colorbar(scatter, ax=axes[1], pad=0.02)
    colorbar.set_label("Active spin orbitals")

    for output_format in formats:
        fig.savefig(output_stem.with_suffix(f".{output_format}"), dpi=dpi)
    plt.close(fig)


def format_stat(row: pd.Series) -> str:
    return (
        f"n={int(row['cases'])}, Pearson r="
        f"{row['pearson_r_vs_log10_advantage']:.3f} "
        f"(p={row['pearson_p_value']:.3g}), Spearman ρ="
        f"{row['spearman_rho_vs_advantage']:.3f} "
        f"(p={row['spearman_p_value']:.3g})"
    )


def write_report(
    path: Path,
    cases: pd.DataFrame,
    correlations: pd.DataFrame,
    coverage: pd.DataFrame,
) -> None:
    matched_raw = correlation_lookup(
        correlations, "matched_heteronuclear", "none"
    )
    matched_size = correlation_lookup(
        correlations, "matched_heteronuclear", "n_qubits"
    )
    full_raw = correlation_lookup(correlations, "full_hgbs5_sweep", "none")
    full_size = correlation_lookup(
        correlations, "full_hgbs5_sweep", "n_qubits"
    )
    full_two_way = correlation_lookup(
        correlations,
        "full_hgbs5_sweep",
        "molecule+n_qubits",
    )
    direct_raw = correlation_lookup(
        correlations,
        "direct_active_tensors_only",
        "none",
    )
    direct_size = correlation_lookup(
        correlations,
        "direct_active_tensors_only",
        "n_qubits",
    )
    virtual_side_sensitivity = correlations[
        (
            correlations["centralization_metric"]
            == "virtual_weighted_freeman_centralization"
        )
        & (correlations["scope"] == "full_hgbs5_sweep")
        & (correlations["fixed_effects"] == "none")
    ].iloc[0]
    graph_cases = int((cases["graph_status"] == "success").sum())
    input_cases = len(cases)
    valid_cases = int(cases["ratio_valid"].sum())
    unavailable = cases[cases["graph_status"] != "success"]

    valid = cases[cases["ratio_valid"]]
    centralization_size = stats.pearsonr(
        valid[PRIMARY_METRIC],
        valid["n_qubits"],
    )
    advantage_size = stats.pearsonr(
        valid["log10_jw_magnitude_to_fermionic_ratio"],
        valid["n_qubits"],
    )

    robust_negative = (
        float(matched_size["pearson_r_vs_log10_advantage"]) < 0.0
        and float(matched_size["pearson_p_value"]) < 0.05
        and float(full_two_way["pearson_r_vs_log10_advantage"]) < 0.0
        and float(full_two_way["pearson_p_value"]) < 0.05
    )
    raw_r = float(full_raw["pearson_r_vs_log10_advantage"])
    raw_p = float(full_raw["pearson_p_value"])
    if robust_negative:
        conclusion = (
            "The matched and fixed-effect tests support the proposed negative "
            "association."
        )
    elif raw_r < 0.0 and raw_p < 0.05:
        conclusion = (
            "The hypothesized negative relationship is not robustly supported. "
            "There is a significant negative aggregate correlation, but it is "
            "absent in the matched heteronuclear cases and disappears after "
            "active-space-size and molecule fixed effects are removed."
        )
    elif raw_r < 0.0:
        conclusion = (
            "The hypothesized negative relationship is not supported. The "
            "full-sweep estimate points weakly in the proposed direction but "
            "is statistically unresolved, and the matched/fixed-effect tests "
            "do not recover a negative association."
        )
    else:
        conclusion = (
            "The hypothesized negative relationship is not supported. The "
            "aggregate estimate is not even in the proposed direction, and "
            "the matched/fixed-effect tests do not show a stable negative "
            "association."
        )

    lines = [
        "# Coefficient-weighted fermionic graph centralization",
        "",
        "## Definition",
        "",
        "Vertices are active spin orbitals, partitioned at the Hartree--Fock "
        "occupied/virtual boundary. For each complete Hermitian fermionic "
        "parent, the absolute coefficient of its lexicographically canonical "
        "component is added to every occupied--virtual orbital pair jointly "
        "touched by that parent. Parents confined to one partition do not "
        "create edges. This matches the occupied-to-virtual incidence network "
        "behind the qualitative hub-dominance hypothesis.",
        "",
        "For one partition `P` with `m` orbitals, strengths `s_i = sum_j "
        "w_ij`, maximum strength `s_max`, and total bipartite edge weight `W`, "
        "its weighted Freeman centralization is",
        "",
        "`C_P^w = sum_{i in P}(s_max - s_i) / ((m - 1) W)`.",
        "",
        "The graph-level statistic is `max(C_occupied^w, C_virtual^w)`: a graph "
        "is hub dominated if either partition concentrates its incident "
        "coefficient weight on one orbital. Partition-wise normalization is "
        "coefficient-scale invariant, equals zero for equal strength within "
        "each partition, equals one for a partition star, and does not mistake "
        "an unequal occupied/virtual node count for hub structure.",
        "",
        "The response is `E_JW-magnitude / E_fermionic`; values above one favor "
        "fermionic-aware ordering. Correlations use its base-10 logarithm, "
        "because the ratios span orders of magnitude.",
        "",
        "## Coverage",
        "",
        f"- Graphs built: {graph_cases}/{input_cases} HGBS-5 cases.",
        f"- Numerically valid graph/ratio pairs: {valid_cases}.",
        f"- Direct active tensors: {int((cases['graph_source'] == 'active_tensor').sum())}.",
        "- Reconstructed exactly in memory from retained full molecular "
        f"pickles: {int((cases['graph_source'] == 'reconstructed_from_full_molecular_pickle').sum())}.",
        f"- Unavailable/error cases: {len(unavailable)}.",
        "",
        "## Correlation results",
        "",
        f"- Matched B-H/Be-H/Li-H, unadjusted: {format_stat(matched_raw)}.",
        "- Matched B-H/Be-H/Li-H, after removing active-space-size fixed "
        f"effects: {format_stat(matched_size)}.",
        f"- Full available sweep, unadjusted: {format_stat(full_raw)}.",
        "- Full available sweep, after removing active-space-size fixed "
        f"effects: {format_stat(full_size)}.",
        "- Full available sweep, after removing molecule and active-space-size "
        f"fixed effects: {format_stat(full_two_way)}.",
        "- Direct active tensors only, unadjusted: "
        f"{format_stat(direct_raw)}.",
        "- Direct active tensors only, after removing active-space-size fixed "
        f"effects: {format_stat(direct_size)}.",
        "",
        "## Interpretation",
        "",
        conclusion,
        "",
        "As a confounding diagnostic, centralization changes with active-space "
        "size "
        f"(Pearson r={centralization_size.statistic:.3f}, "
        f"p={centralization_size.pvalue:.3g}), while log ordering advantage "
        "changes with active-space size "
        f"(r={advantage_size.statistic:.3f}, p={advantage_size.pvalue:.3g}).",
        "",
        "The virtual-side centralization—the most direct version of the "
        "original virtual-hub claim—also gives no resolved full-sweep "
        "association: "
        f"r={virtual_side_sensitivity['pearson_r_vs_log10_advantage']:.3f} "
        f"(p={virtual_side_sensitivity['pearson_p_value']:.3g}). The case table "
        "also reports occupied-side, partition-mean, and non-partition-normalized "
        "variants for sensitivity analysis.",
        "",
        "The fixed-effect rows are the more relevant tests of the mechanistic "
        "claim: they compare cases at matched active-space size, and the "
        "two-way version also removes stable molecule-family differences. "
        "These are observational associations, not a causal identification of "
        "the ordering mechanism.",
        "",
        "## Missing-source limitation",
        "",
    ]
    if unavailable.empty:
        lines.append("All input cases had a recoverable coefficient tensor.")
    else:
        missing_molecules = ", ".join(sorted(unavailable["molecule"].unique()))
        lines.extend(
            [
                "Coefficient tensors and full molecular pickles were not "
                f"retained for: {missing_molecules}. Those rows remain in the "
                "case and coverage tables with an explicit unavailable reason; "
                "they are not silently dropped.",
            ]
        )
    lines.extend(
        [
            "",
            "See `fermionic_graph_case_metrics.csv` for the auditable case "
            "data, `matched_heteronuclear_cases.csv` for the initial matched "
            "comparison, `correlation_summary.csv` for all statistical tests, "
            "and `coverage_summary.csv` for source coverage.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_arguments()
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive.")

    cases = read_cases(args.input, args.basis)
    analyzed = analyze_cases(cases, args.tolerance)
    correlations = build_correlation_table(analyzed)
    coverage = build_coverage_table(analyzed)
    matched = analyzed[analyzed["molecule"].isin(MATCHED_MOLECULES)].copy()

    args.outdir.mkdir(parents=True, exist_ok=True)
    analyzed.to_csv(args.outdir / "fermionic_graph_case_metrics.csv", index=False)
    matched.to_csv(args.outdir / "matched_heteronuclear_cases.csv", index=False)
    correlations.to_csv(args.outdir / "correlation_summary.csv", index=False)
    coverage.to_csv(args.outdir / "coverage_summary.csv", index=False)
    make_figure(
        analyzed,
        correlations,
        args.outdir / "centralization_vs_ordering_advantage",
        args.formats,
        args.dpi,
    )
    write_report(
        args.outdir / "centralization_report.md",
        analyzed,
        correlations,
        coverage,
    )

    print(f"Wrote {args.outdir / 'fermionic_graph_case_metrics.csv'}")
    print(f"Wrote {args.outdir / 'matched_heteronuclear_cases.csv'}")
    print(f"Wrote {args.outdir / 'correlation_summary.csv'}")
    print(f"Wrote {args.outdir / 'coverage_summary.csv'}")
    print(f"Wrote {args.outdir / 'centralization_report.md'}")


if __name__ == "__main__":
    main()
