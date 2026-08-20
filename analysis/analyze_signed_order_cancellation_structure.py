#!/usr/bin/env python3
"""Test signed-parent order structure against fermionic-ordering advantage.

For each available HGBS-5 active-space Hamiltonian, this script constructs a
parent interaction network.  Nodes are complete Hermitian fermionic parents.
Parents i and j receive interaction weight

    w_ij = |c_i c_j| |S_i intersect S_j|,

where c is the canonical signed parent coefficient and S is spin-orbital
support.  Shared support is an inexpensive pre-BCH interaction proxy; it is
not claimed to be an exact fermionic-commutator norm.

The primary diagnostics ask whether the signed-ascending parent sequence puts
large interaction weight locally, especially for opposite-sign parents.  The
same diagnostics are computed in the fermionic magnitude-descending sequence
as an order control.  Results are compared with the JW-magnitude/fermionic
one-minus-overlap ratio and, for the independently benchmarked subset, the
HF-state BCH2 cancellation ratio.
"""

# ruff: noqa: E402  # Select a writable Matplotlib cache before importing it.

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

_MPL_CACHE = Path(tempfile.gettempdir()) / "qhat-matplotlib-cache"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openfermion import get_fermion_operator
from scipy import stats

try:
    from qhat.analysis.analyze_fermionic_graph_centralization import (
        MolecularDataCache,
        load_case_interaction,
        read_cases,
    )
    from qhat.analysis.benchmark_L_sweep_trotter import (
        HermitianFermionTerm,
        build_hermitian_fermion_terms,
        clean_fermion_operator,
    )
except ImportError:
    from analyze_fermionic_graph_centralization import (
        MolecularDataCache,
        load_case_interaction,
        read_cases,
    )
    from benchmark_L_sweep_trotter import (
        HermitianFermionTerm,
        build_hermitian_fermion_terms,
        clean_fermion_operator,
    )


DEFAULT_INPUT = Path(
    "analysis/fermionic_aware_performance/fermionic_aware_case_performance.csv"
)
DEFAULT_OUTDIR = Path("analysis/signed_order_cancellation_structure")
DEFAULT_CANCELLATION_INPUTS = (
    Path("analysis/fermionic_structure_ablation_hgbs5.csv"),
    Path("analysis/fermionic_structure_ablation_extension_hgbs5.csv"),
)
DEFAULT_TOLERANCE = 1.0e-12
MATCHED_MOLECULES = ("B-H", "Be-H", "Li-H")

PRIMARY_FEATURES = (
    "opposite_sign_interaction_weight_fraction",
    "signed_order_weighted_mean_span_fraction",
    "signed_order_weighted_local_5pct_fraction",
    "signed_order_weighted_local_10pct_fraction",
    "signed_order_opposite_sign_local_10pct_interaction_fraction",
    "signed_order_opposite_sign_local_10pct_enrichment",
    "signed_vs_magnitude_weighted_span_improvement",
    "signed_vs_magnitude_weighted_local_10pct_improvement",
    "signed_vs_magnitude_opposite_sign_local_10pct_improvement",
)

PREDICTED_DIRECTIONS = {
    "opposite_sign_interaction_weight_fraction": 1,
    "signed_order_weighted_mean_span_fraction": -1,
    "signed_order_weighted_local_5pct_fraction": 1,
    "signed_order_weighted_local_10pct_fraction": 1,
    "signed_order_opposite_sign_local_10pct_interaction_fraction": 1,
    "signed_order_opposite_sign_local_10pct_enrichment": 1,
    "signed_vs_magnitude_weighted_span_improvement": 1,
    "signed_vs_magnitude_weighted_local_10pct_improvement": 1,
    "signed_vs_magnitude_opposite_sign_local_10pct_improvement": 1,
}

FEATURE_LABELS = {
    "opposite_sign_interaction_weight_fraction": (
        "Opposite-sign share of interaction weight"
    ),
    "signed_order_weighted_mean_span_fraction": (
        "Signed-order weighted mean span"
    ),
    "signed_order_weighted_local_5pct_fraction": (
        "Signed-order interaction weight within 5%"
    ),
    "signed_order_weighted_local_10pct_fraction": (
        "Signed-order interaction weight within 10%"
    ),
    "signed_order_opposite_sign_local_10pct_interaction_fraction": (
        "Opposite-sign local interaction weight (10%)"
    ),
    "signed_order_opposite_sign_local_10pct_enrichment": (
        "Opposite-sign enrichment in 10% window"
    ),
    "signed_vs_magnitude_weighted_span_improvement": (
        "Signed-vs-magnitude span improvement"
    ),
    "signed_vs_magnitude_weighted_local_10pct_improvement": (
        "Signed-vs-magnitude locality improvement (10%)"
    ),
    "signed_vs_magnitude_opposite_sign_local_10pct_improvement": (
        "Signed-vs-magnitude opposite-sign locality improvement"
    ),
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--basis", default="hgbs-5")
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    parser.add_argument(
        "--cancellation-input",
        action="append",
        type=Path,
        dest="cancellation_inputs",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def parent_coefficient(term: HermitianFermionTerm, tolerance: float) -> float:
    """Canonical signed coefficient used by the production parent ordering."""
    canonical_key = min(term.component_keys)
    value = complex(term.operator.terms[canonical_key])
    if abs(value.imag) > max(100.0 * tolerance, 1.0e-10):
        raise ValueError(f"Non-real canonical parent coefficient: {value!r}")
    return float(value.real)


def fermionic_parent_order_indices(
    terms: Sequence[HermitianFermionTerm],
    coefficients: np.ndarray,
    ordering_method: str,
) -> list[int]:
    """Reproduce the deterministic signed/magnitude parent sort keys."""
    lexicographic_keys = [tuple(sorted(term.component_keys)) for term in terms]
    if ordering_method == "signed_ascending":
        return sorted(
            range(len(terms)),
            key=lambda index: (coefficients[index], lexicographic_keys[index]),
        )
    if ordering_method == "magnitude_descending":
        return sorted(
            range(len(terms)),
            key=lambda index: (-abs(coefficients[index]), lexicographic_keys[index]),
        )
    raise ValueError(f"Unsupported ordering method: {ordering_method!r}")


def parent_support(term: HermitianFermionTerm) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                orbital
                for component in term.component_keys
                for orbital, _ in component
            }
        )
    )


def random_order_local_pair_probability(n_nodes: int, cutoff: int) -> float:
    """Probability that two distinct random positions are at most cutoff apart."""
    if n_nodes <= 1:
        return math.nan
    k = min(max(int(cutoff), 0), n_nodes - 1)
    return float(k * (2 * n_nodes - k - 1) / (n_nodes * (n_nodes - 1)))


def random_order_mean_span_fraction(n_nodes: int) -> float:
    """Expected |position_i-position_j|/(N-1) under a random permutation."""
    if n_nodes <= 1:
        return math.nan
    return float((n_nodes + 1) / (3 * (n_nodes - 1)))


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator <= 0.0 or not math.isfinite(denominator):
        return math.nan
    return float(numerator / denominator)


def order_interaction_metrics(
    coefficients: np.ndarray,
    shared_support_counts: np.ndarray,
    order: Sequence[int],
    prefix: str,
) -> dict[str, float]:
    """Compute coefficient/shared-support-weighted locality in one order."""
    coefficients = np.asarray(coefficients, dtype=float)
    shared_support_counts = np.asarray(shared_support_counts)
    order_array = np.asarray(order, dtype=int)
    n_nodes = len(coefficients)
    if shared_support_counts.shape != (n_nodes, n_nodes):
        raise ValueError("shared_support_counts must be N by N.")
    if sorted(order_array.tolist()) != list(range(n_nodes)):
        raise ValueError("order must be a permutation of parent indices.")

    ordered_coefficients = coefficients[order_array]
    ordered_shared = shared_support_counts[np.ix_(order_array, order_array)]
    absolute_coefficients = np.abs(ordered_coefficients)
    signs = np.sign(ordered_coefficients)
    cutoffs = {
        percent: max(1, int(math.ceil(percent * n_nodes / 100.0)))
        for percent in (1, 5, 10)
    }

    total_weight = 0.0
    opposite_weight = 0.0
    weighted_span = 0.0
    opposite_weighted_span = 0.0
    local_weight = {percent: 0.0 for percent in cutoffs}
    local_opposite_weight = {percent: 0.0 for percent in cutoffs}

    for separation in range(1, n_nodes):
        overlap = np.diagonal(ordered_shared, offset=separation).astype(float)
        mask = overlap > 0.0
        if not np.any(mask):
            continue
        weights = (
            absolute_coefficients[:-separation]
            * absolute_coefficients[separation:]
            * overlap
        )
        selected_weight = float(weights[mask].sum())
        opposite = signs[:-separation] != signs[separation:]
        selected_opposite_weight = float(weights[mask & opposite].sum())
        total_weight += selected_weight
        opposite_weight += selected_opposite_weight
        weighted_span += separation * selected_weight
        opposite_weighted_span += separation * selected_opposite_weight
        for percent, cutoff in cutoffs.items():
            if separation <= cutoff:
                local_weight[percent] += selected_weight
                local_opposite_weight[percent] += selected_opposite_weight

    span_denominator = max(1, n_nodes - 1)
    metrics = {
        f"{prefix}_total_interaction_weight": total_weight,
        f"{prefix}_opposite_sign_interaction_weight_fraction": safe_ratio(
            opposite_weight,
            total_weight,
        ),
        f"{prefix}_weighted_mean_span_fraction": safe_ratio(
            weighted_span,
            total_weight * span_denominator,
        ),
        f"{prefix}_opposite_sign_weighted_mean_span_fraction": safe_ratio(
            opposite_weighted_span,
            opposite_weight * span_denominator,
        ),
    }
    random_span = random_order_mean_span_fraction(n_nodes)
    metrics[f"{prefix}_weighted_span_improvement_vs_random"] = (
        random_span - metrics[f"{prefix}_weighted_mean_span_fraction"]
    )

    for percent, cutoff in cutoffs.items():
        overall_local = safe_ratio(local_weight[percent], total_weight)
        opposite_local_total = safe_ratio(
            local_opposite_weight[percent],
            total_weight,
        )
        opposite_local_conditional = safe_ratio(
            local_opposite_weight[percent],
            opposite_weight,
        )
        local_opposite_share = safe_ratio(
            local_opposite_weight[percent],
            local_weight[percent],
        )
        global_opposite_share = safe_ratio(opposite_weight, total_weight)
        random_local = random_order_local_pair_probability(n_nodes, cutoff)
        metrics[f"{prefix}_weighted_local_{percent}pct_fraction"] = overall_local
        metrics[
            f"{prefix}_opposite_sign_local_{percent}pct_interaction_fraction"
        ] = opposite_local_total
        metrics[
            f"{prefix}_opposite_sign_local_{percent}pct_conditional_fraction"
        ] = opposite_local_conditional
        metrics[f"{prefix}_opposite_sign_local_{percent}pct_enrichment"] = safe_ratio(
            local_opposite_share,
            global_opposite_share,
        )
        metrics[f"{prefix}_weighted_local_{percent}pct_enrichment_vs_random"] = (
            safe_ratio(overall_local, random_local)
        )

    return metrics


def build_parent_order_metrics(
    terms: Sequence[HermitianFermionTerm],
    n_qubits: int,
    tolerance: float,
) -> dict[str, float | int]:
    coefficients = np.asarray(
        [parent_coefficient(term, tolerance) for term in terms],
        dtype=float,
    )
    if np.any(np.abs(coefficients) <= tolerance):
        raise ValueError("Hermitian parent with zero canonical coefficient found.")
    incidence = np.zeros((len(terms), n_qubits), dtype=np.uint8)
    for index, term in enumerate(terms):
        incidence[index, list(parent_support(term))] = 1
    shared_support = incidence @ incidence.T
    np.fill_diagonal(shared_support, 0)

    signed_order = fermionic_parent_order_indices(
        terms,
        coefficients,
        "signed_ascending",
    )
    magnitude_order = fermionic_parent_order_indices(
        terms,
        coefficients,
        "magnitude_descending",
    )
    signed = order_interaction_metrics(
        coefficients,
        shared_support,
        signed_order,
        "signed_order",
    )
    magnitude = order_interaction_metrics(
        coefficients,
        shared_support,
        magnitude_order,
        "magnitude_order",
    )

    # The global opposite-sign share is order invariant.  Keep one concise
    # primary name and verify that both order calculations agree.
    signed_opposite = signed[
        "signed_order_opposite_sign_interaction_weight_fraction"
    ]
    magnitude_opposite = magnitude[
        "magnitude_order_opposite_sign_interaction_weight_fraction"
    ]
    if not math.isclose(signed_opposite, magnitude_opposite, abs_tol=1.0e-12):
        raise RuntimeError("Order-invariant interaction weight changed by permutation.")

    signed_span = signed["signed_order_weighted_mean_span_fraction"]
    magnitude_span = magnitude["magnitude_order_weighted_mean_span_fraction"]
    signed_local = signed["signed_order_weighted_local_10pct_fraction"]
    magnitude_local = magnitude["magnitude_order_weighted_local_10pct_fraction"]
    signed_opposite_local = signed[
        "signed_order_opposite_sign_local_10pct_interaction_fraction"
    ]
    magnitude_opposite_local = magnitude[
        "magnitude_order_opposite_sign_local_10pct_interaction_fraction"
    ]

    interaction_edges = int(np.count_nonzero(np.triu(shared_support, k=1)))
    return {
        "number_of_fermionic_parents": len(terms),
        "parent_interaction_edges": interaction_edges,
        "parent_interaction_edge_density": safe_ratio(
            interaction_edges,
            len(terms) * (len(terms) - 1) / 2.0,
        ),
        "opposite_sign_interaction_weight_fraction": signed_opposite,
        **signed,
        **magnitude,
        "signed_vs_magnitude_weighted_span_improvement": (
            magnitude_span - signed_span
        ),
        "signed_vs_magnitude_weighted_local_10pct_improvement": (
            signed_local - magnitude_local
        ),
        "signed_vs_magnitude_opposite_sign_local_10pct_improvement": (
            signed_opposite_local - magnitude_opposite_local
        ),
    }


def analyze_cases(cases: pd.DataFrame, tolerance: float) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    cache = MolecularDataCache()
    for _, source in cases.iterrows():
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
        base = {
            key: source[key]
            for key in (
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
            )
        }
        base.update(
            {
                "ratio_valid": ratio_valid,
                "jw_magnitude_to_fermionic_ratio": ratio,
                "log10_jw_magnitude_to_fermionic_ratio": (
                    math.log10(ratio) if ratio_valid else math.nan
                ),
                "log10_fermionic_infidelity": (
                    math.log10(fermionic_error)
                    if math.isfinite(fermionic_error) and fermionic_error > 0.0
                    else math.nan
                ),
            }
        )
        print(f"{source['case_id']}: building signed parent-interaction metrics")
        try:
            interaction, data_source, reason = load_case_interaction(source, cache)
            if interaction is None:
                rows.append(
                    {
                        **base,
                        "structure_status": "unavailable",
                        "structure_source": data_source,
                        "structure_unavailable_reason": reason,
                    }
                )
                print("  unavailable")
                continue

            fermion_hamiltonian = clean_fermion_operator(
                get_fermion_operator(interaction),
                tolerance,
            )
            terms = build_hermitian_fermion_terms(fermion_hamiltonian, tolerance)
            metrics = build_parent_order_metrics(
                terms,
                int(source["n_qubits"]),
                tolerance,
            )
            rows.append(
                {
                    **base,
                    "structure_status": "success",
                    "structure_source": data_source,
                    "structure_unavailable_reason": "",
                    **metrics,
                }
            )
            print(
                "  "
                f"parents={metrics['number_of_fermionic_parents']}, "
                f"opp={metrics['opposite_sign_interaction_weight_fraction']:.4f}, "
                f"local10={metrics['signed_order_weighted_local_10pct_fraction']:.4f}, "
                f"ratio={ratio if ratio_valid else 'excluded'}"
            )
        except Exception as error:  # Preserve full-sweep auditability.
            rows.append(
                {
                    **base,
                    "structure_status": "error",
                    "structure_source": "error",
                    "structure_unavailable_reason": (
                        f"{type(error).__name__}: {error}"
                    ),
                }
            )
            print(f"  ERROR: {type(error).__name__}: {error}")

    return pd.DataFrame(rows).sort_values(
        ["molecule", "n_qubits", "active_occupied", "active_vacant"]
    ).reset_index(drop=True)


def load_cancellation_measurements(paths: Iterable[Path]) -> pd.DataFrame:
    """Load one deterministic signed-reference row per independent case."""
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        required = {
            "case_id",
            "schedule",
            "status",
            "bch_cancellation_ratio",
            "bch2_hf_state_norm",
            "bch_pair_abs_sum",
        }
        if not required.issubset(frame.columns):
            continue
        selected = frame[
            frame["status"].astype(str).str.lower().isin({"success", "ok"})
            & frame["schedule"].isin(
                {
                    "fermionic_signed_reference",
                    "fermionic_magnitude_reference",
                    "jw_signed_baseline",
                }
            )
        ].copy()
        selected["cancellation_source_csv"] = str(path)
        frames.append(selected)
    if not frames:
        return pd.DataFrame(columns=["case_id"])

    combined = pd.concat(frames, ignore_index=True)
    values = combined.pivot_table(
        index="case_id",
        columns="schedule",
        values="bch_cancellation_ratio",
        aggfunc="first",
    )
    signed_name = "fermionic_signed_reference"
    magnitude_name = "fermionic_magnitude_reference"
    jw_name = "jw_signed_baseline"
    if signed_name not in values:
        return pd.DataFrame(columns=["case_id"])
    output = pd.DataFrame(
        {
            "case_id": values.index,
            "signed_bch_cancellation_ratio": values[signed_name].to_numpy(),
        }
    )
    output["signed_bch_cancellation_strength"] = -np.log10(
        output["signed_bch_cancellation_ratio"]
    )
    if magnitude_name in values:
        magnitude = values[magnitude_name].reindex(values.index).to_numpy()
        output["fermionic_magnitude_bch_cancellation_ratio"] = magnitude
        output["signed_bch_cancellation_advantage_vs_fermionic_magnitude"] = (
            magnitude / output["signed_bch_cancellation_ratio"]
        )
        output[
            "log10_signed_bch_cancellation_advantage_vs_fermionic_magnitude"
        ] = np.log10(
            output[
                "signed_bch_cancellation_advantage_vs_fermionic_magnitude"
            ]
        )
    if jw_name in values:
        jw = values[jw_name].reindex(values.index).to_numpy()
        output["jw_signed_bch_cancellation_ratio"] = jw
        output["signed_bch_cancellation_advantage_vs_jw_signed"] = (
            jw / output["signed_bch_cancellation_ratio"]
        )
    source_by_case = combined.groupby("case_id")[
        "cancellation_source_csv"
    ].first()
    output["cancellation_source_csv"] = output["case_id"].map(source_by_case)
    return output.reset_index(drop=True)


def fixed_effect_residuals(
    frame: pd.DataFrame,
    values: np.ndarray,
    columns: Sequence[str],
) -> tuple[np.ndarray, int]:
    parts = [np.ones((len(frame), 1), dtype=float)]
    for column in columns:
        indicators = pd.get_dummies(frame[column], dtype=float)
        if indicators.shape[1] > 1:
            parts.append(indicators.iloc[:, 1:].to_numpy(dtype=float))
    design = np.column_stack(parts)
    fitted = design @ np.linalg.lstsq(design, values, rcond=None)[0]
    return values - fitted, int(np.linalg.matrix_rank(design) - 1)


def correlation_row(
    frame: pd.DataFrame,
    feature: str,
    outcome: str,
    scope: str,
    fixed_effects: Sequence[str] = (),
) -> dict[str, Any]:
    usable = frame[
        np.isfinite(frame[feature]) & np.isfinite(frame[outcome])
    ].copy()
    x = usable[feature].to_numpy(dtype=float)
    y = usable[outcome].to_numpy(dtype=float)
    n_controls = 0
    if fixed_effects:
        x, n_controls = fixed_effect_residuals(usable, x, fixed_effects)
        y, y_controls = fixed_effect_residuals(usable, y, fixed_effects)
        if n_controls != y_controls:
            raise RuntimeError("Fixed-effect design ranks differ.")

    pearson_r = pearson_p = spearman_rho = spearman_p = math.nan
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
        spearman_p = float(spearman.pvalue) if not fixed_effects else math.nan

    feature_direction = PREDICTED_DIRECTIONS.get(feature, 1)
    outcome_direction = -1 if outcome == "log10_fermionic_infidelity" else 1
    aligned_direction = feature_direction * outcome_direction
    return {
        "scope": scope,
        "outcome": outcome,
        "feature": feature,
        "feature_label": FEATURE_LABELS.get(feature, feature),
        "fixed_effects": "+".join(fixed_effects) if fixed_effects else "none",
        "cases": len(usable),
        "molecules": int(usable["molecule"].nunique()),
        "active_space_sizes": int(usable["n_qubits"].nunique()),
        "n_controls": n_controls,
        "predicted_direction": (
            "positive" if aligned_direction > 0 else "negative"
        ),
        "pearson_r": pearson_r,
        "pearson_p_value": pearson_p,
        "direction_aligned_pearson_r": aligned_direction * pearson_r,
        "spearman_rho": spearman_rho,
        "spearman_p_value": spearman_p,
    }


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    adjusted = pd.Series(np.nan, index=p_values.index, dtype=float)
    finite = p_values[np.isfinite(p_values)].sort_values()
    if finite.empty:
        return adjusted
    m = len(finite)
    raw = finite.to_numpy(dtype=float) * m / np.arange(1, m + 1)
    monotone = np.minimum.accumulate(raw[::-1])[::-1]
    adjusted.loc[finite.index] = np.minimum(monotone, 1.0)
    return adjusted


def build_correlation_table(cases: pd.DataFrame) -> pd.DataFrame:
    usable = cases[
        (cases["structure_status"] == "success") & cases["ratio_valid"]
    ].copy()
    matched = usable[usable["molecule"].isin(MATCHED_MOLECULES)]
    rows: list[dict[str, Any]] = []
    specifications = (
        (matched, "matched_BH_BeH_LiH", ()),
        (matched, "matched_BH_BeH_LiH", ("n_qubits",)),
        (matched, "matched_BH_BeH_LiH", ("molecule", "n_qubits")),
        (usable, "full_hgbs5_sweep", ()),
        (usable, "full_hgbs5_sweep", ("n_qubits",)),
        (usable, "full_hgbs5_sweep", ("molecule", "n_qubits")),
    )
    outcomes = (
        "log10_jw_magnitude_to_fermionic_ratio",
        "log10_fermionic_infidelity",
    )
    for outcome in outcomes:
        for frame, scope, fixed_effects in specifications:
            for feature in PRIMARY_FEATURES:
                rows.append(
                    correlation_row(
                        frame,
                        feature,
                        outcome,
                        scope,
                        fixed_effects,
                    )
                )
    table = pd.DataFrame(rows)
    table["pearson_fdr_q_value"] = table.groupby(
        ["scope", "outcome", "fixed_effects"],
        group_keys=False,
    )["pearson_p_value"].apply(benjamini_hochberg)
    table["supports_predicted_direction_at_q_0_05"] = (
        (table["direction_aligned_pearson_r"] > 0.0)
        & (table["pearson_fdr_q_value"] < 0.05)
    )
    return table


def build_cancellation_correlation_table(cases: pd.DataFrame) -> pd.DataFrame:
    measured = cases[np.isfinite(cases["signed_bch_cancellation_ratio"])].copy()
    outcomes = (
        "signed_bch_cancellation_strength",
        "log10_signed_bch_cancellation_advantage_vs_fermionic_magnitude",
    )
    rows = [
        correlation_row(
            measured,
            feature,
            outcome,
            "independent_bch_subset",
        )
        for outcome in outcomes
        for feature in PRIMARY_FEATURES
    ]
    # This row tests the mechanism-to-performance link directly: stronger
    # independently measured cancellation should accompany larger advantage.
    rows.append(
        correlation_row(
            measured,
            "signed_bch_cancellation_strength",
            "log10_jw_magnitude_to_fermionic_ratio",
            "independent_bch_subset",
        )
    )
    rows.append(
        correlation_row(
            measured,
            "signed_bch_cancellation_strength",
            "log10_fermionic_infidelity",
            "independent_bch_subset",
        )
    )
    table = pd.DataFrame(rows)
    table["pearson_fdr_q_value"] = table.groupby(
        ["scope", "outcome", "fixed_effects"],
        group_keys=False,
    )["pearson_p_value"].apply(benjamini_hochberg)
    return table


def matched_size_summary(cases: pd.DataFrame) -> pd.DataFrame:
    matched = cases[
        cases["molecule"].isin(MATCHED_MOLECULES)
        & (cases["structure_status"] == "success")
        & cases["ratio_valid"]
    ].copy()
    columns = [
        "n_qubits",
        "molecule",
        "case_id",
        "jw_magnitude_to_fermionic_ratio",
        *PRIMARY_FEATURES,
    ]
    return matched[columns].sort_values(["n_qubits", "molecule"])


def favorable_group_summary(cases: pd.DataFrame) -> pd.DataFrame:
    """Compare feature distributions for ratio > 1 and ratio <= 1 cases."""
    usable = cases[
        (cases["structure_status"] == "success") & cases["ratio_valid"]
    ].copy()
    rows: list[dict[str, Any]] = []
    for scope, frame in (
        ("full_hgbs5_sweep", usable),
        (
            "matched_BH_BeH_LiH",
            usable[usable["molecule"].isin(MATCHED_MOLECULES)],
        ),
    ):
        favorable_mask = frame["jw_magnitude_to_fermionic_ratio"] > 1.0
        for feature in PRIMARY_FEATURES:
            favorable = frame.loc[favorable_mask, feature].dropna().to_numpy(float)
            unfavorable = frame.loc[~favorable_mask, feature].dropna().to_numpy(float)
            direction = PREDICTED_DIRECTIONS[feature]
            statistic = p_value = math.nan
            if favorable.size and unfavorable.size:
                test = stats.mannwhitneyu(
                    direction * favorable,
                    direction * unfavorable,
                    alternative="greater",
                )
                statistic = float(test.statistic)
                p_value = float(test.pvalue)
            rows.append(
                {
                    "scope": scope,
                    "feature": feature,
                    "feature_label": FEATURE_LABELS[feature],
                    "predicted_favorable_direction": (
                        "higher" if direction > 0 else "lower"
                    ),
                    "favorable_cases": favorable.size,
                    "unfavorable_or_tied_cases": unfavorable.size,
                    "favorable_median": (
                        float(np.median(favorable)) if favorable.size else math.nan
                    ),
                    "unfavorable_or_tied_median": (
                        float(np.median(unfavorable))
                        if unfavorable.size
                        else math.nan
                    ),
                    "direction_aligned_median_difference": direction
                    * (
                        float(np.median(favorable))
                        - float(np.median(unfavorable))
                    ),
                    "mann_whitney_u": statistic,
                    "one_sided_p_value": p_value,
                }
            )
    table = pd.DataFrame(rows)
    table["fdr_q_value"] = table.groupby("scope", group_keys=False)[
        "one_sided_p_value"
    ].apply(benjamini_hochberg)
    return table


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.2,
        }
    )


def make_plots(cases: pd.DataFrame, correlations: pd.DataFrame, outdir: Path, dpi: int) -> None:
    configure_plot_style()
    valid = cases[
        (cases["structure_status"] == "success") & cases["ratio_valid"]
    ]
    measured = cases[
        cases["ratio_valid"]
        & np.isfinite(cases["signed_bch_cancellation_ratio"])
        & np.isfinite(cases["log10_jw_magnitude_to_fermionic_ratio"])
    ]

    full_raw = correlations[
        (correlations["scope"] == "full_hgbs5_sweep")
        & (correlations["fixed_effects"] == "none")
        & (
            correlations["outcome"]
            == "log10_jw_magnitude_to_fermionic_ratio"
        )
    ].set_index("feature")
    plot_feature = full_raw["direction_aligned_pearson_r"].idxmax()

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    for molecule, group in valid.groupby("molecule", sort=True):
        axes[0].scatter(
            group[plot_feature],
            group["log10_jw_magnitude_to_fermionic_ratio"],
            s=24,
            alpha=0.75,
            label=molecule,
        )
    axes[0].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[0].set_xlabel(FEATURE_LABELS[plot_feature])
    axes[0].set_ylabel(r"$\log_{10}(\epsilon_{JW_{mag}}/\epsilon_{fermionic})$")
    axes[0].set_title("Full HGBS-5 sweep")

    axes[1].scatter(
        measured["signed_bch_cancellation_strength"],
        measured["log10_jw_magnitude_to_fermionic_ratio"],
        s=42,
        color="#5b4b9a",
    )
    annotation_offsets = {
        "B-B": (8, 14),
        "Be-Be": (-10, 14),
        "BeH2": (-44, -14),
        "F-F": (4, 8),
        "H2O": (4, 6),
        "Li-Li": (7, 13),
        "NH3": (4, 6),
        "O-O": (12, 8),
    }
    for row in measured.itertuples():
        offset = annotation_offsets.get(row.molecule, (3, 3))
        axes[1].annotate(
            row.molecule,
            (row.signed_bch_cancellation_strength, row.log10_jw_magnitude_to_fermionic_ratio),
            xytext=offset,
            textcoords="offset points",
            fontsize=7,
        )
    axes[1].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel(r"Measured cancellation strength $-\log_{10}(R_{BCH})$")
    axes[1].set_ylabel(r"$\log_{10}(\epsilon_{JW_{mag}}/\epsilon_{fermionic})$")
    axes[1].set_title("Independent BCH subset")

    heat = full_raw.loc[list(PRIMARY_FEATURES), "direction_aligned_pearson_r"].to_numpy()[:, None]
    image = axes[2].imshow(heat, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
    axes[2].set_xticks([0], ["direction-aligned r"])
    axes[2].set_yticks(
        range(len(PRIMARY_FEATURES)),
        [FEATURE_LABELS[name] for name in PRIMARY_FEATURES],
        fontsize=7,
    )
    for row_index, value in enumerate(heat[:, 0]):
        axes[2].text(0, row_index, f"{value:.2f}", ha="center", va="center", fontsize=7)
    axes[2].set_title("Predicted-direction correlations")
    fig.colorbar(image, ax=axes[2], fraction=0.06, pad=0.04)
    if len(valid["molecule"].unique()) <= 15:
        axes[0].legend(ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "signed_order_cancellation_diagnostics.png", dpi=dpi)
    fig.savefig(outdir / "signed_order_cancellation_diagnostics.pdf")
    plt.close(fig)


def format_stat(row: pd.Series) -> str:
    return (
        f"r={row['pearson_r']:.3f}, p={row['pearson_p_value']:.3g}, "
        f"q={row['pearson_fdr_q_value']:.3g}, n={int(row['cases'])}"
    )


def write_report(
    cases: pd.DataFrame,
    correlations: pd.DataFrame,
    cancellation_correlations: pd.DataFrame,
    outdir: Path,
) -> None:
    success = cases[cases["structure_status"] == "success"]
    valid = success[success["ratio_valid"]]
    unavailable = cases[cases["structure_status"] == "unavailable"]
    measured = cases[np.isfinite(cases["signed_bch_cancellation_ratio"])]
    full_raw = correlations[
        (correlations["scope"] == "full_hgbs5_sweep")
        & (correlations["fixed_effects"] == "none")
        & (
            correlations["outcome"]
            == "log10_jw_magnitude_to_fermionic_ratio"
        )
    ].sort_values("direction_aligned_pearson_r", ascending=False)
    full_adjusted = correlations[
        (correlations["scope"] == "full_hgbs5_sweep")
        & (correlations["fixed_effects"] == "molecule+n_qubits")
        & (
            correlations["outcome"]
            == "log10_jw_magnitude_to_fermionic_ratio"
        )
    ].sort_values("direction_aligned_pearson_r", ascending=False)
    infidelity_adjusted = correlations[
        (correlations["scope"] == "full_hgbs5_sweep")
        & (correlations["fixed_effects"] == "molecule+n_qubits")
        & (correlations["outcome"] == "log10_fermionic_infidelity")
    ].sort_values("direction_aligned_pearson_r", ascending=False)
    matched_size = correlations[
        (correlations["scope"] == "matched_BH_BeH_LiH")
        & (correlations["fixed_effects"] == "n_qubits")
        & (
            correlations["outcome"]
            == "log10_jw_magnitude_to_fermionic_ratio"
        )
    ].sort_values("direction_aligned_pearson_r", ascending=False)
    direct_cancel = cancellation_correlations[
        (
            cancellation_correlations["feature"]
            == "signed_bch_cancellation_strength"
        )
        & (
            cancellation_correlations["outcome"]
            == "log10_jw_magnitude_to_fermionic_ratio"
        )
    ].iloc[0]
    direct_cancel_infidelity = cancellation_correlations[
        (
            cancellation_correlations["feature"]
            == "signed_bch_cancellation_strength"
        )
        & (
            cancellation_correlations["outcome"]
            == "log10_fermionic_infidelity"
        )
    ].iloc[0]
    cancellation_structure_rows = cancellation_correlations[
        cancellation_correlations["outcome"]
        == "signed_bch_cancellation_strength"
    ].sort_values("direction_aligned_pearson_r", ascending=False)
    best_cancellation_structure = cancellation_structure_rows.iloc[0]
    order_sensitive = set(PRIMARY_FEATURES) - {
        "opposite_sign_interaction_weight_fraction"
    }
    adjusted_order_support = full_adjusted[
        full_adjusted["feature"].isin(order_sensitive)
        & (full_adjusted["direction_aligned_pearson_r"] > 0.0)
        & (full_adjusted["pearson_fdr_q_value"] < 0.05)
    ]
    cancellation_structure = cancellation_correlations[
        cancellation_correlations["feature"].isin(order_sensitive)
        & cancellation_correlations["outcome"].isin(
            {
                "signed_bch_cancellation_strength",
                "log10_signed_bch_cancellation_advantage_vs_fermionic_magnitude",
            }
        )
        & (cancellation_correlations["direction_aligned_pearson_r"] > 0.0)
        & (cancellation_correlations["pearson_fdr_q_value"] < 0.05)
    ]
    direct_cancellation_support = (
        direct_cancel["pearson_r"] > 0.0
        and direct_cancel["pearson_p_value"] < 0.05
    )
    strict_support = (
        not adjusted_order_support.empty
        and not cancellation_structure.empty
        and direct_cancellation_support
    )
    verdict = "verified" if strict_support else "not verified"

    lines = [
        "# Signed-order cancellation-structure test",
        "",
        "## Result",
        "",
        f"The revised cancellation-structure hypothesis is **{verdict} as "
        "stated**. The independent BCH result supports its core cancellation "
        "premise, but the proposed order-sensitive shared-support metrics do "
        "not consistently identify that cancellation structure across the "
        "sweep.",
        "",
        "## Definition",
        "",
        "Nodes are complete Hermitian fermionic parents. For parents `i,j`, "
        "the support-overlap interaction weight is "
        "`w_ij = |c_i c_j| * |S_i intersect S_j|`. This coefficient-weighted "
        "shared-support graph is a pre-BCH proxy, not an exact commutator graph. "
        "The signed order is ascending canonical signed coefficient; the control "
        "order is descending coefficient magnitude.",
        "",
        "The primary metrics include the opposite-sign share of interaction "
        "weight, coefficient-weighted mean edge span, weight within 5% and 10% "
        "of the parent sequence, local opposite-sign weight/enrichment, and "
        "signed-minus-magnitude order improvements. Span is normalized by `N-1`.",
        "",
        "## Coverage",
        "",
        f"- Input HGBS-5 cases: {len(cases)}",
        f"- Parent structures built: {len(success)}",
        f"- Valid performance comparisons: {len(valid)}",
        f"- Independent BCH cancellation cases: {len(measured)}",
        f"- Unavailable source Hamiltonians: {len(unavailable)}",
        "",
        "## Performance correlations",
        "",
        "The response is `log10(epsilon_JW-magnitude / epsilon_fermionic)`, so "
        "positive values favor the fermionic-aware ordering. Reported `r` values "
        "below retain each metric's natural sign; the ranking uses the "
        "predeclared favorable direction (lower span, higher other metrics).",
        "",
        "### Full sweep, raw",
        "",
    ]
    for _, row in full_raw.head(5).iterrows():
        lines.append(f"- {row['feature_label']}: {format_stat(row)}")
    lines.extend(["", "### Full sweep, molecule + active-size adjusted", ""])
    for _, row in full_adjusted.head(5).iterrows():
        lines.append(f"- {row['feature_label']}: {format_stat(row)}")
    lines.extend(["", "### Matched B-H / Be-H / Li-H, active-size adjusted", ""])
    for _, row in matched_size.head(5).iterrows():
        lines.append(f"- {row['feature_label']}: {format_stat(row)}")
    lines.extend(
        [
            "",
            "### Absolute fermionic infidelity, molecule + active-size adjusted",
            "",
        ]
    )
    for _, row in infidelity_adjusted.head(5).iterrows():
        lines.append(f"- {row['feature_label']}: {format_stat(row)}")

    lines.extend(
        [
            "",
            "## Independent cancellation check",
            "",
            "For the eight ablation cases, `R_BCH = ||sum_k v_k|| / "
            "sum_k ||v_k||` is the independently measured HF-state BCH2 "
            "cancellation ratio. Smaller `R_BCH` means more destructive "
            "cancellation; the analysis uses `-log10(R_BCH)` as cancellation "
            "strength.",
            "",
            "- Measured cancellation strength vs performance advantage: "
            f"{format_stat(direct_cancel)}",
            "- Measured cancellation strength vs absolute fermionic infidelity: "
            f"{format_stat(direct_cancel_infidelity)}",
            "",
            "None of the nine structural proxies significantly predicts the "
            "measured cancellation outcomes after FDR correction. The strongest "
            "suggestive relation is "
            f"{best_cancellation_structure['feature_label'].lower()} vs absolute "
            "cancellation strength ("
            f"{format_stat(best_cancellation_structure)}), but it is not conclusive.",
            "",
            "Several full-sweep locality results point opposite to the proposed "
            "direction after molecule and active-size adjustment: larger signed "
            "edge span and smaller local interaction weight accompany more "
            "relative advantage. Those same locality metrics do track lower "
            "absolute fermionic infidelity, but not the relative advantage or "
            "measured cancellation consistently. Therefore the successful BCH "
            "cancellation cannot be reduced to these parent-support span/locality "
            "statistics.",
            "",
            "The detailed CSV also tests every structural feature against both "
            "absolute signed-order cancellation strength and the signed-order "
            "cancellation advantage over fermionic magnitude ordering.",
            "",
            "## Interpretation guardrails",
            "",
            "- Shared orbital support is necessary for many interactions but is "
            "not an exact commutator magnitude or sign.",
            "- The independent cancellation subset contains eight deliberately "
            "selected cases, so its correlations are mechanism checks rather "
            "than population estimates.",
            "- Nine related structural metrics are tested; FDR-adjusted q-values "
            "are reported within each scope/model.",
            "- Cases missing both active tensors and reconstructable molecular "
            "pickles remain unavailable and are not imputed.",
            "",
            "## Files",
            "",
            "- `signed_order_case_metrics.csv`: case-level structure, outcomes, "
            "and cancellation joins.",
            "- `signed_order_performance_correlations.csv`: matched and full-sweep "
            "performance tests.",
            "- `signed_order_cancellation_correlations.csv`: independent BCH tests.",
            "- `matched_heteronuclear_by_size.csv`: B-H/Be-H/Li-H cases grouped "
            "by active-space size.",
            "- `favorable_structure_group_comparison.csv`: favorable-vs-other "
            "feature medians and one-sided rank tests.",
            "- `signed_order_cancellation_diagnostics.{png,pdf}`: diagnostic plot.",
            "",
        ]
    )
    (outdir / "signed_order_cancellation_report.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> None:
    args = parse_arguments()
    cancellation_inputs = (
        tuple(args.cancellation_inputs)
        if args.cancellation_inputs
        else DEFAULT_CANCELLATION_INPUTS
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    cases = read_cases(args.input, args.basis)
    metrics = analyze_cases(cases, args.tolerance)
    cancellation = load_cancellation_measurements(cancellation_inputs)
    joined = metrics.merge(cancellation, on="case_id", how="left")
    correlations = build_correlation_table(joined)
    cancellation_correlations = build_cancellation_correlation_table(joined)
    group_summary = favorable_group_summary(joined)

    joined.to_csv(args.outdir / "signed_order_case_metrics.csv", index=False)
    correlations.to_csv(
        args.outdir / "signed_order_performance_correlations.csv",
        index=False,
    )
    cancellation_correlations.to_csv(
        args.outdir / "signed_order_cancellation_correlations.csv",
        index=False,
    )
    matched_size_summary(joined).to_csv(
        args.outdir / "matched_heteronuclear_by_size.csv",
        index=False,
    )
    group_summary.to_csv(
        args.outdir / "favorable_structure_group_comparison.csv",
        index=False,
    )
    make_plots(joined, correlations, args.outdir, args.dpi)
    write_report(joined, correlations, cancellation_correlations, args.outdir)
    print(f"Wrote signed-order cancellation analysis to {args.outdir}")


if __name__ == "__main__":
    main()
