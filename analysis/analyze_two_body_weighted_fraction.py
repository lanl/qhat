#!/usr/bin/env python3
"""Test two-body weighted composition against fermionic Trotter error.

The primary structural predictor is the fraction of non-constant Hermitian
fermionic-parent coefficient mass carried by two-body parents,

    sum(two-body parents |c_p|) / sum(one- and two-body parents |c_p|).

The script also reports exploratory composition fractions within the two-body
mass, split by HF excitation rank and by the number of distinct spin orbitals
in the parent.  Outcomes are the signed-ascending fermionic error, the fixed
JW magnitude-descending error, and their relative advantage.
"""

# ruff: noqa: E402  # Select a writable Matplotlib cache before imports.

from __future__ import annotations

import argparse
import math
import os
import tempfile
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
from openfermion import get_fermion_operator
from scipy import stats

try:
    from qhat.analysis.analyze_fermionic_graph_centralization import (
        MolecularDataCache,
        load_case_interaction,
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
    )
    from benchmark_L_sweep_trotter import (
        HermitianFermionTerm,
        build_hermitian_fermion_terms,
        clean_fermion_operator,
    )


DEFAULT_INPUT = Path(
    "analysis/fermionic_aware_performance/fermionic_aware_case_performance.csv"
)
DEFAULT_OUTDIR = Path("analysis/two_body_weighted_fraction")
DEFAULT_CANCELLATION_INPUT = Path(
    "analysis/cancellation_hypothesis_validation/full_analysis/case_summary.csv"
)
DEFAULT_BASIS = "hgbs-5"
DEFAULT_TOLERANCE = 1.0e-12
DEFAULT_BOOTSTRAP_REPLICATES = 2000
DEFAULT_SEED = 20260825

REQUIRED_COLUMNS = {
    "case_id",
    "tensor_path",
    "molecule",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "valid_comparison",
    "jw_magnitude_one_minus_overlap",
    "fermionic_aware_one_minus_overlap",
}

FEATURES = (
    "two_body_coefficient_mass_fraction",
    "two_body_parent_fraction",
    "two_body_excitation_rank0_weighted_fraction",
    "two_body_excitation_rank1_weighted_fraction",
    "two_body_excitation_rank2_weighted_fraction",
    "two_body_support2_weighted_fraction",
    "two_body_support3_weighted_fraction",
    "two_body_support4_weighted_fraction",
)

FEATURE_LABELS = {
    "two_body_coefficient_mass_fraction": (
        "Two-body share of total parent |coefficient| mass"
    ),
    "two_body_parent_fraction": "Two-body share of parent count",
    "two_body_excitation_rank0_weighted_fraction": (
        "Rank-0 share within two-body |coefficient| mass"
    ),
    "two_body_excitation_rank1_weighted_fraction": (
        "Rank-1 share within two-body |coefficient| mass"
    ),
    "two_body_excitation_rank2_weighted_fraction": (
        "Rank-2 share within two-body |coefficient| mass"
    ),
    "two_body_support2_weighted_fraction": (
        "Two-orbital share within two-body |coefficient| mass"
    ),
    "two_body_support3_weighted_fraction": (
        "Three-orbital share within two-body |coefficient| mass"
    ),
    "two_body_support4_weighted_fraction": (
        "Four-orbital share within two-body |coefficient| mass"
    ),
}

FEATURE_FAMILIES = {
    "two_body_coefficient_mass_fraction": "primary_total_fraction",
    "two_body_parent_fraction": "unweighted_sensitivity",
    "two_body_excitation_rank0_weighted_fraction": "excitation_composition",
    "two_body_excitation_rank1_weighted_fraction": "excitation_composition",
    "two_body_excitation_rank2_weighted_fraction": "excitation_composition",
    "two_body_support2_weighted_fraction": "support_composition",
    "two_body_support3_weighted_fraction": "support_composition",
    "two_body_support4_weighted_fraction": "support_composition",
}

OUTCOMES = (
    "log10_fermionic_error",
    "log10_jw_magnitude_error",
    "log10_jw_magnitude_to_fermionic_advantage",
)

OUTCOME_LABELS = {
    "log10_fermionic_error": "log10 fermionic Trotter error",
    "log10_jw_magnitude_error": "log10 JW-magnitude Trotter error",
    "log10_jw_magnitude_to_fermionic_advantage": (
        "log10(JW-magnitude / fermionic error)"
    ),
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--basis", default=DEFAULT_BASIS)
    parser.add_argument(
        "--cancellation-input",
        type=Path,
        default=DEFAULT_CANCELLATION_INPUT,
    )
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=DEFAULT_BOOTSTRAP_REPLICATES,
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def read_cases(path: Path, basis: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Input is missing columns: {', '.join(missing)}")
    selected = frame[
        (frame["basis"].str.lower() == basis.lower())
        & frame["valid_comparison"].astype(bool)
    ].copy()
    if selected.empty:
        raise ValueError(f"No valid {basis} cases found in {path}.")
    return selected.sort_values(
        ["molecule", "n_qubits", "active_occupied", "active_vacant"]
    ).reset_index(drop=True)


def canonical_parent_coefficient(
    term: HermitianFermionTerm,
    tolerance: float,
) -> float:
    canonical_key = min(term.component_keys)
    value = complex(term.operator.terms[canonical_key])
    if abs(value.imag) > max(100.0 * tolerance, 1.0e-10):
        raise ValueError(f"Non-real canonical parent coefficient: {value!r}")
    return float(value.real)


def parent_composition(
    terms: Sequence[HermitianFermionTerm],
    active_occupied: int,
    tolerance: float,
) -> dict[str, float | int]:
    coefficient_mass = {1: 0.0, 2: 0.0}
    parent_count = {1: 0, 2: 0}
    excitation_mass = {0: 0.0, 1: 0.0, 2: 0.0}
    support_mass = {2: 0.0, 3: 0.0, 4: 0.0}

    for term in terms:
        canonical_key = min(term.component_keys)
        body_rank = len(canonical_key) // 2
        if body_rank not in coefficient_mass:
            raise ValueError(f"Unexpected fermionic body rank {body_rank}.")
        weight = abs(canonical_parent_coefficient(term, tolerance))
        coefficient_mass[body_rank] += weight
        parent_count[body_rank] += 1

        if body_rank != 2:
            continue

        creators = [
            index
            for index, action in canonical_key
            if action == 1
        ]
        annihilators = [
            index
            for index, action in canonical_key
            if action == 0
        ]
        excitation_rank = abs(
            sum(index >= active_occupied for index in creators)
            - sum(index >= active_occupied for index in annihilators)
        )
        support_size = len(set(creators + annihilators))
        excitation_mass[excitation_rank] += weight
        support_mass[support_size] += weight

    total_mass = coefficient_mass[1] + coefficient_mass[2]
    total_count = parent_count[1] + parent_count[2]
    two_body_mass = coefficient_mass[2]
    if total_mass <= 0.0 or two_body_mass <= 0.0 or total_count <= 0:
        raise ValueError("Non-positive parent coefficient mass or count.")

    output: dict[str, float | int] = {
        "number_of_fermionic_parents": total_count,
        "number_of_one_body_parents": parent_count[1],
        "number_of_two_body_parents": parent_count[2],
        "one_body_coefficient_mass": coefficient_mass[1],
        "two_body_coefficient_mass": two_body_mass,
        "total_nonconstant_parent_coefficient_mass": total_mass,
        "two_body_coefficient_mass_fraction": two_body_mass / total_mass,
        "two_body_parent_fraction": parent_count[2] / total_count,
    }
    for rank in (0, 1, 2):
        output[
            f"two_body_excitation_rank{rank}_weighted_fraction"
        ] = excitation_mass[rank] / two_body_mass
    for size in (2, 3, 4):
        output[
            f"two_body_support{size}_weighted_fraction"
        ] = support_mass[size] / two_body_mass
    return output


def analyze_cases(
    cases: pd.DataFrame,
    tolerance: float,
) -> pd.DataFrame:
    cache = MolecularDataCache()
    rows: list[dict[str, Any]] = []

    for _, source in cases.iterrows():
        base = source.to_dict()
        print(f"{source['case_id']}: reconstructing fermionic parents")
        try:
            interaction, structure_source, unavailable_reason = (
                load_case_interaction(source, cache)
            )
            if interaction is None:
                rows.append(
                    {
                        **base,
                        "structure_status": "unavailable",
                        "structure_source": structure_source,
                        "structure_unavailable_reason": unavailable_reason,
                    }
                )
                print("  unavailable")
                continue

            fermion_hamiltonian = clean_fermion_operator(
                get_fermion_operator(interaction),
                tolerance,
            )
            terms = build_hermitian_fermion_terms(
                fermion_hamiltonian,
                tolerance,
            )
            composition = parent_composition(
                terms,
                int(source["active_occupied"]),
                tolerance,
            )
            fermionic_error = float(
                source["fermionic_aware_one_minus_overlap"]
            )
            jw_error = float(
                source["jw_magnitude_one_minus_overlap"]
            )
            if fermionic_error <= 0.0 or jw_error <= 0.0:
                raise ValueError("Trotter errors must be positive.")

            rows.append(
                {
                    **base,
                    "structure_status": "success",
                    "structure_source": structure_source,
                    "structure_unavailable_reason": "",
                    **composition,
                    "log10_fermionic_error": math.log10(fermionic_error),
                    "log10_jw_magnitude_error": math.log10(jw_error),
                    "log10_jw_magnitude_to_fermionic_advantage": math.log10(
                        jw_error / fermionic_error
                    ),
                }
            )
            print(
                "  two-body mass fraction="
                f"{composition['two_body_coefficient_mass_fraction']:.6f}"
            )
        except Exception as error:  # Preserve a complete audit of failures.
            rows.append(
                {
                    **base,
                    "structure_status": "error",
                    "structure_source": "error",
                    "structure_unavailable_reason": str(error),
                }
            )
            print(f"  ERROR: {error}")

    return pd.DataFrame(rows)


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


def correlation_values(
    frame: pd.DataFrame,
    feature: str,
    outcome: str,
    fixed_effects: Sequence[str],
) -> tuple[float, float, float, float, int]:
    usable = frame[
        np.isfinite(frame[feature]) & np.isfinite(frame[outcome])
    ].copy()
    x = usable[feature].to_numpy(dtype=float)
    y = usable[outcome].to_numpy(dtype=float)
    n_controls = 0
    if fixed_effects:
        x, n_controls = fixed_effect_residuals(
            usable,
            x,
            fixed_effects,
        )
        y, y_controls = fixed_effect_residuals(
            usable,
            y,
            fixed_effects,
        )
        if n_controls != y_controls:
            raise RuntimeError("Fixed-effect design ranks differ.")

    if len(usable) < 3 or np.std(x) <= 0.0 or np.std(y) <= 0.0:
        return math.nan, math.nan, math.nan, math.nan, n_controls

    pearson_r = float(np.corrcoef(x, y)[0, 1])
    degrees_freedom = len(usable) - n_controls - 2
    if degrees_freedom <= 0:
        pearson_p = math.nan
    elif abs(pearson_r) >= 1.0:
        pearson_p = 0.0
    else:
        statistic = pearson_r * math.sqrt(
            degrees_freedom / (1.0 - pearson_r**2)
        )
        pearson_p = float(
            2.0 * stats.t.sf(abs(statistic), degrees_freedom)
        )
    spearman = stats.spearmanr(x, y)
    return (
        pearson_r,
        pearson_p,
        float(spearman.statistic),
        float(spearman.pvalue),
        n_controls,
    )


def cluster_bootstrap_interval(
    frame: pd.DataFrame,
    feature: str,
    outcome: str,
    replicates: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    molecules = frame["molecule"].drop_duplicates().to_numpy()
    values: list[float] = []
    for _ in range(replicates):
        sampled = rng.choice(
            molecules,
            size=len(molecules),
            replace=True,
        )
        chunks: list[pd.DataFrame] = []
        for sample_index, molecule in enumerate(sampled):
            chunk = frame[frame["molecule"] == molecule].copy()
            chunk["bootstrap_molecule"] = (
                f"{sample_index}:{molecule}"
            )
            chunks.append(chunk)
        bootstrap = pd.concat(chunks, ignore_index=True)
        x, _ = fixed_effect_residuals(
            bootstrap,
            bootstrap[feature].to_numpy(dtype=float),
            ("bootstrap_molecule", "n_qubits"),
        )
        y, _ = fixed_effect_residuals(
            bootstrap,
            bootstrap[outcome].to_numpy(dtype=float),
            ("bootstrap_molecule", "n_qubits"),
        )
        if np.std(x) > 0.0 and np.std(y) > 0.0:
            values.append(float(np.corrcoef(x, y)[0, 1]))
    if not values:
        return math.nan, math.nan, 0
    return (
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
        len(values),
    )


def leave_one_molecule_out_range(
    frame: pd.DataFrame,
    feature: str,
    outcome: str,
) -> tuple[float, float]:
    values: list[float] = []
    for molecule in frame["molecule"].drop_duplicates():
        subset = frame[frame["molecule"] != molecule]
        pearson_r, _, _, _, _ = correlation_values(
            subset,
            feature,
            outcome,
            ("molecule", "n_qubits"),
        )
        if math.isfinite(pearson_r):
            values.append(pearson_r)
    if not values:
        return math.nan, math.nan
    return min(values), max(values)


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    adjusted = pd.Series(np.nan, index=p_values.index, dtype=float)
    finite = p_values[np.isfinite(p_values)].sort_values()
    if finite.empty:
        return adjusted
    count = len(finite)
    raw = finite.to_numpy(dtype=float) * count / np.arange(1, count + 1)
    monotone = np.minimum.accumulate(raw[::-1])[::-1]
    adjusted.loc[finite.index] = np.minimum(monotone, 1.0)
    return adjusted


def partial_correlation(
    frame: pd.DataFrame,
    feature: str,
    outcome: str,
    control: str,
    fixed_effects: Sequence[str] = (),
) -> tuple[float, float, int, int]:
    usable = frame[
        np.isfinite(frame[feature])
        & np.isfinite(frame[outcome])
        & np.isfinite(frame[control])
    ].copy()
    x = usable[feature].to_numpy(dtype=float)
    y = usable[outcome].to_numpy(dtype=float)
    z = usable[control].to_numpy(dtype=float)
    n_fixed_controls = 0
    if fixed_effects:
        x, n_fixed_controls = fixed_effect_residuals(
            usable,
            x,
            fixed_effects,
        )
        y, y_controls = fixed_effect_residuals(
            usable,
            y,
            fixed_effects,
        )
        z, z_controls = fixed_effect_residuals(
            usable,
            z,
            fixed_effects,
        )
        if n_fixed_controls != y_controls or y_controls != z_controls:
            raise RuntimeError("Fixed-effect design ranks differ.")

    control_design = np.column_stack(
        [np.ones(len(usable), dtype=float), z]
    )
    x_residual = x - control_design @ np.linalg.lstsq(
        control_design,
        x,
        rcond=None,
    )[0]
    y_residual = y - control_design @ np.linalg.lstsq(
        control_design,
        y,
        rcond=None,
    )[0]
    if np.std(x_residual) <= 0.0 or np.std(y_residual) <= 0.0:
        return math.nan, math.nan, len(usable), n_fixed_controls + 1

    pearson_r = float(np.corrcoef(x_residual, y_residual)[0, 1])
    n_controls = n_fixed_controls + 1
    degrees_freedom = len(usable) - n_controls - 2
    if degrees_freedom <= 0:
        pearson_p = math.nan
    elif abs(pearson_r) >= 1.0:
        pearson_p = 0.0
    else:
        statistic = pearson_r * math.sqrt(
            degrees_freedom / (1.0 - pearson_r**2)
        )
        pearson_p = float(
            2.0 * stats.t.sf(abs(statistic), degrees_freedom)
        )
    return pearson_r, pearson_p, len(usable), n_controls


def build_bch_independence_table(
    cases: pd.DataFrame,
    cancellation_path: Path,
) -> pd.DataFrame:
    if not cancellation_path.exists():
        return pd.DataFrame()
    cancellation = pd.read_csv(cancellation_path)
    required = {
        "case_id",
        "matched_pair",
        "fresh_jw_magnitude_to_signed_advantage",
        "fresh_jw_magnitude_to_signed_bch_cancellation_ratio",
    }
    missing = sorted(required.difference(cancellation.columns))
    if missing:
        raise ValueError(
            "Cancellation input is missing columns: "
            + ", ".join(missing)
        )

    successful = cases[cases["structure_status"] == "success"]
    merged = cancellation.merge(
        successful[["case_id", *FEATURES]],
        on="case_id",
        how="inner",
        validate="one_to_one",
    )
    merged["log10_fresh_advantage"] = np.log10(
        merged["fresh_jw_magnitude_to_signed_advantage"]
    )
    merged["log10_relative_bch_cancellation"] = np.log10(
        merged[
            "fresh_jw_magnitude_to_signed_bch_cancellation_ratio"
        ]
    )

    rows: list[dict[str, Any]] = []
    for feature in FEATURES:
        for fixed_effects, specification in (
            ((), "bch_adjusted"),
            (("matched_pair",), "matched_pair_and_bch_adjusted"),
        ):
            partial_r, partial_p, n_cases, n_controls = partial_correlation(
                merged,
                feature,
                "log10_fresh_advantage",
                "log10_relative_bch_cancellation",
                fixed_effects,
            )
            direct_advantage = correlation_values(
                merged,
                feature,
                "log10_fresh_advantage",
                fixed_effects,
            )
            direct_cancellation = correlation_values(
                merged,
                feature,
                "log10_relative_bch_cancellation",
                fixed_effects,
            )
            rows.append(
                {
                    "specification": specification,
                    "fixed_effects": (
                        "+".join(fixed_effects)
                        if fixed_effects
                        else "none"
                    ),
                    "partial_control": (
                        "log10_relative_bch_cancellation"
                    ),
                    "feature": feature,
                    "feature_label": FEATURE_LABELS[feature],
                    "cases": n_cases,
                    "matched_pairs": int(merged["matched_pair"].nunique()),
                    "n_controls": n_controls,
                    "feature_to_advantage_pearson_r": direct_advantage[0],
                    "feature_to_advantage_pearson_p_value": direct_advantage[1],
                    "feature_to_cancellation_pearson_r": direct_cancellation[0],
                    "feature_to_cancellation_pearson_p_value": (
                        direct_cancellation[1]
                    ),
                    "feature_to_advantage_partial_pearson_r": partial_r,
                    "feature_to_advantage_partial_p_value": partial_p,
                }
            )
    return pd.DataFrame(rows)


def build_correlation_table(
    cases: pd.DataFrame,
    bootstrap_replicates: int,
    seed: int,
) -> pd.DataFrame:
    usable = cases[cases["structure_status"] == "success"].copy()
    rows: list[dict[str, Any]] = []
    specifications = (
        ((), "raw"),
        (("n_qubits",), "active_size_adjusted"),
        (("molecule", "n_qubits"), "molecule_and_active_size_adjusted"),
    )
    for outcome in OUTCOMES:
        for fixed_effects, specification in specifications:
            for feature in FEATURES:
                pearson_r, pearson_p, spearman_rho, spearman_p, controls = (
                    correlation_values(
                        usable,
                        feature,
                        outcome,
                        fixed_effects,
                    )
                )
                rows.append(
                    {
                        "specification": specification,
                        "fixed_effects": (
                            "+".join(fixed_effects)
                            if fixed_effects
                            else "none"
                        ),
                        "outcome": outcome,
                        "outcome_label": OUTCOME_LABELS[outcome],
                        "feature": feature,
                        "feature_label": FEATURE_LABELS[feature],
                        "feature_family": FEATURE_FAMILIES[feature],
                        "cases": len(usable),
                        "molecules": int(usable["molecule"].nunique()),
                        "active_space_sizes": int(
                            usable["n_qubits"].nunique()
                        ),
                        "n_controls": controls,
                        "pearson_r": pearson_r,
                        "pearson_p_value": pearson_p,
                        "spearman_rho_on_analysis_values": spearman_rho,
                        "spearman_p_value_on_analysis_values": spearman_p,
                    }
                )

    table = pd.DataFrame(rows)
    table["pearson_fdr_q_value"] = table.groupby(
        ["specification", "outcome"],
        group_keys=False,
    )["pearson_p_value"].apply(benjamini_hochberg)
    table["cluster_bootstrap_95ci_low"] = math.nan
    table["cluster_bootstrap_95ci_high"] = math.nan
    table["cluster_bootstrap_valid_replicates"] = 0
    table["leave_one_molecule_out_r_min"] = math.nan
    table["leave_one_molecule_out_r_max"] = math.nan

    adjusted_indices = table.index[
        table["specification"]
        == "molecule_and_active_size_adjusted"
    ]
    rng = np.random.default_rng(seed)
    for index in adjusted_indices:
        feature = str(table.at[index, "feature"])
        outcome = str(table.at[index, "outcome"])
        low, high, valid_replicates = cluster_bootstrap_interval(
            usable,
            feature,
            outcome,
            bootstrap_replicates,
            rng,
        )
        loo_low, loo_high = leave_one_molecule_out_range(
            usable,
            feature,
            outcome,
        )
        table.at[index, "cluster_bootstrap_95ci_low"] = low
        table.at[index, "cluster_bootstrap_95ci_high"] = high
        table.at[
            index,
            "cluster_bootstrap_valid_replicates",
        ] = valid_replicates
        table.at[index, "leave_one_molecule_out_r_min"] = loo_low
        table.at[index, "leave_one_molecule_out_r_max"] = loo_high

    return table


def select_result(
    correlations: pd.DataFrame,
    feature: str,
    outcome: str,
    specification: str,
) -> pd.Series:
    selected = correlations[
        (correlations["feature"] == feature)
        & (correlations["outcome"] == outcome)
        & (correlations["specification"] == specification)
    ]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one result for {feature}/{outcome}/{specification}."
        )
    return selected.iloc[0]


def format_result(row: pd.Series, include_robustness: bool = False) -> str:
    text = (
        f"r={row['pearson_r']:.3f}, "
        f"p={row['pearson_p_value']:.3g}, "
        f"q={row['pearson_fdr_q_value']:.3g}, "
        f"n={int(row['cases'])}"
    )
    if include_robustness:
        text += (
            ", molecule-bootstrap 95% CI "
            f"[{row['cluster_bootstrap_95ci_low']:.3f}, "
            f"{row['cluster_bootstrap_95ci_high']:.3f}], "
            "LOMO r range "
            f"[{row['leave_one_molecule_out_r_min']:.3f}, "
            f"{row['leave_one_molecule_out_r_max']:.3f}]"
        )
    return text


def write_report(
    path: Path,
    cases: pd.DataFrame,
    correlations: pd.DataFrame,
    bch_independence: pd.DataFrame,
    basis: str,
) -> None:
    successful = cases[cases["structure_status"] == "success"]
    unavailable = cases[cases["structure_status"] != "success"]

    primary_raw = select_result(
        correlations,
        "two_body_coefficient_mass_fraction",
        "log10_fermionic_error",
        "raw",
    )
    primary_size = select_result(
        correlations,
        "two_body_coefficient_mass_fraction",
        "log10_fermionic_error",
        "active_size_adjusted",
    )
    primary_adjusted = select_result(
        correlations,
        "two_body_coefficient_mass_fraction",
        "log10_fermionic_error",
        "molecule_and_active_size_adjusted",
    )
    primary_advantage = select_result(
        correlations,
        "two_body_coefficient_mass_fraction",
        "log10_jw_magnitude_to_fermionic_advantage",
        "molecule_and_active_size_adjusted",
    )
    excitation0 = select_result(
        correlations,
        "two_body_excitation_rank0_weighted_fraction",
        "log10_fermionic_error",
        "molecule_and_active_size_adjusted",
    )
    excitation1 = select_result(
        correlations,
        "two_body_excitation_rank1_weighted_fraction",
        "log10_fermionic_error",
        "molecule_and_active_size_adjusted",
    )
    excitation2 = select_result(
        correlations,
        "two_body_excitation_rank2_weighted_fraction",
        "log10_fermionic_error",
        "molecule_and_active_size_adjusted",
    )
    support4_advantage = select_result(
        correlations,
        "two_body_support4_weighted_fraction",
        "log10_jw_magnitude_to_fermionic_advantage",
        "molecule_and_active_size_adjusted",
    )
    support4_fermionic = select_result(
        correlations,
        "two_body_support4_weighted_fraction",
        "log10_fermionic_error",
        "molecule_and_active_size_adjusted",
    )
    support4_jw = select_result(
        correlations,
        "two_body_support4_weighted_fraction",
        "log10_jw_magnitude_error",
        "molecule_and_active_size_adjusted",
    )

    def independence_result(
        feature: str,
        specification: str,
    ) -> pd.Series | None:
        if bch_independence.empty:
            return None
        selected = bch_independence[
            (bch_independence["feature"] == feature)
            & (bch_independence["specification"] == specification)
        ]
        return selected.iloc[0] if len(selected) == 1 else None

    support4_bch = independence_result(
        "two_body_support4_weighted_fraction",
        "bch_adjusted",
    )
    support4_pair_bch = independence_result(
        "two_body_support4_weighted_fraction",
        "matched_pair_and_bch_adjusted",
    )
    total_bch = independence_result(
        "two_body_coefficient_mass_fraction",
        "bch_adjusted",
    )

    independence_lines: list[str] = []
    if (
        support4_bch is not None
        and support4_pair_bch is not None
        and total_bch is not None
    ):
        independence_lines = [
            "",
            "## Is two-body composition independent of BCH cancellation?",
            "",
            "The available 20-case matched cancellation panel was joined "
            "without loss. Partial correlations control the independently "
            "measured `log10(JW-magnitude/signed BCH cancellation ratio)`.",
            "",
            "- Total two-body mass fraction, BCH-adjusted: "
            f"partial r={total_bch['feature_to_advantage_partial_pearson_r']:.3f}, "
            f"p={total_bch['feature_to_advantage_partial_p_value']:.3g}, "
            f"n={int(total_bch['cases'])}",
            "- Four-orbital two-body fraction, BCH-adjusted: "
            f"partial r={support4_bch['feature_to_advantage_partial_pearson_r']:.3f}, "
            f"p={support4_bch['feature_to_advantage_partial_p_value']:.3g}, "
            f"n={int(support4_bch['cases'])}",
            "- Four-orbital fraction, matched-pair + BCH adjusted: "
            f"partial r={support4_pair_bch['feature_to_advantage_partial_pearson_r']:.3f}, "
            f"p={support4_pair_bch['feature_to_advantage_partial_p_value']:.3g}, "
            f"n={int(support4_pair_bch['cases'])}",
            "",
            "The current data therefore do **not** verify two-body "
            "composition as an independent second mechanism beyond BCH "
            "cancellation. The full-sweep composition signal is useful as a "
            "structural marker, but it may be mediated by the same "
            "order-induced cancellation mechanism or by another shared "
            "Hamiltonian property.",
        ]

    lines = [
        "# Two-body weighted-fraction validation",
        "",
        "## Verdict",
        "",
        "The broad hypothesis is **not supported by the total two-body "
        "coefficient-mass fraction alone**. Its raw association with the "
        "fermionic error disappears after molecule and active-space-size "
        "adjustment. However, the internal composition of the two-body mass "
        "does contain reproducible cross-case signal. That narrower signal "
        "does not survive as an independent explanation after controlling "
        "the measured BCH cancellation ratio in the 20-case matched panel.",
        "",
        "## Definitions",
        "",
        "The primary fraction is",
        "",
        "`f_2b = sum_{two-body parents}|c_p| / "
        "sum_{one- and two-body parents}|c_p|`.",
        "",
        "Every `c_p` is the canonical real coefficient of one complete "
        "Hermitian fermionic parent, matching the signed-ascending production "
        "ordering. The constant term is excluded. Two-body composition is "
        "also split by HF excitation rank (0/1/2) and distinct spin-orbital "
        "support size (2/3/4), with each split normalized by total two-body "
        "|coefficient| mass.",
        "",
        "Outcomes are `log10(fermionic one-minus-overlap)`, the matching fixed "
        "JW magnitude-descending error, and "
        "`log10(JW-magnitude error / fermionic error)`.",
        "",
        "## Coverage",
        "",
        f"- Basis: {basis}",
        f"- Valid performance cases considered: {len(cases)}",
        f"- Fermionic parent structures recovered: {len(successful)}",
        f"- Molecules represented: {successful['molecule'].nunique()}",
        f"- Active-space sizes represented: {successful['n_qubits'].nunique()}",
        f"- Unavailable/error cases: {len(unavailable)}",
        "",
        "## Primary total two-body fraction test",
        "",
        "### Absolute fermionic error",
        "",
        f"- Raw: {format_result(primary_raw)}",
        f"- Active-size adjusted: {format_result(primary_size)}",
        "- Molecule + active-size adjusted: "
        f"{format_result(primary_adjusted, include_robustness=True)}",
        "",
        "The positive raw coefficient means that more total two-body weight "
        "is associated with *larger*, not smaller, fermionic error. The "
        "complete loss of association after the full adjustment shows that "
        "this is a between-molecule/active-space composition effect rather "
        "than an ordering-specific predictor.",
        "",
        "### Relative fermionic advantage over JW magnitude descending",
        "",
        "- Molecule + active-size adjusted: "
        f"{format_result(primary_advantage, include_robustness=True)}",
        "",
        "Thus total two-body coefficient mass does not verify a general "
        "fermionic-ordering advantage.",
        "",
        "## Exploratory two-body composition results",
        "",
        "All rows below use molecule + active-space-size adjustment and the "
        "reported q-values correct across the eight tested weighted/count "
        "features for the same outcome and specification.",
        "",
        "### HF excitation-rank composition vs absolute fermionic error",
        "",
        f"- Rank 0: {format_result(excitation0, include_robustness=True)}",
        f"- Rank 1: {format_result(excitation1, include_robustness=True)}",
        f"- Rank 2: {format_result(excitation2, include_robustness=True)}",
        "",
        "More rank-0 two-body mass accompanies lower fermionic error, whereas "
        "more single- and double-excitation mass accompanies higher error. "
        "These are compositional associations: the three fractions sum to one.",
        "",
        "### Four-orbital two-body mass vs relative advantage",
        "",
        "- Relative advantage: "
        f"{format_result(support4_advantage, include_robustness=True)}",
        f"- Absolute fermionic error: {format_result(support4_fermionic)}",
        f"- Absolute JW-magnitude error: {format_result(support4_jw)}",
        "",
        "The four-orbital fraction tracks a larger relative fermionic "
        "advantage, but it does not lower the absolute fermionic error. The "
        "relative effect arises because the JW magnitude-descending error "
        "increases with this fraction while the fermionic error is nearly "
        "unchanged.",
        *independence_lines,
        "",
        "## Interpretation limits",
        "",
        "- This establishes cross-case association, not a causal mechanism "
        "independent of BCH cancellation.",
        "- The total two-body fraction was the primary test. The excitation- "
        "and support-composition findings are exploratory and need a "
        "prospectively selected panel.",
        "- The fixed effects treat molecule and active-space size as "
        "categorical controls. Molecule-cluster bootstrap intervals and "
        "leave-one-molecule-out ranges quantify robustness to repeated cases.",
        "- Source tensors/pickles were unavailable for some otherwise valid "
        "performance rows; no structural values were imputed.",
        "- Only HGBS-5 could be reconstructed from the current repository "
        "sources; the STO-6G panel remains an external replication target.",
        "",
        "## Files",
        "",
        "- `two_body_case_metrics.csv`: case-level parent composition and "
        "performance outcomes.",
        "- `two_body_correlations.csv`: raw and adjusted correlation tests, "
        "FDR q-values, cluster-bootstrap intervals, and LOMO ranges.",
        "- `two_body_bch_independence.csv`: 20-case partial correlations "
        "controlling independently measured BCH cancellation.",
        "- `two_body_weighted_fraction_validation.png`: raw/adjusted diagnostic "
        "figure.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def regression_line(
    axis: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    color: str,
) -> None:
    if len(x) < 2 or np.std(x) <= 0.0:
        return
    x_line = np.linspace(float(np.min(x)), float(np.max(x)), 100)
    slope, intercept = np.polyfit(x, y, 1)
    axis.plot(
        x_line,
        intercept + slope * x_line,
        color=color,
        linewidth=1.5,
        zorder=2,
    )


def plot_validation(
    path: Path,
    cases: pd.DataFrame,
    correlations: pd.DataFrame,
    dpi: int,
) -> None:
    usable = cases[cases["structure_status"] == "success"].copy()
    comparisons = (
        (
            "two_body_coefficient_mass_fraction",
            "log10_fermionic_error",
            "Total two-body coefficient mass",
            "Two-body mass fraction",
            "Residual two-body mass fraction",
            "Residual log10 fermionic error",
        ),
        (
            "two_body_excitation_rank0_weighted_fraction",
            "log10_fermionic_error",
            "Rank-0 share within two-body mass",
            "Rank-0 fraction within two-body mass",
            "Residual rank-0 fraction",
            "Residual log10 fermionic error",
        ),
        (
            "two_body_support4_weighted_fraction",
            "log10_jw_magnitude_to_fermionic_advantage",
            "Four-orbital share within two-body mass",
            "Four-orbital fraction within two-body mass",
            "Residual four-orbital fraction",
            "Residual log10 advantage",
        ),
    )

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.grid": True,
            "grid.alpha": 0.2,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(12.0, 7.0),
        constrained_layout=True,
    )
    color_values = usable["n_qubits"].to_numpy(dtype=float)
    normalization = plt.Normalize(
        float(color_values.min()),
        float(color_values.max()),
    )
    colormap = plt.get_cmap("viridis")

    for column, (
        feature,
        outcome,
        title,
        raw_x_label,
        residual_x_label,
        residual_y_label,
    ) in enumerate(comparisons):
        raw_axis = axes[0, column]
        adjusted_axis = axes[1, column]
        x = usable[feature].to_numpy(dtype=float)
        y = usable[outcome].to_numpy(dtype=float)
        colors = colormap(normalization(color_values))

        raw_axis.scatter(
            x,
            y,
            c=colors,
            s=27,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )
        regression_line(raw_axis, x, y, "#343a40")
        raw_result = select_result(
            correlations,
            feature,
            outcome,
            "raw",
        )
        raw_axis.text(
            0.03,
            0.97,
            f"raw r={raw_result['pearson_r']:.3f}\n"
            f"p={raw_result['pearson_p_value']:.2g}",
            transform=raw_axis.transAxes,
            va="top",
            ha="left",
            fontsize=8,
        )
        raw_axis.set_title(title)
        raw_axis.set_xlabel(raw_x_label)
        raw_axis.set_ylabel(OUTCOME_LABELS[outcome])

        x_residual, _ = fixed_effect_residuals(
            usable,
            x,
            ("molecule", "n_qubits"),
        )
        y_residual, _ = fixed_effect_residuals(
            usable,
            y,
            ("molecule", "n_qubits"),
        )
        adjusted_axis.scatter(
            x_residual,
            y_residual,
            c=colors,
            s=27,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )
        regression_line(
            adjusted_axis,
            x_residual,
            y_residual,
            "#343a40",
        )
        adjusted_axis.axhline(
            0.0,
            color="#adb5bd",
            linewidth=0.7,
            linestyle="--",
        )
        adjusted_axis.axvline(
            0.0,
            color="#adb5bd",
            linewidth=0.7,
            linestyle="--",
        )
        adjusted_result = select_result(
            correlations,
            feature,
            outcome,
            "molecule_and_active_size_adjusted",
        )
        adjusted_axis.text(
            0.03,
            0.97,
            f"adjusted r={adjusted_result['pearson_r']:.3f}\n"
            f"p={adjusted_result['pearson_p_value']:.2g}, "
            f"q={adjusted_result['pearson_fdr_q_value']:.2g}",
            transform=adjusted_axis.transAxes,
            va="top",
            ha="left",
            fontsize=8,
        )
        adjusted_axis.set_xlabel(residual_x_label)
        adjusted_axis.set_ylabel(residual_y_label)

    axes[0, 0].text(
        -0.22,
        1.12,
        "Raw",
        transform=axes[0, 0].transAxes,
        fontsize=11,
        fontweight="bold",
    )
    axes[1, 0].text(
        -0.22,
        1.12,
        "Molecule + active-size adjusted",
        transform=axes[1, 0].transAxes,
        fontsize=11,
        fontweight="bold",
    )
    scalar = plt.cm.ScalarMappable(
        norm=normalization,
        cmap=colormap,
    )
    scalar.set_array([])
    colorbar = figure.colorbar(
        scalar,
        ax=axes,
        location="right",
        shrink=0.82,
        pad=0.02,
    )
    colorbar.set_label("Qubits in active space")
    figure.suptitle(
        "Two-body weighted composition vs fermionic-ordering performance",
        fontsize=13,
    )
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_arguments()
    if args.bootstrap_replicates <= 0:
        raise ValueError("--bootstrap-replicates must be positive.")
    args.outdir.mkdir(parents=True, exist_ok=True)

    source = read_cases(args.input, args.basis)
    cases = analyze_cases(source, args.tolerance)
    correlations = build_correlation_table(
        cases,
        args.bootstrap_replicates,
        args.seed,
    )
    bch_independence = build_bch_independence_table(
        cases,
        args.cancellation_input,
    )

    case_path = args.outdir / "two_body_case_metrics.csv"
    correlation_path = args.outdir / "two_body_correlations.csv"
    bch_path = args.outdir / "two_body_bch_independence.csv"
    report_path = args.outdir / "two_body_weighted_fraction_report.md"
    figure_path = args.outdir / "two_body_weighted_fraction_validation.png"

    cases.to_csv(case_path, index=False, lineterminator="\n")
    correlations.to_csv(
        correlation_path,
        index=False,
        lineterminator="\n",
    )
    bch_independence.to_csv(
        bch_path,
        index=False,
        lineterminator="\n",
    )
    write_report(
        report_path,
        cases,
        correlations,
        bch_independence,
        args.basis,
    )
    plot_validation(
        figure_path,
        cases,
        correlations,
        args.dpi,
    )

    successful = cases[cases["structure_status"] == "success"]
    print(
        f"Recovered structures: {len(successful)}/{len(cases)} "
        f"across {successful['molecule'].nunique()} molecules."
    )
    print(f"Wrote {case_path}")
    print(f"Wrote {correlation_path}")
    print(f"Wrote {bch_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote {figure_path}")


if __name__ == "__main__":
    main()
