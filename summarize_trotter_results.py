#!/usr/bin/env python3
"""Summarize focused QHAT Trotter benchmark outputs as CSV and Markdown."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


SCHEME_LABELS = {
    "fermionic_commutation_then_JW": "Fermionic coloring → JW",
    "JW_Pauli_string_commutation": "Direct JW Pauli coloring",
}


def fmt(value: object) -> str:
    if isinstance(value, (float, np.floating)):
        if np.isnan(value):
            return "—"
        return f"{value:.6g}"
    return str(value)


def markdown_table(frame: pd.DataFrame, columns: list[str] | None = None) -> str:
    selected = frame if columns is None else frame[columns]
    headers = [str(column) for column in selected.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in selected.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(fmt(value) for value in row) + " |")
    return "\n".join(lines)


def robustness_summary(trials: pd.DataFrame) -> pd.DataFrame:
    keys = [
        "case_id",
        "molecule",
        "bond_length",
        "basis",
        "active_occupied",
        "active_vacant",
        "qubits",
        "scheme",
        "formula_order",
        "steps",
        "total_time",
    ]
    rows: list[dict[str, object]] = []
    for key, group in trials.groupby(keys, sort=False):
        common = dict(zip(keys, key))
        colored = group[group.ordering_kind.eq("colored")]
        random = group[group.ordering_kind.eq("random")]
        if len(colored) != 1 or random.empty:
            raise ValueError(f"Expected one colored row and random rows for {key}")
        default = colored.iloc[0]
        random_errors = random.spectral_error.to_numpy()
        random_leakage = random.particle_number_leakage.to_numpy()
        best_observed = min(float(default.spectral_error), float(random_errors.min()))
        rows.append(
            {
                **common,
                "colored_error": float(default.spectral_error),
                "random_best_observed_error": float(random_errors.min()),
                "random_p10_error": float(np.quantile(random_errors, 0.1)),
                "random_median_error": float(np.median(random_errors)),
                "random_p90_error": float(np.quantile(random_errors, 0.9)),
                "random_worst_error": float(random_errors.max()),
                "best_observed_error": best_observed,
                "colored_over_best_observed": float(default.spectral_error)
                / max(best_observed, np.finfo(float).tiny),
                "colored_over_random_median": float(default.spectral_error)
                / max(float(np.median(random_errors)), np.finfo(float).tiny),
                "fraction_random_better_than_colored": float(
                    np.mean(random_errors < float(default.spectral_error))
                ),
                "colored_particle_number_leakage": float(
                    default.particle_number_leakage
                ),
                "random_best_particle_number_leakage": float(random_leakage.min()),
                "random_median_particle_number_leakage": float(
                    np.median(random_leakage)
                ),
                "random_worst_particle_number_leakage": float(random_leakage.max()),
                "colored_state_infidelity": float(default.state_infidelity),
                "random_median_state_infidelity": float(
                    np.median(random.state_infidelity)
                ),
            }
        )
    summary = pd.DataFrame(rows)
    comparison_keys = ["case_id", "formula_order", "steps", "total_time"]
    summary["global_best_observed_error"] = summary.groupby(comparison_keys)[
        "best_observed_error"
    ].transform("min")
    denominator = summary.global_best_observed_error.clip(lower=np.finfo(float).tiny)
    summary["scheme_best_over_global_best"] = summary.best_observed_error / denominator
    summary["colored_over_global_best"] = summary.colored_error / denominator
    return summary


def build_markdown(
    config: dict[str, object],
    graphs: pd.DataFrame,
    commutators: pd.DataFrame,
    summary: pd.DataFrame,
) -> str:
    lines = [
        "# Focused Fermionic-vs-JW Trotter Benchmark",
        "",
        "## Benchmark definition",
        "",
        (
            f"The benchmark uses {config['random_orderings']} seeded random permutations "
            "of each decomposition's nonidentity vertices. The colored baseline flattens "
            "deterministic largest-degree-first color groups. The same sampled permutations "
            "are reused across first/second order and all step counts."
        ),
        "",
        f"- Evolution time: `{config['total_time']}`",
        f"- Steps: `{config['steps']}`",
        "- Implemented orders: first and second only",
        f"- Random seed: `{config['seed']}`",
        "- ‘Best observed’ means best among the colored baseline and sampled orders; it is not a proven optimum.",
        "",
        "Fermionic tensor monomials are first put into the same descending-index normal-order convention used by `h2_fermionic.ipynb`, and algebraically equivalent permutations are combined before Hermitian pairing. This keeps the comparison independent of OpenFermion without changing its term convention.",
        "",
        "For color block `B_c = sum(H_i for i in color c)`, the report records pairwise `||[B_c,B_d]||₂`, the norm of the ordered commutator sum as a first-order diagnostic, and the sum of left/right nested-commutator norms as a second-order diagnostic. These are coefficient-weighted explanatory proxies, not rigorous rankings of the full finite-time error.",
        "",
        "## Selected cases",
        "",
    ]

    case_columns = [
        "case_id",
        "molecule",
        "bond_length",
        "basis",
        "active_occupied",
        "active_vacant",
        "qubits",
        "selection_reason",
    ]
    cases = graphs.drop_duplicates("case_id")[case_columns]
    lines.extend([markdown_table(cases), ""])

    for case_id in cases.case_id:
        case_graphs = graphs[graphs.case_id.eq(case_id)].copy()
        case_summary = summary[summary.case_id.eq(case_id)].copy()
        title_row = case_graphs.iloc[0]
        lines.extend(
            [
                f"## {case_id}: {title_row['molecule']} at L={title_row['bond_length']}",
                "",
                title_row["selection_reason"],
                "",
                "### Graph and coefficient-weighted diagnostics",
                "",
            ]
        )
        graph_table = case_graphs[
            [
                "scheme",
                "vertices",
                "noncommuting_edges",
                "colors",
                "block_commutator_norm_sum",
                "ordered_commutator_sum_norm",
                "nested_commutator_norm_sum",
            ]
        ].copy()
        graph_table["scheme"] = graph_table.scheme.map(SCHEME_LABELS)
        lines.extend([markdown_table(graph_table), ""])

        schedule_table = case_graphs[
            ["scheme", "color_order_json", "color_groups_json"]
        ].copy()
        schedule_table["scheme"] = schedule_table.scheme.map(SCHEME_LABELS)
        lines.extend(
            [
                "The color order and vertex membership used by the baseline are:",
                "",
                markdown_table(schedule_table),
                "",
            ]
        )

        pair_table = commutators[commutators.case_id.eq(case_id)][
            [
                "scheme",
                "left_color",
                "right_color",
                "commutator_spectral_norm",
                "normalized_commutator_norm",
                "left_nested_commutator_norm",
                "right_nested_commutator_norm",
            ]
        ].copy()
        pair_table["scheme"] = pair_table.scheme.map(SCHEME_LABELS)
        lines.extend(
            [
                "Pairwise color-block diagnostics:",
                "",
                markdown_table(pair_table),
                "",
                "### Colored baseline errors",
                "",
            ]
        )
        error_table = case_summary[
            ["scheme", "formula_order", "steps", "colored_error"]
        ].copy()
        error_table["scheme"] = error_table.scheme.map(SCHEME_LABELS)
        lines.extend([markdown_table(error_table), ""])

        lines.extend(["### Ordering robustness", ""])
        robust_table = case_summary[case_summary.steps.eq(max(config["steps"]))][
            [
                "scheme",
                "formula_order",
                "colored_error",
                "random_best_observed_error",
                "random_median_error",
                "random_p90_error",
                "colored_over_best_observed",
                "scheme_best_over_global_best",
                "colored_over_global_best",
                "fraction_random_better_than_colored",
            ]
        ].copy()
        robust_table["scheme"] = robust_table.scheme.map(SCHEME_LABELS)
        lines.extend([markdown_table(robust_table), ""])

        lines.extend(["### Particle-number leakage", ""])
        leakage_table = case_summary[case_summary.steps.eq(max(config["steps"]))][
            [
                "scheme",
                "formula_order",
                "colored_particle_number_leakage",
                "random_median_particle_number_leakage",
                "random_worst_particle_number_leakage",
            ]
        ].copy()
        leakage_table["scheme"] = leakage_table.scheme.map(SCHEME_LABELS)
        lines.extend([markdown_table(leakage_table), ""])

    lines.extend(
        [
            "## Interpretation",
            "",
            "The B–B STO-6G and HGBS-5 cases have the same molecule, geometry, active space, graph counts, and deterministic color memberships. Nevertheless, the preferred decomposition reverses when the basis changes. This is evidence that coefficient-weighted block commutators—and for second order, nested commutators—capture behavior that graph topology alone misses. The three-case study supports this explanation but does not establish a universal predictor.",
            "",
            "Fermionic vertices preserve particle number individually. Random JW Pauli-string orders can temporarily separate cancellation partners and leak out of the target particle-number sector. The colored JW baseline may suppress that leakage when the relevant strings remain grouped. In these cases the leakage is much smaller than the spectral error, so it is a physical robustness diagnostic rather than the primary explanation for the operator-norm error.",
            "",
            "The random-ordering ensemble also shows whether a colored result is robust or merely one favorable schedule. Compare the colored error with the random median and best observed error rather than treating one deterministic ordering as definitive.",
            "",
            "Canonical normal-order aggregation is itself important: compared with the earlier raw tensor-slot convention, it reduces the selected four-qubit fermionic graph from 23 vertices/32 edges to 13/12 and each selected six-qubit fermionic graph from 55/360 to 28/90. The colored errors remain unchanged in these three cases because the combined permutations occupied the same commuting blocks, while the random-order statistics change because the sampling unit is now a genuine canonical fermionic term. This is another concrete reason that raw graph size is representation-dependent.",
            "",
            "## Scope limits",
            "",
            "- No fourth-order product formula is included.",
            "- No dense all-case 8-qubit sweep is included.",
            "- The random ensemble samples vertex orders, not additional coloring algorithms.",
            "- Best observed errors are sampled results, not mathematical optima.",
            "- Nested-commutator sums are pairwise heuristics, not the complete schedule-weighted BCH error generator.",
            "- Normalized block-commutator values include full block norms and are sensitive to scalar energy shifts; the absolute norms are the primary diagnostics.",
            "- No gate-count or circuit-depth claim is made.",
        ]
    )
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=Path("trotter_benchmark_results"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    results_dir = args.results_dir.resolve()
    trials = pd.read_csv(results_dir / "ordering_trials.csv")
    graphs = pd.read_csv(results_dir / "graph_diagnostics.csv")
    commutators = pd.read_csv(results_dir / "block_commutators.csv")
    config = json.loads((results_dir / "benchmark_config.json").read_text())

    summary = robustness_summary(trials)
    summary.to_csv(results_dir / "robustness_summary.csv", index=False)
    report = build_markdown(config, graphs, commutators, summary)
    (results_dir / "benchmark_summary.md").write_text(report)
    print(f"Wrote {results_dir / 'robustness_summary.csv'}")
    print(f"Wrote {results_dir / 'benchmark_summary.md'}")


if __name__ == "__main__":
    main()
