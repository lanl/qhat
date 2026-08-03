#!/usr/bin/env python3
"""Create the corrected B2 random-group ordering summary figure.

The fixed reference is the signed-coefficient/lexicographic ordering. Only
successful rows with schedule == "random_groups" are used for the JW- and
fermionic-coloring distributions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

QUBITS = (12, 16, 18)
GRAPH_LEVELS = ("jw", "fermionic")
BASELINE_ORDERING = "signed_coefficient_lexicographic"
SCHEDULE = "random_groups"
METRIC = "state_infidelity"


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Plot B2 random-group mean +/- standard deviation against the signed-coefficient baseline."
    )
    parser.add_argument(
        "--robustness-csv",
        type=Path,
        default=script_dir / "b2_coloring_robustness_results.csv",
    )
    parser.add_argument(
        "--baseline-csv",
        type=Path,
        default=script_dir / "b2_signed_coefficient_baseline_results.csv",
    )
    parser.add_argument(
        "--output-png",
        type=Path,
        default=script_dir / "b2_coloring_robustness_figures" / "b2_random_groups_summary.png",
    )
    parser.add_argument(
        "--output-summary-csv",
        type=Path,
        default=script_dir / "b2_random_groups_summary.csv",
    )
    return parser.parse_args()


def load_inputs(robustness_csv: Path, baseline_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not robustness_csv.is_file():
        raise FileNotFoundError(f"Robustness CSV not found: {robustness_csv}")
    if not baseline_csv.is_file():
        raise FileNotFoundError(f"Baseline CSV not found: {baseline_csv}")

    robustness = pd.read_csv(robustness_csv)
    baseline = pd.read_csv(baseline_csv)

    robustness_required = {"status", "n_qubits", "graph_level", "schedule", METRIC}
    baseline_required = {"status", "n_qubits", "ordering", METRIC}
    missing_robustness = robustness_required.difference(robustness.columns)
    missing_baseline = baseline_required.difference(baseline.columns)
    if missing_robustness:
        raise ValueError(f"Robustness CSV is missing columns: {sorted(missing_robustness)}")
    if missing_baseline:
        raise ValueError(f"Baseline CSV is missing columns: {sorted(missing_baseline)}")

    return robustness, baseline


def build_summary(robustness: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    baseline_rows = baseline.loc[
        baseline["status"].eq("success")
        & baseline["ordering"].eq(BASELINE_ORDERING)
        & baseline["n_qubits"].isin(QUBITS),
        ["n_qubits", METRIC],
    ].copy()

    duplicated = baseline_rows["n_qubits"].duplicated(keep=False)
    if duplicated.any():
        duplicate_qubits = sorted(baseline_rows.loc[duplicated, "n_qubits"].unique())
        raise ValueError(f"Expected one signed baseline per qubit count; duplicates found for {duplicate_qubits}")

    baseline_by_qubits = baseline_rows.set_index("n_qubits")[METRIC]
    missing_qubits = sorted(set(QUBITS).difference(baseline_by_qubits.index))
    if missing_qubits:
        raise ValueError(f"Missing signed baseline for qubit counts: {missing_qubits}")

    selected = robustness.loc[
        robustness["status"].eq("success")
        & robustness["schedule"].eq(SCHEDULE)
        & robustness["graph_level"].isin(GRAPH_LEVELS)
        & robustness["n_qubits"].isin(QUBITS),
        ["n_qubits", "graph_level", METRIC],
    ].copy()

    selected[METRIC] = pd.to_numeric(selected[METRIC], errors="coerce")
    if selected[METRIC].isna().any():
        raise ValueError("Selected random_groups rows contain missing or nonnumeric state_infidelity values.")
    if (selected[METRIC] <= 0).any():
        raise ValueError("state_infidelity must be positive for the logarithmic plot.")

    rows: list[dict[str, float | int | str]] = []
    for n_qubits in QUBITS:
        reference = float(baseline_by_qubits.loc[n_qubits])
        for graph_level in GRAPH_LEVELS:
            values = selected.loc[
                selected["n_qubits"].eq(n_qubits)
                & selected["graph_level"].eq(graph_level),
                METRIC,
            ]
            if values.empty:
                raise ValueError(f"No {SCHEDULE} rows for {n_qubits} qubits and graph_level={graph_level!r}.")

            rows.append(
                {
                    "n_qubits": n_qubits,
                    "graph_level": graph_level,
                    "schedule": SCHEDULE,
                    "number_of_samples": int(values.size),
                    "signed_coefficient_lexicographic_baseline": reference,
                    "mean_state_infidelity": float(values.mean()),
                    "median_state_infidelity": float(values.median()),
                    "std_state_infidelity": float(values.std(ddof=1)),
                    "minimum_state_infidelity": float(values.min()),
                    "maximum_state_infidelity": float(values.max()),
                    "number_beating_baseline": int((values < reference).sum()),
                    "fraction_beating_baseline": float((values < reference).mean()),
                }
            )

    return pd.DataFrame(rows)


def plot_summary(summary: pd.DataFrame, output_png: Path) -> None:
    x = np.asarray(QUBITS, dtype=float)
    fig, ax = plt.subplots(figsize=(8.6, 5.4), constrained_layout=True)

    baseline = (
        summary[["n_qubits", "signed_coefficient_lexicographic_baseline"]]
        .drop_duplicates()
        .set_index("n_qubits")
        .loc[list(QUBITS), "signed_coefficient_lexicographic_baseline"]
        .to_numpy(dtype=float)
    )
    ax.plot(
        x,
        baseline,
        marker="D",
        linewidth=2.0,
        linestyle="--",
        label="Signed coefficient + lexicographic baseline",
    )

    display_names = {
        "jw": "JW-coloring random groups",
        "fermionic": "Fermionic-coloring random groups",
    }

    for graph_level in GRAPH_LEVELS:
        group = summary.loc[summary["graph_level"].eq(graph_level)].set_index("n_qubits").loc[list(QUBITS)]
        mean = group["mean_state_infidelity"].to_numpy(dtype=float)
        std = group["std_state_infidelity"].to_numpy(dtype=float)
        lower = np.maximum(mean - std, np.finfo(float).tiny)
        upper = mean + std

        (line,) = ax.plot(
            x,
            mean,
            marker="o",
            linewidth=2.0,
            label=f"{display_names[graph_level]}: mean",
        )
        ax.fill_between(
            x,
            lower,
            upper,
            alpha=0.18,
            color=line.get_color(),
            label=f"{display_names[graph_level]}: mean +/- 1 std",
        )


    ax.set_yscale("log")
    ax.set_xticks(QUBITS)
    ax.set_xlabel("Number of qubits")
    ax.set_ylabel("State infidelity")
    ax.set_title(r"$B_2$ random-group Trotter ordering robustness")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncol=1)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    robustness, baseline = load_inputs(args.robustness_csv, args.baseline_csv)
    summary = build_summary(robustness, baseline)

    args.output_summary_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_summary_csv, index=False)
    plot_summary(summary, args.output_png)

    print(summary.to_string(index=False))
    print(f"\nSaved figure: {args.output_png}")
    print(f"Saved summary: {args.output_summary_csv}")


if __name__ == "__main__":
    main()

