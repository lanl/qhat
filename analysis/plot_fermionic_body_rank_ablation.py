#!/usr/bin/env python3

"""
Plot the 20-case selective 1B-vs-2B fermionic parent-ordering ablation.

For each Hamiltonian:
  - 1B shuffle:
      randomly permute only 1B parent blocks among the original 1B slots
  - 2B shuffle:
      randomly permute only 2B parent blocks among the original 2B slots

The figure reports the median loss relative to the signed fermionic
reference:

    BCH loss =
        R_BCH(shuffle) / R_BCH(reference)

    Trotter-error loss =
        E(shuffle) / E(reference)

Values > 1 mean that the shuffle is worse than the signed reference.

Expected input:
    analysis/fermionic_body_rank_ablation_20case.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCHEDULE_1B = "signed_parent_1b_blocks_randomized"
SCHEDULE_2B = "signed_parent_2b_blocks_randomized"

BCH_COLUMN = "bch_cancellation_ratio_to_signed_reference"
ERROR_COLUMN = "one_minus_overlap_ratio_to_signed_reference"


MOLECULE_LABELS = {
    "B-B": r"B$_2$",
    "Be-Be": r"Be$_2$",
    "F-F": r"F$_2$",
    "Li-Li": r"Li$_2$",
    "B-H": "B-H",
    "Be-H": "Be-H",
    "Li-H": "Li-H",
    "H2O": r"H$_2$O",
    "NH3": r"NH$_3$",
    "CH4": r"CH$_4$",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot selective one-body versus two-body fermionic "
            "parent-ordering ablations."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=Path("analysis/fermionic_body_rank_ablation_20case.csv"),
        help="Full body-rank ablation CSV.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/fermionic_body_rank_ablation.pdf"),
        help="Output figure. PDF is recommended for the manuscript.",
    )

    parser.add_argument(
        "--expected-samples",
        type=int,
        default=20,
        help="Expected randomized samples for each case and body rank.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for raster output.",
    )

    return parser.parse_args()


def check_required_columns(df: pd.DataFrame) -> None:
    required = {
        "case_id",
        "schedule",
        "sample_index",
        "molecule",
        "ao",
        "av",
        "n_qubits",
        BCH_COLUMN,
        ERROR_COLUMN,
    }

    missing = sorted(required.difference(df.columns))

    if missing:
        raise ValueError(
            "Input CSV is missing required columns:\n  "
            + "\n  ".join(missing)
        )


def build_case_label(row: pd.Series) -> str:
    molecule = str(row["molecule"])
    molecule = MOLECULE_LABELS.get(molecule, molecule)

    ao = int(row["ao"])
    av = int(row["av"])
    n_qubits = int(row["n_qubits"])

    return rf"{molecule} {ao}/{av} ({n_qubits}q)"


def load_and_summarize(
    input_path: Path,
    expected_samples: int,
) -> pd.DataFrame:

    df = pd.read_csv(input_path)
    
    # Normalize active-space column names used by this plotting script.
    df = df.rename(
        columns={
            "active_occupied": "ao",
            "active_vacant": "av",
        }
    )

    check_required_columns(df)

    # Keep successful rows if the benchmark CSV has a status column.
    if "status" in df.columns:
        df = df[df["status"] == "success"].copy()

    # Keep only the two selective body-rank interventions.
    df = df[
        df["schedule"].isin(
            [
                SCHEDULE_1B,
                SCHEDULE_2B,
            ]
        )
    ].copy()

    if df.empty:
        raise ValueError(
            "No 1B/2B selective-shuffle rows were found in the input CSV."
        )

    # Make sure the plotted quantities are numeric.
    df[BCH_COLUMN] = pd.to_numeric(df[BCH_COLUMN], errors="coerce")
    df[ERROR_COLUMN] = pd.to_numeric(df[ERROR_COLUMN], errors="coerce")

    df = df[
        np.isfinite(df[BCH_COLUMN])
        & np.isfinite(df[ERROR_COLUMN])
        & (df[BCH_COLUMN] > 0.0)
        & (df[ERROR_COLUMN] > 0.0)
    ].copy()

    # Check that every case has the complete set of random samples.
    counts = (
        df.groupby(["case_id", "schedule"])
        .size()
        .rename("samples")
        .reset_index()
    )

    bad_counts = counts[counts["samples"] != expected_samples]

    if not bad_counts.empty:
        print(
            "\nWARNING: some case/schedule combinations do not have "
            f"{expected_samples} samples:"
        )
        print(bad_counts.to_string(index=False))
        print()

    # Keep only cases for which BOTH the 1B and 2B interventions exist.
    schedules_per_case = (
        df.groupby("case_id")["schedule"]
        .nunique()
    )

    complete_case_ids = schedules_per_case[
        schedules_per_case == 2
    ].index

    df = df[df["case_id"].isin(complete_case_ids)].copy()

    # Median over randomized schedules for each case and intervention.
    summary = (
        df.groupby(
            [
                "case_id",
                "schedule",
                "molecule",
                "ao",
                "av",
                "n_qubits",
            ],
            as_index=False,
        )
        .agg(
            bch_loss=(BCH_COLUMN, "median"),
            error_loss=(ERROR_COLUMN, "median"),
            samples=("sample_index", "count"),
        )
    )

    # Pivot so each case has one 1B value and one 2B value.
    metadata = (
        summary[
            [
                "case_id",
                "molecule",
                "ao",
                "av",
                "n_qubits",
            ]
        ]
        .drop_duplicates("case_id")
        .set_index("case_id")
    )

    bch = summary.pivot(
        index="case_id",
        columns="schedule",
        values="bch_loss",
    )

    error = summary.pivot(
        index="case_id",
        columns="schedule",
        values="error_loss",
    )

    result = metadata.copy()

    result["bch_1b"] = bch[SCHEDULE_1B]
    result["bch_2b"] = bch[SCHEDULE_2B]

    result["error_1b"] = error[SCHEDULE_1B]
    result["error_2b"] = error[SCHEDULE_2B]

    result = result.reset_index()

    # Sort by active-space size, then case ID.
    # This reproduces the ordering used in the presentation figure.
    result = result.sort_values(
        ["n_qubits", "case_id"],
        kind="stable",
    ).reset_index(drop=True)

    result["label"] = result.apply(build_case_label, axis=1)

    return result


def print_summary(df: pd.DataFrame) -> None:
    n_cases = len(df)

    bch_2b_larger = int(
        (df["bch_2b"] > df["bch_1b"]).sum()
    )

    error_2b_larger = int(
        (df["error_2b"] > df["error_1b"]).sum()
    )

    median_bch_1b = float(df["bch_1b"].median())
    median_bch_2b = float(df["bch_2b"].median())

    median_error_1b = float(df["error_1b"].median())
    median_error_2b = float(df["error_2b"].median())

    print()
    print("Selective body-rank ablation summary")
    print("------------------------------------")
    print(f"Complete cases: {n_cases}")

    print(
        "2B BCH loss > 1B BCH loss: "
        f"{bch_2b_larger}/{n_cases}"
    )

    print(
        "2B error loss > 1B error loss: "
        f"{error_2b_larger}/{n_cases}"
    )

    print(
        "Median BCH loss:"
        f"  1B = {median_bch_1b:.3g}x,"
        f"  2B = {median_bch_2b:.3g}x"
    )

    print(
        "Median Trotter-error loss:"
        f"  1B = {median_error_1b:.3g}x,"
        f"  2B = {median_error_2b:.3g}x"
    )

    print()


def nice_log_limits(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]

    if len(values) == 0:
        return 0.5, 2.0

    # Always include the signed-reference value 1.
    minimum = min(float(values.min()), 1.0)
    maximum = max(float(values.max()), 1.0)

    low_power = np.floor(np.log10(minimum))
    high_power = np.ceil(np.log10(maximum))

    lower = 10.0 ** low_power
    upper = 10.0 ** high_power

    # Give a little room on the left of 1x.
    if lower == 1.0:
        lower = 0.6

    return lower, upper


def draw_panel(
    ax: plt.Axes,
    y: np.ndarray,
    one_body: np.ndarray,
    two_body: np.ndarray,
    title: str,
    xlabel: str,
) -> None:

    # Connecting line between the 1B and 2B median for each case.
    for yi, x1, x2 in zip(y, one_body, two_body):
        ax.plot(
            [x1, x2],
            [yi, yi],
            color="0.72",
            linewidth=1.0,
            zorder=1,
        )

    # 1B intervention.
    ax.scatter(
        one_body,
        y,
        marker="o",
        s=28,
        color="#1f77b4",
        edgecolors="none",
        label="1B shuffle",
        zorder=3,
    )

    # 2B intervention.
    ax.scatter(
        two_body,
        y,
        marker="^",
        s=36,
        color="#ff7f0e",
        edgecolors="none",
        label="2B shuffle",
        zorder=4,
    )

    # Signed fermionic reference.
    ax.axvline(
        1.0,
        linestyle="--",
        linewidth=1.0,
        color="#1f4e79",
        zorder=0,
    )

    ax.set_xscale("log")

    limits = nice_log_limits(
        np.concatenate([one_body, two_body])
    )
    ax.set_xlim(*limits)

    ax.set_title(
        title,
        fontsize=12,
        fontweight="bold",
        color="#002b66",
    )

    ax.set_xlabel(xlabel)

    ax.grid(
        axis="x",
        which="both",
        alpha=0.15,
        linewidth=0.7,
    )

    ax.grid(
        axis="y",
        which="major",
        visible=False,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        frameon=False,
        loc="lower right",
        fontsize=8,
    )


def make_plot(
    df: pd.DataFrame,
    output_path: Path,
    dpi: int,
) -> None:

    labels = df["label"].tolist()

    y = np.arange(len(df))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.5, 6.0),
        sharey=True,
    )

    draw_panel(
        axes[0],
        y,
        df["bch_1b"].to_numpy(float),
        df["bch_2b"].to_numpy(float),
        title="BCH cancellation loss",
        xlabel=r"median $C_{\mathrm{shuffle}}/C_{\mathrm{ref}}$",
    )

    draw_panel(
        axes[1],
        y,
        df["error_1b"].to_numpy(float),
        df["error_2b"].to_numpy(float),
        title="Measured Trotter-error loss",
        xlabel=r"median $E_{\mathrm{shuffle}}/E_{\mathrm{ref}}$",
    )

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(
        labels,
        fontsize=8,
    )

    # Put the first case at the top.
    axes[0].invert_yaxis()

    fig.subplots_adjust(
        left=0.18,
        right=0.985,
        bottom=0.12,
        top=0.92,
        wspace=0.07,
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    # Also create a PNG automatically when PDF is requested.
    if output_path.suffix.lower() == ".pdf":
        png_path = output_path.with_suffix(".png")

        fig.savefig(
            png_path,
            dpi=dpi,
            bbox_inches="tight",
        )

        print(f"Saved: {png_path}")

    print(f"Saved: {output_path}")

    plt.close(fig)


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(
            f"Input CSV not found: {args.input}"
        )

    summary = load_and_summarize(
        args.input,
        args.expected_samples,
    )

    print_summary(summary)

    make_plot(
        summary,
        args.output,
        args.dpi,
    )


if __name__ == "__main__":
    main()