#!/usr/bin/env python3
"""Plot the corrected 20-order random-group benchmark summary.

Input
-----
A summary CSV produced by:

    analysis/benchmark_comparable_diatomics.py

Outputs
-------
1. Mean +/- one standard deviation versus active-space qubits.
2. Fraction of random schedules beating the signed-coefficient baseline.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "molecule",
    "n_qubits",
    "raw_jw_one_minus_overlap",
    "signed_baseline_one_minus_overlap",
    "jw_mean_one_minus_overlap",
    "jw_median_one_minus_overlap",
    "jw_std_one_minus_overlap",
    "jw_schedules_beating_baseline",
    "jw_fraction_beating_baseline",
    "fermionic_mean_one_minus_overlap",
    "fermionic_median_one_minus_overlap",
    "fermionic_std_one_minus_overlap",
    "fermionic_schedules_beating_baseline",
    "fermionic_fraction_beating_baseline",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("analysis/b2_20orders_summary.csv"),
        help=(
            "Corrected summary CSV produced by "
            "benchmark_comparable_diatomics.py."
        ),
    )

    parser.add_argument(
        "--molecule",
        default="B-B",
        help="Molecule to plot: B-B, C-C, or N-N.",
    )

    parser.add_argument(
        "--samples",
        type=int,
        default=20,
        help="Number of random-group schedules per coloring method.",
    )

    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("analysis/b2_20orders"),
        help="Output prefix without a filename extension.",
    )

    parser.add_argument(
        "--show-medians",
        action="store_true",
        help="Also display median curves on the error plot.",
    )

    return parser.parse_args()


def load_summary(path: Path, molecule: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Summary CSV not found: {path}")

    all_rows = pd.read_csv(path)

    missing = REQUIRED_COLUMNS.difference(all_rows.columns)
    if missing:
        raise ValueError(
            "Summary CSV is missing required columns: "
            f"{sorted(missing)}"
        )

    df = all_rows.loc[all_rows["molecule"].eq(molecule)].copy()

    if df.empty:
        available = sorted(
            all_rows["molecule"].dropna().astype(str).unique()
        )
        raise ValueError(
            f"No rows found for molecule={molecule!r}. "
            f"Available molecules: {available}"
        )

    numeric_columns = REQUIRED_COLUMNS - {"molecule"}

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="raise")

    if df["n_qubits"].duplicated().any():
        duplicate_qubits = sorted(
            df.loc[
                df["n_qubits"].duplicated(keep=False),
                "n_qubits",
            ].unique()
        )
        raise ValueError(
            "Expected one row per active-space size. "
            f"Duplicate qubit counts found: {duplicate_qubits}. "
            "Check whether multiple bond lengths or basis sets are present."
        )

    positive_columns = [
        "raw_jw_one_minus_overlap",
        "signed_baseline_one_minus_overlap",
        "jw_mean_one_minus_overlap",
        "jw_median_one_minus_overlap",
        "fermionic_mean_one_minus_overlap",
        "fermionic_median_one_minus_overlap",
    ]

    if (df[positive_columns] <= 0.0).any().any():
        raise ValueError(
            "All plotted errors must be positive because the y-axis "
            "uses logarithmic scaling."
        )

    return df.sort_values("n_qubits").reset_index(drop=True)


def molecule_title(molecule: str) -> str:
    labels = {
        "B-B": r"$B_2$",
        "C-C": r"$C_2$",
        "N-N": r"$N_2$",
    }
    return labels.get(molecule, molecule)


def plot_error_summary(
    df: pd.DataFrame,
    *,
    molecule: str,
    samples: int,
    output_path: Path,
    show_medians: bool,
) -> None:
    x = df["n_qubits"].to_numpy(dtype=float)
    tiny = np.finfo(float).tiny

    fig, ax = plt.subplots(figsize=(9.5, 6.2))

    # Fixed references
    ax.plot(
        x,
        df["signed_baseline_one_minus_overlap"],
        marker="D",
        linewidth=2.2,
        label="Signed coefficient + lexicographic baseline",
    )

    ax.plot(
        x,
        df["raw_jw_one_minus_overlap"],
        marker="x",
        linestyle="--",
        linewidth=1.7,
        label="Raw JW insertion order",
    )

    # JW random-group schedules
    jw_mean = df["jw_mean_one_minus_overlap"].to_numpy(dtype=float)
    jw_std = df["jw_std_one_minus_overlap"].to_numpy(dtype=float)

    jw_line = ax.plot(
        x,
        jw_mean,
        marker="o",
        linewidth=2.2,
        label="JW coloring: random-group mean",
    )[0]

    ax.fill_between(
        x,
        np.maximum(jw_mean - jw_std, tiny),
        jw_mean + jw_std,
        alpha=0.18,
        color=jw_line.get_color(),
        label="JW coloring: mean +/- 1 SD",
    )

    # Fermionic random-group schedules
    fermionic_mean = df[
        "fermionic_mean_one_minus_overlap"
    ].to_numpy(dtype=float)

    fermionic_std = df[
        "fermionic_std_one_minus_overlap"
    ].to_numpy(dtype=float)

    fermionic_line = ax.plot(
        x,
        fermionic_mean,
        marker="s",
        linewidth=2.2,
        label="Fermionic coloring: random-group mean",
    )[0]

    ax.fill_between(
        x,
        np.maximum(fermionic_mean - fermionic_std, tiny),
        fermionic_mean + fermionic_std,
        alpha=0.18,
        color=fermionic_line.get_color(),
        label="Fermionic coloring: mean +/- 1 SD",
    )

    # Optional medians
    if show_medians:
        ax.plot(
            x,
            df["jw_median_one_minus_overlap"],
            linestyle=":",
            linewidth=1.7,
            color=jw_line.get_color(),
            label="JW coloring: median",
        )

        ax.plot(
            x,
            df["fermionic_median_one_minus_overlap"],
            linestyle=":",
            linewidth=1.7,
            color=fermionic_line.get_color(),
            label="Fermionic coloring: median",
        )

    ax.set_yscale("log")
    ax.set_xticks(df["n_qubits"].astype(int))

    ax.set_xlabel("Active-space qubits")
    ax.set_ylabel(
        r"$1-|\langle\psi_{\mathrm{exact}}|"
        r"\psi_{\mathrm{Trotter}}\rangle|$"
    )

    ax.set_title(
        f"{molecule_title(molecule)}/STO-6G "
        "Random Color-Group Ordering\n"
        f"First-order Trotter, t = 1, 100 steps, "
        f"{samples} schedules per method"
    )

    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8.5)

    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_win_rates(
    df: pd.DataFrame,
    *,
    molecule: str,
    samples: int,
    output_path: Path,
) -> None:
    x = df["n_qubits"].to_numpy(dtype=float)

    jw_fraction = df[
        "jw_fraction_beating_baseline"
    ].to_numpy(dtype=float)

    fermionic_fraction = df[
        "fermionic_fraction_beating_baseline"
    ].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(9.5, 5.2))

    jw_line = ax.plot(
        x,
        100.0 * jw_fraction,
        marker="o",
        linewidth=2.2,
        label="JW coloring",
    )[0]

    fermionic_line = ax.plot(
        x,
        100.0 * fermionic_fraction,
        marker="s",
        linewidth=2.2,
        label="Fermionic coloring",
    )[0]

    # Display the number of successful schedules, such as 3/20.
    for row in df.itertuples(index=False):
        ax.annotate(
            f"{int(row.jw_schedules_beating_baseline)}/{samples}",
            (
                row.n_qubits,
                100.0 * row.jw_fraction_beating_baseline,
            ),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            color=jw_line.get_color(),
        )

        ax.annotate(
            (
                f"{int(row.fermionic_schedules_beating_baseline)}"
                f"/{samples}"
            ),
            (
                row.n_qubits,
                100.0 * row.fermionic_fraction_beating_baseline,
            ),
            xytext=(0, -14),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            color=fermionic_line.get_color(),
        )

    ax.set_xticks(df["n_qubits"].astype(int))
    ax.set_ylim(-3.0, 103.0)

    ax.set_xlabel("Active-space qubits")
    ax.set_ylabel("Schedules beating signed baseline (%)")

    ax.set_title(
        f"{molecule_title(molecule)} Random-Group Success Rate "
        f"({samples} schedules per method)"
    )

    ax.grid(True, alpha=0.25)
    ax.legend()

    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    if args.samples <= 0:
        raise ValueError("--samples must be positive.")

    df = load_summary(
        path=args.summary_csv,
        molecule=args.molecule,
    )

    error_output = args.output_prefix.with_name(
        args.output_prefix.name + "_error_summary.png"
    )

    win_rate_output = args.output_prefix.with_name(
        args.output_prefix.name + "_win_rates.png"
    )

    plot_error_summary(
        df,
        molecule=args.molecule,
        samples=args.samples,
        output_path=error_output,
        show_medians=args.show_medians,
    )

    plot_win_rates(
        df,
        molecule=args.molecule,
        samples=args.samples,
        output_path=win_rate_output,
    )

    print(df.to_string(index=False))
    print(f"\nSaved error summary: {error_output}")
    print(f"Saved win-rate plot: {win_rate_output}")


if __name__ == "__main__":
    main()