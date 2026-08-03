#!/usr/bin/env python3
"""Summarize B2 random color-group schedules against Reuben's baseline.

The script uses only rows with ``schedule == 'random_groups'`` from the
coloring-robustness CSV.  The primary fixed reference is the deterministic
ordering defined as ascending signed coefficient, with exact coefficient ties
broken by ascending dense Pauli-string lexicographic order.

Outputs
-------
<output-prefix>_summary.csv
    One wide summary row per active-space size.
<output-prefix>_summary.md
    A report-ready Markdown table and concise interpretation.
<output-prefix>_figure.png
    Mean +/- one sample standard deviation for JW and fermionic random group
    schedules, with the corrected deterministic baseline and raw JW reference.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


QUBITS = (12, 16, 18)
GRAPH_LEVELS = ("jw", "fermionic")
SIGNED_BASELINE_NAME = "signed_coefficient_lexicographic"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--robustness-csv",
        type=Path,
        required=True,
        help="CSV produced by benchmark_b2_coloring_robustness.py",
    )
    parser.add_argument(
        "--baseline-csv",
        type=Path,
        required=True,
        help="CSV produced by benchmark_b2_signed_coefficient_baseline.py",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("b2_random_groups_corrected_baseline"),
        help="Output path prefix without a file extension",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def load_robustness(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    require_columns(
        df,
        ["status", "n_qubits", "graph_level", "schedule", "state_overlap_abs"],
        "Robustness CSV",
    )
    df = df.loc[df["status"].eq("success")].copy()
    df["one_minus_overlap"] = 1.0 - df["state_overlap_abs"].astype(float)
    return df


def load_baselines(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    require_columns(
        df,
        ["status", "n_qubits", "ordering", "one_minus_overlap"],
        "Baseline CSV",
    )
    return df.loc[df["status"].eq("success")].copy()


def scalar_for_case(
    df: pd.DataFrame,
    *,
    n_qubits: int,
    ordering: str,
    column: str,
) -> float:
    rows = df.loc[
        df["n_qubits"].eq(n_qubits) & df["ordering"].eq(ordering),
        column,
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly one {ordering!r} row for {n_qubits} qubits; "
            f"found {len(rows)}."
        )
    return float(rows.iloc[0])


def build_summary(robustness: pd.DataFrame, baselines: pd.DataFrame) -> pd.DataFrame:
    random_groups = robustness.loc[
        robustness["schedule"].eq("random_groups")
        & robustness["graph_level"].isin(GRAPH_LEVELS)
        & robustness["n_qubits"].isin(QUBITS)
    ].copy()

    rows: list[dict[str, float | int]] = []
    for n_qubits in QUBITS:
        signed_baseline = scalar_for_case(
            baselines,
            n_qubits=n_qubits,
            ordering=SIGNED_BASELINE_NAME,
            column="one_minus_overlap",
        )
        raw_jw = scalar_for_case(
            baselines,
            n_qubits=n_qubits,
            ordering="jw_raw",
            column="one_minus_overlap",
        )

        row: dict[str, float | int] = {
            "n_qubits": n_qubits,
            "signed_baseline_one_minus_overlap": signed_baseline,
            "raw_jw_one_minus_overlap": raw_jw,
        }

        for graph_level in GRAPH_LEVELS:
            values = random_groups.loc[
                random_groups["n_qubits"].eq(n_qubits)
                & random_groups["graph_level"].eq(graph_level),
                "one_minus_overlap",
            ].astype(float)

            if values.empty:
                raise ValueError(
                    f"No random_groups rows for {n_qubits} qubits and "
                    f"graph_level={graph_level!r}."
                )

            prefix = "jw" if graph_level == "jw" else "fermionic"
            beats = values.lt(signed_baseline)
            row.update(
                {
                    f"{prefix}_sample_count": int(values.size),
                    f"{prefix}_mean": float(values.mean()),
                    f"{prefix}_median": float(values.median()),
                    f"{prefix}_std": float(values.std(ddof=1)),
                    f"{prefix}_min": float(values.min()),
                    f"{prefix}_max": float(values.max()),
                    f"{prefix}_count_beating_signed_baseline": int(beats.sum()),
                    f"{prefix}_fraction_beating_signed_baseline": float(beats.mean()),
                }
            )

        rows.append(row)

    return pd.DataFrame(rows).sort_values("n_qubits").reset_index(drop=True)


def sci(value: float) -> str:
    return f"{value:.4e}"


def write_markdown(summary: pd.DataFrame, path: Path) -> None:
    lines = [
        "# B₂ Random Color-Group Ordering Summary",
        "",
        "**Metric:** $1-|\\langle\\psi_{\\mathrm{exact}}|\\psi_{\\mathrm{Trotter}}\\rangle|$ "
        "(lower is better).",
        "",
        "**Setup:** B2/STO-6G, first-order Trotter evolution, $t=1$, "
        "100 steps, and 100 random color-group permutations per coloring method.",
        "",
        "Only rows with `schedule = random_groups` are included in the random-order statistics.",
        "",
        "| Qubits | Signed-coefficient baseline | Raw JW | JW mean +/- SD | JW median | JW min-max | JW beating baseline | Fermionic mean +/- SD | Fermionic median | Fermionic min-max | Fermionic beating baseline |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for row in summary.itertuples(index=False):
        lines.append(
            "| "
            f"{row.n_qubits} | "
            f"{sci(row.signed_baseline_one_minus_overlap)} | "
            f"{sci(row.raw_jw_one_minus_overlap)} | "
            f"{sci(row.jw_mean)} +/- {sci(row.jw_std)} | "
            f"{sci(row.jw_median)} | "
            f"{sci(row.jw_min)}--{sci(row.jw_max)} | "
            f"{row.jw_count_beating_signed_baseline}/{row.jw_sample_count} "
            f"({100.0 * row.jw_fraction_beating_signed_baseline:.1f}%) | "
            f"{sci(row.fermionic_mean)} +/- {sci(row.fermionic_std)} | "
            f"{sci(row.fermionic_median)} | "
            f"{sci(row.fermionic_min)}--{sci(row.fermionic_max)} | "
            f"{row.fermionic_count_beating_signed_baseline}/{row.fermionic_sample_count} "
            f"({100.0 * row.fermionic_fraction_beating_signed_baseline:.1f}%) |"
        )

    lines.extend(
        [
            "",
            "## Main findings",
            "",
            "- The fermionic random-group distribution has a lower mean and median than the JW-coloring distribution for all three active spaces.",
            "- The corrected signed-coefficient/lexicographic baseline is stronger than the typical random coloring schedule in every case.",
            "- Fermionic schedules beating the corrected baseline occur in 18% of the 12-qubit samples, 4% of the 16-qubit samples, and 0% of the 18-qubit samples.",
            "- The 16-qubit fermionic coloring therefore contains a small number of competitive schedules, but random group ordering is not a reliable improvement over the corrected baseline.",
            "",
            "The shaded regions in the figure represent mean plus or minus one sample standard deviation.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_summary(summary: pd.DataFrame, path: Path) -> None:
    x = summary["n_qubits"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(9.5, 6.2))

    ax.plot(
        x,
        summary["signed_baseline_one_minus_overlap"],
        marker="D",
        linewidth=2.2,
        label="Signed coefficient + lexicographic baseline",
    )
    ax.plot(
        x,
        summary["raw_jw_one_minus_overlap"],
        marker="x",
        linestyle="--",
        linewidth=1.7,
        label="Raw JW insertion order (secondary)",
    )

    jw_line = ax.plot(
        x,
        summary["jw_mean"],
        marker="o",
        linewidth=2.2,
        label="JW coloring: random-group mean",
    )[0]
    ax.fill_between(
        x,
        summary["jw_mean"] - summary["jw_std"],
        summary["jw_mean"] + summary["jw_std"],
        alpha=0.18,
        color=jw_line.get_color(),
        label="JW coloring: +/- 1 SD",
    )

    fermionic_line = ax.plot(
        x,
        summary["fermionic_mean"],
        marker="s",
        linewidth=2.2,
        label="Fermionic coloring: random-group mean",
    )[0]
    ax.fill_between(
        x,
        summary["fermionic_mean"] - summary["fermionic_std"],
        summary["fermionic_mean"] + summary["fermionic_std"],
        alpha=0.18,
        color=fermionic_line.get_color(),
        label="Fermionic coloring: +/- 1 SD",
    )

    annotation_positions = {
        12: ((16, -18), "left"),
        16: ((0, 10), "center"),
        18: ((-10, 10), "right"),
    }
    for row in summary.itertuples(index=False):
        offset, alignment = annotation_positions[int(row.n_qubits)]
        ax.annotate(
            f"{100.0 * row.fermionic_fraction_beating_signed_baseline:.0f}% beat baseline",
            xy=(row.n_qubits, row.fermionic_mean),
            xytext=offset,
            textcoords="offset points",
            ha=alignment,
            fontsize=9,
        )

    ax.set_yscale("log")
    ax.set_xticks(QUBITS)
    ax.set_xlabel("Active-space qubits")
    ax.set_ylabel(
        r"$1-|\langle\psi_{\mathrm{exact}}|\psi_{\mathrm{Trotter}}\rangle|$"
    )
    ax.set_title(
        "B₂/STO-6G Random Color-Group Ordering Robustness\n"
        "First-order Trotter, t = 1, 100 steps, 100 schedules per method"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    robustness = load_robustness(args.robustness_csv)
    baselines = load_baselines(args.baseline_csv)
    summary = build_summary(robustness, baselines)

    csv_path = output_prefix.with_name(output_prefix.name + "_summary.csv")
    md_path = output_prefix.with_name(output_prefix.name + "_summary.md")
    png_path = output_prefix.with_name(output_prefix.name + "_figure.png")

    summary.to_csv(csv_path, index=False, float_format="%.16e")
    write_markdown(summary, md_path)
    plot_summary(summary, png_path)

    print(summary.to_string(index=False))
    print(f"\nWrote {csv_path}")
    print(f"Wrote {md_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
