#!/usr/bin/env python3
"""Compare fermionic-aware ordering with corrected deterministic JW baselines.

The input is the exact case-level join produced by
``analyze_jw_magnitude_baseline.py``.  The fair JW comparator is the smaller
error from signed-coefficient JW and magnitude-descending JW for the same
Hamiltonian and numerical settings.  This script writes auditable case,
molecule, molecule/basis, and condition summaries and two publication-style
figures.
"""

# ruff: noqa: E402  # Configure the headless Matplotlib backend before imports.

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path

_MPL_CACHE = Path(tempfile.gettempdir()) / "qhat-matplotlib-cache"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path("analysis/jw_magnitude_descending_case_comparison.csv")
DEFAULT_OUTDIR = Path("analysis/fermionic_aware_performance")
DEFAULT_ERROR_FLOOR = 1.0e-12

MOLECULE_LABELS = {
    "B-B": r"B$_2$",
    "Be-Be": r"Be$_2$",
    "F-F": r"F$_2$",
    "Li-Li": r"Li$_2$",
    "O-O": r"O$_2$",
    "C-C": r"C$_2$",
    "B-H": "BH",
    "Be-H": "BeH",
    "Li-H": "LiH",
    "BeH2": r"BeH$_2$",
    "H2O": r"H$_2$O",
    "NH3": r"NH$_3$",
    "CH4": r"CH$_4$",
    "N-N": r"N$_2$",
}

BASIS_STYLES = {
    "hgbs-5": {"marker": "o", "label": "HGBS-5"},
    "sto-6g": {"marker": "s", "label": "STO-6G"},
}

OUTCOME_COLORS = {
    "fermionic_win": "#2a9d8f",
    "tie": "#8d99ae",
    "strongest_jw_win": "#e76f51",
}

REQUIRED_COLUMNS = {
    "case_id",
    "tensor_path",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "trotter_steps",
    "evolution_time",
    "coefficient_tolerance",
    "jw_signed_one_minus_overlap",
    "jw_magnitude_one_minus_overlap",
    "fermionic_signed_one_minus_overlap",
    "jw_signed_source_csv",
    "jw_magnitude_source_csv",
    "fermionic_signed_source_csv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--error-floor", type=float, default=DEFAULT_ERROR_FLOOR)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "pdf", "svg"),
        default=["png", "pdf"],
    )
    return parser.parse_args()


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.bbox": "tight",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.2,
            "grid.linewidth": 0.7,
        }
    )


def close_enough(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-6, abs_tol=1.0e-14)


def classify_outcome(fermionic: float, strongest_jw: float) -> str:
    if close_enough(fermionic, strongest_jw):
        return "tie"
    return "fermionic_win" if fermionic < strongest_jw else "strongest_jw_win"


def read_input(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Input is missing columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError(f"Input contains no rows: {path}")
    return frame


def build_case_table(frame: pd.DataFrame, error_floor: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, source in frame.iterrows():
        signed = float(source["jw_signed_one_minus_overlap"])
        magnitude = float(source["jw_magnitude_one_minus_overlap"])
        strongest = min(signed, magnitude)
        if close_enough(signed, magnitude):
            strongest_name = "tie"
        elif signed < magnitude:
            strongest_name = "jw_signed"
        else:
            strongest_name = "jw_magnitude"

        fermionic_raw = source["fermionic_signed_one_minus_overlap"]
        available = pd.notna(fermionic_raw)
        valid = False
        exclusion_reason = ""
        outcome = "fermionic_not_available"
        fermionic = math.nan
        ratio_signed = math.nan
        ratio_magnitude = math.nan
        ratio_strongest = math.nan
        advantage = math.nan
        log_advantage = math.nan

        if not available:
            exclusion_reason = "fermionic_not_available"
        else:
            fermionic = float(fermionic_raw)
            finite_positive = (
                math.isfinite(fermionic)
                and math.isfinite(strongest)
                and fermionic >= 0.0
                and strongest >= 0.0
            )
            if not finite_positive:
                exclusion_reason = "invalid_error"
                outcome = "invalid_error"
            elif min(fermionic, strongest) <= error_floor:
                exclusion_reason = "below_numerical_floor"
                outcome = "below_numerical_floor"
            elif fermionic == 0.0 or signed == 0.0 or magnitude == 0.0:
                exclusion_reason = "zero_error"
                outcome = "below_numerical_floor"
            else:
                valid = True
                ratio_signed = fermionic / signed
                ratio_magnitude = fermionic / magnitude
                ratio_strongest = fermionic / strongest
                advantage = strongest / fermionic
                log_advantage = math.log10(advantage)
                outcome = classify_outcome(fermionic, strongest)

        occupied = int(source["active_occupied"])
        vacant = int(source["active_vacant"])
        n_qubits = int(source["n_qubits"])
        rows.append(
            {
                "case_id": source["case_id"],
                "tensor_path": source["tensor_path"],
                "molecule": source["molecule"],
                "bond_length": source["bond_length"],
                "basis": source["basis"],
                "active_occupied": occupied,
                "active_vacant": vacant,
                "occupied_fraction": occupied / n_qubits,
                "active_space_imbalance": abs(occupied - vacant) / n_qubits,
                "n_qubits": n_qubits,
                "trotter_steps": int(source["trotter_steps"]),
                "evolution_time": float(source["evolution_time"]),
                "coefficient_tolerance": float(source["coefficient_tolerance"]),
                "numerical_error_floor": error_floor,
                "jw_signed_one_minus_overlap": signed,
                "jw_magnitude_one_minus_overlap": magnitude,
                "strongest_jw_one_minus_overlap": strongest,
                "strongest_jw_ordering": strongest_name,
                "fermionic_aware_one_minus_overlap": fermionic,
                "fermionic_to_jw_signed_ratio": ratio_signed,
                "fermionic_to_jw_magnitude_ratio": ratio_magnitude,
                "fermionic_to_strongest_jw_ratio": ratio_strongest,
                "fermionic_advantage_factor": advantage,
                "log10_fermionic_advantage": log_advantage,
                "fermionic_available": available,
                "valid_comparison": valid,
                "outcome": outcome,
                "exclusion_reason": exclusion_reason,
                "jw_signed_source_csv": source["jw_signed_source_csv"],
                "jw_magnitude_source_csv": source["jw_magnitude_source_csv"],
                "fermionic_aware_source_csv": source[
                    "fermionic_signed_source_csv"
                ],
            }
        )

    cases = pd.DataFrame(rows)
    valid = cases[cases["valid_comparison"]]
    ranks = valid["fermionic_advantage_factor"].rank(
        method="min", ascending=False
    )
    cases["global_advantage_rank"] = pd.Series(pd.NA, index=cases.index, dtype="Int64")
    cases.loc[valid.index, "global_advantage_rank"] = ranks.astype("Int64")
    return cases.sort_values(
        [
            "molecule",
            "basis",
            "n_qubits",
            "active_occupied",
            "active_vacant",
            "case_id",
        ]
    ).reset_index(drop=True)


def geometric_mean(values: pd.Series) -> float:
    array = values.to_numpy(dtype=float)
    return float(np.exp(np.mean(np.log(array))))


def summarize_subset(label: str, subset: pd.DataFrame) -> dict[str, object]:
    valid = subset[subset["valid_comparison"]]
    advantage = valid["fermionic_advantage_factor"]
    wins = int((valid["outcome"] == "fermionic_win").sum())
    ties = int((valid["outcome"] == "tie").sum())
    losses = int((valid["outcome"] == "strongest_jw_win").sum())
    return {
        "group": label,
        "cases_total": len(subset),
        "fermionic_available": int(subset["fermionic_available"].sum()),
        "valid_comparisons": len(valid),
        "excluded_below_floor": int(
            (subset["exclusion_reason"] == "below_numerical_floor").sum()
        ),
        "fermionic_not_available": int(
            (subset["exclusion_reason"] == "fermionic_not_available").sum()
        ),
        "fermionic_wins": wins,
        "ties": ties,
        "strongest_jw_wins": losses,
        "fermionic_win_fraction": wins / len(valid) if len(valid) else math.nan,
        "fermionic_nonloss_fraction": (
            (wins + ties) / len(valid) if len(valid) else math.nan
        ),
        "median_fermionic_advantage": (
            float(advantage.median()) if len(valid) else math.nan
        ),
        "geometric_mean_fermionic_advantage": (
            geometric_mean(advantage) if len(valid) else math.nan
        ),
        "q25_fermionic_advantage": (
            float(advantage.quantile(0.25)) if len(valid) else math.nan
        ),
        "q75_fermionic_advantage": (
            float(advantage.quantile(0.75)) if len(valid) else math.nan
        ),
        "maximum_fermionic_advantage": (
            float(advantage.max()) if len(valid) else math.nan
        ),
        "minimum_fermionic_advantage": (
            float(advantage.min()) if len(valid) else math.nan
        ),
        "median_fermionic_to_jw_signed_ratio": (
            float(valid["fermionic_to_jw_signed_ratio"].median())
            if len(valid)
            else math.nan
        ),
        "median_fermionic_to_strongest_jw_ratio": (
            float(valid["fermionic_to_strongest_jw_ratio"].median())
            if len(valid)
            else math.nan
        ),
    }


def build_summaries(cases: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    molecule_rows = [summarize_subset("ALL", cases)]
    for molecule, group in cases.groupby("molecule", sort=True):
        row = summarize_subset(str(molecule), group)
        row["molecule"] = molecule
        molecule_rows.append(row)
    molecule_summary = pd.DataFrame(molecule_rows)
    molecule_summary["molecule"] = molecule_summary["molecule"].fillna("ALL")
    molecule_summary = molecule_summary.drop(columns="group")

    basis_rows: list[dict[str, object]] = []
    for (molecule, basis), group in cases.groupby(
        ["molecule", "basis"], sort=True
    ):
        row = summarize_subset(f"{molecule}/{basis}", group)
        row["molecule"] = molecule
        row["basis"] = basis
        basis_rows.append(row)
    basis_summary = pd.DataFrame(basis_rows).drop(columns="group")
    return molecule_summary, basis_summary


def qubit_bin(n_qubits: int) -> str:
    if n_qubits <= 8:
        return "04-08"
    if n_qubits <= 12:
        return "10-12"
    if n_qubits <= 16:
        return "14-16"
    return "18-20"


def build_condition_summary(cases: pd.DataFrame) -> pd.DataFrame:
    valid = cases[cases["valid_comparison"]].copy()
    valid["qubit_bin"] = valid["n_qubits"].map(qubit_bin)
    rows: list[dict[str, object]] = []
    groupings = {
        "basis": ["basis"],
        "qubit_bin": ["qubit_bin"],
        "basis_and_qubit_bin": ["basis", "qubit_bin"],
    }
    for group_type, columns in groupings.items():
        grouper: str | list[str] = columns[0] if len(columns) == 1 else columns
        for key, group in valid.groupby(grouper, sort=True):
            values = key if isinstance(key, tuple) else (key,)
            summary = summarize_subset("", group)
            summary["group_type"] = group_type
            summary["group_value"] = "/".join(map(str, values))
            rows.append(summary)
    condition = pd.DataFrame(rows)
    return condition[
        ["group_type", "group_value"]
        + [column for column in condition.columns if column not in {"group_type", "group_value", "group"}]
    ]


def save_figure(
    figure: plt.Figure,
    outdir: Path,
    stem: str,
    formats: list[str],
    dpi: int,
) -> list[Path]:
    paths: list[Path] = []
    for extension in formats:
        path = outdir / f"{stem}.{extension}"
        kwargs = {"dpi": dpi} if extension == "png" else {}
        figure.savefig(path, **kwargs)
        paths.append(path)
    return paths


def molecule_order(molecule_summary: pd.DataFrame) -> list[str]:
    summary = molecule_summary[molecule_summary["molecule"] != "ALL"].copy()
    return summary.sort_values(
        ["median_fermionic_advantage", "fermionic_win_fraction"],
        ascending=False,
    )["molecule"].tolist()


def plot_molecule_overview(
    cases: pd.DataFrame,
    molecule_summary: pd.DataFrame,
    outdir: Path,
    formats: list[str],
    dpi: int,
) -> list[Path]:
    valid = cases[cases["valid_comparison"]].copy()
    order = molecule_order(molecule_summary)
    y_lookup = {molecule: index for index, molecule in enumerate(order)}

    figure, (ax_ratio, ax_outcome) = plt.subplots(
        1,
        2,
        figsize=(12.5, 7.2),
        gridspec_kw={"width_ratios": [2.3, 1.15], "wspace": 0.12},
        sharey=True,
    )

    for molecule in order:
        group = valid[valid["molecule"] == molecule]
        y = y_lookup[molecule]
        for basis, basis_group in group.groupby("basis", sort=True):
            style = BASIS_STYLES.get(
                str(basis), {"marker": "^", "label": str(basis)}
            )
            offsets = np.linspace(-0.16, 0.16, len(basis_group))
            for offset, (_, row) in zip(offsets, basis_group.iterrows()):
                ax_ratio.scatter(
                    row["fermionic_advantage_factor"],
                    y + offset,
                    marker=style["marker"],
                    color=OUTCOME_COLORS[str(row["outcome"])],
                    edgecolor="white",
                    linewidth=0.4,
                    s=28,
                    alpha=0.88,
                    zorder=3,
                )
        values = group["fermionic_advantage_factor"]
        q25, median, q75 = values.quantile([0.25, 0.5, 0.75])
        ax_ratio.plot([q25, q75], [y, y], color="black", linewidth=2.1, zorder=4)
        ax_ratio.scatter(
            median,
            y,
            marker="D",
            color="black",
            edgecolor="white",
            linewidth=0.5,
            s=42,
            zorder=5,
        )

        counts = [
            int((group["outcome"] == "fermionic_win").sum()),
            int((group["outcome"] == "tie").sum()),
            int((group["outcome"] == "strongest_jw_win").sum()),
        ]
        left = 0.0
        for outcome, count in zip(
            ("fermionic_win", "tie", "strongest_jw_win"), counts
        ):
            fraction = count / len(group)
            ax_outcome.barh(
                y,
                fraction,
                left=left,
                height=0.62,
                color=OUTCOME_COLORS[outcome],
                edgecolor="white",
                linewidth=0.5,
            )
            left += fraction
        ax_outcome.text(
            1.02,
            y,
            f"{counts[0]}/{len(group)}",
            va="center",
            ha="left",
            fontsize=8,
        )

    labels = [MOLECULE_LABELS.get(molecule, molecule) for molecule in order]
    ax_ratio.set_yticks(range(len(order)), labels)
    ax_ratio.invert_yaxis()
    ax_ratio.set_xscale("log")
    all_values = valid["fermionic_advantage_factor"].to_numpy(dtype=float)
    lower = 10 ** math.floor(math.log10(all_values.min()))
    upper = 10 ** math.ceil(math.log10(all_values.max()))
    ax_ratio.set_xlim(lower, upper)
    ax_ratio.xaxis.set_major_formatter(LogFormatterMathtext())
    ax_ratio.axvline(1.0, color="black", linewidth=1.0, linestyle="--")
    ax_ratio.set_xlabel(
        r"Fermionic advantage $\min(E_{\mathrm{JW,signed}}, "
        r"E_{\mathrm{JW,mag}})/E_{\mathrm{F}}$"
    )
    ax_ratio.set_title("(a) Case distribution and median [IQR]")
    ax_ratio.grid(axis="y", visible=False)

    ax_outcome.set_xlim(0.0, 1.14)
    ax_outcome.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_outcome.set_xticklabels(["0", "25", "50", "75", "100"])
    ax_outcome.set_xlabel("Outcome share (%)")
    ax_outcome.set_title("(b) Wins / ties / losses")
    ax_outcome.grid(axis="y", visible=False)
    ax_outcome.tick_params(axis="y", left=False, labelleft=False)
    ax_outcome.text(
        1.02,
        -0.85,
        "F wins / n",
        ha="left",
        va="bottom",
        fontsize=8,
    )

    basis_handles = [
        Line2D(
            [0],
            [0],
            marker=style["marker"],
            color="none",
            markerfacecolor="#6c757d",
            markeredgecolor="white",
            markersize=6,
            label=style["label"],
        )
        for style in BASIS_STYLES.values()
    ]
    outcome_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=OUTCOME_COLORS[outcome],
            markersize=6,
            label=label,
        )
        for outcome, label in (
            ("fermionic_win", "Fermionic win"),
            ("tie", "Tie"),
            ("strongest_jw_win", "Strongest JW win"),
        )
    ]
    median_handle = Line2D(
        [0],
        [0],
        marker="D",
        color="black",
        markersize=6,
        label="Median (IQR bar)",
    )
    figure.legend(
        handles=basis_handles + outcome_handles + [median_handle],
        loc="lower center",
        ncol=6,
        frameon=False,
        bbox_to_anchor=(0.5, 0.045),
    )
    figure.suptitle(
        "Fermionic-aware ordering vs strongest deterministic JW baseline",
        y=1.01,
        fontsize=13,
    )
    figure.text(
        0.5,
        0.012,
        "Advantage > 1 favors fermionic-aware ordering; comparisons at or "
        "below 1e-12 are excluded.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    figure.subplots_adjust(bottom=0.2)
    paths = save_figure(
        figure,
        outdir,
        "fermionic_aware_molecule_performance",
        formats,
        dpi,
    )
    plt.close(figure)
    return paths


def plot_molecule_scaling(
    cases: pd.DataFrame,
    molecule_summary: pd.DataFrame,
    outdir: Path,
    formats: list[str],
    dpi: int,
) -> list[Path]:
    valid = cases[cases["valid_comparison"]].copy()
    order = molecule_order(molecule_summary)
    n_columns = 4
    n_rows = math.ceil(len(order) / n_columns)
    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(13.2, 3.15 * n_rows),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    minimum = valid["fermionic_advantage_factor"].min()
    maximum = valid["fermionic_advantage_factor"].max()
    y_lower = 10 ** math.floor(math.log10(minimum))
    y_upper = 10 ** math.ceil(math.log10(maximum))

    for index, molecule in enumerate(order):
        axis = axes.flat[index]
        group = valid[valid["molecule"] == molecule]
        for basis, basis_group in group.groupby("basis", sort=True):
            style = BASIS_STYLES.get(
                str(basis), {"marker": "^", "label": str(basis)}
            )
            axis.scatter(
                basis_group["n_qubits"],
                basis_group["fermionic_advantage_factor"],
                marker=style["marker"],
                s=30,
                alpha=0.9,
                label=style["label"],
            )
        summary_row = molecule_summary[
            molecule_summary["molecule"] == molecule
        ].iloc[0]
        axis.axhline(1.0, color="black", linewidth=0.9, linestyle="--")
        axis.set_yscale("log")
        axis.set_ylim(y_lower, y_upper)
        axis.set_xlim(3, 21)
        axis.set_xticks([4, 8, 12, 16, 20])
        axis.set_title(
            f"{MOLECULE_LABELS.get(molecule, molecule)}  "
            f"({int(summary_row['fermionic_wins'])}/"
            f"{int(summary_row['valid_comparisons'])} wins)"
        )
        axis.text(
            0.03,
            0.05,
            f"median={summary_row['median_fermionic_advantage']:.3g}",
            transform=axis.transAxes,
            fontsize=8,
            va="bottom",
        )

    for index in range(len(order), len(axes.flat)):
        axes.flat[index].axis("off")
    for row in range(n_rows):
        axes[row, 0].set_ylabel("Fermionic advantage")
    for column in range(n_columns):
        axes[-1, column].set_xlabel("Qubits in active space")

    handles = [
        Line2D(
            [0],
            [0],
            marker=style["marker"],
            color="none",
            markerfacecolor=f"C{index}",
            markersize=6,
            label=style["label"],
        )
        for index, style in enumerate(BASIS_STYLES.values())
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )
    figure.suptitle(
        "Active-space size dependence of fermionic-aware ordering",
        fontsize=13,
        y=1.0,
    )
    figure.subplots_adjust(bottom=0.08, hspace=0.32, wspace=0.15)
    paths = save_figure(
        figure,
        outdir,
        "fermionic_aware_active_space_performance",
        formats,
        dpi,
    )
    plt.close(figure)
    return paths


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    frame.to_csv(path, index=False, lineterminator="\n")


def main() -> None:
    args = parse_args()
    if args.error_floor < 0.0:
        raise ValueError("--error-floor must be non-negative")
    args.outdir.mkdir(parents=True, exist_ok=True)
    configure_plot_style()

    source = read_input(args.input)
    cases = build_case_table(source, args.error_floor)
    molecule_summary, basis_summary = build_summaries(cases)
    condition_summary = build_condition_summary(cases)

    write_csv(args.outdir / "fermionic_aware_case_performance.csv", cases)
    write_csv(
        args.outdir / "fermionic_aware_molecule_performance.csv",
        molecule_summary,
    )
    write_csv(
        args.outdir / "fermionic_aware_molecule_basis_performance.csv",
        basis_summary,
    )
    write_csv(
        args.outdir / "fermionic_aware_condition_performance.csv",
        condition_summary,
    )
    overview_paths = plot_molecule_overview(
        cases,
        molecule_summary,
        args.outdir,
        args.formats,
        args.dpi,
    )
    scaling_paths = plot_molecule_scaling(
        cases,
        molecule_summary,
        args.outdir,
        args.formats,
        args.dpi,
    )

    overall = molecule_summary[molecule_summary["molecule"] == "ALL"].iloc[0]
    print(f"Input cases: {len(cases)}")
    print(f"Fermionic-aware rows available: {int(cases.fermionic_available.sum())}")
    print(f"Valid comparisons: {int(overall.valid_comparisons)}")
    print(
        "Fermionic wins / ties / strongest-JW wins: "
        f"{int(overall.fermionic_wins)} / {int(overall.ties)} / "
        f"{int(overall.strongest_jw_wins)}"
    )
    print("Figures:")
    for path in overview_paths + scaling_paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
