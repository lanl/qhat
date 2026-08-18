#!/usr/bin/env python3
"""Plot first-order step scaling and the BCH ordering mechanism.

The script uses the existing final-validation CSVs and produces two figures:

1. ``scaling_plot``
   One panel per molecule/case, comparing signed fermionic and signed JW
   schedules across Trotter step counts.  A fitted log-log slope and an
   anchored ``r^-2`` guide show that ordering changes the leading prefactor,
   not the first-order asymptotic convergence rate.

2. ``bch_mechanism_plot``
   The left panel compares the squared leading BCH-norm advantage with the
   observed error advantage across active spaces.  The right panel performs
   the same comparison for the parent-structure ablation schedules.

3. ``random_order_comparison``
   For the four extension cases, deterministic schedules are overlaid on the
   two random-ablation error distributions.

By default, extension CSVs from ``instruction.md`` are included automatically
when they exist and ignored when they do not.  The script never modifies its
input CSVs.

Run from the repository root::

    python analysis/plot_fermionic_scaling_bch.py

Only NumPy, pandas, and Matplotlib are required.
"""

# ruff: noqa: E402  # Matplotlib cache/backend setup must precede its imports.

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable

# Keep Matplotlib usable in restricted and headless environments.  A user-set
# MPLCONFIGDIR still takes precedence.
_MPL_CACHE = Path(tempfile.gettempdir()) / "qhat-matplotlib-cache"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd


REFERENCE = "fermionic_signed_reference"
JW_SIGNED = "jw_signed_baseline"

FERMIONIC_MAGNITUDE = "fermionic_magnitude_reference"
ROUND_ROBIN = "signed_parent_descendants_round_robin"
RANDOM_BLOCKS = "signed_parent_blocks_randomized"
WITHIN_PARENT = "signed_parent_within_randomized"

DEFAULT_STEP_CSVS = (
    Path("analysis/final_validation_step_scaling_hgbs5.csv"),
    Path("analysis/final_validation_step_scaling_extension_hgbs5.csv"),
)
DEFAULT_ACTIVE_CSVS = (
    Path("analysis/final_validation_active_space_hgbs5.csv"),
)
DEFAULT_ABLATION_CSVS = (
    Path("analysis/fermionic_structure_ablation_hgbs5.csv"),
    Path("analysis/fermionic_structure_ablation_extension_hgbs5.csv"),
)

SCHEDULE_LABELS = {
    REFERENCE: "Fermionic signed",
    JW_SIGNED: "JW signed",
    FERMIONIC_MAGNITUDE: "Fermionic magnitude",
    ROUND_ROBIN: "Round-robin descendants",
    RANDOM_BLOCKS: "Random parent blocks",
    WITHIN_PARENT: "Within-parent shuffle",
}

MOLECULE_LABELS = {
    "B-B": r"B$_2$",
    "Be-Be": r"Be$_2$",
    "F-F": r"F$_2$",
    "Li-Li": r"Li$_2$",
    "O-O": r"O$_2$",
    "B-H": "BH",
    "Be-H": "BeH",
    "Li-H": "LiH",
    "BeH2": r"BeH$_2$",
    "H2O": r"H$_2$O",
    "NH3": r"NH$_3$",
    "CH4": r"CH$_4$",
}

SCHEDULE_STYLES = {
    REFERENCE: {"color": "#1f77b4", "marker": "o"},
    JW_SIGNED: {"color": "#ff7f0e", "marker": "s"},
    FERMIONIC_MAGNITUDE: {"color": "#2ca02c", "marker": "^"},
    ROUND_ROBIN: {"color": "#d62728", "marker": "D"},
    RANDOM_BLOCKS: {"color": "#9467bd", "marker": "v"},
    WITHIN_PARENT: {"color": "#8c564b", "marker": "P"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create publication-style first-order scaling and BCH-mechanism "
            "plots from the existing fermionic-ordering validation CSVs."
        )
    )
    parser.add_argument(
        "--step-csv",
        type=Path,
        nargs="+",
        default=list(DEFAULT_STEP_CSVS),
        help=(
            "Step-scaling CSV path(s). Missing default extension files are "
            "ignored."
        ),
    )
    parser.add_argument(
        "--active-space-csv",
        type=Path,
        nargs="+",
        default=list(DEFAULT_ACTIVE_CSVS),
        help="Active-space validation CSV path(s).",
    )
    parser.add_argument(
        "--ablation-csv",
        type=Path,
        nargs="+",
        default=list(DEFAULT_ABLATION_CSVS),
        help=(
            "Structure-ablation CSV path(s). Missing default extension files "
            "are ignored."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("analysis/fermionic_ordering_plots"),
        help="Output directory.",
    )
    parser.add_argument(
        "--error-column",
        choices=("one_minus_overlap", "state_infidelity"),
        default="one_minus_overlap",
        help=(
            "State-error metric used in both figures. Default: "
            "one_minus_overlap."
        ),
    )
    parser.add_argument(
        "--error-floor",
        type=float,
        default=1.0e-13,
        help="Values at or below this floor are excluded. Default: 1e-13.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG resolution. Default: 300.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "pdf", "svg"),
        default=["png", "pdf"],
        help="Figure formats to write. Default: png pdf.",
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
            "grid.alpha": 0.22,
            "grid.linewidth": 0.7,
            "lines.linewidth": 1.8,
            "lines.markersize": 5.5,
        }
    )


def read_success_csvs(paths: Iterable[Path], label: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists():
            print(f"INFO: {label} CSV not found; skipping: {path}")
            continue
        frame = pd.read_csv(path)
        if "status" in frame.columns:
            failed = frame[frame["status"] != "success"]
            if not failed.empty:
                print(
                    f"WARNING: excluding {len(failed)} non-success rows from "
                    f"{path}"
                )
            frame = frame[frame["status"] == "success"].copy()
        frame["source_csv"] = str(path)
        frames.append(frame)

    if not frames:
        raise FileNotFoundError(f"No readable {label} CSV was found.")
    return pd.concat(frames, ignore_index=True, sort=False)


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = sorted(set(columns).difference(df.columns))
    if missing:
        raise ValueError(f"{label} data is missing columns: {', '.join(missing)}")


def molecule_label(raw: object) -> str:
    value = str(raw)
    return MOLECULE_LABELS.get(value, value)


def fit_log_log(
    x_values: Iterable[float],
    y_values: Iterable[float],
) -> tuple[float, float, float, int]:
    x = np.asarray(list(x_values), dtype=float)
    y = np.asarray(list(y_values), dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    x = x[valid]
    y = y[valid]
    if x.size < 3:
        return math.nan, math.nan, math.nan, int(x.size)

    log_x = np.log10(x)
    log_y = np.log10(y)
    slope, intercept = np.polyfit(log_x, log_y, 1)
    fitted = intercept + slope * log_x
    residual = float(np.sum((log_y - fitted) ** 2))
    total = float(np.sum((log_y - np.mean(log_y)) ** 2))
    r_squared = 1.0 - residual / total if total > 0.0 else math.nan
    if np.std(log_x) > 0.0 and np.std(log_y) > 0.0:
        pearson = float(np.corrcoef(log_x, log_y)[0, 1])
    else:
        pearson = math.nan
    return float(slope), pearson, r_squared, int(x.size)


def save_figure(
    fig: plt.Figure,
    outdir: Path,
    stem: str,
    formats: Iterable[str],
    dpi: int,
) -> list[Path]:
    paths: list[Path] = []
    for extension in formats:
        path = outdir / f"{stem}.{extension}"
        kwargs = {"dpi": dpi} if extension == "png" else {}
        fig.savefig(path, **kwargs)
        paths.append(path)
    return paths


def plot_step_scaling(
    df: pd.DataFrame,
    *,
    error_column: str,
    error_floor: float,
    outdir: Path,
    formats: Iterable[str],
    dpi: int,
) -> tuple[pd.DataFrame, list[Path]]:
    required = (
        "case_id",
        "molecule",
        "basis",
        "active_occupied",
        "active_vacant",
        "n_qubits",
        "schedule",
        "trotter_steps",
        error_column,
    )
    require_columns(df, required, "step-scaling")

    data = df[df["schedule"].isin((REFERENCE, JW_SIGNED))].copy()
    data = data[
        np.isfinite(data["trotter_steps"])
        & np.isfinite(data[error_column])
        & (data["trotter_steps"] > 0)
        & (data[error_column] > error_floor)
    ].copy()
    data = data.sort_values(["case_id", "schedule", "trotter_steps"])
    data = data.drop_duplicates(
        subset=["case_id", "schedule", "trotter_steps"], keep="first"
    )
    if data.empty:
        raise ValueError("No valid step-scaling rows remain after filtering.")

    case_ids = list(dict.fromkeys(data["case_id"].astype(str)))
    n_cases = len(case_ids)
    n_columns = 2 if n_cases > 1 else 1
    n_rows = math.ceil(n_cases / n_columns)
    fig, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(6.0 * n_columns, 4.2 * n_rows),
        squeeze=False,
    )

    summary_rows: list[dict[str, object]] = []
    for panel_index, case_id in enumerate(case_ids):
        ax = axes.flat[panel_index]
        case = data[data["case_id"].astype(str) == case_id]
        metadata = case.iloc[0]
        all_steps = np.sort(case["trotter_steps"].unique().astype(float))

        guide_anchor: tuple[float, float] | None = None
        for schedule in (REFERENCE, JW_SIGNED):
            group = case[case["schedule"] == schedule].sort_values(
                "trotter_steps"
            )
            if group.empty:
                continue
            slope, pearson, r_squared, n_points = fit_log_log(
                group["trotter_steps"], group[error_column]
            )
            style = SCHEDULE_STYLES[schedule]
            label = f"{SCHEDULE_LABELS[schedule]}  ($s={slope:.2f}$)"
            ax.loglog(
                group["trotter_steps"],
                group[error_column],
                color=style["color"],
                marker=style["marker"],
                label=label,
            )
            if schedule == JW_SIGNED:
                first = group.iloc[0]
                guide_anchor = (
                    float(first["trotter_steps"]),
                    float(first[error_column]),
                )
            summary_rows.append(
                {
                    "case_id": case_id,
                    "molecule": metadata["molecule"],
                    "basis": metadata["basis"],
                    "active_occupied": int(metadata["active_occupied"]),
                    "active_vacant": int(metadata["active_vacant"]),
                    "n_qubits": int(metadata["n_qubits"]),
                    "schedule": schedule,
                    "error_column": error_column,
                    "n_fit_points": n_points,
                    "log_error_vs_log_steps_slope": slope,
                    "log_log_pearson": pearson,
                    "fit_r_squared": r_squared,
                    "step_min": int(group["trotter_steps"].min()),
                    "step_max": int(group["trotter_steps"].max()),
                    "error_at_step_min": float(group.iloc[0][error_column]),
                    "error_at_step_max": float(group.iloc[-1][error_column]),
                }
            )

        if guide_anchor is not None and all_steps.size >= 2:
            anchor_r, anchor_error = guide_anchor
            guide = anchor_error * (all_steps / anchor_r) ** -2.0
            ax.loglog(
                all_steps,
                guide,
                color="0.35",
                linestyle="--",
                linewidth=1.2,
                label=r"$r^{-2}$ guide",
            )

        active = (
            f"{int(metadata['active_occupied'])}+"
            f"{int(metadata['active_vacant'])}"
        )
        ax.set_title(
            f"{molecule_label(metadata['molecule'])}: "
            f"{int(metadata['n_qubits'])} qubits ({active})"
        )
        ax.set_xlabel("First-order Trotter steps, $r$")
        ax.set_ylabel(
            r"State error, $1-|\langle\psi_{\rm exact}|\psi_r\rangle|$"
            if error_column == "one_minus_overlap"
            else "State infidelity"
        )
        ax.set_xticks(all_steps)
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.legend(loc="best", frameon=False)

    for unused in range(n_cases, n_rows * n_columns):
        axes.flat[unused].set_visible(False)

    fig.suptitle(
        "First-order Trotter scaling: ordering changes the error prefactor",
        fontsize=14,
        y=1.005,
    )
    fig.tight_layout()
    output_paths = save_figure(
        fig, outdir, "scaling_plot", formats=formats, dpi=dpi
    )
    plt.close(fig)

    summary = pd.DataFrame(summary_rows).sort_values(
        ["molecule", "n_qubits", "schedule"]
    )
    summary.to_csv(outdir / "scaling_plot_fit_summary.csv", index=False)
    data.to_csv(outdir / "scaling_plot_data.csv", index=False)
    return summary, output_paths


def build_matched_bch_data(
    active_df: pd.DataFrame,
    ablation_df: pd.DataFrame,
    *,
    error_column: str,
    error_floor: float,
) -> pd.DataFrame:
    required = (
        "case_id",
        "molecule",
        "basis",
        "active_occupied",
        "active_vacant",
        "n_qubits",
        "schedule",
        "bch2_hf_state_norm",
        error_column,
    )
    require_columns(active_df, required, "active-space")
    require_columns(ablation_df, required, "ablation")

    active = active_df.copy()
    active["source_priority"] = 0
    ablation = ablation_df.copy()
    ablation["source_priority"] = 1
    combined = pd.concat([active, ablation], ignore_index=True, sort=False)
    combined = combined[combined["schedule"].isin((REFERENCE, JW_SIGNED))]
    combined = combined.sort_values("source_priority").drop_duplicates(
        subset=["case_id", "schedule"], keep="first"
    )

    key_columns = (
        "case_id",
        "molecule",
        "basis",
        "active_occupied",
        "active_vacant",
        "n_qubits",
    )
    records: list[dict[str, object]] = []
    for key, group in combined.groupby(list(key_columns), sort=True):
        schedules = group.set_index("schedule")
        if REFERENCE not in schedules.index or JW_SIGNED not in schedules.index:
            continue
        reference = schedules.loc[REFERENCE]
        jw = schedules.loc[JW_SIGNED]
        if isinstance(reference, pd.DataFrame):
            reference = reference.iloc[0]
        if isinstance(jw, pd.DataFrame):
            jw = jw.iloc[0]

        reference_error = float(reference[error_column])
        jw_error = float(jw[error_column])
        reference_bch = float(reference["bch2_hf_state_norm"])
        jw_bch = float(jw["bch2_hf_state_norm"])
        values = (reference_error, jw_error, reference_bch, jw_bch)
        if not all(np.isfinite(value) for value in values):
            continue
        if (
            reference_error <= error_floor
            or jw_error <= error_floor
            or reference_bch <= 0.0
            or jw_bch <= 0.0
        ):
            continue

        record = dict(zip(key_columns, key))
        bch_advantage = jw_bch / reference_bch
        record.update(
            {
                "fermionic_error": reference_error,
                "jw_error": jw_error,
                "fermionic_bch_norm": reference_bch,
                "jw_bch_norm": jw_bch,
                "bch_norm_advantage_jw_over_fermionic": bch_advantage,
                "predicted_error_advantage": bch_advantage**2,
                "observed_error_advantage": jw_error / reference_error,
            }
        )
        records.append(record)

    if not records:
        raise ValueError("No matched fermionic/JW BCH comparisons were found.")
    return pd.DataFrame(records).sort_values(["molecule", "n_qubits"])


def build_ablation_bch_data(
    df: pd.DataFrame,
    *,
    error_column: str,
    error_floor: float,
) -> pd.DataFrame:
    required = (
        "case_id",
        "molecule",
        "schedule",
        "sample_index",
        "bch_norm_ratio_to_signed_reference",
        error_column,
    )
    require_columns(df, required, "ablation")

    data = df.copy()
    data = data[
        np.isfinite(data["bch_norm_ratio_to_signed_reference"])
        & np.isfinite(data[error_column])
        & (data["bch_norm_ratio_to_signed_reference"] > 0.0)
        & (data[error_column] > error_floor)
    ].copy()

    references = (
        data[data["schedule"] == REFERENCE]
        .sort_values("source_csv")
        .drop_duplicates("case_id", keep="first")
        .set_index("case_id")[error_column]
    )
    data["reference_error"] = data["case_id"].map(references)
    data = data[
        np.isfinite(data["reference_error"])
        & (data["reference_error"] > error_floor)
    ].copy()
    data["predicted_error_ratio"] = (
        data["bch_norm_ratio_to_signed_reference"] ** 2
    )
    data["observed_error_ratio"] = (
        data[error_column] / data["reference_error"]
    )
    data = data[
        np.isfinite(data["predicted_error_ratio"])
        & np.isfinite(data["observed_error_ratio"])
        & (data["predicted_error_ratio"] > 0.0)
        & (data["observed_error_ratio"] > 0.0)
    ].copy()
    if data.empty:
        raise ValueError("No valid BCH ablation comparisons were found.")
    return data.sort_values(["molecule", "schedule", "sample_index"])


def equal_log_limits(ax: plt.Axes, x: np.ndarray, y: np.ndarray) -> None:
    values = np.concatenate([x, y])
    values = values[np.isfinite(values) & (values > 0.0)]
    if values.size == 0:
        return
    log_min = math.floor(float(np.log10(values.min())))
    log_max = math.ceil(float(np.log10(values.max())))
    if log_min == log_max:
        log_min -= 1
        log_max += 1
    lower = 10.0**log_min
    upper = 10.0**log_max
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.plot(
        [lower, upper],
        [lower, upper],
        color="0.35",
        linestyle="--",
        linewidth=1.2,
        zorder=1,
    )


def annotate_fit(
    ax: plt.Axes,
    x: Iterable[float],
    y: Iterable[float],
    *,
    exclude_reference: bool = False,
) -> tuple[float, float, float, int]:
    x_array = np.asarray(list(x), dtype=float)
    y_array = np.asarray(list(y), dtype=float)
    if exclude_reference:
        non_reference = ~(
            np.isclose(x_array, 1.0, rtol=1.0e-10, atol=1.0e-12)
            & np.isclose(y_array, 1.0, rtol=1.0e-10, atol=1.0e-12)
        )
        x_array = x_array[non_reference]
        y_array = y_array[non_reference]
    slope, pearson, r_squared, n_points = fit_log_log(x_array, y_array)
    text = (
        f"log-log Pearson $r={pearson:.3f}$\n"
        f"fit slope $={slope:.3f}$  ($n={n_points}$)"
    )
    ax.text(
        0.04,
        0.96,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    return slope, pearson, r_squared, n_points


def plot_bch_mechanism(
    matched: pd.DataFrame,
    ablation: pd.DataFrame,
    *,
    outdir: Path,
    formats: Iterable[str],
    dpi: int,
) -> tuple[pd.DataFrame, list[Path]]:
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.7))
    active_ax, ablation_ax = axes

    molecules = sorted(matched["molecule"].astype(str).unique())
    color_map = plt.get_cmap("tab10")
    marker_cycle = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h")
    for index, molecule in enumerate(molecules):
        group = matched[matched["molecule"].astype(str) == molecule]
        color = color_map(index % 10)
        marker = marker_cycle[index % len(marker_cycle)]
        active_ax.scatter(
            group["predicted_error_advantage"],
            group["observed_error_advantage"],
            s=54,
            color=color,
            marker=marker,
            edgecolor="white",
            linewidth=0.6,
            label=molecule_label(molecule),
            zorder=3,
        )
        for _, row in group.iterrows():
            active_ax.annotate(
                f"{int(row['n_qubits'])}q",
                (
                    float(row["predicted_error_advantage"]),
                    float(row["observed_error_advantage"]),
                ),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=7,
                color=color,
            )

    active_x = matched["predicted_error_advantage"].to_numpy(dtype=float)
    active_y = matched["observed_error_advantage"].to_numpy(dtype=float)
    active_ax.set_xscale("log")
    active_ax.set_yscale("log")
    equal_log_limits(active_ax, active_x, active_y)
    active_ax.axvline(1.0, color="0.65", linewidth=0.9)
    active_ax.axhline(1.0, color="0.65", linewidth=0.9)
    active_fit = annotate_fit(active_ax, active_x, active_y)
    active_ax.set_title("(a) Matched fermionic vs JW schedules")
    active_ax.set_xlabel(
        r"BCH-predicted advantage, "
        r"$(\|\Omega_2^{\rm JW}|\psi_0\rangle\|/"
        r"\|\Omega_2^{\rm F}|\psi_0\rangle\|)^2$"
    )
    active_ax.set_ylabel(
        r"Observed advantage, $\epsilon_{\rm JW}/\epsilon_{\rm F}$"
    )
    active_ax.legend(frameon=False, ncol=2, loc="lower right")

    for schedule in (
        REFERENCE,
        FERMIONIC_MAGNITUDE,
        JW_SIGNED,
        ROUND_ROBIN,
        RANDOM_BLOCKS,
        WITHIN_PARENT,
    ):
        group = ablation[ablation["schedule"] == schedule]
        if group.empty:
            continue
        style = SCHEDULE_STYLES[schedule]
        is_random = schedule in (RANDOM_BLOCKS, WITHIN_PARENT)
        ablation_ax.scatter(
            group["predicted_error_ratio"],
            group["observed_error_ratio"],
            s=28 if is_random else 48,
            color=style["color"],
            marker=style["marker"],
            alpha=0.58 if is_random else 0.9,
            edgecolor="none" if is_random else "white",
            linewidth=0.5,
            label=SCHEDULE_LABELS[schedule],
            zorder=2 if is_random else 3,
        )

    ablation_x = ablation["predicted_error_ratio"].to_numpy(dtype=float)
    ablation_y = ablation["observed_error_ratio"].to_numpy(dtype=float)
    ablation_ax.set_xscale("log")
    ablation_ax.set_yscale("log")
    equal_log_limits(ablation_ax, ablation_x, ablation_y)
    ablation_ax.axvline(1.0, color="0.65", linewidth=0.9)
    ablation_ax.axhline(1.0, color="0.65", linewidth=0.9)
    ablation_fit = annotate_fit(
        ablation_ax, ablation_x, ablation_y, exclude_reference=True
    )
    ablation_ax.set_title("(b) Fermionic parent-structure ablations")
    ablation_ax.set_xlabel(
        r"BCH-predicted ratio, "
        r"$(\|\Omega_2|\psi_0\rangle\|/"
        r"\|\Omega_2^{\rm F}|\psi_0\rangle\|)^2$"
    )
    ablation_ax.set_ylabel(
        r"Observed ratio, $\epsilon/\epsilon_{\rm F}$"
    )
    ablation_ax.legend(frameon=False, ncol=2, loc="lower right")

    diagonal_handle = Line2D(
        [0],
        [0],
        color="0.35",
        linestyle="--",
        linewidth=1.2,
        label="BCH prediction: observed = predicted",
    )
    fig.legend(
        handles=[diagonal_handle],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
    )
    fig.suptitle(
        "Leading state-dependent BCH error predicts ordering performance",
        fontsize=14,
        y=1.02,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    output_paths = save_figure(
        fig, outdir, "bch_mechanism_plot", formats=formats, dpi=dpi
    )
    plt.close(fig)

    fit_summary = pd.DataFrame(
        [
            {
                "panel": "matched_fermionic_vs_jw",
                "log_log_fit_slope": active_fit[0],
                "log_log_pearson": active_fit[1],
                "log_log_fit_r_squared": active_fit[2],
                "n_points": active_fit[3],
            },
            {
                "panel": "parent_structure_ablations",
                "log_log_fit_slope": ablation_fit[0],
                "log_log_pearson": ablation_fit[1],
                "log_log_fit_r_squared": ablation_fit[2],
                "n_points": ablation_fit[3],
            },
        ]
    )
    fit_summary.to_csv(outdir / "bch_mechanism_fit_summary.csv", index=False)
    matched.to_csv(outdir / "bch_matched_schedule_plot_data.csv", index=False)
    ablation.to_csv(outdir / "bch_ablation_plot_data.csv", index=False)
    return fit_summary, output_paths


def plot_random_order_comparison(
    ablation: pd.DataFrame,
    *,
    outdir: Path,
    formats: Iterable[str],
    dpi: int,
) -> list[Path]:
    target_molecules = ("BeH2", "H2O", "NH3", "O-O")
    data = ablation[ablation["molecule"].isin(target_molecules)].copy()
    if data.empty:
        print("INFO: no extension cases found; random-order plot skipped.")
        return []

    random_schedules = (RANDOM_BLOCKS, WITHIN_PARENT)
    deterministic_schedules = (
        REFERENCE,
        FERMIONIC_MAGNITUDE,
        JW_SIGNED,
        ROUND_ROBIN,
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), squeeze=False)

    for panel_index, molecule in enumerate(target_molecules):
        ax = axes.flat[panel_index]
        group = data[data["molecule"] == molecule]
        if group.empty:
            ax.set_visible(False)
            continue

        for position, schedule in enumerate(random_schedules, start=1):
            values = group.loc[
                group["schedule"] == schedule, "observed_error_ratio"
            ].to_numpy(dtype=float)
            if values.size == 0:
                continue
            box = ax.boxplot(
                [values],
                positions=[position],
                widths=0.36,
                showfliers=False,
                patch_artist=True,
                medianprops={"color": "0.2", "linewidth": 1.2},
                whiskerprops={"color": "0.45"},
                capprops={"color": "0.45"},
            )
            box["boxes"][0].set_facecolor(
                SCHEDULE_STYLES[schedule]["color"]
            )
            box["boxes"][0].set_alpha(0.18)
            offsets = np.linspace(-0.10, 0.10, values.size)
            ax.scatter(
                position + offsets,
                np.sort(values),
                s=26,
                color=SCHEDULE_STYLES[schedule]["color"],
                marker=SCHEDULE_STYLES[schedule]["marker"],
                alpha=0.72,
                edgecolor="none",
                zorder=3,
            )

        offsets = np.linspace(-0.18, 0.18, len(deterministic_schedules))
        for offset, schedule in zip(offsets, deterministic_schedules):
            values = group.loc[
                group["schedule"] == schedule, "observed_error_ratio"
            ].to_numpy(dtype=float)
            if values.size == 0:
                continue
            style = SCHEDULE_STYLES[schedule]
            ax.scatter(
                3.0 + offset,
                float(values[0]),
                s=64,
                color=style["color"],
                marker=style["marker"],
                edgecolor="white",
                linewidth=0.6,
                zorder=4,
            )

        ax.axhline(1.0, color="0.35", linestyle="--", linewidth=1.0)
        ax.set_yscale("log")
        ax.set_xticks((1, 2, 3))
        ax.set_xticklabels(("Random parent\nblocks", "Within-parent\nshuffle", "Deterministic"))
        ax.set_title(molecule_label(molecule))
        ax.set_ylabel(r"Error ratio, $\epsilon/\epsilon_{\rm F}$")
        ax.grid(axis="x", visible=False)

    handles = [
        Line2D(
            [0],
            [0],
            color="none",
            marker=SCHEDULE_STYLES[schedule]["marker"],
            markerfacecolor=SCHEDULE_STYLES[schedule]["color"],
            markeredgecolor="white",
            markersize=7,
            label=SCHEDULE_LABELS[schedule],
        )
        for schedule in deterministic_schedules
    ]
    fig.legend(
        handles=handles,
        frameon=False,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
    )
    fig.suptitle(
        "Deterministic schedules within random-order error distributions",
        fontsize=14,
        y=1.01,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    output_paths = save_figure(
        fig, outdir, "random_order_comparison", formats=formats, dpi=dpi
    )
    plt.close(fig)
    return output_paths


def main() -> None:
    args = parse_args()
    if args.error_floor <= 0.0:
        raise ValueError("--error-floor must be positive.")
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive.")

    configure_plot_style()
    args.outdir.mkdir(parents=True, exist_ok=True)

    step_df = read_success_csvs(args.step_csv, "step-scaling")
    active_df = read_success_csvs(args.active_space_csv, "active-space")
    ablation_df = read_success_csvs(args.ablation_csv, "ablation")

    scaling_summary, scaling_paths = plot_step_scaling(
        step_df,
        error_column=args.error_column,
        error_floor=args.error_floor,
        outdir=args.outdir,
        formats=args.formats,
        dpi=args.dpi,
    )
    matched = build_matched_bch_data(
        active_df,
        ablation_df,
        error_column=args.error_column,
        error_floor=args.error_floor,
    )
    ablation = build_ablation_bch_data(
        ablation_df,
        error_column=args.error_column,
        error_floor=args.error_floor,
    )
    mechanism_summary, mechanism_paths = plot_bch_mechanism(
        matched,
        ablation,
        outdir=args.outdir,
        formats=args.formats,
        dpi=args.dpi,
    )
    random_paths = plot_random_order_comparison(
        ablation,
        outdir=args.outdir,
        formats=args.formats,
        dpi=args.dpi,
    )

    print()
    print("Scaling fits:")
    print(
        scaling_summary[
            [
                "molecule",
                "n_qubits",
                "schedule",
                "log_error_vs_log_steps_slope",
                "fit_r_squared",
            ]
        ].to_string(index=False)
    )
    print()
    print("BCH mechanism fits:")
    print(mechanism_summary.to_string(index=False))
    print()
    print("Created figures:")
    for path in scaling_paths + mechanism_paths + random_paths:
        print(f"  {path}")
    print()
    print(f"Plot data and fit summaries: {args.outdir}")


if __name__ == "__main__":
    main()
