#!/usr/bin/env python3
"""Generate standalone BCH presentation panels from saved validation CSVs.

This script is deliberately post-processing only.  It never invokes the BCH
or Trotter benchmarks; every plotted value is reconstructed from committed CSV
artifacts.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
from xml.sax.saxutils import escape

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PRIMARY_BASELINE = "jw_magnitude_baseline"
PROPOSED_ORDERING = "fermionic_signed_reference"
RANDOM_PARENT_BLOCKS = "signed_parent_blocks_randomized"
WITHIN_PARENT_SHUFFLE = "signed_parent_within_randomized"
PARENT_ABLATION_SCHEDULES = (
    RANDOM_PARENT_BLOCKS,
    WITHIN_PARENT_SHUFFLE,
)
TROTTER_STEPS = (50, 100, 200)

DEFAULT_CASE_SUMMARY = Path(
    "analysis/cancellation_hypothesis_validation/full_analysis/case_summary.csv"
)
DEFAULT_MATCHED_PAIRS = Path(
    "analysis/cancellation_hypothesis_validation/full_analysis/"
    "matched_pair_summary.csv"
)
DEFAULT_PARENT_ABLATION = Path(
    "analysis/cancellation_hypothesis_validation/pilot_analysis/"
    "primary_aggregated_deltas.csv"
)
DEFAULT_OUTPUT_DIR = Path("analysis/presentation_bch_panels")

EXPECTED_HELDOUT_R = 0.903
EXPECTED_HELDOUT_P = 5.28e-8
EXPECTED_HELDOUT_DIRECTION = (17, 20)
EXPECTED_MATCHED_R = 0.916
EXPECTED_MATCHED_P = 1.96e-4
EXPECTED_MATCHED_DIRECTION = (9, 10)
EXPECTED_RANDOM_MEDIANS = (5.947, 99.12)
EXPECTED_WITHIN_MEDIANS = (1.0, 1.0)

OUTCOME_COLORS = {
    "favorable": "#D55E00",
    "negative_control": "#0072B2",
}
OUTCOME_LABELS = {
    "favorable": "Favorable",
    "negative_control": "Control",
}
CONCORDANCE_COLORS = {
    True: "#009E73",
    False: "#CC79A7",
}
SCHEDULE_COLORS = {
    RANDOM_PARENT_BLOCKS: "#D55E00",
    WITHIN_PARENT_SHUFFLE: "#0072B2",
}
SCHEDULE_MARKERS = {
    RANDOM_PARENT_BLOCKS: "o",
    WITHIN_PARENT_SHUFFLE: "D",
}
SCHEDULE_LABELS = {
    RANDOM_PARENT_BLOCKS: "Random parent blocks",
    WITHIN_PARENT_SHUFFLE: "Within-parent shuffle",
}

HELDOUT_LABEL_OFFSETS = {
    "Li-Li_2.66_hgbs-5_as-006-006": (36, 30),
    "F-F_1.28_hgbs-5_as-012-002": (9, -20),
    "Be-H_1.47_hgbs-5_as-005-013": (32, -40),
    "CH4_s-1.00_hgbs-5_as-008-010": (9, -19),
    "H2O_s-1.00_hgbs-5_as-008-012": (-95, 9),
    "B-B_1.70_hgbs-5_as-010-010": (9, -19),
}
RANDOM_LABEL_OFFSETS = {
    "Be-H": (9, -18),
    "H2O": (9, -19),
    "Li-H": (9, 8),
    "Li-Li": (-66, 7),
}


@dataclass(frozen=True)
class CorrelationSummary:
    pearson_r: float
    p_value: float
    direction_count: int
    direction_total: int


@dataclass(frozen=True)
class PresentationResults:
    heldout: CorrelationSummary
    matched_pairs: CorrelationSummary
    parent_medians: dict[str, tuple[float, float]]
    output_files: tuple[Path, ...]


def repository_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_from_root(root: Path, path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def require_columns(
    frame: pd.DataFrame,
    columns: Iterable[str],
    source: str,
) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def require_finite(frame: pd.DataFrame, columns: Sequence[str], source: str) -> None:
    values = frame[list(columns)].to_numpy(dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{source} contains non-finite plotting values.")


def correlation_summary(x: Sequence[float], y: Sequence[float]) -> CorrelationSummary:
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    if len(x_array) != len(y_array) or len(x_array) < 3:
        raise ValueError("Pearson correlation requires at least three paired values.")
    result = stats.pearsonr(x_array, y_array)
    direction = int(np.count_nonzero(x_array * y_array > 0.0))
    return CorrelationSummary(
        pearson_r=float(result.statistic),
        p_value=float(result.pvalue),
        direction_count=direction,
        direction_total=len(x_array),
    )


def assert_close(
    actual: float,
    expected: float,
    *,
    label: str,
    atol: float = 0.0,
    rtol: float = 0.0,
) -> None:
    if not math.isclose(actual, expected, abs_tol=atol, rel_tol=rtol):
        raise AssertionError(
            f"{label} mismatch: actual={actual:.12g}, expected={expected:.12g}, "
            f"atol={atol:g}, rtol={rtol:g}"
        )


def validate_expected_statistics(
    heldout: CorrelationSummary,
    matched: CorrelationSummary,
    parent_medians: dict[str, tuple[float, float]],
) -> None:
    assert_close(
        heldout.pearson_r,
        EXPECTED_HELDOUT_R,
        label="held-out Pearson r",
        atol=5.0e-4,
    )
    assert_close(
        heldout.p_value,
        EXPECTED_HELDOUT_P,
        label="held-out Pearson p",
        rtol=0.01,
    )
    if (heldout.direction_count, heldout.direction_total) != EXPECTED_HELDOUT_DIRECTION:
        raise AssertionError("Held-out direction agreement is not 17/20.")

    assert_close(
        matched.pearson_r,
        EXPECTED_MATCHED_R,
        label="matched-pair Pearson r",
        atol=5.0e-4,
    )
    assert_close(
        matched.p_value,
        EXPECTED_MATCHED_P,
        label="matched-pair Pearson p",
        rtol=0.01,
    )
    if (matched.direction_count, matched.direction_total) != EXPECTED_MATCHED_DIRECTION:
        raise AssertionError("Matched-pair direction agreement is not 9/10.")

    random_cancellation, random_error = parent_medians[RANDOM_PARENT_BLOCKS]
    within_cancellation, within_error = parent_medians[WITHIN_PARENT_SHUFFLE]
    assert_close(
        random_cancellation,
        EXPECTED_RANDOM_MEDIANS[0],
        label="random-parent median BCH ratio",
        atol=0.001,
    )
    assert_close(
        random_error,
        EXPECTED_RANDOM_MEDIANS[1],
        label="random-parent median error ratio",
        atol=0.01,
    )
    assert_close(
        within_cancellation,
        EXPECTED_WITHIN_MEDIANS[0],
        label="within-parent median BCH ratio",
        atol=1.0e-10,
    )
    assert_close(
        within_error,
        EXPECTED_WITHIN_MEDIANS[1],
        label="within-parent median error ratio",
        atol=2.0e-6,
    )


def compact_case_label(case_id: str, n_qubits: int) -> str:
    molecule = case_id.split("_")[0]
    replacements = {
        "Li-Li": r"Li$_2$",
        "F-F": r"F$_2$",
        "B-B": r"B$_2$",
        "H2O": r"H$_2$O",
        "CH4": r"CH$_4$",
        "Be-H": "BeH",
    }
    return f"{replacements.get(molecule, molecule)} ({n_qubits}q)"


def prepare_heldout_cases(case_summary: pd.DataFrame) -> tuple[pd.DataFrame, CorrelationSummary]:
    source = "case_summary.csv"
    required = (
        "case_id",
        "molecule",
        "expected_outcome",
        "n_qubits",
        "bch_discovery_case",
        "signed_bch_cancellation_ratio",
        "fresh_signed_one_minus_overlap",
        "fresh_jw_magnitude_one_minus_overlap",
        "fresh_jw_magnitude_bch_cancellation_ratio",
        "fresh_jw_magnitude_to_signed_bch_cancellation_ratio",
        "fresh_jw_magnitude_to_signed_advantage",
    )
    require_columns(case_summary, required, source)
    if len(case_summary) != 20:
        raise AssertionError(f"Expected 20 held-out cases, found {len(case_summary)}.")
    discovery = case_summary["bch_discovery_case"].astype(str).str.lower()
    if discovery.isin(("true", "1", "yes")).any():
        raise AssertionError("The presentation panel contains a BCH discovery case.")
    outcomes = set(case_summary["expected_outcome"].astype(str))
    if outcomes != set(OUTCOME_COLORS):
        raise AssertionError(f"Unexpected held-out outcome labels: {sorted(outcomes)}")

    numeric = (
        "signed_bch_cancellation_ratio",
        "fresh_signed_one_minus_overlap",
        "fresh_jw_magnitude_one_minus_overlap",
        "fresh_jw_magnitude_bch_cancellation_ratio",
        "fresh_jw_magnitude_to_signed_bch_cancellation_ratio",
        "fresh_jw_magnitude_to_signed_advantage",
    )
    require_finite(case_summary, numeric, source)
    if not np.all(case_summary[list(numeric)].to_numpy(dtype=float) > 0.0):
        raise ValueError(f"{source} contains non-positive ratios or errors.")

    expected_error_ratio = (
        case_summary["fresh_jw_magnitude_one_minus_overlap"].to_numpy(float)
        / case_summary["fresh_signed_one_minus_overlap"].to_numpy(float)
    )
    expected_bch_ratio = (
        case_summary["fresh_jw_magnitude_bch_cancellation_ratio"].to_numpy(float)
        / case_summary["signed_bch_cancellation_ratio"].to_numpy(float)
    )
    np.testing.assert_allclose(
        case_summary["fresh_jw_magnitude_to_signed_advantage"],
        expected_error_ratio,
        rtol=1.0e-11,
        atol=0.0,
        err_msg="Primary error comparator is not JW magnitude versus fermionic signed.",
    )
    np.testing.assert_allclose(
        case_summary["fresh_jw_magnitude_to_signed_bch_cancellation_ratio"],
        expected_bch_ratio,
        rtol=1.0e-11,
        atol=0.0,
        err_msg="Primary BCH comparator is not JW magnitude versus fermionic signed.",
    )

    data = case_summary[
        ["case_id", "molecule", "expected_outcome", "n_qubits"]
    ].copy()
    data["baseline_schedule"] = PRIMARY_BASELINE
    data["proposed_schedule"] = PROPOSED_ORDERING
    data["bch_ratio_jwmag_to_fermionic_signed"] = expected_bch_ratio
    data["error_ratio_jwmag_to_fermionic_signed"] = expected_error_ratio
    data["log10_relative_cancellation"] = np.log10(expected_bch_ratio)
    data["log10_fresh_advantage"] = np.log10(expected_error_ratio)
    data["direction_concordant"] = (
        data["log10_relative_cancellation"] * data["log10_fresh_advantage"] > 0.0
    )
    data["point_label"] = ""
    for index, row in data.iterrows():
        if row["case_id"] in HELDOUT_LABEL_OFFSETS:
            data.at[index, "point_label"] = compact_case_label(
                str(row["case_id"]), int(row["n_qubits"])
            )

    if set(data["baseline_schedule"]) != {PRIMARY_BASELINE}:
        raise AssertionError("A non-primary baseline entered the held-out panel.")
    if set(data["proposed_schedule"]) != {PROPOSED_ORDERING}:
        raise AssertionError("A non-proposed ordering entered the held-out panel.")

    summary = correlation_summary(
        data["log10_relative_cancellation"],
        data["log10_fresh_advantage"],
    )
    return data.sort_values("case_id").reset_index(drop=True), summary


def prepare_matched_pairs(
    matched_pair_summary: pd.DataFrame,
    heldout: pd.DataFrame,
) -> tuple[pd.DataFrame, CorrelationSummary]:
    source = "matched_pair_summary.csv"
    required = (
        "matched_pair",
        "n_qubits",
        "favorable_case_id",
        "negative_control_case_id",
        "delta_log10_relative_cancellation",
        "delta_log10_fresh_advantage",
        "direction_concordant",
    )
    require_columns(matched_pair_summary, required, source)
    if len(matched_pair_summary) != 10:
        raise AssertionError(
            f"Expected 10 active-space-matched pairs, found {len(matched_pair_summary)}."
        )
    require_finite(
        matched_pair_summary,
        ("delta_log10_relative_cancellation", "delta_log10_fresh_advantage"),
        source,
    )

    heldout_by_case = heldout.set_index("case_id")
    rows: list[dict[str, object]] = []
    for row in matched_pair_summary.itertuples(index=False):
        favorable = heldout_by_case.loc[str(row.favorable_case_id)]
        control = heldout_by_case.loc[str(row.negative_control_case_id)]
        if favorable["expected_outcome"] != "favorable":
            raise AssertionError(f"Pair {row.matched_pair} has the wrong favorable case.")
        if control["expected_outcome"] != "negative_control":
            raise AssertionError(f"Pair {row.matched_pair} has the wrong control case.")
        x = float(
            favorable["log10_relative_cancellation"]
            - control["log10_relative_cancellation"]
        )
        y = float(favorable["log10_fresh_advantage"] - control["log10_fresh_advantage"])
        assert_close(
            float(row.delta_log10_relative_cancellation),
            x,
            label=f"{row.matched_pair} cancellation delta",
            atol=1.0e-11,
        )
        assert_close(
            float(row.delta_log10_fresh_advantage),
            y,
            label=f"{row.matched_pair} error delta",
            atol=1.0e-11,
        )
        concordant = x * y > 0.0
        if bool(row.direction_concordant) != concordant:
            raise AssertionError(f"Pair {row.matched_pair} has a stale direction flag.")
        rows.append(
            {
                "matched_pair": row.matched_pair,
                "n_qubits": int(row.n_qubits),
                "favorable_case_id": row.favorable_case_id,
                "negative_control_case_id": row.negative_control_case_id,
                "baseline_schedule": PRIMARY_BASELINE,
                "proposed_schedule": PROPOSED_ORDERING,
                "delta_log10_relative_cancellation": x,
                "delta_log10_fresh_advantage": y,
                "direction_concordant": concordant,
                "point_label": (
                    f"{row.matched_pair}\nH$_2$O vs F$_2$"
                    if not concordant
                    else ""
                ),
            }
        )

    data = pd.DataFrame(rows).sort_values("matched_pair").reset_index(drop=True)
    if set(data["baseline_schedule"]) != {PRIMARY_BASELINE}:
        raise AssertionError("A non-primary baseline entered the matched-pair panel.")
    if set(data["proposed_schedule"]) != {PROPOSED_ORDERING}:
        raise AssertionError("A non-proposed ordering entered the matched-pair panel.")
    summary = correlation_summary(
        data["delta_log10_relative_cancellation"],
        data["delta_log10_fresh_advantage"],
    )
    return data, summary


def prepare_parent_ablation(
    aggregated_deltas: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, tuple[float, float]]]:
    source = "primary_aggregated_deltas.csv"
    required = (
        "case_id",
        "molecule",
        "schedule",
        "trotter_steps",
        "evolution_time",
        "delta_log10_cancellation",
        "delta_log10_error",
    )
    require_columns(aggregated_deltas, required, source)
    selected = aggregated_deltas[
        aggregated_deltas["schedule"].isin(PARENT_ABLATION_SCHEDULES)
    ].copy()
    if set(selected["schedule"].astype(str)) != set(PARENT_ABLATION_SCHEDULES):
        raise AssertionError("The parent-ablation input is missing a required schedule.")
    if set(selected["trotter_steps"].astype(int)) != set(TROTTER_STEPS):
        raise AssertionError("Parent ablation must use exactly 50/100/200 Trotter steps.")
    if not np.allclose(selected["evolution_time"].to_numpy(float), 1.0):
        raise AssertionError("Parent ablation mixes evolution times.")
    require_finite(
        selected,
        ("delta_log10_cancellation", "delta_log10_error"),
        source,
    )

    for key, group in selected.groupby(["case_id", "schedule"], sort=True):
        steps = tuple(sorted(group["trotter_steps"].astype(int)))
        if steps != TROTTER_STEPS:
            raise AssertionError(f"{key} does not contain exactly 50/100/200 steps: {steps}")

    collapsed = (
        selected.groupby(["case_id", "molecule", "schedule"], as_index=False)
        .agg(
            delta_log10_cancellation=("delta_log10_cancellation", "median"),
            delta_log10_error=("delta_log10_error", "median"),
            step_rows=("trotter_steps", "size"),
        )
        .sort_values(["schedule", "case_id"])
        .reset_index(drop=True)
    )
    if collapsed["case_id"].nunique() != 4 or len(collapsed) != 8:
        raise AssertionError("Expected four pilot cases for each of two ablation schedules.")
    if set(collapsed["schedule"]) != set(PARENT_ABLATION_SCHEDULES):
        raise AssertionError("An unintended schedule entered the parent-ablation panel.")
    if not (collapsed["step_rows"] == 3).all():
        raise AssertionError("Each parent-ablation point must collapse three step rows.")

    collapsed["reference_schedule"] = PROPOSED_ORDERING
    collapsed["steps_collapsed"] = "50/100/200"
    collapsed["cancellation_ratio_to_signed"] = 10.0 ** collapsed[
        "delta_log10_cancellation"
    ]
    collapsed["error_ratio_to_signed"] = 10.0 ** collapsed["delta_log10_error"]
    collapsed["point_label"] = np.where(
        collapsed["schedule"] == RANDOM_PARENT_BLOCKS,
        collapsed["molecule"],
        "",
    )

    medians: dict[str, tuple[float, float]] = {}
    for schedule in PARENT_ABLATION_SCHEDULES:
        group = collapsed[collapsed["schedule"] == schedule]
        medians[schedule] = (
            float(group["cancellation_ratio_to_signed"].median()),
            float(group["error_ratio_to_signed"].median()),
        )
    return collapsed, medians


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.labelsize": 17,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 13,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 1.1,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def finish_axes(ax: plt.Axes, x: Sequence[float], y: Sequence[float]) -> None:
    x_values = np.append(np.asarray(x, dtype=float), 0.0)
    y_values = np.append(np.asarray(y, dtype=float), 0.0)
    x_span = max(float(np.ptp(x_values)), 0.2)
    y_span = max(float(np.ptp(y_values)), 0.2)
    ax.set_xlim(float(x_values.min() - 0.12 * x_span), float(x_values.max() + 0.12 * x_span))
    ax.set_ylim(float(y_values.min() - 0.12 * y_span), float(y_values.max() + 0.12 * y_span))
    ax.axvline(0.0, color="0.35", linewidth=1.25, linestyle="--", zorder=0)
    ax.axhline(0.0, color="0.35", linewidth=1.25, linestyle="--", zorder=0)
    ax.grid(True, color="0.88", linewidth=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=5, width=1.0)


def save_png_pdf(fig: plt.Figure, output_dir: Path, stem: str) -> tuple[Path, Path]:
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    return png, pdf


def format_p_value(value: float) -> str:
    if value < 1.0e-3:
        return f"{value:.2e}"
    return f"{value:.3f}"


def plot_heldout_cases(
    data: pd.DataFrame,
    summary: CorrelationSummary,
    output_dir: Path,
) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for outcome in ("negative_control", "favorable"):
        group = data[data["expected_outcome"] == outcome]
        ax.scatter(
            group["log10_relative_cancellation"],
            group["log10_fresh_advantage"],
            s=100,
            color=OUTCOME_COLORS[outcome],
            edgecolor="white",
            linewidth=0.9,
            label=OUTCOME_LABELS[outcome],
            zorder=3,
        )
    for row in data.itertuples(index=False):
        if not row.point_label:
            continue
        offset = HELDOUT_LABEL_OFFSETS[str(row.case_id)]
        ax.annotate(
            row.point_label,
            (row.log10_relative_cancellation, row.log10_fresh_advantage),
            xytext=offset,
            textcoords="offset points",
            fontsize=11.5,
            color="0.18",
            arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.7},
            zorder=4,
        )
    finish_axes(ax, data["log10_relative_cancellation"], data["log10_fresh_advantage"])
    ax.set_xlabel(
        r"$\log_{10}(R_{\mathrm{BCH,JWmag}}/R_{\mathrm{BCH,F\!\!-signed}})$"
    )
    ax.set_ylabel(
        r"$\log_{10}(\epsilon_{\mathrm{JWmag}}/\epsilon_{\mathrm{F\!\!-signed}})$"
    )
    ax.legend(frameon=False, loc="lower right")
    ax.text(
        0.04,
        0.96,
        (
            f"Pearson $r$ = {summary.pearson_r:.3f}\n"
            f"$p$ = {format_p_value(summary.p_value)}\n"
            f"Direction agreement = {summary.direction_count}/{summary.direction_total}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        bbox={"facecolor": "white", "edgecolor": "0.82", "pad": 5.0},
    )
    fig.tight_layout()
    paths = save_png_pdf(fig, output_dir, "slide2_heldout_cases")
    plt.close(fig)
    return paths


def plot_matched_pairs(
    data: pd.DataFrame,
    summary: CorrelationSummary,
    output_dir: Path,
) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for concordant in (True, False):
        group = data[data["direction_concordant"] == concordant]
        ax.scatter(
            group["delta_log10_relative_cancellation"],
            group["delta_log10_fresh_advantage"],
            s=110,
            marker="o" if concordant else "X",
            color=CONCORDANCE_COLORS[concordant],
            edgecolor="white",
            linewidth=0.9,
            label="Concordant" if concordant else "Discordant",
            zorder=3,
        )
    discordant = data[~data["direction_concordant"]]
    for row in discordant.itertuples(index=False):
        ax.annotate(
            row.point_label,
            (row.delta_log10_relative_cancellation, row.delta_log10_fresh_advantage),
            xytext=(14, 12),
            textcoords="offset points",
            fontsize=12,
            color="0.18",
            arrowprops={"arrowstyle": "-", "color": "0.5", "lw": 0.8},
            zorder=4,
        )
    finish_axes(
        ax,
        data["delta_log10_relative_cancellation"],
        data["delta_log10_fresh_advantage"],
    )
    ax.set_xlabel(r"$\Delta\log_{10}$ relative BCH cancellation")
    ax.set_ylabel(r"$\Delta\log_{10}$ JWmag / F-signed error advantage")
    ax.legend(frameon=False, loc="lower right")
    ax.text(
        0.04,
        0.96,
        (
            f"Pearson $r$ = {summary.pearson_r:.3f}\n"
            f"$p$ = {format_p_value(summary.p_value)}\n"
            f"Direction agreement = {summary.direction_count}/{summary.direction_total}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        bbox={"facecolor": "white", "edgecolor": "0.82", "pad": 5.0},
    )
    fig.tight_layout()
    paths = save_png_pdf(fig, output_dir, "slide2_matched_pairs")
    plt.close(fig)
    return paths


def plot_parent_ablation(
    data: pd.DataFrame,
    medians: dict[str, tuple[float, float]],
    output_dir: Path,
) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    for schedule in PARENT_ABLATION_SCHEDULES:
        group = data[data["schedule"] == schedule]
        ax.scatter(
            group["delta_log10_cancellation"],
            group["delta_log10_error"],
            s=112,
            marker=SCHEDULE_MARKERS[schedule],
            color=SCHEDULE_COLORS[schedule],
            edgecolor="white",
            linewidth=0.9,
            label=SCHEDULE_LABELS[schedule],
            zorder=3,
        )
    random_rows = data[data["schedule"] == RANDOM_PARENT_BLOCKS]
    for row in random_rows.itertuples(index=False):
        ax.annotate(
            str(row.point_label),
            (row.delta_log10_cancellation, row.delta_log10_error),
            xytext=RANDOM_LABEL_OFFSETS[str(row.molecule)],
            textcoords="offset points",
            fontsize=12,
            color="0.18",
            arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.7},
            zorder=4,
        )
    finish_axes(ax, data["delta_log10_cancellation"], data["delta_log10_error"])
    ax.set_xlabel(r"$\Delta\log_{10}$ BCH cancellation ratio")
    ax.set_ylabel(r"$\Delta\log_{10}$ Trotter error")
    ax.legend(frameon=False, loc="upper left")

    random_bch, random_error = medians[RANDOM_PARENT_BLOCKS]
    within_bch, within_error = medians[WITHIN_PARENT_SHUFFLE]
    ax.text(
        0.97,
        0.05,
        (
            "Median across cases\n"
            f"Random blocks: $R_{{\\mathrm{{BCH}}}}$ = {random_bch:.3f}×  |  "
            f"error = {random_error:.2f}×\n"
            f"Within-parent: $R_{{\\mathrm{{BCH}}}}$ = {within_bch:.2f}×  |  "
            f"error = {within_error:.2f}×"
        ),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11.5,
        bbox={"facecolor": "white", "edgecolor": "0.82", "pad": 5.0},
    )
    fig.tight_layout()
    paths = save_png_pdf(fig, output_dir, "slide3_parent_block_ablation")
    plt.close(fig)
    return paths


def svg_text(
    x: float,
    y: float,
    text: str,
    *,
    size: int,
    weight: str = "normal",
    anchor: str = "start",
    fill: str = "#202124",
) -> str:
    return (
        f'<text x="{x}" y="{y}" font-family="DejaVu Sans, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
        f'fill="{fill}">{escape(text)}</text>'
    )


def draw_parent_block(
    x: float,
    y: float,
    parent: str,
    descendants: Sequence[str],
    *,
    shuffled: bool,
) -> list[str]:
    colors = {"P1": "#56B4E9", "P2": "#E69F00", "P3": "#009E73"}
    width = 250
    height = 108
    elements = [
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="8" '
        f'fill="{colors[parent]}" fill-opacity="0.18" stroke="{colors[parent]}" '
        'stroke-width="3"/>',
        svg_text(
            x + 18,
            y + 31,
            f"{parent} ↻" if shuffled else parent,
            size=23,
            weight="bold",
        ),
    ]
    cell_y = y + 48
    cell_width = 61
    for index, descendant in enumerate(descendants):
        cell_x = x + 18 + index * 72
        elements.append(
            f'<rect x="{cell_x}" y="{cell_y}" width="{cell_width}" height="42" '
            'rx="4" fill="#FFFFFF" stroke="#606368" stroke-width="1.5"/>'
        )
        elements.append(
            svg_text(
                cell_x + cell_width / 2,
                cell_y + 28,
                descendant,
                size=18,
                anchor="middle",
            )
        )
    return elements


def write_parent_order_schematic(output_dir: Path) -> tuple[Path, Path]:
    svg_path = output_dir / "slide3_parent_order_schematic.svg"
    data_path = output_dir / "slide3_parent_order_schematic_data.csv"
    width, height = 1450, 610
    label_x = 34
    block_x = (420, 690, 960)
    row_y = (38, 228, 418)
    rows = (
        (
            "Intact signed order",
            ("P1", "P2", "P3"),
            (("d1", "d2", "d3"),) * 3,
            (False, False, False),
            ("Reference", "Parent + descendant order intact"),
        ),
        (
            "Parent-block intervention",
            ("P3", "P1", "P2"),
            (("d1", "d2", "d3"),) * 3,
            (False, False, False),
            ("Whole parent blocks permuted", "Descendant order + contiguity preserved"),
        ),
        (
            "Within-parent control",
            ("P1", "P2", "P3"),
            (("d3", "d1", "d2"), ("d2", "d3", "d1"), ("d2", "d1", "d3")),
            (True, True, True),
            ("Parent-block order fixed", "Descendants shuffled within each block"),
        ),
    )
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
    ]
    data_rows: list[dict[str, object]] = []
    for row_index, (label, parents, descendants, shuffled, notes) in enumerate(rows):
        y = row_y[row_index]
        elements.append(svg_text(label_x, y + 49, label, size=25, weight="bold"))
        for line_index, note in enumerate(notes):
            elements.append(
                svg_text(
                    label_x,
                    y + 78 + 23 * line_index,
                    note,
                    size=15,
                    fill="#55585C",
                )
            )
        for slot, (parent, order, is_shuffled) in enumerate(
            zip(parents, descendants, shuffled), start=1
        ):
            elements.extend(
                draw_parent_block(
                    block_x[slot - 1],
                    y,
                    parent,
                    order,
                    shuffled=is_shuffled,
                )
            )
            data_rows.append(
                {
                    "row": label,
                    "slot": slot,
                    "parent": parent,
                    "descendant_order": "|".join(order),
                    "parent_order_changed": row_index == 1,
                    "descendant_order_changed": row_index == 2,
                    "contiguity_preserved": True,
                }
            )
        if row_index < 2:
            arrow_y = y + 135
            elements.append(
                f'<path d="M 800 {arrow_y - 16} L 800 {arrow_y + 12}" '
                'stroke="#777A7E" stroke-width="3" fill="none"/>'
            )
            elements.append(
                f'<path d="M 791 {arrow_y + 3} L 800 {arrow_y + 14} '
                f'L 809 {arrow_y + 3}" fill="none" stroke="#777A7E" stroke-width="3"/>'
            )
    elements.append("</svg>")
    svg_path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    pd.DataFrame(data_rows).to_csv(data_path, index=False)
    return svg_path, data_path


def write_statistics(
    output_dir: Path,
    heldout: CorrelationSummary,
    matched: CorrelationSummary,
    medians: dict[str, tuple[float, float]],
) -> Path:
    rows = [
        {
            "figure": "slide2_heldout_cases",
            "baseline_schedule": PRIMARY_BASELINE,
            "proposed_schedule": PROPOSED_ORDERING,
            "pearson_r": heldout.pearson_r,
            "pearson_p_value": heldout.p_value,
            "direction_agreement": heldout.direction_count,
            "direction_total": heldout.direction_total,
            "schedule": "",
            "median_bch_ratio": "",
            "median_error_ratio": "",
        },
        {
            "figure": "slide2_matched_pairs",
            "baseline_schedule": PRIMARY_BASELINE,
            "proposed_schedule": PROPOSED_ORDERING,
            "pearson_r": matched.pearson_r,
            "pearson_p_value": matched.p_value,
            "direction_agreement": matched.direction_count,
            "direction_total": matched.direction_total,
            "schedule": "",
            "median_bch_ratio": "",
            "median_error_ratio": "",
        },
    ]
    for schedule in PARENT_ABLATION_SCHEDULES:
        bch_median, error_median = medians[schedule]
        rows.append(
            {
                "figure": "slide3_parent_block_ablation",
                "baseline_schedule": "",
                "proposed_schedule": "",
                "pearson_r": "",
                "pearson_p_value": "",
                "direction_agreement": "",
                "direction_total": "",
                "schedule": schedule,
                "median_bch_ratio": bch_median,
                "median_error_ratio": error_median,
            }
        )
    path = output_dir / "recomputed_statistics.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def generate_all(
    *,
    case_summary_path: Path,
    matched_pair_path: Path,
    parent_ablation_path: Path,
    output_dir: Path,
) -> PresentationResults:
    case_summary = pd.read_csv(case_summary_path)
    matched_pair_summary = pd.read_csv(matched_pair_path)
    parent_ablation = pd.read_csv(parent_ablation_path)

    heldout_data, heldout_stats = prepare_heldout_cases(case_summary)
    matched_data, matched_stats = prepare_matched_pairs(
        matched_pair_summary,
        heldout_data,
    )
    parent_data, parent_medians = prepare_parent_ablation(parent_ablation)
    validate_expected_statistics(heldout_stats, matched_stats, parent_medians)

    output_dir.mkdir(parents=True, exist_ok=True)
    configure_plot_style()
    outputs: list[Path] = []
    outputs.extend(plot_heldout_cases(heldout_data, heldout_stats, output_dir))
    heldout_csv = output_dir / "slide2_heldout_cases_data.csv"
    heldout_data.to_csv(heldout_csv, index=False)
    outputs.append(heldout_csv)

    outputs.extend(plot_matched_pairs(matched_data, matched_stats, output_dir))
    matched_csv = output_dir / "slide2_matched_pairs_data.csv"
    matched_data.to_csv(matched_csv, index=False)
    outputs.append(matched_csv)

    outputs.extend(plot_parent_ablation(parent_data, parent_medians, output_dir))
    parent_csv = output_dir / "slide3_parent_block_ablation_data.csv"
    parent_data.to_csv(parent_csv, index=False)
    outputs.append(parent_csv)

    schematic_svg, schematic_csv = write_parent_order_schematic(output_dir)
    outputs.extend((schematic_svg, schematic_csv))
    outputs.append(
        write_statistics(output_dir, heldout_stats, matched_stats, parent_medians)
    )

    return PresentationResults(
        heldout=heldout_stats,
        matched_pairs=matched_stats,
        parent_medians=parent_medians,
        output_files=tuple(outputs),
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-summary", type=Path, default=DEFAULT_CASE_SUMMARY)
    parser.add_argument("--matched-pairs", type=Path, default=DEFAULT_MATCHED_PAIRS)
    parser.add_argument("--parent-ablation", type=Path, default=DEFAULT_PARENT_ABLATION)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    root = repository_root()
    results = generate_all(
        case_summary_path=resolve_from_root(root, args.case_summary),
        matched_pair_path=resolve_from_root(root, args.matched_pairs),
        parent_ablation_path=resolve_from_root(root, args.parent_ablation),
        output_dir=resolve_from_root(root, args.output_dir),
    )
    print(
        "Held-out cases: "
        f"r={results.heldout.pearson_r:.9f}, "
        f"p={results.heldout.p_value:.9g}, "
        f"direction={results.heldout.direction_count}/"
        f"{results.heldout.direction_total}"
    )
    print(
        "Matched pairs: "
        f"r={results.matched_pairs.pearson_r:.9f}, "
        f"p={results.matched_pairs.p_value:.9g}, "
        f"direction={results.matched_pairs.direction_count}/"
        f"{results.matched_pairs.direction_total}"
    )
    for schedule in PARENT_ABLATION_SCHEDULES:
        bch_median, error_median = results.parent_medians[schedule]
        print(
            f"{schedule}: median R_BCH={bch_median:.9f}x, "
            f"median error={error_median:.9f}x"
        )
    print("Generated files:")
    for path in results.output_files:
        print(f"  {path}")


if __name__ == "__main__":
    main()
