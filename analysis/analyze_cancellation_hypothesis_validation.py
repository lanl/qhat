#!/usr/bin/env python3
"""Analyze within-Hamiltonian tests of the BCH-cancellation hypothesis.

The primary test compares every alternative ordering with the signed fermionic
reference for the same case, time, and Trotter step count:

    delta_cancel = log10(R_BCH,alternative / R_BCH,signed)
    delta_error  = log10(error_alternative / error_signed)

A positive within-case association means that weakening destructive BCH
cancellation accompanies a larger Trotter error. Random schedule replicates
are first reduced to one median per case/schedule/condition, and repeated
Trotter conditions are then reduced to one median per case/schedule for the
primary test. Neither seeds nor step-count repeats can inflate its sample size.
"""

# ruff: noqa: E402  # Configure a writable Matplotlib cache before import.

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
from scipy import stats


REFERENCE = "fermionic_signed_reference"
DEFAULT_INPUT = Path(
    "analysis/cancellation_hypothesis_validation/pilot_ablation_results.csv"
)
DEFAULT_MANIFEST = Path(
    "analysis/cancellation_hypothesis_validation/pilot_panel_manifest.csv"
)
DEFAULT_OUTDIR = Path(
    "analysis/cancellation_hypothesis_validation/pilot_analysis"
)

SCHEDULE_LABELS = {
    "fermionic_magnitude_reference": "Fermionic magnitude",
    "jw_signed_baseline": "JW signed",
    "jw_magnitude_baseline": "JW magnitude",
    "signed_parent_descendants_round_robin": "Round robin",
    "signed_parent_blocks_randomized": "Random parent blocks",
    "signed_parent_within_randomized": "Within-parent shuffle",
}

REQUIRED_RESULT_COLUMNS = {
    "status",
    "case_id",
    "molecule",
    "schedule",
    "sample_index",
    "trotter_steps",
    "evolution_time",
    "bch_cancellation_ratio",
    "one_minus_overlap",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, dest="inputs")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--error-floor", type=float, default=1.0e-13)
    parser.add_argument("--bootstrap-samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260822)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def read_results(paths: Sequence[Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        frame = pd.read_csv(path)
        missing = sorted(REQUIRED_RESULT_COLUMNS.difference(frame.columns))
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        frame["source_csv"] = str(path)
        frames.append(frame)
    if not frames:
        raise ValueError("At least one validation result CSV is required.")
    combined = pd.concat(frames, ignore_index=True)
    combined = combined[combined["status"] == "success"].copy()
    dedupe = [
        "case_id",
        "schedule",
        "sample_index",
        "random_seed",
        "trotter_steps",
        "evolution_time",
    ]
    dedupe = [column for column in dedupe if column in combined.columns]
    return combined.drop_duplicates(dedupe, keep="last").reset_index(drop=True)


def build_within_case_deltas(
    results: pd.DataFrame,
    manifest: pd.DataFrame,
    error_floor: float,
) -> pd.DataFrame:
    keys = ["case_id", "trotter_steps", "evolution_time"]
    references = results[results["schedule"] == REFERENCE].copy()
    duplicate_references = references.duplicated(keys, keep=False)
    if duplicate_references.any():
        duplicates = references.loc[duplicate_references, keys]
        raise ValueError(
            "Expected one signed reference per condition; duplicates found:\n"
            f"{duplicates.to_string(index=False)}"
        )
    if references.empty:
        raise ValueError("No fermionic signed reference rows found.")

    reference_columns = keys + [
        "bch_cancellation_ratio",
        "one_minus_overlap",
    ]
    references = references[reference_columns].rename(
        columns={
            "bch_cancellation_ratio": "reference_bch_cancellation_ratio",
            "one_minus_overlap": "reference_one_minus_overlap",
        }
    )
    alternatives = results[results["schedule"] != REFERENCE].copy()
    deltas = alternatives.merge(references, on=keys, how="inner", validate="many_to_one")
    manifest_columns = [
        "case_id",
        "matched_pair",
        "expected_outcome",
        "fermionic_advantage_factor",
    ]
    if "n_qubits" not in deltas.columns:
        manifest_columns.append("n_qubits")
    deltas = deltas.merge(
        manifest[manifest_columns],
        on="case_id",
        how="left",
        validate="many_to_one",
    )

    positive = (
        np.isfinite(deltas["bch_cancellation_ratio"])
        & np.isfinite(deltas["reference_bch_cancellation_ratio"])
        & (deltas["bch_cancellation_ratio"] > 0.0)
        & (deltas["reference_bch_cancellation_ratio"] > 0.0)
        & np.isfinite(deltas["one_minus_overlap"])
        & np.isfinite(deltas["reference_one_minus_overlap"])
        & (deltas["one_minus_overlap"] > error_floor)
        & (deltas["reference_one_minus_overlap"] > error_floor)
    )
    deltas["delta_valid"] = positive
    deltas["cancellation_ratio_to_signed"] = np.where(
        positive,
        deltas["bch_cancellation_ratio"]
        / deltas["reference_bch_cancellation_ratio"],
        np.nan,
    )
    deltas["error_ratio_to_signed"] = np.where(
        positive,
        deltas["one_minus_overlap"] / deltas["reference_one_minus_overlap"],
        np.nan,
    )
    deltas["delta_log10_cancellation"] = np.where(
        positive,
        np.log10(deltas["cancellation_ratio_to_signed"]),
        np.nan,
    )
    deltas["delta_log10_error"] = np.where(
        positive,
        np.log10(deltas["error_ratio_to_signed"]),
        np.nan,
    )
    deltas["signed_has_stronger_cancellation"] = (
        deltas["delta_log10_cancellation"] > 0.0
    )
    deltas["signed_has_lower_error"] = deltas["delta_log10_error"] > 0.0
    deltas["direction_concordant"] = (
        deltas["delta_log10_cancellation"]
        * deltas["delta_log10_error"]
        > 0.0
    )
    return deltas.sort_values(
        ["case_id", "evolution_time", "trotter_steps", "schedule", "sample_index"]
    ).reset_index(drop=True)


def aggregate_random_replicates(deltas: pd.DataFrame) -> pd.DataFrame:
    valid = deltas[deltas["delta_valid"]].copy()
    group_keys = [
        "case_id",
        "molecule",
        "matched_pair",
        "expected_outcome",
        "fermionic_advantage_factor",
        "n_qubits",
        "schedule",
        "trotter_steps",
        "evolution_time",
    ]
    log_columns = [
        "delta_log10_cancellation",
        "delta_log10_error",
    ]
    aggregated = valid.groupby(group_keys, as_index=False)[log_columns].median()
    aggregated["cancellation_ratio_to_signed"] = 10.0 ** aggregated[
        "delta_log10_cancellation"
    ]
    aggregated["error_ratio_to_signed"] = 10.0 ** aggregated[
        "delta_log10_error"
    ]
    counts = valid.groupby(group_keys, as_index=False).size().rename(
        columns={"size": "replicates"}
    )
    aggregated = aggregated.merge(counts, on=group_keys, validate="one_to_one")
    aggregated["signed_has_stronger_cancellation"] = (
        aggregated["delta_log10_cancellation"] > 0.0
    )
    aggregated["signed_has_lower_error"] = aggregated["delta_log10_error"] > 0.0
    aggregated["direction_concordant"] = (
        aggregated["delta_log10_cancellation"]
        * aggregated["delta_log10_error"]
        > 0.0
    )
    return aggregated.sort_values(
        ["case_id", "evolution_time", "trotter_steps", "schedule"]
    ).reset_index(drop=True)


def collapse_conditions_for_inference(aggregated: pd.DataFrame) -> pd.DataFrame:
    """Return one log-space median per case and alternative schedule."""
    columns = [
        "delta_log10_cancellation",
        "delta_log10_error",
    ]
    collapsed = aggregated.groupby(
        ["case_id", "schedule"],
        as_index=False,
    )[columns].median()
    collapsed["cancellation_ratio_to_signed"] = 10.0 ** collapsed[
        "delta_log10_cancellation"
    ]
    collapsed["error_ratio_to_signed"] = 10.0 ** collapsed[
        "delta_log10_error"
    ]
    return collapsed


def collapse_random_conditions_for_inference(deltas: pd.DataFrame) -> pd.DataFrame:
    """Return one log-space median per case and randomized schedule seed."""
    random_blocks = deltas[
        deltas["delta_valid"]
        & (deltas["schedule"] == "signed_parent_blocks_randomized")
    ].copy()
    seed_columns = ["case_id", "sample_index"]
    if "random_seed" in random_blocks.columns:
        seed_columns.append("random_seed")
    return random_blocks.groupby(seed_columns, as_index=False)[
        ["delta_log10_cancellation", "delta_log10_error"]
    ].median()


def residualize_by_groups(
    frame: pd.DataFrame,
    column: str,
    groups: Sequence[str],
) -> np.ndarray:
    values = frame[column].to_numpy(dtype=float)
    residual = values.copy()
    for _, indices in frame.groupby(list(groups), sort=False).groups.items():
        positions = frame.index.get_indexer(indices)
        residual[positions] -= float(np.mean(values[positions]))
    return residual


def correlation_statistics(x: np.ndarray, y: np.ndarray) -> dict[str, float | int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    result: dict[str, float | int] = {
        "observations": int(len(x)),
        "pearson_r": math.nan,
        "pearson_p_value": math.nan,
        "spearman_rho": math.nan,
        "spearman_p_value": math.nan,
        "ols_slope": math.nan,
        "ols_intercept": math.nan,
    }
    if len(x) < 3 or np.std(x) == 0.0 or np.std(y) == 0.0:
        return result
    pearson = stats.pearsonr(x, y)
    spearman = stats.spearmanr(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    result.update(
        {
            "pearson_r": float(pearson.statistic),
            "pearson_p_value": float(pearson.pvalue),
            "spearman_rho": float(spearman.statistic),
            "spearman_p_value": float(spearman.pvalue),
            "ols_slope": float(slope),
            "ols_intercept": float(intercept),
        }
    )
    return result


def fixed_effect_correlation(frame: pd.DataFrame) -> dict[str, float | int]:
    groups = [
        column
        for column in ["case_id", "trotter_steps", "evolution_time"]
        if column in frame.columns
    ]
    x = residualize_by_groups(frame, "delta_log10_cancellation", groups)
    y = residualize_by_groups(frame, "delta_log10_error", groups)
    return correlation_statistics(x, y)


def bootstrap_case_fixed_effect_r(
    frame: pd.DataFrame,
    samples: int,
    seed: int,
) -> tuple[float, float, int]:
    case_ids = frame["case_id"].unique()
    if len(case_ids) < 2 or samples <= 0:
        return math.nan, math.nan, 0
    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    for _ in range(samples):
        selected = rng.choice(case_ids, size=len(case_ids), replace=True)
        pieces: list[pd.DataFrame] = []
        for draw_index, case_id in enumerate(selected):
            piece = frame[frame["case_id"] == case_id].copy()
            piece["case_id"] = f"{case_id}#bootstrap-{draw_index}"
            pieces.append(piece)
        sample = pd.concat(pieces, ignore_index=True)
        estimate = float(fixed_effect_correlation(sample)["pearson_r"])
        if math.isfinite(estimate):
            estimates.append(estimate)
    if not estimates:
        return math.nan, math.nan, 0
    lower, upper = np.quantile(estimates, [0.025, 0.975])
    return float(lower), float(upper), len(estimates)


def build_schedule_summary(aggregated: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for schedule, group in aggregated.groupby("schedule", sort=True):
        rows.append(
            {
                "schedule": schedule,
                "schedule_label": SCHEDULE_LABELS.get(schedule, schedule),
                "case_conditions": len(group),
                "cases": int(group["case_id"].nunique()),
                "median_cancellation_ratio_to_signed": float(
                    group["cancellation_ratio_to_signed"].median()
                ),
                "median_error_ratio_to_signed": float(
                    group["error_ratio_to_signed"].median()
                ),
                "signed_stronger_cancellation_fraction": float(
                    group["signed_has_stronger_cancellation"].mean()
                ),
                "signed_lower_error_fraction": float(
                    group["signed_has_lower_error"].mean()
                ),
                "direction_concordance_fraction": float(
                    group["direction_concordant"].mean()
                ),
            }
        )
    return pd.DataFrame(rows)


def build_case_summary(
    aggregated: pd.DataFrame,
    results: pd.DataFrame,
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    reference = results[results["schedule"] == REFERENCE].copy()
    reference = reference.groupby("case_id", as_index=False).agg(
        signed_bch_cancellation_ratio=("bch_cancellation_ratio", "median"),
        signed_one_minus_overlap=("one_minus_overlap", "median"),
    )
    deterministic_errors = results[
        results["schedule"].isin(
            {REFERENCE, "jw_signed_baseline", "jw_magnitude_baseline"}
        )
    ].pivot_table(
        index="case_id",
        columns="schedule",
        values="one_minus_overlap",
        aggfunc="median",
    )
    deterministic_cancellation = results[
        results["schedule"].isin(
            {REFERENCE, "jw_signed_baseline", "jw_magnitude_baseline"}
        )
    ].pivot_table(
        index="case_id",
        columns="schedule",
        values="bch_cancellation_ratio",
        aggfunc="median",
    )
    required_schedules = {
        REFERENCE,
        "jw_signed_baseline",
        "jw_magnitude_baseline",
    }
    missing_schedules = sorted(required_schedules.difference(deterministic_errors.columns))
    if missing_schedules:
        raise ValueError(
            "Fresh performance classification requires schedules: "
            f"{missing_schedules}"
        )
    fresh = pd.DataFrame(
        {
            "case_id": deterministic_errors.index,
            "fresh_signed_one_minus_overlap": deterministic_errors[REFERENCE],
            "fresh_jw_signed_one_minus_overlap": deterministic_errors[
                "jw_signed_baseline"
            ],
            "fresh_jw_magnitude_one_minus_overlap": deterministic_errors[
                "jw_magnitude_baseline"
            ],
        }
    ).reset_index(drop=True)
    fresh["fresh_best_jw_one_minus_overlap"] = fresh[
        [
            "fresh_jw_signed_one_minus_overlap",
            "fresh_jw_magnitude_one_minus_overlap",
        ]
    ].min(axis=1)
    fresh["fresh_best_jw_schedule"] = deterministic_errors[
        ["jw_signed_baseline", "jw_magnitude_baseline"]
    ].idxmin(axis=1).to_numpy()
    fresh["fresh_best_jw_bch_cancellation_ratio"] = [
        float(deterministic_cancellation.loc[case_id, schedule])
        for case_id, schedule in zip(
            fresh["case_id"],
            fresh["fresh_best_jw_schedule"],
            strict=True,
        )
    ]
    fresh["fresh_best_jw_to_signed_bch_cancellation_ratio"] = (
        fresh["fresh_best_jw_bch_cancellation_ratio"]
        / fresh["case_id"].map(deterministic_cancellation[REFERENCE])
    )
    fresh["fresh_jw_to_signed_advantage"] = (
        fresh["fresh_best_jw_one_minus_overlap"]
        / fresh["fresh_signed_one_minus_overlap"]
    )
    fresh["fresh_observed_outcome"] = np.where(
        fresh["fresh_jw_to_signed_advantage"] > 1.0,
        "favorable",
        "negative_control",
    )
    rows: list[dict[str, Any]] = []
    for case_id, group in aggregated.groupby("case_id", sort=True):
        stats_row = correlation_statistics(
            group["delta_log10_cancellation"].to_numpy(float),
            group["delta_log10_error"].to_numpy(float),
        )
        rows.append(
            {
                "case_id": case_id,
                "alternative_conditions": len(group),
                "alternatives_with_weaker_cancellation": int(
                    group["signed_has_stronger_cancellation"].sum()
                ),
                "alternatives_with_higher_error": int(
                    group["signed_has_lower_error"].sum()
                ),
                "direction_concordant_alternatives": int(
                    group["direction_concordant"].sum()
                ),
                "within_case_pearson_r": stats_row["pearson_r"],
                "within_case_spearman_rho": stats_row["spearman_rho"],
            }
        )
    summary = pd.DataFrame(rows)
    summary = manifest.merge(summary, on="case_id", how="left", validate="one_to_one")
    summary = summary.merge(reference, on="case_id", how="left", validate="one_to_one")
    summary = summary.merge(fresh, on="case_id", how="left", validate="one_to_one")
    summary["performance_label_reproduced"] = (
        summary["expected_outcome"] == summary["fresh_observed_outcome"]
    )
    if "historical_fermionic_one_minus_overlap" in summary:
        summary["fresh_to_historical_signed_error_ratio"] = (
            summary["fresh_signed_one_minus_overlap"]
            / summary["historical_fermionic_one_minus_overlap"]
        )
    return summary


def build_validation_statistics(
    deltas: pd.DataFrame,
    aggregated: pd.DataFrame,
    case_summary: pd.DataFrame,
    bootstrap_samples: int,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    primary_frame = collapse_conditions_for_inference(aggregated)
    fixed = fixed_effect_correlation(primary_frame)
    lower, upper, valid_bootstrap = bootstrap_case_fixed_effect_r(
        primary_frame,
        bootstrap_samples,
        seed,
    )
    rows.append(
        {
            "test": "within_case_delta_correlation",
            "description": (
                "Case-centered delta log10 BCH cancellation ratio versus delta "
                "log10 one-minus-overlap after collapsing Trotter conditions"
            ),
            **fixed,
            "bootstrap_95ci_lower": lower,
            "bootstrap_95ci_upper": upper,
            "valid_bootstrap_samples": valid_bootstrap,
        }
    )

    deterministic_schedules = {
        "fermionic_magnitude_reference",
        "jw_signed_baseline",
        "jw_magnitude_baseline",
        "signed_parent_descendants_round_robin",
    }
    deterministic = primary_frame[
        primary_frame["schedule"].isin(deterministic_schedules)
    ].copy()
    deterministic_fixed = fixed_effect_correlation(deterministic)
    deterministic_lower, deterministic_upper, deterministic_bootstrap = (
        bootstrap_case_fixed_effect_r(
            deterministic,
            bootstrap_samples,
            seed + 1,
        )
    )
    rows.append(
        {
            "test": "within_case_deterministic_controls",
            "description": (
                "Case-centered delta correlation using only the four deterministic "
                "alternative schedules after collapsing Trotter conditions"
            ),
            **deterministic_fixed,
            "bootstrap_95ci_lower": deterministic_lower,
            "bootstrap_95ci_upper": deterministic_upper,
            "valid_bootstrap_samples": deterministic_bootstrap,
        }
    )

    random_blocks = collapse_random_conditions_for_inference(deltas)
    random_fixed = fixed_effect_correlation(random_blocks)
    random_lower, random_upper, random_bootstrap = bootstrap_case_fixed_effect_r(
        random_blocks,
        bootstrap_samples,
        seed + 2,
    )
    rows.append(
        {
            "test": "within_case_random_parent_blocks",
            "description": (
                "Case/condition-centered delta correlation across raw randomized "
                "parent-block schedules"
            ),
            **random_fixed,
            "bootstrap_95ci_lower": random_lower,
            "bootstrap_95ci_upper": random_upper,
            "valid_bootstrap_samples": random_bootstrap,
        }
    )

    for (steps, evolution_time), condition in aggregated.groupby(
        ["trotter_steps", "evolution_time"],
        sort=True,
    ):
        condition_fixed = fixed_effect_correlation(condition)
        rows.append(
            {
                "test": f"within_case_condition_r{steps:g}_t{evolution_time:g}",
                "description": (
                    "Case-centered delta correlation at "
                    f"r={steps:g}, T={evolution_time:g}"
                ),
                **condition_fixed,
                "bootstrap_95ci_lower": math.nan,
                "bootstrap_95ci_upper": math.nan,
                "valid_bootstrap_samples": 0,
            }
        )

    direct = case_summary[
        np.isfinite(
            case_summary["fresh_best_jw_to_signed_bch_cancellation_ratio"]
        )
        & (case_summary["fresh_best_jw_to_signed_bch_cancellation_ratio"] > 0.0)
        & np.isfinite(case_summary["fresh_jw_to_signed_advantage"])
        & (case_summary["fresh_jw_to_signed_advantage"] > 0.0)
    ]
    direct_stats = correlation_statistics(
        np.log10(
            direct["fresh_best_jw_to_signed_bch_cancellation_ratio"].to_numpy(
                float
            )
        ),
        np.log10(direct["fresh_jw_to_signed_advantage"].to_numpy(float)),
    )
    rows.append(
        {
            "test": "relative_cancellation_vs_fresh_advantage",
            "description": (
                "Best-JW/signed BCH cancellation ratio versus freshly recomputed "
                "best-JW/signed error ratio"
            ),
            **direct_stats,
            "bootstrap_95ci_lower": math.nan,
            "bootstrap_95ci_upper": math.nan,
            "valid_bootstrap_samples": 0,
        }
    )
    absolute = case_summary[
        np.isfinite(case_summary["signed_bch_cancellation_ratio"])
        & (case_summary["signed_bch_cancellation_ratio"] > 0.0)
        & np.isfinite(case_summary["fresh_jw_to_signed_advantage"])
        & (case_summary["fresh_jw_to_signed_advantage"] > 0.0)
    ]
    absolute_stats = correlation_statistics(
        -np.log10(absolute["signed_bch_cancellation_ratio"].to_numpy(float)),
        np.log10(absolute["fresh_jw_to_signed_advantage"].to_numpy(float)),
    )
    rows.append(
        {
            "test": "absolute_signed_cancellation_vs_fresh_advantage",
            "description": (
                "Absolute signed cancellation strength versus freshly recomputed "
                "best-JW/signed error ratio"
            ),
            **absolute_stats,
            "bootstrap_95ci_lower": math.nan,
            "bootstrap_95ci_upper": math.nan,
            "valid_bootstrap_samples": 0,
        }
    )
    return pd.DataFrame(rows)


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.2,
        }
    )


def make_plot(
    aggregated: pd.DataFrame,
    case_summary: pd.DataFrame,
    outdir: Path,
    dpi: int,
) -> None:
    configure_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4))
    for schedule, group in aggregated.groupby("schedule", sort=True):
        axes[0].scatter(
            group["delta_log10_cancellation"],
            group["delta_log10_error"],
            s=38,
            alpha=0.8,
            label=SCHEDULE_LABELS.get(schedule, schedule),
        )
    limits = axes[0].axis()
    low = min(limits[0], limits[2])
    high = max(limits[1], limits[3])
    axes[0].plot([low, high], [low, high], color="black", linestyle=":", linewidth=0.8)
    axes[0].axhline(0.0, color="black", linestyle="--", linewidth=0.7)
    axes[0].axvline(0.0, color="black", linestyle="--", linewidth=0.7)
    axes[0].set_xlabel(r"$\Delta\log_{10} R_{BCH}$ (alternative / signed)")
    axes[0].set_ylabel(r"$\Delta\log_{10}\epsilon$ (alternative / signed)")
    axes[0].set_title("Within-case ordering ablations")
    axes[0].legend(frameon=False, fontsize=7)

    valid = case_summary[
        np.isfinite(
            case_summary["fresh_best_jw_to_signed_bch_cancellation_ratio"]
        )
        & (case_summary["fresh_best_jw_to_signed_bch_cancellation_ratio"] > 0.0)
        & np.isfinite(case_summary["fresh_jw_to_signed_advantage"])
        & (case_summary["fresh_jw_to_signed_advantage"] > 0.0)
    ]
    colors = valid["expected_outcome"].map(
        {"favorable": "#2a9d8f", "negative_control": "#e76f51"}
    )
    x = np.log10(valid["fresh_best_jw_to_signed_bch_cancellation_ratio"])
    y = np.log10(valid["fresh_jw_to_signed_advantage"])
    axes[1].scatter(x, y, c=colors, s=48)
    for row, x_value, y_value in zip(valid.itertuples(), x, y, strict=True):
        axes[1].annotate(
            row.molecule,
            (x_value, y_value),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )
    axes[1].axhline(0.0, color="black", linestyle="--", linewidth=0.7)
    axes[1].axvline(0.0, color="black", linestyle="--", linewidth=0.7)
    axes[1].set_xlabel(r"$\log_{10}(R_{best-JW}/R_{signed})$")
    axes[1].set_ylabel(r"$\log_{10}(\epsilon_{best-JW}/\epsilon_{signed})$")
    axes[1].set_title("BCH-held-out matched cases")
    fig.tight_layout()
    fig.savefig(outdir / "cancellation_validation.png", dpi=dpi)
    fig.savefig(outdir / "cancellation_validation.pdf")
    plt.close(fig)


def format_number(value: object, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{number:.{digits}g}" if math.isfinite(number) else "NA"


def write_report(
    results: pd.DataFrame,
    manifest: pd.DataFrame,
    deltas: pd.DataFrame,
    aggregated: pd.DataFrame,
    schedule_summary: pd.DataFrame,
    case_summary: pd.DataFrame,
    validation_statistics: pd.DataFrame,
    outdir: Path,
) -> None:
    primary = validation_statistics[
        validation_statistics["test"] == "within_case_delta_correlation"
    ].iloc[0]
    deterministic = validation_statistics[
        validation_statistics["test"] == "within_case_deterministic_controls"
    ].iloc[0]
    random_blocks = validation_statistics[
        validation_statistics["test"] == "within_case_random_parent_blocks"
    ].iloc[0]
    condition_rows = validation_statistics[
        validation_statistics["test"].str.startswith("within_case_condition_")
    ]
    direct = validation_statistics[
        validation_statistics["test"] == "relative_cancellation_vs_fresh_advantage"
    ].iloc[0]
    absolute = validation_statistics[
        validation_statistics["test"]
        == "absolute_signed_cancellation_vs_fresh_advantage"
    ].iloc[0]
    concordance = float(aggregated["direction_concordant"].mean())
    signed_better_both = float(
        (
            aggregated["signed_has_stronger_cancellation"]
            & aggregated["signed_has_lower_error"]
        ).mean()
    )
    primary_supports = (
        float(primary["pearson_r"]) > 0.0
        and float(primary["pearson_p_value"]) < 0.05
        and concordance >= 0.7
    )
    direct_supports = (
        float(direct["pearson_r"]) > 0.0
        and float(direct["pearson_p_value"]) < 0.05
    )
    if primary_supports and direct_supports:
        verdict = "supports both the within-Hamiltonian and cross-case claims"
    elif primary_supports:
        verdict = (
            "supports the within-Hamiltonian cancellation link, but does not "
            "yet establish the cross-case claim"
        )
    else:
        verdict = "does not yet establish the proposed cancellation mechanism"
    reproduced = int(case_summary["performance_label_reproduced"].fillna(False).sum())

    lines = [
        "# BCH cancellation hypothesis validation",
        "",
        "## Result",
        "",
        f"This validation panel **{verdict}**. The "
        "primary analysis changes only the Pauli schedule within each fixed "
        "Hamiltonian and compares every alternative with the signed fermionic "
        "reference.",
        "",
        "## Coverage",
        "",
        f"- Manifest cases: {len(manifest)}",
        f"- Successful benchmark rows: {len(results)}",
        f"- Valid raw alternative/reference deltas: {int(deltas['delta_valid'].sum())}",
        f"- Case/schedule/condition medians: {len(aggregated)}",
        f"- Primary case/schedule inference units: {int(primary['observations'])}",
        f"- Historical performance labels reproduced from current tensors: "
        f"{reproduced}/{len(case_summary)}",
        f"- Trotter step counts: {sorted(results['trotter_steps'].unique().tolist())}",
        f"- Evolution times: {sorted(results['evolution_time'].unique().tolist())}",
        "",
        "## Primary within-case test",
        "",
        "The variables are `delta log10 cancellation = log10(R_alt/R_signed)` "
        "and `delta log10 error = log10(epsilon_alt/epsilon_signed)`. A positive "
        "association is the predicted direction.",
        "",
        f"- Case-centered Pearson r: {format_number(primary['pearson_r'])}",
        f"- p-value: {format_number(primary['pearson_p_value'])}",
        "- Case-block bootstrap 95% interval: "
        f"[{format_number(primary['bootstrap_95ci_lower'])}, "
        f"{format_number(primary['bootstrap_95ci_upper'])}]",
        f"- Direction-concordant deltas: {concordance:.1%}",
        f"- Signed stronger-cancellation and lower-error deltas: {signed_better_both:.1%}",
        "- Deterministic-control sensitivity: "
        f"r={format_number(deterministic['pearson_r'])}, "
        f"p={format_number(deterministic['pearson_p_value'])}, "
        f"n={int(deterministic['observations'])}",
        "- Random-parent-block sensitivity (raw seeds, centered within case): "
        f"r={format_number(random_blocks['pearson_r'])}, "
        f"p={format_number(random_blocks['pearson_p_value'])}, "
        f"n={int(random_blocks['observations'])}",
        "- Step-specific Pearson r values: "
        + ", ".join(
            f"{row.test.removeprefix('within_case_condition_')}="
            f"{format_number(row.pearson_r)}"
            for row in condition_rows.itertuples()
        ),
        "",
        "## Held-out cross-case check",
        "",
        "The panel cases were not among the eight BCH cases used to formulate "
        "the hypothesis. Favorable/control labels were used only to construct "
        "the matched panel; the reported advantage below is freshly recomputed "
        "from this run's JW-signed, JW-magnitude, and signed-reference rows.",
        "",
        "- Best-JW/signed cancellation ratio vs best-JW/signed error ratio: "
        f"r={format_number(direct['pearson_r'])}, "
        f"p={format_number(direct['pearson_p_value'])}, n={int(direct['observations'])}",
        f"- Absolute signed cancellation strength (supplemental): "
        f"r={format_number(absolute['pearson_r'])}, "
        f"p={format_number(absolute['pearson_p_value'])}, "
        f"n={int(absolute['observations'])}",
        "",
        "## Schedule medians relative to signed fermionic",
        "",
        "| Schedule | BCH cancellation ratio | Error ratio | Concordance |",
        "|---|---:|---:|---:|",
    ]
    for row in schedule_summary.itertuples():
        lines.append(
            f"| {row.schedule_label} | "
            f"{row.median_cancellation_ratio_to_signed:.4g} | "
            f"{row.median_error_ratio_to_signed:.4g} | "
            f"{row.direction_concordance_fraction:.1%} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation limits",
            "",
            "- This is a BCH-held-out test, not a fully prospective performance "
            "test: favorable/control labels came from the existing Trotter sweep.",
            "- Random seeds are summarized by their median before the primary "
            "correlation; repeated step counts are also collapsed before primary "
            "inference, preventing either from inflating sample size.",
            "- `R_BCH` is evaluated on the HF initial state. A later campaign "
            "should repeat the analysis with propagated local-error vectors over time.",
            "- A pilot panel can establish workflow correctness and effect direction, "
            "but the full 20-case panel is needed for a report-level mechanism claim.",
            "",
            "## Files",
            "",
            "- `within_case_deltas.csv`: every raw alternative/reference comparison.",
            "- `primary_aggregated_deltas.csv`: one median per case/schedule/condition.",
            "- `schedule_summary.csv`: schedule-level ratios and concordance.",
            "- `case_summary.csv`: case-level outcome and cancellation results.",
            "- `validation_statistics.csv`: primary and held-out correlations.",
            "- `cancellation_validation.{png,pdf}`: diagnostic figure.",
            "",
        ]
    )
    (outdir / "validation_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_arguments()
    inputs = args.inputs or [DEFAULT_INPUT]
    if args.error_floor <= 0.0:
        raise ValueError("--error-floor must be positive.")
    if args.bootstrap_samples < 0:
        raise ValueError("--bootstrap-samples cannot be negative.")
    results = read_results(inputs)
    manifest = pd.read_csv(args.manifest)
    manifest_case_ids = set(manifest["case_id"])
    extra_cases = sorted(set(results["case_id"]).difference(manifest_case_ids))
    if extra_cases:
        print(f"WARNING: excluding result rows outside the manifest: {extra_cases}")
    results = results[results["case_id"].isin(manifest_case_ids)].copy()
    missing_cases = sorted(set(manifest["case_id"]).difference(results["case_id"]))
    if missing_cases:
        print(f"WARNING: no successful rows yet for cases: {missing_cases}")
    deltas = build_within_case_deltas(results, manifest, args.error_floor)
    aggregated = aggregate_random_replicates(deltas)
    if aggregated.empty:
        raise ValueError("No valid alternative/reference deltas were available.")
    schedule_summary = build_schedule_summary(aggregated)
    case_summary = build_case_summary(aggregated, results, manifest)
    validation_statistics = build_validation_statistics(
        deltas,
        aggregated,
        case_summary,
        args.bootstrap_samples,
        args.seed,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    deltas.to_csv(args.outdir / "within_case_deltas.csv", index=False)
    aggregated.to_csv(args.outdir / "primary_aggregated_deltas.csv", index=False)
    schedule_summary.to_csv(args.outdir / "schedule_summary.csv", index=False)
    case_summary.to_csv(args.outdir / "case_summary.csv", index=False)
    validation_statistics.to_csv(
        args.outdir / "validation_statistics.csv",
        index=False,
    )
    make_plot(aggregated, case_summary, args.outdir, args.dpi)
    write_report(
        results,
        manifest,
        deltas,
        aggregated,
        schedule_summary,
        case_summary,
        validation_statistics,
        args.outdir,
    )
    print(f"Wrote cancellation validation analysis to {args.outdir}")


if __name__ == "__main__":
    main()
