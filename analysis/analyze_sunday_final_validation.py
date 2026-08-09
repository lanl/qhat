#!/usr/bin/env python3
"""Summarize the final fermionic-ordering validation campaign.

Run from the QHAT repository root after run_sunday_final_validation_hgbs5.sh.
The script is read-only with respect to benchmark CSVs.  It writes compact
summary CSVs, text conclusions, and separate PNG figures under
analysis/final_validation_summary/ by default.

No SciPy or seaborn is required: only numpy, pandas, and matplotlib.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REFERENCE = "fermionic_signed_reference"
FMAG = "fermionic_magnitude_reference"
JW = "jw_signed_baseline"
ROUND = "signed_parent_descendants_round_robin"
BLOCK_RANDOM = "signed_parent_blocks_randomized"
WITHIN_RANDOM = "signed_parent_within_randomized"

SCHEDULE_LABELS = {
    REFERENCE: "Fermionic signed",
    FMAG: "Fermionic magnitude",
    JW: "JW signed",
    ROUND: "Round robin",
    BLOCK_RANDOM: "Random parent blocks",
    WITHIN_RANDOM: "Within-parent shuffle",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ablation",
        type=Path,
        default=Path("analysis/fermionic_structure_ablation_hgbs5.csv"),
    )
    parser.add_argument(
        "--step-scaling",
        type=Path,
        default=Path("analysis/final_validation_step_scaling_hgbs5.csv"),
    )
    parser.add_argument(
        "--active-space",
        type=Path,
        default=Path("analysis/final_validation_active_space_hgbs5.csv"),
    )
    parser.add_argument(
        "--be2-time",
        type=Path,
        default=Path("analysis/final_validation_be2_fixed_dt_time_hgbs5.csv"),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("analysis/final_validation_summary"),
    )
    parser.add_argument(
        "--error-floor",
        type=float,
        default=1.0e-13,
        help="Ignore smaller values in log fits/plots. Default: 1e-13.",
    )
    return parser.parse_args()


def read_success(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "status" in df.columns:
        failed = df[df["status"] != "success"]
        if not failed.empty:
            print(f"WARNING: {path} contains {len(failed)} non-success rows.")
        df = df[df["status"] == "success"].copy()
    return df


def safe_log10(values: pd.Series, floor: float) -> np.ndarray:
    arr = values.to_numpy(dtype=float)
    arr = arr[np.isfinite(arr) & (arr > floor)]
    return np.log10(arr)


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return math.nan
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def fit_log_slope(steps: pd.Series, errors: pd.Series, floor: float) -> tuple[float, float, int]:
    x = steps.to_numpy(dtype=float)
    y = errors.to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > floor)
    if np.count_nonzero(mask) < 3:
        return math.nan, math.nan, int(np.count_nonzero(mask))
    lx = np.log10(x[mask])
    ly = np.log10(y[mask])
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = intercept + slope * lx
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else math.nan
    return float(slope), float(r2), int(np.count_nonzero(mask))


def summarize_ablation(df: pd.DataFrame, outdir: Path) -> list[str]:
    if df.empty:
        return ["Ablation CSV not found; ablation summary skipped."]

    rows: list[dict[str, float | str | int]] = []
    conclusions: list[str] = []

    same_owner_zero = (
        np.allclose(df["same_owner_pair_fraction"].fillna(0.0), 0.0)
        and np.allclose(df["same_owner_abs_weight_fraction"].fillna(0.0), 0.0)
        and np.allclose(df["bch_same_owner_norm"].fillna(0.0), 0.0)
    )
    conclusions.append(
        "Ablation: same-owner noncommuting contribution is zero in every row."
        if same_owner_zero
        else "Ablation: nonzero same-owner commutator contributions were detected."
    )

    for molecule, group in df.groupby("molecule", sort=True):
        ref_rows = group[group["schedule"] == REFERENCE]
        if ref_rows.empty:
            continue
        ref = ref_rows.iloc[0]
        block = group[group["schedule"] == BLOCK_RANDOM]
        within = group[group["schedule"] == WITHIN_RANDOM]

        row: dict[str, float | str | int] = {
            "molecule": molecule,
            "reference_one_minus_overlap": float(ref["one_minus_overlap"]),
            "reference_bch": float(ref["bch2_hf_state_norm"]),
        }
        if not block.empty:
            row.update(
                {
                    "random_block_samples": int(len(block)),
                    "random_block_error_ratio_min": float(block["one_minus_overlap_ratio_to_signed_reference"].min()),
                    "random_block_error_ratio_median": float(block["one_minus_overlap_ratio_to_signed_reference"].median()),
                    "random_block_error_ratio_max": float(block["one_minus_overlap_ratio_to_signed_reference"].max()),
                    "random_block_bch_ratio_median": float(block["bch_norm_ratio_to_signed_reference"].median()),
                    "random_blocks_beating_reference": int((block["one_minus_overlap"] < float(ref["one_minus_overlap"])).sum()),
                }
            )
            mask = (
                np.isfinite(block["bch2_hf_state_norm"].to_numpy(dtype=float))
                & np.isfinite(block["one_minus_overlap"].to_numpy(dtype=float))
                & (block["bch2_hf_state_norm"].to_numpy(dtype=float) > 0.0)
                & (block["one_minus_overlap"].to_numpy(dtype=float) > 0.0)
            )
            if np.count_nonzero(mask) >= 3:
                x = np.log10(block.loc[mask, "bch2_hf_state_norm"].to_numpy(dtype=float))
                y = np.log10(block.loc[mask, "one_minus_overlap"].to_numpy(dtype=float))
                row["random_block_log_bch_error_pearson"] = correlation(x, y)
                row["random_block_loglog_slope"] = float(np.polyfit(x, y, 1)[0])
        if not within.empty:
            row.update(
                {
                    "within_parent_samples": int(len(within)),
                    "within_parent_error_ratio_min": float(within["one_minus_overlap_ratio_to_signed_reference"].min()),
                    "within_parent_error_ratio_median": float(within["one_minus_overlap_ratio_to_signed_reference"].median()),
                    "within_parent_error_ratio_max": float(within["one_minus_overlap_ratio_to_signed_reference"].max()),
                    "within_parent_bch_ratio_min": float(within["bch_norm_ratio_to_signed_reference"].min()),
                    "within_parent_bch_ratio_max": float(within["bch_norm_ratio_to_signed_reference"].max()),
                }
            )
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "ablation_summary.csv", index=False)

    if not summary.empty:
        for _, row in summary.iterrows():
            molecule = str(row["molecule"])
            if "random_block_error_ratio_median" in row and pd.notna(row.get("random_block_error_ratio_median")):
                conclusions.append(
                    f"Ablation {molecule}: randomizing whole parent blocks gives median error ratio "
                    f"{float(row['random_block_error_ratio_median']):.3g} relative to signed fermionic."
                )
            if "within_parent_error_ratio_median" in row and pd.notna(row.get("within_parent_error_ratio_median")):
                conclusions.append(
                    f"Ablation {molecule}: within-parent shuffling gives median error ratio "
                    f"{float(row['within_parent_error_ratio_median']):.6g}."
                )
    return conclusions


def summarize_step_scaling(df: pd.DataFrame, outdir: Path, floor: float) -> list[str]:
    if df.empty:
        return ["Step-scaling CSV not found; step-scaling summary skipped."]

    rows: list[dict[str, float | str | int]] = []
    for (case_id, molecule, schedule), group in df.groupby(
        ["case_id", "molecule", "schedule"], sort=True
    ):
        group = group.sort_values("trotter_steps")
        slope, r2, nfit = fit_log_slope(
            group["trotter_steps"], group["one_minus_overlap"], floor
        )
        rows.append(
            {
                "case_id": case_id,
                "molecule": molecule,
                "schedule": schedule,
                "n_points": int(len(group)),
                "n_fit_points": nfit,
                "step_min": int(group["trotter_steps"].min()),
                "step_max": int(group["trotter_steps"].max()),
                "error_at_min_steps": float(group.iloc[0]["one_minus_overlap"]),
                "error_at_max_steps": float(group.iloc[-1]["one_minus_overlap"]),
                "log_error_vs_log_steps_slope": slope,
                "fit_r2": r2,
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "step_scaling_summary.csv", index=False)

    # One separate figure per molecule/case (no subplots).
    for (case_id, molecule), group in df.groupby(["case_id", "molecule"], sort=True):
        fig, ax = plt.subplots()
        for schedule, sched_group in group.groupby("schedule", sort=True):
            sched_group = sched_group.sort_values("trotter_steps")
            mask = sched_group["one_minus_overlap"] > floor
            if not mask.any():
                continue
            ax.loglog(
                sched_group.loc[mask, "trotter_steps"],
                sched_group.loc[mask, "one_minus_overlap"],
                marker="o",
                label=SCHEDULE_LABELS.get(schedule, schedule),
            )
        ax.set_xlabel("First-order Trotter steps")
        ax.set_ylabel("1 - |overlap|")
        ax.set_title(f"Step scaling: {molecule} ({case_id})")
        ax.legend()
        fig.tight_layout()
        safe_case = str(case_id).replace("/", "_")
        fig.savefig(outdir / f"step_scaling_{safe_case}.png", dpi=200)
        plt.close(fig)

    conclusions: list[str] = []
    for molecule, group in summary.groupby("molecule", sort=True):
        ref = group[group["schedule"] == REFERENCE]
        jw = group[group["schedule"] == JW]
        if not ref.empty:
            r = ref.iloc[0]
            conclusions.append(
                f"Step scaling {molecule}: signed-fermionic fitted slope = "
                f"{float(r['log_error_vs_log_steps_slope']):.3f} "
                f"(R^2={float(r['fit_r2']):.3f}, {int(r['n_fit_points'])} points above floor)."
            )
        if not jw.empty:
            j = jw.iloc[0]
            conclusions.append(
                f"Step scaling {molecule}: JW-signed fitted slope = "
                f"{float(j['log_error_vs_log_steps_slope']):.3f}."
            )
    return conclusions


def summarize_active_space(df: pd.DataFrame, outdir: Path) -> list[str]:
    if df.empty:
        return ["Active-space CSV not found; active-space summary skipped."]

    # Expect one row per case/schedule for t=1, r=100.
    key_cols = [
        "case_id",
        "molecule",
        "basis",
        "active_occupied",
        "active_vacant",
        "n_qubits",
    ]
    compact = df[
        key_cols
        + [
            "schedule",
            "one_minus_overlap",
            "bch2_hf_state_norm",
            "bch_cancellation_ratio",
            "weighted_orientation_flip_fraction_vs_signed_reference",
        ]
    ].copy()

    err = compact.pivot_table(index=key_cols, columns="schedule", values="one_minus_overlap", aggfunc="first")
    bch = compact.pivot_table(index=key_cols, columns="schedule", values="bch2_hf_state_norm", aggfunc="first")
    cancel = compact.pivot_table(index=key_cols, columns="schedule", values="bch_cancellation_ratio", aggfunc="first")

    records: list[dict[str, float | str | int | bool]] = []
    for index in err.index:
        record = dict(zip(key_cols, index))
        if REFERENCE not in err.columns or JW not in err.columns:
            continue
        signed_error = float(err.loc[index, REFERENCE])
        jw_error = float(err.loc[index, JW])
        signed_bch = float(bch.loc[index, REFERENCE])
        jw_bch = float(bch.loc[index, JW])
        record.update(
            {
                "signed_error": signed_error,
                "jw_signed_error": jw_error,
                "signed_advantage_factor_vs_jw": jw_error / signed_error if signed_error > 0.0 else math.nan,
                "signed_beats_jw": bool(signed_error < jw_error),
                "signed_bch": signed_bch,
                "jw_signed_bch": jw_bch,
                "bch_advantage_factor_vs_jw": jw_bch / signed_bch if signed_bch > 0.0 else math.nan,
                "signed_cancellation_ratio": float(cancel.loc[index, REFERENCE]),
                "jw_signed_cancellation_ratio": float(cancel.loc[index, JW]),
            }
        )
        if FMAG in err.columns:
            record["fermionic_magnitude_error_ratio_to_signed"] = float(err.loc[index, FMAG]) / signed_error
        if ROUND in err.columns:
            record["round_robin_error_ratio_to_signed"] = float(err.loc[index, ROUND]) / signed_error
        records.append(record)

    summary = pd.DataFrame(records).sort_values(["molecule", "n_qubits"])
    summary.to_csv(outdir / "active_space_summary.csv", index=False)

    mol_rows: list[dict[str, float | str | int]] = []
    for molecule, group in summary.groupby("molecule", sort=True):
        advantages = group["signed_advantage_factor_vs_jw"].to_numpy(dtype=float)
        mol_rows.append(
            {
                "molecule": molecule,
                "cases": int(len(group)),
                "signed_wins_vs_jw": int(group["signed_beats_jw"].sum()),
                "signed_losses_vs_jw": int((~group["signed_beats_jw"]).sum()),
                "median_signed_advantage_factor_vs_jw": float(np.median(advantages)),
                "max_signed_advantage_factor_vs_jw": float(np.max(advantages)),
            }
        )
    mol_summary = pd.DataFrame(mol_rows)
    mol_summary.to_csv(outdir / "active_space_by_molecule.csv", index=False)

    valid = summary[
        np.isfinite(summary["signed_advantage_factor_vs_jw"])
        & np.isfinite(summary["bch_advantage_factor_vs_jw"])
        & (summary["signed_advantage_factor_vs_jw"] > 0.0)
        & (summary["bch_advantage_factor_vs_jw"] > 0.0)
    ]
    if len(valid) >= 3:
        x = np.log10(valid["bch_advantage_factor_vs_jw"].to_numpy(dtype=float))
        y = np.log10(valid["signed_advantage_factor_vs_jw"].to_numpy(dtype=float))
        pearson = correlation(x, y)
        slope = float(np.polyfit(x, y, 1)[0])
    else:
        pearson = math.nan
        slope = math.nan

    # Separate scatter figure.
    if not valid.empty:
        fig, ax = plt.subplots()
        for molecule, group in valid.groupby("molecule", sort=True):
            ax.scatter(
                group["bch_advantage_factor_vs_jw"],
                group["signed_advantage_factor_vs_jw"],
                label=molecule,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(1.0, linewidth=1.0)
        ax.axvline(1.0, linewidth=1.0)
        ax.set_xlabel("BCH advantage: JW BCH / fermionic-signed BCH")
        ax.set_ylabel("Error advantage: JW error / fermionic-signed error")
        ax.set_title("Does BCH cancellation identify favorable fermionic cases?")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "active_space_bch_vs_error_advantage.png", dpi=200)
        plt.close(fig)

    conclusions = [
        f"Active-space sweep: {int(summary['signed_beats_jw'].sum())}/{len(summary)} cases have lower error with signed fermionic than JW signed.",
        f"Across active spaces, Pearson correlation between log BCH advantage and log error advantage = {pearson:.3f}; fitted slope = {slope:.3f}.",
    ]
    for _, row in mol_summary.iterrows():
        conclusions.append(
            f"Active-space {row['molecule']}: signed fermionic wins {int(row['signed_wins_vs_jw'])}/{int(row['cases'])}; "
            f"median JW/signed error ratio = {float(row['median_signed_advantage_factor_vs_jw']):.3g}."
        )
    return conclusions


def summarize_be2_time(df: pd.DataFrame, outdir: Path) -> list[str]:
    if df.empty:
        return ["Be2 time-scan CSV not found; time-scan summary skipped."]

    group_cols = ["evolution_time", "trotter_steps", "trotter_dt"]
    err = df.pivot_table(index=group_cols, columns="schedule", values="one_minus_overlap", aggfunc="first").reset_index()
    bch = df.pivot_table(index=group_cols, columns="schedule", values="bch2_hf_state_norm", aggfunc="first").reset_index()

    summary = err.copy()
    if REFERENCE in err.columns:
        for schedule in [FMAG, JW, ROUND]:
            if schedule in err.columns:
                summary[f"{schedule}_error_ratio_to_signed"] = err[schedule] / err[REFERENCE]
    schedule_cols = [s for s in [REFERENCE, FMAG, JW, ROUND] if s in err.columns]
    if schedule_cols:
        summary["best_schedule"] = err[schedule_cols].idxmin(axis=1)
        summary["best_error"] = err[schedule_cols].min(axis=1)

    # Add BCH values with suffix.
    for schedule in schedule_cols:
        summary[f"{schedule}_bch"] = bch[schedule]

    summary = summary.sort_values("evolution_time")
    summary.to_csv(outdir / "be2_fixed_dt_time_summary.csv", index=False)

    fig, ax = plt.subplots()
    for schedule in schedule_cols:
        mask = summary[schedule] > 0.0
        ax.loglog(
            summary.loc[mask, "evolution_time"],
            summary.loc[mask, schedule],
            marker="o",
            label=SCHEDULE_LABELS.get(schedule, schedule),
        )
    ax.set_xlabel("Total evolution time")
    ax.set_ylabel("1 - |overlap|")
    ax.set_title("Be2 fixed-dt trajectory test (dt = 0.01)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "be2_fixed_dt_time_scan.png", dpi=200)
    plt.close(fig)

    conclusions: list[str] = []
    for _, row in summary.iterrows():
        conclusions.append(
            f"Be2 t={float(row['evolution_time']):.2f}, r={int(row['trotter_steps'])}: best schedule = "
            f"{SCHEDULE_LABELS.get(str(row['best_schedule']), str(row['best_schedule']))}; "
            f"best error = {float(row['best_error']):.3e}."
        )
    if ROUND in summary.columns and REFERENCE in summary.columns:
        ratio = summary[ROUND] / summary[REFERENCE]
        conclusions.append(
            "Be2 round-robin/signed error ratios across increasing time: "
            + ", ".join(f"{value:.3g}" for value in ratio.to_numpy(dtype=float))
            + "."
        )
    return conclusions


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    ablation = read_success(args.ablation)
    step = read_success(args.step_scaling)
    active = read_success(args.active_space)
    be2 = read_success(args.be2_time)

    conclusions: list[str] = []
    conclusions.extend(summarize_ablation(ablation, args.outdir))
    conclusions.extend(summarize_step_scaling(step, args.outdir, args.error_floor))
    conclusions.extend(summarize_active_space(active, args.outdir))
    conclusions.extend(summarize_be2_time(be2, args.outdir))

    report = args.outdir / "final_validation_findings.txt"
    report.write_text("\n".join(f"- {line}" for line in conclusions) + "\n", encoding="utf-8")

    print()
    print("Final-validation summaries written to:")
    print(f"  {args.outdir}")
    print()
    print("Key generated files:")
    for path in sorted(args.outdir.iterdir()):
        print(f"  {path}")
    print()
    print("Read final_validation_findings.txt first, then inspect the summary CSVs and figures.")


if __name__ == "__main__":
    main()
