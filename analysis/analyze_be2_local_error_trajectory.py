#!/usr/bin/env python3
"""Summarize and plot the Be2 trajectory-aware local-error case study."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REFERENCE = "fermionic_signed_reference"
FMAG = "fermionic_magnitude_reference"
JW = "jw_signed_baseline"
ROUND = "signed_parent_descendants_round_robin"
LABELS = {
    REFERENCE: "Fermionic signed",
    FMAG: "Fermionic magnitude",
    JW: "JW signed",
    ROUND: "Round robin",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("analysis/be2_local_error_trajectory_hgbs5.csv"),
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/be2_local_error_trajectory_summary.csv"),
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=Path("analysis/be2_local_error_trajectory_proxy.png"),
    )
    parser.add_argument(
        "--ablation",
        type=Path,
        default=Path("analysis/fermionic_structure_ablation_hgbs5.csv"),
        help="Optional existing ablation CSV for the t=0 HF-BCH sanity check.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.input)
    if "status" in df.columns:
        df = df[df["status"] == "success"].copy()
    if df.empty:
        raise RuntimeError("No successful trajectory rows found.")

    metric = "trajectory_bch_proxy_2err_over_dt2"
    needed = {"checkpoint_time", "schedule", metric}
    missing = needed - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing columns: {sorted(missing)}")

    rows = []
    for checkpoint, group in df.groupby("checkpoint_time", sort=True):
        values = {
            row.schedule: float(getattr(row, metric))
            for row in group.itertuples(index=False)
        }
        if REFERENCE not in values:
            continue
        ref = values[REFERENCE]
        best_schedule = min(values, key=values.get)
        rows.append(
            {
                "checkpoint_time": float(checkpoint),
                "signed_proxy": ref,
                "fermionic_magnitude_proxy": values.get(FMAG, np.nan),
                "jw_signed_proxy": values.get(JW, np.nan),
                "round_robin_proxy": values.get(ROUND, np.nan),
                "fermionic_magnitude_ratio_to_signed": values.get(FMAG, np.nan) / ref,
                "jw_signed_ratio_to_signed": values.get(JW, np.nan) / ref,
                "round_robin_ratio_to_signed": values.get(ROUND, np.nan) / ref,
                "best_schedule": best_schedule,
                "best_schedule_label": LABELS.get(best_schedule, best_schedule),
            }
        )

    summary = pd.DataFrame(rows).sort_values("checkpoint_time")
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary, index=False)

    fig, ax = plt.subplots()
    for schedule, group in df.groupby("schedule", sort=True):
        group = group.sort_values("checkpoint_time")
        ax.semilogy(
            group["checkpoint_time"],
            group[metric],
            marker="o",
            label=LABELS.get(schedule, schedule),
        )
    ax.set_xlabel("Exact-trajectory checkpoint time")
    ax.set_ylabel(r"Trajectory BCH proxy: $2\|\delta\psi_{local}\|/\Delta t^2$")
    ax.set_title("Be2/HGBS-5/20q: local ordering sensitivity along exact trajectory")
    ax.legend()
    fig.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, dpi=220)
    plt.close(fig)

    rr = summary.dropna(subset=["round_robin_ratio_to_signed"])
    crossing = rr[rr["round_robin_ratio_to_signed"] < 1.0]

    print(summary.to_string(index=False))
    print()
    if not crossing.empty:
        first = crossing.iloc[0]
        print(
            "RESULT: round-robin local error becomes smaller than signed "
            f"fermionic by checkpoint t={first['checkpoint_time']:g}."
        )
        print(
            "Interpretation: the preferred instantaneous ordering changes as "
            "the exact state evolves, directly supporting a trajectory-dependent "
            "explanation of the Be2 finite-time crossover."
        )
    else:
        print(
            "RESULT: round-robin never beats signed fermionic in local error "
            "over the sampled checkpoints."
        )
        print(
            "Interpretation: the t=1 cumulative crossover is then not explained "
            "by local error magnitudes alone; time-integrated vector directions "
            "and cancellation become the leading mechanism to investigate."
        )

    if args.ablation.exists():
        abl = pd.read_csv(args.ablation)
        if "status" in abl.columns:
            abl = abl[abl["status"] == "success"]
        case_id = str(df.iloc[0]["case_id"])
        hf = abl[
            (abl["case_id"] == case_id)
            & (abl["schedule"].isin(df["schedule"].unique()))
        ][["schedule", "bch2_hf_state_norm"]].drop_duplicates("schedule")
        t0 = df[np.isclose(df["checkpoint_time"], 0.0)][
            ["schedule", metric]
        ]
        check = hf.merge(t0, on="schedule", how="inner")
        if not check.empty:
            check["trajectory_proxy_over_hf_bch"] = (
                check[metric] / check["bch2_hf_state_norm"]
            )
            print()
            print("t=0 sanity check: trajectory proxy / existing HF BCH")
            print(check.to_string(index=False))

    print()
    print(f"Wrote: {args.summary}")
    print(f"Wrote: {args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())