#!/usr/bin/env python3
"""Analyze the Be2 cross-time coherence benchmark across HGBS-5 and STO-6G."""

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
        "--summary",
        type=Path,
        default=Path("analysis/be2_cross_time_coherence_summary.csv"),
    )
    parser.add_argument(
        "--trace",
        type=Path,
        default=Path("analysis/be2_cross_time_coherence_trace.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/be2_cross_time_coherence_figures"),
    )
    parser.add_argument(
        "--findings",
        type=Path,
        default=Path("analysis/be2_cross_time_coherence_findings.txt"),
    )
    parser.add_argument(
        "--comparison-output",
        type=Path,
        default=Path("analysis/be2_cross_time_coherence_comparison.csv"),
    )
    return parser.parse_args()


def clean_basis(value: str) -> str:
    return str(value).upper().replace("HGBS-5", "HGBS-5").replace("STO-6G", "STO-6G")


def make_schedule_plot(
    df: pd.DataFrame,
    *,
    basis: str,
    metric: str,
    ylabel: str,
    filename: Path,
    logy: bool,
) -> None:
    group = df[df["basis"] == basis].copy()
    if group.empty:
        return
    fig, ax = plt.subplots()
    for schedule, sg in group.groupby("schedule", sort=True):
        sg = sg.sort_values("trotter_steps")
        if logy:
            ax.loglog(
                sg["trotter_steps"],
                sg[metric],
                marker="o",
                label=LABELS.get(schedule, schedule),
            )
        else:
            ax.semilogx(
                sg["trotter_steps"],
                sg[metric],
                marker="o",
                label=LABELS.get(schedule, schedule),
            )
    ax.set_xlabel("Trotter steps")
    ax.set_ylabel(ylabel)
    ax.set_title(f"Be2/{basis}/20q: cross-time Trotter coherence")
    ax.legend()
    fig.tight_layout()
    filename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filename, dpi=220)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.summary)
    if "status" in df.columns:
        df = df[df["status"] == "success"].copy()
    if df.empty:
        raise RuntimeError("No successful summary rows found.")

    numeric = [
        "trotter_steps",
        "actual_one_minus_overlap",
        "sum_local_defect_norm2",
        "vector_interference_term",
        "vector_interference_fraction",
        "vector_cancellation_ratio",
        "overlap_contribution_cancellation_ratio",
        "real_overlap_loss_cancellation_ratio",
        "overlap_reconstruction_abs_error",
    ]
    for column in numeric:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for basis in sorted(df["basis"].dropna().unique()):
        safe = str(basis).lower().replace("-", "").replace(" ", "_")
        make_schedule_plot(
            df,
            basis=basis,
            metric="actual_one_minus_overlap",
            ylabel=r"Final $1-|\langle\psi_{exact}|\psi_{Trotter}\rangle|$",
            filename=args.output_dir / f"{safe}_final_error_vs_steps.png",
            logy=True,
        )
        make_schedule_plot(
            df,
            basis=basis,
            metric="vector_cancellation_ratio",
            ylabel=r"Vector cancellation ratio $\|\Delta\|^2/\sum_j\|\delta_j\|^2$",
            filename=args.output_dir / f"{safe}_vector_cancellation_vs_steps.png",
            logy=True,
        )
        make_schedule_plot(
            df,
            basis=basis,
            metric="overlap_contribution_cancellation_ratio",
            ylabel=r"Overlap cancellation ratio $|\sum_j a_j|/\sum_j|a_j|$",
            filename=args.output_dir / f"{safe}_overlap_cancellation_vs_steps.png",
            logy=True,
        )

    # Build a compact signed-vs-round-robin comparison for every basis/step count.
    rows = []
    for (basis, steps), group in df.groupby(["basis", "trotter_steps"], sort=True):
        indexed = group.set_index("schedule")
        if REFERENCE not in indexed.index or ROUND not in indexed.index:
            continue
        signed = indexed.loc[REFERENCE]
        rr = indexed.loc[ROUND]
        rows.append(
            {
                "basis": basis,
                "trotter_steps": int(steps),
                "signed_final_error": signed["actual_one_minus_overlap"],
                "round_robin_final_error": rr["actual_one_minus_overlap"],
                "round_robin_final_error_ratio_to_signed": (
                    rr["actual_one_minus_overlap"] / signed["actual_one_minus_overlap"]
                ),
                "signed_sum_local_defect_norm2": signed["sum_local_defect_norm2"],
                "round_robin_sum_local_defect_norm2": rr["sum_local_defect_norm2"],
                "round_robin_local_defect_ratio_to_signed": (
                    rr["sum_local_defect_norm2"] / signed["sum_local_defect_norm2"]
                ),
                "signed_vector_cancellation_ratio": signed["vector_cancellation_ratio"],
                "round_robin_vector_cancellation_ratio": rr["vector_cancellation_ratio"],
                "signed_overlap_cancellation_ratio": signed[
                    "overlap_contribution_cancellation_ratio"
                ],
                "round_robin_overlap_cancellation_ratio": rr[
                    "overlap_contribution_cancellation_ratio"
                ],
                "signed_vector_interference_fraction": signed[
                    "vector_interference_fraction"
                ],
                "round_robin_vector_interference_fraction": rr[
                    "vector_interference_fraction"
                ],
            }
        )
    comparison = pd.DataFrame(rows)
    args.comparison_output.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.comparison_output, index=False)

    lines = []
    lines.append("Be2 cross-time coherence findings")
    lines.append("=" * 38)
    lines.append("")
    lines.append(
        "Interpretation rule: a vector_cancellation_ratio below 1 means the "
        "cross-time interference term is destructive; the smaller the ratio, "
        "the stronger the net cancellation relative to the available local "
        "defect power."
    )
    lines.append("")

    for row in comparison.itertuples(index=False):
        lines.append(
            f"{row.basis}, steps={row.trotter_steps}: "
            f"RR/f-signed final error={row.round_robin_final_error_ratio_to_signed:.6g}; "
            f"RR/f-signed local-defect power={row.round_robin_local_defect_ratio_to_signed:.6g}; "
            f"vector cancellation signed={row.signed_vector_cancellation_ratio:.6g}, "
            f"RR={row.round_robin_vector_cancellation_ratio:.6g}; "
            f"overlap cancellation signed={row.signed_overlap_cancellation_ratio:.6g}, "
            f"RR={row.round_robin_overlap_cancellation_ratio:.6g}."
        )
    lines.append("")

    # Report whether the key mechanism pattern occurs at each basis's 100-step row.
    for basis in sorted(comparison["basis"].unique()) if not comparison.empty else []:
        bg = comparison[comparison["basis"] == basis]
        target = bg.iloc[(bg["trotter_steps"] - 100).abs().argsort()[:1]]
        if target.empty:
            continue
        row = target.iloc[0]
        rr_local_worse = row["round_robin_local_defect_ratio_to_signed"] > 1.0
        rr_final_better = row["round_robin_final_error_ratio_to_signed"] < 1.0
        rr_cancels_more = (
            row["round_robin_vector_cancellation_ratio"]
            < row["signed_vector_cancellation_ratio"]
        )
        if rr_local_worse and rr_final_better and rr_cancels_more:
            lines.append(
                f"{basis}: DIRECT SUPPORT for the coherence mechanism near "
                f"{int(row['trotter_steps'])} steps.  Round robin has larger "
                "summed local defects but smaller final error, together with "
                "stronger destructive cross-time vector cancellation."
            )
        elif rr_local_worse and rr_final_better:
            lines.append(
                f"{basis}: round robin is locally worse but globally better, "
                "but the vector-cancellation diagnostic is not smaller.  "
                "Inspect the overlap-contribution decomposition and trace."
            )
        else:
            lines.append(
                f"{basis}: the 100-step behavior does not show the same "
                "locally-worse/globally-better pattern; this is a useful "
                "basis-dependent control rather than a failure."
            )

    max_recon = float(df["overlap_reconstruction_abs_error"].max())
    lines.append("")
    lines.append(f"Maximum telescoping overlap reconstruction error: {max_recon:.3e}")

    args.findings.parent.mkdir(parents=True, exist_ok=True)
    args.findings.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(comparison.to_string(index=False))
    print()
    print("\n".join(lines[-6:]))
    print()
    print(f"Wrote: {args.comparison_output}")
    print(f"Wrote: {args.findings}")
    print(f"Figures: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
