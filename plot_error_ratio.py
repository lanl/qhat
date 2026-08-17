#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

CSV_PATH = Path("analysis/l_sweep_20case_fermionic_comparison.csv")

# Change these to the two ordering names you want to compare.
METHOD_A = "fermionic_coloring"
METHOD_B = "jw_coloring"

ERROR_COL = "state_infidelity"
R_COL = "trotter_steps"
CASE_COL = "case_id"


# ----------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------

df = pd.read_csv(CSV_PATH)

if "status" in df.columns:
    df = df[df["status"] == "success"].copy()

df = df[
    df["ordering"].isin([METHOD_A, METHOD_B])
].copy()

# Remove unusable values for a log-ratio plot.
df = df[
    np.isfinite(df[ERROR_COL])
    & (df[ERROR_COL] > 0)
    & np.isfinite(df[R_COL])
    & (df[R_COL] > 0)
].copy()


# ----------------------------------------------------------------------
# Put A and B next to each other for the same case and r
# ----------------------------------------------------------------------

wide = (
    df.pivot_table(
        index=[CASE_COL, R_COL],
        columns="ordering",
        values=ERROR_COL,
        aggfunc="first",
    )
    .dropna(subset=[METHOD_A, METHOD_B])
    .reset_index()
)

wide["error_ratio"] = wide[METHOD_A] / wide[METHOD_B]
wide["log10_ratio"] = np.log10(wide["error_ratio"])


# ----------------------------------------------------------------------
# Aggregate across all L-sweep cases
#
# Median in log space:
#
#   10^[ median(log10(E_A / E_B)) ]
#
# This is a natural central value for multiplicative ratios.
# ----------------------------------------------------------------------

summary = (
    wide.groupby(R_COL)["log10_ratio"]
    .agg(
        median="median",
        q25=lambda x: x.quantile(0.25),
        q75=lambda x: x.quantile(0.75),
        count="count",
    )
    .reset_index()
)

summary["ratio_median"] = 10.0 ** summary["median"]
summary["ratio_q25"] = 10.0 ** summary["q25"]
summary["ratio_q75"] = 10.0 ** summary["q75"]


# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 6))

# Individual L-sweep cases
for case_id, case_df in wide.groupby(CASE_COL):
    case_df = case_df.sort_values(R_COL)

    ax.plot(
        case_df[R_COL],
        case_df["error_ratio"],
        marker="o",
        markersize=3,
        linewidth=0.8,
        alpha=0.18,
    )

# Interquartile range across cases
ax.fill_between(
    summary[R_COL],
    summary["ratio_q25"],
    summary["ratio_q75"],
    alpha=0.2,
    label="L-sweep IQR",
)

# Central L-sweep trend
ax.plot(
    summary[R_COL],
    summary["ratio_median"],
    marker="o",
    linewidth=2.5,
    label="Geometric median across cases",
)

# Equal-error reference
ax.axhline(
    1.0,
    linestyle="--",
    linewidth=1.5,
    label="Equal error",
)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Trotter steps $r$")
ax.set_ylabel(
    rf"$E_{{{METHOD_A}}} / E_{{{METHOD_B}}}$"
)

ax.set_title(
    "Trotter-error ratio vs. number of Trotter steps"
)

ax.grid(True, which="both", alpha=0.25)
ax.legend()

fig.tight_layout()

output = Path(
    f"analysis/error_ratio_{METHOD_A}_vs_{METHOD_B}.png"
)
fig.savefig(output, dpi=300)

print(f"Saved: {output}")
print()
print(summary.to_string(index=False))