#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm


CSV_PATH = Path("analysis/l_sweep_trotter_results.csv")

ERROR_COL = "state_infidelity"
CASE_COL = "case_id"
ORDERING_COL = "ordering"
R_COL = "trotter_steps"

FORMULA_ORDER = 1

ORDERINGS = {
    "jw_raw": "JW raw baseline",
    "jw_coloring": "JW coloring",
    "fermionic_coloring": "Fermionic coloring",
}


# ----------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------

df = pd.read_csv(CSV_PATH)

df = df[
    (df["status"] == "success")
    & (df["formula_order"] == FORMULA_ORDER)
    & (df[ORDERING_COL].isin(ORDERINGS))
].copy()

df = df[
    np.isfinite(df[ERROR_COL])
    & (df[ERROR_COL] > 0.0)
    & np.isfinite(df[R_COL])
    & (df[R_COL] > 0)
].copy()


# ----------------------------------------------------------------------
# Fit log(error) = intercept + slope * log(r)
#
# One slope per:
#     case + ordering
# ----------------------------------------------------------------------

rows = []

for (case_id, ordering), group in df.groupby(
    [CASE_COL, ORDERING_COL]
):
    group = (
        group.sort_values(R_COL)
        .drop_duplicates(R_COL)
    )

    # Need at least 3 r values for a meaningful fit.
    if len(group) < 3:
        continue

    r = group[R_COL].to_numpy(dtype=float)
    error = group[ERROR_COL].to_numpy(dtype=float)

    good = (
        np.isfinite(r)
        & np.isfinite(error)
        & (r > 0)
        & (error > 0)
    )

    r = r[good]
    error = error[good]

    if len(r) < 3:
        continue

    slope, intercept = np.polyfit(
        np.log10(r),
        np.log10(error),
        1,
    )

    rows.append(
        {
            "case_id": case_id,
            "ordering": ordering,
            "slope": slope,
            "intercept": intercept,
            "n_points": len(r),
        }
    )


slopes = pd.DataFrame(rows)

print()
print("Slope summary:")
print(
    slopes.groupby("ordering")["slope"]
    .agg(["count", "mean", "median", "std"])
)
print()


# ----------------------------------------------------------------------
# Plot distributions
# ----------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 5.5))

all_slopes = slopes["slope"].to_numpy()

xmin = np.nanpercentile(all_slopes, 1)
xmax = np.nanpercentile(all_slopes, 99)

# Give the graph a little extra room.
padding = 0.15 * (xmax - xmin)
xmin -= padding
xmax += padding

x = np.linspace(xmin, xmax, 500)

for ordering, label in ORDERINGS.items():
    values = slopes.loc[
        slopes["ordering"] == ordering,
        "slope",
    ].to_numpy()

    if len(values) < 2:
        continue

    # Histogram shown as a probability density.
    ax.hist(
        values,
        bins=20,
        density=True,
        alpha=0.20,
        label=f"{label} histogram",
    )

    # Gaussian fitted from empirical mean and standard deviation.
    mu = np.mean(values)
    sigma = np.std(values, ddof=1)

    if sigma > 0:
        gaussian = norm.pdf(x, loc=mu, scale=sigma)

        ax.plot(
            x,
            gaussian,
            linewidth=2,
            label=(
                f"{label}: "
                rf"$\mu={mu:.3f},\ \sigma={sigma:.3f}$"
            ),
        )


# ----------------------------------------------------------------------
# Theoretical prediction
# ----------------------------------------------------------------------

expected_slope = -2 * FORMULA_ORDER

ax.axvline(
    expected_slope,
    linestyle="--",
    linewidth=2,
    label=rf"Theory: slope = ${expected_slope}$",
)

ax.set_xlabel(
    r"Fitted slope $s$ from $E(r)\propto r^s$"
)
ax.set_ylabel("Probability density")

ax.set_title(
    "L-sweep distribution of Trotter-error scaling exponents"
)

ax.grid(True, alpha=0.25)
ax.legend(fontsize=9)

fig.tight_layout()

output = Path(
    "analysis/l_sweep_state_infidelity_slope_distribution.png"
)

fig.savefig(output, dpi=300)

print(f"Saved: {output}")