#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CSV = Path("analysis/nh3_scaling_demo.csv")

df = pd.read_csv(CSV)

df = df[
    (df["status"] == "success")
    & (df["formula_order"] == 1)
].copy()

labels = {
    "jw_raw": "JW raw baseline",
    "jw_coloring": "JW coloring",
    "fermionic_coloring": "Fermionic coloring",
}

fig, ax = plt.subplots(figsize=(7.2, 5.2))

for ordering, label in labels.items():
    g = (
        df[df["ordering"] == ordering]
        .sort_values("trotter_steps")
        .copy()
    )

    r = g["trotter_steps"].to_numpy(dtype=float)
    error = g["state_infidelity"].to_numpy(dtype=float)

    good = np.isfinite(error) & (error > 0.0)
    r = r[good]
    error = error[good]

    slope, intercept = np.polyfit(
        np.log10(r),
        np.log10(error),
        1,
    )

    ax.loglog(
        r,
        error,
        marker="o",
        linewidth=2,
        label=f"{label}: slope = {slope:.3f}",
    )


# Add theoretical r^-2 guide.
raw = (
    df[df["ordering"] == "jw_raw"]
    .sort_values("trotter_steps")
)

r0 = float(raw["trotter_steps"].iloc[0])
e0 = float(raw["state_infidelity"].iloc[0])

r_guide = np.array(
    [
        df["trotter_steps"].min(),
        df["trotter_steps"].max(),
    ],
    dtype=float,
)

guide = e0 * (r0 / r_guide) ** 2

ax.loglog(
    r_guide,
    guide,
    linestyle="--",
    linewidth=1.5,
    label=r"Theory: $r^{-2}$",
)

ax.set_xlabel("Trotter steps $r$")
ax.set_ylabel(
    r"State infidelity "
    r"$1-|\langle\psi_{\rm exact}|\psi_{\rm Trot}\rangle|^2$"
)

ax.set_title(
    "NH$_3$ / STO-6G, active space (2 occ, 4 vac), $T=1$"
)

ax.grid(True, which="both", alpha=0.25)
ax.legend()

fig.tight_layout()

output = Path("analysis/nh3_trotter_scaling.png")
fig.savefig(output, dpi=300)

print(f"Saved: {output}")