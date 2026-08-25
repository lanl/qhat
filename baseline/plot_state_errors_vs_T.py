#!/usr/bin/env python3
r"""
plot_state_errors_vs_T.py
-------------------------
Trotter state error against evolution time, one 3x3 grid per molecule per
active space.

Grid layout, for a single (molecule, basis, mapping, active space):
    columns = Trotter order (1, 2, 4)
    rows    = n_steps per U, i.e. delta_t = 1/n_steps
    x axis  = evolution time T, log scale
    y axis  = state error, log scale
    lines   = bond length L (single atoms have none and plot as one "atom" line)

Companion to plot_state_errors.py, which fixes T and puts qubit count on x.
Here the active space is fixed and T is on x, so the figure answers how error
accumulates with evolution time rather than how it varies across active spaces.

Note that error is NOT expected to be monotonic in T at fixed delta_t: Trotter
error accumulates as a coherent sum of rotations that can partially cancel at
particular evolution times. A dip at one T is a real feature of the dynamics,
not a bad data point.

Takes either one CSV or a <basis>/<mapping>/ directory, in which case every
state_errors_*.csv in it is plotted, giving n_molecules x n_active_spaces
figures.

Usage:
    python plot_state_errors_vs_T.py state_errors/hgbs-5/jw
    python plot_state_errors_vs_T.py state_errors/hgbs-5/jw --qubits 12
    python plot_state_errors_vs_T.py state_errors/hgbs-5/jw/state_errors_Li-Li.csv
    python plot_state_errors_vs_T.py <dir> --orders 2,4 --n-steps 1,2 --out-dir figs
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd


def L_label(v):
    """Bond length as a stable group key; single atoms carry no L."""
    return "atom" if pd.isna(v) else f"L={v:g}"


def one_figure(df, molecule, basis, mapping, active, n_qubits, orders, steps,
               args):
    Ts = sorted(int(t) for t in df["T"].unique())
    keys = sorted(df["Lkey"].unique(),
                  key=lambda s: (s == "atom",
                                 float(s[2:]) if s != "atom" else 0.0))

    cmap = plt.get_cmap(args.cmap)
    colors = {k: cmap(i / max(1, len(keys) - 1)) for i, k in enumerate(keys)}

    nrow, ncol = len(steps), len(orders)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 3.6 * nrow),
                             sharex=True, sharey=True, squeeze=False)

    for r, ns in enumerate(steps):
        for c, order in enumerate(orders):
            ax = axes[r][c]
            sub = df[(df["n_steps"] == ns) & (df["order"] == order)]

            for k in keys:
                s = sub[sub["Lkey"] == k].sort_values("T")
                if s.empty:
                    continue
                ax.plot(s["T"], s["err"], marker="o", ms=4, lw=1.3,
                        color=colors[k], label=k)

            if not args.linx:
                ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.25, lw=0.5)

            if r == 0:
                ax.set_title(f"order {order}", fontsize=12)

            # n_steps is PER U; delta_t = 1/n_steps is the step size and
            # total_steps = T/delta_t is what sets the cost. total_steps
            # therefore varies along x within a panel, so it is given as a
            # range rather than a single number.
            dt = 1.0 / ns
            note = (f"$n_{{steps}}/U$ = {ns}\n$\\Delta t$ = {dt:g}")
            ax.text(0.42, 0.26, note,
                    transform=ax.transAxes, fontsize=9.5, va="top", ha="left",
                    bbox=dict(fc="white", ec="0.6", alpha=0.9,
                              boxstyle="round,pad=0.3"))

            if c == 0:
                ax.set_ylabel("state error\n"
                              r"$\|\psi_{trot}-\psi_{exact}\|$", fontsize=11)

            ax.set_xticks(Ts)
            ax.set_xticklabels([f"{t:,}" for t in Ts], fontsize=11)
            if not args.linx:
                ax.xaxis.set_minor_locator(ticker.NullLocator())
                ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
            if r == nrow - 1:
                ax.set_xlabel("evolution time $T$  (applications of $U$)")

    fig.suptitle(f"Trotter state error vs evolution time   {molecule}   "
                 f"{basis}   {mapping.upper()}   {active}  ({n_qubits}q)",
                 fontsize=13, y=0.94)

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title="bond length", loc="upper center",
                   bbox_to_anchor=(0.5, 0.885), ncol=min(len(labels), 10),
                   fontsize=10, title_fontsize=11, frameon=True)

    top = 1.0 - (0.9 / (2.2 * nrow))
    fig.tight_layout(rect=[0, 0, 1, top])
    return fig


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", help="a <basis>/<mapping>/ directory, or one CSV")
    ap.add_argument("--out-dir", default=None,
                    help="where to write PNGs (default: alongside each CSV)")
    ap.add_argument("--qubits", default=None,
                    help="comma-separated qubit counts (default: all present, "
                         "one figure each)")
    ap.add_argument("--active", default=None,
                    help="comma-separated active spaces, e.g. 004-008 "
                         "(default: all present)")
    ap.add_argument("--T", default=None, dest="T_values",
                    help="comma-separated T (default: all present)")
    ap.add_argument("--orders", default=None,
                    help="comma-separated orders (default: all present)")
    ap.add_argument("--n-steps", default=None, dest="n_steps",
                    help="comma-separated n_steps per U (default: all present)")
    ap.add_argument("--linx", action="store_true",
                    help="linear T axis. Default is log base 2, since T is "
                         "swept in powers of two.")
    ap.add_argument("--cmap", default="viridis")
    args = ap.parse_args()

    src = Path(args.source)
    csvs = sorted(src.glob("state_errors_*.csv")) if src.is_dir() else [src]
    if not csvs:
        raise SystemExit(f"no state_errors_*.csv in {src}")

    n_fig = 0
    for path in csvs:
        df = pd.read_csv(path, comment="#")
        if df.empty:
            print(f"  {path.name}: empty, skipped")
            continue
        df["Lkey"] = df["L"].map(L_label)

        orders = ([int(o) for o in args.orders.split(",")] if args.orders
                  else sorted(int(o) for o in df["order"].unique()))
        steps = ([int(s) for s in args.n_steps.split(",")] if args.n_steps
                 else sorted(int(s) for s in df["n_steps"].unique()))

        df = df[df["order"].isin(orders) & df["n_steps"].isin(steps)]
        if args.T_values:
            df = df[df["T"].isin([int(t) for t in args.T_values.split(",")])]
        if args.qubits:
            df = df[df["n_qubits"].isin(
                [int(q) for q in args.qubits.split(",")])]
        if args.active:
            df = df[df["active"].isin(
                [a.strip() for a in args.active.split(",")])]
        if df.empty:
            print(f"  {path.name}: no rows match the requested filters")
            continue

        molecule = df["molecule"].iloc[0]
        basis = df["basis"].iloc[0]
        mapping = df["mapping"].iloc[0]
        out_dir = Path(args.out_dir) if args.out_dir else path.parent
        out_dir.mkdir(parents=True, exist_ok=True)

        spaces = (df.drop_duplicates("active")
                    .sort_values("n_qubits")[["active", "n_qubits"]]
                    .values.tolist())
        for active, n_qubits in spaces:
            d = df[df["active"] == active]
            n_T = d["T"].nunique()
            if n_T < 2:
                print(f"  {molecule} {active}: only {n_T} T value, skipped")
                continue
            fig = one_figure(d, molecule, basis, mapping, active,
                             int(n_qubits), orders, steps, args)
            out = out_dir / (f"state_errors_vsT_{molecule}_{basis}_{mapping}"
                             f"_as{active}.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)
            n_fig += 1
            print(f"wrote {out}   ({len(d)} rows, {d['Lkey'].nunique()} L, "
                  f"T {sorted(int(t) for t in d['T'].unique())})")

    print(f"\n{n_fig} figures written.")


if __name__ == "__main__":
    main()
