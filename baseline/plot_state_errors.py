#!/usr/bin/env python3
"""
plot_state_errors.py
--------------------
Trotter state error from compute_trotter_state_error.py, one 3x3 grid per
molecule per evolution time T.

Grid layout, for a single (molecule, basis, mapping, T):
    columns = Trotter order (1, 2, 4)
    rows    = n_steps per U, i.e. delta_t = 1/n_steps
    x axis  = qubit count, annotated with the active space
    y axis  = state error, log scale
    lines   = bond length L (single atoms have none and plot as one "atom" line)

T is now a swept axis rather than a per-Hamiltonian lookup, so it is constant
across a figure and named in the title instead of on the x axis. Basis and
mapping are read from the data, not the path.

Takes either one CSV or a <basis>/<mapping>/ directory, in which case every
state_errors_*.csv in it is plotted, giving n_molecules x n_T figures.

Usage:
    python plot_state_errors.py state_errors/hgbs-5/jw
    python plot_state_errors.py state_errors/hgbs-5/jw --T 1024
    python plot_state_errors.py state_errors/hgbs-5/jw/state_errors_Li-Li.csv
    python plot_state_errors.py <dir> --orders 2,4 --n-steps 1,2 --cmap plasma
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


def one_figure(df, molecule, basis, mapping, T, orders, steps, args):
    qubits = sorted(int(q) for q in df["n_qubits"].unique())
    keys = sorted(df["Lkey"].unique(),
                  key=lambda s: (s == "atom",
                                 float(s[2:]) if s != "atom" else 0.0))
    active_of = (df.drop_duplicates("n_qubits")
                   .set_index("n_qubits")["active"].to_dict())

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
                s = sub[sub["Lkey"] == k].sort_values("n_qubits")
                if s.empty:
                    continue
                ax.plot(s["n_qubits"], s["err"], marker="o", ms=4, lw=1.3,
                        color=colors[k], label=k)

            if args.logx:
                ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.25, lw=0.5)

            if r == 0:
                ax.set_title(f"order {order}", fontsize=12)

            # n_steps is PER U; the meaningful step size is delta_t = 1/n_steps
            # and total_steps = T/delta_t sets the cost. dt_norm = delta_t*sum|c|
            # is the DIMENSIONLESS step, and it is only constant across active
            # spaces where the shifted one-norm exceeds 1. Below that,
            # normalization clamps to 2 and dt_norm tracks lambda instead, so
            # the panel's Hamiltonians are not taking the same step in units of
            # ||H||. Showing the spread makes that visible rather than implicit.
            dt = 1.0 / ns
            note = (f"$n_{{steps}}/U$ = {ns}\n$\\Delta t$ = {dt:g}\n"
                    f"steps = {int(T * ns):,}")
            if "dt_norm" in sub.columns and not sub.empty:
                lo, hi = sub["dt_norm"].min(), sub["dt_norm"].max()
                # note += ("\n$\\Delta t\\,\\Sigma|c|$ = "
                #          + (f"{lo:.3f}" if hi - lo < 1e-9
                #             else f"{lo:.3f}-{hi:.3f}"))
            ax.text(0.42, 0.26, note,
                    transform=ax.transAxes, fontsize=9.5, va="top", ha="left",
                    bbox=dict(fc="white", ec="0.6", alpha=0.9,
                              boxstyle="round,pad=0.3"))

            if c == 0:
                ax.set_ylabel("state error\n"
                              r"$\|\psi_{trot}-\psi_{exact}\|$", fontsize=11)

            ax.set_xticks(qubits)
            ax.set_xticklabels([str(q) for q in qubits], fontsize=11)
            if args.logx:
                ax.xaxis.set_minor_locator(ticker.NullLocator())
                ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
            if r == nrow - 1:
                ax.set_xlabel("(qubits, active space)")
                ax.set_xticklabels(
                    [f"{q}\n{active_of.get(q, '')}" for q in qubits],
                    fontsize=10)

    fig.suptitle(f"Trotter state error   {molecule}   {basis}   "
                 f"{mapping.upper()}   T = {T:,}", fontsize=13, y=0.94)

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
    ap.add_argument("--T", default=None, dest="T_values",
                    help="comma-separated T to plot (default: all present)")
    ap.add_argument("--orders", default=None,
                    help="comma-separated orders (default: all present)")
    ap.add_argument("--n-steps", default=None, dest="n_steps",
                    help="comma-separated n_steps per U (default: all present)")
    ap.add_argument("--logx", action="store_true",
                    help="log-scale the qubit axis. Off by default: 4-12 qubits "
                         "is under half a decade and no power law in qubit "
                         "count is expected, so it changes little.")
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
        Ts = ([int(t) for t in args.T_values.split(",")] if args.T_values
              else sorted(int(t) for t in df["T"].unique()))

        df = df[df["order"].isin(orders) & df["n_steps"].isin(steps)]
        if df.empty:
            print(f"  {path.name}: no rows match the requested orders/n_steps")
            continue

        molecule = df["molecule"].iloc[0]
        basis = df["basis"].iloc[0]
        mapping = df["mapping"].iloc[0]
        out_dir = Path(args.out_dir) if args.out_dir else path.parent
        out_dir.mkdir(parents=True, exist_ok=True)

        for T in Ts:
            d = df[df["T"] == T]
            if d.empty:
                continue
            fig = one_figure(d, molecule, basis, mapping, T, orders, steps, args)
            out = out_dir / (f"state_errors_{molecule}_{basis}_{mapping}"
                             f"_T{T}.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)
            n_fig += 1
            msg = (f"wrote {out}   ({len(d)} rows, {d['Lkey'].nunique()} L, "
                   f"qubits {sorted(int(q) for q in d['n_qubits'].unique())})")
            if "dt_norm" in d.columns:
                at1 = d[d["n_steps"] == d["n_steps"].min()]["dt_norm"]
                if not at1.empty and at1.max() - at1.min() > 1e-9:
                    msg += (f"\n    NOTE dt_norm spans {at1.min():.3f}-"
                            f"{at1.max():.3f} at the coarsest step: active "
                            f"spaces are NOT taking the same dimensionless "
                            f"step (one-norm clamp)")
            print(msg)

    print(f"\n{n_fig} figures written.")


if __name__ == "__main__":
    main()
