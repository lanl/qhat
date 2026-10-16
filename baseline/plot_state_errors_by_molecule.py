#!/usr/bin/env python3
"""
plot_state_errors_by_molecule.py
--------------------------------
State error across the whole library, one panel per molecule.

Grid layout, for a single (basis, mapping, T, order, delta_t):
    one panel per molecule, 3 columns
    x axis  = active space (qubit count as the position, active space as label)
    y axis  = state error, log scale
    colour  = bond length as a fraction of equilibrium, L / L_eq

Because L_min = L_min_frac * L_eq, absolute bond lengths are not comparable
between molecules -- Li-Li's 1.60 A and H-H's 0.38 A are the same point in the
sweep. Colour therefore keys on L / L_eq, recovered as
L_min_frac * L / min(L), which is the same scale in every panel. A colourbar
replaces a per-molecule legend for the same reason.

Panels are grouped by period, taken from the heaviest element in the species,
and each period gets its own rows rather than running into the next. With 3
columns the nine library diatomics divide exactly: H-H, He-H, He-He across the
first row, then the six Li/Be/B species across two more. --no-period-rows packs
them continuously instead.

Complements plot_state_errors.py, which fixes the molecule and grids over order
and delta_t. Here order and delta_t are fixed and the grid is molecules, so the
default invocation gives one figure per T -- three figures for three T.

Single atoms have no bond length; they are dropped with a note, since the
L / L_eq colour axis is meaningless for them.

Usage:
    python plot_state_errors_by_molecule.py state_errors/hgbs-5/jw
    python plot_state_errors_by_molecule.py <dir> --order 4 --n-steps 4
    python plot_state_errors_by_molecule.py <dir> --T 512 --out-dir figs
    python plot_state_errors_by_molecule.py <dir> --order 1,2,4 --n-steps 1,2,4
    
    
    python plot_state_errors_by_molecule.py state_errors/hgbs-5/jw --out-dir plots/active_space_vs_state_error_all_molecule_hgbs_jw
    python plot_state_errors_by_molecule.py state_errors/hgbs-5/bk --out-dir plots/active_space_vs_state_error_all_molecule_hgbs_bk
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# the sweep order used by build_config_L_sweep.py
DEFAULT_MOLECULES = ["H-H", "He-H", "He-He", "Li-H", "Li-Li",
                     "Be-H", "Be-Be", "B-H", "B-B"]

ELEMENT_PERIOD = {"H": 1, "He": 1,
                  "Li": 2, "Be": 2, "B": 2, "C": 2, "N": 2, "O": 2,
                  "F": 2, "Ne": 2,
                  "Na": 3, "Mg": 3, "Al": 3, "Si": 3, "P": 3, "S": 3,
                  "Cl": 3, "Ar": 3}


def species_period(name):
    """
    Period of a species, taken as the heaviest constituent: Li-H is period 2
    because the chemistry is set by the Li. Bare symbols (single atoms) work
    unchanged. Unknown elements sort last rather than raising.
    """
    return max(ELEMENT_PERIOD.get(part, 99) for part in name.split("-"))


def period_layout(molecules, ncol):
    """
    Arrange panels so each period occupies whole rows, padding with blanks
    rather than letting one period run into the next. Returns the padded slot
    list (None for a blank) and a list of (period, first_row, last_row) bands.

    For the nine library diatomics this is exact with ncol=3: period 1 is
    H-H, He-H, He-He (one row) and period 2 is the six Li/Be/B species (two).
    """
    from itertools import groupby
    slots, bands = [], []
    for period, grp in groupby(molecules, key=species_period):
        grp = list(grp)
        first_row = len(slots) // ncol
        slots.extend(grp)
        while len(slots) % ncol:
            slots.append(None)
        bands.append((period, first_row, len(slots) // ncol - 1))
    return slots, bands


def load(source):
    """Concatenate every state_errors_*.csv in a directory, or read one CSV."""
    p = Path(source)
    if p.is_dir():
        files = sorted(p.glob("state_errors_*.csv"))
        if not files:
            raise SystemExit(f"no state_errors_*.csv in {p}")
        frames = []
        for f in files:
            d = pd.read_csv(f, comment="#")
            if not d.empty:
                frames.append(d)
        if not frames:
            raise SystemExit(f"every CSV in {p} was empty")
        print(f"read {len(files)} CSVs from {p}")
        return pd.concat(frames, ignore_index=True)
    return pd.read_csv(p, comment="#")


def one_figure(df, molecules, basis, mapping, T, order, ns, args, norm, cmap):
    ncol = args.ncol
    if args.no_period_rows:
        slots, bands = list(molecules), []
    else:
        slots, bands = period_layout(molecules, ncol)
    nrow = int(np.ceil(len(slots) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow),
                             squeeze=False, sharey=args.sharey,
                             layout="constrained")

    for idx, mol in enumerate(slots):
        ax = axes[idx // ncol][idx % ncol]
        if mol is None:
            ax.axis("off")
            continue
        d = df[df["molecule"] == mol]
        if d.empty:
            ax.set_title(f"{mol}  (no data)", fontsize=10, color="0.5")
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        # active space position on x, ordered by qubit count
        spaces = (d.drop_duplicates("active")
                    .sort_values("n_qubits")[["active", "n_qubits"]]
                    .values.tolist())
        xpos = {a: i for i, (a, q) in enumerate(spaces)}

        for frac, g in d.groupby("L_frac"):
            g = g.copy()
            g["x"] = g["active"].map(xpos)
            g = g.sort_values("x")
            ax.plot(g["x"], g["err"], marker="o", ms=4, lw=1.3,
                    color=cmap(norm(frac)))

        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.25, lw=0.5)
        ax.set_title(mol, fontsize=11)
        ax.set_xticks(range(len(spaces)))
        ax.set_xticklabels([f"{a}\n{q}q" for a, q in spaces], fontsize=8)
        if idx % ncol == 0:
            ax.set_ylabel("state error\n"
                          r"$\|\psi_{trot}-\psi_{exact}\|$", fontsize=10)
        if idx // ncol == nrow - 1:
            ax.set_xlabel("active space", fontsize=10)

    for j in range(len(slots), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    # one label per period band, on the left margin of its first row
    for period, first_row, last_row in bands:
        mid = axes[(first_row + last_row) // 2][0]
        mid.annotate(f"period {period}", xy=(-0.28, 0.5),
                     xycoords="axes fraction", rotation=90,
                     va="center", ha="center", fontsize=12, color="0.35")

    fig.suptitle(f"Trotter state error   {basis}   {mapping.upper()}   "
                 f"T = {T:,}   order {order}   "
                 f"$\\Delta t$ = {1.0/ns:g}   ({T*ns:,} steps)",
                 fontsize=13)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
    cbar.set_label(r"bond length  $L / L_{eq}$", fontsize=10)
    return fig


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", help="a <basis>/<mapping>/ directory, or one CSV")
    ap.add_argument("--out-dir", default=None,
                    help="where to write PNGs (default: alongside the source)")
    ap.add_argument("--T", default=None, dest="T_values",
                    help="comma-separated T (default: all present)")
    ap.add_argument("--order", default="2",
                    help="Trotter order; one figure per value (default: 2)")
    ap.add_argument("--n-steps", default="1", dest="n_steps",
                    help="n_steps per U, i.e. delta_t = 1/n_steps; one figure "
                         "per value (default: 1, so delta_t = 1)")
    ap.add_argument("--molecules", default=None,
                    help="comma-separated panel order "
                         "(default: the 9 library diatomics)")
    ap.add_argument("--ncol", type=int, default=3)
    ap.add_argument("--no-period-rows", action="store_true",
                    dest="no_period_rows",
                    help="pack panels continuously instead of giving each "
                         "period its own rows. Grouping by period is on by "
                         "default; with 3 columns the nine library diatomics "
                         "split exactly, 3 in period 1 and 6 in period 2")
    ap.add_argument("--sharey", action="store_true",
                    help="share the y axis across panels, so error magnitudes "
                         "are directly comparable between molecules")
    ap.add_argument("--L-min-frac", type=float, default=0.6, dest="L_min_frac",
                    help="the --L-min-frac used when the library was built; "
                         "sets the L/L_eq scale (default: 0.6)")
    ap.add_argument("--cmap", default="viridis")
    args = ap.parse_args()

    df = load(args.source)

    atoms = df[df["L"].isna()]
    df = df[df["L"].notna()]
    if df.empty:
        raise SystemExit("no rows with a bond length; nothing to plot")

    # L / L_eq, comparable across molecules. L_min = L_min_frac * L_eq, so
    # L / L_eq = L_min_frac * L / L_min, taken per molecule. Rounded to 3 dp
    # because the .dat filenames carry L to 2 dp, and that rounding differs per
    # molecule -- at 4 dp the same sweep position would split into separate
    # groups. Sweep positions are 0.267 apart, so 3 dp resolves them cleanly.
    df = df.copy()
    df["L_frac"] = df.groupby("molecule")["L"].transform(
        lambda s: args.L_min_frac * s / s.min()).round(3)

    if args.molecules:
        molecules = [m.strip() for m in args.molecules.split(",")]
    else:
        present = set(df["molecule"].unique())
        molecules = [m for m in DEFAULT_MOLECULES if m in present]
        molecules += sorted(present - set(DEFAULT_MOLECULES))
    if not args.no_period_rows:
        # stable sort: period first, existing order preserved within a period
        molecules.sort(key=species_period)
    if not molecules:
        raise SystemExit("no molecules to plot")

    Ts = ([int(t) for t in args.T_values.split(",")] if args.T_values
          else sorted(int(t) for t in df["T"].unique()))
    orders = [int(o) for o in args.order.split(",")]
    steps = [int(s) for s in args.n_steps.split(",")]

    bases = sorted(df["basis"].unique())
    maps = sorted(df["mapping"].unique())
    if len(bases) > 1 or len(maps) > 1:
        print(f"WARNING: mixed basis {bases} / mapping {maps}; point this at a "
              f"single <basis>/<mapping>/ directory")
    basis, mapping = "+".join(bases), "+".join(maps)

    src = Path(args.source)
    out_dir = Path(args.out_dir) if args.out_dir else (
        src if src.is_dir() else src.parent)
    out_dir.mkdir(parents=True, exist_ok=True)

    norm = Normalize(vmin=df["L_frac"].min(), vmax=df["L_frac"].max())
    cmap = plt.get_cmap(args.cmap)

    n = 0
    for T in Ts:
        for order in orders:
            for ns in steps:
                d = df[(df["T"] == T) & (df["order"] == order)
                       & (df["n_steps"] == ns)]
                if d.empty:
                    print(f"  no rows for T={T} order={order} n_steps={ns}")
                    continue
                fig = one_figure(d, molecules, basis, mapping, T, order, ns,
                                 args, norm, cmap)
                out = out_dir / (f"state_errors_bymol_{basis}_{mapping}"
                                 f"_T{T}_ord{order}_ns{ns}.png")
                fig.savefig(out, dpi=150, bbox_inches="tight")
                plt.close(fig)
                n += 1
                print(f"wrote {out}   ({len(d)} rows, "
                      f"{d['molecule'].nunique()} molecules, "
                      f"{d['L_frac'].nunique()} bond lengths)")

    print(f"\n{n} figures written.")
    if not atoms.empty:
        print(f"({len(atoms)} single-atom rows dropped: no bond length, so the "
              f"L/L_eq colour axis does not apply. Use plot_state_errors.py "
              f"for those.)")


if __name__ == "__main__":
    main()