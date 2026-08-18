#!/usr/bin/env python3
"""Compare four QHAT Trotter orderings across multiple active spaces.

The four canonical methods are:
  1. jw_raw
  2. jw_magnitude_descending_lexicographic   (the JW baseline)
  3. fermionic_signed_coefficient_lexicographic / fermionic_signed_reference
  4. fermionic_coloring

This script is designed to combine rows from different QHAT result CSVs, e.g.
new deterministic-ordering results plus L-sweep coloring results, while requiring
all plotted rows to use the SAME error metric and Trotter settings.

Typical QHAT-root usage:

  python analysis/plot_four_way_active_spaces.py \
      --csv analysis/deterministic_orderings_hgbs5_results.csv \
            analysis/l_sweep_trotter_state_t1.csv \
      --basis hgbs-5 \
      --steps 100 \
      --time 1.0 \
      --formula-order 1 \
      --metric state_infidelity \
      --out-dir analysis/four_way_active_space_figures

You may pass any number of CSV files. Rows that are not one of the four target
methods are ignored.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Canonical names used internally by the plotter.
JW_RAW = "jw_raw"
JW_MAG = "jw_magnitude_descending"
FERM_AWARE = "fermionic_aware"
FERM_COLOR = "fermionic_coloring"

METHOD_ORDER = [JW_RAW, JW_MAG, FERM_AWARE, FERM_COLOR]

# QHAT names seen in the current L-sweep / deterministic-ordering outputs.
ALIASES = {
    # Raw JW
    "jw_raw": JW_RAW,
    "raw_jw": JW_RAW,

    # Descending-magnitude JW baseline
    "desc_magnitude": JW_MAG,
    "jw_magnitude_descending": JW_MAG,
    "jw_magnitude_descending_lexicographic": JW_MAG,
    "jw_magnitude_reference": JW_MAG,

    # Fermionic-aware signed parent ordering
    "fermionic_signed_reference": FERM_AWARE,
    "fermionic_signed": FERM_AWARE,
    "fermionic_signed_coefficient_lexicographic": FERM_AWARE,

    # Fermionic coloring
    "fermionic_coloring": FERM_COLOR,
    "fermionic_color": FERM_COLOR,
}

LABELS = {
    JW_RAW: "JW raw",
    JW_MAG: "JW descending |c| baseline",
    FERM_AWARE: "Fermionic-aware signed",
    FERM_COLOR: "Fermionic coloring",
}

# Match the simple marker/line vocabulary of plot_ordering_study.py.
MARKERS = {
    JW_RAW: "o",
    JW_MAG: "v",
    FERM_AWARE: "s",
    FERM_COLOR: "^",
}

LINESTYLES = {
    JW_RAW: "-",
    JW_MAG: "--",
    FERM_AWARE: "-.",
    FERM_COLOR: ":",
}


METRIC_LABELS = {
    "state_infidelity": "state infidelity",
    "state_vector_2norm_error": r"state-vector error  $\\|\\psi_{Trot}-\\psi_{exact}\\|$",
    "phase_aligned_state_2norm_error": "phase-aligned state-vector error",
    "one_minus_overlap": r"$1-|\\langle\\psi_{exact}|\\psi_{Trot}\\rangle|$",
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--csv",
        nargs="+",
        type=Path,
        required=True,
        help="One or more detailed QHAT result CSVs. Pass deterministic and coloring outputs together.",
    )
    ap.add_argument("--out-dir", type=Path, default=Path("analysis/four_way_active_space_figures"))
    ap.add_argument(
        "--metric",
        default="state_infidelity",
        choices=(
            "state_infidelity",
            "state_vector_2norm_error",
            "phase_aligned_state_2norm_error",
            "one_minus_overlap",
        ),
        help=(
            "Error metric to compare. state_infidelity is the safest default when mixing "
            "older coloring CSVs with the newer deterministic sweep."
        ),
    )
    ap.add_argument("--basis", default="hgbs-5", help="Basis filter; use 'all' to disable.")
    ap.add_argument("--molecules", default=None, help="Comma-separated molecules; default: all complete molecules.")
    ap.add_argument("--steps", type=int, default=100, help="Trotter steps (default: 100).")
    ap.add_argument("--time", type=float, default=1.0, dest="evolution_time", help="Evolution time (default: 1.0).")
    ap.add_argument("--formula-order", type=int, default=1, choices=(1, 2, 4), help="Trotter formula order.")
    ap.add_argument(
        "--bond-length",
        type=float,
        default=None,
        help="Optional fixed bond length. If omitted, all matching bond lengths are shown within each active-space group.",
    )
    ap.add_argument(
        "--complete-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Plot only Hamiltonian cases having all four methods (default: true).",
    )
    ap.add_argument(
        "--ratio-panel",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add error / JW-descending-magnitude ratio panel (default: true).",
    )
    ap.add_argument("--plot-floor", type=float, default=1e-16, help="Positive floor for log plotting only.")
    ap.add_argument("--dpi", type=int, default=180)
    return ap.parse_args()


def _first_existing(df: pd.DataFrame, names: tuple[str, ...]) -> str | None:
    return next((name for name in names if name in df.columns), None)


def _parse_active_text(value: object) -> tuple[float, float]:
    if pd.isna(value):
        return (math.nan, math.nan)
    text = str(value).strip()
    m = re.fullmatch(r"\s*(\d+)\s*[-+_/,:]\s*(\d+)\s*", text)
    if not m:
        return (math.nan, math.nan)
    return (float(m.group(1)), float(m.group(2)))


def normalize_one_csv(path: Path, metric: str, source_index: int) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path, comment="#")
    if df.empty:
        return pd.DataFrame()

    if "status" in df.columns:
        df = df[df["status"].fillna("success").astype(str).str.lower().eq("success")].copy()

    ordering_col = _first_existing(df, ("ordering", "schedule"))
    if ordering_col is None:
        return pd.DataFrame()

    df["method"] = df[ordering_col].astype(str).map(lambda x: ALIASES.get(x.strip()))
    df = df[df["method"].notna()].copy()
    if df.empty:
        return df

    # Strict metric handling: never silently replace one physical error metric by another.
    if metric in df.columns:
        df["error"] = pd.to_numeric(df[metric], errors="coerce")
    elif metric == "one_minus_overlap" and "state_overlap_abs" in df.columns:
        df["error"] = 1.0 - pd.to_numeric(df["state_overlap_abs"], errors="coerce")
    else:
        methods = ", ".join(sorted(df["method"].dropna().unique()))
        raise ValueError(
            f"{path} contains target method(s) [{methods}] but does not provide metric {metric!r}. "
            "Choose a metric available in every input source; state_infidelity is usually the common choice."
        )

    # Standardize metadata names from deterministic, coloring, and older ordering-study outputs.
    molecule_col = _first_existing(df, ("molecule", "mol"))
    if molecule_col is None:
        raise ValueError(f"{path}: missing molecule column")
    df["molecule_std"] = df[molecule_col].astype(str)

    basis_col = _first_existing(df, ("basis",))
    df["basis_std"] = df[basis_col].astype(str) if basis_col else "unknown"

    if "active_occupied" in df.columns and "active_vacant" in df.columns:
        df["active_occupied_std"] = pd.to_numeric(df["active_occupied"], errors="coerce")
        df["active_vacant_std"] = pd.to_numeric(df["active_vacant"], errors="coerce")
    elif "active" in df.columns:
        parsed = df["active"].map(_parse_active_text)
        df["active_occupied_std"] = parsed.map(lambda t: t[0])
        df["active_vacant_std"] = parsed.map(lambda t: t[1])
    else:
        df["active_occupied_std"] = math.nan
        df["active_vacant_std"] = math.nan

    if "n_qubits" in df.columns:
        df["n_qubits_std"] = pd.to_numeric(df["n_qubits"], errors="coerce")
    else:
        df["n_qubits_std"] = df["active_occupied_std"] + df["active_vacant_std"]

    bond_col = _first_existing(df, ("bond_length", "L", "geometry_scale"))
    df["bond_length_std"] = pd.to_numeric(df[bond_col], errors="coerce") if bond_col else math.nan

    case_col = _first_existing(df, ("case_id",))
    if case_col:
        df["case_id_std"] = df[case_col].astype(str)
    else:
        def make_case(row: pd.Series) -> str:
            return (
                f"{row['molecule_std']}|{row['basis_std']}|"
                f"{row['active_occupied_std']}|{row['active_vacant_std']}|"
                f"{row['bond_length_std']}"
            )
        df["case_id_std"] = df.apply(make_case, axis=1)

    steps_col = _first_existing(df, ("trotter_steps",))
    df["steps_std"] = pd.to_numeric(df[steps_col], errors="coerce") if steps_col else math.nan

    time_col = _first_existing(df, ("evolution_time",))
    df["time_std"] = pd.to_numeric(df[time_col], errors="coerce") if time_col else math.nan

    formula_col = _first_existing(df, ("formula_order",))
    # Current deterministic fermionic-aware benchmark is explicitly first order even though
    # its CSV does not carry formula_order.
    df["formula_order_std"] = (
        pd.to_numeric(df[formula_col], errors="coerce") if formula_col else 1
    )

    df["source_file"] = str(path)
    df["source_index"] = source_index

    keep = [
        "case_id_std", "molecule_std", "basis_std", "active_occupied_std",
        "active_vacant_std", "n_qubits_std", "bond_length_std", "steps_std",
        "time_std", "formula_order_std", "method", "error", "source_file",
        "source_index",
    ]
    return df[keep].dropna(subset=["error"]).copy()


def load_data(paths: list[Path], metric: str) -> pd.DataFrame:
    frames = []
    for i, path in enumerate(paths):
        frame = normalize_one_csv(path, metric, i)
        if not frame.empty:
            frames.append(frame)
            print(f"loaded {path}: {len(frame)} target-method rows")
    if not frames:
        raise SystemExit("No rows for the four target methods were found in the supplied CSVs.")
    return pd.concat(frames, ignore_index=True)


def filter_data(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    out = df.copy()

    if args.basis.lower() != "all":
        out = out[out["basis_std"].str.lower().eq(args.basis.lower())]

    if args.molecules:
        wanted = {m.strip() for m in args.molecules.split(",") if m.strip()}
        out = out[out["molecule_std"].isin(wanted)]

    # Only enforce filters where the standardized metadata exists.
    if out["steps_std"].notna().any():
        out = out[np.isclose(out["steps_std"], args.steps)]
    if out["time_std"].notna().any():
        out = out[np.isclose(out["time_std"], args.evolution_time)]
    if out["formula_order_std"].notna().any():
        out = out[np.isclose(out["formula_order_std"], args.formula_order)]
    if args.bond_length is not None and out["bond_length_std"].notna().any():
        out = out[np.isclose(out["bond_length_std"], args.bond_length)]

    return out.copy()


def collapse_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate copies of the same method/case from overlapping CSVs.

    JW raw is often present in both deterministic and coloring files. If copies differ
    materially, fail rather than silently comparing mismatched runs.
    """
    keys = [
        "case_id_std", "molecule_std", "basis_std", "active_occupied_std",
        "active_vacant_std", "n_qubits_std", "bond_length_std", "steps_std",
        "time_std", "formula_order_std", "method",
    ]
    rows = []
    for _, g in df.groupby(keys, dropna=False, sort=False):
        vals = g["error"].to_numpy(dtype=float)
        finite = vals[np.isfinite(vals)]
        if len(finite) == 0:
            continue
        lo, hi = float(np.min(finite)), float(np.max(finite))
        scale = max(abs(lo), abs(hi), 1e-30)
        if len(finite) > 1 and abs(hi - lo) / scale > 1e-6:
            first = g.iloc[0]
            raise ValueError(
                "Conflicting duplicate results for "
                f"{first['case_id_std']} / {first['method']}: {finite.tolist()}. "
                "The CSVs may use different Trotter settings or Hamiltonian revisions."
            )
        # Preserve the first source in command-line order.
        rows.append(g.sort_values("source_index").iloc[0])
    return pd.DataFrame(rows).reset_index(drop=True)


def active_label(row: pd.Series) -> str:
    occ = row["active_occupied_std"]
    vac = row["active_vacant_std"]
    nq = row["n_qubits_std"]
    if pd.notna(occ) and pd.notna(vac):
        return f"{int(occ)}+{int(vac)}\n{int(nq)}q"
    if pd.notna(nq):
        return f"{int(nq)}q"
    return str(row["case_id_std"])


def case_sort_table(df: pd.DataFrame) -> pd.DataFrame:
    meta_cols = [
        "case_id_std", "molecule_std", "basis_std", "active_occupied_std",
        "active_vacant_std", "n_qubits_std", "bond_length_std",
    ]
    meta = df[meta_cols].drop_duplicates("case_id_std").copy()
    meta["_nq"] = meta["n_qubits_std"].fillna(np.inf)
    meta["_occ"] = meta["active_occupied_std"].fillna(np.inf)
    meta["_vac"] = meta["active_vacant_std"].fillna(np.inf)
    meta["_L"] = meta["bond_length_std"].fillna(-np.inf)
    return meta.sort_values(["_nq", "_occ", "_vac", "_L", "case_id_std"])


def plot_one_molecule(dmol: pd.DataFrame, args: argparse.Namespace, out_dir: Path) -> Path | None:
    mol = str(dmol["molecule_std"].iloc[0])
    cases = case_sort_table(dmol)

    # Keep only complete four-method cases if requested.
    if args.complete_only:
        counts = dmol.groupby("case_id_std")["method"].nunique()
        complete_ids = counts[counts == len(METHOD_ORDER)].index
        cases = cases[cases["case_id_std"].isin(complete_ids)]
        dmol = dmol[dmol["case_id_std"].isin(complete_ids)].copy()

    if cases.empty:
        print(f"skip {mol}: no complete four-method cases after filtering")
        return None

    case_ids = cases["case_id_std"].tolist()
    x = np.arange(len(case_ids), dtype=float)
    labels = [active_label(row) for _, row in cases.iterrows()]

    # If one active-space label appears more than once, append bond length to distinguish cases.
    dup_labels = pd.Series(labels).duplicated(keep=False).to_numpy()
    for i, is_dup in enumerate(dup_labels):
        if is_dup:
            L = cases.iloc[i]["bond_length_std"]
            if pd.notna(L):
                labels[i] += f"\nL={L:g}"

    nrows = 2 if args.ratio_panel else 1
    height = 7.0 if args.ratio_panel else 4.7
    fig, axes = plt.subplots(
        nrows, 1,
        figsize=(max(8.5, 0.85 * len(case_ids) + 4.0), height),
        sharex=True,
        squeeze=False,
        gridspec_kw={"height_ratios": [2.2, 1.25]} if args.ratio_panel else None,
    )
    ax = axes[0][0]

    # Pivot once so every method is aligned to the identical Hamiltonian case.
    pivot = dmol.pivot(index="case_id_std", columns="method", values="error").reindex(case_ids)

    for method in METHOD_ORDER:
        if method not in pivot.columns:
            continue
        y = pivot[method].to_numpy(dtype=float)
        yplot = np.maximum(y, args.plot_floor)
        ax.plot(
            x, yplot,
            marker=MARKERS[method],
            linestyle=LINESTYLES[method],
            linewidth=1.6,
            markersize=6.5,
            label=LABELS[method],
        )

    ax.set_yscale("log")
    ax.set_ylabel(METRIC_LABELS.get(args.metric, args.metric))
    ax.grid(True, axis="y", which="major", alpha=0.3, ls="--", lw=0.7)
    ax.set_title(
        f"{mol}: Trotter ordering across active spaces\n"
        f"{args.basis}, T={args.evolution_time:g}, r={args.steps}, first-order"
    )
    ax.legend(loc="best", frameon=True)

    if args.ratio_panel:
        rax = axes[1][0]
        if JW_MAG not in pivot.columns:
            raise ValueError(f"{mol}: missing JW descending-magnitude baseline")
        baseline = pivot[JW_MAG].to_numpy(dtype=float)
        for method in METHOD_ORDER:
            if method not in pivot.columns:
                continue
            values = pivot[method].to_numpy(dtype=float)
            ratio = np.divide(
                values,
                baseline,
                out=np.full_like(values, np.nan, dtype=float),
                where=baseline > 0,
            )
            rax.plot(
                x,
                np.maximum(ratio, args.plot_floor),
                marker=MARKERS[method],
                linestyle=LINESTYLES[method],
                linewidth=1.4,
                markersize=5.5,
                label=LABELS[method],
            )
        rax.axhline(1.0, color="0.35", lw=1.1, ls="--")
        rax.set_yscale("log")
        rax.set_ylabel("error / JW |c| baseline")
        rax.grid(True, axis="y", which="major", alpha=0.3, ls="--", lw=0.7)

    axes[-1][0].set_xticks(x)
    axes[-1][0].set_xticklabels(labels, rotation=0, fontsize=9)
    axes[-1][0].set_xlabel("active occupied + active virtual orbitals (total qubits)")

    fig.tight_layout()
    safe_mol = re.sub(r"[^A-Za-z0-9_.-]+", "_", mol)
    out = out_dir / f"four_way_active_spaces_{safe_mol}.png"
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return out


def write_summary(df: pd.DataFrame, args: argparse.Namespace, out_dir: Path) -> Path:
    keys = [
        "case_id_std", "molecule_std", "basis_std", "active_occupied_std",
        "active_vacant_std", "n_qubits_std", "bond_length_std",
    ]
    wide = df.pivot_table(index=keys, columns="method", values="error", aggfunc="first").reset_index()
    for method in METHOD_ORDER:
        if method not in wide.columns:
            wide[method] = np.nan

    baseline = wide[JW_MAG]
    for method in (JW_RAW, FERM_AWARE, FERM_COLOR):
        wide[f"{method}_ratio_to_jw_magnitude"] = wide[method] / baseline
        wide[f"jw_magnitude_advantage_over_{method}"] = baseline / wide[method]

    wide["all_four_present"] = wide[METHOD_ORDER].notna().all(axis=1)
    out = out_dir / "four_way_active_space_summary.csv"
    wide.to_csv(out, index=False)
    print(f"wrote {out}")
    return out


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = load_data(args.csv, args.metric)
    df = filter_data(df, args)
    if df.empty:
        raise SystemExit("No target rows remain after the requested basis/Trotter filters.")

    df = collapse_duplicates(df)
    write_summary(df, args, args.out_dir)

    molecules = sorted(df["molecule_std"].unique())
    if args.molecules:
        requested = [m.strip() for m in args.molecules.split(",") if m.strip()]
        molecules = [m for m in requested if m in set(molecules)]

    made = 0
    for mol in molecules:
        dmol = df[df["molecule_std"] == mol].copy()
        if plot_one_molecule(dmol, args, args.out_dir) is not None:
            made += 1

    print(f"\n{made} molecule figure(s) written")
    if made == 0:
        print("Tip: use --no-complete-only to inspect which methods are missing per case.")


if __name__ == "__main__":
    main()
