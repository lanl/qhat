#!/usr/bin/env python3
"""Batch QHAT L-sweep Trotter-error analyzer.

This script extracts the full-library case-study logic from
``qhat_L_sweep_trotter_demo.ipynb`` and applies it to every generated JW
Hamiltonian up to a configurable qubit limit.

For the reference L-sweep library, the default settings select 108 JW cases:

    108 cases x 2 decompositions x 2 Trotter orders x 4 step counts
    = 1,728 CSV rows.

Run this script from the QHAT repository root, for example:

    python run_all_L_sweep_cases.py \
        --library demo_L_sweep_library \
        --output qhat_full_library_coloring_trotter_errors.csv

Dependencies: numpy, pandas, networkx, scipy.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import networkx as nx
import numpy as np
import pandas as pd
from scipy.linalg import expm, norm


# -----------------------------------------------------------------------------
# Dense one-qubit matrices
# -----------------------------------------------------------------------------

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI_MATRIX = {"I": I2, "X": X, "Y": Y, "Z": Z}

NAME_RE = re.compile(r"_as-(\d{3})-(\d{3})_(jw|bk)\.dat$", re.IGNORECASE)

RESULT_COLUMNS = [
    "molecule",
    "L",
    "basis",
    "active_occupied",
    "active_vacant",
    "qubits",
    "scheme",
    "vertices",
    "noncommuting_edges",
    "colors",
    "order",
    "steps",
    "spectral_error",
    "reconstruction_error",
    "file",
]

SCHEME_ORDER = {
    "fermionic_commutation_then_JW": 0,
    "JW_Pauli_string_commutation": 1,
}


@dataclass(frozen=True)
class Case:
    molecule: str
    bond_length: float
    basis: str
    active_occupied: int
    active_vacant: int
    mapping: str
    qubits: int
    path: Path

    @property
    def key(self) -> tuple[str, float, str, int, int]:
        return (
            self.molecule,
            self.bond_length,
            self.basis,
            self.active_occupied,
            self.active_vacant,
        )


# -----------------------------------------------------------------------------
# File inventory and parsing
# -----------------------------------------------------------------------------


def build_inventory(library_root: Path) -> list[Case]:
    """Find every generated QHAT ``_jw.dat`` and ``_bk.dat`` file."""
    library_root = library_root.resolve()
    if not library_root.is_dir():
        raise FileNotFoundError(f"Sweep library does not exist: {library_root}")

    cases: list[Case] = []
    for path in sorted(library_root.rglob("*.dat")):
        match = NAME_RE.search(path.name)
        if match is None:
            continue

        relative = path.relative_to(library_root)
        if len(relative.parts) < 4:
            print(f"WARNING: unexpected library path, skipping: {path}", file=sys.stderr)
            continue

        occupied = int(match.group(1))
        vacant = int(match.group(2))
        try:
            bond_length = float(relative.parts[1])
        except ValueError:
            print(f"WARNING: invalid bond-length directory, skipping: {path}", file=sys.stderr)
            continue

        cases.append(
            Case(
                molecule=relative.parts[0],
                bond_length=bond_length,
                basis=relative.parts[2],
                active_occupied=occupied,
                active_vacant=vacant,
                mapping=match.group(3).lower(),
                qubits=occupied + vacant,
                path=path.resolve(),
            )
        )

    return cases


def tensor_path_for_dat(path: Path) -> Path:
    """Return the tensor file paired with a mapped QHAT ``.dat`` file."""
    replaced = re.sub(
        r"_(?:jw|bk)\.dat$",
        ".tensors.npz",
        str(path),
        flags=re.IGNORECASE,
    )
    return Path(replaced)


def parse_qhat_dat(path: Path) -> tuple[dict[str, str], list[tuple[str, complex]]]:
    """Parse QHAT's default coefficient/Pauli-word ``.dat`` format."""
    metadata: dict[str, str] = {}
    terms: list[tuple[str, complex]] = []

    with path.open("r", encoding="ascii") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith("#"):
                content = line[1:].strip()
                if "=" in content:
                    key, value = content.split("=", maxsplit=1)
                    metadata[key.strip()] = value.strip()
                continue

            fields = line.split()
            if len(fields) != 2:
                raise ValueError(
                    f"Expected '<coefficient> <Pauli word>' in {path}:{line_number}; "
                    f"found {line!r}"
                )

            coefficient = complex(fields[0])
            pauli_word = fields[1].upper()
            invalid = set(pauli_word) - set(PAULI_MATRIX)
            if invalid:
                raise ValueError(
                    f"Invalid Pauli symbols {sorted(invalid)} in {path}:{line_number}"
                )
            terms.append((pauli_word, coefficient))

    if not terms:
        raise ValueError(f"No Pauli terms found in {path}")

    lengths = {len(word) for word, _ in terms}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent Pauli-word lengths in {path}: {sorted(lengths)}")

    return metadata, terms


# -----------------------------------------------------------------------------
# Matrix construction
# -----------------------------------------------------------------------------


def kron_all(items: Iterable[np.ndarray]) -> np.ndarray:
    output = np.array([[1]], dtype=complex)
    for item in items:
        output = np.kron(output, item)
    return output


def pauli_matrix(pauli_word: str) -> np.ndarray:
    return kron_all(PAULI_MATRIX[symbol] for symbol in pauli_word)


def ladder_matrices(n_qubits: int) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Construct dense Jordan-Wigner annihilation and creation matrices."""
    annihilation: list[np.ndarray] = []
    creation: list[np.ndarray] = []

    for mode in range(n_qubits):
        factors_a: list[np.ndarray] = []
        factors_c: list[np.ndarray] = []

        for qubit in range(n_qubits):
            if qubit < mode:
                factors_a.append(Z)
                factors_c.append(Z)
            elif qubit == mode:
                factors_a.append((X + 1j * Y) / 2)
                factors_c.append((X - 1j * Y) / 2)
            else:
                factors_a.append(I2)
                factors_c.append(I2)

        annihilation.append(kron_all(factors_a))
        creation.append(kron_all(factors_c))

    return annihilation, creation


def monomial_matrix(
    key: tuple[tuple[int, int], ...],
    annihilation: Sequence[np.ndarray],
    creation: Sequence[np.ndarray],
) -> np.ndarray:
    dim = annihilation[0].shape[0]
    output = np.eye(dim, dtype=complex)
    for mode, action in key:
        output = output @ (creation[mode] if action else annihilation[mode])
    return output


def dagger_key(key: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    return tuple((mode, 1 - action) for mode, action in reversed(key))


def tensor_monomials(
    tensor_file: Path,
    tolerance: float = 1.0e-12,
) -> tuple[dict[tuple[tuple[int, int], ...], complex], int]:
    """Read the QHAT InteractionOperator tensors as fermionic monomials."""
    with np.load(tensor_file) as data:
        constant = (
            complex(np.asarray(data["constant"]).item())
            if "constant" in data.files
            else 0.0 + 0.0j
        )
        one_body = np.asarray(data["one_body"])
        two_body = np.asarray(data["two_body"])

    if one_body.ndim != 2 or one_body.shape[0] != one_body.shape[1]:
        raise ValueError(f"Invalid one_body shape in {tensor_file}: {one_body.shape}")

    n_qubits = int(one_body.shape[0])
    if two_body.shape != (n_qubits, n_qubits, n_qubits, n_qubits):
        raise ValueError(f"Invalid two_body shape in {tensor_file}: {two_body.shape}")

    coefficients: dict[tuple[tuple[int, int], ...], complex] = {(): constant}

    for p, q in np.argwhere(np.abs(one_body) > tolerance):
        key = ((int(p), 1), (int(q), 0))
        coefficients[key] = complex(one_body[p, q])

    for p, q, r, s in np.argwhere(np.abs(two_body) > tolerance):
        key = (
            (int(p), 1),
            (int(q), 1),
            (int(r), 0),
            (int(s), 0),
        )
        coefficients[key] = complex(two_body[p, q, r, s])

    return coefficients, n_qubits


def hermitian_vertices_from_tensors(
    tensor_file: Path,
    tolerance: float = 1.0e-12,
) -> tuple[list[np.ndarray], list[tuple[tuple, tuple]], int]:
    """Group each monomial with its Hermitian conjugate and map it to JW."""
    coefficients, n_qubits = tensor_monomials(tensor_file, tolerance)
    annihilation, creation = ladder_matrices(n_qubits)

    used: set[tuple[tuple[int, int], ...]] = set()
    vertices: list[np.ndarray] = []
    labels: list[tuple[tuple, tuple]] = []

    for key, coefficient in coefficients.items():
        if key in used:
            continue

        partner = dagger_key(key)
        matrix = (
            monomial_matrix(key, annihilation, creation)
            if key
            else np.eye(2**n_qubits, dtype=complex)
        )

        if partner == key:
            vertex = coefficient * matrix
            used.add(key)
        else:
            partner_coefficient = coefficients.get(partner, 0.0 + 0.0j)
            partner_matrix = monomial_matrix(partner, annihilation, creation)
            vertex = coefficient * matrix + partner_coefficient * partner_matrix
            used.update((key, partner))

        if norm(vertex, 2) > tolerance:
            vertices.append(vertex)
            labels.append((key, partner))

    if not vertices:
        raise ValueError(f"No nonzero fermionic vertices reconstructed from {tensor_file}")

    return vertices, labels, n_qubits


# -----------------------------------------------------------------------------
# Noncommutation graph and product formulas
# -----------------------------------------------------------------------------


def matrix_noncommutation_graph(
    vertices: Sequence[np.ndarray],
    tolerance: float = 1.0e-10,
) -> tuple[nx.Graph, dict[int, int], dict[int, list[int]]]:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(vertices)))

    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            commutator = vertices[i] @ vertices[j] - vertices[j] @ vertices[i]
            if norm(commutator, 2) > tolerance:
                graph.add_edge(i, j)

    colors = nx.coloring.greedy_color(graph, strategy="largest_first")
    color_groups: dict[int, list[int]] = {}
    for node, color in colors.items():
        color_groups.setdefault(color, []).append(node)

    return graph, colors, color_groups


def product_formula_from_groups(
    vertices: Sequence[np.ndarray],
    color_groups: dict[int, list[int]],
    total_time: float,
    nsteps: int,
    order: int,
) -> np.ndarray:
    if order not in (1, 2):
        raise ValueError(f"Only first- and second-order formulas are supported; got {order}")
    if nsteps <= 0:
        raise ValueError(f"nsteps must be positive; got {nsteps}")

    dim = vertices[0].shape[0]
    identity = np.eye(dim, dtype=complex)
    delta_t = total_time / nsteps
    ordered_groups = [color_groups[color] for color in sorted(color_groups)]

    if order == 1:
        schedule = [(group, delta_t) for group in ordered_groups]
    else:
        schedule = [
            *((group, delta_t / 2) for group in ordered_groups),
            *((group, delta_t / 2) for group in reversed(ordered_groups)),
        ]

    one_step = identity.copy()
    for group, tau in schedule:
        for index in group:
            one_step = expm(-1j * tau * vertices[index]) @ one_step

    return np.linalg.matrix_power(one_step, nsteps)


# -----------------------------------------------------------------------------
# Per-case and batch analysis
# -----------------------------------------------------------------------------


def analyze_coloring_case(
    dat_path: Path,
    total_time: float = 1.0,
    orders: Sequence[int] = (1, 2),
    steps_list: Sequence[int] = (1, 2, 5, 10),
    reconstruction_tolerance: float = 1.0e-8,
) -> list[dict[str, object]]:
    dat_path = Path(dat_path)
    tensor_path = tensor_path_for_dat(dat_path)
    if not tensor_path.exists():
        raise FileNotFoundError(f"Missing tensor file: {tensor_path}")

    _, mapped_terms = parse_qhat_dat(dat_path)
    fermionic_vertices, _, n_qubits = hermitian_vertices_from_tensors(tensor_path)

    pauli_word_lengths = {len(word) for word, _ in mapped_terms}
    if pauli_word_lengths != {n_qubits}:
        raise ValueError(
            f"Tensor file has {n_qubits} qubits but {dat_path} has Pauli lengths "
            f"{sorted(pauli_word_lengths)}"
        )

    fermionic_graph, _, fermionic_groups = matrix_noncommutation_graph(
        fermionic_vertices
    )

    pauli_vertices = [coefficient * pauli_matrix(word) for word, coefficient in mapped_terms]
    pauli_graph, _, pauli_groups = matrix_noncommutation_graph(pauli_vertices)

    zero = np.zeros_like(fermionic_vertices[0])
    hamiltonian_from_tensors = sum(fermionic_vertices, zero.copy())
    hamiltonian_from_dat = sum(pauli_vertices, zero.copy())
    reconstruction_error = float(norm(hamiltonian_from_tensors - hamiltonian_from_dat, 2))

    if reconstruction_error > reconstruction_tolerance:
        raise ValueError(
            f"Tensor/JW reconstruction error {reconstruction_error:.3e} exceeds "
            f"{reconstruction_tolerance:.3e}"
        )

    exact_evolution = expm(-1j * total_time * hamiltonian_from_dat)

    output: list[dict[str, object]] = []
    decompositions = (
        (
            "fermionic_commutation_then_JW",
            fermionic_vertices,
            fermionic_groups,
            fermionic_graph,
        ),
        (
            "JW_Pauli_string_commutation",
            pauli_vertices,
            pauli_groups,
            pauli_graph,
        ),
    )

    for scheme, vertices, groups, graph in decompositions:
        for order in orders:
            for steps in steps_list:
                approximate_evolution = product_formula_from_groups(
                    vertices=vertices,
                    color_groups=groups,
                    total_time=total_time,
                    nsteps=steps,
                    order=order,
                )
                output.append(
                    {
                        "scheme": scheme,
                        "vertices": len(vertices),
                        "noncommuting_edges": graph.number_of_edges(),
                        "colors": len(groups),
                        "order": int(order),
                        "steps": int(steps),
                        "spectral_error": float(norm(approximate_evolution - exact_evolution, 2)),
                        "reconstruction_error": reconstruction_error,
                    }
                )

    return output


def sort_results(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.reindex(columns=RESULT_COLUMNS)

    output = frame.copy()
    output["_scheme_order"] = output["scheme"].map(SCHEME_ORDER).fillna(99)
    output = output.sort_values(
        [
            "molecule",
            "L",
            "basis",
            "active_occupied",
            "active_vacant",
            "_scheme_order",
            "order",
            "steps",
        ],
        kind="stable",
    ).drop(columns="_scheme_order")
    return output.reindex(columns=RESULT_COLUMNS).reset_index(drop=True)


def completed_case_keys(
    frame: pd.DataFrame,
    rows_per_case: int,
) -> set[tuple[str, float, str, int, int]]:
    if frame.empty:
        return set()

    key_columns = [
        "molecule",
        "L",
        "basis",
        "active_occupied",
        "active_vacant",
    ]
    counts = frame.groupby(key_columns, dropna=False).size()
    return {
        (str(key[0]), float(key[1]), str(key[2]), int(key[3]), int(key[4]))
        for key, count in counts.items()
        if int(count) == rows_per_case
    }


def save_results(rows: list[dict[str, object]], output_path: Path) -> pd.DataFrame:
    frame = sort_results(pd.DataFrame(rows, columns=RESULT_COLUMNS))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)
    return frame


def run_batch(args: argparse.Namespace) -> int:
    library = args.library.resolve()
    output_path = args.output.resolve()
    skipped_path = args.skipped_output.resolve()

    inventory = build_inventory(library)
    jw_cases = [
        case
        for case in inventory
        if case.mapping == "jw" and case.qubits <= args.max_qubits
    ]

    print(f"Library: {library}")
    print(f"Mapped .dat files found: {len(inventory)}")
    print(f"Eligible JW cases at <= {args.max_qubits} qubits: {len(jw_cases)}")

    if args.expect_cases is not None and len(jw_cases) != args.expect_cases:
        print(
            f"WARNING: expected {args.expect_cases} eligible cases but found {len(jw_cases)}. "
            "The script will process the cases that are present.",
            file=sys.stderr,
        )

    rows_per_case = 2 * len(args.orders) * len(args.steps)
    expected_rows = len(jw_cases) * rows_per_case
    print(
        "Expected result rows: "
        f"{len(jw_cases)} cases x 2 schemes x {len(args.orders)} orders "
        f"x {len(args.steps)} step counts = {expected_rows}"
    )

    existing = pd.DataFrame(columns=RESULT_COLUMNS)
    if args.resume and output_path.exists():
        existing = pd.read_csv(output_path)
        missing_columns = set(RESULT_COLUMNS) - set(existing.columns)
        if missing_columns:
            raise ValueError(
                f"Cannot resume from {output_path}; missing columns: {sorted(missing_columns)}"
            )
        existing = existing.reindex(columns=RESULT_COLUMNS)
        print(f"Loaded {len(existing)} existing rows from {output_path}")

    done_keys = completed_case_keys(existing, rows_per_case)
    result_rows = existing.to_dict(orient="records")
    skipped_rows: list[dict[str, object]] = []

    batch_start = time.perf_counter()
    for case_number, case in enumerate(jw_cases, start=1):
        if case.key in done_keys:
            print(f"[{case_number}/{len(jw_cases)}] SKIP completed {case.key}")
            continue

        # Remove partial rows for this case before recomputing it.
        result_rows = [
            row
            for row in result_rows
            if not (
                str(row["molecule"]) == case.molecule
                and float(row["L"]) == case.bond_length
                and str(row["basis"]) == case.basis
                and int(row["active_occupied"]) == case.active_occupied
                and int(row["active_vacant"]) == case.active_vacant
            )
        ]

        case_start = time.perf_counter()
        print(f"[{case_number}/{len(jw_cases)}] RUN {case.key}")

        try:
            measurements = analyze_coloring_case(
                dat_path=case.path,
                total_time=args.total_time,
                orders=args.orders,
                steps_list=args.steps,
                reconstruction_tolerance=args.reconstruction_tolerance,
            )
            for measurement in measurements:
                result_rows.append(
                    {
                        "molecule": case.molecule,
                        "L": case.bond_length,
                        "basis": case.basis,
                        "active_occupied": case.active_occupied,
                        "active_vacant": case.active_vacant,
                        "qubits": case.qubits,
                        **measurement,
                        "file": str(case.path),
                    }
                )

            save_results(result_rows, output_path)
            elapsed = time.perf_counter() - case_start
            print(
                f"    completed {len(measurements)} rows in {elapsed:.2f} s; "
                f"checkpoint saved"
            )
        except Exception as exc:  # Continue the suite and record the failure.
            elapsed = time.perf_counter() - case_start
            skipped_rows.append(
                {
                    "molecule": case.molecule,
                    "L": case.bond_length,
                    "basis": case.basis,
                    "active_occupied": case.active_occupied,
                    "active_vacant": case.active_vacant,
                    "qubits": case.qubits,
                    "file": str(case.path),
                    "reason": f"{type(exc).__name__}: {exc}",
                    "elapsed_seconds": elapsed,
                }
            )
            pd.DataFrame(skipped_rows).to_csv(skipped_path, index=False)
            print(f"    FAILED after {elapsed:.2f} s: {type(exc).__name__}: {exc}")

    final_results = save_results(result_rows, output_path)
    if skipped_rows:
        pd.DataFrame(skipped_rows).to_csv(skipped_path, index=False)
    elif skipped_path.exists():
        skipped_path.unlink()

    total_elapsed = time.perf_counter() - batch_start
    unique_cases = (
        final_results[
            ["molecule", "L", "basis", "active_occupied", "active_vacant"]
        ]
        .drop_duplicates()
        .shape[0]
        if not final_results.empty
        else 0
    )

    print()
    print(f"Saved results: {output_path}")
    print(f"Rows written: {len(final_results)}")
    print(f"Completed unique cases: {unique_cases}/{len(jw_cases)}")
    print(f"Skipped cases in this run: {len(skipped_rows)}")
    if skipped_rows:
        print(f"Skipped-case log: {skipped_path}")
    print(f"Batch runtime: {total_elapsed:.2f} s")

    if len(final_results) != expected_rows:
        print(
            f"WARNING: expected {expected_rows} rows but the output contains "
            f"{len(final_results)} rows.",
            file=sys.stderr,
        )
        return 1

    return 0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the fermionic-coloring versus JW-Pauli-coloring dense Trotter "
            "comparison over a generated QHAT L-sweep library."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("demo_L_sweep_library"),
        help="Root of the generated L-sweep library (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("qhat_full_library_coloring_trotter_errors.csv"),
        help="Main result CSV (default: %(default)s)",
    )
    parser.add_argument(
        "--skipped-output",
        type=Path,
        default=Path("qhat_full_library_coloring_trotter_skipped.csv"),
        help="Failure log CSV (default: %(default)s)",
    )
    parser.add_argument(
        "--max-qubits",
        type=int,
        default=6,
        help="Largest dense exact-evolution case (default: %(default)s)",
    )
    parser.add_argument(
        "--orders",
        type=int,
        nargs="+",
        default=[1, 2],
        choices=[1, 2],
        help="Product-formula orders (default: 1 2)",
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=[1, 2, 5, 10],
        help="Trotter step counts (default: 1 2 5 10)",
    )
    parser.add_argument(
        "--total-time",
        type=float,
        default=1.0,
        help="Total evolution time (default: %(default)s)",
    )
    parser.add_argument(
        "--reconstruction-tolerance",
        type=float,
        default=1.0e-8,
        help="Maximum tensor-to-JW matrix discrepancy (default: %(default)s)",
    )
    parser.add_argument(
        "--expect-cases",
        type=int,
        default=108,
        help="Warn unless this many eligible cases are found; use 0 to disable",
    )
    parser.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Resume from a partially written result CSV (default: enabled)",
    )

    args = parser.parse_args(argv)
    if args.max_qubits < 1:
        parser.error("--max-qubits must be positive")
    if any(step < 1 for step in args.steps):
        parser.error("every --steps value must be positive")
    if args.total_time <= 0:
        parser.error("--total-time must be positive")
    if args.expect_cases == 0:
        args.expect_cases = None
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        return run_batch(args)
    except KeyboardInterrupt:
        print("Interrupted. The last completed case is already checkpointed.", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
