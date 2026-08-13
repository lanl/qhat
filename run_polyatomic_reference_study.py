#!/usr/bin/env python3
"""Run the deterministic polyatomic reference-geometry Trotter study.

Designed for the LANL/QHAT ``L-sweep`` branch and intended to be run from the
repository root.  It uses only existing repository scripts and does not modify
QHAT source code.

Study definition
----------------
Molecules:
    BeH2, H2O, NH3, CH4
Bases:
    sto-6g, hgbs-5
Geometry:
    reference geometry only (scale s=1.00)
Requested active-space sizes:
    4, 6, 8, 10, 12, 14, 16, 18, 20 qubits

Only active spaces supported by each molecule/basis are benchmarked.
The five deterministic orderings are:
    jw_raw
    signed_coefficient_lexicographic
    jw_magnitude_descending_lexicographic
    fermionic_signed_coefficient_lexicographic
    fermionic_magnitude_descending_lexicographic

Modes
-----
configs:
    Write/rewrite the deterministic config files only.
tensors:
    Generate configs, then generate only missing tensor NPZ files.
dry-run:
    Generate configs/tensors, discover actually available active spaces, then
    ask the missing-only benchmark driver to print what would still be run.
full:
    Generate configs/tensors and run every missing deterministic ordering row.
status:
    Do not generate or benchmark anything; just print currently available
    reference-geometry tensor cases and completion counts in the output CSV.

Restart behavior
----------------
- Config generation is deterministic and inexpensive.
- build_L_sweep_tensors.py skips an existing tensor unless --overwrite is used.
- append_missing_deterministic_orderings.py preserves existing successful rows
  and runs only missing ordering rows for the same Trotter settings.
Therefore rerunning ``full`` after an interruption is the intended recovery.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

MOLECULES = ("BeH2", "H2O", "NH3", "CH4")
BASES = ("sto-6g", "hgbs-5")
REQUESTED_QUBITS = tuple(range(4, 21, 2))
SCALE = 1.0
SCALE_TAG = "s-1.00"
STEPS = 100
EVOLUTION_TIME = 1.0

ORDERINGS = (
    "jw_raw",
    "signed_coefficient_lexicographic",
    "jw_magnitude_descending_lexicographic",
    "fermionic_signed_coefficient_lexicographic",
    "fermionic_magnitude_descending_lexicographic",
)

TENSOR_NAME_PATTERN = re.compile(
    r"_as-(?P<occupied>\d{3})-(?P<vacant>\d{3})(?:_JW)?\.tensors\.npz$",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run all uncompleted reference-geometry polyatomic deterministic-ordering cases."
    )
    parser.add_argument(
        "mode",
        nargs="?",
        choices=("configs", "tensors", "dry-run", "full", "status"),
        default="full",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="QHAT repository root (default: current working directory).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/polyatomic_reference_deterministic_ordering.csv"),
        help="Combined deterministic-ordering CSV.",
    )
    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required repository file not found: {path}")


def run_command(command: list[str], cwd: Path) -> None:
    printable = " ".join(command)
    print("\n$ " + printable, flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def relpath(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def generate_configs(repo_root: Path, generator: Path, library: Path) -> None:
    run_command(
        [
            sys.executable,
            relpath(generator, repo_root),
            "--molecules",
            ",".join(MOLECULES),
            "--bases",
            ",".join(BASES),
            "--scale-values",
            f"{SCALE:.1f}",
            "--active-sizes",
            ",".join(str(q) for q in REQUESTED_QUBITS),
            "--mappings",
            "JW",
            "--library",
            relpath(library, repo_root),
        ],
        cwd=repo_root,
    )


def generate_missing_tensors(
    repo_root: Path,
    tensor_builder: Path,
    library: Path,
    summary_dir: Path,
) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)
    for molecule in MOLECULES:
        subset = library / molecule / SCALE_TAG
        if not subset.exists():
            print(f"\nSKIP tensor generation for {molecule}: {subset} does not exist")
            continue
        summary = summary_dir / f"tensor_generation_{molecule}_{SCALE_TAG}.csv"
        run_command(
            [
                sys.executable,
                relpath(tensor_builder, repo_root),
                "--library",
                relpath(subset, repo_root),
                "--summary",
                relpath(summary, repo_root),
            ],
            cwd=repo_root,
        )


def parse_active_space_from_tensor(path: Path) -> tuple[int, int, int] | None:
    match = TENSOR_NAME_PATTERN.search(path.name)
    if match is None:
        return None
    occupied = int(match.group("occupied"))
    vacant = int(match.group("vacant"))
    return occupied, vacant, occupied + vacant


def discover_available_cases(library: Path) -> dict[tuple[str, str], list[Path]]:
    cases: dict[tuple[str, str], list[Path]] = {}
    for molecule in MOLECULES:
        for basis in BASES:
            directory = library / molecule / SCALE_TAG / basis
            paths: list[Path] = []
            if directory.exists():
                for path in directory.glob("*.tensors.npz"):
                    parsed = parse_active_space_from_tensor(path)
                    if parsed is None:
                        continue
                    _, _, n_qubits = parsed
                    if n_qubits in REQUESTED_QUBITS:
                        paths.append(path)
            paths.sort(key=lambda p: (parse_active_space_from_tensor(p)[2], p.name))
            cases[(molecule, basis)] = paths
    return cases


def validate_unique_sizes(cases: dict[tuple[str, str], list[Path]]) -> None:
    problems: list[str] = []
    for (molecule, basis), paths in cases.items():
        by_size: dict[int, list[Path]] = defaultdict(list)
        for path in paths:
            parsed = parse_active_space_from_tensor(path)
            assert parsed is not None
            by_size[parsed[2]].append(path)
        for n_qubits, matches in sorted(by_size.items()):
            if len(matches) > 1:
                problems.append(
                    f"{molecule}/{basis}/{n_qubits}q has {len(matches)} tensors: "
                    + ", ".join(str(path) for path in matches)
                )
    if problems:
        raise RuntimeError(
            "Ambiguous tensor cases were found. Refusing to benchmark:\n  "
            + "\n  ".join(problems)
        )


def print_inventory(cases: dict[tuple[str, str], list[Path]]) -> None:
    print("\n" + "=" * 88)
    print("Reference-geometry polyatomic tensor inventory")
    print("=" * 88)
    total = 0
    for molecule in MOLECULES:
        for basis in BASES:
            paths = cases[(molecule, basis)]
            sizes = [parse_active_space_from_tensor(path)[2] for path in paths]
            total += len(paths)
            text = " ".join(str(size) for size in sizes) if sizes else "(none)"
            print(f"{molecule:5s}  {basis:7s}  qubits: {text}")
    print(f"Total available molecule/basis/active-space cases: {total}")


def read_completed_rows(output: Path) -> dict[tuple[str, str, str, int], set[str]]:
    completed: dict[tuple[str, str, str, int], set[str]] = defaultdict(set)
    if not output.exists() or output.stat().st_size == 0:
        return completed
    with output.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") != "success":
                continue
            try:
                n_qubits = int(row.get("n_qubits", ""))
            except ValueError:
                continue
            try:
                steps = int(float(row.get("trotter_steps", "")))
                evolution_time = float(row.get("evolution_time", ""))
            except ValueError:
                continue
            if steps != STEPS or abs(evolution_time - EVOLUTION_TIME) > 1.0e-12:
                continue
            key = (
                str(row.get("molecule", "")),
                str(row.get("basis", "")).lower(),
                str(row.get("bond_length", "")),
                n_qubits,
            )
            completed[key].add(str(row.get("ordering", "")))
    return completed


def print_completion_status(
    cases: dict[tuple[str, str], list[Path]],
    output: Path,
) -> None:
    completed = read_completed_rows(output)
    total_cases = 0
    complete_cases = 0
    missing_rows = 0
    print("\n" + "=" * 88)
    print("Deterministic-ordering completion status")
    print("=" * 88)
    for molecule in MOLECULES:
        for basis in BASES:
            for path in cases[(molecule, basis)]:
                parsed = parse_active_space_from_tensor(path)
                assert parsed is not None
                n_qubits = parsed[2]
                total_cases += 1
                key = (molecule, basis.lower(), SCALE_TAG, n_qubits)
                present = completed.get(key, set())
                missing = [name for name in ORDERINGS if name not in present]
                if not missing:
                    complete_cases += 1
                    state = "COMPLETE"
                else:
                    missing_rows += len(missing)
                    state = "missing " + ", ".join(missing)
                print(f"{molecule:5s} {basis:7s} {n_qubits:2d}q: {state}")
    print(
        f"Complete cases: {complete_cases}/{total_cases}; "
        f"missing ordering rows: {missing_rows}"
    )
    print(f"Output CSV: {output}")


def qubit_sizes(paths: Iterable[Path]) -> list[int]:
    sizes: list[int] = []
    for path in paths:
        parsed = parse_active_space_from_tensor(path)
        if parsed is not None:
            sizes.append(parsed[2])
    return sorted(set(sizes))


def benchmark_missing(
    repo_root: Path,
    benchmark_driver: Path,
    library: Path,
    output: Path,
    cases: dict[tuple[str, str], list[Path]],
    dry_run: bool,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    for molecule in MOLECULES:
        for basis in BASES:
            sizes = qubit_sizes(cases[(molecule, basis)])
            if not sizes:
                print(f"\nSKIP benchmark {molecule}/{basis}: no tensors available")
                continue
            command = [
                sys.executable,
                relpath(benchmark_driver, repo_root),
                "--library",
                relpath(library, repo_root),
                "--molecules",
                molecule,
                "--basis",
                basis,
                "--qubits",
                *[str(size) for size in sizes],
                "--bond-length",
                f"{molecule}={SCALE_TAG}",
                "--orderings",
                *ORDERINGS,
                "--steps",
                str(STEPS),
                "--time",
                str(EVOLUTION_TIME),
                "--output",
                relpath(output, repo_root),
                # We explicitly run jw_raw if it is missing, so unrelated
                # baseline CSVs are not needed for this controlled study.
                "--reference-glob",
                "",
            ]
            if dry_run:
                command.append("--dry-run")
            run_command(command, cwd=repo_root)


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()

    generator = repo_root / "hamiltonian_generator/build_config_polyatomic_sweep.py"
    tensor_builder = repo_root / "hamiltonian_generator/build_L_sweep_tensors.py"
    benchmark_driver = repo_root / "analysis/append_missing_deterministic_orderings.py"
    library = repo_root / "hamiltonian_generator/polyatomic_library"
    summary_dir = repo_root / "analysis/polyatomic_reference_study"
    output = args.output
    if not output.is_absolute():
        output = repo_root / output

    require_file(generator)
    require_file(tensor_builder)
    require_file(benchmark_driver)

    print("=" * 88)
    print("Polyatomic deterministic-ordering reference-geometry study")
    print("=" * 88)
    print(f"Mode:                {args.mode}")
    print(f"Repository:          {repo_root}")
    print(f"Molecules:           {' '.join(MOLECULES)}")
    print(f"Bases:               {' '.join(BASES)}")
    print(f"Geometry scale:      {SCALE_TAG}")
    print(f"Requested qubits:    {' '.join(str(q) for q in REQUESTED_QUBITS)}")
    print(f"Trotter:             first order, t={EVOLUTION_TIME}, steps={STEPS}")
    print(f"Combined output:     {output}")

    if args.mode == "status":
        cases = discover_available_cases(library)
        validate_unique_sizes(cases)
        print_inventory(cases)
        print_completion_status(cases, output)
        return

    generate_configs(repo_root, generator, library)
    if args.mode == "configs":
        return

    generate_missing_tensors(repo_root, tensor_builder, library, summary_dir)
    cases = discover_available_cases(library)
    validate_unique_sizes(cases)
    print_inventory(cases)

    if args.mode == "tensors":
        print_completion_status(cases, output)
        return

    benchmark_missing(
        repo_root=repo_root,
        benchmark_driver=benchmark_driver,
        library=library,
        output=output,
        cases=cases,
        dry_run=(args.mode == "dry-run"),
    )
    print_completion_status(cases, output)


if __name__ == "__main__":
    main()

