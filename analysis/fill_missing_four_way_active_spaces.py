#!/usr/bin/env python3
"""Audit and fill missing four-way HGBS-5 active-space comparisons.

The target universe combines successful HGBS-5 deterministic active-space
rows with the H-H and He-He bond/active-space cases in the original coloring
sweep.  For each exact tensor/settings key, this script collects:

1. ``jw_raw``
2. ``jw_magnitude_descending_lexicographic``
3. ``fermionic_signed_coefficient_lexicographic``
4. ``fermionic_coloring``

Existing results are reused from every detailed CSV under ``analysis/``.  With
``--execute``, absent tensors are regenerated from their tracked QHAT configs
and evaluates only the missing rows.  Coloring/raw-JW rows use the matrix-free
fixed-particle/fixed-Sz benchmark; deterministic coefficient orderings use the
matching deterministic benchmark implementation.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "qhat-matplotlib-cache"),
)

try:
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as matrix_free
except ImportError:
    import benchmark_b2_signed_coefficient_baseline as baseline
    import benchmark_b2_active_spaces_matrix_free as matrix_free


REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "analysis"

TARGET_CSVS = (
    ANALYSIS_DIR / "deterministic_ordering_additions_hgbs5.csv",
    ANALYSIS_DIR / "heteronuclear_hgbs-5_deterministic_ordering.csv",
    ANALYSIS_DIR / "polyatomic_reference_deterministic_ordering.csv",
)
LIGHT_DIATOMIC_TARGET_CSV = ANALYSIS_DIR / "l_sweep_trotter_state_t1.csv"
LIGHT_DIATOMIC_MOLECULES = {"H-H", "He-He"}

JW_RAW = "jw_raw"
JW_MAGNITUDE = "jw_magnitude_descending_lexicographic"
FERMIONIC_AWARE = "fermionic_signed_coefficient_lexicographic"
FERMIONIC_COLORING = "fermionic_coloring"
METHODS = (JW_RAW, JW_MAGNITUDE, FERMIONIC_AWARE, FERMIONIC_COLORING)

ALIASES = {
    "jw_raw": JW_RAW,
    "raw_jw": JW_RAW,
    "jw_magnitude_descending_lexicographic": JW_MAGNITUDE,
    "jw_magnitude_descending": JW_MAGNITUDE,
    "jw_magnitude_reference": JW_MAGNITUDE,
    "desc_magnitude": JW_MAGNITUDE,
    "fermionic_signed_coefficient_lexicographic": FERMIONIC_AWARE,
    "fermionic_signed_reference": FERMIONIC_AWARE,
    "fermionic_signed": FERMIONIC_AWARE,
    "fermionic_coloring": FERMIONIC_COLORING,
    "fermionic_color": FERMIONIC_COLORING,
}


@dataclass(frozen=True)
class TargetCase:
    case_id: str
    tensor_path: str
    molecule: str
    bond_length: str
    basis: str
    active_occupied: int
    active_vacant: int
    n_qubits: int
    trotter_steps: int
    evolution_time: float
    coefficient_tolerance: float

    @property
    def key(self) -> tuple[str, int, float, float]:
        return (
            self.tensor_path,
            self.trotter_steps,
            self.evolution_time,
            self.coefficient_tolerance,
        )

    @property
    def tensor_abspath(self) -> Path:
        return REPO_ROOT / self.tensor_path


@dataclass(frozen=True)
class Result:
    error: float
    source_csv: str
    ordering_original: str
    source_priority: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basis", default="hgbs-5")
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--time", type=float, default=1.0, dest="evolution_time")
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/four_way_missing_fermionic_coloring_hgbs5.csv"),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path("analysis/four_way_active_space_missing_audit.csv"),
    )
    parser.add_argument(
        "--complete-output",
        type=Path,
        default=Path("analysis/four_way_active_space_complete_hgbs5.csv"),
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Generate absent tensors and calculate missing coloring rows.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Execute only the first N missing rows; useful for a smoke test.",
    )
    parser.add_argument(
        "--molecules",
        nargs="+",
        default=None,
        help="Optional exact molecule filter for execution and auditing.",
    )
    parser.add_argument("--parallel-threshold", type=int, default=2**16)
    parser.add_argument("--no-spin-sector", action="store_true")
    return parser.parse_args()


def resolve_output(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def close_enough(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-6, abs_tol=1.0e-14)


def cross_engine_close(left: float, right: float) -> bool:
    """Allow tiny dense-vs-sector roundoff without hiding real conflicts."""
    return math.isclose(left, right, rel_tol=5.0e-3, abs_tol=1.0e-12)


def source_priority(path: Path, method: str) -> int:
    if path.name == "four_way_missing_fermionic_coloring_hgbs5.csv":
        return 0
    if path in TARGET_CSVS:
        return 1
    if method == FERMIONIC_COLORING and path.name in {
        "l_sweep_trotter_state_t1.csv",
        "polyatomic_trotter_results.csv",
    }:
        return 2
    return 3


def read_targets(args: argparse.Namespace) -> list[TargetCase]:
    targets: dict[tuple[str, int, float, float], TargetCase] = {}
    wanted_molecules = set(args.molecules or [])
    target_specs = [
        (path, FERMIONIC_AWARE, None)
        for path in TARGET_CSVS
    ]
    target_specs.append(
        (
            LIGHT_DIATOMIC_TARGET_CSV,
            FERMIONIC_COLORING,
            LIGHT_DIATOMIC_MOLECULES,
        )
    )
    for path, defining_ordering, allowed_molecules in target_specs:
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if row.get("status") != "success":
                    continue
                if row.get("ordering") != defining_ordering:
                    continue
                if allowed_molecules and row.get("molecule") not in allowed_molecules:
                    continue
                if row.get("basis", "").lower() != args.basis.lower():
                    continue
                if wanted_molecules and row.get("molecule") not in wanted_molecules:
                    continue
                if int(row["trotter_steps"]) != args.steps:
                    continue
                if int(row.get("formula_order") or 1) != 1:
                    continue
                if not close_enough(float(row["evolution_time"]), args.evolution_time):
                    continue
                tolerance = float(row.get("coefficient_tolerance") or 1.0e-12)
                if not close_enough(tolerance, args.tolerance):
                    continue

                target = TargetCase(
                    case_id=row["case_id"],
                    tensor_path=normalized_tensor_path(row["tensor_path"]),
                    molecule=row["molecule"],
                    bond_length=row["bond_length"],
                    basis=row["basis"],
                    active_occupied=int(row["active_occupied"]),
                    active_vacant=int(row["active_vacant"]),
                    n_qubits=int(row["n_qubits"]),
                    trotter_steps=int(row["trotter_steps"]),
                    evolution_time=float(row["evolution_time"]),
                    coefficient_tolerance=tolerance,
                )
                existing = targets.get(target.key)
                if existing is not None and existing != target:
                    raise ValueError(f"Conflicting target metadata for {target.key}")
                targets[target.key] = target

    return sorted(
        targets.values(),
        key=lambda case: (
            case.molecule,
            case.n_qubits,
            case.active_occupied,
            case.active_vacant,
            case.tensor_path,
        ),
    )


def detailed_csv_paths(excluded: set[Path]) -> Iterable[Path]:
    for path in sorted(ANALYSIS_DIR.rglob("*.csv")):
        if path.resolve() not in excluded:
            yield path


def normalized_tensor_path(value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        try:
            return str(path.relative_to(REPO_ROOT))
        except ValueError:
            return str(path)
    return str(path)


def scan_results(
    targets: list[TargetCase],
    excluded: set[Path],
) -> dict[tuple[str, int, float, float], dict[str, Result]]:
    target_keys = {target.key for target in targets}
    results: dict[tuple[str, int, float, float], dict[str, Result]] = {
        key: {} for key in target_keys
    }
    for path in detailed_csv_paths(excluded):
        try:
            stream = path.open(newline="", encoding="utf-8")
        except OSError:
            continue
        with stream:
            reader = csv.DictReader(stream)
            fields = set(reader.fieldnames or [])
            required = {
                "status",
                "tensor_path",
                "ordering",
                "trotter_steps",
                "evolution_time",
                "state_infidelity",
            }
            if not required.issubset(fields):
                continue
            for row in reader:
                if row.get("status") != "success":
                    continue
                original = row.get("ordering", "")
                method = ALIASES.get(original)
                if method is None:
                    continue
                try:
                    key = (
                        normalized_tensor_path(row["tensor_path"]),
                        int(row["trotter_steps"]),
                        float(row["evolution_time"]),
                        float(row.get("coefficient_tolerance") or 1.0e-12),
                    )
                    error = float(row["state_infidelity"])
                    formula_order = int(row.get("formula_order") or 1)
                except (KeyError, TypeError, ValueError):
                    continue
                if key not in target_keys or formula_order != 1:
                    continue
                if not math.isfinite(error) or error < 0.0:
                    continue

                relative_source = str(path.relative_to(REPO_ROOT))
                candidate = Result(
                    error,
                    relative_source,
                    original,
                    source_priority(path, method),
                )
                existing = results[key].get(method)
                if existing is not None:
                    if not cross_engine_close(existing.error, candidate.error):
                        raise ValueError(
                            "Conflicting duplicate result for "
                            f"{key} / {method}: {existing.error} from "
                            f"{existing.source_csv} vs {candidate.error} from "
                            f"{candidate.source_csv}"
                        )
                    if candidate.source_priority < existing.source_priority:
                        results[key][method] = candidate
                    continue
                results[key][method] = candidate
    return results


def config_path_for(case: TargetCase) -> Path:
    tensor = case.tensor_abspath
    filename = tensor.name.removesuffix(".tensors.npz") + "_JW.config"
    return tensor.with_name(filename)


def patched_config(case: TargetCase) -> str:
    config_path = config_path_for(case)
    if not config_path.exists():
        raise FileNotFoundError(f"Missing config for {case.case_id}: {config_path}")
    text = config_path.read_text(encoding="utf-8")
    case_prefix = case.case_id.split("_as-", maxsplit=1)[0]
    file_stub = case.tensor_abspath.parent / case_prefix
    logfile = case.tensor_abspath.with_name(case.case_id + "_JW.log")
    text, stub_count = re.subn(
        r"(?m)^general\.file_stub\s*=.*$",
        f"general.file_stub = {str(file_stub)!r}",
        text,
    )
    text, log_count = re.subn(
        r"(?m)^general\.logfile\s*=.*$",
        f"general.logfile = {str(logfile)!r}",
        text,
    )
    if stub_count != 1 or log_count != 1:
        raise ValueError(
            f"Could not patch exactly one file_stub/logfile in {config_path}"
        )
    return text


def generate_tensor(case: TargetCase) -> None:
    if case.tensor_abspath.exists():
        return
    case.tensor_abspath.parent.mkdir(parents=True, exist_ok=True)
    config_text = patched_config(case)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".config",
            prefix=f"{case.case_id}-",
            dir=tempfile.gettempdir(),
            encoding="utf-8",
            delete=False,
        ) as stream:
            stream.write(config_text)
            temporary_path = Path(stream.name)

        environment = os.environ.copy()
        pythonpath = str(REPO_ROOT.parent)
        if environment.get("PYTHONPATH"):
            pythonpath += os.pathsep + environment["PYTHONPATH"]
        environment["PYTHONPATH"] = pythonpath
        subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "hamiltonian_generator" / "hamgen.py"),
                str(temporary_path),
            ],
            cwd=REPO_ROOT,
            env=environment,
            check=True,
        )
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    if not case.tensor_abspath.exists():
        raise FileNotFoundError(
            f"hamgen completed without creating {case.tensor_abspath}"
        )


def write_failure(
    writer: csv.DictWriter,
    stream: Any,
    case: TargetCase,
    error: Exception,
) -> None:
    row = matrix_free.blank_row()
    row.update(
        {
            "status": "failed",
            "error_message": f"{type(error).__name__}: {error}",
            "case_id": case.case_id,
            "tensor_path": case.tensor_path,
            "molecule": case.molecule,
            "bond_length": case.bond_length,
            "basis": case.basis,
            "active_occupied": case.active_occupied,
            "active_vacant": case.active_vacant,
            "n_qubits": case.n_qubits,
            "coefficient_tolerance": case.coefficient_tolerance,
        }
    )
    writer.writerow(row)
    stream.flush()


def execute_missing(
    targets: list[TargetCase],
    results: dict[tuple[str, int, float, float], dict[str, Result]],
    args: argparse.Namespace,
    output: Path,
) -> None:
    missing = [target for target in targets if len(results[target.key]) < len(METHODS)]
    if args.limit is not None:
        missing = missing[: args.limit]
    if not missing:
        print("No missing four-way rows require execution.")
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    completed, raw_infidelities = matrix_free.load_resume_data(output)
    write_header = not output.exists() or output.stat().st_size == 0
    matrix_free.warm_up_numba()

    with output.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=matrix_free.FIELDNAMES,
            lineterminator="\n",
        )
        if write_header:
            writer.writeheader()
            stream.flush()

        for index, case in enumerate(missing, start=1):
            print()
            print("=" * 88)
            print(f"[{index}/{len(missing)}] {case.case_id}")
            present = results[case.key]
            missing_methods = [method for method in METHODS if method not in present]
            print("  Missing methods: " + ", ".join(missing_methods))
            raw = present.get(JW_RAW)
            raw_key = (
                case.case_id,
                1,
                case.trotter_steps,
                case.evolution_time,
            )
            if raw is not None:
                raw_infidelities[raw_key] = raw.error
                completed.add(
                    (
                        case.case_id,
                        JW_RAW,
                        1,
                        case.trotter_steps,
                        case.evolution_time,
                    )
                )
            try:
                print("  Generating missing tensor..." if not case.tensor_abspath.exists() else "  Reusing tensor...")
                generate_tensor(case)

                matrix_free_methods = [
                    method
                    for method in (JW_RAW, FERMIONIC_COLORING)
                    if method in missing_methods
                ]
                if matrix_free_methods:
                    requested = [JW_RAW]
                    if FERMIONIC_COLORING in matrix_free_methods:
                        requested.append(FERMIONIC_COLORING)
                    matrix_free.benchmark_case(
                        tensor_path=Path(case.tensor_path),
                        requested_orderings=requested,
                        formula_orders=[1],
                        steps_list=[case.trotter_steps],
                        evolution_time=case.evolution_time,
                        tolerance=case.coefficient_tolerance,
                        spin_preserving=not args.no_spin_sector,
                        parallel_threshold=args.parallel_threshold,
                        completed=completed,
                        raw_infidelities=raw_infidelities,
                        writer=writer,
                        output_stream=stream,
                    )

                deterministic_methods = [
                    method
                    for method in (JW_MAGNITUDE, FERMIONIC_AWARE)
                    if method in missing_methods
                ]
                if deterministic_methods:
                    raw_reference = (
                        {"state_infidelity": raw.error}
                        if raw is not None
                        else {}
                    )
                    deterministic_rows = baseline.benchmark_case(
                        tensor_path=Path(case.tensor_path),
                        args=args,
                        ordering_names=deterministic_methods,
                        raw_reference=raw_reference,
                    )
                    for deterministic_row in deterministic_rows:
                        output_row = matrix_free.blank_row()
                        output_row.update(
                            {
                                field: deterministic_row[field]
                                for field in matrix_free.FIELDNAMES
                                if field in deterministic_row
                            }
                        )
                        output_row.update(
                            {
                                "case_id": case.case_id,
                                "tensor_path": case.tensor_path,
                                "formula_order": 1,
                            }
                        )
                        writer.writerow(output_row)
                        stream.flush()
            except Exception as error:
                print(f"FAILED {case.case_id}: {type(error).__name__}: {error}")
                write_failure(writer, stream, case, error)


AUDIT_FIELDS = [
    "case_id",
    "tensor_path",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "trotter_steps",
    "evolution_time",
    "coefficient_tolerance",
    "jw_raw_present",
    "jw_magnitude_present",
    "fermionic_aware_present",
    "fermionic_coloring_present",
    "missing_methods",
    "all_four_present",
]

COMPLETE_FIELDS = [
    "status",
    "case_id",
    "tensor_path",
    "molecule",
    "bond_length",
    "basis",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "ordering",
    "formula_order",
    "trotter_steps",
    "evolution_time",
    "state_infidelity",
    "coefficient_tolerance",
    "source_csv",
    "ordering_original",
]


def write_audit_and_complete(
    targets: list[TargetCase],
    results: dict[tuple[str, int, float, float], dict[str, Result]],
    audit_output: Path,
    complete_output: Path,
) -> tuple[int, int]:
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    audit_rows: list[dict[str, Any]] = []
    complete_rows: list[dict[str, Any]] = []
    complete_count = 0
    for case in targets:
        methods = results[case.key]
        missing = [method for method in METHODS if method not in methods]
        if not missing:
            complete_count += 1
        audit_rows.append(
            {
                "case_id": case.case_id,
                "tensor_path": case.tensor_path,
                "molecule": case.molecule,
                "bond_length": case.bond_length,
                "basis": case.basis,
                "active_occupied": case.active_occupied,
                "active_vacant": case.active_vacant,
                "n_qubits": case.n_qubits,
                "trotter_steps": case.trotter_steps,
                "evolution_time": case.evolution_time,
                "coefficient_tolerance": case.coefficient_tolerance,
                "jw_raw_present": JW_RAW in methods,
                "jw_magnitude_present": JW_MAGNITUDE in methods,
                "fermionic_aware_present": FERMIONIC_AWARE in methods,
                "fermionic_coloring_present": FERMIONIC_COLORING in methods,
                "missing_methods": ";".join(missing),
                "all_four_present": not missing,
            }
        )
        for method in METHODS:
            result = methods.get(method)
            if result is None:
                continue
            complete_rows.append(
                {
                    "status": "success",
                    "case_id": case.case_id,
                    "tensor_path": case.tensor_path,
                    "molecule": case.molecule,
                    "bond_length": case.bond_length,
                    "basis": case.basis,
                    "active_occupied": case.active_occupied,
                    "active_vacant": case.active_vacant,
                    "n_qubits": case.n_qubits,
                    "ordering": method,
                    "formula_order": 1,
                    "trotter_steps": case.trotter_steps,
                    "evolution_time": case.evolution_time,
                    "state_infidelity": result.error,
                    "coefficient_tolerance": case.coefficient_tolerance,
                    "source_csv": result.source_csv,
                    "ordering_original": result.ordering_original,
                }
            )

    with audit_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=AUDIT_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)
    with complete_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=COMPLETE_FIELDS,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(complete_rows)
    return complete_count, len(complete_rows)


def main() -> None:
    args = parse_args()
    if args.steps <= 0 or args.evolution_time <= 0.0 or args.tolerance <= 0.0:
        raise ValueError("Steps, time, and tolerance must be positive.")
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be positive.")

    output = resolve_output(args.output)
    audit_output = resolve_output(args.audit_output)
    complete_output = resolve_output(args.complete_output)
    excluded = {audit_output.resolve(), complete_output.resolve()}

    targets = read_targets(args)
    results = scan_results(targets, excluded)
    method_counts = {
        method: sum(method in results[target.key] for target in targets)
        for method in METHODS
    }
    print(f"Target exact cases: {len(targets)}")
    for method in METHODS:
        print(
            f"Existing {method}: {method_counts[method]}/{len(targets)}; "
            f"missing {len(targets) - method_counts[method]}"
        )

    if args.execute:
        execute_missing(targets, results, args, output)
        results = scan_results(targets, excluded)

    complete_count, long_rows = write_audit_and_complete(
        targets,
        results,
        audit_output,
        complete_output,
    )
    remaining = len(targets) - complete_count
    print(f"Complete four-way cases: {complete_count}/{len(targets)}")
    print(f"Remaining incomplete cases: {remaining}")
    print(f"Consolidated long rows: {long_rows}/{len(targets) * len(METHODS)}")
    print(f"Audit: {audit_output.relative_to(REPO_ROOT)}")
    print(f"Complete data: {complete_output.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
