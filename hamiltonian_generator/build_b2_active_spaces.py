#!/usr/bin/env python3
"""Create the nine B2/STO-6G active-space configurations and tensor files.

Place this file in ``hamiltonian_generator/`` on the QHAT ``L-sweep`` branch,
or run it from the repository root after installing QHAT with ``pip install -e .``.

The generated active spaces are the exact sequence used by the current
``build_config_L_sweep.py`` 40/60 heuristic for B2:

    (2, 2), (2, 4), (4, 4), (4, 6), (4, 8),
    (6, 8), (6, 10), (8, 10), (10, 10)

The two numbers are active occupied and active vacant spin orbitals, so their
sum is the qubit count.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from qhat.common.logging_utils import configure_logging
from qhat.hamiltonian_generator.hamgen import get_ham2
from qhat.hamiltonian_generator.hamgen_types import (
    GeneralConfigurationUser,
    HamiltonianConfiguration,
    State,
)


ACTIVE_SPACES: tuple[tuple[int, int], ...] = (
    (2, 2),
    (2, 4),
    (4, 4),
    (4, 6),
    (4, 8),
    (6, 8),
    (6, 10),
    (8, 10),
    (10, 10),
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the nine B2/STO-6G/JW active-space configs and "
            "their QHAT tensor NPZ files."
        )
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/b2_active_space_library"),
        help="Output library root.",
    )
    parser.add_argument(
        "--bond-length",
        type=float,
        default=None,
        help=(
            "Fixed B-B separation in Angstrom. By default, use the same "
            "Pyykko covalent-radius estimate as build_config_L_sweep.py "
            "(2 * 0.01 * B.covalent_radius_pyykko)."
        ),
    )
    parser.add_argument(
        "--max-qubits",
        type=int,
        default=20,
        help="Generate only active spaces at or below this size.",
    )
    parser.add_argument(
        "--configs-only",
        action="store_true",
        help="Write config files without running get_ham2 or saving tensors.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate tensor files that already exist.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/b2_active_space_tensor_summary.csv"),
        help="CSV summary path.",
    )
    return parser.parse_args()


def resolve_bond_length(requested_bond_length: float | None) -> float:
    """Use the L-sweep branch's Pyykko covalent-radius estimate by default."""
    if requested_bond_length is not None:
        return requested_bond_length

    import mendeleev

    boron = mendeleev.element("B")
    radius_pm = boron.covalent_radius_pyykko
    if radius_pm is None:
        raise ValueError(
            "Mendeleev did not provide B.covalent_radius_pyykko; "
            "pass --bond-length explicitly."
        )
    return round(2.0 * 0.01 * float(radius_pm), 4)


def validate_arguments(args: argparse.Namespace) -> None:
    if args.bond_length <= 0.0:
        raise ValueError("--bond-length must be positive.")
    if args.max_qubits < 4 or args.max_qubits > 20:
        raise ValueError("--max-qubits must be between 4 and 20.")


def write_configuration(
    library_root: Path,
    bond_length: float,
    active_occupied: int,
    active_vacant: int,
) -> Path:
    """Write one B2/STO-6G/JW QHAT configuration file."""
    molecule = "B-B"
    basis = "sto-6g"
    mapping = "JW"
    bond_text = f"{bond_length:.2f}"
    stub = f"{molecule}_{bond_text}_{basis}"
    extended = (
        f"{stub}_as-{active_occupied:03d}-{active_vacant:03d}_{mapping}"
    )

    output_directory = library_root / molecule / bond_text / basis
    output_directory.mkdir(parents=True, exist_ok=True)

    absolute_stub = str((output_directory.resolve() / stub))
    config_path = output_directory / f"{extended}.config"
    half_length = 0.5 * bond_length

    lines = [
        "general.print_verbose()",
        (
            'general.logfile = '
            f'"{absolute_stub}_as-{active_occupied:03d}-'
            f'{active_vacant:03d}_{mapping}.log"'
        ),
        f'general.file_stub = "{absolute_stub}"',
        'general.file_format = "default"',
        f"L = {bond_length!r}",
        f'hamiltonian.add_atom("B", {-half_length!r}, 0.0, 0.0)',
        f'hamiltonian.add_atom("B", {half_length!r}, 0.0, 0.0)',
        'hamiltonian.basis = "sto-6g"',
        f"hamiltonian.num_active_occupied = {active_occupied}",
        f"hamiltonian.num_active_vacant = {active_vacant}",
        'hamiltonian.f2q_mapping = "JW"',
        "",
    ]
    config_path.write_text("\n".join(lines), encoding="utf-8")
    return config_path


def load_configuration_file(config_path: Path) -> State:
    """Load one QHAT Python configuration without invoking hamgen.py CLI."""
    config_script = config_path.read_text(encoding="utf-8").rstrip("\n")
    general = GeneralConfigurationUser()
    hamiltonian = HamiltonianConfiguration()
    namespace: dict[str, Any] = {
        "general": general,
        "hamiltonian": hamiltonian,
    }
    exec(compile(config_script, str(config_path), "exec"), namespace)
    return State(config_script, general, hamiltonian)


def tensor_path_for_state(state: State) -> Path:
    """Return the ``*.tensors.npz`` path corresponding to the active space."""
    return Path(state.filename_ham2()).with_suffix(".tensors.npz")


def save_tensors(tensor_path: Path, active_hamiltonian: Any) -> None:
    """Save the active-space InteractionOperator arrays."""
    import numpy as np

    tensor_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        tensor_path,
        constant=active_hamiltonian.constant,
        one_body=active_hamiltonian.one_body_tensor,
        two_body=active_hamiltonian.two_body_tensor,
    )


def main() -> None:
    args = parse_arguments()
    args.bond_length = resolve_bond_length(args.bond_length)
    validate_arguments(args)

    configure_logging(level="verbose", logfile=None)
    print(f"Using fixed B-B bond length: {args.bond_length:.4f} Angstrom")

    selected_spaces = [
        pair for pair in ACTIVE_SPACES if sum(pair) <= args.max_qubits
    ]
    if not selected_spaces:
        raise ValueError("No active spaces satisfy --max-qubits.")

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "config_path",
        "tensor_path",
        "bond_length",
        "basis",
        "mapping",
        "active_occupied",
        "active_vacant",
        "n_qubits",
        "status",
        "error",
    ]

    successful = 0
    skipped = 0
    failed = 0

    with args.summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()

        for case_index, (active_occupied, active_vacant) in enumerate(
            selected_spaces,
            start=1,
        ):
            n_qubits = active_occupied + active_vacant
            print()
            print("=" * 78)
            print(
                f"[{case_index}/{len(selected_spaces)}] "
                f"B2/STO-6G as-{active_occupied:03d}-{active_vacant:03d} "
                f"({n_qubits} qubits)"
            )

            config_path = write_configuration(
                args.library,
                args.bond_length,
                active_occupied,
                active_vacant,
            )
            row = {
                "config_path": str(config_path),
                "tensor_path": "",
                "bond_length": f"{args.bond_length:.8f}",
                "basis": "sto-6g",
                "mapping": "JW",
                "active_occupied": active_occupied,
                "active_vacant": active_vacant,
                "n_qubits": n_qubits,
                "status": "",
                "error": "",
            }

            print(f"Config: {config_path}")

            if args.configs_only:
                row["status"] = "config_only"
                writer.writerow(row)
                stream.flush()
                successful += 1
                continue

            try:
                state = load_configuration_file(config_path)
                tensor_path = tensor_path_for_state(state)
                row["tensor_path"] = str(tensor_path)

                if tensor_path.exists() and not args.overwrite:
                    row["status"] = "skipped"
                    print(f"Tensor already exists: {tensor_path}")
                    skipped += 1
                else:
                    active_hamiltonian = get_ham2(state)
                    actual_n_qubits = int(active_hamiltonian.n_qubits)
                    if actual_n_qubits != n_qubits:
                        raise ValueError(
                            "Generated tensor size does not match the requested "
                            f"active space: {actual_n_qubits} != {n_qubits}."
                        )
                    save_tensors(tensor_path, active_hamiltonian)
                    row["status"] = "success"
                    print(f"Saved tensor: {tensor_path}")
                    successful += 1
            except Exception as error:  # Keep the full sweep resumable.
                row["status"] = "failed"
                row["error"] = f"{type(error).__name__}: {error}"
                print(f"FAILED: {row['error']}")
                failed += 1

            writer.writerow(row)
            stream.flush()

    print()
    print("=" * 78)
    print(f"Active spaces selected: {len(selected_spaces)}")
    print(f"Successful/config-only: {successful}")
    print(f"Already existed:        {skipped}")
    print(f"Failed:                  {failed}")
    print(f"Summary:                 {args.summary}")
    print("=" * 78)


if __name__ == "__main__":
    main()
