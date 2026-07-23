#!/usr/bin/env python3
"""
Generate active-space tensor NPZ files for all L-sweep configurations.

This script stops before the Jordan-Wigner transformation. It creates:

    *_as-XXX-XXX.tensors.npz

containing:

    constant
    one_body
    two_body
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import numpy as np

from qhat.hamiltonian_generator.hamgen import get_ham2
from qhat.hamiltonian_generator.hamgen_types import (
    GeneralConfigurationUser,
    HamiltonianConfiguration,
    State,
)
from qhat.common.logging_utils import configure_logging


def load_configuration_file(config_path: Path) -> State:
    """
    Load one QHAT configuration file without using command-line
    arguments from hamgen.py.
    """
    config_script = config_path.read_text(encoding="utf-8").rstrip("\n")

    general = GeneralConfigurationUser()
    hamiltonian = HamiltonianConfiguration()

    namespace: dict[str, Any] = {
        "general": general,
        "hamiltonian": hamiltonian,
    }

    exec(
        compile(config_script, str(config_path), "exec"),
        namespace,
    )

    return State(
        config_script,
        general,
        hamiltonian,
    )


def get_tensor_path(state: State) -> Path:
    """
    Convert:

        ..._as-002-002.pickle

    into:

        ..._as-002-002.tensors.npz
    """
    ham2_path = Path(state.filename_ham2())
    return ham2_path.with_suffix(".tensors.npz")


def save_tensors(
    tensor_path: Path,
    active_hamiltonian: Any,
) -> None:
    """Save the active-space InteractionOperator tensors."""
    tensor_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    np.savez_compressed(
        tensor_path,
        constant=active_hamiltonian.constant,
        one_body=active_hamiltonian.one_body_tensor,
        two_body=active_hamiltonian.two_body_tensor,
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate active-space tensor NPZ files from "
            "L-sweep configuration files."
        )
    )

    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/library"),
        help="Directory containing the L-sweep configurations.",
    )

    parser.add_argument(
        "--summary",
        type=Path,
        default=Path(
            "hamiltonian_generator/l_sweep_tensor_summary.csv"
        ),
        help="CSV file recording successful and failed cases.",
    )

    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N configurations.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate an NPZ file even when it already exists.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    configure_logging(
        level="verbose",
        logfile=None,
    )

    config_paths = sorted(
        args.library.rglob("*_JW.config")
    )

    total_discovered = len(config_paths)

    print(f"JW configurations discovered: {total_discovered}")

    if total_discovered != 156:
        print(
            "Warning: expected 156 configurations, "
            f"but found {total_discovered}."
        )

    if args.limit is not None:
        config_paths = config_paths[: args.limit]
        print(f"Testing only the first {len(config_paths)} cases.")

    args.summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "config_path",
        "tensor_path",
        "molecule",
        "bond_length",
        "basis",
        "active_occupied",
        "active_vacant",
        "n_qubits",
        "status",
        "error",
    ]

    successful = 0
    skipped = 0
    failed = 0

    with args.summary.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as summary_file:
        writer = csv.DictWriter(
            summary_file,
            fieldnames=fieldnames,
        )
        writer.writeheader()

        for index, config_path in enumerate(
            config_paths,
            start=1,
        ):
            print()
            print(
                f"[{index}/{len(config_paths)}] "
                f"{config_path}"
            )

            row = {
                "config_path": str(config_path),
                "tensor_path": "",
                "molecule": "",
                "bond_length": "",
                "basis": "",
                "active_occupied": "",
                "active_vacant": "",
                "n_qubits": "",
                "status": "",
                "error": "",
            }

            try:
                state = load_configuration_file(
                    config_path
                )

                geometry = (
                    state.config_hamiltonian.geometry()
                )

                row["molecule"] = "-".join(
                    atom[0] for atom in geometry
                )
                row["basis"] = (
                    state.config_hamiltonian.basis
                )
                row["active_occupied"] = (
                    state.config_hamiltonian
                    .num_active_occupied
                )
                row["active_vacant"] = (
                    state.config_hamiltonian
                    .num_active_vacant
                )

                tensor_path = get_tensor_path(state)
                row["tensor_path"] = str(tensor_path)

                if (
                    tensor_path.exists()
                    and not args.overwrite
                ):
                    row["status"] = "skipped"
                    row["n_qubits"] = (
                        row["active_occupied"]
                        + row["active_vacant"]
                    )

                    print(
                        "Tensor file already exists; "
                        "skipping."
                    )
                    skipped += 1

                else:
                    # Build or load the active-space
                    # InteractionOperator.
                    active_hamiltonian = get_ham2(state)

                    # Explicitly save the NPZ. This also handles
                    # cases where get_ham2 loaded an existing
                    # pickle but the NPZ was absent.
                    save_tensors(
                        tensor_path,
                        active_hamiltonian,
                    )

                    row["bond_length"] = (
                        active_hamiltonian.separation
                    )
                    row["n_qubits"] = (
                        active_hamiltonian.n_qubits
                    )
                    row["status"] = "success"

                    print(f"Saved: {tensor_path}")
                    successful += 1

            except Exception as error:
                row["status"] = "failed"
                row["error"] = (
                    f"{type(error).__name__}: {error}"
                )

                print(f"FAILED: {row['error']}")
                failed += 1

            writer.writerow(row)
            summary_file.flush()

    print()
    print("=" * 60)
    print(f"Configurations processed: {len(config_paths)}")
    print(f"Successful:              {successful}")
    print(f"Already existed:         {skipped}")
    print(f"Failed:                  {failed}")
    print(f"Summary CSV:             {args.summary}")
    print("=" * 60)


if __name__ == "__main__":
    main()