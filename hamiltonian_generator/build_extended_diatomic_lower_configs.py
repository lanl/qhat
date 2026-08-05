#!/usr/bin/env python3
"""Create lower-qubit Li2, Be2, and F2 STO-6G configurations."""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

LIBRARY = (
    REPO_ROOT
    / "hamiltonian_generator"
    / "extended_diatomic_20q_library"
)

BASIS = "sto-6g"

CASES = {
    "Li-Li": {
        "atom": "Li",
        "bond_length": "2.66",
        # Preserve six active occupied orbitals.
        "active_spaces": [
            (6, 2),   # 8 qubits
            (6, 4),   # 10 qubits
            (6, 6),   # 12 qubits
            (6, 8),   # 14 qubits
            (6, 10),  # 16 qubits
        ],
    },
    "Be-Be": {
        "atom": "Be",
        "bond_length": "2.04",
        # Preserve eight active occupied orbitals.
        "active_spaces": [
            (8, 2),  # 10 qubits
            (8, 4),  # 12 qubits
            (8, 6),  # 14 qubits
            (8, 8),  # 16 qubits
        ],
    },
    "F-F": {
        "atom": "F",
        "bond_length": "1.28",
        # Preserve two active vacant orbitals.
        "active_spaces": [
            (4, 2),   # 6 qubits
            (6, 2),   # 8 qubits
            (8, 2),   # 10 qubits
            (10, 2),  # 12 qubits
            (12, 2),  # 14 qubits
            (14, 2),  # 16 qubits
        ],
    },
}


def write_config(
    *,
    molecule: str,
    atom: str,
    bond_length: str,
    active_occupied: int,
    active_vacant: int,
) -> Path:
    output_directory = (
        LIBRARY
        / molecule
        / bond_length
        / BASIS
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    active_suffix = (
        f"as-{active_occupied:03d}-{active_vacant:03d}"
    )
    base_name = f"{molecule}_{bond_length}_{BASIS}"

    config_path = (
        output_directory
        / f"{base_name}_{active_suffix}_JW.config"
    )
    logfile = (
        output_directory
        / f"{base_name}_{active_suffix}_JW.log"
    )
    file_stub = output_directory / base_name

    contents = "\n".join(
        [
            "general.print_verbose()",
            f"general.logfile = {str(logfile)!r}",
            f"general.file_stub = {str(file_stub)!r}",
            'general.file_format = "default"',
            f"L = {bond_length}",
            (
                f'hamiltonian.add_atom("{atom}", '
                "-0.5 * L, 0.0, 0.0)"
            ),
            (
                f'hamiltonian.add_atom("{atom}", '
                " 0.5 * L, 0.0, 0.0)"
            ),
            f'hamiltonian.basis = "{BASIS}"',
            (
                "hamiltonian.num_active_occupied = "
                f"{active_occupied}"
            ),
            (
                "hamiltonian.num_active_vacant = "
                f"{active_vacant}"
            ),
            'hamiltonian.f2q_mapping = "JW"',
            "",
        ]
    )

    if config_path.exists():
        print(f"Already exists: {config_path}")
    else:
        config_path.write_text(contents, encoding="utf-8")
        print(f"Created:        {config_path}")

    return config_path


def main() -> None:
    count = 0

    for molecule, settings in CASES.items():
        for active_occupied, active_vacant in settings[
            "active_spaces"
        ]:
            write_config(
                molecule=molecule,
                atom=settings["atom"],
                bond_length=settings["bond_length"],
                active_occupied=active_occupied,
                active_vacant=active_vacant,
            )
            count += 1

    print()
    print(f"Configurations checked: {count}")
    print(f"Library: {LIBRARY}")


if __name__ == "__main__":
    main()