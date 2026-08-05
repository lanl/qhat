#!/usr/bin/env python3
"""Create HGBS-5 basis-transfer configurations for Li2, Be2, F2, and B2."""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

OUTPUT_ROOT = (
    REPO_ROOT
    / "hamiltonian_generator"
    / "hgbs5_diatomic_basis_transfer"
)

BASIS = "hgbs-5"

CASES = {
    "Li-Li": {
        "atom": "Li",
        "bond_length": "2.66",
        "active_spaces": [
            ("lower", 6, 2),    # 8 qubits
            ("lower", 6, 4),    # 10
            ("lower", 6, 6),    # 12
            ("screen", 6, 8),   # 14
            ("night", 6, 10),   # 16
            ("night", 6, 12),   # 18
            ("night", 6, 14),   # 20
        ],
    },
    "Be-Be": {
        "atom": "Be",
        "bond_length": "2.04",
        "active_spaces": [
            ("lower", 8, 2),    # 10 qubits
            ("lower", 8, 4),    # 12
            ("screen", 8, 6),   # 14
            ("night", 8, 8),    # 16
            ("night", 8, 10),   # 18
            ("night", 8, 12),   # 20
        ],
    },
    "F-F": {
        "atom": "F",
        "bond_length": "1.28",
        "active_spaces": [
            ("lower", 4, 2),    # 6 qubits
            ("lower", 6, 2),    # 8
            ("lower", 8, 2),    # 10
            ("lower", 10, 2),   # 12
            ("screen", 12, 2),  # 14
            ("night", 14, 2),   # 16
            ("night", 16, 2),   # 18
            ("night", 18, 2),   # 20
        ],
    },
    "B-B": {
        "atom": "B",
        "bond_length": "1.70",
        "active_spaces": [
            ("lower", 2, 2),    # 4 qubits
            ("lower", 2, 4),    # 6
            ("lower", 4, 4),    # 8
            ("lower", 4, 6),    # 10
            ("lower", 4, 8),    # 12
            ("screen", 6, 8),   # 14
            ("night", 6, 10),   # 16
            ("night", 8, 10),   # 18
            ("night", 10, 10),  # 20
        ],
    },
}


def write_config(
    *,
    phase: str,
    molecule: str,
    atom: str,
    bond_length: str,
    active_occupied: int,
    active_vacant: int,
) -> Path:
    output_directory = (
        OUTPUT_ROOT
        / phase
        / molecule
        / bond_length
        / BASIS
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    base_name = f"{molecule}_{bond_length}_{BASIS}"
    active_suffix = (
        f"as-{active_occupied:03d}-{active_vacant:03d}"
    )

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
            f"general.logfile = {str(logfile.resolve())!r}",
            f"general.file_stub = {str(file_stub.resolve())!r}",
            'general.file_format = "default"',
            f"L = {float(bond_length)}",
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
        for phase, occupied, vacant in settings["active_spaces"]:
            n_qubits = occupied + vacant

            print(
                f"{phase:6s} {molecule:5s}: "
                f"occupied={occupied:2d}, "
                f"vacant={vacant:2d}, "
                f"qubits={n_qubits:2d}"
            )

            write_config(
                phase=phase,
                molecule=molecule,
                atom=settings["atom"],
                bond_length=settings["bond_length"],
                active_occupied=occupied,
                active_vacant=vacant,
            )
            count += 1

    print()
    print(f"Configurations checked: {count}")
    print(f"Output root: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()