#!/usr/bin/env python3
"""Create controlled fixed-active-occupied B2/STO-6G configurations."""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

OUTPUT_ROOT = (
    REPO_ROOT
    / "hamiltonian_generator"
    / "b2_controlled_active_space_libraries"
)

MOLECULE = "B-B"
ATOM = "B"
BOND_LENGTH_TEXT = "1.70"
BOND_LENGTH = 1.70
BASIS = "sto-6g"

SERIES = {
    "occ-002": {
        "active_occupied": 2,
        "active_vacant": [2, 4, 6, 8, 10],
    },
    "occ-004": {
        "active_occupied": 4,
        "active_vacant": [2, 4, 6, 8, 10],
    },
    "occ-006": {
        "active_occupied": 6,
        "active_vacant": [2, 4, 6, 8, 10],
    },
}


def write_config(
    series_name: str,
    active_occupied: int,
    active_vacant: int,
) -> Path:
    output_directory = (
        OUTPUT_ROOT
        / series_name
        / MOLECULE
        / BOND_LENGTH_TEXT
        / BASIS
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    base_name = f"{MOLECULE}_{BOND_LENGTH_TEXT}_{BASIS}"
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

    half_length = BOND_LENGTH / 2.0

    contents = "\n".join(
        [
            "general.print_verbose()",
            f"general.logfile = {str(logfile.relative_to(REPO_ROOT))!r}",
            f"general.file_stub = {str(file_stub.relative_to(REPO_ROOT))!r}",
            'general.file_format = "default"',
            f"L = {BOND_LENGTH}",
            (
                f'hamiltonian.add_atom("{ATOM}", '
                f"-{half_length}, 0.0, 0.0)"
            ),
            (
                f'hamiltonian.add_atom("{ATOM}", '
                f"{half_length}, 0.0, 0.0)"
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

    for series_name, settings in SERIES.items():
        active_occupied = settings["active_occupied"]

        for active_vacant in settings["active_vacant"]:
            n_qubits = active_occupied + active_vacant

            print(
                f"{series_name}: occupied={active_occupied}, "
                f"vacant={active_vacant}, qubits={n_qubits}"
            )

            write_config(
                series_name=series_name,
                active_occupied=active_occupied,
                active_vacant=active_vacant,
            )
            count += 1

    print()
    print(f"Configurations checked: {count}")
    print(f"Output root: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
