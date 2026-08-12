#!/usr/bin/env python3
"""Generate missing C2/N2 QHAT configs and active-space tensor files.

This script is designed to extend the existing C2/N2 case study to both
STO-6G and HGBS-5 and to fill the even active-space sizes from 4 through
20 qubits (or a user-selected subset).

The active occupied/vacant sequence follows the same 40/60 heuristic used by
``build_config_L_sweep.py``.  Existing configs and tensor files are preserved
by default, so the script is safe to rerun when only some qubit sizes are
missing.

Default fixed geometries use the same L-sweep convention as the other
homonuclear diatomic studies: the sum of the two Pyykko covalent radii.
With the current mendeleev data this is approximately C-C = 1.50 Angstrom
and N-N = 1.42 Angstrom.  Use ``--bond-length`` to override either value.

Run from the QHAT repository root after ``python -m pip install -e .``.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any, Sequence

import basis_set_exchange
import mendeleev

from qhat.common.logging_utils import configure_logging
from qhat.hamiltonian_generator.hamgen import get_ham2
from qhat.hamiltonian_generator.hamgen_types import (
    GeneralConfigurationUser,
    HamiltonianConfiguration,
    State,
)


MOLECULES: dict[str, tuple[str, str]] = {
    "C-C": ("C", "C"),
    "N-N": ("N", "N"),
}
DEFAULT_BASES = ("sto-6g", "hgbs-5")
DEFAULT_QUBITS = tuple(range(4, 21, 2))
MAPPING = "JW"


SUMMARY_FIELDS = [
    "molecule",
    "bond_length",
    "basis",
    "total_electrons",
    "total_spin_orbitals",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "config_path",
    "tensor_path",
    "config_status",
    "tensor_status",
    "error",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/comparable_diatomic_library"),
        help="Output library root.",
    )
    parser.add_argument(
        "--molecules",
        nargs="+",
        default=list(MOLECULES),
        choices=sorted(MOLECULES),
        help="Molecules to generate. Default: C-C N-N.",
    )
    parser.add_argument(
        "--bases",
        nargs="+",
        default=list(DEFAULT_BASES),
        help="Basis sets to generate. Default: sto-6g hgbs-5.",
    )
    parser.add_argument(
        "--qubits",
        nargs="+",
        type=int,
        default=list(DEFAULT_QUBITS),
        help="Requested even active-space sizes. Default: 4 6 ... 20.",
    )
    parser.add_argument(
        "--bond-length",
        action="append",
        default=[],
        metavar="MOLECULE=VALUE",
        help=(
            "Override a fixed bond length in Angstrom, for example "
            "--bond-length C-C=1.50 --bond-length N-N=1.42. "
            "May be supplied more than once."
        ),
    )
    parser.add_argument(
        "--configs-only",
        action="store_true",
        help="Write/verify configs but do not call get_ham2 or save tensors.",
    )
    parser.add_argument(
        "--overwrite-configs",
        action="store_true",
        help="Rewrite config files that already exist.",
    )
    parser.add_argument(
        "--overwrite-tensors",
        action="store_true",
        help="Regenerate tensor files that already exist.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/c2_n2_case_study/tensor_generation_summary.csv"),
        help="CSV generation summary.",
    )
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    if not args.qubits:
        raise ValueError("--qubits cannot be empty.")
    if any(q <= 0 or q % 2 for q in args.qubits):
        raise ValueError("--qubits must contain positive even integers.")
    if len(set(args.qubits)) != len(args.qubits):
        raise ValueError("--qubits contains duplicate values.")
    if not args.bases:
        raise ValueError("--bases cannot be empty.")


def parse_bond_length_overrides(values: Sequence[str]) -> dict[str, float]:
    overrides: dict[str, float] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(
                "--bond-length must use MOLECULE=VALUE, e.g. C-C=1.50."
            )
        molecule, text = value.split("=", 1)
        molecule = molecule.strip()
        if molecule not in MOLECULES:
            raise ValueError(
                f"Unknown molecule in --bond-length: {molecule!r}. "
                f"Choose from {sorted(MOLECULES)}."
            )
        length = float(text)
        if length <= 0.0:
            raise ValueError("Bond lengths must be positive.")
        overrides[molecule] = length
    return overrides


def pyykko_bond_length(molecule: str) -> float:
    atom1, atom2 = MOLECULES[molecule]
    radii_angstrom: list[float] = []
    for symbol in (atom1, atom2):
        element = mendeleev.element(symbol)
        radius_pm = element.covalent_radius_pyykko
        if radius_pm is None:
            raise ValueError(
                f"No Pyykko covalent radius for {symbol}; "
                "pass --bond-length explicitly."
            )
        radii_angstrom.append(0.01 * float(radius_pm))
    return round(sum(radii_angstrom), 4)


def resolved_bond_lengths(overrides: dict[str, float]) -> dict[str, float]:
    return {
        molecule: overrides.get(molecule, pyykko_bond_length(molecule))
        for molecule in MOLECULES
    }


def norb_for_shell(shell: str) -> int:
    # Spin-orbital counts, matching build_config_L_sweep.py.
    return {"s": 2, "p": 6, "d": 10, "f": 14, "g": 18}[shell]


def atomic_spin_orbital_count(basis: str, atomic_number: int) -> int:
    """Return the number of spin orbitals supplied by one atomic basis."""
    try:
        basis_string = basis_set_exchange.get_basis(
            basis,
            elements=atomic_number,
            fmt="nwchem",
        )
    except (KeyError, ValueError) as error:
        raise ValueError(
            f"Basis {basis!r} is unavailable for atomic number {atomic_number}."
        ) from error

    orbital_count = 0
    for line in basis_string.splitlines():
        if not line.startswith("#BASIS SET:"):
            continue
        tokens = line.split("->")
        if len(tokens) != 2:
            continue
        orbitals = tokens[1][2:-1]
        for orbital in orbitals.split(","):
            orbital = orbital.strip()
            digit_count = 0
            for character in orbital:
                if character.isdigit():
                    digit_count += 1
                else:
                    break
            if digit_count == 0:
                continue
            count = int(orbital[:digit_count])
            shell = orbital[digit_count:]
            if shell in {"s", "p", "d", "f", "g"}:
                orbital_count += count * norb_for_shell(shell)

    if orbital_count <= 0:
        raise ValueError(
            f"Could not infer a spin-orbital count for basis={basis}, Z={atomic_number}."
        )
    return orbital_count


def molecular_counts(molecule: str, basis: str) -> tuple[int, int]:
    atom1, atom2 = MOLECULES[molecule]
    z1 = int(mendeleev.element(atom1).atomic_number)
    z2 = int(mendeleev.element(atom2).atomic_number)
    total_electrons = z1 + z2
    total_spin_orbitals = (
        atomic_spin_orbital_count(basis, z1)
        + atomic_spin_orbital_count(basis, z2)
    )
    return total_electrons, total_spin_orbitals


def active_space_sequence(
    total_electrons: int,
    total_spin_orbitals: int,
    max_active: int,
) -> list[tuple[int, int]]:
    """Reproduce the L-sweep 40/60 occupied/vacant active-space progression."""
    ratio_ideal = 0.4
    total_vacancies = total_spin_orbitals - total_electrons
    if total_electrons <= 0 or total_vacancies <= 0:
        raise ValueError(
            "The molecular basis must provide occupied and vacant spin orbitals."
        )

    n_act_occ = 2 - total_electrons % 2
    n_act_vac = n_act_occ
    spaces: list[tuple[int, int]] = []

    while n_act_occ + n_act_vac <= max_active:
        if n_act_occ <= total_electrons and n_act_vac <= total_vacancies:
            spaces.append((n_act_occ, n_act_vac))

        ratio_md = (n_act_occ + 1) / (n_act_occ + n_act_vac + 2)
        preferred = "vacant" if ratio_md >= ratio_ideal else "occupied"
        fallback = "occupied" if preferred == "vacant" else "vacant"

        changed = False
        for choice in (preferred, fallback):
            if choice == "vacant" and n_act_vac + 2 <= total_vacancies:
                n_act_vac += 2
                changed = True
                break
            if choice == "occupied" and n_act_occ + 2 <= total_electrons:
                n_act_occ += 2
                changed = True
                break
        if not changed:
            break

    return spaces


def write_configuration(
    library_root: Path,
    molecule: str,
    basis: str,
    bond_length: float,
    active_occupied: int,
    active_vacant: int,
) -> Path:
    atom1, atom2 = MOLECULES[molecule]
    bond_text = f"{bond_length:.2f}"
    stub = f"{molecule}_{bond_text}_{basis}"
    extended = (
        f"{stub}_as-{active_occupied:03d}-{active_vacant:03d}_{MAPPING}"
    )
    output_directory = library_root / molecule / bond_text / basis
    output_directory.mkdir(parents=True, exist_ok=True)

    absolute_stub = str(output_directory.resolve() / stub)
    config_path = output_directory / f"{extended}.config"
    half_length = 0.5 * bond_length

    lines = [
        "general.print_verbose()",
        (
            'general.logfile = '
            f'"{absolute_stub}_as-{active_occupied:03d}-'
            f'{active_vacant:03d}_{MAPPING}.log"'
        ),
        f'general.file_stub = "{absolute_stub}"',
        'general.file_format = "default"',
        f"L = {bond_length!r}",
        f'hamiltonian.add_atom("{atom1}", {-half_length!r}, 0.0, 0.0)',
        f'hamiltonian.add_atom("{atom2}", {half_length!r}, 0.0, 0.0)',
        f'hamiltonian.basis = "{basis}"',
        f"hamiltonian.num_active_occupied = {active_occupied}",
        f"hamiltonian.num_active_vacant = {active_vacant}",
        f'hamiltonian.f2q_mapping = "{MAPPING}"',
        "",
    ]
    config_path.write_text("\n".join(lines), encoding="utf-8")
    return config_path


def expected_config_path(
    library_root: Path,
    molecule: str,
    basis: str,
    bond_length: float,
    active_occupied: int,
    active_vacant: int,
) -> Path:
    bond_text = f"{bond_length:.2f}"
    stub = f"{molecule}_{bond_text}_{basis}"
    extended = (
        f"{stub}_as-{active_occupied:03d}-{active_vacant:03d}_{MAPPING}"
    )
    return library_root / molecule / bond_text / basis / f"{extended}.config"


def load_configuration_file(config_path: Path) -> State:
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
    return Path(state.filename_ham2()).with_suffix(".tensors.npz")


def save_tensors(tensor_path: Path, active_hamiltonian: Any) -> None:
    import numpy as np

    tensor_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        tensor_path,
        constant=active_hamiltonian.constant,
        one_body=active_hamiltonian.one_body_tensor,
        two_body=active_hamiltonian.two_body_tensor,
    )


def generate_case(
    *,
    args: argparse.Namespace,
    molecule: str,
    basis: str,
    bond_length: float,
    total_electrons: int,
    total_spin_orbitals: int,
    active_occupied: int,
    active_vacant: int,
) -> dict[str, Any]:
    n_qubits = active_occupied + active_vacant
    config_path = expected_config_path(
        args.library,
        molecule,
        basis,
        bond_length,
        active_occupied,
        active_vacant,
    )

    row: dict[str, Any] = {
        "molecule": molecule,
        "bond_length": f"{bond_length:.8f}",
        "basis": basis,
        "total_electrons": total_electrons,
        "total_spin_orbitals": total_spin_orbitals,
        "active_occupied": active_occupied,
        "active_vacant": active_vacant,
        "n_qubits": n_qubits,
        "config_path": str(config_path),
        "tensor_path": "",
        "config_status": "",
        "tensor_status": "",
        "error": "",
    }

    try:
        if config_path.exists() and not args.overwrite_configs:
            row["config_status"] = "existing"
        else:
            config_path = write_configuration(
                args.library,
                molecule,
                basis,
                bond_length,
                active_occupied,
                active_vacant,
            )
            row["config_path"] = str(config_path)
            row["config_status"] = "written"

        if args.configs_only:
            row["tensor_status"] = "not_requested"
            return row

        state = load_configuration_file(config_path)
        tensor_path = tensor_path_for_state(state)
        row["tensor_path"] = str(tensor_path)

        if tensor_path.exists() and not args.overwrite_tensors:
            row["tensor_status"] = "existing"
            return row

        active_hamiltonian = get_ham2(state)
        actual_n_qubits = int(active_hamiltonian.n_qubits)
        if actual_n_qubits != n_qubits:
            raise ValueError(
                "Generated tensor size does not match requested active space: "
                f"{actual_n_qubits} != {n_qubits}."
            )
        save_tensors(tensor_path, active_hamiltonian)
        row["tensor_status"] = "written"
    except Exception as error:
        row["tensor_status"] = "failed"
        row["error"] = f"{type(error).__name__}: {error}"

    return row


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    bases = [basis.lower() for basis in args.bases]
    requested_qubits = sorted(set(args.qubits))
    overrides = parse_bond_length_overrides(args.bond_length)
    bond_lengths = resolved_bond_lengths(overrides)

    configure_logging(level="verbose", logfile=None)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    print("C2/N2 active-space generation")
    print(f"Library: {args.library}")
    print(f"Molecules: {args.molecules}")
    print(f"Bases: {bases}")
    print(f"Requested qubits: {requested_qubits}")
    for molecule in args.molecules:
        print(f"{molecule} bond length: {bond_lengths[molecule]:.4f} Angstrom")

    rows: list[dict[str, Any]] = []

    for molecule in args.molecules:
        bond_length = bond_lengths[molecule]
        for basis in bases:
            total_electrons, total_spin_orbitals = molecular_counts(
                molecule,
                basis,
            )
            sequence = active_space_sequence(
                total_electrons,
                total_spin_orbitals,
                max(requested_qubits),
            )
            by_qubits = {
                occupied + vacant: (occupied, vacant)
                for occupied, vacant in sequence
            }

            missing_from_sequence = [
                q for q in requested_qubits if q not in by_qubits
            ]
            if missing_from_sequence:
                raise ValueError(
                    f"Requested active spaces are unavailable for {molecule}/{basis}: "
                    f"{missing_from_sequence}. Available: {sorted(by_qubits)}"
                )

            print()
            print("=" * 88)
            print(
                f"{molecule} / {basis}: electrons={total_electrons}, "
                f"spin_orbitals={total_spin_orbitals}"
            )
            print("=" * 88)

            for q in requested_qubits:
                active_occupied, active_vacant = by_qubits[q]
                print(
                    f"{q:2d}q -> active occupied={active_occupied:2d}, "
                    f"active vacant={active_vacant:2d}"
                )
                row = generate_case(
                    args=args,
                    molecule=molecule,
                    basis=basis,
                    bond_length=bond_length,
                    total_electrons=total_electrons,
                    total_spin_orbitals=total_spin_orbitals,
                    active_occupied=active_occupied,
                    active_vacant=active_vacant,
                )
                rows.append(row)
                if row["error"]:
                    print(f"  FAILED: {row['error']}")
                else:
                    print(
                        f"  config={row['config_status']}, "
                        f"tensor={row['tensor_status']}"
                    )

    with args.summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    failures = [row for row in rows if row["error"]]
    print()
    print("=" * 88)
    print(f"Cases considered: {len(rows)}")
    print(f"Failures: {len(failures)}")
    print(f"Summary: {args.summary}")
    print("=" * 88)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
