#!/usr/bin/env python3
"""Extend B-H, Be-H, and Li-H active-space tensors on the QHAT L-sweep branch.

The script preserves the bond lengths used by the existing heteronuclear
4/6/8-qubit study, follows the same 40/60 active occupied/vacant progression
as build_config_L_sweep.py, and generates every requested active space that is
available in the chosen basis. Unsupported sizes are reported and skipped
instead of aborting the whole multi-molecule run.

Run from the QHAT repository root after `python -m pip install -e .`.
"""

from __future__ import annotations

import argparse
import csv
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
    "B-H": ("B", "H"),
    "Be-H": ("Be", "H"),
    "Li-H": ("Li", "H"),
}

# Exact L values stored in the existing L-sweep configs.  The directory names
# are rounded to two decimals: 1.29, 1.47, and 1.81, respectively.
DEFAULT_BOND_LENGTHS: dict[str, float] = {
    "B-H": 1.287,
    "Be-H": 1.474,
    "Li-H": 1.815,
}

DEFAULT_BASES = ("sto-6g", "hgbs-5")
DEFAULT_QUBITS = tuple(range(10, 21, 2))
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
        default=Path("hamiltonian_generator/library"),
        help="Library root. Default: hamiltonian_generator/library.",
    )
    parser.add_argument(
        "--molecules",
        nargs="+",
        default=list(MOLECULES),
        choices=sorted(MOLECULES),
        help="Molecules to generate. Default: B-H Be-H Li-H.",
    )
    parser.add_argument(
        "--bases",
        nargs="+",
        default=list(DEFAULT_BASES),
        help="Basis sets. Default: sto-6g hgbs-5.",
    )
    parser.add_argument(
        "--qubits",
        nargs="+",
        type=int,
        default=list(DEFAULT_QUBITS),
        help="Requested even active-space sizes. Default: 10 12 ... 20.",
    )
    parser.add_argument(
        "--bond-length",
        action="append",
        default=[],
        metavar="MOLECULE=VALUE",
        help=(
            "Override a fixed bond length in Angstrom. May be repeated. "
            "Defaults preserve the existing 4/6/8-qubit heteronuclear study."
        ),
    )
    parser.add_argument(
        "--configs-only",
        action="store_true",
        help="Write/verify configs but do not generate tensors.",
    )
    parser.add_argument(
        "--overwrite-configs",
        action="store_true",
        help="Rewrite configs that already exist.",
    )
    parser.add_argument(
        "--overwrite-tensors",
        action="store_true",
        help="Regenerate tensors that already exist.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("analysis/heteronuclear_extended_tensor_summary.csv"),
        help="Generation summary CSV.",
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
                "--bond-length must use MOLECULE=VALUE, e.g. B-H=1.287."
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


def resolved_bond_lengths(overrides: dict[str, float]) -> dict[str, float]:
    return {
        molecule: overrides.get(molecule, DEFAULT_BOND_LENGTHS[molecule])
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
            f"Could not infer spin-orbital count for basis={basis}, "
            f"Z={atomic_number}."
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
    """Reproduce the L-sweep 40/60 occupied/vacant progression."""
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
    extended = f"{stub}_as-{active_occupied:03d}-{active_vacant:03d}_{MAPPING}"
    return library_root / molecule / bond_text / basis / f"{extended}.config"


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
    extended = f"{stub}_as-{active_occupied:03d}-{active_vacant:03d}_{MAPPING}"
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


def unsupported_row(
    *,
    molecule: str,
    basis: str,
    bond_length: float,
    total_electrons: int,
    total_spin_orbitals: int,
    n_qubits: int,
    available: Sequence[int],
) -> dict[str, Any]:
    return {
        "molecule": molecule,
        "bond_length": f"{bond_length:.8f}",
        "basis": basis,
        "total_electrons": total_electrons,
        "total_spin_orbitals": total_spin_orbitals,
        "active_occupied": "",
        "active_vacant": "",
        "n_qubits": n_qubits,
        "config_path": "",
        "tensor_path": "",
        "config_status": "unsupported",
        "tensor_status": "unsupported",
        "error": f"Unavailable under 40/60 progression; available sizes: {list(available)}",
    }


def write_summary(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_arguments()
    validate_arguments(args)
    bases = [basis.lower() for basis in args.bases]
    requested_qubits = sorted(set(args.qubits))
    overrides = parse_bond_length_overrides(args.bond_length)
    bond_lengths = resolved_bond_lengths(overrides)

    configure_logging(level="verbose", logfile=None)
    print("Heteronuclear active-space extension")
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
            total_electrons, total_spin_orbitals = molecular_counts(molecule, basis)
            sequence = active_space_sequence(
                total_electrons,
                total_spin_orbitals,
                max(requested_qubits),
            )
            by_qubits = {
                occupied + vacant: (occupied, vacant)
                for occupied, vacant in sequence
            }
            available = sorted(by_qubits)

            print()
            print("=" * 88)
            print(
                f"{molecule} / {basis}: electrons={total_electrons}, "
                f"spin_orbitals={total_spin_orbitals}, available={available}"
            )
            print("=" * 88)

            for q in requested_qubits:
                if q not in by_qubits:
                    print(f"{q:2d}q -> SKIP (not available in this basis)")
                    rows.append(
                        unsupported_row(
                            molecule=molecule,
                            basis=basis,
                            bond_length=bond_length,
                            total_electrons=total_electrons,
                            total_spin_orbitals=total_spin_orbitals,
                            n_qubits=q,
                            available=available,
                        )
                    )
                    continue

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
                print(
                    f"      config={row['config_status']} "
                    f"tensor={row['tensor_status']}"
                )
                if row["error"]:
                    print(f"      ERROR: {row['error']}")

    write_summary(args.summary, rows)
    failed = sum(row["tensor_status"] == "failed" for row in rows)
    unsupported = sum(row["tensor_status"] == "unsupported" for row in rows)
    ready = sum(row["tensor_status"] in {"written", "existing"} for row in rows)
    print()
    print(f"Summary: {args.summary}")
    print(f"Ready tensor cases: {ready}")
    print(f"Unsupported requested cases: {unsupported}")
    print(f"Failed cases: {failed}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

