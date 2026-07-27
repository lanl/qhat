#!/usr/bin/env python3
"""Generate QHAT Hamiltonian-generator configs for small polyatomic molecules.

This is a companion to ``build_config_L_sweep.py``.  Instead of placing two
atoms at +/-L/2, it starts from a reference Cartesian geometry and uniformly
scales all coordinates by a dimensionless factor.

The built-in neutral molecules are:
    BeH2   3 atoms, linear
    H2O    3 atoms, bent
    NH3    4 atoms, trigonal pyramidal
    CH4    5 atoms, tetrahedral

Examples
--------
Generate the recommended first-pass STO-6G configs only::

    python build_config_polyatomic_sweep.py \
        --molecules BeH2,H2O,NH3,CH4 \
        --bases sto-6g \
        --scale-values 0.8,1.0,1.5 \
        --active-sizes 4,6,8 \
        --library polyatomic_library

Generate and immediately run both STO-6G and HGBS-5 configs::

    python build_config_polyatomic_sweep.py \
        --molecules BeH2,H2O,NH3,CH4 \
        --bases sto-6g,hgbs-5 \
        --scale-values 0.8,1.0,1.5 \
        --active-sizes 4,6,8 \
        --library polyatomic_library \
        --hamgen-dir . \
        --run

Notes
-----
* Active-space sizes are numbers of spin orbitals, hence numbers of qubits.
* The active occupied/vacant progression mirrors build_config_L_sweep.py:
  4 qubits -> (2 occupied, 2 vacant)
  6 qubits -> (2 occupied, 4 vacant)
  8 qubits -> (4 occupied, 4 vacant)
  whenever the molecule and basis have enough orbitals.
* QHAT's current hamgen.py constructs neutral molecules.  Therefore, the
  built-in cases here are all neutral; charged species are intentionally not
  included.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import basis_set_exchange
import mendeleev


# =================================================================================================
# Molecule definitions
# =================================================================================================

Atom = tuple[str, float, float, float]


@dataclass(frozen=True)
class MoleculeSpec:
    """Reference geometry and metadata for one neutral molecule."""

    name: str
    geometry: tuple[Atom, ...]
    description: str


def water_geometry() -> tuple[Atom, ...]:
    """H2O with O at the origin, r(O-H)=0.9584 A and angle H-O-H=104.45 deg."""
    r_oh = 0.9584
    half_angle = math.radians(104.45 / 2.0)
    x = r_oh * math.sin(half_angle)
    z = r_oh * math.cos(half_angle)
    return (
        ("O", 0.0, 0.0, 0.0),
        ("H", +x, 0.0, z),
        ("H", -x, 0.0, z),
    )


def ammonia_geometry() -> tuple[Atom, ...]:
    """NH3 with r(N-H)=1.0124 A and angle H-N-H=106.7 deg."""
    r_nh = 1.0124
    hnh_angle = math.radians(106.7)

    # Three equivalent H vectors separated by 120 degrees in azimuth.
    # For polar angle alpha, v_i dot v_j = 1.5*cos(alpha)^2 - 0.5.
    cos_alpha_squared = (math.cos(hnh_angle) + 0.5) / 1.5
    if cos_alpha_squared < 0.0:
        raise ValueError("Invalid NH3 reference angle")
    cos_alpha = math.sqrt(cos_alpha_squared)
    sin_alpha = math.sqrt(1.0 - cos_alpha_squared)
    radial_xy = r_nh * sin_alpha
    z = r_nh * cos_alpha

    atoms: list[Atom] = [("N", 0.0, 0.0, 0.0)]
    for azimuth_deg in (0.0, 120.0, 240.0):
        phi = math.radians(azimuth_deg)
        atoms.append(("H", radial_xy * math.cos(phi), radial_xy * math.sin(phi), z))
    return tuple(atoms)


def methane_geometry() -> tuple[Atom, ...]:
    """CH4 with C at the origin and r(C-H)=1.086 A."""
    r_ch = 1.086
    a = r_ch / math.sqrt(3.0)
    return (
        ("C", 0.0, 0.0, 0.0),
        ("H", +a, +a, +a),
        ("H", +a, -a, -a),
        ("H", -a, +a, -a),
        ("H", -a, -a, +a),
    )


MOLECULES: dict[str, MoleculeSpec] = {
    "BEH2": MoleculeSpec(
        name="BeH2",
        geometry=(
            ("Be", 0.0, 0.0, 0.0),
            ("H", -1.3264, 0.0, 0.0),
            ("H", +1.3264, 0.0, 0.0),
        ),
        description="linear BeH2",
    ),
    "H2O": MoleculeSpec(
        name="H2O",
        geometry=water_geometry(),
        description="bent H2O",
    ),
    "NH3": MoleculeSpec(
        name="NH3",
        geometry=ammonia_geometry(),
        description="trigonal-pyramidal NH3",
    ),
    "CH4": MoleculeSpec(
        name="CH4",
        geometry=methane_geometry(),
        description="tetrahedral CH4",
    ),
}


# =================================================================================================
# Basis and orbital-count helpers
# =================================================================================================


def spin_orbitals_for_shell(shell: str) -> int:
    """Return the number of spin orbitals represented by one shell."""
    counts = {"s": 2, "p": 6, "d": 10, "f": 14, "g": 18}
    try:
        return counts[shell]
    except KeyError as exc:
        raise ValueError(f"Unsupported shell label: {shell!r}") from exc


def get_atomic_spin_orbital_count(basis: str, element: str) -> int | None:
    """Mirror the orbital-count logic used by build_config_L_sweep.py."""
    atomic_number = mendeleev.element(element).atomic_number
    try:
        basis_string = basis_set_exchange.get_basis(
            basis,
            elements=atomic_number,
            fmt="nwchem",
        )
    except (KeyError, ValueError):
        return None

    orbital_count = 0
    for line in basis_string.splitlines():
        if not line.startswith("#BASIS SET:"):
            continue

        tokens = line.split("->")
        if len(tokens) != 2:
            continue
        orbitals = tokens[1].strip().strip("[]")

        for orbital in orbitals.split(","):
            orbital = orbital.strip()
            digit_count = 0
            for char in orbital:
                if char.isdigit():
                    digit_count += 1
                else:
                    break
            if digit_count == 0:
                continue

            count = int(orbital[:digit_count])
            shell = orbital[digit_count:].strip().lower()
            if shell not in {"s", "p", "d", "f", "g"}:
                continue
            orbital_count += count * spin_orbitals_for_shell(shell)

    return orbital_count if orbital_count > 0 else None


def molecule_electron_count(spec: MoleculeSpec) -> int:
    """Total electrons for a neutral molecule."""
    return sum(mendeleev.element(atom[0]).atomic_number for atom in spec.geometry)


def molecule_spin_orbital_count(spec: MoleculeSpec, basis: str) -> int | None:
    """Total spin-orbital count for the molecule and basis."""
    total = 0
    for element, *_ in spec.geometry:
        count = get_atomic_spin_orbital_count(basis, element)
        if count is None:
            return None
        total += count
    return total


# =================================================================================================
# Active-space construction
# =================================================================================================


def active_space_sequence(
    total_electrons: int,
    total_spin_orbitals: int,
    max_active: int,
) -> list[tuple[int, int]]:
    """Build the same occupied/vacant progression as build_config_L_sweep.py."""
    ratio_ideal = 0.4
    total_vacancies = total_spin_orbitals - total_electrons

    if total_electrons <= 0:
        raise ValueError("The molecule must have at least one electron")
    if total_vacancies <= 0:
        raise ValueError("The basis must provide at least one vacant spin orbital")

    n_act_occ = 2 - total_electrons % 2
    n_act_vac = n_act_occ
    spaces: list[tuple[int, int]] = []

    while n_act_occ + n_act_vac <= max_active:
        if n_act_occ <= total_electrons and n_act_vac <= total_vacancies:
            spaces.append((n_act_occ, n_act_vac))

        ratio_md = (n_act_occ + 1) / (n_act_occ + n_act_vac + 2)

        if ratio_md >= ratio_ideal:
            preferred = "vacant"
        else:
            preferred = "occupied"

        changed = False
        for choice in (preferred, "occupied" if preferred == "vacant" else "vacant"):
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


# =================================================================================================
# Config writing
# =================================================================================================


def scaled_geometry(spec: MoleculeSpec, scale: float) -> tuple[Atom, ...]:
    if scale <= 0.0:
        raise ValueError("Geometry scale factors must be positive")
    return tuple(
        (element, scale * x, scale * y, scale * z)
        for element, x, y, z in spec.geometry
    )


def write_config(
    spec: MoleculeSpec,
    basis: str,
    scale: float,
    occupied: int,
    vacant: int,
    mapping: str,
    library_root: str | Path,
) -> Path:
    """Write one QHAT .config file and return its path."""
    scale_tag = f"s-{scale:.2f}"
    stub = f"{spec.name}_{scale_tag}_{basis}"
    extended = f"{stub}_as-{occupied:03d}-{vacant:03d}_{mapping}"

    path = Path(library_root) / spec.name / scale_tag / basis
    path.mkdir(parents=True, exist_ok=True)

    # Deliberately omit active-space and mapping tags from file_stub.  QHAT's
    # State class appends those tags to active-space and qubit-Hamiltonian files,
    # while all active spaces can reuse the same Hartree-Fock pickle.
    abs_stub = str((path.resolve() / stub))
    filename = path / f"{extended}.config"

    geometry = scaled_geometry(spec, scale)

    with filename.open("w", encoding="utf-8") as fout:
        print(f"# Molecule: {spec.name} ({spec.description})", file=fout)
        print(f"# Uniform geometry scale relative to the built-in reference: {scale}", file=fout)
        print("general.print_verbose()", file=fout)
        print(f'general.logfile = "{abs_stub}_as-{occupied:03d}-{vacant:03d}_{mapping}.log"', file=fout)
        print(f'general.file_stub = "{abs_stub}"', file=fout)
        print('general.file_format = "default"', file=fout)

        for element, x, y, z in geometry:
            print(
                f'hamiltonian.add_atom("{element}", '
                f"{x:.12f}, {y:.12f}, {z:.12f})",
                file=fout,
            )

        print(f'hamiltonian.basis = "{basis}"', file=fout)
        print(f"hamiltonian.num_active_occupied = {occupied}", file=fout)
        print(f"hamiltonian.num_active_vacant = {vacant}", file=fout)
        print(f'hamiltonian.f2q_mapping = "{mapping}"', file=fout)

    return filename


# =================================================================================================
# Generation and execution
# =================================================================================================


def parse_csv_strings(value: str) -> list[str]:
    values = [item.strip() for item in value.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one comma-separated value")
    return values


def parse_csv_floats(value: str) -> list[float]:
    try:
        values = [float(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Scale values must be numbers") from exc
    if not values or any(item <= 0.0 for item in values):
        raise argparse.ArgumentTypeError("Scale values must be positive")
    return values


def parse_csv_even_ints(value: str) -> list[int]:
    try:
        values = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Active-space sizes must be integers") from exc
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one active-space size")
    if any(item < 2 or item % 2 != 0 for item in values):
        raise argparse.ArgumentTypeError("Active-space sizes must be positive even integers")
    return sorted(set(values))


def resolve_molecules(names: Sequence[str]) -> list[MoleculeSpec]:
    result: list[MoleculeSpec] = []
    for name in names:
        key = name.replace("-", "").upper()
        try:
            result.append(MOLECULES[key])
        except KeyError as exc:
            available = ", ".join(spec.name for spec in MOLECULES.values())
            raise ValueError(f"Unknown molecule {name!r}. Available: {available}") from exc
    return result


def generate_configs(
    molecules: Sequence[MoleculeSpec],
    bases: Sequence[str],
    scale_values: Sequence[float],
    requested_active_sizes: Sequence[int],
    mappings: Sequence[str],
    library_root: str | Path,
) -> list[Path]:
    config_files: list[Path] = []
    max_active = max(requested_active_sizes)
    requested = set(requested_active_sizes)

    print("Generating polyatomic QHAT configurations")
    print("=" * 80)

    for spec in molecules:
        total_electrons = molecule_electron_count(spec)
        print(f"\n{spec.name}: {spec.description}; electrons={total_electrons}")

        for basis in bases:
            total_spin_orbitals = molecule_spin_orbital_count(spec, basis)
            if total_spin_orbitals is None:
                print(f"  SKIP basis={basis}: basis unavailable for at least one element")
                continue

            spaces = [
                pair
                for pair in active_space_sequence(
                    total_electrons,
                    total_spin_orbitals,
                    max_active,
                )
                if sum(pair) in requested
            ]

            missing = sorted(requested - {sum(pair) for pair in spaces})
            print(
                f"  basis={basis:7s} total_spin_orbitals={total_spin_orbitals:3d} "
                f"active_spaces={spaces}"
            )
            if missing:
                print(f"    WARNING: unavailable requested active sizes: {missing}")

            for scale in scale_values:
                for occupied, vacant in spaces:
                    for mapping in mappings:
                        filename = write_config(
                            spec=spec,
                            basis=basis,
                            scale=scale,
                            occupied=occupied,
                            vacant=vacant,
                            mapping=mapping,
                            library_root=library_root,
                        )
                        config_files.append(filename)
                        print(
                            f"    {spec.name:5s} scale={scale:4.2f} {basis:7s} "
                            f"as={occupied:03d}-{vacant:03d} {mapping:2s} -> {filename}"
                        )

    print("\n" + "=" * 80)
    print(f"Total configs written: {len(config_files)}")
    return config_files


def run_hamgen(config_files: Iterable[Path], hamgen_dir: str | Path) -> None:
    config_files = list(config_files)
    hamgen_path = Path(hamgen_dir) / "hamgen.py"
    if not hamgen_path.exists():
        raise FileNotFoundError(f"hamgen.py not found at {hamgen_path.resolve()}")

    print(f"\nRunning hamgen.py on {len(config_files)} configs...")
    failures: list[tuple[Path, int]] = []

    for index, config_file in enumerate(config_files, start=1):
        print(f"[{index}/{len(config_files)}] {config_file}")
        result = subprocess.run(
            [sys.executable, str(hamgen_path), str(config_file)],
            check=False,
        )
        if result.returncode != 0:
            failures.append((config_file, result.returncode))
            print(f"  WARNING: exited with code {result.returncode}")

    if failures:
        print("\nHamiltonian generation completed with failures:")
        for config_file, return_code in failures:
            print(f"  code={return_code:3d}  {config_file}")
        raise SystemExit(1)

    print("\nHamiltonian generation completed successfully.")


# =================================================================================================
# Command line
# =================================================================================================


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate QHAT Hamiltonian-generator configs for BeH2, H2O, NH3, "
            "and CH4 by uniformly scaling reference Cartesian geometries."
        )
    )
    parser.add_argument(
        "--molecules",
        type=parse_csv_strings,
        default=parse_csv_strings("BeH2,H2O,NH3,CH4"),
        help="Comma-separated molecules: BeH2,H2O,NH3,CH4",
    )
    parser.add_argument(
        "--bases",
        type=parse_csv_strings,
        default=parse_csv_strings("sto-6g"),
        help="Comma-separated basis sets (default: sto-6g)",
    )
    parser.add_argument(
        "--scale-values",
        type=parse_csv_floats,
        default=parse_csv_floats("0.8,1.0,1.5"),
        help="Comma-separated uniform geometry scales (default: 0.8,1.0,1.5)",
    )
    parser.add_argument(
        "--active-sizes",
        type=parse_csv_even_ints,
        default=parse_csv_even_ints("4,6,8"),
        help="Comma-separated active-space sizes/qubits (default: 4,6,8)",
    )
    parser.add_argument(
        "--mappings",
        type=parse_csv_strings,
        default=parse_csv_strings("JW"),
        help="Comma-separated mappings (default: JW)",
    )
    parser.add_argument(
        "--library",
        default="polyatomic_library",
        help="Output library root (default: polyatomic_library)",
    )
    parser.add_argument(
        "--hamgen-dir",
        default=".",
        help="Directory containing hamgen.py (default: current directory)",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run hamgen.py after writing all configs",
    )
    parser.add_argument(
        "--list-molecules",
        action="store_true",
        help="List built-in molecules and exit",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.list_molecules:
        for spec in MOLECULES.values():
            print(f"{spec.name:5s}  {len(spec.geometry)} atoms  {spec.description}")
        return

    try:
        molecules = resolve_molecules(args.molecules)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    bases = [basis.lower() for basis in args.bases]
    mappings = [mapping.upper() for mapping in args.mappings]
    unsupported_mappings = sorted(set(mappings) - {"JW", "BK"})
    if unsupported_mappings:
        raise SystemExit(f"Unsupported mappings: {unsupported_mappings}; use JW and/or BK")

    config_files = generate_configs(
        molecules=molecules,
        bases=bases,
        scale_values=args.scale_values,
        requested_active_sizes=args.active_sizes,
        mappings=mappings,
        library_root=args.library,
    )

    if not config_files:
        raise SystemExit("No configs were generated. Check molecule and basis availability.")

    if args.run:
        run_hamgen(config_files, args.hamgen_dir)
    else:
        print("\nConfigs only were generated. Add --run to execute hamgen.py.")


if __name__ == "__main__":
    main()