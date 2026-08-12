#!/usr/bin/env python3
"""Generate and benchmark an O2 triplet active-space ordering study on QHAT L-sweep.

Run this file from the QHAT repository root after installing QHAT editable:

    python -m pip install -e .
    python analysis/run_o2_fermionic_ordering_study.py smoke
    python analysis/run_o2_fermionic_ordering_study.py full

The full mode covers STO-6G and HGBS-5, even active spaces 4..20 qubits,
and the five deterministic orderings already defined by
analysis/benchmark_b2_signed_coefficient_baseline.py:

    jw_raw
    signed_coefficient_lexicographic
    jw_magnitude_descending_lexicographic
    fermionic_signed_coefficient_lexicographic
    fermionic_magnitude_descending_lexicographic

O2 requires special handling because its X-state reference is triplet.  The
current L-sweep Hamiltonian generator infers multiplicity only from electron
parity, which would make neutral even-electron O2 a singlet.  This driver
therefore performs an explicit multiplicity-3 ROHF calculation, applies QHAT's
existing active-space machinery, records the actual ROHF active determinant,
and uses that determinant for exact-sector and Trotter comparisons.

Resume behavior is ON by default:
  * existing configs are preserved;
  * existing triplet-aware tensors are preserved;
  * successful (case_id, ordering, steps, time) result rows are skipped;
  * each completed ordering row is flushed immediately.

"Full" is exhaustive over the requested basis/active-space/order-method grid;
it is not an exhaustive factorial permutation of every Pauli term.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import pickle
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

import basis_set_exchange
import mendeleev
import numpy as np
from scipy.sparse.linalg import expm_multiply

from openfermion import FermionOperator, MolecularData, get_fermion_operator, jordan_wigner
try:
    from openfermion import get_number_preserving_sparse_operator
except ImportError:  # Compatibility with older OpenFermion exports.
    from openfermion.linalg import get_number_preserving_sparse_operator
from openfermionpyscf import PyscfMolecularData
from openfermionpyscf._run_pyscf import compute_integrals, compute_scf, prepare_pyscf_molecule

try:
    from qhat.common.logging_utils import configure_logging
    from qhat.hamiltonian_generator.hamgen import apply_active_space
    from qhat.hamiltonian_generator.hamgen_types import (
        GeneralConfigurationUser,
        HamiltonianConfiguration,
        State,
    )
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as matrixfree
    from qhat.analysis import benchmark_b2_coloring_robustness as robustness
    from qhat.analysis.benchmark_L_sweep_trotter import (
        DEFAULT_TOLERANCE,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
    )
except ImportError:
    # Fallback when this file is placed directly in analysis/ and QHAT's package
    # path has not been installed in editable mode.
    from benchmark_b2_signed_coefficient_baseline import (  # type: ignore
        ORDERING_NAMES as _ORDERING_NAMES,
    )
    import benchmark_b2_signed_coefficient_baseline as baseline  # type: ignore
    import benchmark_b2_active_spaces_matrix_free as matrixfree  # type: ignore
    import benchmark_b2_coloring_robustness as robustness  # type: ignore
    from benchmark_L_sweep_trotter import (  # type: ignore
        DEFAULT_TOLERANCE,
        build_pauli_noncommutation_graph,
        clean_fermion_operator,
        load_interaction_operator,
        parse_case_metadata,
        real_coefficient,
    )
    from qhat.common.logging_utils import configure_logging
    from qhat.hamiltonian_generator.hamgen import apply_active_space
    from qhat.hamiltonian_generator.hamgen_types import (
        GeneralConfigurationUser,
        HamiltonianConfiguration,
        State,
    )


MOLECULE = "O-O"
ELEMENTS = ("O", "O")
MULTIPLICITY = 3
CHARGE = 0
MAPPING = "JW"
DEFAULT_BASES = ("sto-6g", "hgbs-5")
DEFAULT_QUBITS = tuple(range(4, 21, 2))

EXTRA_FIELDS = [
    "multiplicity",
    "charge",
    "active_alpha_electrons",
    "active_beta_electrons",
    "reference_determinant_bitstring",
]
RESULT_FIELDS = list(baseline.FIELDNAMES) + EXTRA_FIELDS

GENERATION_FIELDS = [
    "molecule",
    "bond_length",
    "basis",
    "multiplicity",
    "charge",
    "total_electrons",
    "total_spin_orbitals",
    "active_occupied",
    "active_vacant",
    "n_qubits",
    "active_alpha_electrons",
    "active_beta_electrons",
    "reference_determinant_bitstring",
    "config_path",
    "ham1_pickle_path",
    "ham2_pickle_path",
    "tensor_path",
    "config_status",
    "tensor_status",
    "error",
]

SUMMARY_FIELDS = [
    "case_id",
    "molecule",
    "bond_length",
    "basis",
    "n_qubits",
    "active_occupied",
    "active_vacant",
    "active_alpha_electrons",
    "active_beta_electrons",
    "jw_raw_state_infidelity",
    "signed_coefficient_lexicographic_ratio_to_raw",
    "jw_magnitude_descending_lexicographic_ratio_to_raw",
    "fermionic_signed_coefficient_lexicographic_ratio_to_raw",
    "fermionic_magnitude_descending_lexicographic_ratio_to_raw",
    "best_ordering",
    "best_state_infidelity_ratio_to_raw",
]


# -----------------------------------------------------------------------------
# CLI and active-space configuration
# -----------------------------------------------------------------------------


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        choices=("smoke", "full"),
        default="full",
        help="smoke = STO-6G/4q only; full = STO-6G + HGBS-5, 4..20q.",
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=Path("hamiltonian_generator/o2_active_space_library"),
        help="Output library root for O2 configs, pickles, and tensors.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/o2_fermionic_ordering_study"),
        help="Directory for benchmark CSVs and summaries.",
    )
    parser.add_argument(
        "--bases",
        nargs="+",
        default=None,
        help="Override basis list. Example: --bases sto-6g hgbs-5",
    )
    parser.add_argument(
        "--qubits",
        nargs="+",
        type=int,
        default=None,
        help="Override active-space qubits. Example: --qubits 4 6 8 10",
    )
    parser.add_argument(
        "--bond-length",
        type=float,
        default=None,
        help=(
            "O-O separation in Angstrom. Default: sum of two Pyykko covalent "
            "radii, matching the L-sweep diatomic convention."
        ),
    )
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--time", type=float, default=1.0, dest="evolution_time")
    parser.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    parser.add_argument("--parallel-threshold", type=int, default=2**16)
    parser.add_argument(
        "--no-spin-sector",
        action="store_true",
        help="Use fixed-particle exact sector only instead of fixed particle + Sz.",
    )
    parser.add_argument(
        "--overwrite-configs",
        action="store_true",
        help="Rewrite existing O2 config files.",
    )
    parser.add_argument(
        "--overwrite-tensors",
        action="store_true",
        help="Regenerate O2 ROHF active-space tensors and reference metadata.",
    )
    parser.add_argument(
        "--fresh-results",
        action="store_true",
        help="Overwrite benchmark CSV instead of resuming successful rows.",
    )
    return parser.parse_args()


def resolved_requested_grid(args: argparse.Namespace) -> tuple[list[str], list[int]]:
    if args.bases is not None:
        bases = [str(value).lower() for value in args.bases]
    elif args.mode == "smoke":
        bases = ["sto-6g"]
    else:
        bases = list(DEFAULT_BASES)

    if args.qubits is not None:
        qubits = sorted(set(int(value) for value in args.qubits))
    elif args.mode == "smoke":
        qubits = [4]
    else:
        qubits = list(DEFAULT_QUBITS)

    if not bases:
        raise ValueError("At least one basis is required.")
    if not qubits or any(q <= 0 or q % 2 for q in qubits):
        raise ValueError("--qubits must contain positive even integers.")
    if args.steps <= 0:
        raise ValueError("--steps must be positive.")
    if args.evolution_time <= 0.0:
        raise ValueError("--time must be positive.")
    if args.tolerance <= 0.0:
        raise ValueError("--tolerance must be positive.")
    return bases, qubits


def pyykko_o2_bond_length() -> float:
    radius_pm = mendeleev.element("O").covalent_radius_pyykko
    if radius_pm is None:
        raise ValueError("No Pyykko covalent radius for O; pass --bond-length.")
    return round(2.0 * 0.01 * float(radius_pm), 4)


def norb_for_shell(shell: str) -> int:
    return {"s": 2, "p": 6, "d": 10, "f": 14, "g": 18}[shell]


def atomic_spin_orbital_count(basis: str, atomic_number: int) -> int:
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
            f"Could not infer spin-orbital count for basis={basis}, Z={atomic_number}."
        )
    return orbital_count


def molecular_counts(basis: str) -> tuple[int, int]:
    z = int(mendeleev.element("O").atomic_number)
    total_electrons = 2 * z - CHARGE
    total_spin_orbitals = 2 * atomic_spin_orbital_count(basis, z)
    return total_electrons, total_spin_orbitals


def active_space_sequence(
    total_electrons: int,
    total_spin_orbitals: int,
    max_active: int,
) -> list[tuple[int, int]]:
    """Copy the L-sweep 40/60 occupied/vacant progression."""
    ratio_ideal = 0.4
    total_vacancies = total_spin_orbitals - total_electrons
    if total_electrons <= 0 or total_vacancies <= 0:
        raise ValueError("Basis must provide occupied and vacant spin orbitals.")

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


# -----------------------------------------------------------------------------
# Configs and explicit triplet ROHF generation
# -----------------------------------------------------------------------------


def case_paths(
    library: Path,
    basis: str,
    bond_length: float,
    active_occupied: int,
    active_vacant: int,
) -> dict[str, Path]:
    bond_text = f"{bond_length:.2f}"
    stub_name = f"{MOLECULE}_{bond_text}_{basis}"
    directory = library / MOLECULE / bond_text / basis
    extended = f"{stub_name}_as-{active_occupied:03d}-{active_vacant:03d}_{MAPPING}"
    file_stub = directory / stub_name
    ham2_pickle = Path(
        f"{file_stub}_as-{active_occupied:03d}-{active_vacant:03d}.pickle"
    )
    return {
        "directory": directory,
        "file_stub": file_stub,
        "config": directory / f"{extended}.config",
        "ham1_pickle": Path(f"{file_stub}.pickle"),
        "ham2_pickle": ham2_pickle,
        "tensor": ham2_pickle.with_suffix(".tensors.npz"),
    }


def write_configuration(
    paths: dict[str, Path],
    basis: str,
    bond_length: float,
    active_occupied: int,
    active_vacant: int,
) -> None:
    paths["directory"].mkdir(parents=True, exist_ok=True)
    half_length = 0.5 * bond_length
    absolute_stub = str(paths["file_stub"].resolve())
    lines = [
        "# O2 triplet configuration generated by run_o2_fermionic_ordering_study.py",
        "# NOTE: current L-sweep hamgen.py does not consume multiplicity/charge;",
        "# this study's driver performs the explicit multiplicity-3 ROHF calculation.",
        "general.print_verbose()",
        f'general.logfile = "{absolute_stub}_as-{active_occupied:03d}-{active_vacant:03d}_{MAPPING}.log"',
        f'general.file_stub = "{absolute_stub}"',
        'general.file_format = "default"',
        f"L = {bond_length!r}",
        f'hamiltonian.add_atom("O", {-half_length!r}, 0.0, 0.0)',
        f'hamiltonian.add_atom("O", {half_length!r}, 0.0, 0.0)',
        f'hamiltonian.basis = "{basis}"',
        f"hamiltonian.num_active_occupied = {active_occupied}",
        f"hamiltonian.num_active_vacant = {active_vacant}",
        f'hamiltonian.f2q_mapping = "{MAPPING}"',
        f"hamiltonian.multiplicity = {MULTIPLICITY}",
        f"hamiltonian.charge = {CHARGE}",
        "",
    ]
    paths["config"].write_text("\n".join(lines), encoding="utf-8")


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


def run_o2_rohf(basis: str, bond_length: float) -> tuple[PyscfMolecularData, np.ndarray]:
    """Run neutral O2 multiplicity-3 ROHF and populate OpenFermion integrals."""
    geometry = [
        ("O", (-0.5 * bond_length, 0.0, 0.0)),
        ("O", (0.5 * bond_length, 0.0, 0.0)),
    ]
    molecule = MolecularData(
        geometry=geometry,
        basis=basis,
        multiplicity=MULTIPLICITY,
        charge=CHARGE,
        description=f"O2_triplet_{bond_length:.2f}",
    )
    pyscf_molecule = prepare_pyscf_molecule(molecule)
    molecule.n_orbitals = int(pyscf_molecule.nao_nr())
    molecule.n_qubits = 2 * molecule.n_orbitals
    molecule.nuclear_repulsion = float(pyscf_molecule.energy_nuc())

    print(f"  Running O2 multiplicity-{MULTIPLICITY} ROHF for {basis} ...", flush=True)
    start = time.perf_counter()
    pyscf_scf = compute_scf(pyscf_molecule)
    pyscf_scf.verbose = 0
    pyscf_scf.run()
    if not bool(getattr(pyscf_scf, "converged", True)):
        raise RuntimeError(f"PySCF ROHF did not converge for O2/{basis}.")
    molecule.hf_energy = float(pyscf_scf.e_tot)
    molecule._pyscf_data = {"mol": pyscf_molecule, "scf": pyscf_scf}
    molecule.canonical_orbitals = pyscf_scf.mo_coeff.astype(float)
    molecule.orbital_energies = pyscf_scf.mo_energy.astype(float)
    one_body, two_body = compute_integrals(pyscf_molecule, pyscf_scf)
    molecule.one_body_integrals = one_body
    molecule.two_body_integrals = two_body
    molecule.overlap_integrals = pyscf_scf.get_ovlp()
    molecule.basis = basis
    molecule.separation = bond_length
    molecule.hf_time = time.perf_counter() - start

    pyscf_data = PyscfMolecularData.__new__(PyscfMolecularData)
    pyscf_data.__dict__.update(molecule.__dict__)
    mo_occ = np.asarray(pyscf_scf.mo_occ, dtype=float)

    singly = int(np.count_nonzero(np.isclose(mo_occ, 1.0, atol=1e-7)))
    if singly != 2:
        raise ValueError(
            f"Expected two singly occupied ROHF orbitals for triplet O2, found {singly}: "
            f"mo_occ={mo_occ.tolist()}"
        )
    print(
        f"  ROHF converged: E_HF={pyscf_data.hf_energy:.12f}, "
        f"spatial orbitals={pyscf_data.n_orbitals}, singly occupied={singly}",
        flush=True,
    )
    return pyscf_data, mo_occ


def active_reference_from_rohf(
    mo_occ: np.ndarray,
    total_electrons: int,
    active_occupied: int,
    active_vacant: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return active spatial occupations and alpha/beta spin-orbital determinant."""
    n_frozen_occupied_spin = total_electrons - active_occupied
    if n_frozen_occupied_spin < 0 or n_frozen_occupied_spin % 2:
        raise ValueError(
            "O2 active space must freeze an even number of occupied spin orbitals."
        )
    first_active_spatial = n_frozen_occupied_spin // 2
    active_spatial_count = (active_occupied + active_vacant) // 2
    last_active_spatial = first_active_spatial + active_spatial_count
    active_mo_occ = np.asarray(
        mo_occ[first_active_spatial:last_active_spatial],
        dtype=float,
    )
    if active_mo_occ.size != active_spatial_count:
        raise ValueError(
            f"Active ROHF occupation slice has {active_mo_occ.size} orbitals, "
            f"expected {active_spatial_count}."
        )

    reference: list[bool] = []
    for occupation in active_mo_occ:
        if np.isclose(occupation, 2.0, atol=1e-7):
            reference.extend((True, True))
        elif np.isclose(occupation, 1.0, atol=1e-7):
            # ROHF triplet convention: singly occupied orbitals carry alpha spin.
            reference.extend((True, False))
        elif np.isclose(occupation, 0.0, atol=1e-7):
            reference.extend((False, False))
        else:
            raise ValueError(
                "Unexpected ROHF orbital occupation in active space: "
                f"{occupation!r}; full active mo_occ={active_mo_occ.tolist()}"
            )

    determinant = np.asarray(reference, dtype=bool)
    if determinant.size != active_occupied + active_vacant:
        raise RuntimeError("Reference determinant has the wrong qubit count.")
    if int(np.sum(determinant)) != active_occupied:
        raise ValueError(
            f"Reference has {int(np.sum(determinant))} electrons but active_occupied="
            f"{active_occupied}. active_mo_occ={active_mo_occ.tolist()}"
        )
    alpha = int(np.sum(determinant[::2]))
    beta = int(np.sum(determinant[1::2]))
    if alpha - beta != MULTIPLICITY - 1:
        raise ValueError(
            f"Active O2 determinant has Nalpha-Nbeta={alpha-beta}, expected "
            f"{MULTIPLICITY-1} for multiplicity {MULTIPLICITY}."
        )
    return active_mo_occ, determinant


def save_triplet_tensors(
    tensor_path: Path,
    active_hamiltonian: Any,
    reference: np.ndarray,
    active_mo_occ: np.ndarray,
) -> None:
    tensor_path.parent.mkdir(parents=True, exist_ok=True)
    alpha = int(np.sum(reference[::2]))
    beta = int(np.sum(reference[1::2]))
    np.savez_compressed(
        tensor_path,
        constant=active_hamiltonian.constant,
        one_body=active_hamiltonian.one_body_tensor,
        two_body=active_hamiltonian.two_body_tensor,
        reference_determinant=reference.astype(np.int8),
        active_mo_occ=np.asarray(active_mo_occ, dtype=float),
        multiplicity=np.asarray(MULTIPLICITY, dtype=np.int64),
        charge=np.asarray(CHARGE, dtype=np.int64),
        active_alpha=np.asarray(alpha, dtype=np.int64),
        active_beta=np.asarray(beta, dtype=np.int64),
    )


def load_triplet_reference(tensor_path: Path, n_qubits: int) -> dict[str, Any]:
    with np.load(tensor_path, allow_pickle=False) as data:
        required = {
            "reference_determinant",
            "multiplicity",
            "charge",
            "active_alpha",
            "active_beta",
        }
        missing = required.difference(data.files)
        if missing:
            raise ValueError(
                f"{tensor_path} is missing O2 triplet metadata {sorted(missing)}. "
                "Regenerate it with --overwrite-tensors."
            )
        reference = np.asarray(data["reference_determinant"], dtype=np.int8).astype(bool)
        multiplicity = int(np.asarray(data["multiplicity"]).item())
        charge = int(np.asarray(data["charge"]).item())
        active_alpha = int(np.asarray(data["active_alpha"]).item())
        active_beta = int(np.asarray(data["active_beta"]).item())
    if reference.size != n_qubits:
        raise ValueError(
            f"Reference determinant length {reference.size} != n_qubits {n_qubits}."
        )
    if multiplicity != MULTIPLICITY or charge != CHARGE:
        raise ValueError(
            f"Unexpected O2 tensor metadata multiplicity={multiplicity}, charge={charge}."
        )
    if active_alpha != int(np.sum(reference[::2])):
        raise ValueError("Stored active_alpha does not match reference determinant.")
    if active_beta != int(np.sum(reference[1::2])):
        raise ValueError("Stored active_beta does not match reference determinant.")
    return {
        "reference": reference,
        "multiplicity": multiplicity,
        "charge": charge,
        "active_alpha": active_alpha,
        "active_beta": active_beta,
    }


def bitstring(reference: np.ndarray) -> str:
    return "".join("1" if value else "0" for value in reference.tolist())


def generate_o2_cases(
    args: argparse.Namespace,
    bases: Sequence[str],
    qubits: Sequence[int],
    bond_length: float,
) -> list[Path]:
    args.library.mkdir(parents=True, exist_ok=True)
    generation_rows: list[dict[str, Any]] = []
    tensor_paths: list[Path] = []

    for basis in bases:
        total_electrons, total_spin_orbitals = molecular_counts(basis)
        sequence = active_space_sequence(
            total_electrons,
            total_spin_orbitals,
            max(qubits),
        )
        by_qubits = {
            occupied + vacant: (occupied, vacant)
            for occupied, vacant in sequence
        }
        unavailable = [q for q in qubits if q not in by_qubits]
        if unavailable:
            raise ValueError(
                f"Requested O2/{basis} active spaces unavailable: {unavailable}. "
                f"Available: {sorted(by_qubits)}"
            )

        case_info: list[tuple[int, int, int, dict[str, Path]]] = []
        for q in qubits:
            active_occupied, active_vacant = by_qubits[q]
            paths = case_paths(
                args.library,
                basis,
                bond_length,
                active_occupied,
                active_vacant,
            )
            config_status = "existing"
            if args.overwrite_configs or not paths["config"].exists():
                write_configuration(
                    paths,
                    basis,
                    bond_length,
                    active_occupied,
                    active_vacant,
                )
                config_status = "written"
            case_info.append((q, active_occupied, active_vacant, paths))
            generation_rows.append(
                {
                    "molecule": MOLECULE,
                    "bond_length": f"{bond_length:.8f}",
                    "basis": basis,
                    "multiplicity": MULTIPLICITY,
                    "charge": CHARGE,
                    "total_electrons": total_electrons,
                    "total_spin_orbitals": total_spin_orbitals,
                    "active_occupied": active_occupied,
                    "active_vacant": active_vacant,
                    "n_qubits": q,
                    "active_alpha_electrons": "",
                    "active_beta_electrons": "",
                    "reference_determinant_bitstring": "",
                    "config_path": str(paths["config"]),
                    "ham1_pickle_path": str(paths["ham1_pickle"]),
                    "ham2_pickle_path": str(paths["ham2_pickle"]),
                    "tensor_path": str(paths["tensor"]),
                    "config_status": config_status,
                    "tensor_status": "",
                    "error": "",
                }
            )

        need_generation = args.overwrite_tensors or any(
            not paths["tensor"].exists()
            for _, _, _, paths in case_info
        )
        if not need_generation:
            # Also verify that existing tensors came from this triplet-aware driver.
            try:
                for q, _, _, paths in case_info:
                    load_triplet_reference(paths["tensor"], q)
            except ValueError:
                need_generation = True

        full_rohf: PyscfMolecularData | None = None
        full_mo_occ: np.ndarray | None = None
        if need_generation:
            full_rohf, full_mo_occ = run_o2_rohf(basis, bond_length)
            if int(full_rohf.n_electrons) != total_electrons:
                raise ValueError(
                    f"ROHF electron count {full_rohf.n_electrons} != expected {total_electrons}."
                )
            if int(full_rohf.n_qubits) != total_spin_orbitals:
                raise ValueError(
                    f"ROHF spin-orbital count {full_rohf.n_qubits} != expected "
                    f"{total_spin_orbitals}."
                )
            # The full ROHF object is shared by all active-space configs for this basis.
            with case_info[0][3]["ham1_pickle"].open("wb") as stream:
                pickle.dump(full_rohf, stream)

        for row, (q, active_occupied, active_vacant, paths) in zip(
            generation_rows[-len(case_info):],
            case_info,
        ):
            try:
                regenerate = args.overwrite_tensors or not paths["tensor"].exists()
                if not regenerate:
                    try:
                        ref_data = load_triplet_reference(paths["tensor"], q)
                    except ValueError:
                        regenerate = True

                if regenerate:
                    if full_rohf is None or full_mo_occ is None:
                        full_rohf, full_mo_occ = run_o2_rohf(basis, bond_length)
                        with paths["ham1_pickle"].open("wb") as stream:
                            pickle.dump(full_rohf, stream)
                    state = load_configuration_file(paths["config"])
                    active_hamiltonian = apply_active_space(state, full_rohf)
                    if int(active_hamiltonian.n_qubits) != q:
                        raise ValueError(
                            f"Generated active Hamiltonian has {active_hamiltonian.n_qubits} "
                            f"qubits, expected {q}."
                        )
                    active_mo_occ, reference = active_reference_from_rohf(
                        full_mo_occ,
                        total_electrons,
                        active_occupied,
                        active_vacant,
                    )
                    with paths["ham2_pickle"].open("wb") as stream:
                        pickle.dump(active_hamiltonian, stream)
                    save_triplet_tensors(
                        paths["tensor"],
                        active_hamiltonian,
                        reference,
                        active_mo_occ,
                    )
                    row["tensor_status"] = "written"
                    ref_data = {
                        "reference": reference,
                        "active_alpha": int(np.sum(reference[::2])),
                        "active_beta": int(np.sum(reference[1::2])),
                    }
                else:
                    row["tensor_status"] = "existing"

                reference = np.asarray(ref_data["reference"], dtype=bool)
                row["active_alpha_electrons"] = int(ref_data["active_alpha"])
                row["active_beta_electrons"] = int(ref_data["active_beta"])
                row["reference_determinant_bitstring"] = bitstring(reference)
                tensor_paths.append(paths["tensor"])
                print(
                    f"  {basis:7s} {q:2d}q -> occ={active_occupied:2d}, "
                    f"vac={active_vacant:2d}, ref={bitstring(reference)}, "
                    f"tensor={row['tensor_status']}",
                    flush=True,
                )
            except Exception as error:
                row["tensor_status"] = "failed"
                row["error"] = f"{type(error).__name__}: {error}"
                print(f"  FAILED O2/{basis}/{q}q: {row['error']}", flush=True)

    generation_summary = args.output_dir / "o2_tensor_generation_summary.csv"
    generation_summary.parent.mkdir(parents=True, exist_ok=True)
    with generation_summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=GENERATION_FIELDS)
        writer.writeheader()
        writer.writerows(generation_rows)

    failures = [row for row in generation_rows if row["error"]]
    if failures:
        raise RuntimeError(
            f"O2 generation had {len(failures)} failed case(s). See {generation_summary}."
        )
    return sorted(
        tensor_paths,
        key=lambda path: (path.parent.name, parse_case_metadata(path, load_interaction_operator(path)[1]).n_qubits),
    )


# -----------------------------------------------------------------------------
# Generalized open-shell reference utilities
# -----------------------------------------------------------------------------


def number_sector_basis_indices_from_reference(reference: np.ndarray) -> np.ndarray:
    reference = np.asarray(reference, dtype=bool)
    occupied = np.where(reference)[0]
    unoccupied = np.where(~reference)[0]
    indices: list[int] = []
    maximum_order = min(len(occupied), len(unoccupied))
    for order in range(maximum_order + 1):
        for occupied_removed, unoccupied_added in itertools.product(
            itertools.combinations(occupied, order),
            itertools.combinations(unoccupied, order),
        ):
            determinant = reference.copy()
            determinant[list(occupied_removed)] = False
            determinant[list(unoccupied_added)] = True
            indices.append(matrixfree.determinant_to_index(determinant))
    return np.asarray(indices, dtype=np.int64)


def spin_sector_basis_indices_from_reference(reference: np.ndarray) -> np.ndarray:
    reference = np.asarray(reference, dtype=bool)
    n_qubits = int(reference.size)
    if n_qubits % 2:
        raise ValueError("Spin-preserving reference requires an even qubit count.")

    occupied_alpha = np.where(reference[::2])[0] * 2
    unoccupied_alpha = np.where(~reference[::2])[0] * 2
    occupied_beta = np.where(reference[1::2])[0] * 2 + 1
    unoccupied_beta = np.where(~reference[1::2])[0] * 2 + 1
    alpha_maximum = min(len(occupied_alpha), len(unoccupied_alpha))
    beta_maximum = min(len(occupied_beta), len(unoccupied_beta))

    indices: list[int] = []
    for total_order in range(alpha_maximum + beta_maximum + 1):
        for alpha_order in range(alpha_maximum + 1):
            beta_order = total_order - alpha_order
            if beta_order < 0 or beta_order > beta_maximum:
                continue
            products = itertools.product(
                itertools.combinations(occupied_alpha, alpha_order),
                itertools.combinations(unoccupied_alpha, alpha_order),
                itertools.combinations(occupied_beta, beta_order),
                itertools.combinations(unoccupied_beta, beta_order),
            )
            for alpha_removed, alpha_added, beta_removed, beta_added in products:
                determinant = reference.copy()
                determinant[list(alpha_removed)] = False
                determinant[list(alpha_added)] = True
                determinant[list(beta_removed)] = False
                determinant[list(beta_added)] = True
                indices.append(matrixfree.determinant_to_index(determinant))
    return np.asarray(indices, dtype=np.int64)


def build_reference_state(reference: np.ndarray) -> np.ndarray:
    reference = np.asarray(reference, dtype=bool)
    state = np.zeros(2 ** int(reference.size), dtype=np.complex128)
    state[matrixfree.determinant_to_index(reference)] = 1.0
    return state


def exact_reference_state_from_reference(
    fermion_hamiltonian: FermionOperator,
    reference: np.ndarray,
    evolution_time: float,
    tolerance: float,
    spin_preserving: bool,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    reference = np.asarray(reference, dtype=bool)
    n_qubits = int(reference.size)
    n_electrons = int(np.sum(reference))
    fermion_without_identity = matrixfree.remove_fermionic_identity(
        fermion_hamiltonian,
        tolerance,
    )
    build_start = time.perf_counter()
    sparse_hamiltonian = get_number_preserving_sparse_operator(
        fermion_without_identity,
        num_qubits=n_qubits,
        num_electrons=n_electrons,
        spin_preserving=spin_preserving,
        reference_determinant=reference,
        excitation_level=None,
    )
    exact_build_time = time.perf_counter() - build_start
    if spin_preserving:
        basis_indices = spin_sector_basis_indices_from_reference(reference)
    else:
        basis_indices = number_sector_basis_indices_from_reference(reference)
    if sparse_hamiltonian.shape[0] != basis_indices.size:
        raise ValueError(
            "OpenFermion sector dimension does not match reconstructed open-shell "
            f"basis order: {sparse_hamiltonian.shape[0]} != {basis_indices.size}."
        )
    initial_sector_state = np.zeros(sparse_hamiltonian.shape[0], dtype=np.complex128)
    initial_sector_state[0] = 1.0
    evolution_start = time.perf_counter()
    exact_sector_state = expm_multiply(
        (-1j * evolution_time) * sparse_hamiltonian,
        initial_sector_state,
    )
    exact_evolution_time = time.perf_counter() - evolution_start
    return (
        np.asarray(exact_sector_state, dtype=np.complex128),
        basis_indices,
        exact_build_time,
        exact_evolution_time,
    )


def build_open_shell_bch_evaluator(
    pauli_keys: Sequence[Any],
    coefficients: dict[Any, complex],
    pauli_graph: Any,
    reference: np.ndarray,
    tolerance: float,
) -> Any:
    """Generalize QHAT's HFCommutatorEvaluator to an arbitrary determinant."""
    n_qubits = int(reference.size)
    reference_index = matrixfree.determinant_to_index(np.asarray(reference, dtype=bool))
    real_coefficients = np.asarray(
        [real_coefficient(coefficients[key], tolerance) for key in pauli_keys],
        dtype=np.float64,
    )
    left_values: list[int] = []
    right_values: list[int] = []
    targets: list[int] = []
    amplitudes: list[complex] = []
    for first, second in pauli_graph.edges:
        left = min(int(first), int(second))
        right = max(int(first), int(second))
        product_phase, product_key = robustness.multiply_pauli_keys(
            pauli_keys[left], pauli_keys[right]
        )
        target, pauli_phase = robustness.apply_pauli_to_basis_index(
            product_key,
            reference_index,
            n_qubits,
        )
        amplitude = (
            2.0
            * real_coefficients[left]
            * real_coefficients[right]
            * product_phase
            * pauli_phase
        )
        if abs(amplitude) <= tolerance:
            continue
        left_values.append(left)
        right_values.append(right)
        targets.append(target)
        amplitudes.append(amplitude)

    if not amplitudes:
        return robustness.HFCommutatorEvaluator(
            number_of_pauli_terms=len(pauli_keys),
            left_indices=np.empty(0, dtype=np.int32),
            right_indices=np.empty(0, dtype=np.int32),
            target_bins=np.empty(0, dtype=np.int32),
            pair_amplitudes=np.empty(0, dtype=np.complex128),
            number_of_target_bins=0,
        )
    _, target_bins = np.unique(np.asarray(targets, dtype=np.int64), return_inverse=True)
    return robustness.HFCommutatorEvaluator(
        number_of_pauli_terms=len(pauli_keys),
        left_indices=np.asarray(left_values, dtype=np.int32),
        right_indices=np.asarray(right_values, dtype=np.int32),
        target_bins=np.asarray(target_bins, dtype=np.int32),
        pair_amplitudes=np.asarray(amplitudes, dtype=np.complex128),
        number_of_target_bins=int(np.max(target_bins)) + 1,
    )


# -----------------------------------------------------------------------------
# Deterministic ordering benchmark with triplet reference
# -----------------------------------------------------------------------------


def result_key(row: dict[str, str]) -> tuple[str, str, int, float] | None:
    try:
        return (
            row["case_id"],
            row["ordering"],
            int(row["trotter_steps"]),
            float(row["evolution_time"]),
        )
    except (KeyError, TypeError, ValueError):
        return None


def load_completed_results(output: Path) -> set[tuple[str, str, int, float]]:
    completed: set[tuple[str, str, int, float]] = set()
    if not output.exists():
        return completed
    with output.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") != "success":
                continue
            key = result_key(row)
            if key is not None:
                completed.add(key)
    return completed


def load_raw_references(output: Path) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    if not output.exists():
        return result
    with output.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") != "success" or row.get("ordering") != "jw_raw":
                continue
            try:
                result[row["case_id"]] = {
                    "one_minus_overlap": float(row["one_minus_overlap"]),
                    "state_infidelity": float(row["state_infidelity"]),
                    "bch2_hf_state_norm": float(row["bch2_hf_state_norm"]),
                }
            except (KeyError, TypeError, ValueError):
                continue
    return result


def benchmark_triplet_case(
    tensor_path: Path,
    args: argparse.Namespace,
    ordering_names: Sequence[str],
    raw_reference: dict[str, float] | None,
) -> list[dict[str, Any]]:
    interaction, n_qubits = load_interaction_operator(tensor_path)
    metadata = parse_case_metadata(tensor_path, n_qubits)
    ref_data = load_triplet_reference(tensor_path, n_qubits)
    reference = np.asarray(ref_data["reference"], dtype=bool)
    n_electrons = int(np.sum(reference))
    if n_electrons != metadata.active_occupied:
        raise ValueError(
            f"Stored O2 reference has {n_electrons} electrons, but filename says "
            f"active_occupied={metadata.active_occupied}."
        )

    fermion_hamiltonian = clean_fermion_operator(
        get_fermion_operator(interaction),
        args.tolerance,
    )
    full_jw_hamiltonian = jordan_wigner(fermion_hamiltonian)
    full_jw_hamiltonian.compress(abs_tol=args.tolerance)
    final_coefficients = {
        key: coefficient
        for key, coefficient in full_jw_hamiltonian.terms.items()
        if key != () and abs(coefficient) > args.tolerance
    }
    raw_pauli_keys = list(final_coefficients)
    if not raw_pauli_keys:
        raise ValueError("Identity-free JW Hamiltonian has no Pauli terms.")

    print("    Building the five existing deterministic ordering definitions...", flush=True)
    orderings = baseline.build_deterministic_orderings(
        fermion_hamiltonian=fermion_hamiltonian,
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        tolerance=args.tolerance,
    )
    raw_index_by_key = {key: index for index, key in enumerate(raw_pauli_keys)}
    ordering_indices = {
        name: [raw_index_by_key[key] for key in keys]
        for name, keys in orderings.items()
    }

    print("    Building open-shell BCH evaluator...", flush=True)
    jw_graph, _ = build_pauli_noncommutation_graph(raw_pauli_keys)
    bch_evaluator = build_open_shell_bch_evaluator(
        pauli_keys=raw_pauli_keys,
        coefficients=final_coefficients,
        pauli_graph=jw_graph,
        reference=reference,
        tolerance=args.tolerance,
    )

    print(
        f"    Building exact sector reference: N={n_electrons}, "
        f"Nalpha={ref_data['active_alpha']}, Nbeta={ref_data['active_beta']}...",
        flush=True,
    )
    (
        exact_sector_state,
        exact_basis_indices,
        exact_build_time,
        exact_evolution_time,
    ) = exact_reference_state_from_reference(
        fermion_hamiltonian=fermion_hamiltonian,
        reference=reference,
        evolution_time=args.evolution_time,
        tolerance=args.tolerance,
        spin_preserving=not args.no_spin_sector,
    )
    number_basis_indices = number_sector_basis_indices_from_reference(reference)
    initial_state = build_reference_state(reference)
    compiled_raw_terms = matrixfree.compile_ordered_terms(
        raw_pauli_keys,
        final_coefficients,
        n_qubits,
        args.tolerance,
    )

    raw_reference = dict(raw_reference or {})
    raw_one_minus_overlap = raw_reference.get("one_minus_overlap")
    raw_infidelity = raw_reference.get("state_infidelity")
    raw_bch_norm = raw_reference.get("bch2_hf_state_norm")
    rows: list[dict[str, Any]] = []

    for ordering_name in ordering_names:
        pauli_keys = orderings[ordering_name]
        pauli_order_indices = ordering_indices[ordering_name]
        print(f"    Running {ordering_name} ...", flush=True)
        row = {field: "" for field in RESULT_FIELDS}
        baseline.add_metadata(row, metadata, tensor_path)
        row.update(
            {
                "status": "success",
                "number_of_pauli_terms": len(raw_pauli_keys),
                "ordering": ordering_name,
                "ordering_definition": baseline.ordering_definition(ordering_name),
                "pauli_order_hash": baseline.hash_pauli_order(pauli_keys, n_qubits),
                "first_pauli_string": baseline.dense_pauli_string(pauli_keys[0], n_qubits),
                "first_coefficient": real_coefficient(
                    final_coefficients[pauli_keys[0]], args.tolerance
                ),
                "last_pauli_string": baseline.dense_pauli_string(pauli_keys[-1], n_qubits),
                "last_coefficient": real_coefficient(
                    final_coefficients[pauli_keys[-1]], args.tolerance
                ),
                "trotter_steps": args.steps,
                "trotter_dt": args.evolution_time / args.steps,
                "evolution_time": args.evolution_time,
                "exact_sector_dimension": exact_sector_state.size,
                "exact_build_time_seconds": exact_build_time,
                "exact_evolution_time_seconds": exact_evolution_time,
                "coefficient_tolerance": args.tolerance,
                "multiplicity": MULTIPLICITY,
                "charge": CHARGE,
                "active_alpha_electrons": ref_data["active_alpha"],
                "active_beta_electrons": ref_data["active_beta"],
                "reference_determinant_bitstring": bitstring(reference),
            }
        )

        bch_norm = bch_evaluator.evaluate(pauli_order_indices)
        row["bch2_hf_state_norm"] = bch_norm
        ordered_terms = [compiled_raw_terms[index] for index in pauli_order_indices]
        trotter_start = time.perf_counter()
        approximate_state, nominal_exponential_count = matrixfree.evolve_trotter_state(
            initial_state=initial_state,
            terms=ordered_terms,
            formula_order=1,
            trotter_steps=args.steps,
            evolution_time=args.evolution_time,
            parallel_threshold=args.parallel_threshold,
        )
        row["trotter_runtime_seconds"] = time.perf_counter() - trotter_start
        row["nominal_exponential_count"] = nominal_exponential_count
        metrics = matrixfree.compare_states(
            exact_sector_state=exact_sector_state,
            exact_basis_indices=exact_basis_indices,
            approximate_full_state=approximate_state,
            number_basis_indices=number_basis_indices,
        )
        row.update(metrics)
        overlap = float(metrics["state_overlap_abs"])
        one_minus_overlap = 1.0 - overlap
        infidelity = float(metrics["state_infidelity"])
        row["one_minus_overlap"] = one_minus_overlap

        if ordering_name == "jw_raw":
            raw_one_minus_overlap = one_minus_overlap
            raw_infidelity = infidelity
            raw_bch_norm = bch_norm
            row["one_minus_overlap_ratio_to_raw_jw"] = 1.0
            row["state_infidelity_ratio_to_raw_jw"] = 1.0
            row["bch_squared_ratio_to_raw_jw"] = 1.0
        else:
            if raw_one_minus_overlap is not None:
                row["one_minus_overlap_ratio_to_raw_jw"] = baseline.safe_ratio(
                    one_minus_overlap,
                    float(raw_one_minus_overlap),
                )
            if raw_infidelity is not None:
                row["state_infidelity_ratio_to_raw_jw"] = baseline.safe_ratio(
                    infidelity,
                    float(raw_infidelity),
                )
            if raw_bch_norm is not None:
                bch_ratio = baseline.safe_ratio(bch_norm, float(raw_bch_norm))
                row["bch_squared_ratio_to_raw_jw"] = (
                    bch_ratio**2 if isinstance(bch_ratio, float) else ""
                )
        rows.append(row)
        print(
            f"      1-overlap={one_minus_overlap:.6e}  "
            f"infidelity={infidelity:.6e}  "
            f"ratio/raw={row['state_infidelity_ratio_to_raw_jw']}  "
            f"BCH_ref={bch_norm:.6e}",
            flush=True,
        )
    return rows


def run_benchmarks(
    args: argparse.Namespace,
    tensor_paths: Sequence[Path],
) -> Path:
    output = args.output_dir / "o2_deterministic_ordering_results.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.fresh_results and output.exists():
        output.unlink()

    completed = load_completed_results(output)
    raw_references = load_raw_references(output)
    append = output.exists() and output.stat().st_size > 0
    mode = "a" if append else "w"

    with output.open(mode, newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_FIELDS)
        if not append:
            writer.writeheader()
            stream.flush()

        for case_index, tensor_path in enumerate(tensor_paths, start=1):
            interaction, n_qubits = load_interaction_operator(tensor_path)
            del interaction
            metadata = parse_case_metadata(tensor_path, n_qubits)
            pending = [
                ordering
                for ordering in baseline.ORDERING_NAMES
                if (metadata.case_id, ordering, args.steps, args.evolution_time)
                not in completed
            ]
            print()
            print("-" * 100)
            print(f"[{case_index}/{len(tensor_paths)}] {metadata.case_id}")
            if not pending:
                print("    SKIP: all five orderings already successful.", flush=True)
                continue

            # Ratios for a resumed case need the previously completed raw row.
            raw_reference = raw_references.get(metadata.case_id)
            if "jw_raw" not in pending and raw_reference is None:
                # Incomplete/corrupt old CSV: recompute raw so later ratios remain valid.
                pending = ["jw_raw"] + pending

            try:
                rows = benchmark_triplet_case(
                    tensor_path=tensor_path,
                    args=args,
                    ordering_names=pending,
                    raw_reference=raw_reference,
                )
                for row in rows:
                    writer.writerow(row)
                    stream.flush()  # protect every completed ordering against shutdown
                    key = result_key({key: str(value) for key, value in row.items()})
                    if key is not None:
                        completed.add(key)
                    if row["ordering"] == "jw_raw":
                        raw_references[metadata.case_id] = {
                            "one_minus_overlap": float(row["one_minus_overlap"]),
                            "state_infidelity": float(row["state_infidelity"]),
                            "bch2_hf_state_norm": float(row["bch2_hf_state_norm"]),
                        }
            except Exception as error:
                print(f"    FAILED: {type(error).__name__}: {error}", flush=True)
                failure = {field: "" for field in RESULT_FIELDS}
                failure.update(
                    {
                        "status": "failed",
                        "error_message": f"{type(error).__name__}: {error}",
                        "case_id": metadata.case_id,
                        "tensor_path": str(tensor_path),
                        "molecule": metadata.molecule,
                        "bond_length": metadata.bond_length,
                        "basis": metadata.basis,
                        "active_occupied": metadata.active_occupied,
                        "active_vacant": metadata.active_vacant,
                        "n_qubits": metadata.n_qubits,
                        "trotter_steps": args.steps,
                        "trotter_dt": args.evolution_time / args.steps,
                        "evolution_time": args.evolution_time,
                        "coefficient_tolerance": args.tolerance,
                        "multiplicity": MULTIPLICITY,
                        "charge": CHARGE,
                    }
                )
                writer.writerow(failure)
                stream.flush()
    return output


# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------


def write_summary(result_path: Path, output_path: Path) -> None:
    successful: dict[tuple[str, str], dict[str, str]] = {}
    with result_path.open("r", newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row.get("status") != "success":
                continue
            # Latest successful duplicate, if any, wins.
            successful[(row["case_id"], row["ordering"])] = row

    by_case: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for (case_id, ordering), row in successful.items():
        by_case[case_id][ordering] = row

    rows: list[dict[str, Any]] = []
    for case_id, methods in by_case.items():
        if "jw_raw" not in methods:
            continue
        raw = methods["jw_raw"]
        raw_error = float(raw["state_infidelity"])
        summary = {field: "" for field in SUMMARY_FIELDS}
        summary.update(
            {
                "case_id": case_id,
                "molecule": raw["molecule"],
                "bond_length": raw["bond_length"],
                "basis": raw["basis"],
                "n_qubits": int(raw["n_qubits"]),
                "active_occupied": int(raw["active_occupied"]),
                "active_vacant": int(raw["active_vacant"]),
                "active_alpha_electrons": int(raw["active_alpha_electrons"]),
                "active_beta_electrons": int(raw["active_beta_electrons"]),
                "jw_raw_state_infidelity": raw_error,
            }
        )
        best_name = "jw_raw"
        best_ratio = 1.0
        for ordering in baseline.ORDERING_NAMES:
            if ordering == "jw_raw" or ordering not in methods:
                continue
            ratio_text = methods[ordering].get("state_infidelity_ratio_to_raw_jw", "")
            try:
                ratio = float(ratio_text)
            except (TypeError, ValueError):
                continue
            summary[f"{ordering}_ratio_to_raw"] = ratio
            if ratio < best_ratio:
                best_ratio = ratio
                best_name = ordering
        summary["best_ordering"] = best_name
        summary["best_state_infidelity_ratio_to_raw"] = best_ratio
        rows.append(summary)

    rows.sort(key=lambda row: (str(row["basis"]), int(row["n_qubits"])))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    if rows:
        print()
        print("Compact O2 summary (state-infidelity ratio to raw JW)")
        print("basis     q   JW-signed   JW-|c|   ferm-signed   ferm-|c|   best")
        for row in rows:
            def fmt(name: str) -> str:
                value = row.get(name, "")
                return f"{float(value):9.4f}" if value != "" else "        -"
            print(
                f"{str(row['basis']):8s} {int(row['n_qubits']):2d} "
                f"{fmt('signed_coefficient_lexicographic_ratio_to_raw')} "
                f"{fmt('jw_magnitude_descending_lexicographic_ratio_to_raw')} "
                f"{fmt('fermionic_signed_coefficient_lexicographic_ratio_to_raw')} "
                f"{fmt('fermionic_magnitude_descending_lexicographic_ratio_to_raw')} "
                f"{row['best_ordering']}"
            )


def main() -> None:
    args = parse_arguments()
    bases, qubits = resolved_requested_grid(args)
    bond_length = (
        float(args.bond_length)
        if args.bond_length is not None
        else pyykko_o2_bond_length()
    )
    if bond_length <= 0.0:
        raise ValueError("--bond-length must be positive.")

    configure_logging(level="verbose", logfile=None)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("O2 triplet deterministic fermionic-aware ordering study")
    print("=" * 100)
    print(f"Mode:              {args.mode}")
    print(f"Molecule:          O2 / {MOLECULE}")
    print(f"Multiplicity:      {MULTIPLICITY}")
    print(f"Charge:            {CHARGE}")
    print(f"Bond length:       {bond_length:.4f} Angstrom")
    print(f"Bases:             {bases}")
    print(f"Active qubits:     {qubits}")
    print(f"Trotter:           first order, t={args.evolution_time}, steps={args.steps}")
    print(f"Library:           {args.library}")
    print(f"Output directory:  {args.output_dir}")
    print("Resume:            ON" if not args.fresh_results else "Resume:            OFF")
    print()
    print("Phase 1/2: generating/verifying O2 triplet configs and tensors")
    tensor_paths = generate_o2_cases(args, bases, qubits, bond_length)

    print()
    print("Phase 2/2: deterministic ordering benchmark")
    matrixfree.warm_up_numba()
    result_path = run_benchmarks(args, tensor_paths)
    summary_path = args.output_dir / "o2_deterministic_ordering_summary.csv"
    write_summary(result_path, summary_path)

    print()
    print("=" * 100)
    print("Complete")
    print(f"Detailed results:   {result_path}")
    print(f"Summary:            {summary_path}")
    print(f"Generation summary: {args.output_dir / 'o2_tensor_generation_summary.csv'}")
    print("=" * 100)


if __name__ == "__main__":
    main()