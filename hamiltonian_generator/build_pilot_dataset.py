"""
Build and optionally run a small Hamiltonian dataset pilot.

The existing build_config.py creates a broad library sweep.  This script is
more targeted: it emits configs that include both small active spaces and
20-30 qubit active spaces for representative atomic numbers up to Z=30.
"""

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys


DEFAULT_ELEMENTS = [1, 3, 4, 6, 8, 12, 14, 20, 26, 30]
DEFAULT_TARGET_QUBITS = [2, 4, 6, 8, 10, 12, 14, 16, 20, 22, 24, 26, 28, 30]
DEFAULT_MAPPINGS = ["JW", "BK"]
DEFAULT_FAMILIES = ["homonuclear", "hydride"]
MAPPING_ORDER = {"JW": 0, "BK": 1}


def parse_csv_ints(value):
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def parse_csv_strings(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def norb_for_shell(shell):
    if shell == "s":
        return 2
    if shell == "p":
        return 6
    if shell == "d":
        return 10
    if shell == "f":
        return 14
    if shell == "g":
        return 18
    raise KeyError(f"Invalid shell '{shell}'")


def get_orbital_count(basis_set_exchange, basis, atomic_number):
    try:
        basis_string = basis_set_exchange.get_basis(basis, elements=atomic_number, fmt="nwchem")
    except KeyError:
        return None

    orbital_count = 0
    for line in basis_string.splitlines():
        if not line.startswith("#BASIS SET:"):
            continue
        orbitals = line.split("->")[1][2:-1]
        for orbital in orbitals.split(","):
            idx = 0
            for char in orbital:
                if char.isdigit():
                    idx += 1
                else:
                    break
            count = int(orbital[:idx])
            shell = orbital[idx:]
            orbital_count += count * norb_for_shell(shell)
    return orbital_count


def load_elements(atomic_numbers, basis):
    try:
        import basis_set_exchange
        import mendeleev
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Missing Hamiltonian generator dependencies. Install the project dependencies first, "
            "for example with: python3 -m pip install -e ."
        ) from exc

    elements = {}
    for atomic_number in atomic_numbers:
        element = mendeleev.element(atomic_number)
        radius = 0.01 * element.atomic_radius if element.atomic_radius is not None else None
        orbital_count = get_orbital_count(basis_set_exchange, basis, atomic_number)
        if orbital_count is None:
            continue
        elements[atomic_number] = {
            "atomic_number": atomic_number,
            "symbol": element.symbol,
            "radius": radius,
            "orbitals": orbital_count,
        }
    return elements


def choose_active_space(total_electrons, total_orbitals, target_qubits):
    total_vacancies = total_orbitals - total_electrons
    if target_qubits > total_orbitals or total_vacancies <= 0:
        return None

    candidates = []
    for active_occupied in range(1, min(total_electrons, target_qubits - 1) + 1):
        active_vacant = target_qubits - active_occupied
        if active_vacant <= 0 or active_vacant > total_vacancies:
            continue
        if (total_electrons - active_occupied) % 2 != 0:
            continue
        if target_qubits % 2 != 0:
            continue
        ratio_error = abs((active_occupied / target_qubits) - 0.4)
        candidates.append((ratio_error, active_occupied, active_vacant))

    if not candidates:
        return None

    _, active_occupied, active_vacant = min(candidates)
    return active_occupied, active_vacant


def molecule_specs(elements, families):
    hydrogen = elements.get(1)
    for element in elements.values():
        if "homonuclear" in families:
            yield "homonuclear", element, element
        if "hydride" in families and element["atomic_number"] != 1 and hydrogen is not None:
            yield "hydride", element, hydrogen


def molecule_spacing(element1, element2):
    if element1["radius"] is not None and element2["radius"] is not None:
        return element1["radius"] + element2["radius"]
    return 1.0


def config_text(row, write_initial_state):
    flag = "True" if write_initial_state else "False"
    return "\n".join([
        "general.print_verbose()",
        f'general.logfile = "{row["name"]}.log"',
        f'general.file_stub = "{row["stub"]}"',
        'general.file_format = "default"',
        f"general.write_initial_state = {flag}",
        f"L = {row['spacing']:.8f}",
        f'hamiltonian.add_atom("{row["atom1"]}", -0.5 * L, 0.0, 0.0)',
        f'hamiltonian.add_atom("{row["atom2"]}",  0.5 * L, 0.0, 0.0)',
        f'hamiltonian.basis = "{row["basis"]}"',
        f'hamiltonian.num_active_occupied = {row["active_occupied"]}',
        f'hamiltonian.num_active_vacant = {row["active_vacant"]}',
        f'hamiltonian.f2q_mapping = "{row["mapping"]}"',
        "",
    ])


def build_rows(args):
    elements = load_elements(args.elements, args.basis)
    rows = []
    seen = set()

    for family, element1, element2 in molecule_specs(elements, args.families):
        total_electrons = element1["atomic_number"] + element2["atomic_number"]
        total_orbitals = element1["orbitals"] + element2["orbitals"]
        spacing = molecule_spacing(element1, element2)
        for target_qubits in args.target_qubits:
            active_space = choose_active_space(total_electrons, total_orbitals, target_qubits)
            if active_space is None:
                continue
            active_occupied, active_vacant = active_space
            for mapping in args.mappings:
                stub = f"{element1['symbol']}-{element2['symbol']}_{spacing:4.2f}_{args.basis}"
                name = f"{stub}_as-{active_occupied:03d}-{active_vacant:03d}_{mapping.lower()}"
                key = (name, family)
                if key in seen:
                    continue
                seen.add(key)
                rows.append({
                    "name": name,
                    "family": family,
                    "basis": args.basis,
                    "mapping": mapping,
                    "z1": element1["atomic_number"],
                    "z2": element2["atomic_number"],
                    "atom1": element1["symbol"],
                    "atom2": element2["symbol"],
                    "spacing": spacing,
                    "target_qubits": target_qubits,
                    "active_occupied": active_occupied,
                    "active_vacant": active_vacant,
                    "total_electrons": total_electrons,
                    "total_orbitals": total_orbitals,
                    "stub": stub,
                })

    rows.sort(key=lambda row: (
        row["target_qubits"],
        row["z1"],
        row["z2"],
        row["family"],
        row["mapping"],
    ))
    return limit_rows(rows, args.max_configs)


def limit_rows(rows, max_configs):
    if not max_configs or len(rows) <= max_configs:
        return rows

    grouped = {}
    for row in rows:
        case_key = (
            row["family"],
            row["basis"],
            row["z1"],
            row["z2"],
            row["target_qubits"],
            row["active_occupied"],
            row["active_vacant"],
            row["spacing"],
        )
        grouped.setdefault(row["target_qubits"], {}).setdefault(case_key, []).append(row)

    for target_group in grouped.values():
        for case_rows in target_group.values():
            case_rows.sort(key=lambda row: MAPPING_ORDER.get(row["mapping"], 99))

    selected = []
    while len(selected) < max_configs and grouped:
        for target_qubits in sorted(grouped):
            case_groups = grouped[target_qubits]
            if case_groups:
                case_key = sorted(case_groups)[0]
                case_rows = case_groups.pop(case_key)
                if len(selected) + len(case_rows) > max_configs:
                    return selected
                selected.extend(case_rows)
                if len(selected) == max_configs:
                    break
        grouped = {key: value for key, value in grouped.items() if value}

    return selected


def write_configs(rows, output_dir, write_initial_state):
    output_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        config_dir = (
            output_dir
            / f"{row['z1']:03d}-{row['z2']:03d}_{row['atom1']}-{row['atom2']}"
            / f"{row['target_qubits']:02d}q"
            / row["basis"]
        )
        config_dir.mkdir(parents=True, exist_ok=True)
        config_path = config_dir / f"{row['name']}.config"
        config_path.write_text(config_text(row, write_initial_state), encoding="utf-8")
        row["config_path"] = str(config_path)
        row["output_path"] = str(config_dir / f"{row['name']}.dat")


def write_manifest(rows, manifest_path):
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "name",
        "family",
        "basis",
        "mapping",
        "z1",
        "z2",
        "atom1",
        "atom2",
        "target_qubits",
        "active_occupied",
        "active_vacant",
        "total_electrons",
        "total_orbitals",
        "spacing",
        "config_path",
        "output_path",
        "status",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def run_configs(rows, repo_root, keep_going):
    hamgen = repo_root / "hamiltonian_generator" / "hamgen.py"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root.parent) + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("MPLCONFIGDIR", "/tmp/qhat_matplotlib")

    for row in rows:
        config_path = Path(row["config_path"]).resolve()
        result = subprocess.run(
            [sys.executable, str(hamgen), str(config_path)],
            cwd=config_path.parent,
            env=env,
            check=False,
        )
        row["status"] = "generated" if result.returncode == 0 else f"failed:{result.returncode}"
        if result.returncode != 0 and not keep_going:
            raise SystemExit(result.returncode)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="library_pilot")
    parser.add_argument("--basis", default="sto-6g")
    parser.add_argument("--elements", type=parse_csv_ints, default=DEFAULT_ELEMENTS)
    parser.add_argument("--target-qubits", type=parse_csv_ints, default=DEFAULT_TARGET_QUBITS)
    parser.add_argument("--mappings", type=parse_csv_strings, default=DEFAULT_MAPPINGS)
    parser.add_argument("--families", type=parse_csv_strings, default=DEFAULT_FAMILIES)
    parser.add_argument("--max-configs", type=int, default=40)
    parser.add_argument("--write-initial-state", action="store_true")
    parser.add_argument("--run", action="store_true", help="Run hamgen.py for each generated config.")
    parser.add_argument("--keep-going", action="store_true", help="Continue after failed hamgen runs.")
    return parser.parse_args()


def main():
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    output_dir = (repo_root / args.output_dir).resolve()
    manifest_path = output_dir / "manifest.csv"

    rows = build_rows(args)
    write_configs(rows, output_dir, args.write_initial_state)
    for row in rows:
        row["status"] = "config-only"

    if args.run:
        run_configs(rows, repo_root, args.keep_going)

    write_manifest(rows, manifest_path)
    print(f"Wrote {len(rows)} configs to {output_dir}")
    print(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    main()
