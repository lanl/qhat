# build_config_L_sweep.py
# -----------------------
# A companion script to build_config.py that replaces the two fixed
# interatomic separations (1 Å and sum-of-radii) with a user-defined
# sweep over a range of bond lengths L.
#
# Changes from build_config.py:
#   - spacing_loop() replaced with L_sweep_loop() using numpy.linspace
#   - L_eq estimated from covalent_radius_pyykko (from mendeleev)
#     (build_config.py uses atomic_radius for its proportional spacing)
#   - L_min = L_min_frac * L_eq, L_max = L_max_frac * L_eq
#   - file_stub uses absolute path so outputs save alongside config files
#   - active_space_loop() adds hard break at max_active qubits
#   - System filter: defaults to the 9 library molecules plus atoms H--Ne;
#     pass a positional argument to override, e.g. "Li-Li,Be-Be" or "C,O"
#   - Monatomic systems use one atom at the origin and have no L sweep
#   - --run flag to call hamgen.py on all configs after writing
#
# Everything else (element_loop, basis_loop, configuration_loop,
# active_space_loop, mapping_loop, write_config) mirrors build_config.py.
#
# Usage:
#   # all 9 library molecules plus atoms H--Ne
#   python build_config_L_sweep.py --L-steps 50 --max-active 12 --run
#
#   # specific molecules
#   python build_config_L_sweep.py Li-Li,Be-Be --L-steps 50 --run
#
#   # all elements up to period 6 (same scope as build_config.py)
#   python build_config_L_sweep.py --all-molecules --max-period 6 --L-steps 50 --run

import argparse
import basis_set_exchange
import math
import mendeleev
import numpy
import subprocess
import sys
from pathlib import Path


# =================================================================================================
# Default system list: existing molecules plus atoms H through Ne
# =================================================================================================

LIBRARY_MOLECULES = {
    ("H",), ("He",), ("Li",), ("Be",), ("B",),
    ("C",), ("N",), ("O",), ("F",), ("Ne",),
    ("H",  "H" ), ("He", "H" ), ("He", "He"),
    ("Li", "H" ), ("Li", "Li"),
    ("Be", "H" ), ("Be", "Be"),
    ("B",  "H" ), ("B",  "B" ),
}

MAX_ATOMIC_NUMBER = 10


# =================================================================================================

class Count:
    def __init__(self):
        self.count = {"sto-6g": 0, "hgbs-5": 0}


# =================================================================================================

def write_config(basis, Z1, atom1, Z2, atom2, L, occupied, vacant, mapping,
                 library_root="library"):
    if atom2 is None:
        stub = f"{atom1}_atom_{basis}"
        path = Path(library_root) / "atoms" / atom1 / basis
    else:
        stub = f"{atom1}-{atom2}_{L:4.2f}_{basis}"
        path = Path(library_root) / f"{atom1}-{atom2}" / f"{L:.2f}" / basis
    extended = f"{stub}_as-{occupied:03d}-{vacant:03d}_{mapping}"
    path.mkdir(parents=True, exist_ok=True)
    abs_stub = str(path.resolve() / stub)

    filename = path / f"{extended}.config"
    with open(filename, 'w') as fout:
        print(f'general.print_verbose()',                                              file=fout)
        print(f'general.logfile = "{abs_stub}_as-{occupied:03d}-{vacant:03d}_{mapping}.log"', file=fout)
        print(f'general.file_stub = "{abs_stub}"',                                    file=fout)
        print(f'general.file_format = "default"',                                     file=fout)
        if atom2 is None:
            print(f'hamiltonian.add_atom("{atom1}", 0.0, 0.0, 0.0)',                  file=fout)
        else:
            print(f'L = {L}',                                                          file=fout)
            print(f'hamiltonian.add_atom("{atom1}", -0.5 * L, 0.0, 0.0)',             file=fout)
            print(f'hamiltonian.add_atom("{atom2}",  0.5 * L, 0.0, 0.0)',             file=fout)
        print(f'hamiltonian.basis = "{basis}"',                                        file=fout)
        print(f'hamiltonian.num_active_occupied = {occupied}',                         file=fout)
        print(f'hamiltonian.num_active_vacant = {vacant}',                             file=fout)
        print(f'hamiltonian.f2q_mapping = "{mapping}"',                                file=fout)
    return filename


# =================================================================================================

def norb_for_shell(c):
    return {"s": 2, "p": 6, "d": 10, "f": 14, "g": 18}[c]


def get_orbital_count(basis, element):
    try:
        basis_string = basis_set_exchange.get_basis(basis, elements=element, fmt="nwchem")
    except KeyError:
        return None
    orbital_count = 0
    for line in basis_string.split('\n'):
        if line[:11] == "#BASIS SET:":
            tokens   = line.split('->')
            orbitals = tokens[1][2:-1]
            for orbital in orbitals.split(','):
                idx = 0
                for c in orbital:
                    if c in "0123456789":
                        idx += 1
                count = int(orbital[:idx])
                shell = orbital[idx:]
                if shell not in {"s", "p", "d", "f", "g"}:
                    continue
                orbital_count += count * norb_for_shell(shell)
    return orbital_count


def load_elements():
    elements = []
    for element in mendeleev.get_all_elements():
        radius     = 0.01 * element.atomic_radius \
                     if element.atomic_radius is not None else None
        cov_radius = 0.01 * element.covalent_radius_pyykko \
                     if element.covalent_radius_pyykko is not None else radius
        elements.append({
            "atomic number" : element.atomic_number,
            "name"          : element.name,
            "symbol"        : element.symbol,
            "group"         : element.group.symbol if element.group is not None else '',
            "period"        : element.period,
            "radius"        : radius,
            "cov_radius"    : cov_radius,
            "sto-6g"        : get_orbital_count("sto-6g", element.atomic_number),
            "hgbs-5"        : get_orbital_count("hgbs-5", element.atomic_number),
        })
    return elements


# =================================================================================================

def do_the_thing(elements, element1, basis, configuration, element2, L,
                 total_orbitals, n_act_occ, n_act_vac, mapping, count, indent,
                 library_root, config_files):
    count.count[basis] += 1
    sym1 = element1["symbol"]
    sym2 = element2["symbol"] if element2 is not None else None
    Z1   = element1["atomic number"]
    Z2   = element2["atomic number"] if element2 is not None else None
    system = f"{sym1:>2s}-{sym2:<2s}" if sym2 is not None else f"{sym1:>2s} atom"
    spacing = f"L={L:5.2f}Å" if L is not None else " " * 8
    print('  ' * indent + '  '.join([
        f"{count.count['sto-6g']:07d}",
        f"{count.count['hgbs-5']:07d}",
        f"{configuration[:7]:7s}",
        system,
        spacing,
        f"{basis:6s}",
        f"{mapping:2s}",
        f"{total_orbitals:3}",
        f"{n_act_occ:3}",
        f"{n_act_vac:3}",
    ]))
    fname = write_config(basis, Z1, sym1, Z2, sym2, L, n_act_occ, n_act_vac,
                         mapping, library_root)
    config_files.append(fname)


# =================================================================================================

def mapping_loop(elements, element1, basis, configuration, element2, L,
                 total_orbitals, n_act_occ, n_act_vac, count, indent,
                 library_root, config_files):
    for mapping in ["JW", "BK"]:
        do_the_thing(elements, element1, basis, configuration, element2, L,
                     total_orbitals, n_act_occ, n_act_vac, mapping, count, indent,
                     library_root, config_files)


# =================================================================================================

def active_space_loop(elements, element1, basis, configuration, element2, L,
                      count, indent, library_root, config_files, max_active):
    ratio_ideal     = 0.4
    total_electrons = element1["atomic number"]
    total_orbitals  = element1[basis]
    if element2 is not None:
        total_electrons += element2["atomic number"]
        total_orbitals  += element2[basis]
    total_vacancies = total_orbitals - total_electrons
    n_act_occ       = 2 - total_electrons % 2
    n_act_vac       = n_act_occ
    if total_vacancies < n_act_vac:
        system = element1["symbol"] if element2 is None else \
                 f"{element1['symbol']}-{element2['symbol']}"
        print('  ' * indent +
              f"SKIP: {system} {basis} has {total_vacancies} vacant spin "
              f"orbitals; the first active space requires {n_act_vac}.")
        return
    active_lo       = n_act_occ + n_act_vac
    active_hi       = total_orbitals + 1

    for active_size in range(active_lo, active_hi, 2):
        assert n_act_occ + n_act_vac == active_size
        assert (total_electrons - n_act_occ) % 2 == 0

        # Hard stop at max_active — unlike build_config.py which keeps
        # powers-of-2 and full-space sizes beyond max_active
        if active_size > max_active:
            break

        clause1 = active_size <= max_active
        clause2 = 2 ** int(math.log2(active_size)) == active_size
        clause3 = active_size == (total_orbitals - total_orbitals % 2)

        if clause1 or clause2 or clause3:
            mapping_loop(elements, element1, basis, configuration, element2, L,
                         total_orbitals, n_act_occ, n_act_vac, count, indent,
                         library_root, config_files)

        ratio_md = (n_act_occ + 1) / (n_act_occ + n_act_vac + 2)
        if ratio_md >= ratio_ideal:
            if n_act_vac + 2 <= total_vacancies:
                n_act_vac += 2
            else:
                n_act_occ += 2
        else:
            if n_act_occ + 2 <= total_electrons:
                n_act_occ += 2
            else:
                n_act_vac += 2


# =================================================================================================

def L_sweep_loop(elements, element1, basis, configuration, element2, count, indent,
                 L_steps, L_min_frac, L_max_frac, library_root, config_files, max_active):
    """
    Replaces spacing_loop() from build_config.py.
    Estimates L_eq from sum of Pyykko covalent radii, then sweeps
    L_min = L_min_frac * L_eq .. L_max = L_max_frac * L_eq.
    """
    r1 = element1["cov_radius"]
    r2 = element2["cov_radius"]
    if r1 is None or r2 is None:
        print('  ' * indent +
              f"SKIP: no covalent radius for {element1['symbol']} or {element2['symbol']}")
        return

    L_eq     = r1 + r2
    L_min    = round(L_min_frac * L_eq, 4)
    L_max    = round(L_max_frac * L_eq, 4)
    L_values = numpy.linspace(L_min, L_max, L_steps).round(4)

    print('  ' * indent +
          f"L_SWEEP: L_eq={L_eq:.2f}  L_min={L_min:.2f}  L_max={L_max:.2f}  steps={L_steps}")

    for L in L_values:
        active_space_loop(elements, element1, basis, configuration, element2,
                          float(L), count, indent + 1, library_root, config_files,
                          max_active)


# =================================================================================================

def configuration_loop(elements, element1, basis, count, indent,
                        L_steps, L_min_frac, L_max_frac, library_root, config_files,
                        max_active, molecules=None):
    for configuration in ["hydride", "homonuclear"]:
        if configuration == "homonuclear":
            element2 = element1
        elif configuration == "hydride":
            if element1["atomic number"] == 1:
                continue
            element2 = elements[0]
            assert element2["name"] == "Hydrogen"
        else:
            raise NotImplementedError

        if molecules is not None:
            sym1 = element1["symbol"]
            sym2 = element2["symbol"]
            if (sym1, sym2) not in molecules and (sym2, sym1) not in molecules:
                continue

        print('  ' * indent + f"CONFIGURATION_LOOP: {configuration}")
        L_sweep_loop(elements, element1, basis, configuration, element2,
                     count, indent + 1, L_steps, L_min_frac, L_max_frac,
                     library_root, config_files, max_active)


# =================================================================================================

def basis_loop(elements, element1, count, indent,
               L_steps, L_min_frac, L_max_frac, library_root, config_files,
               max_active, molecules=None):
    for basis in ["sto-6g", "hgbs-5"]:
        if element1[basis] is None:
            continue
        print('  ' * indent + f"BASIS_LOOP: {basis}")
        configuration_loop(elements, element1, basis, count, indent + 1,
                           L_steps, L_min_frac, L_max_frac, library_root,
                           config_files, max_active, molecules=molecules)


# =================================================================================================

def element_loop(elements, count, indent, max_period,
                 L_steps, L_min_frac, L_max_frac, library_root, config_files,
                 max_active, molecules=None):
    for element1 in elements:
        if element1["period"] > max_period:
            continue
        if molecules is not None:
            sym = element1["symbol"]
            if not any(sym in pair for pair in molecules):
                continue
        print('  ' * indent + f"ELEMENT_LOOP: {element1['name']}")
        basis_loop(elements, element1, count, indent + 1,
                   L_steps, L_min_frac, L_max_frac, library_root, config_files,
                   max_active, molecules=molecules)


# =================================================================================================

def atom_loop(elements, count, indent, max_period, library_root, config_files,
              max_active, molecules=None):
    for element in elements:
        if element["atomic number"] > MAX_ATOMIC_NUMBER:
            continue
        if element["period"] > max_period:
            continue
        if molecules is not None and (element["symbol"],) not in molecules:
            continue
        print('  ' * indent + f"ATOM: {element['name']}")
        for basis in ["sto-6g", "hgbs-5"]:
            if element[basis] is None:
                continue
            active_space_loop(
                elements, element, basis, "atom", None, None,
                count, indent + 1, library_root, config_files, max_active
            )


# =================================================================================================

def run_hamgen(config_files, hamgen_dir):
    hamgen_path = Path(hamgen_dir) / "hamgen.py"
    if not hamgen_path.exists():
        print(f"\nERROR: hamgen.py not found at '{hamgen_path}'.")
        sys.exit(1)
    print(f"\nRunning hamgen.py on {len(config_files)} configs...")
    for i, cfg in enumerate(config_files, 1):
        print(f"[{i}/{len(config_files)}] {cfg.name}")
        result = subprocess.run(
            [sys.executable, str(hamgen_path), str(cfg)],
            capture_output=False,
        )
        if result.returncode != 0:
            print(f"  WARNING: exited with code {result.returncode}")
        else:
            print(f"  Done.")


# =================================================================================================

def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate hamgen configs sweeping bond length L. "
                    "Companion script to build_config.py."
    )
    ap.add_argument("molecules",       type=str, nargs="?", default=None,
                    help="Comma-separated systems, e.g. Li-Li,Be-Be or C,O "
                         "(default: 9 molecules plus atoms H through Ne)")
    ap.add_argument("--all-molecules", action="store_true", dest="all_molecules",
                    help="Generate for all elements up to --max-period (no filter)")
    ap.add_argument("--L-steps",       type=int,   default=10,  dest="L_steps",
                    help="Number of bond length points (default: 10)")
    ap.add_argument("--L-min-frac",    type=float, default=0.6, dest="L_min_frac",
                    help="L_min = frac * L_eq (default: 0.6)")
    ap.add_argument("--L-max-frac",    type=float, default=3.0, dest="L_max_frac",
                    help="L_max = frac * L_eq (default: 3.0)")
    ap.add_argument("--max-active",    type=int,   default=12,  dest="max_active",
                    help="Max active space in qubits (default: 12)")
    ap.add_argument("--max-period",    type=int,   default=6,   dest="max_period",
                    help="Max element period (default: 6)")
    ap.add_argument("--library",       type=str,   default="library",
                    help="Root directory for library (default: library)")
    ap.add_argument("--hamgen-dir",    type=str,   default=".", dest="hamgen_dir",
                    help="Directory containing hamgen.py (default: .)")
    ap.add_argument("--run",           action="store_true",
                    help="Run hamgen.py on all configs after generating them")
    return ap.parse_args()


# =================================================================================================

def main():
    args         = parse_args()
    elements     = load_elements()
    config_files = []

    if args.all_molecules:
        molecules = None
    elif args.molecules:
        molecules = set()
        for m in args.molecules.split(","):
            parts = tuple(p for p in m.strip().split("-") if p)
            if len(parts) not in {1, 2}:
                print(f"ERROR: invalid system specification {m!r}.")
                sys.exit(1)
            molecules.add(parts)
    else:
        molecules = LIBRARY_MOLECULES

    count = Count()
    element_loop(
        elements, count, indent=0,
        max_period   = args.max_period,
        L_steps      = args.L_steps,
        L_min_frac   = args.L_min_frac,
        L_max_frac   = args.L_max_frac,
        library_root = args.library,
        config_files = config_files,
        max_active   = args.max_active,
        molecules    = molecules,
    )
    atom_loop(
        elements, count, indent=0,
        max_period   = args.max_period,
        library_root = args.library,
        config_files = config_files,
        max_active   = args.max_active,
        molecules    = molecules,
    )

    print(f"\nTotal configs written: {len(config_files)}")

    if args.run:
        run_hamgen(config_files, args.hamgen_dir)
    else:
        print("\nTo run hamgen on all configs:")
        print("  python build_config_L_sweep.py --run")


if __name__ == "__main__":
    main()