import basis_set_exchange
import math
import mendeleev
import os
import pprint

# =================================================================================================

class Count:
    count = {"sto-6g": 0, "hgbs-5": 0}

# =================================================================================================

def write_config(basis, Z1, atom1, Z2, atom2, spacing, occupied, vacant, mapping):
    stub = f"{atom1}-{atom2}_{spacing:4.2f}_{basis}"
    extended = f"{stub}_as-{occupied:03d}-{vacant:03d}_{mapping}"
    path = f"library/{Z1:03d}-{Z2:03d}_{atom1}-{atom2}/{spacing:4.2f}/{basis}"

    os.makedirs(path, exist_ok=True)
    filename = path + "/" + f"{extended}.config"
    with open(filename, 'w') as fout:
        print(f'general.print_verbose()', file=fout)
        print(f'general.logfile = "{extended}.log"', file=fout)
        print(f'general.file_stub = "{stub}"', file=fout)
        print(f'general.file_format = "default"', file=fout)
        print(f'L = {spacing}', file=fout)
        print(f'hamiltonian.add_atom("{atom1}", -0.5 * L, 0.0, 0.0)', file=fout)
        print(f'hamiltonian.add_atom("{atom2}",  0.5 * L, 0.0, 0.0)', file=fout)
        print(f'hamiltonian.basis = "{basis}"', file=fout)
        print(f'hamiltonian.num_active_occupied = {occupied}', file=fout)
        print(f'hamiltonian.num_active_vacant = {vacant}', file=fout)
        print(f'hamiltonian.f2q_mapping = "{mapping}"', file=fout)

# =================================================================================================

def norb_for_shell(c):
    if c == 's':
        return 2
    elif c == 'p':
        return 6
    elif c == 'd':
        return 10
    elif c == 'f':
        return 14
    elif c == 'g':
        return 18
    else:
        raise KeyError(f"Invalid shell '{c}'")

# =================================================================================================

def get_orbital_count(basis, element):
    try:
        basis_string = basis_set_exchange.get_basis(basis, elements=element, fmt="nwchem")
    except KeyError:
        return None
    orbital_count = 0
    for line in basis_string.split('\n'):
        if line[:11] == "#BASIS SET:":
            tokens = line.split('->')
            orbitals = tokens[1][2:-1]
            orbital_list = orbitals.split(',')
            for orbital in orbital_list:
                idx = 0
                for c in orbital:
                    if c in "0123456789":
                        idx += 1
                count = int(orbital[:idx])
                shell = orbital[idx:]
                shell_norb = norb_for_shell(shell)
                norb = count * shell_norb
                orbital_count += norb
    return orbital_count


# =================================================================================================
def load_elements():
    we = 0
    for element in mendeleev.get_all_elements():
        we = max(we, len(element.name))
    wr = 6
    elements = list()
    for element in mendeleev.get_all_elements():
        radius = 0.01 * element.atomic_radius if element.atomic_radius is not None else None
        sto6g_count = get_orbital_count("sto-6g", element.atomic_number)
        hgbs5_count = get_orbital_count("hgbs-5", element.atomic_number)
        print('   '.join([
            f"{element.atomic_number:3}",
            f"{element.name:{we}}",
            f"{element.symbol:2}",
            f"{radius:{wr}.{wr-2}f}" if radius is not None else f"{' ':{wr}}",
            f"{sto6g_count:3}" if sto6g_count is not None else f"{' ':3}",
            f"{hgbs5_count:3}" if hgbs5_count is not None else f"{' ':3}",
            f"{element.period:1}",
            f"{element.group.symbol:5}" if element.group is not None else f"{' ':5}",
            ]))
        elements.append({
            "atomic number": element.atomic_number,
            "name": element.name,
            "symbol": element.symbol,
            "group": element.group.symbol if element.group is not None else '',
            "period": element.period,
            "radius": radius,
            "sto-6g": sto6g_count,
            "hgbs-5": hgbs5_count,
            })
    return elements

# =================================================================================================

def do_the_thing(elements, element1, basis, configuration, element2, spacing, isep, total_orbitals,
                 n_act_occ, n_act_vac, mapping, count, indent):
    count.count[basis] += 1
    Z1 = element1["atomic number"]
    Z2 = element2["atomic number"]
    sym1 = element1["symbol"]
    sym2 = element2["symbol"]
    message = '  ' * indent + '  '.join([
        f"{count.count['sto-6g']:07d}",
        f"{count.count['hgbs-5']:07d}",
        f"{configuration[:7]:7s}",
        f"{sym1:>2s}-{sym2:<2s}",
        f"{spacing[:3]:3s}",
        f"{isep:4.2f}Å",
        f"{basis:6s}",
        f"{mapping:2s}",
        f"{total_orbitals:3}",
        f"{n_act_occ:3}",
        f"{n_act_vac:3}",
        f"{n_act_occ/(n_act_occ+n_act_vac):4.2f}",
        ])
    print(message)
    write_config(basis, Z1, sym1, Z2, sym2, isep, n_act_occ, n_act_vac, mapping)

# =================================================================================================

def mapping_loop(elements, element1, basis, configuration, element2, spacing, isep, total_orbitals,
                 n_act_occ, n_act_vac, count, indent):
    # Loop over fermion-to-qubit
    for mapping in ["JW", "BK"]:
        print('  ' * indent + f"MAPPING_LOOP: {mapping}")
        do_the_thing(elements, element1, basis, configuration, element2, spacing, isep,
                     total_orbitals, n_act_occ, n_act_vac, mapping, count, indent+1)

# =================================================================================================

def active_space_loop(elements, element1, basis, configuration, element2, spacing, isep, count,
                      indent):
    # Loop over active space
    # -- constraints:
    #    -- the number of frozen, occupied orbitals must be even (from pySCF)
    #    -- the number of active orbitals must be even (from pySCF)
    # -- We want to balance the occupied vs vacant active orbitals so that about
    #    40% of the active space is occupied.  This number is a little arbitrary,
    #    but sufficient for now.
    ratio_ideal = 0.4
    total_electrons = element1["atomic number"] + element2["atomic number"]
    total_orbitals = element1[basis] + element2[basis]
    total_vacancies = total_orbitals - total_electrons
    n_act_occ = 2 - total_electrons % 2
    n_act_vac = n_act_occ
    active_lo = n_act_occ + n_act_vac
    active_hi = total_orbitals + 1
    for active_size in range(active_lo, active_hi, 2):
        assert n_act_occ + n_act_vac == active_size
        assert (total_electrons - n_act_occ) % 2 == 0
        clause1 = active_size <= 16
        clause2 = 2**int(math.log2(active_size)) == active_size
        clause3 = active_size == (total_orbitals - total_orbitals % 2)

        if clause1 or clause2 or clause3:
            print('  ' * indent + f"ACTIVE_SPACE_LOOP: {active_size}")
            mapping_loop(elements, element1, basis, configuration, element2, spacing, isep,
                         total_orbitals, n_act_occ, n_act_vac, count, indent+1)

        # Increase the number of active orbitals for the next iteration
        ratio_md = (n_act_occ + 1) / (n_act_occ + n_act_vac + 2)
        if ratio_md >= ratio_ideal:
            # increase number of active, vacant orbitals (if possible)
            if n_act_vac + 2 <= total_vacancies:
                n_act_vac += 2
            else:
                n_act_occ += 2
        else:
            # increase number of active, occupied orbitals (if possible)
            if n_act_occ + 2 <= total_electrons:
                n_act_occ += 2
            else:
                n_act_vac += 2

# =================================================================================================

def spacing_loop(elements, element1, basis, configuration, element2, count, indent):
    # Loop over spacing
    for spacing in ["fixed", "proportional"]:
        if spacing == "fixed":
            isep = 1.0
        elif spacing == "proportional":
            if element1["radius"] is None:
                # We don't have an atomic radius, so we can't compute the
                # "proportional" interatomic separation
                continue
            isep = element1["radius"] + element2["radius"]
        else:
            raise NotImplementedError("invalid configuration")
        print('  ' * indent + f"SPACING_LOOP: {spacing}")
        active_space_loop(elements, element1, basis, configuration, element2, spacing, isep,
                          count, indent+1)

# =================================================================================================

def configuration_loop(elements, element1, basis, count, indent):
    # Loop over configuration
    for configuration in ["hydride", "homonuclear"]:
        if configuration == "homonuclear":
            element2 = element1
        elif configuration == "hydride":
            if element1["atomic number"] == 1:
                # hydrogen hydride is the same as homonuclear hydrogen: don't duplicate
                continue
            element2 = elements[0]
            assert element2["name"] == "Hydrogen", \
                    f"element 2 should be Hydrogen but instead is {element2['name']}"
        else:
            raise NotImplementedError("invalid configuration")
        print('  ' * indent + f"CONFIGURATION_LOOP: {configuration}")
        spacing_loop(elements, element1, basis, configuration, element2, count, indent+1)

# =================================================================================================

def basis_loop(elements, element1, count, indent):
    # Loop over each basis set
    for basis in ["sto-6g", "hgbs-5"]:
        if element1[basis] is None:
            # This element doesn't exist in this basis set
            continue
        print('  ' * indent + f"BASIS_LOOP: {basis}")
        configuration_loop(elements, element1, basis, count, indent+1)

# =================================================================================================

def element_loop(elements, count, indent):
    # Loop over each element
    for element1 in elements:
        if element1["period"] > 6:
            # TODO: Have to check rules about data for period 7
            continue
        print('  ' * indent + f"ELEMENT_LOOP: {element1['name']}")
        basis_loop(elements, element1, count, indent+1)

# =================================================================================================

def main():
    elements = load_elements()

    element_loop(elements, Count(), 0)

if __name__ == "__main__":
    main()
