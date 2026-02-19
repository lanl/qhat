from hamgen_types import GeneralConfiguration, GeneralConfigurationUser, \
        HamiltonianConfiguration, State

import argparse
import math
import mendeleev
import numpy as np
from openfermion import InteractionOperator, MolecularData, jordan_wigner, bravyi_kitaev
from openfermionpyscf import PyscfMolecularData
from openfermionpyscf._run_pyscf import compute_integrals, compute_scf, prepare_pyscf_molecule
import pickle
import pprint
import time

# -------------------------------------------------------------------------------------------------

def load_configuration() -> State:

    # Set up and read command-line arguments
    parser = argparse.ArgumentParser()
    default_config = "config.py"
    parser.add_argument(
            "configuration_file",
            nargs='?',
            default=default_config,
            help=f"Name of the configuration file; defaults to \"{default_config}\"")
    args = parser.parse_args()

    # Read the configuration file
    with open(args.configuration_file, 'r') as fin:
        config_script = fin.read()
    if len(config_script) > 0 and config_script[-1] == '\n':
        config_script = config_script[:-1]

    # Execute the configuration file
    general = GeneralConfigurationUser()
    hamiltonian = HamiltonianConfiguration()
    exec(config_script)

    # Build the state (does some post-processing of user configuration)
    state = State(config_script, general, hamiltonian)

    state.log("\n".join([
        f"Contents of configuration file \"{args.configuration_file}\":",
        config_script
        ]))

    return state

# -------------------------------------------------------------------------------------------------

def define_active_orbitals(state, molecule: PyscfMolecularData) -> tuple[int,int]:
    assert(state.config_hamiltonian.num_active_occupied is not None)
    assert(state.config_hamiltonian.num_active_vacant is not None)

    # Total spin (single-occupancy) orbitals (active or frozen)
    n_tot_orb_so = molecule.n_qubits
    n_tot_occ_so = molecule.n_electrons
    n_tot_vac_so = n_tot_orb_so - n_tot_occ_so

    # Active spin (single-occupancy) orbitals
    n_act_occ_so = state.config_hamiltonian.num_active_occupied
    n_act_vac_so = state.config_hamiltonian.num_active_vacant
    n_act_orb_so = n_act_occ_so + n_act_vac_so

    # Frozen spin (single-occupancy) orbitals
    n_frz_orb_so = n_tot_orb_so - n_act_orb_so
    n_frz_occ_so = n_tot_occ_so - n_act_occ_so
    n_frz_vac_so = n_tot_vac_so - n_act_vac_so

    # Informational printout
    msg = [f'single-occupancy (spin) orbital counts:', ]
    msg.append(f'   {"":12}  {"total":6}  {"active":6}  {"frozen":6}')
    msg.append(f'   {"num orbitals":12}  {n_tot_orb_so:6}  {n_act_orb_so:6}  {n_frz_orb_so:6}')
    msg.append(f'   {"num occupied":12}  {n_tot_occ_so:6}  {n_act_occ_so:6}  {n_frz_occ_so:6}')
    msg.append(f'   {"num vacant"  :12}  {n_tot_vac_so:6}  {n_act_vac_so:6}  {n_frz_vac_so:6}')
    state.log('\n'.join(msg))

    # Verify total orbitals
    assert n_tot_occ_so > 0, "Must have at least one occupied orbital"
    assert n_tot_vac_so > 0, "Must have at least one vacant orbital"
    assert n_tot_occ_so + n_tot_vac_so == n_tot_orb_so, "Inconsistent orbital counts"

    # Verify active orbitals
    assert n_act_occ_so > 0, "Must have at least one active occupied orbital"
    assert n_act_vac_so > 0, "Must have at least one active vacant orbital"
    assert n_act_occ_so + n_act_vac_so == n_act_orb_so, "Inconsistent active orbital counts"

    # Verify frozen orbitals
    assert n_frz_occ_so >= 0, "Must have a non-negative number of frozen occupied orbitals"
    assert n_frz_vac_so >= 0, "Must have a non-negative number of frozen vacant orbitals"
    assert n_frz_occ_so + n_frz_vac_so == n_frz_orb_so, "Inconsistent frozen orbital counts"

    # pySCF/OpenFermion works in terms of spatial (double-occupancy) orbitals when listing
    # frozen-occupied and active-all, so the corresponding spin (single-occupancy) counts must be
    # even
    assert n_frz_occ_so % 2 == 0, "Number of frozen occupied orbitals must be even"
    assert n_act_orb_so % 2 == 0, "Number of active orbitals must be even"

    # Indices (single-occupancy/spin)
    idx_frz_occ_so =                  n_frz_occ_so
    idx_act_occ_so = idx_frz_occ_so + n_act_occ_so
    idx_act_vac_so = idx_act_occ_so + n_act_vac_so
    idx_frz_vac_so = idx_act_vac_so + n_frz_vac_so

    w = math.floor(math.log10(n_tot_orb_so)) + 1

    msg = [f'single-occupancy (spin) orbital indices:', ]
    msg.append(f'   {"frozen":6} {"occupied":8}  [{0:{w}}, {idx_frz_occ_so:{w}})')
    msg.append(f'   {"active":6} {"occupied":8}  [{idx_frz_occ_so:{w}}, {idx_act_occ_so:{w}})')
    msg.append(f'   {"active":6} {"vacant":8}  [{idx_act_occ_so:{w}}, {idx_act_vac_so:{w}})')
    msg.append(f'   {"frozen":6} {"vacant":8}  [{idx_act_vac_so:{w}}, {idx_frz_vac_so:{w}})')
    state.log_verbose('\n'.join(msg))

    # Indices (double-occupancy/spatial)
    n_frz_occ_do = n_frz_occ_so // 2
    n_act_orb_do = n_act_orb_so // 2

    idx_act_occ = n_frz_occ_do
    idx_frz_vac = idx_act_occ + n_act_orb_do

    msg = [f'double-occupancy (spatial) orbital indices:', ]
    msg.append(f'   {"frozen":6} {"occupied":8}  [{0:{w}}, {idx_act_occ:{w}})')
    msg.append(f'   {"active":6} {"total":8}  [{idx_act_occ:{w}}, {idx_frz_vac:{w}})')
    state.log_verbose('\n'.join(msg))

    metadata = {
            "n_frz_occ_so": n_frz_occ_so,
            "n_act_occ_so": n_act_occ_so,
            "n_act_vac_so": n_act_vac_so,
            "n_frz_vac_so": n_frz_vac_so,
            }

    return idx_act_occ, idx_frz_vac, metadata

# -------------------------------------------------------------------------------------------------

# Copied from openfermionpyscf.run_pyscf then modified
def my_run_pyscf(state, molecule, verbose=False):
    # Prepare pyscf molecule.
    pyscf_molecule = prepare_pyscf_molecule(molecule)
    molecule.n_orbitals = int(pyscf_molecule.nao_nr())
    molecule.n_qubits = 2 * molecule.n_orbitals
    molecule.nuclear_repulsion = float(pyscf_molecule.energy_nuc())

    # Run SCF.
    pyscf_scf = compute_scf(pyscf_molecule)
    pyscf_scf.verbose = 0
    #pyscf_scf.frozen = 8    # TODO: Get the right number, right now just testing
    # run() passes through various types of machinery to eventually call SCF.scf() (assuming RHF,
    # which is a subclass of SCF)
    pyscf_scf.run()
    molecule.hf_energy = float(pyscf_scf.e_tot)
    if verbose:
        print('Hartree-Fock energy for {} ({} electrons) is {}.'.format(
            molecule.name, molecule.n_electrons, molecule.hf_energy))

    # Hold pyscf data in molecule. They are required to compute density
    # matrices and other quantities.
    molecule._pyscf_data = pyscf_data = {}
    pyscf_data['mol'] = pyscf_molecule
    pyscf_data['scf'] = pyscf_scf

    # Populate fields.
    molecule.canonical_orbitals = pyscf_scf.mo_coeff.astype(float)
    molecule.orbital_energies = pyscf_scf.mo_energy.astype(float)

    # Get integrals.
    one_body_integrals, two_body_integrals = compute_integrals(
        pyscf_molecule, pyscf_scf)
    molecule.one_body_integrals = one_body_integrals
    molecule.two_body_integrals = two_body_integrals
    molecule.overlap_integrals = pyscf_scf.get_ovlp()

    # Return updated molecule instance.
    pyscf_molecular_data = PyscfMolecularData.__new__(PyscfMolecularData)
    pyscf_molecular_data.__dict__.update(molecule.__dict__)
    pyscf_molecular_data.save()
    return pyscf_molecular_data

# -------------------------------------------------------------------------------------------------

def compute_Hartree_Fock(state):

    tstart = time.time()

    # -- hamiltonian is an InteractionOperator (defined in openfermion)
    # -- type of `molecule_data` is `MolecularData` (defined in openfermion)
    state.log_verbose("Construct instance of MolecularData.")
    join_str = '\n' + 11 * ' '
    geometry = state.config_hamiltonian.geometry()
    electron_count = 0
    for atom in geometry:
        electron_count += mendeleev.element(atom[0]).atomic_number
    print(f"electron_count = {electron_count}")
    geom_str = join_str.join(str(k) for k in geometry)
    state.log_verbose(f"geometry: {join_str}{geom_str}")
    molecule_data = MolecularData(
            geometry,
            state.config_hamiltonian.basis,
            1 + electron_count % 2,
            0, # TODO: is the charge really always 0?
            )

    # Generate a PyscfMolecularData
    # -- defined by openfermionpyscf, inherits from openfermion.chem.MolecularData
    state.log_verbose("Construct instance of PyscfMolecularData.")
    ham1_HartreeFock = my_run_pyscf(state, molecule_data)

    tstop = time.time()

    # Annotate with metadata
    ham1_HartreeFock.basis = state.config_hamiltonian.basis
    if len(geometry) == 2:
        c0 = geometry[0][1]
        c1 = geometry[1][1]
        ham1_HartreeFock.separation = math.sqrt(sum((c0[n] - c1[n])**2 for n in range(3)))
    else:
        ham1_HartreeFock.separation = None
    ham1_HartreeFock.hf_time = tstop - tstart

    state.log(f"Hartree-Fock energy from run_pyscf = {ham1_HartreeFock.hf_energy}")

    return ham1_HartreeFock

# -------------------------------------------------------------------------------------------------

def apply_active_space(state, ham1_HartreeFock):
    tstart = time.time()
    state.log_verbose("Define active space.")
    # -- the most useful statistics for users are in terms of single-occupancy (spin) orbitals
    # -- openfermion defines active space in terms of double-occupancy (spatial) orbitals
    idx_act_occ, idx_frz_vac, asmeta = define_active_orbitals(state, ham1_HartreeFock)

    # Build an InteractionOperator (from openfermion) with one- and two-body integrals computed
    state.log_verbose("Build instance of InteractionOperator.")
    ham2_ActiveSpace = ham1_HartreeFock.get_molecular_hamiltonian(
        occupied_indices=range(idx_act_occ),
        active_indices=range(idx_act_occ, idx_frz_vac),
    )
    tstop = time.time()
    # Pass along metadata that needs to be preserved
    ham2_ActiveSpace.hf_energy = ham1_HartreeFock.hf_energy
    ham2_ActiveSpace.asmeta = asmeta
    ham2_ActiveSpace.basis = ham1_HartreeFock.basis
    ham2_ActiveSpace.separation = ham1_HartreeFock.separation
    ham2_ActiveSpace.hf_time = ham1_HartreeFock.hf_time
    ham2_ActiveSpace.as_time = tstop - tstart
    return ham2_ActiveSpace

# -------------------------------------------------------------------------------------------------

def map_fermions_to_qubits(state, ham2_ActiveSpace):
    # Transform fermionic operators to qubit operators
    tstart = time.time()
    ham3_Fermion2Qubit = None
    mapping_name = state.config_hamiltonian.fermion_to_qubit_name()
    if mapping_name == "jordan-wigner":
        ham3_Fermion2Qubit = jordan_wigner(ham2_ActiveSpace)
    elif mapping_name == "bravyi-kitaev":
        ham3_Fermion2Qubit = bravyi_kitaev(ham2_ActiveSpace)
    else:
        mapping = state.config_hamiltonian.f2q_mapping()
        raise NotImplementedError(f"invalid fermion-to-qubit mapping \"{mapping}\"")
    tstop = time.time()
    # Propagate previously-computed metadata that we need to preserve across all run modes
    ham3_Fermion2Qubit.hf_energy = ham2_ActiveSpace.hf_energy
    ham3_Fermion2Qubit.asmeta = ham2_ActiveSpace.asmeta
    ham3_Fermion2Qubit.basis = ham2_ActiveSpace.basis
    ham3_Fermion2Qubit.separation = ham2_ActiveSpace.separation
    ham3_Fermion2Qubit.f2q_mapping = mapping_name
    ham3_Fermion2Qubit.hf_time = ham2_ActiveSpace.hf_time
    ham3_Fermion2Qubit.as_time = ham2_ActiveSpace.as_time
    ham3_Fermion2Qubit.f2q_time = tstop - tstart
    # Return result
    return ham3_Fermion2Qubit

# -------------------------------------------------------------------------------------------------

def write_data(state, pauli_sum, num_qubits):
    encoding = 'ascii'
    if state.config_general.file_format == "hamlib":
        # Build the string to match the HamLib format
        hamlib_str = bytearray()
        tail = "]".encode(encoding)
        space = " ".encode(encoding)[0]
        for string, coefficient in pauli_sum.items():
            w1 = 20
            w2 = w1 + (6 if coefficient >= 0 else 7)
            term = bytearray(f"({coefficient:{w2}.{w1}e}+0j) [", encoding)
            for qubit in string:
                term.extend(f"{qubit[1]}{qubit[0]} ".encode(encoding))
            if term[-1] == space:
                del term[-1:] # replace the extra trailing space with a closing brace
            term.extend(tail)
            hamlib_str.extend(term + " +\n".encode(encoding))
        del hamlib_str[-3:]
        # TODO: This needs to be saved to the appropriate data file (HDF5)
        with open(state.filename_ham3(), 'w') as f:
            for key, value in state.metadata.items():
                print(f"# {key} = {value}", file=f)
            print(hamlib_str.decode(), file=f)
    else:
        with open(state.filename_ham3(), 'w') as f:
            for key, value in state.metadata.items():
                print(f"# {key} = {value}", file=f)
            for string, coefficient in pauli_sum.items():
                s = bytearray('I' * num_qubits, encoding)
                for term in string:
                    s[term[0]] = term[1].encode(encoding)[0]
                print(f"{coefficient:24.17e} {s.decode()}", file=f)

# -------------------------------------------------------------------------------------------------

def get_ham1(state):
    ham1_filename = state.filename_ham1()
    try:
        with open(ham1_filename, 'rb') as file:
            ham1_HartreeFock = pickle.load(file)
    # TODO: Python errors get a little weird if you have an exception inside an exception.  So
    #       instead the exception clause should _only_ flag that we're going to recompute versus
    #       load, and the actual recomputing should be outside of the except clause.
    except FileNotFoundError as err:
        state.log(f"Could not load \"{ham1_filename}\".  Recomputing from the beginning.")
        # Compute ham1_ActiveSpace from scratch
        state.log("Perform Hartree-Fock calculation.")
        ham1_HartreeFock = compute_Hartree_Fock(state)
        state.log(f"Pickle to \"{ham1_filename}\" file.")
        # Save ham2_ActiveSpace for later re-use
        with open(ham1_filename, 'wb') as ham1_file:
            pickle.dump(ham1_HartreeFock, ham1_file)
        return ham1_HartreeFock
    else:
        state.log(' '.join([
            f"Loaded \"{ham1_filename}\".",
            "Continuing from after the Hartree-Fock calculation."]))
        return ham1_HartreeFock

# -------------------------------------------------------------------------------------------------

def get_ham2(state):
    ham2_filename = state.filename_ham2()
    state.log(f"Trying to load \"{ham2_filename}\".")
    try:
        with open(ham2_filename, 'rb') as file:
            ham2_ActiveSpace = pickle.load(file)
    # TODO: Python errors get a little weird if you have an exception inside an exception.  So
    #       instead the exception clause should _only_ flag that we're going to recompute versus
    #       load, and the actual recomputing should be outside of the except clause.
    except FileNotFoundError as err:
        state.log(' '.join([
            f"Could not load \"{ham2_filename}\".",
            f"Trying to load \"{state.filename_ham1()}\"."]))
        # Get ham1_HartreeFock (by loading or by recomputing, depending on data availability)
        ham1_HartreeFock = get_ham1(state)
        # Recompute ham2_ActiveSpace from ham1_HartreeFock
        state.log("Apply active space.")
        ham2_ActiveSpace = apply_active_space(state, ham1_HartreeFock)
        state.log(f"Pickle to \"{ham2_filename}\" file.")
        # Save ham2_ActiveSpace for later re-use
        with open(ham2_filename, 'wb') as ham2_file:
            pickle.dump(ham2_ActiveSpace, ham2_file)
        # Save the one-body and two-body tensors using numpy's `save` function.  These files can be
        # loaded using numpy's `load` function.
        t_filename = ham2_filename[:ham2_filename.rfind('.')] + ".tensors.npz"
        np.savez_compressed(t_filename,
                            one_body=ham2_ActiveSpace.one_body_tensor,
                            two_body=ham2_ActiveSpace.two_body_tensor)
        return ham2_ActiveSpace
    else:
        state.log(' '.join([
            f"Loaded \"{ham2_filename}\".",
            "Continuing from after the active space is applied."]))
        return ham2_ActiveSpace

# -------------------------------------------------------------------------------------------------

# TODO: Annotate all relevant quantities with units

def compute_metadata(state, ham3_Fermion2Qubit):
    # number of orbitals
    s_frz_occ = "number of frozen, occupied, single-occupancy orbitals"
    s_act_occ = "number of active, occupied, single-occupancy orbitals"
    s_act_vac = "number of active, vacant, single-occupancy orbitals"
    s_frz_vac = "number of frozen, vacant, single-occupancy orbitals"
    n_frz_occ = ham3_Fermion2Qubit.asmeta["n_frz_occ_so"]
    n_act_occ = ham3_Fermion2Qubit.asmeta["n_act_occ_so"]
    n_act_vac = ham3_Fermion2Qubit.asmeta["n_act_vac_so"]
    n_frz_vac = ham3_Fermion2Qubit.asmeta["n_frz_vac_so"]
    state.metadata[s_frz_occ] = n_frz_occ
    state.metadata[s_act_occ] = n_act_occ
    state.metadata[s_act_vac] = n_act_vac
    state.metadata[s_frz_vac] = n_frz_vac
    # size of active space
    state.metadata["size of active space"] = n_act_occ + n_act_vac
    # number of qubits
    state.metadata["number of qubits"] = n_act_occ + n_act_vac
    # basis set
    state.metadata["basis set"] = ham3_Fermion2Qubit.basis
    # interatomic spacing (in angstroms)
    state.metadata["interatomic separation (angstroms)"] = ham3_Fermion2Qubit.separation
    # fermion-to-qubit transformation
    state.metadata["fermion-to-qubit operator mapping"] = ham3_Fermion2Qubit.f2q_mapping
    # number of terms in sum of Pauli strings
    state.metadata["number of terms in sum of Pauli strings"] = len(ham3_Fermion2Qubit.terms)
    # one-norm of sum of Pauli strings
    one_norm = sum(abs(coefficient) for coefficient in ham3_Fermion2Qubit.terms.values())
    state.metadata["one-norm of sum of Pauli strings (Hartrees)"] = one_norm
    # eigenvalue bounds (see Trotter Workflow document for where this comes from)
    E0 = ham3_Fermion2Qubit.terms[tuple()]
    dE = one_norm - abs(E0)
    state.metadata["eigenvalue lower bound from shifted one-norm (Hartrees)"] = E0 - dE
    state.metadata["eigenvalue upper bound from shifted one-norm (Hartrees)"] = E0 + dE
    # Hartree-Fock energy
    state.metadata["Hartree-Fock energy (Hartrees)"] = ham3_Fermion2Qubit.hf_energy
    # pySCF runtime
    state.metadata["runtime for Hartree-Fock calculation (seconds)"] = ham3_Fermion2Qubit.hf_time
    state.metadata["runtime for active space calculation (seconds)"] = ham3_Fermion2Qubit.as_time
    state.metadata["runtime for fermion-to-qubit calculation (seconds)"] = ham3_Fermion2Qubit.f2q_time
    state.log("metadata:\n" + pprint.pformat(state.metadata, indent=3))

# -------------------------------------------------------------------------------------------------

def run():
    state = load_configuration()

    # TODO: If the final result files already exist, exit as no-op
    # TODO: Add the ability to save and reload ham3 similar to ham2 and ham1

    # Get ham2_ActiveSpace (by loading or by recomputing, depending on data availability)
    ham2_ActiveSpace = get_ham2(state)

    # Compute ham3_Fermion2Qubit from ham2_ActiveSpace
    state.log("Map fermionic operators to qubit operators.")
    ham3_Fermion2Qubit = map_fermions_to_qubits(state, ham2_ActiveSpace)

    state.log("Generate sum of Pauli strings.")
    pauli_sum = ham3_Fermion2Qubit.terms

    compute_metadata(state, ham3_Fermion2Qubit)

    ham3_filename = state.filename_ham3()
    state.log(f"Save sum of Pauli strings to data file \"{ham3_filename}\".")
    write_data(state, pauli_sum, ham2_ActiveSpace.n_qubits)

    state.log("Hamiltonian generation complete.")

# -------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    run()
