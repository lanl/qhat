import time
import numpy as np
import math
from pyLIQTR.PhaseEstimation.pe import PhaseEstimation
from pyLIQTR.utils.resource_analysis import estimate_resources

from openfermion.transforms import jordan_wigner

# --------------------------------------------------
# Helper functions for Pauli string arithmetic
# --------------------------------------------------
def pauli_dict(pauli_key):
    """Convert a tuple of (qubit, operator) pairs into a dictionary."""
    return {q: op for (q, op) in pauli_key}

def anticommute(key1, key2):
    """
    Return True if the Pauli strings (given by key1 and key2)
    anticommute. They anticommute if the number of qubits
    on which they both act nontrivially with different operators is odd.
    """
    d1 = pauli_dict(key1)
    d2 = pauli_dict(key2)
    count = 0
    for q in d1:
        if q in d2 and d1[q] != d2[q]:
            count += 1
    return (count % 2 == 1)

def pauli_product(p1, p2):
    """
    Multiply two single-qubit Pauli operators (ignoring phase).
    """
    if p1 == 'I':
        return p2
    if p2 == 'I':
        return p1
    if p1 == p2:
        return 'I'
    # For different Paulis, the unique remaining operator is returned.
    if {p1, p2} == {'X', 'Y'}:
        return 'Z'
    if {p1, p2} == {'X', 'Z'}:
        return 'Y'
    if {p1, p2} == {'Y', 'Z'}:
        return 'X'
    return 'I'

def pauli_product_key(key1, key2):
    """
    Given two Pauli strings (each represented as a tuple of (qubit, op)),
    return the Pauli string (as a sorted tuple) corresponding to their product.
    The product is computed qubit-by-qubit (ignoring any phase).
    """
    d1 = pauli_dict(key1)
    d2 = pauli_dict(key2)
    prod = {}
    for q in set(d1.keys()).union(d2.keys()):
        op1 = d1.get(q, 'I')
        op2 = d2.get(q, 'I')
        prod[q] = pauli_product(op1, op2)
    # Remove qubits where the product is identity.
    prod = {q: op for q, op in prod.items() if op != 'I'}
    return tuple(sorted(prod.items()))

# --------------------------------------------------
# Functions to compute (nested) commutator norms directly
# --------------------------------------------------
def compute_commutator(op1, op2):
    """
    Given two QubitOperator terms op1 and op2 (each with one Pauli string),
    compute their commutator [op1, op2].

    Returns:
      - (None, 0.0) if they commute;
      - ((prod_key, coeff), norm) if they anticommute,
        where prod_key is the resulting Pauli string (as a tuple)
        and coeff = 2*c1*c2 (with c1, c2 the term coefficients).
    The norm is computed as 2*|c1*c2|.
    """
    key1 = list(op1.terms.keys())[0]
    key2 = list(op2.terms.keys())[0]
    c1 = list(op1.terms.values())[0]
    c2 = list(op2.terms.values())[0]
    if not anticommute(key1, key2):
        return None, 0.0
    prod_key = pauli_product_key(key1, key2)
    coeff = 2 * c1 * c2
    norm_val = 2 * abs(c1 * c2)
    return (prod_key, coeff), norm_val

def compute_nested_commutator_norm(H_outer, inner_operator):
    """
    Given an outer QubitOperator H_outer and an inner operator (a tuple (prod_key, coeff))
    representing the inner commutator [H_j, H_k], compute the norm of the nested commutator:
        [H_outer, inner_operator].
    If the Pauli string of H_outer anticommutes with inner_operator's pauli key,
    the norm is 2*|c_outer * (inner coeff)|; otherwise, it is zero.
    """
    if inner_operator is None:
        return 0.0
    key_outer = list(H_outer.terms.keys())[0]
    c_outer = list(H_outer.terms.values())[0]
    if not anticommute(key_outer, inner_operator[0]):
        return 0.0
    return 2 * abs(c_outer * inner_operator[1])

# --------------------------------------------------
# Main Monte Carlo estimation
# --------------------------------------------------
def trotter_error_estimator(pauli_terms, time_limit, batch_size=100):
    """
    Monte Carlo estimation of Trotter error coefficients.

    Reference: Childs et al., "Theory of Trotter Error" (arXiv:1912.08854v3)

    For first-order formulas (Equation 145):
      Error ≤ (t²/2) * C1, where C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||

    For second-order formulas (Equation 152):
      Error ≤ (t³/12) * C21 + (t³/24) * C22, where
      C21 = Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]||  (all triples with k < i and k < j)
      C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||       (all pairs with k < j)

    Returns:
      (C1, C2) where C2 = C21/12 + C22/24
    """
    N = len(pauli_terms)

    # Dictionary to cache computed pair commutators.
    # Keys are sorted tuples (i, j) and values are (operator, norm).
    comm_dict = {}

    # ---------------------------
    # Estimate C1 = sum_{i<j} ||[H_i, H_j]||
    # ---------------------------
    C1_sum = 0.0
    samples_C1 = 0
    start_time = time.time()
    while time.time() - start_time < time_limit/3:
        for _ in range(batch_size):
            # Sample a pair (i, j) with i != j.
            i = np.random.randint(0, N)
            j = np.random.randint(0, N - 1)
            if j >= i:
                j += 1
            op, norm_val = compute_commutator(pauli_terms[i], pauli_terms[j])
            # Cache the computed commutator.
            key = tuple(sorted((i, j)))
            comm_dict[key] = (op, norm_val)
            C1_sum += norm_val
            samples_C1 += 1
        if time.time() - start_time >= time_limit:
            break
    total_C1 = N * (N - 1)/2  # Total number of distinct pairs (i < j)
    C1_est = C1_sum * (total_C1 / samples_C1) if samples_C1 > 0 else 0.0
    print(f"  C1: {samples_C1} samples")

    # ---------------------------
    # Estimate C21 = sum_{k<j, k<i} ||[H_i, [H_j, H_k]]||
    # ---------------------------
    C21_sum = 0.0
    samples_C21 = 0
    start_time = time.time()
    while time.time() - start_time < time_limit/3:
        for _ in range(batch_size):
            # Choose k in [0, N-2] and then sample two distinct indices i and j from {k+1, ..., N-1}
            k = np.random.randint(0, N - 1)
            valid = np.arange(k + 1, N)
            if len(valid) < 2:
                continue
            i, j = np.random.choice(valid, size=2, replace=False)
            # Inner commutator: [H_j, H_k]
            key_inner = tuple(sorted((j, k)))
            if key_inner in comm_dict:
                inner_op, _ = comm_dict[key_inner]
            else:
                inner_op, _ = compute_commutator(pauli_terms[j], pauli_terms[k])
                comm_dict[key_inner] = (inner_op, 0 if inner_op is None else 2 * abs(list(pauli_terms[j].terms.values())[0] *
                                                                                      list(pauli_terms[k].terms.values())[0]))
            # Compute nested commutator norm: [H_i, [H_j, H_k]]
            nested_norm = compute_nested_commutator_norm(pauli_terms[i], inner_op)
            C21_sum += nested_norm
            samples_C21 += 1
        if time.time() - start_time >= time_limit:
            break
    # Total number of valid triples: for each k, there are C(N-k-1, 2) choices of (i,j).
    total_C21 = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
    C21_est = C21_sum * (total_C21 / samples_C21) if samples_C21 > 0 else 0.0
    print(f"  C21: {samples_C21} samples")

    # ---------------------------
    # Estimate C22 = sum_{k<j} ||[H_k, [H_k, H_j]]||
    # ---------------------------
    C22_sum = 0.0
    samples_C22 = 0
    start_time = time.time()
    while time.time() - start_time < time_limit/3:
        for _ in range(batch_size):
            # Sample a pair (k, j) with k < j.
            k = np.random.randint(0, N - 1)
            j = np.random.randint(k + 1, N)
            key_inner = tuple(sorted((k, j)))
            if key_inner in comm_dict:
                inner_op, _ = comm_dict[key_inner]
            else:
                inner_op, _ = compute_commutator(pauli_terms[k], pauli_terms[j])
                comm_dict[key_inner] = (inner_op, 0 if inner_op is None else 2 * abs(list(pauli_terms[k].terms.values())[0] *
                                                                                      list(pauli_terms[j].terms.values())[0]))
            # Compute nested commutator norm: [H_k, [H_k, H_j]]
            nested_norm = compute_nested_commutator_norm(pauli_terms[k], inner_op)
            C22_sum += nested_norm
            samples_C22 += 1
        if time.time() - start_time >= time_limit:
            break
    total_C22 = N * (N - 1) / 2  # Total number of (k, j) pairs with k < j
    C22_est = C22_sum * (total_C22 / samples_C22) if samples_C22 > 0 else 0.0
    print(f"  C22: {samples_C22} samples")

    # ---------------------------
    # Final output
    # ---------------------------
    return C1_est/2, C21_est/12 + C22_est/24

def build_active_space(molecule, n_active_electrons_per_atom, n_active_unocc_orbitals_per_atom):
    n_atoms = len(molecule.geometry)
    n_electrons_total = molecule.n_electrons

    n_active_electrons = n_active_electrons_per_atom * n_atoms
    n_active_DE_occ_orbs = n_active_electrons//2
    n_active_DE_unocc_orbs = (n_active_unocc_orbitals_per_atom * n_atoms)//2

    n_DE_occupied = n_electrons_total//2

    frozen_occupied = n_DE_occupied - n_active_DE_occ_orbs
    occupied_indices = range(frozen_occupied)

    active_occ_start = frozen_occupied
    active_occ_stop  = n_DE_occupied
    active_vir_start = n_DE_occupied
    active_vir_stop  = n_DE_occupied + n_active_DE_unocc_orbs

    active_indices = list(range(active_occ_start, active_occ_stop)) \
                + list(range(active_vir_start, active_vir_stop))

    molecular_hamiltonian = molecule.get_molecular_hamiltonian(
        occupied_indices=occupied_indices,
        active_indices=active_indices
    )

    n_active_unocc_total = n_active_unocc_orbitals_per_atom * n_atoms

    init_state = [0]*n_active_unocc_total + [1]*n_active_electrons

    return molecular_hamiltonian, init_state


def trotter_resource_estimator(active_hamiltonian, init_state, trot_ord):
    # things we care less about right now
    trotterize = True
    ev_time = 1 # this has no effect on the resource estimate
    trot_num = 1 # just do a single trotter step
    precision_order = 1 # and a single bit of precision
    # NOTE: skipping phase_offset for now

    gse_args = {
        'mol_ham' : active_hamiltonian,
        'trotterize': trotterize,
        'ev_time': ev_time,
        'trot_ord': trot_ord,
        'trot_num': trot_num,
    }

    gse_inst = PhaseEstimation(
        precision_order=precision_order,
        init_state=init_state,
        kwargs=gse_args
    )

    gse_inst.generate_circuit()
    gse_circuit = gse_inst.pe_circuit

    # calculate resource estimate with pyliqtr/qualtran

    re_qualtran = estimate_resources(gse_circuit)

    return re_qualtran['T']

def generate_resource_estimate(molecule,
                               n_active_electrons_per_atom,
                               n_active_unocc_orbitals_per_atom,
                               epsilon_per_atom,
                               output_qubits,
                               trotter_order,
                               trotter_error_runtime):

    active_space_size = (n_active_electrons_per_atom + n_active_unocc_orbitals_per_atom)*molecule.n_atoms

    print("\n-----------------------------")
    print("Starting resource estimation with active space of size", active_space_size)

    # --------------------------------------
    # Step 1: Restrict to the Active Space
    # --------------------------------------
    print("\nStep 1: Restricting to the active space...")
    start_time = time.time()

    active_hamiltonian, init_state = build_active_space(
        molecule, n_active_electrons_per_atom, n_active_unocc_orbitals_per_atom
    )

    end_time = time.time()
    print(f"    ...finished in {end_time - start_time:.2f} seconds.")

    # ------------------------------------------------------
    # Step 2: Convert to Pauli strings via JW Transform
    # ------------------------------------------------------
    print("\nStep 2: Applying Jordan-Wigner transformation...")
    start_time = time.time()

    H = list(jordan_wigner(active_hamiltonian))

    end_time = time.time()
    print(f"    ...finished in {end_time - start_time:.2f} seconds, number of Pauli strings: {len(H)}")

    # ----------------------------------------------
    # Step 3: Estimate Trotter Error Coefficients
    # ----------------------------------------------
    print("\nStep 3: Estimating Trotter error coefficients...")
    start_time = time.time()

    _, c2 = trotter_error_estimator(H, trotter_error_runtime)

    end_time = time.time()
    print(f"    ...finished in {end_time - start_time:.2f} seconds, estimated coefficient: C = {c2:.6f}")

    # ----------------------------------------------
    # Step 4: Estimate Trotter T Count
    # ----------------------------------------------
    print("\nStep 4: Estimating Trotter T count...")
    start_time = time.time()

    trotter_Tcount = trotter_resource_estimator(active_hamiltonian, init_state, trotter_order)

    end_time = time.time()
    print(f"    ...finished in {end_time - start_time:.2f} seconds, estimated T count per step: {trotter_Tcount:.2e}")

    # ----------------------------------------------
    # Step 5: Compute Final Resource Estimation
    # ----------------------------------------------

    start_time = time.time()

    epsilon = epsilon_per_atom * molecule.n_atoms
    m = output_qubits

    t = np.pi / ((epsilon / 2) * (2 ** (m - 1)))  # for U = e^{-iHt}
    trotter_steps = (2 * np.pi * (2**m - 1) * c2 * (t**2) / (epsilon / 2)) ** 0.5
    total_Tcount = (2**m - 1) * trotter_steps * trotter_Tcount
    total_qubits = active_space_size + output_qubits

    # ----------------------------------------------
    # Final Output
    # ----------------------------------------------

    print("\nFinal result:")
    print(f"    {total_qubits} qubits, {total_Tcount:.2e} T gates")

    return total_qubits, total_Tcount
