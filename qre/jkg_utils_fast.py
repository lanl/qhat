"""
Optimized version of jkg_utils.py for faster Monte Carlo estimation of nested commutators.

Key optimizations:
1. Binary encoding of Pauli strings (2 bits per qubit: 00=I, 01=X, 10=Y, 11=Z)
2. NumPy vectorization for batch operations
3. Numba JIT compilation for critical functions
4. Bitwise operations for anticommutation checks
"""

import time
import numpy as np
import math
from numba import njit, prange

# --------------------------------------------------
# Binary encoding for Pauli strings
# --------------------------------------------------
# Each Pauli operator is encoded as 2 bits:
#   I = 00 (0)
#   X = 01 (1)
#   Y = 10 (2)
#   Z = 11 (3)

PAULI_I = 0
PAULI_X = 1
PAULI_Y = 2
PAULI_Z = 3

def encode_pauli_string(pauli_key, n_qubits):
    """
    Encode a Pauli string (tuple of (qubit, op) pairs) into two integer arrays.

    Args:
        pauli_key: tuple of (qubit, operator) pairs
        n_qubits: total number of qubits

    Returns:
        x_bits, z_bits: Two bit arrays representing the Pauli string
        (X = (1,0), Y = (1,1), Z = (0,1), I = (0,0))
    """
    x_bits = 0
    z_bits = 0

    for qubit, op in pauli_key:
        if op == 'X':
            x_bits |= (1 << qubit)
        elif op == 'Y':
            x_bits |= (1 << qubit)
            z_bits |= (1 << qubit)
        elif op == 'Z':
            z_bits |= (1 << qubit)

    return x_bits, z_bits


@njit
def pauli_anticommute(x1, z1, x2, z2):
    """
    Check if two Pauli strings anticommute using bitwise operations.

    Two Pauli strings anticommute if the number of positions where
    they both have non-identity, non-equal operators is odd.

    Args:
        x1, z1: X and Z bits for first Pauli string
        x2, z2: X and Z bits for second Pauli string

    Returns:
        True if they anticommute, False otherwise
    """
    # Operators differ if either x-bits or z-bits differ
    # Both are non-identity if (x|z) != 0
    diff = ((x1 ^ x2) | (z1 ^ z2)) & (x1 | z1) & (x2 | z2)

    # Count bits set in diff
    count = 0
    while diff:
        count += diff & 1
        diff >>= 1

    return (count % 2) == 1


@njit
def compute_commutator_norm_fast(x1, z1, c1, x2, z2, c2):
    """
    Compute the norm of commutator [P1, P2] where P1, P2 are Pauli strings.

    Returns 2*|c1*c2| if they anticommute, 0 otherwise.
    """
    if not pauli_anticommute(x1, z1, x2, z2):
        return 0.0
    return 2.0 * abs(c1 * c2)


@njit
def pauli_product_bits(x1, z1, x2, z2):
    """
    Compute the product of two Pauli strings using bitwise operations.
    Returns (x_prod, z_prod) representing the product Pauli string.
    Ignores phase factors.
    """
    # X * X = I, Z * Z = I, X * Z = Y, Z * X = Y, etc.
    x_prod = x1 ^ x2
    z_prod = z1 ^ z2

    return x_prod, z_prod


@njit
def compute_nested_commutator_norm_fast(x_outer, z_outer, c_outer,
                                        x_inner, z_inner, c_inner):
    """
    Compute norm of nested commutator [P_outer, P_inner].

    Returns 2*|c_outer * c_inner| if they anticommute, 0 otherwise.
    """
    if not pauli_anticommute(x_outer, z_outer, x_inner, z_inner):
        return 0.0
    return 2.0 * abs(c_outer * c_inner)


@njit(parallel=True)
def batch_compute_C1(x_bits, z_bits, coeffs, indices, N):
    """
    Compute C1 samples in parallel for multiple random pairs.

    Args:
        x_bits, z_bits: Arrays of X and Z bits for all Pauli strings
        coeffs: Array of coefficients for all Pauli strings
        indices: Array of (i, j) pairs to sample
        N: number of Pauli strings

    Returns:
        Array of norm values for each pair
    """
    n_samples = len(indices)
    norms = np.zeros(n_samples)

    for idx in prange(n_samples):
        i = indices[idx, 0]
        j = indices[idx, 1]
        norms[idx] = compute_commutator_norm_fast(
            x_bits[i], z_bits[i], coeffs[i],
            x_bits[j], z_bits[j], coeffs[j]
        )

    return norms


@njit(parallel=True)
def batch_compute_C21(x_bits, z_bits, coeffs, indices, N):
    """
    Compute C21 samples in parallel for multiple random triples (i, j, k).

    Returns:
        Array of nested commutator norms for [H_i, [H_j, H_k]]
    """
    n_samples = len(indices)
    norms = np.zeros(n_samples)

    for idx in prange(n_samples):
        i = indices[idx, 0]
        j = indices[idx, 1]
        k = indices[idx, 2]

        # First compute inner commutator [H_j, H_k]
        if not pauli_anticommute(x_bits[j], z_bits[j], x_bits[k], z_bits[k]):
            norms[idx] = 0.0
            continue

        # Product of H_j and H_k
        x_inner, z_inner = pauli_product_bits(x_bits[j], z_bits[j],
                                              x_bits[k], z_bits[k])
        c_inner = 2.0 * coeffs[j] * coeffs[k]

        # Now compute nested commutator norm
        norms[idx] = compute_nested_commutator_norm_fast(
            x_bits[i], z_bits[i], coeffs[i],
            x_inner, z_inner, c_inner
        )

    return norms


@njit(parallel=True)
def batch_compute_C22(x_bits, z_bits, coeffs, indices, N):
    """
    Compute C22 samples in parallel for multiple random pairs (k, j).

    Returns:
        Array of nested commutator norms for [H_k, [H_k, H_j]]
    """
    n_samples = len(indices)
    norms = np.zeros(n_samples)

    for idx in prange(n_samples):
        k = indices[idx, 0]
        j = indices[idx, 1]

        # First compute inner commutator [H_k, H_j]
        if not pauli_anticommute(x_bits[k], z_bits[k], x_bits[j], z_bits[j]):
            norms[idx] = 0.0
            continue

        # Product of H_k and H_j
        x_inner, z_inner = pauli_product_bits(x_bits[k], z_bits[k],
                                              x_bits[j], z_bits[j])
        c_inner = 2.0 * coeffs[k] * coeffs[j]

        # Now compute nested commutator norm with H_k
        norms[idx] = compute_nested_commutator_norm_fast(
            x_bits[k], z_bits[k], coeffs[k],
            x_inner, z_inner, c_inner
        )

    return norms


def preprocess_pauli_terms(pauli_terms):
    """
    Convert QubitOperator terms to efficient binary representation.

    Args:
        pauli_terms: List of QubitOperator terms

    Returns:
        x_bits, z_bits, coeffs, n_qubits
    """
    N = len(pauli_terms)

    # Find maximum qubit index
    n_qubits = 0
    for term in pauli_terms:
        key = list(term.terms.keys())[0]
        if key:  # non-empty
            max_qubit = max(q for q, _ in key)
            n_qubits = max(n_qubits, max_qubit + 1)

    # Encode all Pauli strings
    x_bits = np.zeros(N, dtype=np.int64)
    z_bits = np.zeros(N, dtype=np.int64)
    coeffs = np.zeros(N, dtype=np.complex128)

    for i, term in enumerate(pauli_terms):
        key = list(term.terms.keys())[0]
        coeff = list(term.terms.values())[0]

        x_bits[i], z_bits[i] = encode_pauli_string(key, n_qubits)
        coeffs[i] = coeff

    return x_bits, z_bits, coeffs, n_qubits


def trotter_error_estimator_fast(pauli_terms, time_limit, batch_size=10000):
    """
    Fast Monte Carlo estimation of nested commutator norms.

    This is an optimized version of trotter_error_estimator that uses:
    - Binary encoding for Pauli strings
    - Numba JIT compilation
    - Vectorized batch operations
    - Parallel processing

    Args:
        pauli_terms: List of QubitOperator terms
        time_limit: Total time limit in seconds
        batch_size: Number of samples per batch (larger = better parallelization)

    Returns:
        (C1_est, C2_est): Estimated first and second order commutator norms
    """
    N = len(pauli_terms)

    print(f"Preprocessing {N} Pauli terms...")
    start_prep = time.time()
    x_bits, z_bits, coeffs, n_qubits = preprocess_pauli_terms(pauli_terms)
    print(f"  Preprocessing done in {time.time() - start_prep:.3f}s ({n_qubits} qubits)")

    # ---------------------------
    # Estimate C1 = sum_{i<j} ||[H_i, H_j]||
    # ---------------------------
    print(f"Estimating C1 with batch_size={batch_size}...")
    C1_sum = 0.0
    samples_C1 = 0
    start_time = time.time()

    while time.time() - start_time < time_limit / 3:
        # Generate random pairs in batch
        i_vals = np.random.randint(0, N, size=batch_size)
        j_vals = np.random.randint(0, N - 1, size=batch_size)
        j_vals = np.where(j_vals >= i_vals, j_vals + 1, j_vals)

        indices = np.column_stack([i_vals, j_vals])

        # Compute norms in parallel
        norms = batch_compute_C1(x_bits, z_bits, coeffs, indices, N)

        C1_sum += np.sum(norms)
        samples_C1 += batch_size

        if time.time() - start_time >= time_limit / 3:
            break

    total_C1 = N * (N - 1) / 2
    C1_est = C1_sum * (total_C1 / samples_C1) if samples_C1 > 0 else 0.0
    print(f"  C1 estimation: {samples_C1} samples in {time.time() - start_time:.3f}s")
    print(f"  C1 estimate: {C1_est:.6f}")

    # ---------------------------
    # Estimate C21 = sum_{k<j, k<i} ||[H_i, [H_j, H_k]]||
    # ---------------------------
    print(f"Estimating C21 with batch_size={batch_size}...")
    C21_sum = 0.0
    samples_C21 = 0
    start_time = time.time()

    while time.time() - start_time < time_limit / 3:
        # Generate random triples (i, j, k) with k < i and k < j
        k_vals = np.random.randint(0, max(1, N - 1), size=batch_size)

        # For each k, sample two distinct indices from {k+1, ..., N-1}
        valid_samples = []
        for k in k_vals:
            valid_range = N - k - 1
            if valid_range < 2:
                continue
            # Sample two distinct indices
            choices = np.random.choice(valid_range, size=2, replace=False)
            i, j = k + 1 + choices[0], k + 1 + choices[1]
            valid_samples.append([i, j, k])

        if not valid_samples:
            continue

        indices = np.array(valid_samples, dtype=np.int64)

        # Compute nested norms in parallel
        norms = batch_compute_C21(x_bits, z_bits, coeffs, indices, N)

        C21_sum += np.sum(norms)
        samples_C21 += len(indices)

        if time.time() - start_time >= time_limit / 3:
            break

    total_C21 = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
    C21_est = C21_sum * (total_C21 / samples_C21) if samples_C21 > 0 else 0.0
    print(f"  C21 estimation: {samples_C21} samples in {time.time() - start_time:.3f}s")
    print(f"  C21 estimate: {C21_est:.6f}")

    # ---------------------------
    # Estimate C22 = sum_{k<j} ||[H_k, [H_k, H_j]]||
    # ---------------------------
    print(f"Estimating C22 with batch_size={batch_size}...")
    C22_sum = 0.0
    samples_C22 = 0
    start_time = time.time()

    while time.time() - start_time < time_limit / 3:
        # Generate random pairs (k, j) with k < j
        k_vals = np.random.randint(0, max(1, N - 1), size=batch_size)
        j_vals = np.array([np.random.randint(k + 1, N) for k in k_vals], dtype=np.int64)

        indices = np.column_stack([k_vals, j_vals])

        # Compute nested norms in parallel
        norms = batch_compute_C22(x_bits, z_bits, coeffs, indices, N)

        C22_sum += np.sum(norms)
        samples_C22 += batch_size

        if time.time() - start_time >= time_limit / 3:
            break

    total_C22 = N * (N - 1) / 2
    C22_est = C22_sum * (total_C22 / samples_C22) if samples_C22 > 0 else 0.0
    print(f"  C22 estimation: {samples_C22} samples in {time.time() - start_time:.3f}s")
    print(f"  C22 estimate: {C22_est:.6f}")

    # ---------------------------
    # Final output
    # ---------------------------
    return C1_est, C21_est / 12 + C22_est / 24


def generate_resource_estimate_fast(molecule,
                                    n_active_electrons_per_atom,
                                    n_active_unocc_orbitals_per_atom,
                                    epsilon_per_atom,
                                    output_qubits,
                                    trotter_order,
                                    trotter_error_runtime):
    """
    Fast version of generate_resource_estimate using optimized commutator estimation.

    This function is a drop-in replacement for generate_resource_estimate
    from jkg_utils.py, but uses the optimized trotter_error_estimator_fast.
    """
    from openfermion.transforms import jordan_wigner
    from .jkg_utils import build_active_space, trotter_resource_estimator

    active_space_size = (n_active_electrons_per_atom + n_active_unocc_orbitals_per_atom) * molecule.n_atoms

    print("\n-----------------------------")
    print("Starting FAST resource estimation with active space of size", active_space_size)

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
    # Step 3: Estimate Trotter Error Coefficients (FAST VERSION)
    # ----------------------------------------------
    print("\nStep 3: Estimating Trotter error coefficients (FAST)...")
    start_time = time.time()

    _, c2 = trotter_error_estimator_fast(H, trotter_error_runtime)

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
