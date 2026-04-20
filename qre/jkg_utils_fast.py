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
# Configuration for exact computation
# --------------------------------------------------
# Time budget for exact computation feasibility check (seconds)
EXACT_COMPUTATION_TIME_BUDGET = 60

# Memory limit for exact computation tracking (MB)
EXACT_COMPUTATION_MEMORY_LIMIT_MB = 100

# --------------------------------------------------
# Configuration for convergence monitoring
# --------------------------------------------------
# Enable convergence monitoring and reporting
ENABLE_CONVERGENCE_MONITORING = True

# Relative standard error threshold for early termination (e.g., 0.01 = 1%)
# Set to 0 to disable early termination based on convergence
CONVERGENCE_THRESHOLD = 0.01  # Terminate when SE/mean < 1%

# Minimum number of samples before checking convergence
MIN_SAMPLES_FOR_CONVERGENCE = 100000

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


def should_use_exact_tracking(N, time_budget=None, memory_limit_mb=None):
    """
    Decide whether to use exact computation tracking based on time and memory constraints.

    Args:
        N: Number of Pauli terms
        time_budget: Time budget in seconds (defaults to EXACT_COMPUTATION_TIME_BUDGET)
        memory_limit_mb: Memory limit in MB (defaults to EXACT_COMPUTATION_MEMORY_LIMIT_MB)

    Returns:
        bool: True if exact tracking is feasible and beneficial
    """
    if time_budget is None:
        time_budget = EXACT_COMPUTATION_TIME_BUDGET
    if memory_limit_mb is None:
        memory_limit_mb = EXACT_COMPUTATION_MEMORY_LIMIT_MB

    # Compute combination counts
    c1_count = N * (N - 1) // 2

    # For C21, use approximation for large N to avoid expensive computation
    if N < 500:
        c21_count = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
    else:
        c21_count = N**3 // 6  # Approximation

    c22_count = N * (N - 1) // 2

    # Memory check (12 bytes per triple for set overhead, 8 per pair)
    estimated_memory_mb = (c1_count * 8 + c21_count * 12 + c22_count * 8) / (1024**2)
    if estimated_memory_mb > memory_limit_mb:
        return False

    # Time check (coupon collector: N*ln(N) samples needed on average)
    # Use throughput estimates: 15M/s for C1, 8M/s for C21, 11M/s for C22
    c1_time = c1_count * np.log(max(2, c1_count)) / 15e6 if c1_count > 0 else 0
    c21_time = c21_count * np.log(max(2, c21_count)) / 8e6 if c21_count > 0 else 0
    c22_time = c22_count * np.log(max(2, c22_count)) / 11e6 if c22_count > 0 else 0
    estimated_time = c1_time + c21_time + c22_time

    if estimated_time > time_budget:
        return False

    return True


def trotter_error_estimator_fast(pauli_terms, time_limit, batch_size=10000):
    """
    Fast Monte Carlo estimation of nested commutator norms with exact computation for small systems.

    Reference: Childs et al., "Theory of Trotter Error" (arXiv:1912.08854v3)

    This is an optimized version of trotter_error_estimator that uses:
    - Binary encoding for Pauli strings
    - Numba JIT compilation
    - Vectorized batch operations
    - Parallel processing
    - Early termination with exact results when all combinations are sampled

    For first-order formulas (Equation 145):
      Error ≤ (t²/2) * C1, where C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||

    For second-order formulas (Equation 152):
      Error ≤ (t³/12) * C21 + (t³/24) * C22, where
      C21 = Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]||  (all triples with k < i and k < j)
      C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||       (all pairs with k < j)

    Args:
        pauli_terms: List of QubitOperator terms
        time_limit: Total time limit in seconds
        batch_size: Number of samples per batch (larger = better parallelization)

    Returns:
        (C1_est, C2_est): Estimated first and second order commutator norms
        where C1_est = C1/2 and C2_est = C21/12 + C22/24
    """
    N = len(pauli_terms)

    print(f"Preprocessing {N} Pauli terms...")
    start_prep = time.time()
    x_bits, z_bits, coeffs, n_qubits = preprocess_pauli_terms(pauli_terms)
    print(f"  Preprocessing done in {time.time() - start_prep:.3f}s ({n_qubits} qubits)")

    # Check if exact computation is feasible
    use_exact = should_use_exact_tracking(N, time_limit)
    if use_exact:
        print(f"  Exact computation is feasible for N={N} - enabling tracking")
        seen_c1 = set()
        seen_c21 = set()
        seen_c22 = set()
        c1_values = {}
        c21_values = {}
        c22_values = {}

        # Compute total combinations for completion check
        total_c1 = N * (N - 1) // 2
        total_c21 = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
        total_c22 = N * (N - 1) // 2
    else:
        print(f"  Using Monte Carlo estimation (N={N} too large for exact computation)")

    # Warmup: trigger Numba JIT compilation before timing
    print(f"Warming up Numba JIT compilation...")
    if N >= 2:
        dummy_indices = np.array([[0, 1]], dtype=np.int64)
        batch_compute_C1(x_bits, z_bits, coeffs, dummy_indices, N)
        if N >= 3:
            dummy_indices_3 = np.array([[0, 1, 2]], dtype=np.int64)
            batch_compute_C21(x_bits, z_bits, coeffs, dummy_indices_3, N)
            batch_compute_C22(x_bits, z_bits, coeffs, dummy_indices[:, [0, 1]], N)
    print(f"  Warmup complete")

    # ---------------------------
    # Estimate C1 = sum_{i<j} ||[H_i, H_j]||
    # ---------------------------
    print(f"Estimating C1 with batch_size={batch_size}...")
    C1_sum = 0.0
    C1_sum_sq = 0.0  # For variance calculation
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
        C1_sum_sq += np.sum(norms**2)
        samples_C1 += batch_size

        # Check convergence if enabled
        if ENABLE_CONVERGENCE_MONITORING and samples_C1 >= MIN_SAMPLES_FOR_CONVERGENCE:
            total_C1_count = N * (N - 1) / 2
            C1_est_current = C1_sum * (total_C1_count / samples_C1)
            # Standard error: SE = std(samples) * sqrt(N_total / N_samples) / sqrt(N_samples)
            sample_mean = C1_sum / samples_C1
            sample_var = (C1_sum_sq / samples_C1) - sample_mean**2
            if sample_var > 0:
                sample_std = np.sqrt(sample_var)
                se = sample_std * np.sqrt(total_C1_count / samples_C1) / np.sqrt(samples_C1)
                rel_se = se / C1_est_current if C1_est_current > 0 else float('inf')

                # Check if converged
                if CONVERGENCE_THRESHOLD > 0 and rel_se < CONVERGENCE_THRESHOLD:
                    print(f"  C1 CONVERGED: SE/mean = {rel_se:.4f} < {CONVERGENCE_THRESHOLD:.4f}")
                    print(f"  C1 estimate: {C1_est_current:.6f} ± {se:.6f} (95% CI)")
                    C1_est = C1_est_current
                    break

        # Track exact values if enabled
        if use_exact:
            for idx in range(len(indices)):
                i, j = int(indices[idx, 0]), int(indices[idx, 1])
                pair = (min(i, j), max(i, j))
                if pair not in seen_c1:
                    seen_c1.add(pair)
                    c1_values[pair] = norms[idx]

            # Check if we've seen all pairs
            if len(seen_c1) == total_c1:
                C1_exact = sum(c1_values.values())
                print(f"  C1 EXACT: All {total_c1} pairs sampled in {time.time() - start_time:.3f}s")
                print(f"  C1 exact value: {C1_exact:.6f}")
                C1_est = C1_exact
                break

        if time.time() - start_time >= time_limit / 3:
            break

    if not use_exact or len(seen_c1) < total_c1:
        total_C1 = N * (N - 1) / 2
        C1_est = C1_sum * (total_C1 / samples_C1) if samples_C1 > 0 else 0.0
        print(f"  C1 estimation: {samples_C1} samples in {time.time() - start_time:.3f}s")

        # Report standard error if convergence monitoring enabled
        if ENABLE_CONVERGENCE_MONITORING and samples_C1 > 0:
            sample_mean = C1_sum / samples_C1
            sample_var = (C1_sum_sq / samples_C1) - sample_mean**2
            if sample_var > 0:
                sample_std = np.sqrt(sample_var)
                se = sample_std * np.sqrt(total_C1 / samples_C1) / np.sqrt(samples_C1)
                rel_se = se / C1_est if C1_est > 0 else float('inf')
                print(f"  C1 estimate: {C1_est:.6f} ± {se:.6f} (rel. SE: {rel_se:.4f})")
            else:
                print(f"  C1 estimate: {C1_est:.6f}")
        else:
            print(f"  C1 estimate: {C1_est:.6f}")

        if use_exact:
            print(f"    (Sampled {len(seen_c1)}/{total_c1} unique pairs, {100*len(seen_c1)/total_c1:.1f}% coverage)")

    # ---------------------------
    # Estimate C21 = sum_{k<j, k<i} ||[H_i, [H_j, H_k]]||
    # ---------------------------
    print(f"Estimating C21 with batch_size={batch_size}...")
    C21_sum = 0.0
    C21_sum_sq = 0.0  # For variance calculation
    samples_C21 = 0
    start_time = time.time()

    while time.time() - start_time < time_limit / 3:
        # Generate random triples (i, j, k) with k < i and k < j
        # Only sample k values that allow at least 2 choices above them
        k_vals = np.random.randint(0, max(1, N - 2), size=batch_size)

        # For each k, sample two distinct indices from {k+1, ..., N-1}
        # Vectorized approach: sample i and j independently, then ensure distinctness
        valid_ranges = N - k_vals - 1

        # Sample first index i from [k+1, N)
        i_offsets = np.random.randint(0, valid_ranges, dtype=np.int64)
        i_vals = k_vals + 1 + i_offsets

        # Sample second index j from [k+1, N) \ {i}
        # Use rejection sampling trick: sample from [0, valid_range-1), then adjust
        j_offsets = np.random.randint(0, valid_ranges - 1, dtype=np.int64)
        j_offsets = np.where(j_offsets >= i_offsets, j_offsets + 1, j_offsets)
        j_vals = k_vals + 1 + j_offsets

        indices = np.column_stack([i_vals, j_vals, k_vals])

        # Compute nested norms in parallel
        norms = batch_compute_C21(x_bits, z_bits, coeffs, indices, N)

        C21_sum += np.sum(norms)
        C21_sum_sq += np.sum(norms**2)
        samples_C21 += len(indices)

        # Check convergence if enabled
        if ENABLE_CONVERGENCE_MONITORING and samples_C21 >= MIN_SAMPLES_FOR_CONVERGENCE:
            total_C21_count = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
            C21_est_current = C21_sum * (total_C21_count / samples_C21)
            sample_mean = C21_sum / samples_C21
            sample_var = (C21_sum_sq / samples_C21) - sample_mean**2
            if sample_var > 0:
                sample_std = np.sqrt(sample_var)
                se = sample_std * np.sqrt(total_C21_count / samples_C21) / np.sqrt(samples_C21)
                rel_se = se / C21_est_current if C21_est_current > 0 else float('inf')

                if CONVERGENCE_THRESHOLD > 0 and rel_se < CONVERGENCE_THRESHOLD:
                    print(f"  C21 CONVERGED: SE/mean = {rel_se:.4f} < {CONVERGENCE_THRESHOLD:.4f}")
                    print(f"  C21 estimate: {C21_est_current:.6f} ± {se:.6f} (95% CI)")
                    C21_est = C21_est_current
                    break

        # Track exact values if enabled
        if use_exact:
            for idx in range(len(indices)):
                i, j, k = int(indices[idx, 0]), int(indices[idx, 1]), int(indices[idx, 2])
                # Canonical form: k is smallest, then sort i and j
                triple = (k, min(i, j), max(i, j))
                if triple not in seen_c21:
                    seen_c21.add(triple)
                    c21_values[triple] = norms[idx]

            # Check if we've seen all triples
            if len(seen_c21) == total_c21:
                C21_exact = sum(c21_values.values())
                print(f"  C21 EXACT: All {total_c21} triples sampled in {time.time() - start_time:.3f}s")
                print(f"  C21 exact value: {C21_exact:.6f}")
                C21_est = C21_exact
                break

        if time.time() - start_time >= time_limit / 3:
            break

    if not use_exact or len(seen_c21) < total_c21:
        total_C21 = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
        C21_est = C21_sum * (total_C21 / samples_C21) if samples_C21 > 0 else 0.0
        print(f"  C21 estimation: {samples_C21} samples in {time.time() - start_time:.3f}s")

        # Report standard error if convergence monitoring enabled
        if ENABLE_CONVERGENCE_MONITORING and samples_C21 > 0:
            sample_mean = C21_sum / samples_C21
            sample_var = (C21_sum_sq / samples_C21) - sample_mean**2
            if sample_var > 0:
                sample_std = np.sqrt(sample_var)
                se = sample_std * np.sqrt(total_C21 / samples_C21) / np.sqrt(samples_C21)
                rel_se = se / C21_est if C21_est > 0 else float('inf')
                print(f"  C21 estimate: {C21_est:.6f} ± {se:.6f} (rel. SE: {rel_se:.4f})")
            else:
                print(f"  C21 estimate: {C21_est:.6f}")
        else:
            print(f"  C21 estimate: {C21_est:.6f}")

        if use_exact:
            print(f"    (Sampled {len(seen_c21)}/{total_c21} unique triples, {100*len(seen_c21)/total_c21:.1f}% coverage)")

    # ---------------------------
    # Estimate C22 = sum_{k<j} ||[H_k, [H_k, H_j]]||
    # ---------------------------
    print(f"Estimating C22 with batch_size={batch_size}...")
    C22_sum = 0.0
    C22_sum_sq = 0.0  # For variance calculation
    samples_C22 = 0
    start_time = time.time()

    while time.time() - start_time < time_limit / 3:
        # Generate random pairs (k, j) with k < j
        k_vals = np.random.randint(0, max(1, N - 1), size=batch_size)
        # Vectorized: for each k, sample j from [k+1, N)
        j_vals = k_vals + 1 + np.random.randint(0, N - 1 - k_vals, dtype=np.int64)

        indices = np.column_stack([k_vals, j_vals])

        # Compute nested norms in parallel
        norms = batch_compute_C22(x_bits, z_bits, coeffs, indices, N)

        C22_sum += np.sum(norms)
        C22_sum_sq += np.sum(norms**2)
        samples_C22 += batch_size

        # Check convergence if enabled
        if ENABLE_CONVERGENCE_MONITORING and samples_C22 >= MIN_SAMPLES_FOR_CONVERGENCE:
            total_C22_count = N * (N - 1) / 2
            C22_est_current = C22_sum * (total_C22_count / samples_C22)
            sample_mean = C22_sum / samples_C22
            sample_var = (C22_sum_sq / samples_C22) - sample_mean**2
            if sample_var > 0:
                sample_std = np.sqrt(sample_var)
                se = sample_std * np.sqrt(total_C22_count / samples_C22) / np.sqrt(samples_C22)
                rel_se = se / C22_est_current if C22_est_current > 0 else float('inf')

                if CONVERGENCE_THRESHOLD > 0 and rel_se < CONVERGENCE_THRESHOLD:
                    print(f"  C22 CONVERGED: SE/mean = {rel_se:.4f} < {CONVERGENCE_THRESHOLD:.4f}")
                    print(f"  C22 estimate: {C22_est_current:.6f} ± {se:.6f} (95% CI)")
                    C22_est = C22_est_current
                    break

        # Track exact values if enabled
        if use_exact:
            for idx in range(len(indices)):
                k, j = int(indices[idx, 0]), int(indices[idx, 1])
                pair = (min(k, j), max(k, j))
                if pair not in seen_c22:
                    seen_c22.add(pair)
                    c22_values[pair] = norms[idx]

            # Check if we've seen all pairs
            if len(seen_c22) == total_c22:
                C22_exact = sum(c22_values.values())
                print(f"  C22 EXACT: All {total_c22} pairs sampled in {time.time() - start_time:.3f}s")
                print(f"  C22 exact value: {C22_exact:.6f}")
                C22_est = C22_exact
                break

        if time.time() - start_time >= time_limit / 3:
            break

    if not use_exact or len(seen_c22) < total_c22:
        total_C22 = N * (N - 1) / 2
        C22_est = C22_sum * (total_C22 / samples_C22) if samples_C22 > 0 else 0.0
        print(f"  C22 estimation: {samples_C22} samples in {time.time() - start_time:.3f}s")

        # Report standard error if convergence monitoring enabled
        if ENABLE_CONVERGENCE_MONITORING and samples_C22 > 0:
            sample_mean = C22_sum / samples_C22
            sample_var = (C22_sum_sq / samples_C22) - sample_mean**2
            if sample_var > 0:
                sample_std = np.sqrt(sample_var)
                se = sample_std * np.sqrt(total_C22 / samples_C22) / np.sqrt(samples_C22)
                rel_se = se / C22_est if C22_est > 0 else float('inf')
                print(f"  C22 estimate: {C22_est:.6f} ± {se:.6f} (rel. SE: {rel_se:.4f})")
            else:
                print(f"  C22 estimate: {C22_est:.6f}")
        else:
            print(f"  C22 estimate: {C22_est:.6f}")

        if use_exact:
            print(f"    (Sampled {len(seen_c22)}/{total_c22} unique pairs, {100*len(seen_c22)/total_c22:.1f}% coverage)")

    # ---------------------------
    # Final output
    # ---------------------------
    # Check if we achieved exact computation
    if use_exact and len(seen_c1) == total_c1 and len(seen_c21) == total_c21 and len(seen_c22) == total_c22:
        print("\n" + "="*70)
        print("✅ EXACT COMPUTATION ACHIEVED")
        print(f"   All {total_c1} C1 pairs sampled")
        print(f"   All {total_c21} C21 triples sampled")
        print(f"   All {total_c22} C22 pairs sampled")
        print("="*70)

    # Return C1 and C2 as defined in Childs et al. (arXiv:1912.08854v3)
    # C1 is divided by 2 as per the convention (see VERIFICATION_OF_USER_FIX.md)
    # See Equations 145 and 152 for first- and second-order formulas respectively
    return C1_est / 2, C21_est / 12 + C22_est / 24


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
