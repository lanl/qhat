"""
Fast implementation for Trotter error coefficient estimation.

Based on: Childs et al., "Theory of Trotter Error" (arXiv:1912.08854v3)

Key optimizations:
1. Binary encoding of Pauli strings (X=(1,0), Y=(1,1), Z=(0,1), I=(0,0))
2. NumPy vectorization for batch operations
3. Numba JIT compilation with parallel processing (prange)
4. Bitwise operations for anticommutation checks
5. Exact computation with early termination for small systems

Performance: 100-150x faster than the reference implementation in trotter_coefficients.py

For validation and testing, see trotter_coefficients.py for the reference implementation.
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

# Batch size for exact computation - larger than Monte Carlo for efficiency
# Should provide progress updates every few seconds for large systems
EXACT_COMPUTATION_BATCH_SIZE = 1_000_000  # ~125ms per batch at 8M/sec throughput

# Progress reporting configuration
INITIAL_PROGRESS_INTERVAL = 2.0      # First progress report after 2 seconds
PROGRESS_INTERVAL_GROWTH = 1.5       # Exponential growth factor (1.5x each time)
MAX_PROGRESS_INTERVAL = 300.0        # Cap at 5 minutes between reports

# --------------------------------------------------
# Throughput configuration for feasibility estimation
# --------------------------------------------------
# These values are empirically derived from benchmarking and can be overridden
# if you want to calibrate for specific hardware.
#
# HOW TO RECALIBRATE:
# -------------------
# Method 1: Use the calibration script (RECOMMENDED)
#   cd analysis
#   python3.11 calibrate_throughput.py --runs 3 --N 50
#   # This runs 3 calibration runs, measures throughput, and outputs update code
#
# Method 2: Manual calibration
#   1. Run exact computation on a moderately-sized system (N=50 works well)
#   2. Look at the output for lines like:
#        "C1 EXACT: 1225 pairs computed in 0.082s (15.0M pairs/sec)"
#        "C21 EXACT: 19600 triples computed in 2.450s (8.0M triples/sec)"
#        "C22 EXACT: 1225 pairs computed in 0.111s (11.0M pairs/sec)"
#   3. The throughput values are shown in parentheses (M = millions per second)
#   4. Run multiple times and average to account for system load variations
#
#   Example manual calibration:
#     from openfermion import QubitOperator
#     import numpy as np
#     np.random.seed(42)
#     # Generate N=50 random terms
#     terms = []
#     for i in range(50):
#         qubit = i % 20
#         pauli = np.random.choice(['X', 'Y', 'Z'])
#         coeff = np.random.uniform(0.5, 2.0)
#         terms.append(QubitOperator(f'{pauli}{qubit}', coeff))
#     # Run exact computation and observe throughput in output
#     c1, c2 = trotter_error_estimator_fast(terms, 60, mode='exact')
#
# HOW TO UPDATE THE CONFIGURATION:
# --------------------------------
# Option 1: Modify this file directly (affects all future runs)
#   Change the values in THROUGHPUT_CONFIG below
#
# Option 2: Override at runtime (for specific calculations)
#   from trotter_coefficients_fast import THROUGHPUT_CONFIG
#   my_config = THROUGHPUT_CONFIG.copy()
#   my_config['c1_samples_per_sec'] = 20e6   # Your measured value
#   my_config['c21_samples_per_sec'] = 10e6  # Your measured value
#   my_config['c22_samples_per_sec'] = 14e6  # Your measured value
#   # Use in feasibility check:
#   should_use_exact_tracking(N=100, throughput_config=my_config)
#
# Option 3: Set environment-specific defaults in your application
#   import trotter_coefficients_fast as tcf
#   tcf.THROUGHPUT_CONFIG['c1_samples_per_sec'] = 20e6  # Modifies global default
#
# WHY THESE VALUES MATTER:
# ------------------------
# The throughput values are used to estimate whether exact computation will
# complete within the time budget. If your hardware is faster/slower than
# the defaults, the feasibility checks may be too conservative/optimistic.
# Accurate calibration ensures optimal mode selection (exact vs Monte Carlo).
#
THROUGHPUT_CONFIG = {
    'c1_samples_per_sec': 15e6,   # Pair commutators: ||[H_i, H_j]||
    'c21_samples_per_sec': 8e6,   # Nested triple commutators: ||[H_i, [H_j, H_k]]||
    'c22_samples_per_sec': 11e6,  # Double commutators: ||[H_k, [H_k, H_j]]||
}

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


def should_use_exact_tracking(N, time_budget=None, memory_limit_mb=None, throughput_config=None):
    """
    Decide whether to use exact computation tracking based on time and memory constraints.

    Args:
        N: Number of Pauli terms
        time_budget: Time budget in seconds (defaults to EXACT_COMPUTATION_TIME_BUDGET)
        memory_limit_mb: Memory limit in MB (defaults to EXACT_COMPUTATION_MEMORY_LIMIT_MB)
        throughput_config: Dictionary with throughput rates in samples/sec
            {'c1_samples_per_sec': float, 'c21_samples_per_sec': float, 'c22_samples_per_sec': float}
            If None, uses THROUGHPUT_CONFIG defaults (15M/8M/11M samples/sec)

    Returns:
        bool: True if exact tracking is feasible and beneficial
    """
    if time_budget is None:
        time_budget = EXACT_COMPUTATION_TIME_BUDGET
    if memory_limit_mb is None:
        memory_limit_mb = EXACT_COMPUTATION_MEMORY_LIMIT_MB
    if throughput_config is None:
        throughput_config = THROUGHPUT_CONFIG

    # Compute combination counts
    c1_count = N * (N - 1) // 2
    c21_count = N * (N - 1) * (N - 2) // 6  # C(N, 3) - number of ways to choose 3 items from N
    c22_count = N * (N - 1) // 2

    # Memory check (12 bytes per triple for set overhead, 8 per pair)
    estimated_memory_mb = (c1_count * 8 + c21_count * 12 + c22_count * 8) / (1024**2)
    if estimated_memory_mb > memory_limit_mb:
        return False

    # Time check using coupon collector problem:
    # Expected samples to see all N items = N * H_N ≈ N * (ln(N) + γ)
    # where H_N is the N-th harmonic number and γ ≈ 0.5772 is Euler-Mascheroni constant
    # This is ~11% more accurate than the simple N*ln(N) approximation
    EULER_MASCHERONI = 0.5772156649

    # Use throughput estimates from config
    c1_time = c1_count * (np.log(max(2, c1_count)) + EULER_MASCHERONI) / throughput_config['c1_samples_per_sec'] if c1_count > 0 else 0
    c21_time = c21_count * (np.log(max(2, c21_count)) + EULER_MASCHERONI) / throughput_config['c21_samples_per_sec'] if c21_count > 0 else 0
    c22_time = c22_count * (np.log(max(2, c22_count)) + EULER_MASCHERONI) / throughput_config['c22_samples_per_sec'] if c22_count > 0 else 0
    estimated_time = c1_time + c21_time + c22_time

    if estimated_time > time_budget:
        return False

    return True


def generate_all_pairs(N):
    """
    Generate all unique pairs (i, j) with i < j.

    Args:
        N: Number of terms

    Returns:
        Array of shape (N*(N-1)//2, 2) with all pairs
    """
    count = N * (N - 1) // 2
    indices = np.empty((count, 2), dtype=np.int64)
    idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            indices[idx, 0] = i
            indices[idx, 1] = j
            idx += 1
    return indices


def generate_all_triples(N):
    """
    Generate all unique triples (i, j, k) with k < i and k < j.

    Args:
        N: Number of terms

    Returns:
        Array of shape (total_count, 3) with all triples
    """
    # Count total first - this is C(N, 3)
    total_count = N * (N - 1) * (N - 2) // 6
    indices = np.empty((total_count, 3), dtype=np.int64)

    idx = 0
    for k in range(N - 1):
        for i in range(k + 1, N):
            for j in range(i + 1, N):
                indices[idx, 0] = i
                indices[idx, 1] = j
                indices[idx, 2] = k
                idx += 1
    return indices


def _compute_exact_generic(x_bits, z_bits, coeffs, N, batch_size, config_general,
                           label, index_generator, batch_compute_fn, unit_name):
    """
    Generic function for exact computation of commutator norms.

    This function provides the common implementation for C1, C21, and C22 exact computation,
    reducing code duplication while maintaining performance.

    Args:
        x_bits, z_bits: Arrays of X and Z bits for all Pauli strings
        coeffs: Array of coefficients for all Pauli strings
        N: Number of Pauli strings
        batch_size: Batch size for computation. If None, uses EXACT_COMPUTATION_BATCH_SIZE.
        config_general: GeneralConfiguration object for logging
        label: Label for logging (e.g., "C1", "C21", "C22")
        index_generator: Function to generate all indices, called with (N)
        batch_compute_fn: Function to compute batch norms, called with (x_bits, z_bits, coeffs, batch, N)
        unit_name: Name of units for progress reporting (e.g., "pairs", "triples")

    Returns:
        Exact sum of commutator norms
    """
    if batch_size is None:
        batch_size = EXACT_COMPUTATION_BATCH_SIZE

    config_general.log_verbose(f"  Computing {label} exactly (deterministic enumeration)...")
    start_time = time.time()
    last_progress_time = start_time
    progress_interval = INITIAL_PROGRESS_INTERVAL

    # Generate all indices
    all_indices = index_generator(N)
    total_count = len(all_indices)

    sum_value = 0.0
    # Process in batches for efficiency and progress reporting
    for start_idx in range(0, total_count, batch_size):
        end_idx = min(start_idx + batch_size, total_count)
        batch = all_indices[start_idx:end_idx]
        norms = batch_compute_fn(x_bits, z_bits, coeffs, batch, N)
        sum_value += np.sum(norms)

        # Adaptive progress reporting: frequent at first, then backs off exponentially
        current_time = time.time()
        if current_time - last_progress_time >= progress_interval:
            percent = 100.0 * end_idx / total_count
            elapsed = current_time - start_time
            rate = end_idx / elapsed if elapsed > 0 else 0
            eta = (total_count - end_idx) / rate if rate > 0 else 0
            config_general.log_verbose(f"    Progress: {end_idx:,}/{total_count:,} {unit_name} ({percent:.1f}%) - "
                  f"{rate/1e6:.1f}M {unit_name}/sec - ETA {eta:.1f}s")
            last_progress_time = current_time
            # Increase interval for next report (exponential backoff)
            progress_interval = min(progress_interval * PROGRESS_INTERVAL_GROWTH,
                                   MAX_PROGRESS_INTERVAL)

    elapsed = time.time() - start_time
    rate = total_count / elapsed if elapsed > 0 else 0
    config_general.log_verbose(f"  {label} EXACT: {total_count:,} {unit_name} computed in {elapsed:.3f}s ({rate/1e6:.1f}M {unit_name}/sec)")
    config_general.log_verbose(f"  {label} exact value: {sum_value:.6f}")

    return sum_value


def compute_C1_exact(x_bits, z_bits, coeffs, N, batch_size, config_general):
    """
    Compute C1 exactly by enumerating all pairs deterministically.

    Args:
        batch_size: Batch size for computation. If None, uses EXACT_COMPUTATION_BATCH_SIZE.
                   Larger batches are more efficient but give less frequent progress updates.
        config_general: GeneralConfiguration object for logging

    Returns:
        Exact value of C1 = Σᵢ<ⱼ ||[Hᵢ, Hⱼ]||
    """
    return _compute_exact_generic(
        x_bits, z_bits, coeffs, N, batch_size, config_general,
        label="C1",
        index_generator=generate_all_pairs,
        batch_compute_fn=batch_compute_C1,
        unit_name="pairs"
    )


def compute_C21_exact(x_bits, z_bits, coeffs, N, batch_size, config_general):
    """
    Compute C21 exactly by enumerating all triples deterministically.

    Args:
        batch_size: Batch size for computation. If None, uses EXACT_COMPUTATION_BATCH_SIZE.
        config_general: GeneralConfiguration object for logging

    Returns:
        Exact value of C21 = Σₖ<ᵢ,ₖ<ⱼ ||[Hᵢ, [Hⱼ, Hₖ]]||
    """
    return _compute_exact_generic(
        x_bits, z_bits, coeffs, N, batch_size, config_general,
        label="C21",
        index_generator=generate_all_triples,
        batch_compute_fn=batch_compute_C21,
        unit_name="triples"
    )


def compute_C22_exact(x_bits, z_bits, coeffs, N, batch_size, config_general):
    """
    Compute C22 exactly by enumerating all pairs deterministically.

    Args:
        batch_size: Batch size for computation. If None, uses EXACT_COMPUTATION_BATCH_SIZE.
        config_general: GeneralConfiguration object for logging

    Returns:
        Exact value of C22 = Σₖ<ⱼ ||[Hₖ, [Hₖ, Hⱼ]]||
    """
    return _compute_exact_generic(
        x_bits, z_bits, coeffs, N, batch_size, config_general,
        label="C22",
        index_generator=generate_all_pairs,
        batch_compute_fn=batch_compute_C22,
        unit_name="pairs"
    )


def trotter_error_estimator_fast(pauli_terms, time_limit, config_general,
                                  batch_size=10000, mode='monte_carlo', auto_exact=False):
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
        config_general: GeneralConfiguration object with log_verbose method for logging
        batch_size: Number of samples per batch (larger = better parallelization)
        mode: Computation mode, one of:
            'monte_carlo' (default): Monte Carlo sampling with tracking for early exit
            'exact': Force deterministic exact computation (raises error if infeasible)
        auto_exact: If True and mode='monte_carlo', automatically switch to exact
                   computation when feasible (default False for backward compatibility)

    Returns:
        (C1_est, C2_est): First and second order commutator norms
        where C1_est = C1/2 and C2_est = C21/12 + C22/24

    Raises:
        ValueError: If mode='exact' but exact computation is not feasible

    Notes:
        - Default behavior (mode='monte_carlo', auto_exact=False) is fully backward compatible
        - Set auto_exact=True for automatic optimization on small systems
        - Set mode='exact' to force exact computation (will error if infeasible)
    """
    N = len(pauli_terms)

    # Preprocessing (shared by all paths)
    config_general.log_verbose(f"Preprocessing {N} Pauli terms...")
    start_prep = time.time()
    x_bits, z_bits, coeffs, n_qubits = preprocess_pauli_terms(pauli_terms)
    config_general.log_verbose(f"  Preprocessing done in {time.time() - start_prep:.3f}s ({n_qubits} qubits)")

    # Validate mode parameter
    if mode not in ['monte_carlo', 'exact']:
        raise ValueError(f"Invalid mode '{mode}'. Must be 'monte_carlo' or 'exact'.")

    # Feasibility check (used by multiple paths)
    is_exact_feasible = should_use_exact_tracking(N, time_limit)

    # Decide which computation path to take
    if mode == 'exact':
        # Path B: User explicitly requested exact computation
        if not is_exact_feasible:
            raise ValueError(
                f"Exact computation not feasible for N={N} terms. "
                f"Estimated time or memory exceeds limits. "
                f"Use mode='monte_carlo' instead."
            )
        config_general.log_verbose(f"  Using EXACT computation (user requested, N={N})")
        return _compute_exact_path(x_bits, z_bits, coeffs, N, batch_size, config_general)

    else:  # mode == 'monte_carlo'
        # Path A/C: Monte Carlo, with optional auto-switch
        if auto_exact and is_exact_feasible:
            # Auto-switch to exact computation
            config_general.log_verbose(f"  Auto-switching to EXACT computation (N={N}, feasible within limits)")
            return _compute_exact_path(x_bits, z_bits, coeffs, N, batch_size, config_general)
        else:
            # Use Monte Carlo
            if is_exact_feasible:
                config_general.log_verbose(f"  Using MONTE CARLO with tracking (N={N})")
            else:
                config_general.log_verbose(f"  Using MONTE CARLO sampling (N={N}, too large for exact)")
            use_tracking = is_exact_feasible  # Enable tracking if feasible
            return _compute_monte_carlo_path(
                x_bits, z_bits, coeffs, N, time_limit, batch_size,
                use_tracking, is_exact_feasible, config_general
            )


def _compute_exact_path(x_bits, z_bits, coeffs, N, batch_size, config_general):
    """
    Path B: Deterministic exact computation.

    Args:
        config_general: GeneralConfiguration object for logging
    """
    # Warmup JIT
    config_general.log_verbose(f"Warming up Numba JIT compilation...")
    if N >= 2:
        dummy_indices = np.array([[0, 1]], dtype=np.int64)
        batch_compute_C1(x_bits, z_bits, coeffs, dummy_indices, N)
        if N >= 3:
            dummy_indices_3 = np.array([[0, 1, 2]], dtype=np.int64)
            batch_compute_C21(x_bits, z_bits, coeffs, dummy_indices_3, N)
            batch_compute_C22(x_bits, z_bits, coeffs, dummy_indices[:, [0, 1]], N)
    config_general.log_verbose(f"  Warmup complete")

    # Compute exactly
    C1_exact = compute_C1_exact(x_bits, z_bits, coeffs, N, batch_size, config_general)
    C21_exact = compute_C21_exact(x_bits, z_bits, coeffs, N, batch_size, config_general)
    C22_exact = compute_C22_exact(x_bits, z_bits, coeffs, N, batch_size, config_general)

    # Report
    config_general.log_verbose("\n" + "="*70)
    config_general.log_verbose("✅ EXACT COMPUTATION COMPLETED")
    config_general.log_verbose(f"   All combinations enumerated deterministically")
    config_general.log_verbose("="*70)

    return C1_exact / 2, C21_exact / 12 + C22_exact / 24


def _compute_monte_carlo_path(x_bits, z_bits, coeffs, N, time_limit, batch_size,
                               use_tracking, is_feasible, config_general):
    """
    Path A: Monte Carlo sampling with optional tracking for early exit.

    This is the current implementation, extracted into a helper function.

    Args:
        config_general: GeneralConfiguration object for logging
    """
    # Setup tracking if enabled
    if use_tracking:
        config_general.log_verbose(f"  Exact computation is feasible for N={N} - enabling tracking")
        seen_c1 = set()
        seen_c21 = set()
        seen_c22 = set()
        c1_values = {}
        c21_values = {}
        c22_values = {}

        # Compute total combinations for completion check
        total_c1 = N * (N - 1) // 2
        total_c21 = N * (N - 1) * (N - 2) // 6  # C(N, 3)
        total_c22 = N * (N - 1) // 2

    # Warmup: trigger Numba JIT compilation before timing
    config_general.log_verbose(f"Warming up Numba JIT compilation...")
    if N >= 2:
        dummy_indices = np.array([[0, 1]], dtype=np.int64)
        batch_compute_C1(x_bits, z_bits, coeffs, dummy_indices, N)
        if N >= 3:
            dummy_indices_3 = np.array([[0, 1, 2]], dtype=np.int64)
            batch_compute_C21(x_bits, z_bits, coeffs, dummy_indices_3, N)
            batch_compute_C22(x_bits, z_bits, coeffs, dummy_indices[:, [0, 1]], N)
    config_general.log_verbose(f"  Warmup complete")

    # ---------------------------
    # Estimate C1 = sum_{i<j} ||[H_i, H_j]||
    # ---------------------------
    config_general.log_verbose(f"Estimating C1 with batch_size={batch_size}...")
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

        # Track exact values if enabled
        if use_tracking:
            for idx in range(len(indices)):
                i, j = int(indices[idx, 0]), int(indices[idx, 1])
                pair = (min(i, j), max(i, j))
                if pair not in seen_c1:
                    seen_c1.add(pair)
                    c1_values[pair] = norms[idx]

            # Check if we've seen all pairs
            if len(seen_c1) == total_c1:
                C1_exact = sum(c1_values.values())
                config_general.log_verbose(f"  C1 EXACT: All {total_c1} pairs sampled in {time.time() - start_time:.3f}s")
                config_general.log_verbose(f"  C1 exact value: {C1_exact:.6f}")
                C1_est = C1_exact
                break

        if time.time() - start_time >= time_limit / 3:
            break

    if not use_tracking or len(seen_c1) < total_c1:
        total_C1 = N * (N - 1) / 2
        C1_est = C1_sum * (total_C1 / samples_C1) if samples_C1 > 0 else 0.0
        config_general.log_verbose(f"  C1 estimation: {samples_C1} samples in {time.time() - start_time:.3f}s")
        config_general.log_verbose(f"  C1 estimate: {C1_est:.6f}")

        if use_tracking:
            config_general.log_verbose(f"    (Sampled {len(seen_c1)}/{total_c1} unique pairs, {100*len(seen_c1)/total_c1:.1f}% coverage)")

    # ---------------------------
    # Estimate C21 = sum_{k<j, k<i} ||[H_i, [H_j, H_k]]||
    # ---------------------------
    config_general.log_verbose(f"Estimating C21 with batch_size={batch_size}...")
    C21_sum = 0.0
    samples_C21 = 0
    start_time = time.time()

    # Need at least 3 terms for triples
    if N < 3:
        C21_est = 0.0
        config_general.log_verbose(f"  C21 estimation: N={N} too small for triples, C21=0")
    else:
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
            samples_C21 += len(indices)

            # Track exact values if enabled
            if use_tracking:
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
                    config_general.log_verbose(f"  C21 EXACT: All {total_c21} triples sampled in {time.time() - start_time:.3f}s")
                    config_general.log_verbose(f"  C21 exact value: {C21_exact:.6f}")
                    C21_est = C21_exact
                    break

            if time.time() - start_time >= time_limit / 3:
                break

        if not use_tracking or len(seen_c21) < total_c21:
            total_C21 = N * (N - 1) * (N - 2) // 6  # C(N, 3)
            C21_est = C21_sum * (total_C21 / samples_C21) if samples_C21 > 0 else 0.0
            config_general.log_verbose(f"  C21 estimation: {samples_C21} samples in {time.time() - start_time:.3f}s")
            config_general.log_verbose(f"  C21 estimate: {C21_est:.6f}")

            if use_tracking:
                config_general.log_verbose(f"    (Sampled {len(seen_c21)}/{total_c21} unique triples, {100*len(seen_c21)/total_c21:.1f}% coverage)")

    # ---------------------------
    # Estimate C22 = sum_{k<j} ||[H_k, [H_k, H_j]]||
    # ---------------------------
    config_general.log_verbose(f"Estimating C22 with batch_size={batch_size}...")
    C22_sum = 0.0
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
        samples_C22 += batch_size

        # Track exact values if enabled
        if use_tracking:
            for idx in range(len(indices)):
                k, j = int(indices[idx, 0]), int(indices[idx, 1])
                pair = (min(k, j), max(k, j))
                if pair not in seen_c22:
                    seen_c22.add(pair)
                    c22_values[pair] = norms[idx]

            # Check if we've seen all pairs
            if len(seen_c22) == total_c22:
                C22_exact = sum(c22_values.values())
                config_general.log_verbose(f"  C22 EXACT: All {total_c22} pairs sampled in {time.time() - start_time:.3f}s")
                config_general.log_verbose(f"  C22 exact value: {C22_exact:.6f}")
                C22_est = C22_exact
                break

        if time.time() - start_time >= time_limit / 3:
            break

    if not use_tracking or len(seen_c22) < total_c22:
        total_C22 = N * (N - 1) / 2
        C22_est = C22_sum * (total_C22 / samples_C22) if samples_C22 > 0 else 0.0
        config_general.log_verbose(f"  C22 estimation: {samples_C22} samples in {time.time() - start_time:.3f}s")
        config_general.log_verbose(f"  C22 estimate: {C22_est:.6f}")

        if use_tracking:
            config_general.log_verbose(f"    (Sampled {len(seen_c22)}/{total_c22} unique pairs, {100*len(seen_c22)/total_c22:.1f}% coverage)")

    # ---------------------------
    # Final output
    # ---------------------------
    # Check if we achieved exact computation
    if use_tracking and len(seen_c1) == total_c1 and len(seen_c21) == total_c21 and len(seen_c22) == total_c22:
        config_general.log_verbose("\n" + "="*70)
        config_general.log_verbose("✅ EXACT COMPUTATION ACHIEVED")
        config_general.log_verbose(f"   All {total_c1} C1 pairs sampled")
        config_general.log_verbose(f"   All {total_c21} C21 triples sampled")
        config_general.log_verbose(f"   All {total_c22} C22 pairs sampled")
        config_general.log_verbose("="*70)

    # Return C1 and C2 as defined in Childs et al. (arXiv:1912.08854v3)
    # C1 is divided by 2 as per the convention (see VERIFICATION_OF_USER_FIX.md)
    # See Equations 145 and 152 for first- and second-order formulas respectively
    return C1_est / 2, C21_est / 12 + C22_est / 24

