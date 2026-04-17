"""
Example script showing how to use both original and optimized implementations
for resource estimation.

This demonstrates drop-in replacement of the slow trotter_error_estimator
with the fast version.
"""

import sys

# Example usage - replace these with actual molecule data
def example_usage():
    """
    Example showing how to switch between original and fast implementations.
    """

    print("=" * 70)
    print("EXAMPLE: Using both implementations")
    print("=" * 70)

    # Import the original module
    from jkg_utils import generate_resource_estimate

    # Import the fast module
    from jkg_utils_fast import generate_resource_estimate_fast

    print("\nThis example shows how to use both versions:")
    print()
    print("1. ORIGINAL VERSION:")
    print("   from jkg_utils import generate_resource_estimate")
    print("   qubits, tcount = generate_resource_estimate(")
    print("       molecule, n_active_electrons_per_atom,")
    print("       n_active_unocc_orbitals_per_atom, epsilon_per_atom,")
    print("       output_qubits, trotter_order, trotter_error_runtime")
    print("   )")
    print()
    print("2. FAST VERSION (drop-in replacement):")
    print("   from jkg_utils_fast import generate_resource_estimate_fast")
    print("   qubits, tcount = generate_resource_estimate_fast(")
    print("       molecule, n_active_electrons_per_atom,")
    print("       n_active_unocc_orbitals_per_atom, epsilon_per_atom,")
    print("       output_qubits, trotter_order, trotter_error_runtime")
    print("   )")
    print()
    print("3. OR, use just the fast commutator estimator:")
    print("   from openfermion.transforms import jordan_wigner")
    print("   from jkg_utils import build_active_space")
    print("   from jkg_utils_fast import trotter_error_estimator_fast")
    print()
    print("   # Get Hamiltonian and convert to Pauli strings")
    print("   active_hamiltonian, init_state = build_active_space(...)")
    print("   H = list(jordan_wigner(active_hamiltonian))")
    print()
    print("   # Use fast estimator")
    print("   C1, C2 = trotter_error_estimator_fast(H, time_limit=10.0)")
    print()

    print("\n" + "=" * 70)
    print("KEY OPTIMIZATIONS IN THE FAST VERSION:")
    print("=" * 70)
    print()
    print("1. Binary encoding of Pauli strings")
    print("   - Each Pauli operator encoded in 2 bits")
    print("   - Enables fast bitwise operations")
    print()
    print("2. Numba JIT compilation")
    print("   - Near-C performance for critical functions")
    print("   - Automatic parallelization with prange")
    print()
    print("3. Batch processing with NumPy")
    print("   - Generate many random samples at once")
    print("   - Vectorized operations")
    print()
    print("4. Tunable batch_size parameter")
    print("   - Default: 10,000 (vs 100 in original)")
    print("   - Larger batches = better parallelization")
    print("   - Adjust based on problem size")
    print()

    print("\n" + "=" * 70)
    print("EXPECTED SPEEDUPS:")
    print("=" * 70)
    print()
    print("Problem size          | Expected speedup")
    print("-" * 70)
    print("Small (50-100 terms)  | 5-20x")
    print("Medium (200-500)      | 20-50x")
    print("Large (1000+)         | 50-200x")
    print()
    print("Note: Speedup depends on:")
    print("  - Number of CPU cores available")
    print("  - Batch size (tune for your problem)")
    print("  - Sparsity of Pauli strings")
    print()


def show_code_comparison():
    """
    Show side-by-side comparison of key functions.
    """
    print("\n" + "=" * 70)
    print("CODE COMPARISON: Anticommutation check")
    print("=" * 70)
    print()
    print("ORIGINAL (Python with dictionaries):")
    print("-" * 70)
    print("""
def anticommute(key1, key2):
    d1 = pauli_dict(key1)  # Convert to dict
    d2 = pauli_dict(key2)  # Convert to dict
    count = 0
    for q in d1:           # Iterate over qubits
        if q in d2 and d1[q] != d2[q]:
            count += 1
    return (count % 2 == 1)
""")
    print()
    print("FAST VERSION (Numba JIT with bitwise ops):")
    print("-" * 70)
    print("""
@njit  # Compiled to machine code!
def pauli_anticommute(x1, z1, x2, z2):
    # Bitwise XOR and AND operations - single CPU instruction!
    diff = (x1 ^ x2) & (x2 | z2) & (x1 | z1)

    # Fast bit counting
    count = 0
    while diff:
        count += diff & 1
        diff >>= 1

    return (count % 2) == 1
""")
    print()
    print("Why is this faster?")
    print("  - No dictionary operations")
    print("  - Compiled to native machine code")
    print("  - Bitwise operations (single CPU cycles)")
    print("  - Better CPU cache utilization")
    print()


if __name__ == "__main__":
    example_usage()
    show_code_comparison()

    print("\n" + "=" * 70)
    print("To run benchmarks, execute:")
    print("  python benchmark_commutators.py")
    print("=" * 70)
    print()
