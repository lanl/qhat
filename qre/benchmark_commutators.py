"""
Benchmark script to compare original vs. optimized Monte Carlo commutator estimation.

Usage:
    python benchmark_commutators.py
"""

import time
import numpy as np
from openfermion import QubitOperator

# Import both versions
from jkg_utils import trotter_error_estimator
from jkg_utils_fast import trotter_error_estimator_fast


def generate_random_pauli_terms(n_terms, n_qubits, max_weight=4):
    """
    Generate random Pauli terms for testing.

    Args:
        n_terms: Number of Pauli terms to generate
        n_qubits: Number of qubits
        max_weight: Maximum weight (number of non-identity operators) per term

    Returns:
        List of QubitOperator terms
    """
    pauli_ops = ['X', 'Y', 'Z']
    terms = []

    for _ in range(n_terms):
        # Random weight between 1 and max_weight
        weight = np.random.randint(1, max_weight + 1)

        # Random qubits (without replacement)
        qubits = np.random.choice(n_qubits, size=weight, replace=False)

        # Random Pauli operators
        ops = np.random.choice(pauli_ops, size=weight)

        # Build the Pauli string
        pauli_string = ' '.join([f'{op}{q}' for q, op in zip(qubits, ops)])

        # Random coefficient
        coeff = np.random.randn() + 1j * np.random.randn()

        terms.append(QubitOperator(pauli_string, coeff))

    return terms


def run_benchmark(n_terms, n_qubits, time_limit, batch_sizes=None):
    """
    Run benchmark comparing original and fast implementations.

    Args:
        n_terms: Number of Pauli terms
        n_qubits: Number of qubits
        time_limit: Time limit for estimation (seconds)
        batch_sizes: List of batch sizes to test for fast version
    """
    if batch_sizes is None:
        batch_sizes = [100, 1000, 10000]

    print("=" * 70)
    print(f"BENCHMARK: {n_terms} Pauli terms, {n_qubits} qubits, {time_limit}s time limit")
    print("=" * 70)

    # Generate random Pauli terms
    print("\nGenerating random Pauli terms...")
    pauli_terms = generate_random_pauli_terms(n_terms, n_qubits)
    print(f"Generated {len(pauli_terms)} terms")

    # Run original implementation
    print("\n" + "-" * 70)
    print("ORIGINAL IMPLEMENTATION (batch_size=100)")
    print("-" * 70)
    start = time.time()
    C1_orig, C2_orig = trotter_error_estimator(pauli_terms, time_limit, batch_size=100)
    time_orig = time.time() - start

    print(f"\nResults:")
    print(f"  C1 = {C1_orig:.6f}")
    print(f"  C2 = {C2_orig:.6f}")
    print(f"  Total time: {time_orig:.3f}s")

    # Run optimized implementation with different batch sizes
    results_fast = []

    for batch_size in batch_sizes:
        print("\n" + "-" * 70)
        print(f"OPTIMIZED IMPLEMENTATION (batch_size={batch_size})")
        print("-" * 70)
        start = time.time()
        C1_fast, C2_fast = trotter_error_estimator_fast(pauli_terms, time_limit, batch_size=batch_size)
        time_fast = time.time() - start

        print(f"\nResults:")
        print(f"  C1 = {C1_fast:.6f}")
        print(f"  C2 = {C2_fast:.6f}")
        print(f"  Total time: {time_fast:.3f}s")
        print(f"  Speedup: {time_orig / time_fast:.2f}x")
        print(f"  C1 relative difference: {abs(C1_fast - C1_orig) / abs(C1_orig) * 100:.2f}%")
        print(f"  C2 relative difference: {abs(C2_fast - C2_orig) / abs(C2_orig) * 100:.2f}%")

        results_fast.append({
            'batch_size': batch_size,
            'C1': C1_fast,
            'C2': C2_fast,
            'time': time_fast,
            'speedup': time_orig / time_fast
        })

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nOriginal: {time_orig:.3f}s")
    print(f"\nOptimized versions:")
    for result in results_fast:
        print(f"  batch_size={result['batch_size']:5d}: "
              f"{result['time']:6.3f}s (speedup: {result['speedup']:5.2f}x)")

    best_result = max(results_fast, key=lambda r: r['speedup'])
    print(f"\nBest speedup: {best_result['speedup']:.2f}x with batch_size={best_result['batch_size']}")

    return {
        'original': {'C1': C1_orig, 'C2': C2_orig, 'time': time_orig},
        'fast': results_fast
    }


if __name__ == "__main__":
    # Small test
    print("\n\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 20 + "SMALL TEST (50 terms)" + " " * 27 + "║")
    print("╚" + "═" * 68 + "╝")
    run_benchmark(n_terms=50, n_qubits=10, time_limit=3.0, batch_sizes=[100, 1000])

    # Medium test
    print("\n\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 19 + "MEDIUM TEST (200 terms)" + " " * 26 + "║")
    print("╚" + "═" * 68 + "╝")
    run_benchmark(n_terms=200, n_qubits=20, time_limit=6.0, batch_sizes=[1000, 5000, 10000])

    # Large test
    print("\n\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 19 + "LARGE TEST (1000 terms)" + " " * 25 + "║")
    print("╚" + "═" * 68 + "╝")
    run_benchmark(n_terms=1000, n_qubits=50, time_limit=12.0, batch_sizes=[5000, 10000, 20000])

    print("\n\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
