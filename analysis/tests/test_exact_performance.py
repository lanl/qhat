"""
Test performance characteristics of exact vs Monte Carlo computation.
"""

import time
import numpy as np
from openfermion import QubitOperator
import sys
import os

from qhat.analysis.trotter_coefficients_fast import trotter_error_estimator_fast



def generate_random_hamiltonian(N, n_qubits=50):
    """Generate random Hamiltonian for testing."""
    terms = []
    for i in range(N):
        qubit = np.random.randint(0, n_qubits)
        pauli = np.random.choice(['X', 'Y', 'Z'])
        coeff = np.random.uniform(0.5, 2.0)
        terms.append(QubitOperator(f'{pauli}{qubit}', coeff))
    return terms


def test_exact_faster_for_small_N():
    """Exact should be faster than Monte Carlo for small N."""
    print("\n" + "="*70)
    print("Performance Test: Exact vs Monte Carlo for N=50")
    print("="*70)

    np.random.seed(42)
    terms = generate_random_hamiltonian(N=50)

    # Time Monte Carlo (without auto-exact to ensure it doesn't switch)
    print("\nTiming Monte Carlo (60s budget)...")
    start = time.time()
    c1_mc, c2_mc = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=False)
    time_mc = time.time() - start

    # Time exact
    print("\nTiming Exact...")
    start = time.time()
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 60, mode='exact')
    time_ex = time.time() - start

    print("\n" + "-"*70)
    print("Results:")
    print(f"  Monte Carlo: {time_mc:.2f}s → C1={c1_mc:.6f}, C2={c2_mc:.6f}")
    print(f"  Exact:       {time_ex:.2f}s → C1={c1_ex:.6f}, C2={c2_ex:.6f}")
    print(f"  Speedup:     {time_mc/time_ex:.1f}×")
    print("-"*70)

    # Exact should be faster
    if time_ex < time_mc:
        print("✓ Exact is faster than Monte Carlo for N=50")
    else:
        print(f"⚠ Warning: Exact ({time_ex:.2f}s) was slower than MC ({time_mc:.2f}s)")

    # Results should be close (exact is deterministic, MC is stochastic)
    rel_diff_c1 = abs(c1_ex - c1_mc) / c1_ex if c1_ex > 0 else 0
    rel_diff_c2 = abs(c2_ex - c2_mc) / c2_ex if c2_ex > 0 else 0
    print(f"\nRelative differences: C1={100*rel_diff_c1:.1f}%, C2={100*rel_diff_c2:.1f}%")


def test_auto_exact_switches():
    """Auto-exact should switch to exact for small N."""
    print("\n" + "="*70)
    print("Performance Test: Auto-exact switching for N=30")
    print("="*70)

    np.random.seed(43)
    terms = generate_random_hamiltonian(N=30)

    # Time auto-exact (should switch to exact)
    print("\nTiming Auto-exact (should switch)...")
    start = time.time()
    c1_auto, c2_auto = trotter_error_estimator_fast(terms, 60, mode='monte_carlo', auto_exact=True)
    time_auto = time.time() - start

    # Time exact directly
    print("\nTiming Exact directly...")
    start = time.time()
    c1_ex, c2_ex = trotter_error_estimator_fast(terms, 60, mode='exact')
    time_ex = time.time() - start

    print("\n" + "-"*70)
    print("Results:")
    print(f"  Auto-exact: {time_auto:.2f}s → C1={c1_auto:.6f}, C2={c2_auto:.6f}")
    print(f"  Exact:      {time_ex:.2f}s → C1={c1_ex:.6f}, C2={c2_ex:.6f}")
    print(f"  Time ratio: {time_auto/time_ex:.2f}")
    print("-"*70)

    # Auto-exact should give same answer as exact (both deterministic)
    assert abs(c1_auto - c1_ex) < 1e-10, "Auto-exact C1 differs from exact"
    assert abs(c2_auto - c2_ex) < 1e-10, "Auto-exact C2 differs from exact"
    print("✓ Auto-exact gives identical results to exact")

    # Times should be similar (within 2× due to preprocessing overhead)
    if time_auto < 2 * time_ex:
        print("✓ Auto-exact time overhead is acceptable")
    else:
        print(f"⚠ Warning: Auto-exact overhead is large ({time_auto/time_ex:.1f}×)")


def test_scaling_with_N():
    """Test how computation time scales with N."""
    print("\n" + "="*70)
    print("Performance Test: Scaling with N")
    print("="*70)

    test_sizes = [10, 20, 30, 40, 50]
    times = []

    print(f"\n{'N':>5} | {'Pairs':>10} | {'Triples':>12} | {'Time (s)':>10} | {'Rate (M/s)':>11}")
    print("-"*70)

    for N in test_sizes:
        np.random.seed(100 + N)
        terms = generate_random_hamiltonian(N)

        start = time.time()
        c1, c2 = trotter_error_estimator_fast(terms, 60, mode='exact')
        elapsed = time.time() - start

        n_pairs = N * (N - 1) // 2
        n_triples = sum(range(N-1)) if N >= 3 else 0
        total_ops = n_pairs * 2 + n_triples  # C1, C22 use pairs, C21 uses triples
        rate = total_ops / elapsed / 1e6 if elapsed > 0 else 0

        times.append(elapsed)
        print(f"{N:5} | {n_pairs:10,} | {n_triples:12,} | {elapsed:10.3f} | {rate:11.2f}")

    print("-"*70)
    print("✓ Scaling test completed")


if __name__ == "__main__":
    print("="*70)
    print("PERFORMANCE TESTS")
    print("="*70)

    test_exact_faster_for_small_N()
    test_auto_exact_switches()
    test_scaling_with_N()

    print("\n" + "="*70)
    print("✅ ALL PERFORMANCE TESTS COMPLETED!")
    print("="*70)
