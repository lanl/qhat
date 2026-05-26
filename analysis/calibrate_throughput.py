#!/usr/bin/env python3.11
"""
Calibration script for measuring throughput on your hardware.

This script runs exact computation benchmarks and reports the measured
throughput values for C1, C21, and C22 computations. Use these values
to update THROUGHPUT_CONFIG in trotter_coefficients_fast.py.

Usage:
    python3.11 calibrate_throughput.py [--runs N] [--N TERMS]

Options:
    --runs N      Number of calibration runs to average (default: 3)
    --N TERMS     Number of Pauli terms to use (default: 50)
"""

import sys
import argparse
import numpy as np

from openfermion import QubitOperator

from qhat.analysis.trotter_coefficients_fast import trotter_error_estimator_fast, THROUGHPUT_CONFIG


def generate_random_hamiltonian(N, n_qubits=20, seed=None):
    """Generate random Hamiltonian for calibration."""
    if seed is not None:
        np.random.seed(seed)

    terms = []
    for i in range(N):
        qubit = i % n_qubits
        pauli = np.random.choice(['X', 'Y', 'Z'])
        coeff = np.random.uniform(0.5, 2.0)
        terms.append(QubitOperator(f'{pauli}{qubit}', coeff))
    return terms


def parse_throughput_from_output(output_lines):
    """
    Parse throughput values from computation output.

    Looks for lines like:
      "C1 EXACT: 1225 pairs computed in 0.082s (15.0M pairs/sec)"
    """
    throughputs = {}

    for line in output_lines:
        if 'C1 EXACT:' in line and 'pairs/sec' in line:
            # Extract throughput: look for pattern like "(15.0M pairs/sec)"
            try:
                parts = line.split('(')[-1].split('M')[0]
                throughputs['c1'] = float(parts) * 1e6
            except (ValueError, IndexError):
                pass

        elif 'C21 EXACT:' in line and 'triples/sec' in line:
            try:
                parts = line.split('(')[-1].split('M')[0]
                throughputs['c21'] = float(parts) * 1e6
            except (ValueError, IndexError):
                pass

        elif 'C22 EXACT:' in line and 'pairs/sec' in line:
            try:
                parts = line.split('(')[-1].split('M')[0]
                throughputs['c22'] = float(parts) * 1e6
            except (ValueError, IndexError):
                pass

    return throughputs


def run_calibration(N=50, num_runs=3):
    """
    Run calibration benchmarks.

    Args:
        N: Number of Pauli terms (default 50 gives good measurement)
        num_runs: Number of runs to average

    Returns:
        dict: Average throughput values
    """
    print("="*70)
    print("THROUGHPUT CALIBRATION")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Number of Pauli terms: {N}")
    print(f"  Number of runs: {num_runs}")
    print(f"  Current defaults in THROUGHPUT_CONFIG:")
    print(f"    c1_samples_per_sec:  {THROUGHPUT_CONFIG['c1_samples_per_sec']/1e6:.1f}M/s")
    print(f"    c21_samples_per_sec: {THROUGHPUT_CONFIG['c21_samples_per_sec']/1e6:.1f}M/s")
    print(f"    c22_samples_per_sec: {THROUGHPUT_CONFIG['c22_samples_per_sec']/1e6:.1f}M/s")
    print()

    all_throughputs = {'c1': [], 'c21': [], 'c22': []}

    for run in range(num_runs):
        print("-"*70)
        print(f"Run {run + 1}/{num_runs}")
        print("-"*70)

        # Generate fresh random Hamiltonian for each run (different seed)
        terms = generate_random_hamiltonian(N, seed=42 + run)

        # Capture output by redirecting stdout
        import io
        import contextlib

        output_buffer = io.StringIO()
        with contextlib.redirect_stdout(output_buffer):
            c1, c2 = trotter_error_estimator_fast(terms, 60, mode='exact')

        output = output_buffer.getvalue()
        print(output)  # Print to user's console

        # Parse throughput values
        throughputs = parse_throughput_from_output(output.split('\n'))

        for key in ['c1', 'c21', 'c22']:
            if key in throughputs:
                all_throughputs[key].append(throughputs[key])

        print(f"Run {run + 1} results: C1={c1:.3f}, C2={c2:.6f}")
        if throughputs:
            print(f"Measured throughputs:")
            for key, val in throughputs.items():
                print(f"  {key}: {val/1e6:.1f}M/s")

    # Compute averages
    print()
    print("="*70)
    print("CALIBRATION RESULTS")
    print("="*70)

    averages = {}
    for key in ['c1', 'c21', 'c22']:
        if all_throughputs[key]:
            avg = np.mean(all_throughputs[key])
            std = np.std(all_throughputs[key])
            averages[key] = avg
            print(f"\n{key.upper()} throughput:")
            print(f"  Average: {avg/1e6:.1f}M samples/sec")
            print(f"  Std dev: {std/1e6:.1f}M samples/sec")
            print(f"  Values:  {[f'{v/1e6:.1f}M' for v in all_throughputs[key]]}")

    # Generate update code
    print()
    print("-"*70)
    print("TO UPDATE THROUGHPUT_CONFIG:")
    print("-"*70)
    print("\nOption 1: Modify trotter_coefficients_fast.py directly")
    print("Replace THROUGHPUT_CONFIG with:\n")
    print("THROUGHPUT_CONFIG = {")
    if 'c1' in averages:
        print(f"    'c1_samples_per_sec': {averages['c1']:.0f},   # Your hardware: {averages['c1']/1e6:.1f}M/s")
    if 'c21' in averages:
        print(f"    'c21_samples_per_sec': {averages['c21']:.0f},  # Your hardware: {averages['c21']/1e6:.1f}M/s")
    if 'c22' in averages:
        print(f"    'c22_samples_per_sec': {averages['c22']:.0f},  # Your hardware: {averages['c22']/1e6:.1f}M/s")
    print("}")

    print("\nOption 2: Override at runtime")
    print("Add to your code:\n")
    print("from trotter_coefficients_fast import THROUGHPUT_CONFIG")
    print("my_config = THROUGHPUT_CONFIG.copy()")
    if 'c1' in averages:
        print(f"my_config['c1_samples_per_sec'] = {averages['c1']:.0f}")
    if 'c21' in averages:
        print(f"my_config['c21_samples_per_sec'] = {averages['c21']:.0f}")
    if 'c22' in averages:
        print(f"my_config['c22_samples_per_sec'] = {averages['c22']:.0f}")
    print("# Then pass my_config to should_use_exact_tracking()")

    print("\n" + "="*70)
    print("✅ CALIBRATION COMPLETE")
    print("="*70)

    return averages


def main():
    parser = argparse.ArgumentParser(
        description='Calibrate throughput for your hardware',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--runs', type=int, default=3,
                        help='Number of calibration runs to average (default: 3)')
    parser.add_argument('--N', type=int, default=50,
                        help='Number of Pauli terms to use (default: 50)')

    args = parser.parse_args()

    if args.N < 10:
        print("Warning: N < 10 may not give accurate throughput measurements")
        print("Recommended: N >= 30 for reliable results")
        print()

    if args.N > 100:
        print(f"Warning: N={args.N} may take a long time")
        print("Recommended: N=30-70 for calibration")
        print()

    run_calibration(N=args.N, num_runs=args.runs)


if __name__ == "__main__":
    main()
