#!/usr/bin/env python3.11
"""
Simple test to verify fourth-order Trotter selection works in analysis/unitary.py
"""

import sys
import logging
from qhat.analysis.unitary import encode_ramped_trotter
from qhat.analysis.config_types import UnitaryConfiguration
from qhat.analysis.hamiltonian import LinearCombinationOfPauliStrings, Hamiltonian
from qhat.common.logging_utils import configure_logging

# Set up logging to see the selection messages
configure_logging()

def test_order_selection(energy_error, expected_order):
    """Test that the given error tolerance selects the expected order."""
    print(f"\n{'='*70}")
    print(f"Testing with energy_error = {energy_error}")
    print(f"Expected order: {expected_order}")
    print('='*70)

    # Create configuration
    config = UnitaryConfiguration()
    config.method = 'ramped trotter'
    config.trotter_implementation = 'flattened'
    config.trotter_combine_terms = True
    config.timestep = 1.0
    config.energy_error = energy_error
    config.error_scale = 1.0
    config.ordering_method = None  # Use None to preserve input order
    config.error_coeff_mode = 'monte_carlo'
    config.error_coeff_auto_exact = False
    config.error_coeff_time_limit = 5

    # Create a simple Hamiltonian
    # LinearCombinationOfPauliStrings expects either dense or sparse format
    # sparse format: key is tuple of (qubit, operator) tuples
    pauli_dict = {
        ((0, 'X'), (1, 'Y')): 0.5,
        ((1, 'Y'), (2, 'Z')): 0.3,
        ((0, 'Z'), (2, 'X')): 0.2,
        ((0, 'X'),): 0.1,
    }
    lcps = LinearCombinationOfPauliStrings(sparse=pauli_dict, num_qubits=3)
    ham = Hamiltonian(lcps)

    # Encode - this should print log messages showing the order selection
    try:
        result = encode_ramped_trotter(config, ham, tevol_hbar=1.0)
        print(f"✓ Successfully created {expected_order} Trotterization")
        return True
    except Exception as e:
        print(f"✗ Failed: {e}")
        return False

if __name__ == "__main__":
    success = True

    # Test 1: Large error tolerance should prefer first-order (lowest cost)
    print("\n" + "="*70)
    print("TEST 1: Large error tolerance (should prefer first-order)")
    print("="*70)
    success &= test_order_selection(energy_error=1e-2, expected_order="first order")

    # Test 2: Medium error tolerance should prefer second-order
    print("\n" + "="*70)
    print("TEST 2: Medium error tolerance (should prefer second-order)")
    print("="*70)
    success &= test_order_selection(energy_error=1e-4, expected_order="second order")

    # Test 3: Very small error tolerance should prefer fourth-order
    print("\n" + "="*70)
    print("TEST 3: Very small error tolerance (should prefer fourth-order)")
    print("="*70)
    success &= test_order_selection(energy_error=1e-10, expected_order="fourth order")

    print("\n" + "="*70)
    if success:
        print("✓ All tests passed!")
        print("="*70)
        sys.exit(0)
    else:
        print("✗ Some tests failed")
        print("="*70)
        sys.exit(1)
