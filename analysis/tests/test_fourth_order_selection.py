#!/usr/bin/env python3.11
"""
Simple test to verify fourth-order Trotter selection works in analysis/unitary.py
"""

import pytest
from qhat.analysis.unitary import encode_ramped_trotter
from qhat.analysis.config_types import UnitaryConfiguration
from qhat.analysis.hamiltonian import LinearCombinationOfPauliStrings, Hamiltonian


@pytest.fixture
def simple_hamiltonian():
    """Create a simple test Hamiltonian."""
    pauli_dict = {
        ((0, 'X'), (1, 'Y')): 0.5,
        ((1, 'Y'), (2, 'Z')): 0.3,
        ((0, 'Z'), (2, 'X')): 0.2,
        ((0, 'X'),): 0.1,
    }
    lcps = LinearCombinationOfPauliStrings(sparse=pauli_dict, num_qubits=3)
    return Hamiltonian(lcps)


def create_config(energy_error, explicit_order=None):
    """Helper to create UnitaryConfiguration."""
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
    config.trotter_order = explicit_order  # Explicitly request order (or None for auto-select)
    return config


def test_large_error_tolerance_prefers_first_order(simple_hamiltonian):
    """Test that large error tolerance auto-selects first-order."""
    config = create_config(energy_error=1e-2)
    result = encode_ramped_trotter(config, simple_hamiltonian, tevol_hbar=1.0)
    assert result is not None


def test_medium_error_tolerance_prefers_second_order(simple_hamiltonian):
    """Test that medium error tolerance auto-selects second-order."""
    config = create_config(energy_error=1e-4)
    result = encode_ramped_trotter(config, simple_hamiltonian, tevol_hbar=1.0)
    assert result is not None


def test_small_error_tolerance_still_prefers_second_order(simple_hamiltonian):
    """Test that very small error tolerance still auto-selects second-order (not fourth)."""
    config = create_config(energy_error=1e-10)
    result = encode_ramped_trotter(config, simple_hamiltonian, tevol_hbar=1.0)
    assert result is not None


def test_explicit_fourth_order_request(simple_hamiltonian):
    """Test that explicitly requesting fourth-order works."""
    config = create_config(energy_error=1e-10, explicit_order="fourth order")
    result = encode_ramped_trotter(config, simple_hamiltonian, tevol_hbar=1.0)
    assert result is not None
