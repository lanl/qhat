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


def test_original_implementation_is_deprecated(simple_hamiltonian, caplog):
    config = create_config(energy_error=1e-2)
    config.trotter_implementation = 'original'
    with pytest.warns(DeprecationWarning):
        encode_ramped_trotter(config, simple_hamiltonian, tevol_hbar=1.0)
    assert "trotter_implementation='original' is deprecated" in caplog.text


# -------------------------------------------------------------------------------------------------
# "original" with a controlled algorithm (plan D5)
# -------------------------------------------------------------------------------------------------

CONTROLLED_METHODS = ["QPE: QHAT flexible", "QPE: qualtran textbook", "QPE: qualtran qubitization",
                      "QPE: pyLIQTR qubitized", "controlled time evolution",
                      "Controlled Time Evolution"]


def algorithm_config(method):
    from qhat.analysis.config_types import AlgorithmConfiguration
    config = AlgorithmConfiguration()
    config.method = method
    return config


@pytest.mark.parametrize("method", CONTROLLED_METHODS)
def test_original_with_controlled_algorithm_raises(method):
    from qhat.analysis.algorithm import check_unitary_supports_algorithm
    config = create_config(energy_error=1e-2)
    config.trotter_implementation = 'original'
    with pytest.raises(ValueError, match='trotter_implementation = "flattened"'):
        check_unitary_supports_algorithm(config, algorithm_config(method))


@pytest.mark.parametrize("implementation, method", [
    ('original', "time evolution"),
    *[('flattened', method) for method in CONTROLLED_METHODS + ["time evolution"]]])
def test_allowed_combinations_do_not_raise(implementation, method):
    from qhat.analysis.algorithm import check_unitary_supports_algorithm
    config = create_config(energy_error=1e-2)
    config.trotter_implementation = implementation
    check_unitary_supports_algorithm(config, algorithm_config(method))


def test_original_check_ignores_other_unitaries():
    from qhat.analysis.algorithm import check_unitary_supports_algorithm
    config = UnitaryConfiguration()
    config.encode_pauli_lcu(energy_error=1e-3)
    check_unitary_supports_algorithm(config, algorithm_config("QPE: qualtran qubitization"))


def test_driver_rejects_original_with_qpe(monkeypatch, tmp_path):
    """The driver checks before loading the Hamiltonian (the file named here doesn't exist)."""
    import sys
    from qhat.analysis import driver
    config_file = tmp_path / "config.py"
    config_file.write_text(
        f'general.output_directory = "{tmp_path}"\n'
        'hamiltonian.load_pauli_strings("does_not_exist.json")\n'
        'unitary.encode_ramped_trotter(energy_error=1e-2, trotter_implementation="original")\n'
        'algorithm.method = "QPE: QHAT flexible"\n'
        'algorithm.num_phase_qubits = 3\n')
    monkeypatch.setattr(sys, 'argv', ['driver.py', str(config_file)])
    with pytest.raises(ValueError, match='trotter_implementation = "flattened"'):
        driver.run()
