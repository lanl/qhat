"""
Tests for exact numerical simulation functionality.

Tests cover:
- Basic exact simulation with known operators
- Multiple input states
- Output file naming (_exact_final suffix)
- Comparison with approximate simulation
- Matrix-free operators
- Norm preservation
"""

import numpy as np
import pytest
import tempfile
import os
from pathlib import Path

from qhat.analysis.analysis import exact_numerical_simulation
from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.file_io import save_state, load_state


# =================================================================================================
# Unit Tests: Basic Exact Simulation
# =================================================================================================

def test_exact_simulation_identity():
    """Test exact simulation with identity operator."""
    # 2x2 identity
    identity = np.eye(2, dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    # Mock hamiltonian (not actually used in simulation)
    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            # Save initial state
            save_state('test_state.npy', state)

            # Run simulation
            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=identity
            )

            # Verify results
            assert 'exact_simulations' in results
            assert len(results['exact_simulations']) == 1
            assert results['exact_simulations'][0]['input_file'] == 'test_state.npy'
            assert results['exact_simulations'][0]['output_file'] == 'test_state_exact_final.npy'

            # Load and verify output
            final_state = load_state('test_state_exact_final.npy')
            np.testing.assert_array_almost_equal(final_state, state)

        finally:
            os.chdir(original_dir)


def test_exact_simulation_pauli_z():
    """Test exact simulation with Pauli Z operator."""
    # Pauli Z = [[1, 0], [0, -1]]
    pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=pauli_z
            )

            # Load output: Z|0⟩ = |0⟩
            final_state = load_state('test_state_exact_final.npy')
            expected = np.array([1.0, 0.0], dtype=complex)
            np.testing.assert_array_almost_equal(final_state, expected)

        finally:
            os.chdir(original_dir)


def test_exact_simulation_pauli_x():
    """Test exact simulation with Pauli X operator."""
    # Pauli X = [[0, 1], [1, 0]]
    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=pauli_x
            )

            # Load output: X|0⟩ = |1⟩
            final_state = load_state('test_state_exact_final.npy')
            expected = np.array([0.0, 1.0], dtype=complex)
            np.testing.assert_array_almost_equal(final_state, expected)

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: Multiple States
# =================================================================================================

def test_exact_simulation_multiple_states():
    """Test exact simulation with multiple input states."""
    identity = np.eye(2, dtype=complex)
    state1 = np.array([1.0, 0.0], dtype=complex)
    state2 = np.array([0.0, 1.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = ['state1.npy', 'state2.npy']

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('state1.npy', state1)
            save_state('state2.npy', state2)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=identity
            )

            assert len(results['exact_simulations']) == 2
            assert results['exact_simulations'][0]['input_file'] == 'state1.npy'
            assert results['exact_simulations'][0]['output_file'] == 'state1_exact_final.npy'
            assert results['exact_simulations'][1]['input_file'] == 'state2.npy'
            assert results['exact_simulations'][1]['output_file'] == 'state2_exact_final.npy'

            # Verify both outputs exist
            assert os.path.exists('state1_exact_final.npy')
            assert os.path.exists('state2_exact_final.npy')

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: Output File Naming
# =================================================================================================

def test_exact_simulation_output_naming():
    """Test that output files have _exact_final suffix."""
    identity = np.eye(2, dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'my_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('my_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=identity
            )

            # Verify naming: my_state.npy -> my_state_exact_final.npy
            assert results['exact_simulations'][0]['output_file'] == 'my_state_exact_final.npy'
            assert os.path.exists('my_state_exact_final.npy')
            assert not os.path.exists('my_state_final.npy')  # Should not use regular suffix

        finally:
            os.chdir(original_dir)


def test_exact_simulation_preserves_path():
    """Test that output files preserve directory path."""
    identity = np.eye(2, dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'subdir/my_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            os.makedirs('subdir', exist_ok=True)
            save_state('subdir/my_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=identity
            )

            # Verify path preserved
            assert results['exact_simulations'][0]['output_file'] == 'subdir/my_state_exact_final.npy'
            assert os.path.exists('subdir/my_state_exact_final.npy')

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: Norm Preservation
# =================================================================================================

def test_exact_simulation_norm_preservation():
    """Test that unitary evolution preserves state norm."""
    # Random unitary (2x2)
    # Use Pauli X as a simple unitary
    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    state = np.array([0.6, 0.8], dtype=complex)  # Normalized

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=pauli_x
            )

            # Verify norm preserved
            assert abs(results['exact_simulations'][0]['output_norm'] - 1.0) < 1e-10

            # Load and verify
            final_state = load_state('test_state_exact_final.npy')
            assert abs(np.linalg.norm(final_state) - 1.0) < 1e-10

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Unit Tests: Matrix-Free Operators
# =================================================================================================

def test_exact_simulation_matrix_free():
    """Test exact simulation with matrix-free operator."""
    from qhat.analysis.matrix_operations import PauliStringOperator

    # Create matrix-free identity operator
    pauli_dict = {'II': 1.0}  # 2-qubit identity
    matrix_free_op = PauliStringOperator(pauli_dict, num_qubits=2)

    state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=matrix_free_op
            )

            # Load output
            final_state = load_state('test_state_exact_final.npy')
            np.testing.assert_array_almost_equal(final_state, state)

        finally:
            os.chdir(original_dir)


def test_exact_simulation_matrix_free_pauli_x():
    """Test exact simulation with matrix-free Pauli X."""
    from qhat.analysis.matrix_operations import PauliStringOperator

    # Create matrix-free Pauli X operator
    pauli_dict = {'X': 1.0}
    matrix_free_op = PauliStringOperator(pauli_dict, num_qubits=1)

    state = np.array([1.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=matrix_free_op
            )

            # Load output: X|0⟩ = |1⟩
            final_state = load_state('test_state_exact_final.npy')
            expected = np.array([0.0, 1.0], dtype=complex)
            np.testing.assert_array_almost_equal(final_state, expected)

        finally:
            os.chdir(original_dir)


# =================================================================================================
# Error Handling Tests
# =================================================================================================

def test_exact_simulation_missing_state_file():
    """Test error when state file doesn't exist."""
    identity = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'nonexistent.npy'

    class MockHamiltonian:
        pass

    with pytest.raises(Exception):  # Will be FileNotFoundError or similar
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                exact_numerical_simulation(
                    config,
                    hamiltonian=MockHamiltonian(),
                    exact_matrix=identity
                )
            finally:
                os.chdir(original_dir)


def test_exact_simulation_dimension_mismatch():
    """Test error when state dimension doesn't match matrix."""
    matrix_2x2 = np.eye(2, dtype=complex)
    state_4d = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    class MockHamiltonian:
        pass

    with pytest.raises(ValueError, match="Dimension mismatch"):
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                save_state('test_state.npy', state_4d)

                exact_numerical_simulation(
                    config,
                    hamiltonian=MockHamiltonian(),
                    exact_matrix=matrix_2x2
                )
            finally:
                os.chdir(original_dir)


def test_exact_simulation_invalid_input_type():
    """Test error when exact_simulation_inputs is not string or list."""
    identity = np.eye(2, dtype=complex)

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 123  # Invalid type

    class MockHamiltonian:
        pass

    with pytest.raises(ValueError, match="must be a string or list of strings"):
        exact_numerical_simulation(
            config,
            hamiltonian=MockHamiltonian(),
            exact_matrix=identity
        )


# =================================================================================================
# Integration Test: Comparison with Known Result
# =================================================================================================

def test_exact_simulation_known_result():
    """Test exact simulation against analytically known result."""
    # Use Hadamard operator: H = 1/sqrt(2) * [[1, 1], [1, -1]]
    hadamard = (1.0 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)  # |0⟩

    config = AnalysisConfiguration()
    config.exact_simulation_inputs = 'test_state.npy'

    class MockHamiltonian:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        original_dir = os.getcwd()
        os.chdir(tmpdir)
        try:
            save_state('test_state.npy', state)

            results = exact_numerical_simulation(
                config,
                hamiltonian=MockHamiltonian(),
                exact_matrix=hadamard
            )

            # Load output: H|0⟩ = 1/sqrt(2) (|0⟩ + |1⟩)
            final_state = load_state('test_state_exact_final.npy')
            expected = (1.0 / np.sqrt(2)) * np.array([1.0, 1.0], dtype=complex)
            np.testing.assert_array_almost_equal(final_state, expected)

        finally:
            os.chdir(original_dir)
