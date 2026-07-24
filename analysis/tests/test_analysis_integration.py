"""
Integration tests for the analysis framework with numerical simulation.

Tests the full workflow including:
- Matrix caching when multiple analyses need it
- Integration with analyze_algorithm
- Configuration validation
- End-to-end numerical simulation workflow
"""

import numpy as np
import pytest
import tempfile
import os

from qhat.analysis.config_types import (
    AnalysisConfiguration,
    GeneralConfigurationUser,
    GeneralConfiguration,
)
from qhat.analysis.analysis import analyze_algorithm


# =============================================================================
# Mock Classes for Testing
# =============================================================================

class MockAlgorithm:
    """Mock algorithm with tensor_contract method."""

    def __init__(self, matrix, compute_count_tracker=None):
        """
        Initialize with a unitary matrix.

        Args:
            matrix: The unitary matrix to return
            compute_count_tracker: Optional list to track compute calls
        """
        self._matrix = matrix
        self._compute_count = compute_count_tracker

    def tensor_contract(self):
        """Return the matrix and increment counter if tracker provided."""
        if self._compute_count is not None:
            self._compute_count[0] += 1
        return self._matrix


# =============================================================================
# Configuration Validation Tests
# =============================================================================

def test_analyze_algorithm_requires_at_least_one_analysis():
    """Test that analyze_algorithm requires at least one analysis."""
    config_analysis = AnalysisConfiguration()
    # Don't set any analysis options

    unitary = np.eye(4, dtype=complex)
    algorithm = MockAlgorithm(unitary)

    with pytest.raises(ValueError, match="No analyses requested"):
        analyze_algorithm(config_analysis, algorithm)

    print("✓ Test passed: analyze_algorithm_requires_at_least_one_analysis")


def test_analyze_algorithm_accepts_numerical_simulation():
    """Test that analyze_algorithm accepts numerical_simulation_inputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_analysis = AnalysisConfiguration()

        # Create input state
        input_path = os.path.join(tmpdir, "state.npy")
        state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(input_path, state)

        config_analysis.numerical_simulation_inputs = input_path

        unitary = np.eye(4, dtype=complex)
        algorithm = MockAlgorithm(unitary)

        # Should not raise
        results = analyze_algorithm(config_analysis, algorithm)

        assert 'numerical_simulation' in results
        print("✓ Test passed: analyze_algorithm_accepts_numerical_simulation")


# =============================================================================
# Matrix Caching Tests
# =============================================================================

def test_matrix_computed_once_for_single_analysis():
    """Test that matrix is computed exactly once when only one analysis needs it."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_analysis = AnalysisConfiguration()

        # Create input state
        input_path = os.path.join(tmpdir, "state.npy")
        state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(input_path, state)

        config_analysis.numerical_simulation_inputs = input_path

        # Track computation count
        compute_count = [0]
        unitary = np.eye(4, dtype=complex)
        algorithm = MockAlgorithm(unitary, compute_count)

        analyze_algorithm(config_analysis, algorithm)

        # Should compute exactly once
        assert compute_count[0] == 1, f"Expected 1 computation, got {compute_count[0]}"
        print("✓ Test passed: matrix_computed_once_for_single_analysis")


def test_matrix_computed_once_for_multiple_analyses():
    """Test that matrix is computed once and cached when multiple analyses need it."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_analysis = AnalysisConfiguration()

        # Configure both matrix output and numerical simulation
        matrix_output_path = os.path.join(tmpdir, "matrix.npz")
        config_analysis.algorithm_matrix_output_file = matrix_output_path

        input_path = os.path.join(tmpdir, "state.npy")
        state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(input_path, state)
        config_analysis.numerical_simulation_inputs = input_path

        # Track computation count
        compute_count = [0]
        unitary = np.eye(4, dtype=complex)
        algorithm = MockAlgorithm(unitary, compute_count)

        results = analyze_algorithm(config_analysis, algorithm)

        # Should compute exactly once even though two analyses use it
        assert compute_count[0] == 1, f"Expected 1 computation, got {compute_count[0]}"

        # Verify both analyses ran
        assert 'matrix_output' in results
        assert 'numerical_simulation' in results

        # Verify both outputs exist
        assert os.path.exists(matrix_output_path)
        assert os.path.exists(os.path.join(tmpdir, "state_final.npy"))

        print("✓ Test passed: matrix_computed_once_for_multiple_analyses")


def test_matrix_not_computed_if_not_needed():
    """Test that matrix is not computed if no analysis needs it."""
    config_analysis = AnalysisConfiguration()

    # Only set resource_estimator (doesn't need matrix)
    config_analysis.resource_estimator = "pyliqtr"

    # Track computation count
    compute_count = [0]
    unitary = np.eye(4, dtype=complex)
    algorithm = MockAlgorithm(unitary, compute_count)

    try:
        # This will fail because we don't have a proper algorithm for resource estimation,
        # but we can check that tensor_contract was never called
        analyze_algorithm(config_analysis, algorithm)
    except Exception:
        pass  # Expected to fail, we just want to check compute count

    # Should not compute at all
    assert compute_count[0] == 0, f"Expected 0 computations, got {compute_count[0]}"
    print("✓ Test passed: matrix_not_computed_if_not_needed")


# =============================================================================
# End-to-End Integration Tests
# =============================================================================

def test_end_to_end_single_state():
    """Test complete workflow with a single input state."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_analysis = AnalysisConfiguration()

        # Create a known unitary (Pauli X on first qubit)
        unitary = np.array([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0]
        ], dtype=complex)

        # Create input state |00⟩
        input_path = os.path.join(tmpdir, "input.npy")
        input_state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(input_path, input_state)

        config_analysis.numerical_simulation_inputs = input_path

        algorithm = MockAlgorithm(unitary)
        results = analyze_algorithm(config_analysis, algorithm)

        # Verify results structure
        assert 'numerical_simulation' in results
        assert len(results['numerical_simulation']['simulations']) == 1

        sim_result = results['numerical_simulation']['simulations'][0]
        assert sim_result['input_file'] == input_path
        assert sim_result['output_file'] == os.path.join(tmpdir, "input_final.npy")

        # Verify output state is correct (|10⟩ after X gate)
        output_state = np.load(sim_result['output_file'])
        expected = np.array([0.0, 1.0, 0.0, 0.0], dtype=complex)
        assert np.allclose(output_state, expected)

        print("✓ Test passed: end_to_end_single_state")


def test_end_to_end_multiple_states():
    """Test complete workflow with multiple input states."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_analysis = AnalysisConfiguration()

        # Create Hadamard gate on 2 qubits
        H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
        unitary = np.kron(H, H)

        # Create multiple input states
        input_paths = []
        for i in range(3):
            path = os.path.join(tmpdir, f"state_{i}.npy")
            state = np.zeros(4, dtype=complex)
            state[i] = 1.0
            np.save(path, state)
            input_paths.append(path)

        config_analysis.numerical_simulation_inputs = input_paths

        algorithm = MockAlgorithm(unitary)
        results = analyze_algorithm(config_analysis, algorithm)

        # Verify all simulations completed
        assert len(results['numerical_simulation']['simulations']) == 3

        # Verify all output files exist and have correct norm
        for i, sim_result in enumerate(results['numerical_simulation']['simulations']):
            output_path = os.path.join(tmpdir, f"state_{i}_final.npy")
            assert os.path.exists(output_path)

            output_state = np.load(output_path)
            assert np.abs(np.linalg.norm(output_state) - 1.0) < 1e-10

        print("✓ Test passed: end_to_end_multiple_states")


def test_combined_matrix_output_and_simulation():
    """Test that matrix output and numerical simulation work together."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_analysis = AnalysisConfiguration()

        # Create unitary
        unitary = np.eye(4, dtype=complex) * np.exp(1j * np.pi / 4)  # Global phase

        # Configure both analyses
        matrix_path = os.path.join(tmpdir, "matrix.npz")
        config_analysis.algorithm_matrix_output_file = matrix_path

        input_path = os.path.join(tmpdir, "state.npy")
        state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(input_path, state)
        config_analysis.numerical_simulation_inputs = input_path

        # Track computation count
        compute_count = [0]
        algorithm = MockAlgorithm(unitary, compute_count)

        results = analyze_algorithm(config_analysis, algorithm)

        # Verify both outputs
        assert 'matrix_output' in results
        assert 'numerical_simulation' in results

        # Matrix computed exactly once
        assert compute_count[0] == 1

        # Matrix file created
        assert os.path.exists(matrix_path)
        matrix_data = np.load(matrix_path)
        assert np.allclose(matrix_data['matrix'], unitary)

        # Simulation output created
        output_path = os.path.join(tmpdir, "state_final.npy")
        assert os.path.exists(output_path)

        # Verify simulation used the same matrix
        output_state = np.load(output_path)
        expected = unitary @ state
        assert np.allclose(output_state, expected)

        print("✓ Test passed: combined_matrix_output_and_simulation")


def test_large_system_performance():
    """Test that simulation works efficiently with larger systems."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config_analysis = AnalysisConfiguration()

        # Create 8-qubit system (dimension 256)
        n_qubits = 8
        dim = 2**n_qubits

        # Use sparse unitary (just permutation) for efficiency
        unitary = np.eye(dim, dtype=complex)
        # Swap first and last basis states
        unitary[[0, -1]] = unitary[[-1, 0]]

        # Create input state
        input_path = os.path.join(tmpdir, "state.npy")
        state = np.zeros(dim, dtype=complex)
        state[0] = 1.0
        np.save(input_path, state)

        config_analysis.numerical_simulation_inputs = input_path

        algorithm = MockAlgorithm(unitary)
        results = analyze_algorithm(config_analysis, algorithm)

        # Verify output
        output_path = os.path.join(tmpdir, "state_final.npy")
        output_state = np.load(output_path)

        # Should have swapped to last basis state
        expected = np.zeros(dim, dtype=complex)
        expected[-1] = 1.0
        assert np.allclose(output_state, expected)

        print("✓ Test passed: large_system_performance")


# =============================================================================
# Run all tests
# =============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("Running Analysis Integration Tests")
    print("="*70 + "\n")

    print("Configuration Validation Tests:")
    test_analyze_algorithm_requires_at_least_one_analysis()
    test_analyze_algorithm_accepts_numerical_simulation()

    print("\nMatrix Caching Tests:")
    test_matrix_computed_once_for_single_analysis()
    test_matrix_computed_once_for_multiple_analyses()
    test_matrix_not_computed_if_not_needed()

    print("\nEnd-to-End Integration Tests:")
    test_end_to_end_single_state()
    test_end_to_end_multiple_states()
    test_combined_matrix_output_and_simulation()
    test_large_system_performance()

    print("\n" + "="*70)
    print("All integration tests passed! ✓")
    print("="*70)
