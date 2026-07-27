"""
Tests for numerical simulation analysis.

Includes:
- State vector I/O (loading and saving)
- Dimension validation
- Matrix-vector multiplication correctness
- Multiple input state handling
- Error handling for invalid inputs
- Integration with the analysis framework
"""

import numpy as np
import pytest
import tempfile
import os
from pathlib import Path

from qhat.analysis.analysis import numerical_simulation
from qhat.analysis.file_io import (
    load_state,
    save_state,
    _get_state_format_from_extension,
)
from qhat.analysis.config_types import AnalysisConfiguration


# =============================================================================
# State I/O Tests
# =============================================================================

def test_load_state_numpy_basic():
    """Test loading a basic NumPy state vector."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a simple 2-qubit state
        path = os.path.join(tmpdir, "test_state.npy")
        state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(path, state)

        # Load it back
        loaded = load_state(path)

        assert loaded.shape == (4,)
        assert np.allclose(loaded, state)
        assert loaded.dtype == complex
        print("✓ Test passed: load_state_numpy_basic")


def test_load_state_converts_real_to_complex():
    """Test that real arrays are converted to complex."""
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test_state.npy")
        # Save as real array
        state = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
        np.save(path, state)

        # Load should convert to complex
        loaded = load_state(path)

        assert loaded.dtype == complex
        print("✓ Test passed: load_state_converts_real_to_complex")


def test_save_and_load_roundtrip():
    """Test that save/load roundtrip preserves data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test_state.npy")

        # Create a complex state
        original = np.array([0.6+0.8j, 0.0, 0.0, 0.0], dtype=complex)

        # Save and load
        save_state(path, original)
        loaded = load_state(path)

        assert np.allclose(loaded, original)
        print("✓ Test passed: save_and_load_roundtrip")


def test_get_state_format_from_extension():
    """Test format detection from file extensions."""
    assert _get_state_format_from_extension("state.npy") == "numpy"
    assert _get_state_format_from_extension("/path/to/state.npy") == "numpy"
    assert _get_state_format_from_extension("state.NPY") == "numpy"  # case insensitive
    print("✓ Test passed: get_state_format_from_extension")


def test_invalid_extension_raises_error():
    """Test that invalid extensions raise ValueError."""
    with pytest.raises(ValueError, match="Cannot determine state format"):
        _get_state_format_from_extension("state.txt")
    print("✓ Test passed: invalid_extension_raises_error")


def test_load_nonexistent_file_raises_error():
    """Test that loading nonexistent file raises error."""
    with pytest.raises(FileNotFoundError):
        load_state("/nonexistent/path/state.npy")
    print("✓ Test passed: load_nonexistent_file_raises_error")


def test_load_wrong_dimension_not_power_of_2():
    """Test that loading state with non-power-of-2 dimension raises error."""
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test_state.npy")
        # Create state with dimension 5 (not a power of 2)
        state = np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(path, state)

        with pytest.raises(ValueError, match="not a power of 2"):
            load_state(path)
    print("✓ Test passed: load_wrong_dimension_not_power_of_2")


def test_load_multidimensional_array_raises_error():
    """Test that loading multidimensional array raises error."""
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test_state.npy")
        # Create 2D array
        state = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
        np.save(path, state)

        with pytest.raises(ValueError, match="must be 1-dimensional"):
            load_state(path)
    print("✓ Test passed: load_multidimensional_array_raises_error")


# =============================================================================
# Mock Algorithm for Testing
# =============================================================================

class MockAlgorithm:
    """Mock algorithm object for testing numerical simulation."""

    def __init__(self, matrix):
        """Initialize with a specific unitary matrix."""
        self._matrix = matrix

    def tensor_contract(self):
        """Return the stored matrix."""
        return self._matrix


# =============================================================================
# Numerical Simulation Function Tests
# =============================================================================

def test_numerical_simulation_single_input():
    """Test numerical simulation with a single input state."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a simple 2-qubit identity matrix
        dim = 4
        unitary = np.eye(dim, dtype=complex)

        # Create input state
        input_path = os.path.join(tmpdir, "input.npy")
        input_state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(input_path, input_state)

        # Run simulation with precomputed matrix
        algorithm = MockAlgorithm(unitary)
        result = numerical_simulation(
            input_path,
            algorithm,
            unitary_matrix=unitary
        )

        # Verify result structure
        assert 'simulations' in result
        assert len(result['simulations']) == 1
        assert result['simulations'][0]['input_file'] == input_path

        # Verify output file was created
        output_path = os.path.join(tmpdir, "input_final.npy")
        assert os.path.exists(output_path)

        # Verify output state (should be unchanged with identity)
        output_state = np.load(output_path)
        assert np.allclose(output_state, input_state)
        print("✓ Test passed: numerical_simulation_single_input")


def test_numerical_simulation_multiple_inputs():
    """Test numerical simulation with multiple input states."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a simple 2-qubit unitary (bit flip on first qubit)
        dim = 4
        unitary = np.array([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0]
        ], dtype=complex)

        # Create multiple input states
        input_paths = []
        expected_outputs = []

        for i in range(3):
            path = os.path.join(tmpdir, f"input_{i}.npy")
            state = np.zeros(dim, dtype=complex)
            state[i] = 1.0
            np.save(path, state)
            input_paths.append(path)
            # Compute expected output manually
            expected_outputs.append(unitary @ state)

        # Run simulation
        algorithm = MockAlgorithm(unitary)
        result = numerical_simulation(
            input_paths,
            algorithm,
            unitary_matrix=unitary
        )

        # Verify results
        assert len(result['simulations']) == 3

        for i, sim_result in enumerate(result['simulations']):
            output_path = os.path.join(tmpdir, f"input_{i}_final.npy")
            assert os.path.exists(output_path)

            output_state = np.load(output_path)
            assert np.allclose(output_state, expected_outputs[i])
            assert np.abs(sim_result['output_norm'] - 1.0) < 1e-10

        print("✓ Test passed: numerical_simulation_multiple_inputs")


def test_numerical_simulation_preserves_unitarity():
    """Test that simulation preserves norm (unitarity)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a non-trivial unitary (Hadamard on 2 qubits, applied to each)
        H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
        unitary = np.kron(H, H)  # 4x4 unitary

        # Create a normalized input state
        input_path = os.path.join(tmpdir, "input.npy")
        input_state = np.array([0.5, 0.5, 0.5, 0.5], dtype=complex)
        np.save(input_path, input_state)

        # Run simulation
        algorithm = MockAlgorithm(unitary)
        result = numerical_simulation(
            input_path,
            algorithm,
            unitary_matrix=unitary
        )

        # Check that output norm is preserved
        assert np.abs(result['simulations'][0]['output_norm'] - 1.0) < 1e-10

        # Verify actual output
        output_path = os.path.join(tmpdir, "input_final.npy")
        output_state = np.load(output_path)
        assert np.abs(np.linalg.norm(output_state) - 1.0) < 1e-10

        print("✓ Test passed: numerical_simulation_preserves_unitarity")


def test_numerical_simulation_computes_correct_result():
    """Test that simulation computes correct matrix-vector product."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a known unitary and state
        # Pauli X on first qubit of 2-qubit system
        dim = 4
        unitary = np.array([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0]
        ], dtype=complex)

        # Input: |00⟩
        input_path = os.path.join(tmpdir, "input.npy")
        input_state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        np.save(input_path, input_state)

        # Expected output: |10⟩ (first qubit flipped)
        expected_output = np.array([0.0, 1.0, 0.0, 0.0], dtype=complex)

        # Run simulation
        algorithm = MockAlgorithm(unitary)
        numerical_simulation(
            input_path,
            algorithm,
            unitary_matrix=unitary
        )

        # Verify output
        output_path = os.path.join(tmpdir, "input_final.npy")
        output_state = np.load(output_path)
        assert np.allclose(output_state, expected_output)

        print("✓ Test passed: numerical_simulation_computes_correct_result")


def test_numerical_simulation_dimension_mismatch_raises_error():
    """Test that dimension mismatch between state and matrix raises error."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create 2-qubit unitary
        unitary = np.eye(4, dtype=complex)

        # Create 3-qubit state (dimension mismatch)
        input_path = os.path.join(tmpdir, "input.npy")
        input_state = np.zeros(8, dtype=complex)
        input_state[0] = 1.0
        np.save(input_path, input_state)

        # Should raise dimension mismatch error
        algorithm = MockAlgorithm(unitary)
        with pytest.raises(ValueError, match="Dimension mismatch"):
            numerical_simulation(
                input_path,
                algorithm,
                unitary_matrix=unitary
            )

        print("✓ Test passed: numerical_simulation_dimension_mismatch_raises_error")


def test_compute_unitary_matrix_helper():
    """Test the compute_unitary_matrix helper function."""
    from qhat.analysis.matrix_eigendecomposition import compute_unitary_matrix

    # Create unitary
    unitary = np.eye(4, dtype=complex)
    algorithm = MockAlgorithm(unitary)

    # Compute matrix
    result = compute_unitary_matrix(algorithm)

    # Verify it matches
    assert np.allclose(result, unitary)
    assert result.shape == (4, 4)

    print("✓ Test passed: compute_unitary_matrix_helper")


def test_compute_unitary_matrix_missing_tensor_contract():
    """Test error handling when algorithm lacks tensor_contract method."""
    from qhat.analysis.matrix_eigendecomposition import compute_unitary_matrix

    # Create algorithm without tensor_contract method
    class BadAlgorithm:
        pass

    algorithm = BadAlgorithm()

    with pytest.raises(AttributeError, match="does not have a 'tensor_contract.*' method"):
        compute_unitary_matrix(algorithm)

    print("✓ Test passed: compute_unitary_matrix_missing_tensor_contract")


def test_numerical_simulation_invalid_input_type():
    """Test error handling for invalid input type."""
    unitary = np.eye(4, dtype=complex)
    algorithm = MockAlgorithm(unitary)

    with pytest.raises(ValueError, match="must be a string or list of strings"):
        numerical_simulation(
            123,  # Invalid: not string or list
            algorithm,
            unitary_matrix=unitary
        )

    print("✓ Test passed: numerical_simulation_invalid_input_type")


# =============================================================================
# Integration Tests
# =============================================================================

def test_output_filename_generation():
    """Test that output filenames are generated correctly."""
    with tempfile.TemporaryDirectory() as tmpdir:
        unitary = np.eye(4, dtype=complex)

        # Test different path formats
        test_cases = [
            ("state.npy", "state_final.npy"),
            ("my_state.npy", "my_state_final.npy"),
            (f"{tmpdir}/deep/path/state.npy", f"{tmpdir}/deep/path/state_final.npy"),
        ]

        for input_name, expected_output in test_cases:
            # Create directory if needed
            input_path = input_name if not input_name.startswith(tmpdir) else input_name
            if os.path.dirname(input_path):
                os.makedirs(os.path.dirname(input_path), exist_ok=True)

            # Create input file
            input_state = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
            np.save(input_path, input_state)

            # Run simulation
            algorithm = MockAlgorithm(unitary)
            result = numerical_simulation(
                input_path,
                algorithm,
                unitary_matrix=unitary
            )

            # Verify output filename
            assert result['simulations'][0]['output_file'] == expected_output
            assert os.path.exists(expected_output)

        print("✓ Test passed: output_filename_generation")


def test_complex_superposition_state():
    """Test with a complex superposition state."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a rotation unitary
        theta = np.pi / 4
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        # Single-qubit rotation extended to 2 qubits (acts on first qubit)
        R = np.array([
            [cos_t, -sin_t, 0, 0],
            [sin_t, cos_t, 0, 0],
            [0, 0, cos_t, -sin_t],
            [0, 0, sin_t, cos_t]
        ], dtype=complex)

        # Create superposition input state
        input_path = os.path.join(tmpdir, "input.npy")
        input_state = np.array([1/np.sqrt(2), 0, 1/np.sqrt(2), 0], dtype=complex)
        np.save(input_path, input_state)

        # Run simulation
        algorithm = MockAlgorithm(R)
        result = numerical_simulation(
            input_path,
            algorithm,
            unitary_matrix=R
        )

        # Verify norm is preserved
        assert np.abs(result['simulations'][0]['output_norm'] - 1.0) < 1e-10

        # Verify output matches expected
        output_path = os.path.join(tmpdir, "input_final.npy")
        output_state = np.load(output_path)
        expected = R @ input_state
        assert np.allclose(output_state, expected)

        print("✓ Test passed: complex_superposition_state")


# =============================================================================
# Run all tests
# =============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("Running Numerical Simulation Tests")
    print("="*70 + "\n")

    # State I/O tests
    print("State I/O Tests:")
    test_load_state_numpy_basic()
    test_load_state_converts_real_to_complex()
    test_save_and_load_roundtrip()
    test_get_state_format_from_extension()
    test_invalid_extension_raises_error()
    test_load_nonexistent_file_raises_error()
    test_load_wrong_dimension_not_power_of_2()
    test_load_multidimensional_array_raises_error()

    print("\nNumerical Simulation Tests:")
    test_numerical_simulation_single_input()
    test_numerical_simulation_multiple_inputs()
    test_numerical_simulation_preserves_unitarity()
    test_numerical_simulation_computes_correct_result()
    test_numerical_simulation_dimension_mismatch_raises_error()
    test_compute_unitary_matrix_helper()
    test_compute_unitary_matrix_missing_tensor_contract()
    test_numerical_simulation_invalid_input_type()

    print("\nIntegration Tests:")
    test_output_filename_generation()
    test_complex_superposition_state()

    print("\n" + "="*70)
    print("All tests passed! ✓")
    print("="*70)
