"""
Tests for HamLib HDF5 file loading functionality.
"""

import os
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pytest

from qhat.qre.qre_hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings, load_hamlib_hdf5
from qhat.qre.qre_types import GeneralConfiguration, HamiltonianConfiguration


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self):
        self.verbose = False

    def log(self, msg):
        """Mock logging."""
        if self.verbose:
            print(msg)

    def log_verbose(self, msg):
        """Mock verbose logging."""
        if self.verbose:
            print(f"[VERBOSE] {msg}")


def create_test_hamlib_hdf5(filename, include_metadata=True):
    """
    Create a test HamLib HDF5 file with a simple Hamiltonian.

    Creates a 2-qubit Hamiltonian: 1.5 * X0 Z1 - 0.5 * Y0 Y1 + 2.0 * I
    """
    # HamLib format string
    hamlib_string = "(1.5+0j) [X0 Z1] +\n(-0.5+0j) [Y0 Y1] +\n(2.0+0j) []"

    with h5py.File(filename, 'w') as f:
        dataset = f.create_dataset('test_hamiltonian', data=hamlib_string.encode('utf-8'))

        # Add metadata (HamLib v1.1+)
        if include_metadata:
            dataset.attrs['nqubits'] = 2
            dataset.attrs['terms'] = 3
            dataset.attrs['one_norm'] = 4.0  # |1.5| + |-0.5| + |2.0|


def test_load_simple_hamlib_file():
    """Test loading a simple HamLib HDF5 file."""
    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Create test file
        create_test_hamlib_hdf5(tmp_path)

        # Create mock configuration
        config_general = MockConfig()
        config_hamiltonian = MockConfig()
        config_hamiltonian.filename = tmp_path
        config_hamiltonian.hdf5_key = 'test_hamiltonian'

        # Load the Hamiltonian
        H = load_hamlib_hdf5(config_general, config_hamiltonian)

        # Verify it's a Hamiltonian object
        assert isinstance(H, Hamiltonian)

        # Get Pauli strings
        pauli_dict = H.get_all_pauli_strings(return_as="tuples")

        # Expected terms
        expected = {
            ((0, 'X'), (1, 'Z')): 1.5,
            ((0, 'Y'), (1, 'Y')): -0.5,
            tuple(): 2.0  # Identity term
        }

        # Check number of terms
        assert len(pauli_dict) == 3, f"Expected 3 terms, got {len(pauli_dict)}"

        # Check each term
        for pauli, expected_coef in expected.items():
            assert pauli in pauli_dict, f"Missing Pauli term: {pauli}"
            actual_coef = pauli_dict[pauli]
            assert abs(actual_coef - expected_coef) < 1e-10, \
                f"Coefficient mismatch for {pauli}: expected {expected_coef}, got {actual_coef}"

        # Check number of qubits
        assert H.num_qubits() == 2

        print("✓ Simple HamLib file loaded successfully")

    finally:
        # Clean up
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_load_hamlib_autodetect_key():
    """Test auto-detection of HDF5 key when not specified."""
    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Create test file
        create_test_hamlib_hdf5(tmp_path)

        # Create mock configuration WITHOUT specifying hdf5_key
        config_general = MockConfig()
        config_hamiltonian = MockConfig()
        config_hamiltonian.filename = tmp_path

        # Load the Hamiltonian (should auto-detect the key)
        H = load_hamlib_hdf5(config_general, config_hamiltonian)

        # Verify it loaded correctly
        assert isinstance(H, Hamiltonian)
        assert H.num_qubits() == 2

        print("✓ Auto-detection of HDF5 key works")

    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_load_hamlib_with_metadata():
    """Test that metadata is read correctly from HamLib v1.1+ files."""
    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Create test file with metadata
        create_test_hamlib_hdf5(tmp_path, include_metadata=True)

        # Verify metadata was written
        with h5py.File(tmp_path, 'r') as f:
            dataset = f['test_hamiltonian']
            assert 'nqubits' in dataset.attrs
            assert dataset.attrs['nqubits'] == 2

        # Load the file
        config_general = MockConfig()
        config_hamiltonian = MockConfig()
        config_hamiltonian.filename = tmp_path
        config_hamiltonian.hdf5_key = 'test_hamiltonian'

        H = load_hamlib_hdf5(config_general, config_hamiltonian)

        # Verify loaded correctly
        assert isinstance(H, Hamiltonian)
        assert H.num_qubits() == 2

        print("✓ Metadata handling works")

    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_identity_only_hamiltonian():
    """Test a Hamiltonian with only an identity term."""
    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Create file with only identity
        hamlib_string = "(3.14+0j) []"

        with h5py.File(tmp_path, 'w') as f:
            f.create_dataset('identity_ham', data=hamlib_string.encode('utf-8'))

        # Load it
        config_general = MockConfig()
        config_hamiltonian = MockConfig()
        config_hamiltonian.filename = tmp_path
        config_hamiltonian.hdf5_key = 'identity_ham'

        H = load_hamlib_hdf5(config_general, config_hamiltonian)

        # Check
        pauli_dict = H.get_all_pauli_strings(return_as="tuples")
        assert len(pauli_dict) == 1
        assert tuple() in pauli_dict
        assert abs(pauli_dict[tuple()] - 3.14) < 1e-10

        print("✓ Identity-only Hamiltonian works")

    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_real_hamlib_file():
    """Test loading a real HamLib file (H2 molecule with STO-6G basis)."""
    # Path to the real HamLib file
    test_file = Path(__file__).parent / "hamlib_h2_sto6g.hdf5"

    if not test_file.exists():
        print("⚠ Skipping real HamLib file test (file not found)")
        return

    # Load using Jordan-Wigner mapping (ham_JW)
    config_general = MockConfig()
    config_hamiltonian = MockConfig()
    config_hamiltonian.filename = str(test_file)
    config_hamiltonian.hdf5_key = 'ham_JW'

    H = load_hamlib_hdf5(config_general, config_hamiltonian)

    # Verify basic properties
    assert isinstance(H, Hamiltonian)
    assert H.num_qubits() == 4, f"Expected 4 qubits, got {H.num_qubits()}"

    # Check metadata matches
    with h5py.File(test_file, 'r') as f:
        dataset = f['ham_JW']
        expected_terms = dataset.attrs['terms']
        expected_nqubits = dataset.attrs['nqubits']
        expected_one_norm = dataset.attrs['one_norm']

    pauli_dict = H.get_all_pauli_strings(return_as="tuples")
    assert len(pauli_dict) == expected_terms, \
        f"Expected {expected_terms} terms, got {len(pauli_dict)}"

    # Verify some known terms from the file
    # Identity term: -0.3371448491762382
    assert tuple() in pauli_dict
    assert abs(pauli_dict[tuple()] - (-0.3371448491762382)) < 1e-10

    # Check Z0 term: 0.1389276591112677
    z0_term = ((0, 'Z'),)
    assert z0_term in pauli_dict
    assert abs(pauli_dict[z0_term] - 0.1389276591112677) < 1e-10

    print("✓ Real HamLib file (H2 STO-6G) loaded successfully")


if __name__ == "__main__":
    print("Running HamLib loader tests...\n")

    test_load_simple_hamlib_file()
    test_load_hamlib_autodetect_key()
    test_load_hamlib_with_metadata()
    test_identity_only_hamiltonian()
    test_real_hamlib_file()

    print("\n✓ All tests passed!")
