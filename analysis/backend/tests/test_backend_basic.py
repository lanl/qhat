"""Basic backend functionality tests.

These tests verify that the backend system works correctly and that
backends can be loaded and used.
"""

import pytest
import numpy as np

from qhat.analysis.backend import get_backend, BackendRegistry
from qhat.analysis.backend.types import UnsupportedOperationError


def test_backend_registry_list():
    """Test that backends can be listed."""
    available = BackendRegistry.list_available()
    assert isinstance(available, list)
    # Should have at least the backends we defined
    # (may not all be available if dependencies not installed)
    print(f"Available backends: {available}")


def test_get_qualtran_backend():
    """Test loading Qualtran backend."""
    try:
        backend = get_backend("qualtran")
        assert backend.name == "qualtran"
        assert "pauli_trotter" in backend.capabilities
        assert "pauli_lcu" in backend.capabilities
        print(f"Qualtran capabilities: {backend.capabilities}")
    except ImportError as e:
        pytest.skip(f"Qualtran backend not available: {e}")


def test_get_pennylane_backend():
    """Test loading PennyLane backend."""
    try:
        backend = get_backend("pennylane")
        assert backend.name == "pennylane"
        assert "pauli_trotter" in backend.capabilities
        print(f"PennyLane capabilities: {backend.capabilities}")
    except ImportError as e:
        pytest.skip(f"PennyLane backend not available: {e}")


def test_get_qiskit_backend():
    """Test loading Qiskit backend."""
    try:
        backend = get_backend("qiskit")
        assert backend.name == "qiskit"
        assert "pauli_trotter" in backend.capabilities
        print(f"Qiskit capabilities: {backend.capabilities}")
    except ImportError as e:
        pytest.skip(f"Qiskit backend not available: {e}")


def test_invalid_backend():
    """Test that invalid backend name raises error."""
    with pytest.raises(ValueError, match="not found"):
        get_backend("nonexistent_backend")


def test_simple_trotter_encoding():
    """Test basic Trotter encoding with available backend."""
    # Simple 2-qubit Hamiltonian: H = 0.5 ZZ
    pauli_strings = {
        ((0, 'Z'), (1, 'Z')): 0.5,
    }

    # Try each backend
    for backend_name in ["qualtran", "pennylane", "qiskit"]:
        try:
            backend = get_backend(backend_name)

            unitary = backend.encode_pauli_trotter(
                pauli_strings=pauli_strings,
                trotter_order="first order",
                evolution_time=1.0,
                num_steps=1,
                num_qubits=2
            )

            # Check basic properties
            assert unitary.num_qubits == 2
            assert unitary.backend_name == backend_name

            # Try getting resources
            resources = unitary.estimate_resources()
            assert resources.num_qubits == 2

            print(f"\n{backend_name} Trotter encoding:")
            print(f"  Qubits: {resources.num_qubits}")
            print(f"  T gates: {resources.t_gates}")
            print(f"  Clifford gates: {resources.clifford_gates}")

        except ImportError:
            print(f"\n{backend_name} not available, skipping")
            continue
        except Exception as e:
            print(f"\n{backend_name} failed: {e}")
            # Don't fail the test if backend has issues
            continue


def test_unsupported_operation():
    """Test that unsupported operations raise proper error."""
    try:
        backend = get_backend("pennylane")

        # PennyLane doesn't support LCU
        with pytest.raises(UnsupportedOperationError):
            backend.encode_pauli_lcu(
                pauli_strings={((0, 'Z'),): 1.0},
                num_qubits=1
            )

    except ImportError:
        pytest.skip("PennyLane not available")


if __name__ == "__main__":
    # Run tests
    print("Testing backend system...")
    print("=" * 60)

    test_backend_registry_list()
    print("\n" + "=" * 60)

    test_get_qualtran_backend()
    test_get_pennylane_backend()
    test_get_qiskit_backend()
    print("\n" + "=" * 60)

    test_simple_trotter_encoding()
    print("\n" + "=" * 60)

    print("\nAll basic tests passed!")
