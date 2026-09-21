"""Example: Comparing different backends in QHAT.

This example demonstrates how to use QHAT with different quantum computing
backends (Qualtran, PennyLane, Qiskit) and compare their results.
"""

import logging
import numpy as np
from qhat.analysis.backend import get_backend, BackendRegistry
from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def create_simple_hamiltonian():
    """Create a simple 3-qubit Hamiltonian for testing.

    H = 0.5 * Z0*Z1 + 0.3 * X0*X1 + 0.2 * Y1*Y2
    """
    pauli_dict = {
        ((0, 'Z'), (1, 'Z')): 0.5,
        ((0, 'X'), (1, 'X')): 0.3,
        ((1, 'Y'), (2, 'Y')): 0.2,
    }

    return Hamiltonian(LinearCombinationOfPauliStrings(
        sparse=pauli_dict,
        num_qubits=3
    ))


def main():
    """Compare backends for Trotterization."""

    print("=" * 70)
    print("QHAT Backend Comparison Example")
    print("=" * 70)

    # Create Hamiltonian
    hamiltonian = create_simple_hamiltonian()
    pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")

    print(f"\nHamiltonian: {len(pauli_strings)} Pauli terms on {hamiltonian.num_qubits()} qubits")

    # List available backends
    available = BackendRegistry.list_available()
    print(f"\nAvailable backends: {available}")

    # Parameters for Trotterization
    evolution_time = 1.0
    num_steps = 10
    trotter_order = "second order"

    results = {}

    # Try each backend
    for backend_name in ["qualtran", "pennylane", "qiskit"]:
        print(f"\n{'-' * 70}")
        print(f"Testing {backend_name.upper()} backend")
        print(f"{'-' * 70}")

        try:
            # Get backend
            backend = get_backend(backend_name)
            print(f"✓ Loaded {backend.name} backend")
            print(f"  Capabilities: {sorted(backend.capabilities)}")

            # Encode Hamiltonian
            print(f"\nEncoding Hamiltonian with {trotter_order} Trotterization...")
            unitary = backend.encode_pauli_trotter(
                pauli_strings=pauli_strings,
                trotter_order=trotter_order,
                evolution_time=evolution_time,
                num_steps=num_steps,
                num_qubits=hamiltonian.num_qubits()
            )

            print(f"✓ Created unitary: {unitary}")

            # Get resources
            print(f"\nResource Estimation:")
            resources = unitary.estimate_resources()
            print(f"  Qubits: {resources.num_qubits}")
            print(f"  T gates: {resources.t_gates}")
            print(f"  Clifford gates: {resources.clifford_gates}")
            print(f"  Rotation gates: {resources.rotation_gates}")
            if resources.depth:
                print(f"  Depth: {resources.depth}")
            if resources.two_qubit_gates:
                print(f"  Two-qubit gates: {resources.two_qubit_gates}")

            # Try to get matrix (for small systems)
            if hamiltonian.num_qubits() <= 10:
                print(f"\nMatrix Representation:")
                try:
                    matrix = unitary.to_matrix()
                    print(f"  Matrix shape: {matrix.shape}")
                    print(f"  Matrix norm: {np.linalg.norm(matrix):.6f}")

                    # Check unitarity: U† U = I
                    identity_check = matrix.conj().T @ matrix
                    error = np.linalg.norm(identity_check - np.eye(matrix.shape[0]))
                    print(f"  Unitarity check (||U†U - I||): {error:.2e}")

                    results[backend_name] = {
                        'resources': resources,
                        'matrix': matrix
                    }
                except Exception as e:
                    print(f"  Matrix conversion failed: {e}")
                    results[backend_name] = {
                        'resources': resources,
                        'matrix': None
                    }
            else:
                results[backend_name] = {
                    'resources': resources,
                    'matrix': None
                }

        except ImportError as e:
            print(f"✗ {backend_name} not available: {e}")
        except Exception as e:
            print(f"✗ {backend_name} failed: {e}")
            import traceback
            traceback.print_exc()

    # Compare results
    if len(results) >= 2:
        print(f"\n{'=' * 70}")
        print("COMPARISON")
        print(f"{'=' * 70}")

        # Compare resource estimates
        print("\nResource Estimates:")
        print(f"{'Backend':<12} {'T Gates':<12} {'Clifford':<12} {'Rotations':<12}")
        print(f"{'-' * 50}")
        for backend_name, data in results.items():
            res = data['resources']
            print(f"{backend_name:<12} {res.t_gates:<12} {res.clifford_gates:<12} {res.rotation_gates:<12}")

        # Compare matrices (if available)
        matrices = {name: data['matrix'] for name, data in results.items() if data['matrix'] is not None}

        if len(matrices) >= 2:
            print("\nMatrix Comparison:")
            names = list(matrices.keys())
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    name1, name2 = names[i], names[j]
                    mat1, mat2 = matrices[name1], matrices[name2]

                    # Compute Frobenius norm of difference
                    diff_norm = np.linalg.norm(mat1 - mat2, 'fro')
                    rel_diff = diff_norm / np.linalg.norm(mat1, 'fro')

                    print(f"  {name1} vs {name2}:")
                    print(f"    Absolute difference: {diff_norm:.2e}")
                    print(f"    Relative difference: {rel_diff:.2e}")

                    if rel_diff < 1e-10:
                        print(f"    ✓ Matrices are numerically identical")
                    elif rel_diff < 1e-6:
                        print(f"    ✓ Matrices are very close")
                    else:
                        print(f"    ⚠ Matrices differ significantly")

    print(f"\n{'=' * 70}")
    print("Example complete!")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
