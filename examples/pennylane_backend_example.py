"""Example: Using QHAT with PennyLane backend.

This example shows how to use PennyLane as the backend for QHAT operations,
demonstrating PennyLane-specific features and capabilities.
"""

import logging
import numpy as np
from qhat.analysis.backend import get_backend
from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """Demonstrate PennyLane backend usage."""

    print("=" * 70)
    print("QHAT with PennyLane Backend Example")
    print("=" * 70)

    try:
        # Get PennyLane backend
        backend = get_backend("pennylane", device="default.qubit")
        print(f"\n✓ Loaded {backend.name} backend")
        print(f"  Capabilities: {sorted(backend.capabilities)}")

        # Create a 2-qubit Hamiltonian
        # H = 0.5 * ZZ + 0.3 * XX
        pauli_dict = {
            ((0, 'Z'), (1, 'Z')): 0.5,
            ((0, 'X'), (1, 'X')): 0.3,
        }

        hamiltonian = Hamiltonian(LinearCombinationOfPauliStrings(
            sparse=pauli_dict,
            num_qubits=2
        ))

        print(f"\nHamiltonian: {len(pauli_dict)} Pauli terms on 2 qubits")

        # Encode with first-order Trotter
        print("\n" + "-" * 70)
        print("First-Order Trotterization")
        print("-" * 70)

        pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")

        unitary1 = backend.encode_pauli_trotter(
            pauli_strings=pauli_strings,
            trotter_order="first order",
            evolution_time=1.0,
            num_steps=5,
            num_qubits=2
        )

        print(f"Created unitary: {unitary1}")

        resources1 = unitary1.estimate_resources()
        print(f"\nResources (1st order, 5 steps):")
        print(f"  T gates: {resources1.t_gates}")
        print(f"  Rotation gates: {resources1.rotation_gates}")
        print(f"  Operations: {resources1.backend_specific.get('pennylane_num_operations', 'N/A')}")

        # Encode with second-order Trotter
        print("\n" + "-" * 70)
        print("Second-Order Trotterization")
        print("-" * 70)

        unitary2 = backend.encode_pauli_trotter(
            pauli_strings=pauli_strings,
            trotter_order="second order",
            evolution_time=1.0,
            num_steps=5,
            num_qubits=2
        )

        resources2 = unitary2.estimate_resources()
        print(f"Resources (2nd order, 5 steps):")
        print(f"  T gates: {resources2.t_gates}")
        print(f"  Rotation gates: {resources2.rotation_gates}")
        print(f"  Operations: {resources2.backend_specific.get('pennylane_num_operations', 'N/A')}")

        # Compare resource usage
        print("\n" + "-" * 70)
        print("Comparison")
        print("-" * 70)
        print(f"2nd order uses {resources2.rotation_gates / resources1.rotation_gates:.1f}x "
              f"more rotation gates than 1st order")

        # Get matrices and compare
        print("\n" + "-" * 70)
        print("Matrix Comparison")
        print("-" * 70)

        mat1 = unitary1.to_matrix()
        mat2 = unitary2.to_matrix()

        print(f"1st order matrix norm: {np.linalg.norm(mat1):.6f}")
        print(f"2nd order matrix norm: {np.linalg.norm(mat2):.6f}")

        # Check unitarity
        def check_unitarity(mat, name):
            identity = mat.conj().T @ mat
            error = np.linalg.norm(identity - np.eye(mat.shape[0]))
            print(f"{name} unitarity error: {error:.2e}")

        check_unitarity(mat1, "1st order")
        check_unitarity(mat2, "2nd order")

        # Difference between formulas
        diff = np.linalg.norm(mat1 - mat2, 'fro')
        print(f"\nDifference between 1st and 2nd order: {diff:.4f}")
        print("(Higher-order formulas should be more accurate for same step count)")

        # Test controlled operations
        print("\n" + "-" * 70)
        print("Controlled Operations")
        print("-" * 70)

        controlled_unitary = unitary1.controlled(num_controls=1)
        print(f"Created controlled unitary: {controlled_unitary}")
        print(f"  Original qubits: {unitary1.num_qubits}")
        print(f"  Controlled qubits: {controlled_unitary.num_qubits}")

        # Test adjoint
        print("\n" + "-" * 70)
        print("Adjoint Operation")
        print("-" * 70)

        adjoint_unitary = unitary1.adjoint()
        print(f"Created adjoint: {adjoint_unitary}")

        # Verify U† U = I
        mat_adj = adjoint_unitary.to_matrix()
        product = mat_adj @ mat1
        identity_error = np.linalg.norm(product - np.eye(4))
        print(f"U† U identity error: {identity_error:.2e}")
        if identity_error < 1e-10:
            print("✓ Adjoint is correct")

        print("\n" + "=" * 70)
        print("PennyLane backend example complete!")
        print("=" * 70)

    except ImportError as e:
        print(f"\n✗ PennyLane not available: {e}")
        print("\nTo install PennyLane:")
        print("  pip install pennylane")
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
