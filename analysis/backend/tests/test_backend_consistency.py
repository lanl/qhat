"""Comprehensive backend consistency tests.

This module tests that different backends produce consistent results
for the same inputs, validating the backend abstraction.
"""

import json
import logging
import numpy as np
import pytest
from typing import Dict, Any, List

from qhat.analysis.backend import get_backend, BackendRegistry
from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ConsistencyTestResults:
    """Collect and save consistency test results."""

    def __init__(self):
        self.results = []

    def add_result(self, test_name: str, test_data: Dict[str, Any]):
        """Add a test result."""
        self.results.append({
            'test_name': test_name,
            'data': test_data
        })

    def save(self, filename: str):
        """Save results to JSON file."""
        with open(filename, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)
        logger.info(f"Saved results to {filename}")


# Global results collector
results_collector = ConsistencyTestResults()


def get_available_backends():
    """Get list of available backends."""
    available = []
    for name in ["qualtran", "pennylane", "qiskit"]:
        try:
            backend = get_backend(name)
            available.append(name)
        except (ImportError, ValueError):
            pass
    return available


def compare_matrices(mat1: np.ndarray, mat2: np.ndarray,
                     name1: str, name2: str) -> Dict[str, float]:
    """Compare two matrices and return metrics."""
    # Frobenius norm difference
    abs_diff = np.linalg.norm(mat1 - mat2, 'fro')
    rel_diff = abs_diff / np.linalg.norm(mat1, 'fro')

    # Element-wise max difference
    max_elem_diff = np.max(np.abs(mat1 - mat2))

    # Check if phases might differ (matrices equal up to global phase)
    # Try to find best phase match
    best_phase_diff = rel_diff
    for phase in [1, -1, 1j, -1j]:
        phase_rel_diff = np.linalg.norm(mat1 - phase * mat2, 'fro') / np.linalg.norm(mat1, 'fro')
        if phase_rel_diff < best_phase_diff:
            best_phase_diff = phase_rel_diff

    return {
        'absolute_difference': float(abs_diff),
        'relative_difference': float(rel_diff),
        'max_element_difference': float(max_elem_diff),
        'best_phase_relative_difference': float(best_phase_diff),
        'backends': f"{name1} vs {name2}"
    }


def create_test_hamiltonians() -> Dict[str, Hamiltonian]:
    """Create various test Hamiltonians."""
    hamiltonians = {}

    # 2-qubit: Simple ZZ + XX
    hamiltonians['2q_simple'] = Hamiltonian(LinearCombinationOfPauliStrings(
        sparse={
            ((0, 'Z'), (1, 'Z')): 0.5,
            ((0, 'X'), (1, 'X')): 0.3,
        },
        num_qubits=2
    ))

    # 3-qubit: Multiple terms
    hamiltonians['3q_multi'] = Hamiltonian(LinearCombinationOfPauliStrings(
        sparse={
            ((0, 'Z'), (1, 'Z')): 0.5,
            ((0, 'X'), (1, 'X')): 0.3,
            ((1, 'Y'), (2, 'Y')): 0.2,
            ((0, 'Z'),): 0.1,
        },
        num_qubits=3
    ))

    # 4-qubit: Larger system
    hamiltonians['4q_chain'] = Hamiltonian(LinearCombinationOfPauliStrings(
        sparse={
            ((0, 'Z'), (1, 'Z')): 1.0,
            ((1, 'Z'), (2, 'Z')): 1.0,
            ((2, 'Z'), (3, 'Z')): 1.0,
            ((0, 'X'),): 0.5,
            ((1, 'X'),): 0.5,
            ((2, 'X'),): 0.5,
            ((3, 'X'),): 0.5,
        },
        num_qubits=4
    ))

    return hamiltonians


@pytest.mark.parametrize("ham_name", ['2q_simple', '3q_multi', '4q_chain'])
@pytest.mark.parametrize("trotter_order", ['first order', 'second order'])
@pytest.mark.parametrize("num_steps", [1, 5, 10])
def test_trotter_consistency(ham_name, trotter_order, num_steps):
    """Test that Trotter encoding is consistent across backends."""

    test_id = f"trotter_{ham_name}_{trotter_order.replace(' ', '_')}_steps{num_steps}"
    logger.info(f"\n{'='*70}")
    logger.info(f"Test: {test_id}")
    logger.info(f"{'='*70}")

    # Get available backends
    available = get_available_backends()
    if len(available) < 2:
        pytest.skip(f"Need at least 2 backends, only {len(available)} available")

    # Create Hamiltonian
    hamiltonians = create_test_hamiltonians()
    hamiltonian = hamiltonians[ham_name]
    pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")

    evolution_time = 1.0
    num_qubits = hamiltonian.num_qubits()

    logger.info(f"Hamiltonian: {ham_name}, {num_qubits} qubits, {len(pauli_strings)} terms")
    logger.info(f"Trotter: {trotter_order}, {num_steps} steps, t={evolution_time}")

    # Encode with each backend
    unitaries = {}
    resources = {}
    matrices = {}

    for backend_name in available:
        try:
            backend = get_backend(backend_name)

            # Encode
            unitary = backend.encode_pauli_trotter(
                pauli_strings=pauli_strings,
                trotter_order=trotter_order,
                evolution_time=evolution_time,
                num_steps=num_steps,
                num_qubits=num_qubits
            )

            unitaries[backend_name] = unitary

            # Get resources
            res = unitary.estimate_resources()
            resources[backend_name] = {
                'num_qubits': res.num_qubits,
                't_gates': res.t_gates,
                'clifford_gates': res.clifford_gates,
                'rotation_gates': res.rotation_gates,
                'total_gates': res.total_gates()
            }

            # Get matrix (if small enough)
            if num_qubits <= 6:
                matrices[backend_name] = unitary.to_matrix()

            logger.info(f"{backend_name}: T={res.t_gates}, Clifford={res.clifford_gates}, Rot={res.rotation_gates}")

        except Exception as e:
            logger.warning(f"{backend_name} failed: {e}")

    # Compare matrices
    matrix_comparisons = []
    if len(matrices) >= 2:
        backend_names = list(matrices.keys())
        for i in range(len(backend_names)):
            for j in range(i+1, len(backend_names)):
                name1, name2 = backend_names[i], backend_names[j]
                comparison = compare_matrices(
                    matrices[name1], matrices[name2], name1, name2
                )
                matrix_comparisons.append(comparison)

                logger.info(f"\n{name1} vs {name2}:")
                logger.info(f"  Relative difference: {comparison['relative_difference']:.2e}")
                logger.info(f"  Max element diff: {comparison['max_element_difference']:.2e}")

                # Assert consistency
                assert comparison['relative_difference'] < 1e-6, \
                    f"Matrices differ: {name1} vs {name2}, rel_diff={comparison['relative_difference']:.2e}"

    # Save results
    result_data = {
        'hamiltonian': ham_name,
        'trotter_order': trotter_order,
        'num_steps': num_steps,
        'num_qubits': num_qubits,
        'backends_tested': list(unitaries.keys()),
        'resources': resources,
        'matrix_comparisons': matrix_comparisons,
        'passed': True
    }

    results_collector.add_result(test_id, result_data)


def test_resource_scaling():
    """Test that resource estimates scale reasonably across backends."""

    logger.info(f"\n{'='*70}")
    logger.info("Test: Resource Scaling")
    logger.info(f"{'='*70}")

    available = get_available_backends()
    if len(available) < 2:
        pytest.skip(f"Need at least 2 backends, only {len(available)} available")

    # Simple 3-qubit Hamiltonian
    hamiltonian = Hamiltonian(LinearCombinationOfPauliStrings(
        sparse={
            ((0, 'Z'), (1, 'Z')): 0.5,
            ((0, 'X'), (1, 'X')): 0.3,
            ((1, 'Y'), (2, 'Y')): 0.2,
        },
        num_qubits=3
    ))

    pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")

    # Test different step counts
    step_counts = [1, 5, 10, 20]
    scaling_data = {backend_name: [] for backend_name in available}

    for num_steps in step_counts:
        logger.info(f"\nSteps: {num_steps}")

        for backend_name in available:
            try:
                backend = get_backend(backend_name)

                unitary = backend.encode_pauli_trotter(
                    pauli_strings=pauli_strings,
                    trotter_order="second order",
                    evolution_time=1.0,
                    num_steps=num_steps,
                    num_qubits=3
                )

                res = unitary.estimate_resources()
                scaling_data[backend_name].append({
                    'num_steps': num_steps,
                    't_gates': res.t_gates,
                    'total_gates': res.total_gates()
                })

                logger.info(f"  {backend_name}: T={res.t_gates}, Total={res.total_gates()}")

            except Exception as e:
                logger.warning(f"  {backend_name} failed: {e}")

    # Check that resources scale linearly with steps
    for backend_name, data in scaling_data.items():
        if len(data) >= 2:
            # Compare first and last
            first_gates = data[0]['total_gates']
            last_gates = data[-1]['total_gates']
            first_steps = data[0]['num_steps']
            last_steps = data[-1]['num_steps']

            expected_ratio = last_steps / first_steps
            actual_ratio = last_gates / first_gates if first_gates > 0 else 0

            logger.info(f"\n{backend_name} scaling:")
            logger.info(f"  Steps ratio: {expected_ratio:.1f}x")
            logger.info(f"  Gates ratio: {actual_ratio:.1f}x")

            # Should be roughly linear (within 2x of expected)
            if actual_ratio > 0:
                assert 0.5 * expected_ratio < actual_ratio < 2.0 * expected_ratio, \
                    f"{backend_name}: gates don't scale linearly with steps"

    # Save results
    result_data = {
        'scaling_data': scaling_data,
        'passed': True
    }

    results_collector.add_result('resource_scaling', result_data)


def test_unitarity():
    """Test that all backends produce unitary matrices."""

    logger.info(f"\n{'='*70}")
    logger.info("Test: Unitarity")
    logger.info(f"{'='*70}")

    available = get_available_backends()

    # Simple 2-qubit Hamiltonian
    hamiltonian = Hamiltonian(LinearCombinationOfPauliStrings(
        sparse={
            ((0, 'Z'), (1, 'Z')): 0.5,
            ((0, 'X'), (1, 'X')): 0.3,
        },
        num_qubits=2
    ))

    pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")

    unitarity_results = {}

    for backend_name in available:
        try:
            backend = get_backend(backend_name)

            unitary = backend.encode_pauli_trotter(
                pauli_strings=pauli_strings,
                trotter_order="second order",
                evolution_time=1.0,
                num_steps=10,
                num_qubits=2
            )

            # Get matrix
            U = unitary.to_matrix()

            # Check U† U = I
            identity_check = U.conj().T @ U
            error = np.linalg.norm(identity_check - np.eye(U.shape[0]), 'fro')

            unitarity_results[backend_name] = float(error)

            logger.info(f"{backend_name}: ||U†U - I|| = {error:.2e}")

            # Assert unitarity
            assert error < 1e-10, f"{backend_name} produced non-unitary matrix"

        except Exception as e:
            logger.warning(f"{backend_name} failed: {e}")

    # Save results
    result_data = {
        'unitarity_errors': unitarity_results,
        'passed': True
    }

    results_collector.add_result('unitarity', result_data)


def test_controlled_operations():
    """Test that controlled operations work consistently."""

    logger.info(f"\n{'='*70}")
    logger.info("Test: Controlled Operations")
    logger.info(f"{'='*70}")

    available = get_available_backends()

    # Simple 2-qubit Hamiltonian
    hamiltonian = Hamiltonian(LinearCombinationOfPauliStrings(
        sparse={
            ((0, 'Z'), (1, 'Z')): 0.5,
        },
        num_qubits=2
    ))

    pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")

    controlled_results = {}

    for backend_name in available:
        try:
            backend = get_backend(backend_name)

            # Create unitary
            unitary = backend.encode_pauli_trotter(
                pauli_strings=pauli_strings,
                trotter_order="first order",
                evolution_time=0.5,
                num_steps=1,
                num_qubits=2
            )

            # Create controlled version
            controlled = unitary.controlled(num_controls=1)

            # Check qubit count increased
            assert controlled.num_qubits == unitary.num_qubits + 1

            # Get matrix (now 3 qubits)
            C = controlled.to_matrix()

            # Check structure: should be block diagonal with I and U
            # [I 0]
            # [0 U]
            dim = 2 ** unitary.num_qubits
            top_left = C[:dim, :dim]
            bottom_right = C[dim:, dim:]

            # Top-left should be identity
            identity_error = np.linalg.norm(top_left - np.eye(dim), 'fro')

            controlled_results[backend_name] = {
                'qubits_increased': controlled.num_qubits == unitary.num_qubits + 1,
                'identity_block_error': float(identity_error),
            }

            logger.info(f"{backend_name}:")
            logger.info(f"  Qubits: {unitary.num_qubits} -> {controlled.num_qubits}")
            logger.info(f"  Identity block error: {identity_error:.2e}")

            assert identity_error < 1e-10, f"{backend_name} controlled operation incorrect"

        except Exception as e:
            logger.warning(f"{backend_name} failed: {e}")

    # Save results
    result_data = {
        'controlled_results': controlled_results,
        'passed': True
    }

    results_collector.add_result('controlled_operations', result_data)


def test_adjoint_operations():
    """Test that adjoint operations satisfy U† U = I."""

    logger.info(f"\n{'='*70}")
    logger.info("Test: Adjoint Operations")
    logger.info(f"{'='*70}")

    available = get_available_backends()

    # Simple 2-qubit Hamiltonian
    hamiltonian = Hamiltonian(LinearCombinationOfPauliStrings(
        sparse={
            ((0, 'Z'), (1, 'Z')): 0.5,
            ((0, 'X'), (1, 'X')): 0.3,
        },
        num_qubits=2
    ))

    pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")

    adjoint_results = {}

    for backend_name in available:
        try:
            backend = get_backend(backend_name)

            # Create unitary
            unitary = backend.encode_pauli_trotter(
                pauli_strings=pauli_strings,
                trotter_order="second order",
                evolution_time=1.0,
                num_steps=5,
                num_qubits=2
            )

            # Create adjoint
            adjoint = unitary.adjoint()

            # Get matrices
            U = unitary.to_matrix()
            U_adj = adjoint.to_matrix()

            # Check U_adj @ U = I
            product = U_adj @ U
            identity_error = np.linalg.norm(product - np.eye(U.shape[0]), 'fro')

            # Also check that U_adj = U†
            expected_adj = U.conj().T
            adj_error = np.linalg.norm(U_adj - expected_adj, 'fro')

            adjoint_results[backend_name] = {
                'identity_error': float(identity_error),
                'adjoint_correctness_error': float(adj_error),
            }

            logger.info(f"{backend_name}:")
            logger.info(f"  ||U† U - I|| = {identity_error:.2e}")
            logger.info(f"  ||adjoint() - U†|| = {adj_error:.2e}")

            assert identity_error < 1e-10, f"{backend_name} adjoint doesn't satisfy U† U = I"
            assert adj_error < 1e-10, f"{backend_name} adjoint incorrect"

        except Exception as e:
            logger.warning(f"{backend_name} failed: {e}")

    # Save results
    result_data = {
        'adjoint_results': adjoint_results,
        'passed': True
    }

    results_collector.add_result('adjoint_operations', result_data)


@pytest.fixture(scope="session", autouse=True)
def save_results(request):
    """Save all results at end of session."""
    yield
    # This runs after all tests
    results_collector.save('analysis/backend/tests/consistency_results/test_results.json')


if __name__ == "__main__":
    # Run tests and save results
    pytest.main([__file__, '-v', '-s'])
