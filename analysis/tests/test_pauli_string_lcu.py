"""
Test construction of Pauli String LCU block encodings.
"""

import numpy as np
import pytest

from pyLIQTR.ProblemInstances.ProblemInstance import ProblemInstance

from qhat.analysis.config_types import (
    GeneralConfiguration,
    GeneralConfigurationUser,
)
from qhat.analysis.hamiltonian import (
    Hamiltonian,
    LinearCombinationOfPauliStrings,
)
from qhat.analysis.unitary import PauliStringLCU, PyLIQTRPauliStringLCU


def test_pauli_string_lcu1():
    """Test converting a small, simple Hamiltonian to Pauli String LCU
       block encoding."""
    # Create a simple 1-qubit Hamiltonian
    pauli_data = {
        'X': 1.0,
        'Z': 1.0,
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=1, dense=pauli_data)
    H = Hamiltonian(lcps)

    # Convert Hamiltonian to matrix
    H_matrix = H.to_matrix(memory_threshold_gb=1.0)

    # Create unitary operator and convert to matrix
    unitaryop = PauliStringLCU(H, 'AS', probability_eps=0.5)
    unitarymx = unitaryop.tensor_contract()

    # Verify that the upper corner of our unitary is equal to the
    # original Hamiltonian matrix, scaled
    np.testing.assert_array_almost_equal(unitarymx[0:2,0:2], 0.5 * H_matrix)


def test_pauli_string_lcu2():
    """Test converting a small, slightly more complex Hamiltonian to
       Pauli String LCU block encoding."""
    # Create a simple 2-qubit Hamiltonian
    pauli_data = {
        'XX': 1.0,
        'YZ': 1.0,
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    H = Hamiltonian(lcps)

    # Convert Hamiltonian to matrix
    H_matrix = H.to_matrix(memory_threshold_gb=1.0)

    # Create unitary operator and convert to matrix
    unitaryop = PauliStringLCU(H, 'AS', probability_eps=0.5)
    unitarymx = unitaryop.tensor_contract()

    # Verify that the upper corner of our unitary is equal to the
    # original Hamiltonian matrix, scaled
    np.testing.assert_array_almost_equal(unitarymx[0:4,0:4], 0.5 * H_matrix)


# Define a wrapper class to make a QHAT PauliString Hamiltonian behave
# like a PyLIQTR ProblemInstance
class PauliStringInstance(Hamiltonian, ProblemInstance):
    def __init__(self, hamiltonian):
        super().__init__(hamiltonian)

    def __str__(self):
        return str(self.get_all_pauli_strings(return_as='strings'))

    def n_terms(self, **kwargs):
        return len(self.get_all_pauli_strings())

    def n_qubits(self):
        return self.num_qubits()

    def yield_PauliLCU_Info(self, do_pad=0, return_as='strings'):
        for t in self.get_all_pauli_strings(return_as=return_as).items():
            yield t


def test_pauli_string_lcu_pyliqtr1():
    """Test converting a small, simple Hamiltonian to Pauli String LCU
       block encoding."""
    # Create a simple 1-qubit Hamiltonian
    pauli_data = {
        'X': 1.0,
        'Z': 1.0,
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=1, dense=pauli_data)
    inst = PauliStringInstance(lcps)

    # Convert Hamiltonian to matrix
    H_matrix = inst.to_matrix(memory_threshold_gb=1.0)

    # Create unitary operator and convert to matrix
    unitaryop = PauliStringLCU(inst, 'AS', probability_eps=0.5)
    unitarymx = unitaryop.tensor_contract()

    # Verify that the upper corner of our unitary is equal to the
    # original Hamiltonian matrix, scaled
    np.testing.assert_array_almost_equal(unitarymx[0:2,0:2], 0.5 * H_matrix)


def test_pauli_string_lcu_pyliqtr2():
    """Test converting a small, slightly more complex Hamiltonian to
       Pauli String LCU block encoding."""
    # Create a simple 2-qubit Hamiltonian
    pauli_data = {
        'XX': 1.0,
        'YZ': 1.0,
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    inst = PauliStringInstance(lcps)

    # Convert Hamiltonian to matrix
    H_matrix = inst.to_matrix(memory_threshold_gb=1.0)

    # Create unitary operator and convert to matrix
    unitaryop = PyLIQTRPauliStringLCU(inst, 'AS', probability_eps=0.5)
    unitarymx = unitaryop.tensor_contract()

    # Verify that the upper corner of our unitary is equal to the
    # original Hamiltonian matrix, scaled
    np.testing.assert_array_almost_equal(unitarymx[0:4,0:4], 0.5 * H_matrix)


