"""
Test construction of Pauli String LCU block encodings.
"""

import numpy as np
import pytest

from qhat.analysis.config_types import (
    GeneralConfiguration,
    GeneralConfigurationUser,
)
from qhat.analysis.hamiltonian import (
    Hamiltonian,
    LinearCombinationOfPauliStrings,
)
from qhat.analysis.unitary import PauliStringLCU

def test_pauli_string_lcu_simple():
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


def test_pauli_string_lcu_harder():
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


