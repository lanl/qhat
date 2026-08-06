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
    unitaryop = PauliStringLCU(H, 'AS', probability_eps=0.25)
    unitarymx = unitaryop.tensor_contract()

    # Contract unitary tensor in the selection register
    # (highest-order bit)
    h = len(unitarymx) // 2
    unitarytc = 0.5 * (
            unitarymx[0:h,0:h] + unitarymx[0:h,h:] + 
            unitarymx[h:,0:h] + unitarymx[h:,h:]
            )

    # Verify that the upper corner of our unitary is equal to the
    # original Hamiltonian matrix, scaled
    np.testing.assert_array_almost_equal(unitarytc[0:4,0:4], 0.5 * H_matrix)


def test_pauli_string_lcu_harder():
    """Test converting a small, slightly more complex Hamiltonian to
       Pauli String LCU block encoding."""
    # Create a more complex 2-qubit Hamiltonian
    # Note:  we are cheating slightly with this.  After taking square
    # roots and normalizing, the relative probabilities for the three
    # terms should be (0.5, 0.5, sqrt(2)/2).  But since our selection
    # register has only two-bit precision, these probabilities will get
    # rounded to (0.25, 0.25, 0.5), which are proportional to the
    # original coefficients.
    pauli_data = {
        'II': 1.0,
        'XX': 1.0,
        'YZ': 2.0,
    }
    lcps = LinearCombinationOfPauliStrings(num_qubits=2, dense=pauli_data)
    H = Hamiltonian(lcps)

    # Convert Hamiltonian to matrix
    H_matrix = H.to_matrix(memory_threshold_gb=1.0)

    # Create unitary operator and convert to matrix
    unitaryop = PauliStringLCU(H, 'AS', probability_eps=0.25)
    unitarymx = unitaryop.tensor_contract()

    # Contract unitary tensor in the selection register
    # (highest-order 2 bits)
    q = len(unitarymx) // 4
    unitarytc = np.array(q * [q * [0.+0.j]])
    for r in range(4):
        for c in range(4):
            unitarytc[:,:] += unitarymx[r*q:(r+1)*q,c*q:(c+1)*q]
    unitarytc = 0.25 * unitarytc

    # Verify that the upper corner of our unitary is equal to the
    # original Hamiltonian matrix, scaled
    np.testing.assert_array_almost_equal(unitarytc[0:4,0:4], 0.25 * H_matrix)


