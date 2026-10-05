"""
Tests for FCIDump file loading functionality.
"""

import numpy as np
import pytest
from pathlib import Path

from openfermion import InteractionOperator

from qhat.analysis.hamiltonian import Hamiltonian, load_fci


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self):
        self.verbose = False


def test_toy_fci_file():
    """Test loading a toy FCIDump file (simple data for testing)."""
    # Path to the toy FCIDump file
    test_file = Path(__file__).parent / "fcidump.toy"

    if not test_file.exists():
        print("⚠ Skipping toy FCIDump file test (file not found)")
        return

    config_hamiltonian = MockConfig()
    config_hamiltonian.filename = str(test_file)
    config_hamiltonian.fermion_to_qubit_transform = None

    H = load_fci(config_hamiltonian)
    assert isinstance(H, Hamiltonian)
    intop = H.get_core_operator()
    assert isinstance(intop, InteractionOperator)

    # check scalars
    n = H.num_qubits()
    assert n == 3
    ecore = intop.constant
    assert ecore == 99999.0

    # check one-body coefficients
    h1 = intop.one_body_tensor
    assert h1.shape == (3,3)
    for i in range(3):
        for j in range(3):
            # put indices in canonical order, convert to decimal number
            ip = i + 1
            jp = j + 1
            d = (max(ip, jp), min(ip, jp))
            val = 10. * d[0] + d[1]
            assert val == h1[i, j]

    # check two-body coefficients
    h2 = intop.two_body_tensor
    assert h2.shape == (3,3,3,3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    # put indices in canonical order, convert to decimal number
                    ip = i + 1
                    jp = j + 1
                    kp = k + 1
                    lp = l + 1
                    ij = (max(ip, jp), min(ip, jp))
                    kl = (max(kp, lp), min(kp, lp))
                    d = np.empty(4)
                    d[[0,1]], d[[2,3]] = (max(ij, kl), min(ij, kl))
                    val = 1000. * d[0] + 100. * d[1] + 10. * d[2] + d[3]
                    assert val == h2[i, j, k, l]

    print("✓ Toy FCI file loaded successfully")


def test_real_fci_file():
    """Test loading a real FCIDump file (from PyLIQTR example)."""
    # Path to the real FCIDump file
    test_file = Path(__file__).parent / "fcidump.32_2ru_III_3pl"

    if not test_file.exists():
        print("⚠ Skipping real FCIDump file test (file not found)")
        return

    config_hamiltonian = MockConfig()
    config_hamiltonian.filename = str(test_file)
    config_hamiltonian.fermion_to_qubit_transform = None

    H = load_fci(config_hamiltonian)
    assert isinstance(H, Hamiltonian)
    intop = H.get_core_operator()
    assert isinstance(intop, InteractionOperator)

    # check scalars
    n = H.num_qubits()
    assert n == 7
    ecore = intop.constant
    assert ecore == -2885.478832556892

    # check one-body coefficients
    h1 = intop.one_body_tensor
    assert h1.shape == (7,7)
    # spot-check a few values
    # note:  this code is 0-based, but indices in file are 1-based!
    assert h1[0,0] == -12.71913914619708
    assert h1[0,1] == -0.0007824310522580448
    assert h1[6,6] == -3.202720690968762
    # check that all values are filled (nonzero)
    assert np.min(np.fabs(h1)) > 0.0

    # check two-body coefficients
    h2 = intop.two_body_tensor
    assert h2.shape == (7,7,7,7)
    # spot-check a few values
    # note:  this code is 0-based, but indices in file are 1-based!
    assert h2[0,0,0,0] == 0.6441135284888211
    assert h2[0,1,2,3] == 0.003581801198784521
    assert h2[6,6,6,6] == 0.2684082267201577
    # check that all values are filled (nonzero)
    assert np.min(np.fabs(h2)) > 0.0

    print("✓ Real FCI file loaded successfully")


if __name__ == "__main__":
    print("Running FCI loader tests...\n")

    test_toy_fci_file()
    test_real_fci_file()

    print("\n✓ All tests passed!")
