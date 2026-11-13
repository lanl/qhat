"""
Tests for FCIDump file loading functionality.
"""

import numpy as np
import pytest
from pathlib import Path

from openfermion import InteractionOperator

from qhat.analysis.hamiltonian import Hamiltonian, load_fci, load_numpy


class MockConfig:
    """Mock configuration object for testing."""
    def __init__(self):
        self.verbose = False


def test_toy_fci_file():
    """Test loading a toy FCIDump file (simple data for testing)."""
    # Path to the toy FCIDump file
    test_file = Path(__file__).parent / "toy.fci"
    assert test_file.exists()

    config_hamiltonian = MockConfig()
    config_hamiltonian.filename = str(test_file)
    config_hamiltonian.fermion_to_qubit_transform = None

    H = load_fci(config_hamiltonian)
    assert isinstance(H, Hamiltonian)
    intop = H.get_core_operator()
    assert isinstance(intop, InteractionOperator)

    # check scalars
    n = H.num_qubits()
    assert n == 6
    ecore = intop.constant
    assert ecore == 99999.0

    # check one-body coefficients
    h1 = intop.one_body_tensor
    assert h1.shape == (6,6)
    for i in range(3):
        for j in range(3):
            # put indices in canonical order, convert to decimal number
            ip = i + 1
            jp = j + 1
            d = (max(ip, jp), min(ip, jp))
            val = 10. * d[0] + d[1]
            assert val == h1[2*i  , 2*j  ]
            assert val == h1[2*i+1, 2*j+1]

    # check two-body coefficients
    h2 = intop.two_body_tensor
    assert h2.shape == (6,6,6,6)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    # put indices in canonical order, convert to decimal number
                    ip = i + 1
                    jp = j + 1
                    kp = k + 1
                    lp = l + 1
                    il = (max(ip, lp), min(ip, lp))
                    jk = (max(jp, kp), min(jp, kp))
                    d = np.empty(4)
                    d[[0,1]], d[[2,3]] = (max(il, jk), min(il, jk))
                    val = 1000. * d[0] + 100. * d[1] + 10. * d[2] + d[3]
                    assert val == h2[2*i  , 2*j  , 2*k  , 2*l  ]
                    assert val == h2[2*i+1, 2*j+1, 2*k+1, 2*l+1]
                    assert val == h2[2*i  , 2*j+1, 2*k+1, 2*l  ]
                    assert val == h2[2*i+1, 2*j  , 2*k  , 2*l+1]

    print("✓ Toy FCI file loaded successfully")


def check_32_2ru_III_3pl_operator(H, filetype):
#    """Check correctness of operator from file "32_2ru_III_3pl.*"
    assert isinstance(H, Hamiltonian)
    intop = H.get_core_operator()
    assert isinstance(intop, InteractionOperator)

    # check scalars
    n = H.num_qubits()
    assert n == 14
    ecore = intop.constant
    assert ecore == -2885.478832556892

    # check one-body coefficients
    h1 = intop.one_body_tensor
    assert h1.shape == (14,14)
    # spot-check a few values
    # note:  this code is 0-based, but indices in file are 1-based!
    assert h1[0,0] == -12.71913914619708
    assert h1[1,3] == -0.0007824310522580448
    assert h1[13,13] == -3.202720690968762
    # check that all values that should be filled are filled (nonzero)
    assert np.min(np.fabs(h1[0:n:2,0:n:2])) > 0.0
    assert np.min(np.fabs(h1[1:n:2,1:n:2])) > 0.0

    # check two-body coefficients
    h2 = intop.two_body_tensor
    assert h2.shape == (14,14,14,14)
    # spot-check a few values
    # note:  this code is 0-based, but indices in file are 1-based!
    assert h2[0,0,0,0] == 0.32205676424441054
    assert h2[1,2,4,7] == -0.006295966101239615
    assert h2[13,13,13,13] == 0.13420411336007884
    # check that all values that should be filled are filled (nonzero)
    # indices w/ same spin
    assert np.min(np.fabs(h2[0:n:2,0:n:2,0:n:2,0:n:2])) > 0.0
    assert np.min(np.fabs(h2[1:n:2,1:n:2,1:n:2,1:n:2])) > 0.0
    # indices w/ mixed spin
    assert np.min(np.fabs(h2[0:n:2,1:n:2,1:n:2,0:n:2])) > 0.0
    assert np.min(np.fabs(h2[1:n:2,0:n:2,0:n:2,1:n:2])) > 0.0

    print(f"✓ {filetype} loaded successfully")


def test_real_fci_file():
    """Test loading a real FCIDump file (from PyLIQTR example)."""
    # Path to the real FCIDump file
    test_file = Path(__file__).parent / "32_2ru_III_3pl.fci"
    assert test_file.exists()

    config_hamiltonian = MockConfig()
    config_hamiltonian.filename = str(test_file)
    config_hamiltonian.fermion_to_qubit_transform = None

    H = load_fci(config_hamiltonian)
    check_32_2ru_III_3pl_operator(H, "Real FCI file")


def test_numpy_file_openfermion():
    """Test loading a numpy file in openfermion convention."""
    # Path to the numpy file
    test_file = Path(__file__).parent.parent / "examples/Be-H_1.30_sto-6g_as-003-003.tensors.npz"
    assert test_file.exists()

    config_hamiltonian = MockConfig()
    config_hamiltonian.filename = str(test_file)
    config_hamiltonian.fermion_to_qubit_transform = None
    config_hamiltonian.integral_convention = "openfermion"

    H = load_numpy(config_hamiltonian)
    assert isinstance(H, Hamiltonian)
    intop = H.get_core_operator()
    assert isinstance(intop, InteractionOperator)

    # check scalars
    n = H.num_qubits()
    assert n == 6
    ecore = intop.constant
    assert ecore == 0.0

    # check one-body coefficients
    h1 = intop.one_body_tensor
    assert h1.shape == (6,6)
    # spot-check a few values
    assert h1[0,0] == -1.3055432829745055
    assert h1[1,3] == 0.03509463345405264
    assert h1[5,5] == -0.8518347869753975

    # check two-body coefficients
    h2 = intop.two_body_tensor
    assert h2.shape == (6,6,6,6)
    # spot-check a few values
    assert h2[0,0,0,0] == 0.2872954760337895
    assert h2[0,1,1,2] == -0.027478153443205564
    assert h2[5,5,5,5] == 0.2250582818237094

    print("✓ Numpy openfermion file loaded successfully")


def test_numpy_file_chemist_spatial():
    """Test loading a numpy file in chemist spatial convention."""
    # Path to the numpy file
    test_file = Path(__file__).parent / "32_2ru_III_3pl.npz"
    assert test_file.exists()

    config_hamiltonian = MockConfig()
    config_hamiltonian.filename = str(test_file)
    config_hamiltonian.fermion_to_qubit_transform = None
    config_hamiltonian.integral_convention = "chemist_spatial"

    H = load_numpy(config_hamiltonian)
    check_32_2ru_III_3pl_operator(H, "Numpy chemist spatial file")


if __name__ == "__main__":
    print("Running file loader tests...\n")

    test_toy_fci_file()
    test_real_fci_file()
    test_numpy_file_openfermion()
    test_numpy_file_chemist_spatial()

    print("\n✓ All tests passed!")
