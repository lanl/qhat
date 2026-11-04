import numpy as np
from scipy.linalg import expm

from trotter_isomorphism.pauli_keys import pauli_key_to_matrix
from trotter_isomorphism.product_formula import (
    dense_trotter_unitary,
    exact_unitary,
    first_order_trotter_unitary,
    second_order_trotter_unitary,
)


def test_first_order_wrapper_matches_dense_trotter_unitary():
    terms = [(((0, "X"),), 0.25), (((0, "Z"),), -0.5)]

    wrapped = first_order_trotter_unitary(terms, t=0.7, r=3, n_qubits=1)
    direct = dense_trotter_unitary(terms, time=0.7, num_steps=3, n_qubits=1, method=1)

    assert np.allclose(wrapped, direct)


def test_second_order_single_step_has_symmetric_order():
    terms = [(((0, "X"),), 0.25), (((0, "Z"),), -0.5)]
    dt = 0.7

    X = pauli_key_to_matrix(((0, "X"),), 1)
    Z = pauli_key_to_matrix(((0, "Z"),), 1)
    expected = (
        expm(-1j * dt * 0.5 * 0.25 * X)
        @ expm(-1j * dt * 1.0 * -0.5 * Z)
        @ expm(-1j * dt * 0.5 * 0.25 * X)
    )

    actual = second_order_trotter_unitary(terms, t=dt, r=1, n_qubits=1)

    assert np.allclose(actual, expected)


def test_empty_term_list_returns_identity():
    actual = dense_trotter_unitary([], time=1.0, num_steps=4, n_qubits=2, method=2)

    assert np.allclose(actual, np.eye(4))


def test_single_term_trotter_matches_exact_unitary():
    terms = [(((0, "Y"),), 0.3)]
    Y = pauli_key_to_matrix(((0, "Y"),), 1)

    actual = dense_trotter_unitary(terms, time=0.6, num_steps=5, n_qubits=1, method=4)
    expected = exact_unitary(0.3 * Y, time=0.6)

    assert np.allclose(actual, expected)
