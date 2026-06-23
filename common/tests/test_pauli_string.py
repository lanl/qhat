"""Tests for the representation-independent PauliString value type."""

import pytest

from qhat.common.pauli_string import PauliString


def test_construct_from_dense():
    pauli = PauliString.from_dense("IIXYIZ")

    assert pauli.num_qubits == 6
    assert len(pauli) == 6
    assert pauli.to_dense() == "IIXYIZ"
    assert pauli.to_sparse() == ((2, "X"), (3, "Y"), (5, "Z"))
    assert pauli.to_sparse_dict() == {2: "X", 3: "Y", 5: "Z"}


def test_dense_and_sparse_inputs_create_equal_values():
    dense = PauliString.from_dense("IIXYIZ")
    sparse = PauliString.from_sparse(
        {5: "Z", 2: "X", 3: "Y"},
        num_qubits=6,
    )

    assert dense == sparse
    assert hash(dense) == hash(sparse)


def test_pauli_string_can_be_a_dictionary_key():
    pauli = PauliString.from_dense("XII")
    hamiltonian = {pauli: 0.5}

    assert hamiltonian[PauliString.from_sparse(((0, "X"),), 3)] == 0.5


def test_string_representations():
    pauli = PauliString.from_dense("IIXYIZ")

    assert str(pauli) == "IIXYIZ"
    assert repr(pauli) == "PauliString.from_dense('IIXYIZ')"


@pytest.mark.parametrize(
    "dense, sparse",
    [
        ("IIII", ()),
        ("I", ()),
        ("XIII", ((0, "X"),)),
        ("IIZI", ((2, "Z"),)),
        ("XIZI", ((0, "X"), (2, "Z"))),
        ("XYZ", ((0, "X"), (1, "Y"), (2, "Z"))),
    ],
)
def test_dense_to_sparse_conversions(dense, sparse):
    assert PauliString.from_dense(dense).to_sparse() == sparse


@pytest.mark.parametrize(
    "sparse, num_qubits, dense",
    [
        ((), 4, "IIII"),
        ((), 1, "I"),
        (((0, "X"),), 4, "XIII"),
        (((2, "Z"),), 4, "IIZI"),
        (((3, "Y"),), 5, "IIIYI"),
        (((0, "X"), (2, "Z")), 4, "XIZI"),
        (((1, "Y"), (3, "Z")), 5, "IYIZI"),
    ],
)
def test_sparse_to_dense_conversions(sparse, num_qubits, dense):
    assert PauliString.from_sparse(sparse, num_qubits).to_dense() == dense


@pytest.mark.parametrize("dense", ["IIII", "XIII", "XYZI", "XYZIXYZ"])
def test_round_trip_conversion_preserves_dense_string(dense):
    pauli = PauliString.from_dense(dense)

    assert PauliString.from_sparse(pauli.to_sparse(), len(pauli)).to_dense() == dense


def test_sparse_input_is_sorted_and_drops_identity_operators():
    pauli = PauliString.from_sparse(
        ((3, "Z"), (0, "I"), (1, "X")),
        num_qubits=4,
    )

    assert pauli.to_sparse() == ((1, "X"), (3, "Z"))
    assert pauli.to_dense() == "IXIZ"


def test_identity_pauli_string():
    pauli = PauliString.from_sparse({}, num_qubits=4)

    assert pauli.to_dense() == "IIII"
    assert pauli.to_sparse() == ()


def test_zero_qubit_identity_pauli_string():
    dense = PauliString.from_dense("")
    sparse = PauliString.from_sparse((), num_qubits=0)

    assert dense == sparse
    assert dense.to_dense() == ""
    assert dense.to_sparse() == ()


@pytest.mark.parametrize("value", ["IXA", "xyz"])
def test_invalid_dense_input(value):
    with pytest.raises(ValueError):
        PauliString.from_dense(value)


def test_invalid_dense_input_identifies_invalid_operator():
    with pytest.raises(ValueError, match="Invalid Pauli operator 'Q'"):
        PauliString.from_dense("XQZI")


def test_num_qubits_cannot_be_negative():
    with pytest.raises(ValueError, match="non-negative"):
        PauliString.from_sparse((), num_qubits=-1)


@pytest.mark.parametrize(
    "operators, error, message",
    [
        (((-1, "X"),), ValueError, "outside the range"),
        (((3, "X"),), ValueError, "outside the range"),
        (((0, "X"), (0, "Y")), ValueError, "Duplicate"),
        (((0, "A"),), ValueError, "Invalid Pauli operator"),
        (((0.5, "X"),), TypeError, "indices must be integers"),
    ],
)
def test_invalid_sparse_input(operators, error, message):
    with pytest.raises(error, match=message):
        PauliString.from_sparse(operators, num_qubits=3)
