"""
Tests for the Hamiltonian term types in `common/terms/`.

`TestTermContract` runs against every term type registered in `TERM_EXAMPLES`.  A new term type
should add an entry there (see the `BaseTerm` docstring); `test_every_term_class_is_registered`
fails if it does not.  `TestSinglePauliTerm` holds the tests that are specific to `SinglePauliTerm`.
"""

import importlib
import inspect
import numbers
import pkgutil

import attrs
import numpy as np
import pytest
from scipy.linalg import expm

from qhat.common import terms as terms_package
from qhat.common.pauli_string import PauliString
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import pauli_string_to_matrix
from qhat.common.terms.base_term import BaseTerm
from qhat.common.terms.single_pauli_term import SinglePauliTerm


# Maps each term class to a function that returns a list of at least two terms that differ from
# each other.  The function is called several times, and separately built terms must compare equal.
TERM_EXAMPLES = {
    SinglePauliTerm: lambda: [
        SinglePauliTerm("XZI", 0.5),
        SinglePauliTerm("XZI", -0.5),
        SinglePauliTerm("IYY", 0.5),
        SinglePauliTerm("III", 0.25),
    ],
}


def _attribute_names(term):
    if attrs.has(type(term)):
        return [field.name for field in attrs.fields(type(term))]
    return list(vars(term))


def _term_class_keys():
    """(module, qualified name) of every concrete term class defined under `qhat.common.terms`.

    Keyed by name rather than by class, because attrs may leave an unslotted copy of a class
    registered as a subclass until it is garbage collected.
    """
    for module_info in pkgutil.iter_modules(terms_package.__path__):
        importlib.import_module(f"{terms_package.__name__}.{module_info.name}")
    keys, pending = set(), [BaseTerm]
    while pending:
        for subclass in pending.pop().__subclasses__():
            pending.append(subclass)
            if (subclass.__module__.startswith(terms_package.__name__ + ".")
                    and not inspect.isabstract(subclass)):
                keys.add((subclass.__module__, subclass.__qualname__))
    return keys


def test_every_term_class_is_registered():
    registered = {(cls.__module__, cls.__qualname__) for cls in TERM_EXAMPLES}
    missing = _term_class_keys() - registered
    assert not missing, f"Add examples to TERM_EXAMPLES in {__file__} for: {sorted(missing)}"


class TestTermContract:
    """What every `BaseTerm` must provide, checked for each class in `TERM_EXAMPLES`."""

    @pytest.fixture(params=list(TERM_EXAMPLES), ids=lambda cls: cls.__name__)
    def make_examples(self, request):
        return TERM_EXAMPLES[request.param]

    def test_examples_are_distinct_terms(self, make_examples):
        examples = make_examples()
        assert len(examples) >= 2
        assert all(isinstance(term, BaseTerm) for term in examples)

    def test_equal_when_built_separately(self, make_examples):
        for term, same in zip(make_examples(), make_examples()):
            assert term is not same
            assert term == same
            assert hash(term) == hash(same)

    def test_different_terms_are_not_equal(self, make_examples):
        examples = make_examples()
        for i, term in enumerate(examples):
            for other in examples[i + 1:]:
                assert term != other

    def test_not_equal_to_other_types(self, make_examples):
        for term in make_examples():
            assert term != object()
            assert term != str(term)
            assert term != None  # noqa: E711

    def test_usable_in_sets_and_dicts(self, make_examples):
        examples = make_examples()
        assert len(set(examples + make_examples())) == len(examples)
        lookup = {term: i for i, term in enumerate(examples)}
        assert [lookup[term] for term in make_examples()] == list(range(len(examples)))

    def test_attributes_cannot_be_reassigned(self, make_examples):
        for term in make_examples():
            names = _attribute_names(term)
            assert names
            for name in names:
                with pytest.raises(AttributeError):
                    setattr(term, name, getattr(term, name))
            with pytest.raises(AttributeError):
                term.some_new_attribute = 1

    def test_interface(self, make_examples):
        for term in make_examples():
            num_qubits = term.num_qubits()
            assert isinstance(num_qubits, int) and num_qubits >= 1
            pauli_strings = term.to_pauli_strings()
            assert pauli_strings
            for pauli, coefficient in pauli_strings.items():
                assert isinstance(pauli, PauliString) and pauli.num_qubits == num_qubits
                assert isinstance(coefficient, numbers.Real)
            assert isinstance(str(term), str) and str(term)
            assert isinstance(repr(term), str) and repr(term)

    def test_sign_controlled_adds_a_qubit(self, make_examples):
        for term in make_examples():
            signed = term.sign_controlled()
            assert isinstance(signed, BaseTerm)
            assert signed.num_qubits() == term.num_qubits() + 1

    def test_evolution_operator(self, make_examples):
        for term in make_examples():
            n = term.num_qubits()
            op = term.generate_evolution_operator(0.3, 1.1)
            registers = {register.name: register for register in op.signature}
            assert registers["q"].total_bits() == n
            assert "ctrl" not in registers
            assert hasattr(op, "control_phase")
            unitary = op.tensor_contract()
            assert unitary.shape == (2**n, 2**n)
            np.testing.assert_allclose(unitary @ unitary.conj().T, np.eye(2**n), atol=1e-10)

    def test_adjoint_evolution_operator(self, make_examples):
        for term in make_examples():
            unitary = term.generate_evolution_operator(0.3, 1.1).tensor_contract()
            adjoint = term.generate_evolution_operator(0.3, 1.1, is_adjoint=True).tensor_contract()
            np.testing.assert_allclose(adjoint, unitary.conj().T, atol=1e-10)

    def test_controlled_evolution_operator(self, make_examples):
        for term in make_examples():
            n = term.num_qubits()
            op = term.generate_evolution_operator(0.3, 1.1, ctrl_values=(1,))
            registers = {register.name: register for register in op.signature}
            assert registers["q"].total_bits() == n
            assert registers["ctrl"].total_bits() == 1
            assert op.tensor_contract().shape == (2 ** (n + 1), 2 ** (n + 1))


class TestSinglePauliTerm:
    def test_input_formats_give_equal_terms(self):
        dense = SinglePauliTerm("XZI", 0.5)
        from_object = SinglePauliTerm(PauliString.from_dense("XZI"), 0.5)
        sparse = SinglePauliTerm(((0, "X"), (1, "Z")), 0.5, num_qubits=3)
        reordered = SinglePauliTerm(((1, "Z"), (0, "X")), 0.5, num_qubits=3)
        assert dense == from_object == sparse == reordered
        assert len({hash(term) for term in (dense, from_object, sparse, reordered)}) == 1

    def test_sparse_requires_num_qubits(self):
        with pytest.raises(ValueError, match="num_qubits is required"):
            SinglePauliTerm(((0, "X"),), 0.5)

    def test_identity_term(self):
        term = SinglePauliTerm((), 2.0, num_qubits=3)
        assert term.dense_string() == "III"
        assert term.num_qubits() == 3

    @pytest.mark.parametrize("raw, expected", [
        (2, 2.0),
        (np.float64(0.37), 0.37),
        (np.complex128(0.37 + 1e-19j), 0.37),
        (np.array(0.5), 0.5),
    ])
    def test_coefficient_is_stored_as_float(self, raw, expected):
        term = SinglePauliTerm("XZ", raw)
        assert type(term.coefficient) is float
        assert term.coefficient == expected

    @pytest.mark.parametrize("raw, error, match", [
        (0.25j, ValueError, "must be real"),
        (float("nan"), ValueError, "must be finite"),
        (float("inf"), ValueError, "must be finite"),
        (True, TypeError, "real number"),
        ("0.5", TypeError, "real number"),
    ])
    def test_invalid_coefficient_is_rejected(self, raw, error, match):
        with pytest.raises(error, match=match):
            SinglePauliTerm("XZ", raw)

    def test_equality_depends_on_pauli_string_and_coefficient(self):
        term = SinglePauliTerm("XZ", 0.5)
        assert term == SinglePauliTerm("XZ", 0.5)
        assert term != SinglePauliTerm("XZ", 0.25)
        assert term != SinglePauliTerm("ZX", 0.5)
        assert term != ("XZ", 0.5)

    def test_accessors(self):
        term = SinglePauliTerm("XIZ", -0.25)
        assert term.num_qubits() == 3
        assert term.dense_string() == "XIZ"
        assert term.to_pauli_strings() == {PauliString.from_dense("XIZ"): -0.25}
        assert repr(term) == "SinglePauliTerm(XIZ, -0.25)"
        assert "XIZ" in str(term)

    def test_sign_controlled(self):
        assert SinglePauliTerm("XIZ", -0.25).sign_controlled() == SinglePauliTerm("ZXIZ", 0.25)

    def test_evolution_operator_arguments(self):
        term = SinglePauliTerm("XIZ", -0.25)
        op = term.generate_evolution_operator(
            0.3, 1.1, is_adjoint=True, ctrl_values=(1,), controlled_rotation_method="ancilla")
        assert op == PauliStringEvolution(
            "XIZ", -0.25, 0.3, 1.1, is_adjoint=True, ctrl_values=(1,),
            controlled_rotation_method="ancilla")

    def test_evolution_operator_matrix(self):
        term = SinglePauliTerm("XZ", -0.25)
        unitary = term.generate_evolution_operator(0.3, 1.1).tensor_contract()
        expected = expm(-1j * -0.25 * 0.3 / 1.1 * pauli_string_to_matrix("XZ"))
        np.testing.assert_allclose(unitary, expected, atol=1e-12)
