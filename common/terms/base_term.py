from abc import ABC, abstractmethod
from typing import Tuple

from qhat.common import controls
from qhat.common.pauli_string import PauliString

class BaseTerm(ABC):
    """The abstract base class that defines the generic interface for Hamiltonian terms.

    Terms must be immutable value objects with `__eq__` and `__hash__`, because bloqs such as
    `Trotterization` hold them and need to compare and hash equal when built from equal terms.
    The simplest way to meet this is to make each term class `@attrs.frozen` (see
    `SinglePauliTerm`, which keeps a custom `__init__` with `init=False` and `__attrs_init__`);
    its fields must themselves be immutable and hashable (e.g. tuples, not lists).

    Adding a new term type:
    -- Define the class in a module under `common/terms/`, which is where the tests look for
       term classes.
    -- Register it in `TERM_EXAMPLES` in `common/tests/test_terms.py`: add an entry that maps the
       class to a function returning a list of at least two different terms (the function is
       called several times, and separately built terms must compare equal).  That runs the
       generic tests in `TestTermContract` on your class, which check that:
          -- terms built separately from the same values are equal and hash the same, different
             terms are not equal, and terms work as set members and dict keys;
          -- attributes cannot be reassigned or added;
          -- `num_qubits`, `to_pauli_strings`, `str` and `repr` return sensible values;
          -- `sign_controlled` returns a term with one more qubit;
          -- `generate_evolution_operator` returns a bloq with a `q` register of `num_qubits()`
             qubits, a `control_phase`, a unitary matrix, an adjoint that is the conjugate
             transpose, and a `ctrl` register when controlled.
       `test_every_term_class_is_registered` fails if a term class under `common/terms/` has no
       entry.
    -- Put tests for behavior specific to your class in its own class in that file, like
       `TestSinglePauliTerm`.
    """

    @abstractmethod
    def num_qubits(self) -> int:
        pass

    # TODO: What should the sparse / matrix-free version look like?
    #       -- to_sparse_operator: create a new thing that you can do
    #          > sparse_op = term.to_sparse_operator()
    #          > result = sparse_op @ vector
    #       -- __matmul__: the term itself can do matrix-free operations
    #          > result = term @ vector
    #       This is a future design decision, because we can do a lot already with a dense matrix,
    #       but that doesn't scale well so matrix-free operations will be needed to scale certain
    #       analyses.  So for now we explicitly build a dense matrix (for the "simple" but
    #       poorly-scaling versions of these analyses), and later on we'll address matrix-free
    #       operations.
    @abstractmethod
    def to_dense_matrix(self):
        """
        Generates a dense matrix that implements this term.
        """
        pass

    @abstractmethod
    def generate_evolution_operator(
        self,
        t: float,
        hbar: float = 1.0,
        *,
        is_adjoint: bool = False,
        ctrl_values: Tuple[int, ...] = (),
        multi_control_method: str = controls.DEFAULT_MULTI_CONTROL_METHOD,
        controlled_rotation_method: str = controls.DEFAULT_CONTROLLED_ROTATION_METHOD,
        omit_control_phase: bool = False,
    ):
        """The bloq for exp(-i * term * t / hbar), or its adjoint, under the given controls.

        The bloq has a `q` register of `num_qubits()` qubits and, if `ctrl_values` is not empty, a
        `ctrl` register (see `common/controls.py`). `Trotterization` also needs the bloq to have
        a `control_phase` attribute, and to support `tensor_contract()`.
        """
        pass

    # TODO: I'm not sure if this sign_controlled method is the right approach for QHAT, but
    #       discussing the issue with QHAT suggests that in a number of cases we care about the
    #       sign-controlled approach may be less expensive than a controlled operation.  It's just
    #       a matter of deciding if this method is the right approach or if we should consider some
    #       other way to implement the idea of a sign-controlled variant.
    @abstractmethod
    def sign_controlled(self):
        """
        Return the sign-controlled version of this Hermitian term T.

        The returned term is
            R_s = -Z_c ⊗ R,
        where Z_c is the Pauli Z operator on qubit c and R is an operator discussed below.  Let
            V = exp(-i * R * t / hbar).
        Then evolution under R_s is
            exp(-i R_s t / hbar) = |0><0| ⊗ exp(+i R t / hbar) + |1><1| ⊗ exp(-i R t / hbar)
                                 = |0><0| ⊗ V† + |1><1| ⊗ V
                                 = S[V].
        In applications such as QPE, S[V] can replace a controlled unitary C[U] when V**2 = U,
        because
            S[V] = (I ⊗ V†) C[U],
        where C[U] is controlled on qubit c.  Thus S[V] and C[U] have the same relative evolution
        between the control-|0> and control-|1> branches.  They give the same phase-register
        measurement distribution in QPE, although S[V] also applies the target-only operation V†.

        This avoids directly implementing C[U]; instead, it evolves under the enlarged generator
        -Z_c ⊗ R.  Whether this is cheaper depends on whether that enlarged generator can be
        simulated efficiently.

        For U = exp(-i T t / hbar), an exact continuous-time choice is R = T / 2.  For an
        approximated circuit (such as a product formula), construct V so that its implemented
        circuit satisfies V**2 = U.
        """
        pass

    @abstractmethod
    def to_pauli_strings(self) -> dict[PauliString, float]:
        # TODO: This is a temporary bridge method to be removed once Pauli-based analyses
        #       (eigenvalue bounds, nested-commutator bounds, etc) are updated to use the term
        #       interface.
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

