from abc import ABC, abstractmethod

class BaseTerm(ABC):
    """The abstract base class that defines the generic interface for Hamiltonian terms."""

    @abstractmethod
    def num_qubits(self) -> int:
        pass

    @abstractmethod
    def generate_evolution_operator(self, t: float, hbar: float = 1.0):
        """The bloq for exp(-i * scale * term * t / hbar)."""
        pass

    # TODO: I'm not sure if this sign_controlled method is the right approach for QHAT, but
    #       discussing the issue with QHAT suggests that in a number of cases we care about the
    #       sign-controlled approach may be less expensive than a controlled operation.  It's just
    #       a matter of deciding if this method is the right approach or if we should consider some
    #       other way to implement the idea of a sign-controlled variant.
    @abstractmethod
    def sign_controlled(self) -> BaseTerm:
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
    def __str__(self) -> str:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

