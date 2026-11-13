from typing import Optional

from qhat.common.hermiticity import require_real_scalar
from qhat.common.pauli_string import PauliString
# TODO: Should I move the evolution Bloqs into the terms/ directory?
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.terms.base_term import BaseTerm

class SinglePauliTerm(BaseTerm):
    """A single Pauli string as a term in a Hamiltonian."""

    pauli_string: PauliString
    coefficient: float

    def __init__(self, pauli_string, coefficient, num_qubits: Optional[int] = None):
        """`pauli_string` is a PauliString, a dense string, or sparse operators (which need
        `num_qubits`)."""
        # Verify and store Pauli string
        if isinstance(pauli_string, PauliString):
            self.pauli_string = pauli_string
        elif isinstance(pauli_string, str):
            self.pauli_string = PauliString.from_dense(pauli_string)
        else:
            if num_qubits is None:
                raise ValueError("num_qubits is required for a sparse Pauli string.")
            # PauliString.from_sparse will flag any invalid formats for us
            self.pauli_string = PauliString.from_sparse(pauli_string, num_qubits)
        # Ensure coefficient is real (else Hamiltonian cannot be Hermitian) and force to a real
        # data type before storing
        self.coefficient = require_real_scalar(f"coefficient of {self.pauli_string}", coefficient)

    def num_qubits(self) -> int:
        return self.pauli_string.num_qubits

    def generate_evolution_operator(self, t: float, hbar: float = 1.0, *, scale: float = 1.0,
                                    **evolution_options):
        return PauliStringEvolution(self.dense_string(), self.coefficient * scale, t, hbar,
                                    **evolution_options)

    def dense_string(self) -> str:
        return self.pauli_string.to_dense()

    def has_rotations(self) -> bool:
        return set(self.dense_string()) != {'I'}

    def sign_controlled(self) -> 'SinglePauliTerm':
        return SinglePauliTerm('Z' + self.dense_string(), -self.coefficient)

    # Value semantics: Trotterization is an attrs bloq, so equal terms must compare and hash equal.
    def __eq__(self, other) -> bool:
        if not isinstance(other, SinglePauliTerm):
            return NotImplemented
        return self.pauli_string == other.pauli_string and self.coefficient == other.coefficient

    def __hash__(self) -> int:
        return hash((self.pauli_string, self.coefficient))

    def __str__(self) -> str:
        return f"{self.coefficient:10.3e} {self.dense_string()}"

    def __repr__(self) -> str:
        return f"SinglePauliTerm({self.dense_string()}, {self.coefficient})"

