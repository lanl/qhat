from qhat.common.hermiticity import require_real_scalar
from qhat.common.pauli_string import PauliString
from qhat.common.pauli_utils import validate_pauli_string

class SinglePauliTerm:
    """A single Pauli string as a term in a Hamiltonian."""

    pauli_string: PauliString
    coefficient: float

    def __init__(self, pauli_string, coefficient):
        # Verify and store Pauli string
        # TODO: Should we allow any representation that PauliString understands?
        self.pauli_string = PauliString.from_dense(pauli_string)
        # Ensure coefficient is real (else Hamiltonian cannot be Hermitian) and force to a real
        # data type before storing
        self.coefficient = require_real_scalar(f"coefficient of {pauli_string}", coefficient)

    def num_qubits(self) -> int:
        return self.pauli_string.num_qubits

    def dense_string(self) -> str:
        return self.pauli_string.to_dense()

    def __str__(self) -> str:
        return f"{self.coefficient:10.3e} {self.dense_string()}"

    def __repr__(self) -> str:
        return f"SinglePauliTerm({self.dense_string()}, {self.coefficient})"
