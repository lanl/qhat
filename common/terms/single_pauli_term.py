from typing import Optional, Tuple

from qhat.common import controls
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
        return PauliStringEvolution(
            self.dense_string(), self.coefficient, t, hbar,
            is_adjoint=is_adjoint,
            ctrl_values=ctrl_values,
            multi_control_method=multi_control_method,
            controlled_rotation_method=controlled_rotation_method,
            omit_control_phase=omit_control_phase,
        )

    def sign_controlled(self):
        return SinglePauliTerm('Z' + self.dense_string(), -self.coefficient)

    def to_pauli_strings(self) -> dict[PauliString, float]:
        return {self.pauli_string : self.coefficient}

    def dense_string(self) -> str:
        return self.pauli_string.to_dense()

    def __str__(self) -> str:
        return f"{self.coefficient:10.3e} {self.dense_string()}"

    def __repr__(self) -> str:
        return f"SinglePauliTerm({self.dense_string()}, {self.coefficient})"

