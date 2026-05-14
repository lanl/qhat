"""
Pauli String Time Evolution Bloq

This module implements a Qualtran Bloq for time evolution under a single Pauli string Hamiltonian:
    U = exp(-i * coefficient * PauliString * t / ℏ)

where:
    - PauliString is a tensor product of Pauli operators (I, X, Y, Z)
    - coefficient is a scalar weight
    - t is the evolution time
    - ℏ is the reduced Planck constant (defaults to 1 in natural units)
"""

import attrs
import numpy as np
from functools import cached_property
from typing import Dict, TYPE_CHECKING
from scipy.linalg import expm

from cirq import LineQubit, PauliStringPhasorGate, DensePauliString as CirqDensePauliString
from qualtran import Bloq, BloqBuilder, Signature, SoquetT, QBit, Register, Side
from qualtran.cirq_interop import CirqGateAsBloq
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity
from common.pauli_utils import validate_pauli_string

if TYPE_CHECKING:
    import quimb.tensor as qtn


@attrs.frozen
class PauliStringEvolution(Bloq):
    """Time evolution operator for a single Pauli string Hamiltonian.

    Implements the unitary operator:
        U = exp(-i * coefficient * PauliString * time / hbar)

    Attributes:
        pauli_string: String of Pauli operators, e.g., "IXXYIIZ". Each character must be
            one of 'I', 'X', 'Y', or 'Z'.
        coefficient: Scalar coefficient multiplying the Pauli string in the Hamiltonian.
        time: Evolution time parameter.
        hbar: Reduced Planck constant (default 1.0 for natural units).

    Example:
        >>> pse = PauliStringEvolution(pauli_string="XYZ", coefficient=0.5, time=1.0)
        >>> U = pse.tensor_contract()  # Get the unitary matrix
        >>> from qualtran.cirq_interop.t_complexity_protocol import t_complexity
        >>> tc = t_complexity(pse)  # Get resource estimates
    """

    pauli_string: str
    coefficient: float
    time: float
    hbar: float = 1.0

    def __attrs_post_init__(self):
        """Validate the Pauli string contains only valid characters."""
        validate_pauli_string(self.pauli_string)

    @cached_property
    def signature(self) -> Signature:
        """Return the signature with a register of qubits matching the Pauli string length."""
        # Use shape=(n,) to create an array of n individual qubits, not a single n-bit register
        return Signature([Register('q', QBit(), shape=(self.num_qubits,), side=Side.THRU)])

    @cached_property
    def num_qubits(self) -> int:
        """Number of qubits the operator acts on."""
        return len(self.pauli_string)

    def _get_pauli_matrix(self) -> np.ndarray:
        """Get the matrix representation of the Pauli string."""
        # Build up the full Pauli matrix via tensor products
        pauli_matrices = {
            'I': np.array([[1, 0], [0, 1]], dtype=complex),
            'X': np.array([[0, 1], [1, 0]], dtype=complex),
            'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
            'Z': np.array([[1, 0], [0, -1]], dtype=complex),
        }

        result = pauli_matrices[self.pauli_string[0]]
        for char in self.pauli_string[1:]:
            result = np.kron(result, pauli_matrices[char])

        return result

    def _compute_unitary_matrix(self) -> np.ndarray:
        """Compute the unitary matrix for this Pauli string evolution."""
        # Build the Pauli string matrix
        P = self._get_pauli_matrix()

        # Compute the evolution operator: exp(-i * coefficient * P * time / hbar)
        angle = -1j * self.coefficient * self.time / self.hbar
        U = expm(angle * P)

        return U

    def add_my_tensors(
        self,
        tn: 'qtn.TensorNetwork',
        tag: object,
        *,
        incoming: Dict[str, 'SoquetT'],
        outgoing: Dict[str, 'SoquetT'],
    ):
        """Add tensors representing this bloq to the tensor network."""
        import quimb.tensor as qtn

        # Get the unitary matrix for this operation
        unitary = self._compute_unitary_matrix()

        # With shape=(n,), incoming['q'] is an array of Soquets
        import_idx = [incoming['q'][i] for i in range(self.num_qubits)]
        export_idx = [outgoing['q'][i] for i in range(self.num_qubits)]

        # Add the unitary as a tensor to the network
        tn.add(
            qtn.Tensor(
                data=unitary.reshape([2] * (2 * self.num_qubits)),
                inds=tuple(export_idx + import_idx),
                tags=[self.pretty_name(), tag],
            )
        )

    def _build_cirq_gate(self):
        """Build the corresponding Cirq gate for resource estimation."""
        # Create a Cirq DensePauliString directly from the Pauli string
        # We use coefficient=1.0 as the base, and scale via the power parameter below
        cirq_pauli = CirqDensePauliString(self.pauli_string, coefficient=1.0)

        # PauliStringPhasorGate implements exp(-i * exponent * P / 2) where P has eigenvalues ±1
        # We want exp(-i * coefficient * time / hbar * P)
        # So we need: exponent / 2 = coefficient * time / hbar
        # Therefore: exponent = 2 * coefficient * time / hbar
        # And power = exponent / π (since PauliString**power gives exponent = π * power)
        power = 2.0 * self.coefficient * self.time / (self.hbar * np.pi)

        # Create the PauliStringPhasorGate from the base Pauli string
        exp_gate = PauliStringPhasorGate(cirq_pauli) ** power

        return exp_gate

    def build_composite_bloq(self, bb: BloqBuilder, **soqs: SoquetT) -> Dict[str, SoquetT]:
        """Decompose into a Cirq gate for resource estimation.

        This decomposition is needed for pyLIQTR and controlled versions to work properly.
        """
        # Get the Cirq gate representation
        cirq_gate = self._build_cirq_gate()

        # Wrap it as a Bloq
        gate_bloq = CirqGateAsBloq(cirq_gate)

        # Add it to the bloq builder
        return bb.add_d(gate_bloq, **soqs)

    def _t_complexity_(self) -> TComplexity:
        """Return the T-complexity for resource estimation."""
        # Build the internal Cirq gate representation and get its complexity
        cirq_gate = self._build_cirq_gate()
        bloq = CirqGateAsBloq(cirq_gate)
        return t_complexity(bloq)

    def __str__(self) -> str:
        return f"exp(-i * {self.coefficient} * {self.pauli_string} * {self.time} / {self.hbar})"

    def __repr__(self) -> str:
        return (f"PauliStringEvolution(pauli_string='{self.pauli_string}', "
                f"coefficient={self.coefficient}, time={self.time}, hbar={self.hbar})")
