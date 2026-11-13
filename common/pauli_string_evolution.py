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

from functools import cached_property
from typing import TYPE_CHECKING, Dict, Tuple

import attrs
import numpy as np
from scipy.linalg import expm

from numpy.typing import NDArray
from qualtran import Bloq, BloqBuilder, QBit, Register, Side, Signature, Soquet, SoquetT
from qualtran.bloqs.basic_gates import CNOT, GlobalPhase, Hadamard, Rz, SGate
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity

from qhat.common.hermiticity import require_real_scalar
from qhat.common.pauli_utils import validate_pauli_string

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
        """Validate the Pauli string, and store coefficient, time, and hbar as floats."""
        validate_pauli_string(self.pauli_string)
        for name in ("coefficient", "time", "hbar"):
            value = require_real_scalar(f"PauliStringEvolution {name}", getattr(self, name))
            object.__setattr__(self, name, value)

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

    @cached_property
    def _active_qubits(self) -> Tuple[int, ...]:
        """Indices of the qubits where the Pauli string is not the identity."""
        return tuple(i for i, p in enumerate(self.pauli_string) if p != 'I')

    @cached_property
    def _rotation_angle(self) -> float:
        """Angle theta of Rz(theta) = exp(-i theta Z / 2) that implements the parity rotation."""
        return 2.0 * self.coefficient * self.time / self.hbar

    def build_composite_bloq(self, bb: BloqBuilder, q: NDArray[Soquet]) -> Dict[str, SoquetT]:
        """Decompose exactly (global phase included) into basis changes, CNOTs and one Rz.

        Only the non-identity qubits are touched. Each is rotated so that its Pauli becomes Z
        (X: H; Y: S^dagger then H), a CNOT chain collects the parity on the last of them, Rz rotates
        by the parity, and the chain and basis changes are undone. An all-identity string is
        exp(-i * coefficient * time / hbar) times the identity, i.e. a `GlobalPhase`.
        """
        q = np.array(q)
        active = self._active_qubits
        if not active:
            bb.add(GlobalPhase(np.exp(-1j * self.coefficient * self.time / self.hbar)))
            return {'q': q}

        for i in active:
            if self.pauli_string[i] == 'Y':
                q[i] = bb.add(SGate(is_adjoint=True), q=q[i])
            if self.pauli_string[i] in 'XY':
                q[i] = bb.add(Hadamard(), q=q[i])
        for i, j in zip(active[:-1], active[1:]):
            q[i], q[j] = bb.add(CNOT(), ctrl=q[i], target=q[j])
        q[active[-1]] = bb.add(Rz(self._rotation_angle), q=q[active[-1]])
        for i, j in reversed(list(zip(active[:-1], active[1:]))):
            q[i], q[j] = bb.add(CNOT(), ctrl=q[i], target=q[j])
        for i in active:
            if self.pauli_string[i] in 'XY':
                q[i] = bb.add(Hadamard(), q=q[i])
            if self.pauli_string[i] == 'Y':
                q[i] = bb.add(SGate(), q=q[i])
        return {'q': q}

    def _t_complexity_(self) -> TComplexity:
        """T-complexity of the decomposition, counted by hand (no decomposition needed)."""
        active = self._active_qubits
        if not active:
            return TComplexity()
        num_basis = sum(2 if self.pauli_string[i] == 'X' else 4
                        for i in active if self.pauli_string[i] != 'Z')
        num_cnots = 2 * (len(active) - 1)
        return TComplexity(clifford=num_basis + num_cnots) + t_complexity(Rz(self._rotation_angle))

    def __str__(self) -> str:
        return f"exp(-i * {self.coefficient} * {self.pauli_string} * {self.time} / {self.hbar})"

    def __repr__(self) -> str:
        return (f"PauliStringEvolution(pauli_string='{self.pauli_string}', "
                f"coefficient={self.coefficient}, time={self.time}, hbar={self.hbar})")
