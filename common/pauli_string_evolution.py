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

import numbers
from collections import Counter
from functools import cached_property
from typing import TYPE_CHECKING, Dict, Optional, Set, Tuple

import attrs
import numpy as np
from scipy.linalg import expm

from numpy.typing import NDArray
from qualtran import Bloq, BloqBuilder, CtrlSpec, QBit, Register, Side, Signature, Soquet, SoquetT
from qualtran.bloqs.basic_gates import CNOT, GlobalPhase, Hadamard, Rz, SGate
from qualtran.cirq_interop.t_complexity_protocol import TComplexity
from qualtran.resource_counting import QubitCount

from qhat.common import controls
from qhat.common.hermiticity import require_real_scalar
from qhat.common.pauli_utils import validate_pauli_string

if TYPE_CHECKING:
    import quimb.tensor as qtn
    from qualtran import AddControlledT
    from qualtran.resource_counting import BloqCountT, CostKey, SympySymbolAllocator


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
        is_adjoint: If True, the operator is U^dagger = exp(+i * coefficient * PauliString * time
            / hbar). `adjoint()` flips it; `time` is left as given (a negative `time` is a
            different object that happens to be analytically equal).
        ctrl_values: One entry per control qubit, 1 (act on |1>) or 0 (act on |0>); empty for
            no control. `.controlled()` prepends one, e.g. `.controlled(CtrlSpec(cvs=0))` gives
            `(0,)`. The controls are the `ctrl` register, before `q`.
        multi_control_method: How 2 or more controls are built, "ancilla" (default) or "direct";
            see `common/controls.py`. Ignored with 0 or 1 controls.

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
    is_adjoint: bool = False
    ctrl_values: Tuple[int, ...] = ()
    multi_control_method: str = controls.DEFAULT_MULTI_CONTROL_METHOD

    def __attrs_post_init__(self):
        """Validate the Pauli string and controls, and store coefficient, time, and hbar as floats."""
        validate_pauli_string(self.pauli_string)
        for name in ("coefficient", "time", "hbar"):
            value = require_real_scalar(f"PauliStringEvolution {name}", getattr(self, name))
            object.__setattr__(self, name, value)
        object.__setattr__(self, "ctrl_values",
                           controls.validate_ctrl_values("PauliStringEvolution", self.ctrl_values))
        controls.validate_multi_control_method("PauliStringEvolution", self.multi_control_method)

    @cached_property
    def signature(self) -> Signature:
        """The `ctrl` register (if controlled), then one qubit per Pauli string position."""
        # Use shape=(n,) to create an array of n individual qubits, not a single n-bit register
        return Signature(controls.ctrl_registers(self.ctrl_values)
                         + [Register('q', QBit(), shape=(self.num_qubits,), side=Side.THRU)])

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

        # Compute the evolution operator: exp(-i * coefficient * P * time / hbar), conjugate
        # transposed for the adjoint
        angle = -1j * self.coefficient * self.time / self.hbar
        U = expm(angle * P)

        U = U.conj().T if self.is_adjoint else U
        return controls.controlled_matrix(U, self.ctrl_values)

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

        # With shape=(n,), incoming['q'] is an array of Soquets; controls come first
        num_ctrls = len(self.ctrl_values)
        import_idx = ([incoming['ctrl'][i] for i in range(num_ctrls)]
                      + [incoming['q'][i] for i in range(self.num_qubits)])
        export_idx = ([outgoing['ctrl'][i] for i in range(num_ctrls)]
                      + [outgoing['q'][i] for i in range(self.num_qubits)])

        # Add the unitary as a tensor to the network
        tn.add(
            qtn.Tensor(
                data=unitary.reshape([2] * (2 * (num_ctrls + self.num_qubits))),
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
        """Angle theta of Rz(theta) = exp(-i theta Z / 2) that implements the parity rotation.

        The adjoint negates it; the basis changes and the CNOT chain are their own inverse around
        the rotation, so they are unchanged.
        """
        sign = -1.0 if self.is_adjoint else 1.0
        return sign * 2.0 * self.coefficient * self.time / self.hbar

    def build_composite_bloq(
        self, bb: BloqBuilder, q: NDArray[Soquet], ctrl: Optional[NDArray[Soquet]] = None
    ) -> Dict[str, SoquetT]:
        """Decompose exactly (global phase included) into basis changes, CNOTs and Rz rotations.

        Only the non-identity qubits are touched. Each is rotated so that its Pauli becomes Z
        (X: H; Y: S^dagger then H), a CNOT chain collects the parity on the last of them, Rz rotates
        by the parity, and the chain and basis changes are undone. An all-identity string is
        exp(-i * coefficient * time / hbar) times the identity, i.e. a `GlobalPhase`.

        Under control only the rotation is controlled (see `_add_rotation`); the basis changes and
        the CNOT chain cancel when the rotation is skipped. With 2 or more controls and the
        "ancilla" method, the controls are first ANDed into one ancilla, which controls this bloq.
        """
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.add_with_and_ancilla(bb, self, ctrl=ctrl, q=q)
        q = np.array(q)
        active = self._active_qubits
        if not active:
            ctrl = self._add_identity_phase(bb, ctrl)
            return self._outputs(q, ctrl)

        for i in active:
            if self.pauli_string[i] == 'Y':
                q[i] = bb.add(SGate(is_adjoint=True), q=q[i])
            if self.pauli_string[i] in 'XY':
                q[i] = bb.add(Hadamard(), q=q[i])
        for i, j in zip(active[:-1], active[1:]):
            q[i], q[j] = bb.add(CNOT(), ctrl=q[i], target=q[j])
        q[active[-1]], ctrl = self._add_rotation(bb, q[active[-1]], ctrl)
        for i, j in reversed(list(zip(active[:-1], active[1:]))):
            q[i], q[j] = bb.add(CNOT(), ctrl=q[i], target=q[j])
        for i in active:
            if self.pauli_string[i] in 'XY':
                q[i] = bb.add(Hadamard(), q=q[i])
            if self.pauli_string[i] == 'Y':
                q[i] = bb.add(SGate(), q=q[i])
        return self._outputs(q, ctrl)

    def _outputs(self, q: NDArray[Soquet], ctrl: Optional[NDArray[Soquet]]) -> Dict[str, SoquetT]:
        return {'q': q} if ctrl is None else {'ctrl': ctrl, 'q': q}

    def _add_rotation(self, bb: BloqBuilder, target: Soquet, ctrl: Optional[NDArray[Soquet]]):
        """Add Rz(theta) on the parity qubit `target`, controlled on `ctrl` if there is one.

        With one control the rotation is Rz(theta/2) . CNOT . Rz(-+theta/2) . CNOT on the target (2
        rotations, 2 CNOTs, exact): if the control is |1> the CNOTs flip the second rotation's
        sign, so the two add up to Rz(theta); if |0> they cancel. For a control on |0> the second
        angle is +theta/2 instead, which swaps the two cases at no extra cost. With k controls
        ("direct") it is the same idea over all 2^k sign patterns (`controls.add_controlled_rz`).

        Returns the target and control soquets.
        """
        theta = self._rotation_angle
        if ctrl is None:
            return bb.add(Rz(theta), q=target), None
        target, ctrl = controls.add_controlled_rz(bb, theta, target, ctrl, self.ctrl_values)
        return target, np.array(ctrl)

    def _add_identity_phase(self, bb: BloqBuilder, ctrl: Optional[NDArray[Soquet]]):
        """Add exp(-i phi), phi = theta/2, the evolution of an all-identity string.

        Uncontrolled it is a `GlobalPhase`. Under one control it is a phase on the control qubit:
        diag(1, e^{-i phi}) = e^{-i phi/2} Rz(-phi) for a control on |1>, and
        diag(e^{-i phi}, 1) = e^{-i phi/2} Rz(phi) for a control on |0>. With k controls
        ("direct") see `controls.add_controlled_phase`.

        Returns the control soquets.
        """
        phi = self._rotation_angle / 2
        if ctrl is None:
            bb.add(GlobalPhase(np.exp(-1j * phi)))
            return None
        return np.array(controls.add_controlled_phase(bb, phi, ctrl, self.ctrl_values))

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        """The decomposition's gates and their counts, without building the decomposition.

        Qualtran's default builds the decomposition, which is slow when a Hamiltonian has
        thousands of distinct terms. The gates are exactly those `build_composite_bloq` adds,
        angles included. Example: `IXIZI` is 2 H, 2 CNOT and 1 Rz.
        """
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.and_ancilla_call_graph(self)
        active = self._active_qubits
        if not active:
            return set(controls.controlled_phase_call_graph(
                self._rotation_angle / 2, self.ctrl_values).items())
        counts = Counter()
        for i in active:
            if self.pauli_string[i] == 'Y':
                counts[SGate(is_adjoint=True)] += 1
                counts[SGate()] += 1
            if self.pauli_string[i] in 'XY':
                counts[Hadamard()] += 2
        if len(active) > 1:
            counts[CNOT()] += 2 * (len(active) - 1)
        counts += controls.controlled_rz_call_graph(self._rotation_angle, self.ctrl_values)
        return set(counts.items())

    def _t_complexity_(self) -> TComplexity:
        """T-complexity of the decomposition, counted by hand (no decomposition needed).

        Every Rz counts as 1 rotation, whatever its angle (plan D14), so the count doesn't depend
        on the coefficient or time. Uncontrolled: the basis changes and CNOT chain (Cliffords) and
        1 rotation; an all-identity string is free. With k controls ("direct", or k = 1): the same
        with the rotation replaced by 2^k rotations and 2^k CNOTs; an all-identity string is
        2^k - 1 rotations and 2^k - 2 CNOTs. With the "ancilla" method: the AND and its uncompute
        plus the singly-controlled count.
        """
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.and_ancilla_t_complexity(self)
        active = self._active_qubits
        num_ctrls = len(self.ctrl_values)
        if not active:
            return controls.controlled_phase_t_complexity(num_ctrls)
        num_basis = sum(2 if self.pauli_string[i] == 'X' else 4
                        for i in active if self.pauli_string[i] != 'Z')
        num_cnots = 2 * (len(active) - 1)
        return TComplexity(clifford=num_basis + num_cnots) + controls.controlled_rz_t_complexity(
            num_ctrls)

    def my_static_costs(self, cost_key: 'CostKey'):
        """`QubitCount` without decomposing: qubits + controls + AND ancillas (if any)."""
        if isinstance(cost_key, QubitCount):
            return controls.peak_qubits(self.num_qubits, self.ctrl_values, self.multi_control_method)
        return NotImplemented

    def get_ctrl_system(
        self, ctrl_spec: Optional[CtrlSpec] = None
    ) -> Tuple['PauliStringEvolution', 'AddControlledT']:
        """Controlled version: the same bloq with the new control value prepended to `ctrl_values`.

        Only single-qubit controls on |1> or |0> are supported; others raise ValueError.
        """
        return controls.get_ctrl_system(self, ctrl_spec)

    def adjoint(self) -> 'PauliStringEvolution':
        """Return U^dagger as a PauliStringEvolution (flips `is_adjoint`; `time` is unchanged)."""
        return attrs.evolve(self, is_adjoint=not self.is_adjoint)

    def __pow__(self, power: int) -> 'PauliStringEvolution':
        """U**k is the same evolution for k times the time; U**(-k) is (U**k).adjoint().

        Raises ValueError for 0 and for non-integer powers.
        """
        if isinstance(power, bool) or not isinstance(power, numbers.Integral) or power == 0:
            raise ValueError(f"Power must be a non-zero integer, got {power}")
        return attrs.evolve(
            self, time=self.time * abs(int(power)), is_adjoint=self.is_adjoint != (power < 0))

    def __str__(self) -> str:
        text = f"exp(-i * {self.coefficient} * {self.pauli_string} * {self.time} / {self.hbar})"
        return controls.ctrl_str(self.ctrl_values, f"[{text}]†" if self.is_adjoint else text)

    def __repr__(self) -> str:
        adjoint = ", is_adjoint=True" if self.is_adjoint else ""
        ctrl = controls.ctrl_repr(self.ctrl_values, self.multi_control_method)
        return (f"PauliStringEvolution(pauli_string='{self.pauli_string}', "
                f"coefficient={self.coefficient}, time={self.time}, hbar={self.hbar}{adjoint}"
                f"{ctrl})")
