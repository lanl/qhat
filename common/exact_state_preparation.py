"""
Exact preparation of an arbitrary state from |0...0>.

`ExactStatePreparation(amplitudes)` builds the state with a tree of uniformly controlled rotations
(Mottonen et al., quant-ph/0407010; Shende, Bullock & Markov, quant-ph/0406176). Its gate counts
depend only on the number of qubits and on whether the amplitudes are real, so they can be found
without computing the amplitudes (see `state_preparation_call_graph`).
"""

from collections import Counter
from functools import cached_property
from typing import Dict, List, Set, Tuple, TYPE_CHECKING

import attrs
import numpy as np

from qualtran import Bloq, BloqBuilder, QAny, Register, Signature, SoquetT
from qualtran.bloqs.basic_gates import CNOT, GlobalPhase, Ry, Rz

if TYPE_CHECKING:
    import quimb.tensor as qtn
    from qualtran.resource_counting import BloqCountT, SympySymbolAllocator


def state_preparation_call_graph(
    num_qubits: int, real: bool, ssa: 'SympySymbolAllocator'
) -> Set['BloqCountT']:
    """Gate counts of `ExactStatePreparation` on `num_qubits` qubits, with generic angles.

    Real amplitudes: 2^n - 1 Ry and 2^n - 2 CNOT. Complex amplitudes add 2^n - 1 Rz, 2^n - 2
    CNOT and one global phase. E.g. 3 qubits, real: {(Ry, 7), (CNOT, 6)}.
    """
    n = num_qubits
    counts = Counter({Ry(ssa.new_symbol('theta')): 2**n - 1, CNOT(): 2**n - 2})
    if not real:
        counts.update({Rz(ssa.new_symbol('phi')): 2**n - 1, CNOT(): 2**n - 2,
                       GlobalPhase(ssa.new_symbol('z')): 1})
    return {(bloq, count) for bloq, count in counts.items() if count > 0}


def _gray(i: int) -> int:
    return i ^ (i >> 1)


def _multiplexor_angles(alpha: np.ndarray) -> np.ndarray:
    """Rotation angles theta for the Gray-code circuit of a uniformly controlled rotation.

    `alpha[c]` is the angle wanted when the k control qubits hold c. The circuit applies
    R(theta[0]), CNOT, R(theta[1]), CNOT, ..., and the CNOTs flip the sign of theta[i] for control
    value c when popcount(c & gray(i)) is odd, so alpha = M theta with
    M[c, i] = (-1)^popcount(c & gray(i)), and M^T M = 2^k I.
    """
    size = len(alpha)
    c = np.arange(size)[:, None]
    g = np.array([_gray(i) for i in range(size)])[None, :]
    parity = np.vectorize(lambda x: bin(x).count('1') % 2)(c & g)
    return (1 - 2 * parity).T @ alpha / size


def _rotation_matrices(gate, alpha: np.ndarray) -> np.ndarray:
    """2x2 matrices of `gate(a)` for each angle a in `alpha` (Qualtran's Ry and Rz conventions)."""
    c, s = np.cos(alpha / 2), np.sin(alpha / 2)
    if gate is Ry:
        return np.stack([np.stack([c, -s], -1), np.stack([s, c], -1)], -2).astype(complex)
    zero = np.zeros_like(alpha)
    return np.stack([np.stack([c - 1j * s, zero], -1), np.stack([zero, c + 1j * s], -1)], -2)


def _cnot_control_bit(i: int, k: int) -> int:
    """Bit of the control value (0 = least significant) that CNOT i of a k-control multiplexor
    is controlled on: the bit where gray(i) and gray(i + 1) differ, with gray(2^k) taken as 0."""
    nxt = _gray(i + 1) if i + 1 < 2**k else 0
    return (_gray(i) ^ nxt).bit_length() - 1


@attrs.frozen
class ExactStatePreparation(Bloq):
    """Unitary mapping |0...0> to `sum_x amplitudes[x] |x>`, exactly (up to floating point).

    Qubit q[0] is the most significant bit of x. Level k of the circuit rotates q[k], controlled
    on q[0..k-1], so that each prefix of x gets the right total probability: Ry rotations for the
    magnitudes, then (for complex amplitudes) Rz rotations and a global phase for the phases.
    For real amplitudes the last Ry level also sets the signs, so no Rz is needed.

    Args:
        amplitudes: 2^n amplitudes of unit norm (stored as a tuple so the bloq is hashable).

    Registers:
        q: An n-qubit THRU register.
    """

    amplitudes: Tuple[complex, ...] = attrs.field(converter=tuple)

    def __attrs_post_init__(self):
        n = len(self.amplitudes)
        if n < 2 or n & (n - 1):
            raise ValueError(f"Need 2^n amplitudes with n >= 1, got {n}.")
        norm = np.linalg.norm(self.amplitudes)
        if abs(norm - 1) > 1e-10:
            raise ValueError(f"Amplitudes must have unit norm, got norm {norm}.")

    @cached_property
    def num_qubits(self) -> int:
        return len(self.amplitudes).bit_length() - 1

    @cached_property
    def is_real(self) -> bool:
        return not np.any(np.imag(self.amplitudes))

    @cached_property
    def signature(self) -> Signature:
        return Signature([Register('q', QAny(self.num_qubits))])

    @cached_property
    def ry_angles(self) -> List[np.ndarray]:
        """`ry_angles[k][c]`: angle of the Ry on q[k] when q[0..k-1] hold c (2^k angles).

        Each angle sends the probability of prefix c to its two children: Ry(alpha)|0> =
        cos(alpha/2)|0> + sin(alpha/2)|1>. Example: amplitudes (0.6, 0, 0, 0.8) give
        ry_angles = [[2 atan2(0.8, 0.6)], [0, pi]].
        """
        a = np.asarray(self.amplitudes)
        weights = [np.abs(a)**2]  # weights[j]: probability of each prefix of n - j bits
        for _ in range(self.num_qubits - 1):
            weights.append(weights[-1].reshape(-1, 2).sum(axis=1))
        weights.reverse()  # now weights[k] has 2^(k+1) entries
        angles = [2 * np.arctan2(*np.sqrt(w.reshape(-1, 2)).T[::-1]) for w in weights]
        if self.is_real:
            pairs = np.real(a).reshape(-1, 2)
            angles[-1] = 2 * np.arctan2(pairs[:, 1], pairs[:, 0])
        return angles

    @cached_property
    def rz_angles_and_global_phase(self) -> Tuple[List[np.ndarray], float]:
        """`(rz, omega)`: `rz[k][c]` is the angle of the Rz on q[k] when q[0..k-1] hold c, and
        e^{i omega} is the leftover global phase. Empty `rz` and omega = 0 for real amplitudes.

        Rz(beta) = diag(e^{-i beta/2}, e^{i beta/2}) splits a common phase between two children.
        Example: amplitudes (1, 1j)/sqrt(2) give rz = [[pi/2]] and omega = pi/4.
        """
        if self.is_real:
            return [], 0.0
        phases = np.angle(self.amplitudes)
        rz = []
        for _ in range(self.num_qubits):
            pairs = phases.reshape(-1, 2)
            rz.append(pairs[:, 1] - pairs[:, 0])
            phases = pairs.mean(axis=1)
        rz.reverse()
        return rz, float(phases[0])

    def _multiplexors(self):
        """(gate class, target qubit k, angles alpha[c]) for every level, in circuit order."""
        levels = [(Ry, k, alpha) for k, alpha in enumerate(self.ry_angles)]
        rz, _ = self.rz_angles_and_global_phase
        levels += [(Rz, k, alpha) for k, alpha in enumerate(rz)]
        return levels

    def build_composite_bloq(self, bb: BloqBuilder, q: SoquetT) -> Dict[str, SoquetT]:
        qs = bb.split(q)
        for gate, k, alpha in self._multiplexors():
            theta = _multiplexor_angles(alpha)
            for i, angle in enumerate(theta):
                qs[k] = bb.add(gate(float(angle)), q=qs[k])
                if k > 0:
                    ctrl = k - 1 - _cnot_control_bit(i, k)
                    qs[ctrl], qs[k] = bb.add(CNOT(), ctrl=qs[ctrl], target=qs[k])
        if not self.is_real:
            bb.add(GlobalPhase(np.exp(1j * self.rz_angles_and_global_phase[1])))
        return {'q': bb.join(qs, dtype=QAny(self.num_qubits))}

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        return state_preparation_call_graph(self.num_qubits, self.is_real, ssa)

    def tensor_contract(self) -> np.ndarray:
        """Dense unitary, as the product of the levels' block-diagonal matrices (no tensor
        network). Level k acts as sum_c |c><c| (x) R(alpha[c]) (x) I on (q[0..k-1], q[k], rest).
        """
        n = self.num_qubits
        u = np.eye(2**n, dtype=np.complex128)
        for gate, k, alpha in self._multiplexors():
            blocks = _rotation_matrices(gate, np.asarray(alpha, dtype=float))
            u = np.einsum('cab,cbrj->carj', blocks, u.reshape(2**k, 2, 2**(n - k - 1), 2**n))
            u = u.reshape(2**n, 2**n)
        return np.exp(1j * self.rz_angles_and_global_phase[1]) * u

    def add_my_tensors(
        self,
        tn: 'qtn.TensorNetwork',
        tag: object,
        *,
        incoming: Dict[str, SoquetT],
        outgoing: Dict[str, SoquetT],
    ):
        """Add `tensor_contract()` as one dense tensor, also when inside a larger bloq."""
        import quimb.tensor as qtn

        from qualtran._infra.composite_bloq import _flatten_soquet_collection
        from qualtran.simulation.tensor._tensor_data_manipulation import (
            tensor_shape_from_signature,
        )

        data = self.tensor_contract().reshape(tensor_shape_from_signature(self.signature))
        in_ind = _flatten_soquet_collection(incoming[reg.name] for reg in self.signature.lefts())
        out_ind = _flatten_soquet_collection(outgoing[reg.name] for reg in self.signature.rights())
        tn.add(qtn.Tensor(data=data, inds=out_ind + in_ind, tags=[self.pretty_name(), tag]))

    def __str__(self) -> str:
        return f'ExactStatePreparation[{self.num_qubits}]'
