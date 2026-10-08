"""
Quantum phase estimation with a pluggable window state, unitary, and inverse QFT.

`FlexibleQPE` is modeled on Qualtran's `TextbookQPE`. Its fields (the unitary, the phase-register
size, then phase-register components that default to textbook QPE) follow Qualtran 0.4.0, except
that `m_bits` and `ctrl_state_prep` are renamed `num_ancilla_qubits` and `ancilla_prep`, and
`ancilla_prep` is a window state, as in newer Qualtran releases. Its registers are `qpe_reg` plus
the unitary's registers. It differs from Qualtran 0.4.0's `TextbookQPE` in that it:

- decomposes natively as a Bloq, so it supports `decompose_bloq()` and `tensor_contract()` for
  unitaries that are plain Bloqs (0.4.0's cirq-based `cirq.pow` path fails for those);
- fast-forwards U^(2^j) via the unitary's `__pow__` when available;
- uses the structure of QPE to speed up tensor contraction, call graphs, and qubit counts.

See `qpe_window_state.py` for the one intentional deviation from upstream (THRU `qpe_reg`).
"""

from collections import Counter
from functools import cached_property
from typing import Callable, Dict, Optional, Set, Tuple, Type, TYPE_CHECKING

import attrs
import numpy as np

from qualtran import Bloq, BloqBuilder, CtrlSpec, GateWithRegisters, Register, Signature, SoquetT
from qualtran.bloqs.basic_gates import Power
from qualtran.bloqs.qft.qft_text_book import QFTTextBook
from qualtran.resource_counting import get_cost_value, QubitCount
from qualtran.symbolics import is_symbolic, SymbolicInt

from qhat.common.qpe_window_state import QPEWindowStateBase, RectangularWindowState

if TYPE_CHECKING:
    import quimb.tensor as qtn
    from qualtran.resource_counting import BloqCountT, CostKey, SympySymbolAllocator


@attrs.frozen
class FlexibleQPE(GateWithRegisters):
    r"""Phase estimation of `unitary` (Nielsen & Chuang, Sec. 5.2) with pluggable components.

    ```
           ┌─────────┐                              ┌─────────┐
      |0> -│         │-----------------------@------│         │---M--- [m1]:highest bit
           │         │                       |      │         │
      |0> -│         │-----------------@-----+------│         │---M--- [m2]
           │ Ancilla │                 |     |      │ QFT_inv │
      |0> -│  Prep   │-----------@-----+-----+------│         │---M--- [m3]
           │         │           |     |     |      │         │
      |0> -│         │-----@-----+-----+-----+------│         │---M--- [m4]:lowest bit
           └─────────┘     |     |     |     |      └─────────┘
    |Psi> -----------------U-----U^2---U^4---U^8---------------------- |Psi>
    ```

    U^(2^j) is `unitary ** 2**j` if the unitary defines `__pow__` (fast-forwarding), and
    otherwise `Power(unitary, 2**j)` (2^j repetitions). A unitary's `__pow__` must return a bloq
    implementing exactly U^k; e.g. `Trotterization.__pow__` scales both the step count and time.
    `tensor_contract` relies on this: it squares U's matrix instead of contracting each U^(2^j).

    Args:
        unitary: Bloq (THRU registers only) whose eigenphases are estimated.
        num_ancilla_qubits: Number of qubits in the phase (ancilla) register (Qualtran's
            `m_bits`).
        ancilla_prep: Window state prepared on the phase register. Defaults to
            `RectangularWindowState(num_ancilla_qubits)` (a Hadamard on each phase qubit).
        qft_inv: Inverse QFT on the phase register. Defaults to
            `QFTTextBook(num_ancilla_qubits, with_reverse=True).adjoint()`.

    Registers:
        qpe_reg: Phase register of `num_ancilla_qubits` qubits, of type `QFxp(m, m)` with
            `m = num_ancilla_qubits`; must start in |0...0>.
        target registers: All registers of `unitary.signature` (`num_state_qubits` qubits).
    """

    unitary: Bloq
    num_ancilla_qubits: SymbolicInt
    ancilla_prep: QPEWindowStateBase = attrs.field()
    qft_inv: Bloq = attrs.field()

    @ancilla_prep.default
    def _default_ancilla_prep(self):
        return RectangularWindowState(self.num_ancilla_qubits)

    @qft_inv.default
    def _default_inverse_qft(self):
        return QFTTextBook(self.num_ancilla_qubits, with_reverse=True).adjoint()

    def __attrs_post_init__(self):
        if not is_symbolic(self.num_ancilla_qubits):
            if self.ancilla_prep.m_bits != self.num_ancilla_qubits:
                raise ValueError(
                    f"ancilla_prep acts on {self.ancilla_prep.m_bits} qubits but the phase "
                    f"register has {self.num_ancilla_qubits}."
                )
            if self.qft_inv.signature.n_qubits() != self.num_ancilla_qubits:
                raise ValueError(
                    f"qft_inv acts on {self.qft_inv.signature.n_qubits()} qubits but the phase "
                    f"register has {self.num_ancilla_qubits}."
                )
        if any(reg.name == 'qpe_reg' for reg in self.unitary.signature):
            raise ValueError("The unitary may not have a register named 'qpe_reg'.")

    @cached_property
    def num_state_qubits(self) -> int:
        """Number of qubits the unitary acts on (the state whose eigenphase is estimated)."""
        return sum(reg.total_bits() for reg in self.target_registers)

    @cached_property
    def num_total_qubits(self) -> SymbolicInt:
        """Total width of the registers, `num_ancilla_qubits + num_state_qubits`.

        This excludes any scratch qubits allocated inside the unitary; for the peak number of
        qubits in use, see `get_cost_value(qpe, QubitCount())`.
        """
        return self.num_ancilla_qubits + self.num_state_qubits

    @cached_property
    def target_registers(self) -> Tuple[Register, ...]:
        return tuple(self.unitary.signature)

    @cached_property
    def phase_registers(self) -> Tuple[Register, ...]:
        return tuple(self.ancilla_prep.signature)

    @cached_property
    def signature(self) -> Signature:
        return Signature([*self.phase_registers, *self.target_registers])

    def __str__(self) -> str:
        return f'FlexibleQPE[{self.num_ancilla_qubits}]'

    @classmethod
    def from_requirements(
        cls,
        unitary: Bloq,
        phase_error: float,
        probability_of_failure: float,
        *,
        ancilla_prep: Type[QPEWindowStateBase] = RectangularWindowState,
        qft_inv: Optional[Callable[[int], Bloq]] = None,
    ) -> 'FlexibleQPE':
        """Build a QPE meeting Pr[|phase estimate error| > phase_error] <= probability_of_failure.

        Phases are in turns ([0, 1)). Sizing is delegated to the `ancilla_prep` class.
        `qft_inv`, if given, maps the phase-register size to an inverse-QFT bloq.
        """
        prep = ancilla_prep.from_requirements(phase_error, probability_of_failure)
        return cls._from_window_state(unitary, prep, qft_inv)

    @classmethod
    def from_num_phase_qubits(
        cls,
        unitary: Bloq,
        num_phase_qubits: int,
        *,
        ancilla_prep: Type[QPEWindowStateBase] = RectangularWindowState,
        qft_inv: Optional[Callable[[int], Bloq]] = None,
    ) -> 'FlexibleQPE':
        """Build a QPE with exactly `num_phase_qubits` phase qubits."""
        prep = ancilla_prep.from_num_phase_qubits(num_phase_qubits)
        return cls._from_window_state(unitary, prep, qft_inv)

    @classmethod
    def _from_window_state(cls, unitary, prep, qft_inv):
        if qft_inv is None:
            return cls(unitary, prep.m_bits, prep)
        return cls(unitary, prep.m_bits, prep, qft_inv(prep.m_bits))

    def unitary_power(self, k: int) -> Bloq:
        """Return a bloq for U^k, fast-forwarded via `unitary.__pow__` if available."""
        if k == 1:
            return self.unitary
        if hasattr(type(self.unitary), '__pow__'):
            result = self.unitary**k
            if result is not NotImplemented:
                return result
        return Power(self.unitary, k)

    def controlled_power(self, j: int) -> Bloq:
        """Return the bloq for controlled U^(2^j)."""
        return self.unitary_power(2**j).controlled()

    def build_composite_bloq(
        self, bb: BloqBuilder, qpe_reg: SoquetT, **target_soqs: SoquetT
    ) -> Dict[str, SoquetT]:
        m = self.num_ancilla_qubits
        if is_symbolic(m):
            raise NotImplementedError(
                f"Cannot decompose {self} with symbolic num_ancilla_qubits.")
        qpe_reg = bb.add(self.ancilla_prep, qpe_reg=qpe_reg)
        qs = bb.split(qpe_reg)
        target_names = [reg.name for reg in self.target_registers]
        for j in range(m):
            # bb.split is big-endian, so the last qubit is the least significant bit.
            _, add_controlled = self.unitary_power(2**j).get_ctrl_system(CtrlSpec())
            (qs[m - 1 - j],), out_soqs = add_controlled(bb, [qs[m - 1 - j]], target_soqs)
            target_soqs = dict(zip(target_names, out_soqs))
        qpe_reg = bb.join(qs, dtype=self.phase_registers[0].dtype)
        qpe_reg = bb.add(self.qft_inv, q=qpe_reg)
        return {'qpe_reg': qpe_reg, **target_soqs}

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        m = self.num_ancilla_qubits
        if is_symbolic(m):
            # Assumes the unitary is not fast-forwardable.
            return {
                (self.ancilla_prep, 1),
                (self.unitary.controlled(), 2**m - 1),
                (self.qft_inv, 1),
            }
        counts = Counter([self.ancilla_prep, self.qft_inv])
        counts.update(self.controlled_power(j) for j in range(m))
        return set(counts.items())

    def my_static_costs(self, cost_key: 'CostKey'):
        m = self.num_ancilla_qubits
        if not isinstance(cost_key, QubitCount) or is_symbolic(m):
            return NotImplemented
        n = self.num_state_qubits
        widths = [
            get_cost_value(self.ancilla_prep, cost_key) + n,
            get_cost_value(self.qft_inv, cost_key) + n,
        ]
        for j in range(m):
            u_k = self.unitary_power(2**j)
            if isinstance(u_k, Power):
                # Same width as a single controlled U; avoids decomposing 2^j copies.
                u_k = self.unitary
            widths.append(get_cost_value(u_k.controlled(), cost_key) + m - 1)
        return max(widths)

    def _phase_blocks(self) -> np.ndarray:
        """Return W with W[k] = U^k for k = 0, 1, ..., 2^m - 1, where m = `num_ancilla_qubits`.

        `a` holds the powers of two, [U, U^2, U^4, ...]: the m controlled gates in the circuit.
        `W` holds every power, [U^0, U^1, U^2, U^3, ...]: when the phase register is in |k>, the
        gates for the set bits of k fire, applying U^k to the target. `tensor_contract` needs
        every W[k] because the phase register holds a superposition of all |k>.

        U is contracted once and `a` is built by squaring, rather than contracting each
        U^(2^j) bloq, which assumes `unitary_power(k)` is exactly U^k.

        Example (m = 2): a = [U, U^2] and W = [U^0, U^1, U^2, U^3], with U^3 = U^2 @ U.
        """
        m = self.num_ancilla_qubits
        a = [self.unitary.tensor_contract()]
        for _ in range(1, m):
            a.append(a[-1] @ a[-1])
        d = a[0].shape[0]
        w = np.empty((2**m, d, d), dtype=np.complex128)
        w[0] = np.eye(d)
        for k in range(1, 2**m):
            h = k.bit_length() - 1
            w[k] = a[h] @ w[k - 2**h]
        return w

    def tensor_contract(self) -> np.ndarray:
        """Return the dense matrix of this QPE, built from P = `ancilla_prep`,
        W = `_phase_blocks()`, and Q = `qft_inv` instead of contracting the decomposition.

        The output is filled one phase-register input column at a time, so peak memory is close
        to the size of the output (16 * 4^(m + n) bytes for m phase and n target qubits).
        """
        p = self.ancilla_prep.tensor_contract()
        q = self.qft_inv.tensor_contract()
        w = self._phase_blocks()
        num_phases, d, _ = w.shape
        w_flat = w.reshape(num_phases, d * d)
        out = np.empty((num_phases, d, num_phases, d), dtype=np.complex128)
        for b in range(num_phases):
            out[:, :, b, :] = (q @ (p[:, b, None] * w_flat)).reshape(num_phases, d, d)
        return out.reshape(num_phases * d, num_phases * d)

    def add_my_tensors(
        self,
        tn: 'qtn.TensorNetwork',
        tag: object,
        *,
        incoming: Dict[str, SoquetT],
        outgoing: Dict[str, SoquetT],
    ):
        """Add `tensor_contract()` as a single dense tensor, so that a FlexibleQPE inside a
        larger bloq also uses it instead of being decomposed."""
        import quimb.tensor as qtn

        from qualtran._infra.composite_bloq import _flatten_soquet_collection
        from qualtran.simulation.tensor._tensor_data_manipulation import (
            tensor_shape_from_signature,
        )

        data = self.tensor_contract().reshape(tensor_shape_from_signature(self.signature))

        in_ind = _flatten_soquet_collection(incoming[reg.name] for reg in self.signature.lefts())
        out_ind = _flatten_soquet_collection(outgoing[reg.name] for reg in self.signature.rights())
        tn.add(qtn.Tensor(data=data, inds=out_ind + in_ind, tags=[self.pretty_name(), tag]))
