"""
Quantum phase estimation with a pluggable window state, unitary, and inverse QFT.

`FlexibleQPE` is modeled on Qualtran's `TextbookQPE`. Its fields (`unitary`, `m_bits`, then
phase-register components that default to textbook QPE) follow Qualtran 0.4.0, except that
`ctrl_state_prep` is renamed `ancilla_prep` and is a window state, as in newer Qualtran releases.
Its registers are `qpe_reg` plus the unitary's registers. It differs from Qualtran 0.4.0's
`TextbookQPE` in that it:

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

    Args:
        unitary: Bloq (THRU registers only) whose eigenphases are estimated.
        m_bits: Number of qubits in the phase register.
        ancilla_prep: Window state prepared on the phase register. Defaults to
            `RectangularWindowState(m_bits)` (a Hadamard on each phase qubit).
        qft_inv: Inverse QFT on the phase register. Defaults to
            `QFTTextBook(m_bits, with_reverse=True).adjoint()`.

    Registers:
        qpe_reg: Phase register of type `QFxp(m_bits, m_bits)`; must start in |0...0>.
        target registers: All registers of `unitary.signature`.
    """

    unitary: Bloq
    m_bits: SymbolicInt
    ancilla_prep: QPEWindowStateBase = attrs.field()
    qft_inv: Bloq = attrs.field()

    @ancilla_prep.default
    def _default_ancilla_prep(self):
        return RectangularWindowState(self.m_bits)

    @qft_inv.default
    def _default_inverse_qft(self):
        return QFTTextBook(self.m_bits, with_reverse=True).adjoint()

    def __attrs_post_init__(self):
        if not is_symbolic(self.m_bits):
            if self.ancilla_prep.m_bits != self.m_bits:
                raise ValueError(
                    f"ancilla_prep acts on {self.ancilla_prep.m_bits} qubits but the phase "
                    f"register has {self.m_bits}."
                )
            if self.qft_inv.signature.n_qubits() != self.m_bits:
                raise ValueError(
                    f"qft_inv acts on {self.qft_inv.signature.n_qubits()} qubits but the phase "
                    f"register has {self.m_bits}."
                )
        if any(reg.name == 'qpe_reg' for reg in self.unitary.signature):
            raise ValueError("The unitary may not have a register named 'qpe_reg'.")

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
        return f'FlexibleQPE[{self.m_bits}]'

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
        if is_symbolic(self.m_bits):
            raise NotImplementedError(f"Cannot decompose {self} with symbolic m_bits.")
        qpe_reg = bb.add(self.ancilla_prep, qpe_reg=qpe_reg)
        qs = bb.split(qpe_reg)
        target_names = [reg.name for reg in self.target_registers]
        for j in range(self.m_bits):
            # bb.split is big-endian, so the last qubit is the least significant bit.
            _, add_controlled = self.unitary_power(2**j).get_ctrl_system(CtrlSpec())
            (qs[self.m_bits - 1 - j],), out_soqs = add_controlled(
                bb, [qs[self.m_bits - 1 - j]], target_soqs
            )
            target_soqs = dict(zip(target_names, out_soqs))
        qpe_reg = bb.join(qs, dtype=self.phase_registers[0].dtype)
        qpe_reg = bb.add(self.qft_inv, q=qpe_reg)
        return {'qpe_reg': qpe_reg, **target_soqs}

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        if is_symbolic(self.m_bits):
            # Assumes the unitary is not fast-forwardable.
            return {
                (self.ancilla_prep, 1),
                (self.unitary.controlled(), 2**self.m_bits - 1),
                (self.qft_inv, 1),
            }
        counts = Counter([self.ancilla_prep, self.qft_inv])
        counts.update(self.controlled_power(j) for j in range(self.m_bits))
        return set(counts.items())

    def my_static_costs(self, cost_key: 'CostKey'):
        if not isinstance(cost_key, QubitCount) or is_symbolic(self.m_bits):
            return NotImplemented
        n = sum(reg.total_bits() for reg in self.target_registers)
        widths = [
            get_cost_value(self.ancilla_prep, cost_key) + n,
            get_cost_value(self.qft_inv, cost_key) + n,
        ]
        for j in range(self.m_bits):
            u_k = self.unitary_power(2**j)
            if isinstance(u_k, Power):
                # Same width as a single controlled U; avoids decomposing 2^j copies.
                u_k = self.unitary
            widths.append(get_cost_value(u_k.controlled(), cost_key) + self.m_bits - 1)
        return max(widths)

    def _phase_blocks(self) -> np.ndarray:
        """Return W with W[k] = (the matrix applied to the target when qpe_reg = |k>)."""
        a = [self.unitary.tensor_contract()]
        for j in range(1, self.m_bits):
            u_k = self.unitary_power(2**j)
            a.append(a[-1] @ a[-1] if isinstance(u_k, Power) else u_k.tensor_contract())
        d = a[0].shape[0]
        w = np.empty((2**self.m_bits, d, d), dtype=np.complex128)
        w[0] = np.eye(d)
        for k in range(1, 2**self.m_bits):
            h = k.bit_length() - 1
            w[k] = a[h] @ w[k - 2**h]
        return w

    def add_my_tensors(
        self,
        tn: 'qtn.TensorNetwork',
        tag: object,
        *,
        incoming: Dict[str, SoquetT],
        outgoing: Dict[str, SoquetT],
    ):
        """Add a single dense tensor built from the block-diagonal structure of QPE.

        In the phase basis, the controlled powers act as U^k on the target when qpe_reg = |k>,
        so the full matrix is (Q ⊗ I) · diag_k(W_k) · (P ⊗ I).
        """
        import quimb.tensor as qtn

        from qualtran._infra.composite_bloq import _flatten_soquet_collection
        from qualtran.simulation.tensor._tensor_data_manipulation import (
            tensor_shape_from_signature,
        )

        p = self.ancilla_prep.tensor_contract()
        q = self.qft_inv.tensor_contract()
        w = self._phase_blocks()
        data = np.einsum('ak,kij,kb->aibj', q, w, p, optimize=True)
        data = data.reshape(tensor_shape_from_signature(self.signature))

        in_ind = _flatten_soquet_collection(incoming[reg.name] for reg in self.signature.lefts())
        out_ind = _flatten_soquet_collection(outgoing[reg.name] for reg in self.signature.rights())
        tn.add(qtn.Tensor(data=data, inds=out_ind + in_ind, tags=[self.pretty_name(), tag]))
