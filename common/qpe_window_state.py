"""
Window states for the control (phase) register of quantum phase estimation.

Backport of `qualtran.bloqs.phase_estimation.qpe_window_state` from Qualtran releases newer than
0.4.0, so that QHAT's QPE framework matches the newer interface.

Migration note: upstream window states use a RIGHT (allocating) `qpe_reg` register. Here
`qpe_reg` is THRU, so a window state is a unitary that maps |0...0> to the window state and the
full QPE bloq has a square matrix, as with Qualtran 0.4.0's `TextbookQPE`.
"""

import abc
from functools import cached_property
from typing import Dict

import attrs

from qualtran import Bloq, BloqBuilder, QDType, QFxp, Register, Side, Signature, SoquetT
from qualtran.bloqs.basic_gates import Hadamard, OnEach
from qualtran.symbolics import ceil, log2, pi, SymbolicFloat, SymbolicInt


@attrs.frozen
class QPEWindowStateBase(Bloq, metaclass=abc.ABCMeta):
    """Base class for window states prepared on the QPE control register `qpe_reg`."""

    @cached_property
    def m_qdtype(self) -> QDType:
        return QFxp(self.m_bits, self.m_bits)

    @cached_property
    def m_register(self) -> Register:
        return Register('qpe_reg', self.m_qdtype, side=Side.THRU)

    @property
    @abc.abstractmethod
    def m_bits(self) -> SymbolicInt: ...


@attrs.frozen
class RectangularWindowState(QPEWindowStateBase):
    """Window state used in textbook QPE: a Hadamard on every qubit of the control register.

    Args:
        bitsize: Size of the control register to prepare the window state on.

    Registers:
        qpe_reg: A `bitsize`-qubit THRU register of type `QFxp(bitsize, bitsize)`.
    """

    bitsize: SymbolicInt

    @property
    def m_bits(self) -> SymbolicInt:
        return self.bitsize

    @cached_property
    def signature(self) -> Signature:
        return Signature([self.m_register])

    @classmethod
    def from_precision_and_delta(cls, precision: SymbolicInt, delta: SymbolicFloat):
        r"""Estimate $\varphi$ to `precision` bits with probability of failure at most $\delta$.

        Uses Eq. 5.35 of Nielsen & Chuang: `m = n + ceil(log2(2 + 1/(2*delta)))`.
        """
        return cls(precision + ceil(log2(2 + 1 / (2 * delta))))

    @classmethod
    def from_standard_deviation_eps(cls, eps: SymbolicFloat):
        r"""Bound the standard deviation of the estimated phase $\phi$ by $\epsilon$.

        Textbook QPE has standard deviation at most $\pi / \sqrt{2^m}$, so
        `m = ceil(2*log2(pi/eps))`.
        """
        return cls(ceil(2 * log2(pi(eps) / eps)))

    def build_composite_bloq(self, bb: BloqBuilder, qpe_reg: SoquetT) -> Dict[str, SoquetT]:
        qpe_reg = bb.add(OnEach(self.m_bits, Hadamard()), q=qpe_reg)
        return {'qpe_reg': qpe_reg}
