"""
Window states for the control (phase) register of quantum phase estimation.

Backport of `qualtran.bloqs.phase_estimation.qpe_window_state` from Qualtran releases newer than
0.4.0, so that QHAT's QPE framework matches the newer interface.

Migration note: upstream window states use a RIGHT (allocating) `qpe_reg` register. Here
`qpe_reg` is THRU, so a window state is a unitary that maps |0...0> to the window state and the
full QPE bloq has a square matrix, as with Qualtran 0.4.0's `TextbookQPE`.

Requirements convention: every window state's `from_requirements(phase_error,
probability_of_failure)` returns a window state that, used in QPE with the textbook inverse QFT,
guarantees Pr[|estimated phase - true phase| > phase_error] <= probability_of_failure. Phases are
measured in turns, i.e. in [0, 1), and distances are taken modulo 1.
"""

import abc
import math
from functools import cached_property
from typing import Dict

import attrs

from qualtran import Bloq, BloqBuilder, QDType, QFxp, Register, Side, Signature, SoquetT
from qualtran.bloqs.basic_gates import Hadamard, OnEach
from qualtran.symbolics import ceil, log2, pi, SymbolicFloat, SymbolicInt


def precision_bits(phase_error: float) -> int:
    """Smallest n with 2**-n <= phase_error (tolerant of floating-point round-off)."""
    if not 0 < phase_error < 1:
        raise ValueError(f"phase_error must be in (0, 1), got {phase_error}.")
    return max(0, math.ceil(math.log2(1 / phase_error) - 1e-9))


def textbook_confidence_bits(probability_of_failure: float) -> int:
    """Extra bits for textbook QPE to succeed with probability >= 1 - probability_of_failure.

    Nielsen & Chuang Eq. 5.35: ceil(log2(2 + 1/(2*delta))).
    """
    if not 0 < probability_of_failure < 1:
        raise ValueError(
            f"probability_of_failure must be in (0, 1), got {probability_of_failure}.")
    return math.ceil(math.log2(2 + 1 / (2 * probability_of_failure)))


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

    @classmethod
    @abc.abstractmethod
    def from_requirements(
        cls, phase_error: float, probability_of_failure: float
    ) -> 'QPEWindowStateBase':
        """Return the window state meeting the requirements (see the module docstring)."""

    @classmethod
    @abc.abstractmethod
    def from_num_phase_qubits(cls, num_phase_qubits: int) -> 'QPEWindowStateBase':
        """Return the window state on exactly `num_phase_qubits` qubits."""


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
    def from_requirements(cls, phase_error: float, probability_of_failure: float):
        """Precision bits for `phase_error` plus textbook confidence bits (N&C Eq. 5.35)."""
        return cls(precision_bits(phase_error)
                   + textbook_confidence_bits(probability_of_failure))

    @classmethod
    def from_num_phase_qubits(cls, num_phase_qubits: int):
        return cls(num_phase_qubits)

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
