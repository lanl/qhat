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
from qualtran.symbolics import ceil, is_symbolic, log2, SymbolicFloat, SymbolicInt


def precision_bits(phase_error: float) -> int:
    """Smallest n with 2**-n <= phase_error (tolerant of floating-point round-off)."""
    if not 0 < phase_error < 1:
        raise ValueError(f"phase_error must be in (0, 1), got {phase_error}.")
    # Round-off can turn an intended 2**-n into e.g. 0.12499999999999999, making log2 a few ulps
    # above n so that ceil adds a spurious bit; subtracting 1e-12 absorbs this.
    return max(0, math.ceil(math.log2(1 / phase_error) - 1e-12))


@attrs.frozen
class QPEWindowStateBase(Bloq, metaclass=abc.ABCMeta):
    """Base class for window states prepared on the QPE phase register `qpe_reg`.

    Subclasses define `m_bits` and `signature` (normally `Signature([self.m_register])`).

    `m_bits`, `m_qdtype` and `m_register` match Qualtran 0.7.0's `QPEWindowStateBase`, except
    that `m_register` is THRU here (see the module docstring). `from_requirements` and
    `from_num_phase_qubits` are QHAT additions, used by `FlexibleQPE`'s constructors of the same
    names.
    """

    @cached_property
    def m_qdtype(self) -> QDType:
        """Data type of `qpe_reg`: `QFxp(m_bits, m_bits)`, an unsigned fixed-point number with
        all `m_bits` bits after the binary point.

        The register value therefore reads directly as a phase in turns: with `m_bits = 3`, the
        bits 101 mean 0.101 in binary, i.e. 5/8 of a turn.
        """
        return QFxp(self.m_bits, self.m_bits)

    @cached_property
    def m_register(self) -> Register:
        """The phase register: named `qpe_reg`, of type `m_qdtype`, and THRU (Qualtran 0.7.0 uses
        RIGHT). `FlexibleQPE` takes its `qpe_reg` register from the window state's signature.
        """
        return Register('qpe_reg', self.m_qdtype, side=Side.THRU)

    @property
    @abc.abstractmethod
    def m_bits(self) -> SymbolicInt:
        """Number of qubits in the phase register; this is the "m" in `m_qdtype` and
        `m_register`. `FlexibleQPE` requires `num_ancilla_qubits == m_bits`.
        """

    @classmethod
    @abc.abstractmethod
    def from_requirements(
        cls, phase_error: float, probability_of_failure: float
    ) -> 'QPEWindowStateBase':
        """Return a window state of this class for which QPE with the textbook inverse QFT has
        Pr[|estimated phase - true phase| > phase_error] <= probability_of_failure.

        `phase_error` is in turns. E.g. `RectangularWindowState.from_requirements(1/8, 0.1)` has
        `m_bits = 6`: 3 bits for the precision plus 3 for the confidence.
        """

    @classmethod
    @abc.abstractmethod
    def from_num_phase_qubits(cls, num_phase_qubits: int) -> 'QPEWindowStateBase':
        """Return a window state of this class with `m_bits == num_phase_qubits`.

        Subclasses that have other parameters (e.g. a Kaiser window's `alpha`) choose them here.
        """


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
        """`from_precision_and_delta`, with the precision in bits computed from `phase_error`."""
        return cls.from_precision_and_delta(precision_bits(phase_error), probability_of_failure)

    @classmethod
    def from_num_phase_qubits(cls, num_phase_qubits: int):
        return cls(num_phase_qubits)

    @classmethod
    def from_precision_and_delta(cls, precision: SymbolicInt, delta: SymbolicFloat):
        r"""Estimate $\varphi$ to `precision` bits with probability of failure at most $\delta$.

        Uses Eq. 5.35 of Nielsen & Chuang: `m = n + ceil(log2(2 + 1/(2*delta)))`. `delta` may be
        a number, which must be in (0, 1), or a sympy expression.
        """
        if not is_symbolic(delta) and not 0 < delta < 1:
            raise ValueError(f"probability of failure must be in (0, 1), got {delta}.")
        return cls(precision + ceil(log2(2 + 1 / (2 * delta))))

    def build_composite_bloq(self, bb: BloqBuilder, qpe_reg: SoquetT) -> Dict[str, SoquetT]:
        qpe_reg = bb.add(OnEach(self.m_bits, Hadamard()), q=qpe_reg)
        return {'qpe_reg': qpe_reg}
