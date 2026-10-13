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
import numbers
from functools import cached_property
from typing import Dict, Optional

import attrs

from qualtran import Bloq, BloqBuilder, QDType, QFxp, Register, Side, Signature, SoquetT
from qualtran.bloqs.basic_gates import Hadamard, OnEach
from qualtran.symbolics import ceil, is_symbolic, log2, SymbolicFloat, SymbolicInt

from qhat.common.hermiticity import require_real_scalar


def precision_bits(phase_error: float) -> int:
    """Smallest n with 2**-n <= phase_error (tolerant of floating-point round-off)."""
    phase_error = require_real_scalar("phase_error", phase_error)
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

        Subclasses that have other parameters (e.g. a Kaiser window's `alpha`) choose them here,
        or raise a ValueError if the total alone is not enough (e.g. the tapered windows).
        """

    @classmethod
    def from_precision_and_confidence(
        cls, num_precision_qubits: int, num_confidence_qubits: int
    ) -> 'QPEWindowStateBase':
        """Return a window state with `num_precision_qubits` + `num_confidence_qubits` qubits.

        By default only the total matters (`from_num_phase_qubits`); the tapered windows use the
        split (see `tapered_window_state.py`).
        """
        return cls.from_num_phase_qubits(num_precision_qubits + num_confidence_qubits)

    @classmethod
    def from_qubit_counts(
        cls,
        num_phase_qubits: Optional[int] = None,
        num_precision_qubits: Optional[int] = None,
        num_confidence_qubits: Optional[int] = None,
    ) -> 'QPEWindowStateBase':
        """Return a window state from explicit qubit counts, where
        num_phase_qubits = num_precision_qubits + num_confidence_qubits.

        Give `num_phase_qubits` alone (`from_num_phase_qubits`), or any two of the three counts
        (`from_precision_and_confidence`). E.g. `from_qubit_counts(6, num_confidence_qubits=2)`
        is `from_precision_and_confidence(4, 2)`.
        """
        given = {'num_phase_qubits': num_phase_qubits,
                 'num_precision_qubits': num_precision_qubits,
                 'num_confidence_qubits': num_confidence_qubits}
        given = {name: value for name, value in given.items() if value is not None}
        for name, value in given.items():
            if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer, got {value!r}.")
        if set(given) == {'num_phase_qubits'}:
            return cls.from_num_phase_qubits(num_phase_qubits)
        if len(given) == 3:
            if num_precision_qubits + num_confidence_qubits != num_phase_qubits:
                raise ValueError(
                    f"num_phase_qubits ({num_phase_qubits}) must equal num_precision_qubits "
                    f"({num_precision_qubits}) + num_confidence_qubits "
                    f"({num_confidence_qubits}).")
        elif len(given) != 2:
            raise ValueError(
                f"Cannot size {cls.__name__} from {sorted(given) or 'no qubit counts'}. Give two "
                "of num_phase_qubits, num_precision_qubits and num_confidence_qubits; some "
                "window states (e.g. RectangularWindowState) also accept num_phase_qubits alone.")
        elif num_phase_qubits is not None:
            other = 'num_precision_qubits' if num_confidence_qubits is None else \
                'num_confidence_qubits'
            if given[other] > num_phase_qubits:
                raise ValueError(f"{other} ({given[other]}) exceeds num_phase_qubits "
                                 f"({num_phase_qubits}).")
            if num_precision_qubits is None:
                num_precision_qubits = num_phase_qubits - num_confidence_qubits
            else:
                num_confidence_qubits = num_phase_qubits - num_precision_qubits
        return cls.from_precision_and_confidence(num_precision_qubits, num_confidence_qubits)


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
