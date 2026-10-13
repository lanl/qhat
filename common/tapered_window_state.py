"""
Tapered (DPSS-family) window states for QPE, from Patel, O'Gorman, Gidney et al.,
"Optimal coherent quantum phase estimation via tapering", arXiv:2403.18927.

A tapered window state is `sum_n phi[n] |n>` on p = l + m phase qubits (N = 2^p), with a smooth
real taper `phi` in place of the rectangular window. Its size has two parts:

- l = `num_precision_qubits`: success means the estimate is within delta = 2^-(l+1) of the true
  phase (in turns, modulo 1).
- m = `num_confidence_qubits`: sets the taper's bandwidth (K = 2^(m-1) - 1) and so how often the
  estimate misses delta.

Example: l = 3, m = 3 uses 6 phase qubits and gives delta = 1/16 with probability of failure about
1e-8; the rectangular window needs 26 extra qubits for that.

Classes:

- `DPSSWindowState`: the exact discrete prolate spheroidal sequence (paper Sec. IV C), prepared
  with `ExactStatePreparation` on all p qubits (2^p - 1 rotations).
"""

import abc
import itertools
import math
from functools import cached_property
from typing import ClassVar, Dict, Set, TYPE_CHECKING

import attrs
import numpy as np

from qualtran import BloqBuilder, Signature, SoquetT
from qualtran.resource_counting import QubitCount

from qhat.common.exact_state_preparation import (
    ExactStatePreparation,
    state_preparation_call_graph,
)
from qhat.common.qpe_window_state import precision_bits, QPEWindowStateBase

if TYPE_CHECKING:
    from qualtran.resource_counting import BloqCountT, CostKey, SympySymbolAllocator


def half_bandwidth(num_confidence_qubits: int) -> int:
    """K = 2^(m-1) - 1: the taper keeps the 2K + 1 frequencies -K..K. E.g. m = 3 gives K = 3."""
    return 2**(num_confidence_qubits - 1) - 1


def dpss_average_failure_bound(num_precision_qubits: int, num_confidence_qubits: int) -> float:
    """Theorem 1 of the paper: bound on the DPSS failure probability, averaged over the offset of
    the true phase from the N-point grid. E.g. l = 3, m = 5 gives about 8.6e-9.
    """
    n = 2**(num_precision_qubits + num_confidence_qubits)
    k = half_bandwidth(num_confidence_qubits)
    scale = 2 / math.pi**2
    return min(8 * math.exp(-(2 * k - 1) / (scale * math.log(4 * n))),
               10 * math.exp(-(2 * k - 6) / (scale * math.log(100 * k + 75))))


@attrs.frozen
class TaperedWindowStateBase(QPEWindowStateBase):
    """Base for window states `sum_n taper()[n] |n>` on l + m qubits (see the module docstring).

    Subclasses define `taper()`, the circuit, and `failure_bound_factor` (used by
    `confidence_qubits`).

    Args:
        num_precision_qubits: l >= 0; success means |estimate - phase| <= 2^-(l+1).
        num_confidence_qubits: m >= 1; the taper keeps frequencies -K..K, K = 2^(m-1) - 1.

    Registers:
        qpe_reg: An (l + m)-qubit THRU register of type `QFxp(l + m, l + m)`.
    """

    num_precision_qubits: int
    num_confidence_qubits: int

    failure_bound_factor: ClassVar[float] = 1.0

    def __attrs_post_init__(self):
        if self.num_precision_qubits < 0:
            raise ValueError(
                f"num_precision_qubits must be >= 0, got {self.num_precision_qubits}.")
        if self.num_confidence_qubits < 1:
            raise ValueError(
                f"num_confidence_qubits must be >= 1, got {self.num_confidence_qubits}.")

    @property
    def m_bits(self) -> int:
        return self.num_precision_qubits + self.num_confidence_qubits

    @cached_property
    def signature(self) -> Signature:
        return Signature([self.m_register])

    @property
    def half_bandwidth(self) -> int:
        return half_bandwidth(self.num_confidence_qubits)

    @abc.abstractmethod
    def taper(self) -> np.ndarray:
        """The 2^(l+m) real amplitudes phi[n] of the window state (unit norm, sum > 0); n is the
        integer value of `qpe_reg`, and the taper is symmetric about n = (N - 1)/2."""

    @classmethod
    def confidence_qubits(cls, num_precision_qubits: int, probability_of_failure: float) -> int:
        """Number of confidence qubits m for `from_requirements`.

        The smallest m with `failure_bound_factor * dpss_average_failure_bound(l, m)` <=
        `probability_of_failure`. Theorem 1 bounds the average over phases; tests check that the
        worst case also meets the target. The factor is 1 for the DPSS. E.g. l = 3,
        p_fail = 1e-8: m = 5.
        """
        return next(m for m in itertools.count(1)
                    if cls.failure_bound_factor * dpss_average_failure_bound(
                        num_precision_qubits, m) <= probability_of_failure)

    @classmethod
    def from_requirements(cls, phase_error: float, probability_of_failure: float):
        """l from delta = 2^-(l+1) <= phase_error (one bit fewer than the rectangular window,
        whose guarantee is |error| < 2^-l), then m from `confidence_qubits`.
        """
        if not 0 < probability_of_failure < 1:
            raise ValueError(
                f"probability of failure must be in (0, 1), got {probability_of_failure}.")
        l = precision_bits(phase_error) - 1
        return cls(l, cls.confidence_qubits(l, probability_of_failure))

    @classmethod
    def from_num_phase_qubits(cls, num_phase_qubits: int):
        raise ValueError(
            f"{cls.__name__} cannot be sized from num_phase_qubits alone: the taper depends on "
            "how the phase qubits split into precision qubits l and confidence qubits m. Give "
            "num_precision_qubits and num_confidence_qubits (either one may be replaced by "
            "num_phase_qubits = l + m).")

    @classmethod
    def from_precision_and_confidence(cls, num_precision_qubits: int,
                                      num_confidence_qubits: int):
        return cls(num_precision_qubits, num_confidence_qubits)

    def my_static_costs(self, cost_key: 'CostKey'):
        # No scratch qubits; avoids decomposing, which would compute the taper.
        if isinstance(cost_key, QubitCount):
            return self.m_bits
        return NotImplemented


def dpss_taper(num_precision_qubits: int, num_confidence_qubits: int) -> np.ndarray:
    """The DPSS of length N = 2^(l+m) and half-bandwidth W = (2K+1)/(2N), with unit norm and
    positive sum (scipy's `dpss(N, NW=K + 1/2)`)."""
    from scipy.signal.windows import dpss

    n = 2**(num_precision_qubits + num_confidence_qubits)
    phi = dpss(n, half_bandwidth(num_confidence_qubits) + 0.5, norm=2)
    return phi if phi.sum() > 0 else -phi


@attrs.frozen
class DPSSWindowState(TaperedWindowStateBase):
    """Window state with the exact DPSS taper (paper Sec. IV C), the taper with the smallest
    average failure probability. Prepared with `ExactStatePreparation` on all l + m qubits,
    i.e. 2^(l+m) - 1 Ry rotations and 2^(l+m) - 2 CNOTs.

    Args and registers: see `TaperedWindowStateBase`.
    """

    def taper(self) -> np.ndarray:
        return dpss_taper(self.num_precision_qubits, self.num_confidence_qubits)

    def build_composite_bloq(self, bb: BloqBuilder, qpe_reg: SoquetT) -> Dict[str, SoquetT]:
        return {'qpe_reg': bb.add(ExactStatePreparation(self.taper()), q=qpe_reg)}

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        # From the circuit's structure, so the taper is not computed.
        return state_preparation_call_graph(self.m_bits, True, ssa)
