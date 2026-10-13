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

Classes (config names in `qhat.analysis.algorithm.ANCILLA_PREPS`):

- `DPSSWindowState` ("dpss"): the exact discrete prolate spheroidal sequence (paper Sec. IV C),
  prepared with `ExactStatePreparation` on all p qubits (2^p - 1 rotations).
- `TruncatedDPSSWindowState` ("truncated dpss"): the DPSS with only its 2K + 1 lowest frequencies
  kept (Sec. VII A), prepared with the paper's App. A circuit (2^m - 1 + m rotations plus a QFT).
- `BandwidthOptimizedWindowState` ("bandwidth-optimized"): the best taper on those same 2K + 1
  frequencies (Sec. VII B), with the same circuit.
"""

import abc
import itertools
import math
from collections import Counter
from fractions import Fraction
from functools import cached_property
from typing import ClassVar, Dict, Set, Tuple, TYPE_CHECKING

import attrs
import numpy as np

from qualtran import BloqBuilder, QAny, Signature, SoquetT
from qualtran.bloqs.basic_gates import CNOT, ZPowGate
from qualtran.bloqs.qft.qft_text_book import QFTTextBook
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
        worst case also meets the target. The factor is 1 for the DPSS and 4 for the band-limited
        tapers (paper eq. 107). E.g. l = 3, p_fail = 1e-8: m = 5 for the DPSS, 6 for the others.
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


@attrs.frozen
class BandLimitedWindowStateBase(TaperedWindowStateBase):
    """Base for tapers made of the 2K + 1 frequencies -K..K only, prepared with the paper's
    App. A circuit:

    1. `ExactStatePreparation` of the real frequency amplitudes on the m low qubits;
       frequency k sits at k mod 2^m (two's complement) and 2^(m-1) is left empty;
    2. a Z rotation on each of those m qubits, giving frequency k the phase e^{-i pi k (N-1)/N}
       that centres the taper at n = (N - 1)/2 (otherwise it wraps around n = 0 and fails);
    3. CNOTs from the m-register's top bit to the l high qubits (sign extension to p bits);
    4. `QFTTextBook(l + m, with_reverse=True)`.

    Cost: 2^m - 1 Ry, m Z rotations, 2^m - 2 + l CNOTs and one QFT; no scratch qubits.
    Subclasses define `frequency_amplitudes()`.

    Args and registers: see `TaperedWindowStateBase`.
    """

    # Eq. 107: truncating to 2K + 1 frequencies at most quadruples the DPSS failure probability.
    failure_bound_factor: ClassVar[float] = 4.0

    @abc.abstractmethod
    def frequency_amplitudes(self) -> np.ndarray:
        """The K + 1 real amplitudes r[0..K] of the centred taper; frequency -k also gets r[k].
        Normalized as r[0]^2 + 2 sum_{k>0} r[k]^2 = 1, with r[0] > 0.

        The taper is phi[n] = (r[0] + 2 sum_{k>0} r[k] cos(2 pi k (n - s) / N)) / sqrt(N) with
        s = (N - 1)/2. E.g. m = 1 (K = 0) gives r = [1], the rectangular window.
        """

    def taper(self) -> np.ndarray:
        n = 2**self.m_bits
        r = self.frequency_amplitudes()
        k = np.arange(1, len(r))
        cosines = np.cos(2 * np.pi * np.outer(np.arange(n) - (n - 1) / 2, k) / n)
        return (r[0] + 2 * cosines @ r[1:]) / np.sqrt(n)

    def frequency_state(self) -> np.ndarray:
        """The 2^m amplitudes prepared in step 1: r[|k|] at index k mod 2^m for k = -K..K, and 0
        at index 2^(m-1). E.g. m = 2: [r0, r1, 0, r1]."""
        m = self.num_confidence_qubits
        r = self.frequency_amplitudes()
        state = np.zeros(2**m)
        state[:len(r)] = r
        if len(r) > 1:
            state[-(len(r) - 1):] = r[:0:-1]
        return state

    def phase_ramp_exponents(self) -> Tuple[float, ...]:
        """ZPowGate exponents for step 2, one per m-register qubit (most significant first).

        Bit j of the m-bit value has weight c_j (-2^(m-1) for the sign bit, 2^j otherwise) and
        gets ZPowGate(-c_j (N - 1)/N), reduced to (-1, 1]. E.g. l = 1, m = 2: (-1/4, -7/8).
        """
        m, n = self.num_confidence_qubits, 2**self.m_bits
        weights = [-2**(m - 1)] + [2**j for j in reversed(range(m - 1))]
        exponents = []
        for c in weights:
            t = Fraction(-c * (n - 1), n) % 2
            exponents.append(float(t - 2 if t > 1 else t))
        return tuple(exponents)

    def build_composite_bloq(self, bb: BloqBuilder, qpe_reg: SoquetT) -> Dict[str, SoquetT]:
        l, m = self.num_precision_qubits, self.num_confidence_qubits
        qs = bb.split(qpe_reg)  # big-endian: qs[0] is the most significant bit
        low = bb.join(qs[l:], dtype=QAny(m))
        low = bb.add(ExactStatePreparation(self.frequency_state()), q=low)
        low_qs = bb.split(low)
        for j, exponent in enumerate(self.phase_ramp_exponents()):
            low_qs[j] = bb.add(ZPowGate(exponent=exponent), q=low_qs[j])
        qs[l:] = low_qs
        for i in range(l):
            qs[l], qs[i] = bb.add(CNOT(), ctrl=qs[l], target=qs[i])
        qpe_reg = bb.join(qs, dtype=self.m_qdtype)
        qpe_reg = bb.add(QFTTextBook(self.m_bits, with_reverse=True), q=qpe_reg)
        return {'qpe_reg': qpe_reg}

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        # From the circuit's structure, so the taper is not computed.
        counts = Counter(dict(state_preparation_call_graph(self.num_confidence_qubits, True, ssa)))
        counts.update(ZPowGate(exponent=t) for t in self.phase_ramp_exponents())
        counts[CNOT()] += self.num_precision_qubits
        counts[QFTTextBook(self.m_bits, with_reverse=True)] += 1
        return {(bloq, count) for bloq, count in counts.items() if count > 0}


@attrs.frozen
class TruncatedDPSSWindowState(BandLimitedWindowStateBase):
    """The DPSS taper with only its 2K + 1 lowest frequencies kept, renormalized (paper eq. 92).
    Failure probability at most 4x the DPSS's (eq. 107), for far fewer gates when l is large.

    Args and registers: see `TaperedWindowStateBase`.
    """

    def frequency_amplitudes(self) -> np.ndarray:
        n, k_max = 2**self.m_bits, self.half_bandwidth
        k = np.arange(k_max + 1)
        spectrum = np.fft.fft(dpss_taper(self.num_precision_qubits, self.num_confidence_qubits))
        # Undo the linear phase of a taper centred at (N - 1)/2; what remains is real.
        r = np.real(spectrum[:k_max + 1] * np.exp(1j * np.pi * k * (n - 1) / n))
        return r / np.sqrt(r[0]**2 + 2 * np.sum(r[1:]**2))


@attrs.frozen
class BandwidthOptimizedWindowState(BandLimitedWindowStateBase):
    """The taper on frequencies -K..K with the smallest average failure probability (paper
    Sec. VII B): at least as good as `TruncatedDPSSWindowState`, with the same circuit.

    Computed over symmetric tapers (the optimum is symmetric) as the smallest right singular
    vector of G, where (G v)^2 summed over quadrature points is the out-of-band energy
    int_{W < |f| < 1/2} |phihat(f)|^2 df of the taper with coefficients v. An SVD resolves this
    energy down to about 1e-30; an eigensolver for the in-band energy 1 - energy stops near 1e-16.

    Args and registers: see `TaperedWindowStateBase`.
    """

    def frequency_amplitudes(self) -> np.ndarray:
        n, k_max = 2**self.m_bits, self.half_bandwidth
        nodes, weights = np.polynomial.legendre.leggauss(8)
        # u = N f, on half-unit panels covering the out-of-band frequencies [K + 1/2, N/2].
        panel_starts = np.arange(2 * k_max + 1, n) / 2
        r_factor = np.zeros((0, k_max + 1))
        for chunk in np.array_split(panel_starts, max(1, len(panel_starts) // 4096)):
            u = (chunk[:, None] + (nodes + 1) / 4).ravel()
            g = np.sqrt(np.tile(weights, len(chunk)))[:, None] * self._spectra(u)
            r_factor = np.linalg.qr(np.vstack([r_factor, g]), mode='r')
        v = np.linalg.svd(r_factor)[2][-1]
        r = np.concatenate([v[:1], v[1:] / np.sqrt(2)])
        return r if r[0] > 0 else -r

    def _spectra(self, u: np.ndarray) -> np.ndarray:
        """|phihat(u/N)| (up to a common factor) of the orthonormal symmetric basis functions
        1/sqrt(N) and sqrt(2/N) cos(2 pi k (n - s)/N), one column per k = 0..K.

        Uses D(x) = sum_n cos(2 pi (n - s) x) = sin(pi N x)/sin(pi x), with sin(pi N x) computed
        from the fractional part of u so it stays accurate for large N.
        """
        n, k_max = 2**self.m_bits, self.half_bandwidth
        whole = np.floor(u)
        sin_pi_u = np.where(whole % 2, -1.0, 1.0) * np.sin(np.pi * (u - whole))

        def dirichlet(shift):  # D((u + shift)/N)
            x = (u + shift) / n
            return (-1)**shift * sin_pi_u / np.sin(np.pi * np.minimum(x, 1 - x))

        columns = [dirichlet(0)]
        columns += [(dirichlet(k) + dirichlet(-k)) / np.sqrt(2) for k in range(1, k_max + 1)]
        return np.array(columns).T
