"""
Trotterization Bloq for time evolution under non-commuting Hamiltonians.

This module implements Trotterization with ramped coefficients, expanding all steps and
ramps into a single flat sequence of Pauli string evolutions. Adjacent terms with the
same Pauli string are combined for efficiency.
"""

import logging
import numbers
from collections import Counter
from functools import cached_property
from math import cbrt
from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Set, Tuple, Union

import attrs
import numpy as np

from qualtran import Bloq, BloqBuilder, CtrlSpec, QBit, Register, Side, Signature, SoquetT
from qualtran.resource_counting import QubitCount

if TYPE_CHECKING:
    import quimb.tensor as qtn
    from qualtran import AddControlledT
    from qualtran.resource_counting import BloqCountT, CostKey, SympySymbolAllocator
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity

from qhat.common import controls
from qhat.common.commuting_pauli_string_evolution import CommutingPauliStringEvolution
from qhat.common.hermiticity import require_real_scalar
from qhat.common.pauli_utils import validate_pauli_string

# Get logger for this module
logger = logging.getLogger(__name__)


import math
import numpy as np

def Sn(n):
    if n == 2:
        return np.array((0.5, 0.5))
    else:
        # note: usual notation is 1/(k+1) but that uses k as the source power and we use n as the
        #       target power so k = n - 2 and k + 1 = n - 1
        a = 1 / (4 - math.pow(4, 1/(n-1)))
        b = 1 - 4 * a
        base = Sn(n - 2)
        outer = base * a
        inner = base * b
        return np.concatenate((outer, outer, inner, outer, outer))


def suzuki_fiveterm_recursion(n):
    assert isinstance(n, int)
    assert n >= 2
    assert n % 2 == 0
    return tuple(Sn(n))


def get_trotterization_coefficients(method):
    """Get coefficients for a named Trotterization method.

    Args:
        method: Either a string/int naming a method, or a sequence of coefficients.
            Supported names:
            - 1, "first order": First-order (1 cycle)
            - 2, "second order", "verlet", "leapfrog": Second-order Verlet/Leapfrog (1 cycle)
            - 3, "third order", "ruth 1983": Third-order Ruth (1983) (5 cycles)
            - "symmetrized ruth 1983": Fourth-order symmetrized Ruth (10 cycles)
            - 4, "fourth order", "suzuki5", "suzuki 5", "suzuki 1990", "suzuki 1990 5-term":
              Fourth-order Suzuki (1990) 5-term - DEFAULT (5 cycles)
            - "bm4", "blanes and moan 2002 fourth order", "blanes moan 4":
              Fourth-order Blanes & Moan (2002) - recommended by Ostmeyer (6 cycles)
            - "opt4", "ostmeyer 2023 fourth order", "optimised fourth order":
              Fourth-order Ostmeyer (2023) - high theoretical efficiency (5 cycles)
            - 6, "sixth order", "bm6", "blanes and moan 2002 sixth order", "blanes moan 6":
              Sixth-order Blanes & Moan (2002) (10 cycles)
            - 8, "eighth order", "morales", "morales 2022", "morales 2025":
              Eighth-order Morales et al. (2022) (17 cycles)

    Returns:
        Tuple of coefficients for the ramped Trotterization

    Raises:
        ValueError: If method is unknown

    Example:
        >>> coeffs = get_trotterization_coefficients("second order")
        >>> coeffs
        (0.5, 0.5)
    """
    if isinstance(method, str) or isinstance(method, int):
        m = method
        if isinstance(method, str):
            m = method.lower()

        if m in (1, "first order"):
            # First-order method (1 cycle)
            return (1.0,)
        elif m in (2, "second order", "verlet", "leapfrog"):
            # Second-order Verlet/Leapfrog method (1 cycle)
            return (0.5, 0.5)
        elif m in (3, "third order", "ruth 1983"):
            # Third-order Ruth (1983) (5 cycles)
            return (7./24., 3./8., 3./8., -25./24., 1.0)
        elif m == "symmetrized ruth 1983":
            # Symmetrized version of Ruth (1983) that raises it to fourth order (10 cycles)
            # Less efficient than Suzuki 1990 5-term
            return (7./48., 3./16., 3./16., -25./48., 0.5, 0.5, -25./48., 3./16., 3./16., 7./48.)
        elif m in (4, "fourth order", "suzuki5", "suzuki 5", "suzuki 1990", "suzuki 1990 5-term"):
            # Fourth-order Suzuki (1990) 5-term recursion (5 cycles) - DEFAULT
            # Apply the recursion relation from Suzuki (1990) to the standard second-order method
            # to get a fourth-order method; also discussed in Hatano and Suzuki (2005) and
            # Ostmeyer (2023). Ostmeyer 2023 recommends Suzuki 5-term recursion for
            # low-precision Hamiltonians with more than two terms in the Hamiltonian.
            return suzuki_fiveterm_recursion(4)
        elif m in ("bm4", "blanes and moan 2002 fourth order", "blanes moan 4"):
            # Fourth-order Blanes & Moan (2002) from equation 46 of Ostmeyer (2023) (6 cycles)
            # Ostmeyer 2023 recommends this as a reasonable default and notes it is optimal when
            # there are exactly two terms in the Hamiltonian. Has favorable error accumulation.
            # From paper equation 46:
            a1 =  0.07920369643119569
            b1 =  0.209515106613362
            a2 =  0.353172906049774
            b2 = -0.143851773179818
            a3 = -0.04206508035771955
            b3 =  0.5 - b1 - b2
            a4 =  1.0 - 2.0 * (a1 + a2 + a3)
            # By symmetry: a5 = a3, a6 = a2, a7 = a1 and b4 = b3, b5 = b2, b6 = b1

            # Convert (a,b) to (c,d) using equations 4-7:
            c1 = a1
            d1 = b1 - c1
            c2 = a2 - d1
            d2 = b2 - c2
            c3 = a3 - d2
            d3 = b3 - c3
            c4 = a4 - d3
            d4 = b3 - c4  # b4 = b3 by symmetry
            c5 = a3 - d4  # a5 = a3 by symmetry
            d5 = b2 - c5  # b5 = b2 by symmetry
            c6 = a2 - d5  # a6 = a2 by symmetry
            d6 = b1 - c6  # b6 = b1 by symmetry

            return (c1, d1, c2, d2, c3, d3, c4, d4, c5, d5, c6, d6)
        elif m in ("opt4", "ostmeyer 2023 fourth order", "optimised fourth order"):
            # Fourth-order Ostmeyer (2023) from equation 40 (5 cycles)
            # This is the theoretically most efficient q = 5 cycles fourth-order scheme.
            # However, empirically it has unfavorable error accumulation over time and is
            # usually outperformed by Blanes & Moan's scheme (46).
            # From paper equation 40:
            a1 =  0.09257547473195787
            b1 =  0.2540996315529392
            a2 =  0.4627160310210738
            b2 = -0.1676517240119692
            a3 =  0.5 - (a1 + a2)
            b3 =  1.0 - 2.0 * (b1 + b2)
            # By symmetry: a4 = a2, a5 = a1 and b4 = b2, b5 = b1

            # Convert (a,b) to (c,d) using equations 4-7:
            c1 = a1
            d1 = b1 - c1
            c2 = a2 - d1
            d2 = b2 - c2
            c3 = a3 - d2
            d3 = b3 - c3
            c4 = a2 - d3  # a4 = a2 by symmetry
            d4 = b2 - c4  # b4 = b2 by symmetry
            c5 = a1 - d4  # a5 = a1 by symmetry
            d5 = b1 - c5  # b5 = b1 by symmetry

            return (c1, d1, c2, d2, c3, d3, c4, d4, c5, d5)
        elif m in (6, "sixth order", "bm6", "blanes and moan 2002 sixth order", "blanes moan 6"):
            # Sixth-order Blanes & Moan (2002) from equation 53 of Ostmeyer (2023) (10 cycles)
            # Ostmeyer 2023 recommends this for high precision (relative error < 10^(-4)).
            # Most efficient known order n = 6 scheme.
            # From paper equation 53:
            a1 =  0.0502627644003922
            a2 =  0.413514300428344
            a3 =  0.04507988979439977
            a4 = -0.188054853819569
            a5 =  0.54196067845078
            a6 =  1.0 - 2.0 * (a1 + a2 + a3 + a4 + a5)

            b1 =  0.148816447901042
            b2 = -0.132385865767784
            b3 =  0.067307604692185
            b4 =  0.432666402578175
            b5 =  0.5 - (b1 + b2 + b3 + b4)
            # By symmetry: b6 = b5, b7 = b4, b8 = b3, b9 = b2, b10 = b1

            # Convert (a,b) to (c,d) using equations 4-7:
            c1 = a1
            d1 = b1 - c1
            c2 = a2 - d1
            d2 = b2 - c2
            c3 = a3 - d2
            d3 = b3 - c3
            c4 = a4 - d3
            d4 = b4 - c4
            c5 = a5 - d4
            d5 = b5 - c5
            c6 = a6 - d5
            d6 = b5 - c6  # b6 = b5 by symmetry
            c7 = a5 - d6  # a7 = a5 by symmetry
            d7 = b4 - c7  # b7 = b4 by symmetry
            c8 = a4 - d7  # a8 = a4 by symmetry
            d8 = b3 - c8  # b8 = b3 by symmetry
            c9 = a3 - d8  # a9 = a3 by symmetry
            d9 = b2 - c9  # b9 = b2 by symmetry
            c10 = a2 - d9  # a10 = a2 by symmetry
            d10 = b1 - c10  # b10 = b1 by symmetry

            return (c1, d1, c2, d2, c3, d3, c4, d4, c5, d5,
                    c6, d6, c7, d7, c8, d8, c9, d9, c10, d10)
        elif m in (8, "eighth order", "morales", "morales 2022", "morales 2025"):
            # Eighth-order Morales et al. (2022) from Ostmeyer (2023) (17 cycles)
            # Recommended by Ostmeyer (2023); paper modified and updated on arXiv as Morales et al. (2025)
            b1 =  0.12783360986284110837857554950443
            b2 =  0.56148845266356446893590729572808
            b3 = -0.38400573301491401473462588779099
            b4 =  0.15982762208609923217390166127256
            b5 = -0.40049110428180105319963667975074
            b6 =  0.18669648149540687549831902999911
            b7 =  0.26020394234904150277316667709864
            b8 =  0.29137384767986663096528500968049
            k = 0.5 - (b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8)
            return (b1/2, b1/2,
                    b2/2, b2/2,
                    b3/2, b3/2,
                    b4/2, b4/2,
                    b5/2, b5/2,
                    b6/2, b6/2,
                    b7/2, b7/2,
                    b8/2, b8/2,
                    k, k,
                    b8/2, b8/2,
                    b7/2, b7/2,
                    b6/2, b6/2,
                    b5/2, b5/2,
                    b4/2, b4/2,
                    b3/2, b3/2,
                    b2/2, b2/2,
                    b1/2, b1/2)
        elif m[:m.rstrip().rfind(' ')] in ("suzuki recursion", "suzuki 5-term recursion",
                                           "suzuki 1990 5-term recursion"):
            order = int(m[m.rstrip().rfind(' '):].rstrip())
            return suzuki_fiveterm_recursion(order)
        else:
            raise ValueError(f"Unknown Trotter method \"{method}\".")
    else:
        # If it's not a string naming a method, we assume it's a sequence of coefficients
        return tuple(method)


@attrs.frozen
class Trotterization(Bloq):
    """Time evolution using Trotterization with ramped coefficients.

    Implements the approximate time evolution operator U ≈ exp(-i H t / ℏ) for a
    Hamiltonian H = sum_i (h_i * P_i) where P_i are Pauli strings that may not commute.

    The Trotterization uses ramped coefficients (e.g., second-order, fourth-order methods)
    to achieve higher accuracy. All steps and ramps are expanded into a single flat
    sequence of operations. By default, adjacent identical terms are combined for efficiency.

    Attributes:
        pauli_terms: Sequence of (pauli_string, h_i) tuples defining the Hamiltonian.
            Example: [("XY", 0.5), ("ZI", 0.3), ("IZ", 0.2)]
        coefficients: Sequence of coefficients for the ramped method.
            Example: (0.5, 0.5) for second-order, (s2, s2, k, s2, s2) for fourth-order
        time: Total evolution time.
        num_steps: Number of Trotterization steps (time is divided by this).
        hbar: Reduced Planck constant (default 1.0).
        combine_terms: If True (default), combine adjacent identical terms. If False, keep all terms separate.
        tensor_contraction_method: Optional method to force for tensor contraction.
            None or "auto" = auto-select based on step count
            "incremental" = force O(n) incremental contraction
            "structured" = force O(log n) structured contraction (raises error if pattern not detected)
            "qualtran" = contract the decomposition gate by gate with Qualtran (reference)
        is_adjoint: If True, the operator is U^dagger. `adjoint()` flips it and leaves `time`,
            the coefficients and the term order as given; the decomposition and the tensor apply
            the reversal themselves.
        ctrl_values: Controls, as in `PauliStringEvolution`: one entry per control qubit, 1 (act
            on |1>) or 0 (act on |0>); empty for no control. `.controlled()` prepends one. Every
            term of the decomposition gets the same controls.
        multi_control_method: How 2 or more controls are built, "ancilla" (default) or "direct";
            see `common/controls.py`. With "ancilla" the controls are ANDed once for the whole
            Trotterization, and each term gets one control. Ignored with 0 or 1 controls.

    Example:
        >>> # Second-order Trotterization using named method
        >>> from trotterization import get_trotterization_coefficients
        >>> terms = [("XY", 0.5), ("YZ", 0.3)]
        >>> trotter = Trotterization.from_method(
        ...     pauli_terms=terms,
        ...     method="second order",
        ...     time=1.0,
        ...     num_steps=10
        ... )
        >>> tc = t_complexity(trotter)

        >>> # Or with explicit coefficients
        >>> trotter = Trotterization(
        ...     pauli_terms=terms,
        ...     coefficients=(0.5, 0.5),
        ...     time=1.0,
        ...     num_steps=10,
        ...     hbar=1.0
        ... )
    """

    pauli_terms: Tuple[Tuple[str, float], ...]
    coefficients: Tuple[float, ...]
    time: float
    num_steps: int
    hbar: float = 1.0
    combine_terms: bool = True
    _prologue: Tuple[Tuple[str, float], ...] = ()
    _epilogue: Tuple[Tuple[str, float], ...] = ()
    _repeat_core: Tuple[Tuple[str, float], ...] = ()
    _repeat_bridge: Tuple[Tuple[str, float], ...] = ()
    _symmetric_bookends: bool = True
    tensor_contraction_method: Optional[str] = None
    is_adjoint: bool = False
    ctrl_values: Tuple[int, ...] = ()
    multi_control_method: str = controls.DEFAULT_MULTI_CONTROL_METHOD

    def __attrs_post_init__(self):
        """Validate inputs."""
        if len(self.pauli_terms) == 0:
            raise ValueError("Must provide at least one Pauli term.")
        if self.num_steps <= 0:
            raise ValueError(f"num_steps must be positive, got {self.num_steps}")
        if len(self.coefficients) == 0:
            raise ValueError("Must provide at least one coefficient.")

        # Check all Pauli strings have the same length
        pauli_strings = [term[0] for term in self.pauli_terms]
        lengths = [len(ps) for ps in pauli_strings]
        if len(set(lengths)) > 1:
            raise ValueError(
                f"All Pauli strings must have the same length. Got lengths: {lengths}"
            )

        # Validate each term, and store coefficients, time, and hbar as floats
        pauli_terms = []
        for pauli_string, coeff in self.pauli_terms:
            validate_pauli_string(pauli_string)
            coeff = require_real_scalar(f"coefficient of {pauli_string}", coeff)
            pauli_terms.append((pauli_string, coeff))
        object.__setattr__(self, "pauli_terms", tuple(pauli_terms))
        object.__setattr__(self, "coefficients", tuple(
            require_real_scalar("Trotterization ramp coefficient", c) for c in self.coefficients))
        object.__setattr__(self, "time", require_real_scalar("Trotterization time", self.time))
        object.__setattr__(self, "hbar", require_real_scalar("Trotterization hbar", self.hbar))
        object.__setattr__(self, "ctrl_values",
                           controls.validate_ctrl_values("Trotterization", self.ctrl_values))
        controls.validate_multi_control_method("Trotterization", self.multi_control_method)

        # Initialize pattern basics
        object.__setattr__(self, "_prologue", ())
        object.__setattr__(self, "_epilogue", ())
        object.__setattr__(self, "_symmetric_bookends", True)
        object.__setattr__(self, "_repeat_bridge", ())

        # Initialize pattern core: a single step
        repeat_core = []
        ascending = True
        num_terms = len(self.pauli_terms)
        for coeff in self.coefficients:
            # Select indices for ramp direction
            if ascending:
                indices = range(num_terms) # 0, 1, 2, ..., n-1
            else:
                indices = range(num_terms - 1, -1, -1) # n-1, n-2, n-3, ..., 0
            # Append entire ramp
            for idx in indices:
                repeat_core.append((idx, coeff))
            # Alternate ramp direction
            ascending = not ascending
        object.__setattr__(self, "_repeat_core", tuple(repeat_core))

        # Combine terms
        if self.combine_terms:

            # Combine within a step
            combined = []
            current_idx, current_coeff = self._repeat_core[0]
            for idx, coeff in self._repeat_core[1:]:
                if idx == current_idx:
                    # Same term: combine
                    current_coeff += coeff
                else:
                    # Different term: save current and start new
                    combined.append((current_idx, current_coeff))
                    current_idx = idx
                    current_coeff = coeff
            # Don't forget the last term
            combined.append((current_idx, current_coeff))
            # Overwrite repeat core with combined version
            object.__setattr__(self, "_repeat_core", tuple(combined))

            # Combine across steps
            if len(self._repeat_core) == 1:
                # The point of taking multiple Trotter steps is to lower the error that comes from
                # the Trotterization process.  That error arises from having multiple terms in the
                # Hamiltonian that do not commute with each other.  If you only have a single term
                # in your Hamiltonian, then no error can arise from non-commutation.  This is a
                # silly application of multi-step Trotterization, so something odd is already
                # happening.  But the best resolution is to combine all the steps together.
                # The coefficient is unchanged: with a single step, dt = time / num_steps grows
                # by the factor num_steps instead.
                object.__setattr__(self, "num_steps", 1)
            else:
                # Two terms will never combine across steps, because they would have already been
                # combined within a step.  At this point you have at least three terms, or you have
                # two terms that are different and therefore cannot combine across steps.
                head_idx, head_coeff = self._repeat_core[0]
                tail_idx, tail_coeff = self._repeat_core[-1]
                if head_idx == tail_idx:
                    merged_idx = head_idx
                    merged_coeff = head_coeff + tail_coeff
                    object.__setattr__(self, "_repeat_core", self._repeat_core[1:-1])
                    object.__setattr__(self, "_repeat_bridge", ((merged_idx, merged_coeff),))
                    object.__setattr__(self, "_prologue", ((head_idx, head_coeff),))
                    object.__setattr__(self, "_epilogue", ((tail_idx, tail_coeff),))
                    object.__setattr__(self, "_symmetric_bookends", head_coeff == tail_coeff)

        def trm_str(sequence):
            n = len(sequence)
            return f"{n} {'term' if n == 1 else 'terms'}"
        def rpt_str(n):
            return f"{n} {'repetition' if n == 1 else 'repetitions'}"
        logger.debug("\n".join(["",
            "Trotter reptition pattern:",
            f"-- prologue ({trm_str(self._prologue)})",
            f"-- repeat core ({trm_str(self._repeat_core)} x {rpt_str(self.num_steps)})",
            f"-- repeat bridge ({trm_str(self._repeat_bridge)} x {rpt_str(self.num_steps-1)})",
            f"-- epilogue ({trm_str(self._epilogue)})",
            ]))


    @classmethod
    def from_method(
        cls,
        pauli_terms: Sequence[Tuple[str, float]],
        method,
        time: float,
        num_steps: int,
        hbar: float = 1.0,
        combine_terms: bool = True,
        tensor_contraction_method: Optional[str] = None,
        multi_control_method: str = controls.DEFAULT_MULTI_CONTROL_METHOD,
    ):
        """Create Trotterization using a named method.

        Args:
            pauli_terms: Sequence of (pauli_string, h_i) tuples
            method: Method name (e.g., "second order", "fourth order") or coefficients
            time: Total evolution time
            num_steps: Number of Trotterization steps
            hbar: Reduced Planck constant (default 1.0)
            combine_terms: If True (default), combine adjacent identical terms. If False, keep all terms separate.
            tensor_contraction_method: Optional method to force for tensor contraction.
                None or "auto" = auto-select based on step count
                "incremental" = force O(n) incremental contraction
                "structured" = force O(log n) structured contraction (raises error if pattern not detected)
            multi_control_method: How 2 or more controls are built, "ancilla" or "direct".

        Returns:
            Trotterization instance

        Example:
            >>> trotter = Trotterization.from_method(
            ...     pauli_terms=[("X", 1.0), ("Y", 1.0)],
            ...     method="second order",
            ...     time=1.0,
            ...     num_steps=10
            ... )
        """
        coefficients = get_trotterization_coefficients(method)
        return cls(
            pauli_terms=tuple(pauli_terms),
            coefficients=coefficients,
            time=time,
            num_steps=num_steps,
            hbar=hbar,
            combine_terms=combine_terms,
            tensor_contraction_method=tensor_contraction_method,
            multi_control_method=multi_control_method,
        )

    @cached_property
    def num_qubits(self) -> int:
        """Number of qubits the operator acts on."""
        return len(self.pauli_terms[0][0])

    @cached_property
    def signature(self) -> Signature:
        """The `ctrl` register (if controlled), then a register of qubits."""
        return Signature(controls.ctrl_registers(self.ctrl_values)
                         + [Register('q', QBit(), shape=(self.num_qubits,), side=Side.THRU)])

    @cached_property
    def expanded_sequence(self) -> List[Tuple[int, float]]:
        """Expanded sequence of (term_index, coefficient) pairs.

        This expands all steps and ramps into a flat list. If combine_terms is True,
        adjacent identical terms are combined for efficiency.
        """
        sequence = list()
        sequence.extend(self._prologue)
        for _ in range(self.num_steps - 1):
            sequence.extend(self._repeat_core)
            sequence.extend(self._repeat_bridge)
        sequence.extend(self._repeat_core)
        sequence.extend(self._epilogue)
        return sequence

    def build_composite_bloq(self, bb: BloqBuilder, **soqs: SoquetT) -> Dict[str, SoquetT]:
        """Decompose into a sequence of CommutingPauliStringEvolution bloqs.

        Uses CommutingPauliStringEvolution to enable future grouping of commuting terms.
        Currently, each instance contains a single term, but the infrastructure is in place
        for future optimization where multiple commuting terms can be grouped together.

        Uses the expanded sequence with combined terms for efficiency. The adjoint applies the
        sequence in reverse order, with each term adjointed. Under control every term gets the
        same controls. With 2 or more controls and the "ancilla" method, the controls are first
        ANDed into one ancilla, which controls a singly-controlled Trotterization.
        """
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.add_with_and_ancilla(bb, self, **soqs)

        # Apply each term in the expanded sequence (reversed for the adjoint)
        sequence = reversed(self.expanded_sequence) if self.is_adjoint else self.expanded_sequence
        for term_idx, coeff in sequence:
            soqs = bb.add_d(self._term_evolution(term_idx, coeff), **soqs)

        return soqs

    def _term_evolution(self, term_idx: int, coeff: float) -> CommutingPauliStringEvolution:
        """The bloq the decomposition applies for one (term_index, coefficient) pair.

        The effective evolution is exp(-i * h_i * P_i * coeff * dt / hbar), with the adjoint flag
        and controls of this Trotterization. Currently wraps a single term, but enables future
        grouping of commuting terms. `build_call_graph` uses the same bloqs, so its callees equal
        the decomposition's.
        """
        pauli_string, h_i = self.pauli_terms[term_idx]
        return CommutingPauliStringEvolution(
            pauli_terms=((pauli_string, h_i * coeff),),
            time=self.time / self.num_steps,
            hbar=self.hbar,
            is_adjoint=self.is_adjoint,
            ctrl_values=self.ctrl_values,
            multi_control_method=self.multi_control_method,
        )

    def _occurrences(self) -> Counter:
        """How often each (term_index, coefficient) pair appears in `expanded_sequence`, without
        expanding it: the prologue and epilogue once, the core num_steps times and the bridge
        num_steps - 1 times.

        Example: fourth order with 2 terms and 3 steps has 6 distinct pairs, whatever the steps.
        """
        counts = Counter()
        for segment, repeats in ((self._prologue, 1),
                                 (self._repeat_core, self.num_steps),
                                 (self._repeat_bridge, self.num_steps - 1),
                                 (self._epilogue, 1)):
            if repeats == 0:
                continue  # the bridge of a 1-step Trotterization; a 0 count would be a callee
            for pair in segment:
                counts[pair] += repeats
        return counts

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        """The decomposition's callees and their counts, without building the decomposition.

        Qualtran's default decomposes the whole circuit, which grows with num_steps (e.g. 25 s
        for `QECGatesCost` of a 500-step QPE). Here the work grows only with the length of the
        repeat pattern. The callees are exactly the decomposition's bloqs: the same coefficients,
        time step, adjoint flag and controls. With the "ancilla" method and 2 or more controls,
        they are the AND, its uncompute, and this bloq with one control.
        """
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.and_ancilla_call_graph(self)
        counts = Counter()
        for (term_idx, coeff), n in self._occurrences().items():
            counts[self._term_evolution(term_idx, coeff)] += n
        return set(counts.items())

    def _t_complexity_(self) -> TComplexity:
        """Return T-complexity by counting term occurrences.

        Since T-complexity is independent of time and coefficients, we count how many
        times each term appears and compute the complexity once per unique term.

        Uses CommutingPauliStringEvolution to match the decomposition used in
        build_composite_bloq, with the same controls. The "ancilla" method adds the AND ladder
        once and counts the terms with one control.
        """
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.and_ancilla_t_complexity(self)

        # Count occurrences of each term in the expanded sequence
        term_counts = Counter()
        for (term_idx, coeff), n in self._occurrences().items():
            term_counts[term_idx] += n

        # Compute complexity once per unique term and multiply by occurrence count
        total_complexity = TComplexity()

        for term_idx, count in term_counts.items():
            pauli_string, h_i = self.pauli_terms[term_idx]
            # Coefficient and time don't affect T-complexity, so use dummy values
            # Use CommutingPauliStringEvolution to match build_composite_bloq
            cpse = CommutingPauliStringEvolution(
                pauli_terms=((pauli_string, 1.0),),  # Dummy coefficient
                time=1.0,         # Dummy value - doesn't affect complexity
                hbar=1.0,
                ctrl_values=self.ctrl_values,
                multi_control_method=self.multi_control_method,
            )
            term_complexity = t_complexity(cpse)
            total_complexity += count * term_complexity

        return total_complexity

    def add_my_tensors(
        self,
        tn: 'qtn.TensorNetwork',
        tag: object,
        *,
        incoming: Dict[str, SoquetT],
        outgoing: Dict[str, SoquetT],
    ):
        """Add `tensor_contract()` as a single dense tensor, so that a Trotterization inside a
        larger bloq (e.g. `Power` or a user-built `CompositeBloq`) also uses the optimized
        contraction instead of being decomposed."""
        import quimb.tensor as qtn

        n = self.num_qubits
        num_ctrls = len(self.ctrl_values)
        data = self.tensor_contract().reshape([2] * (2 * (num_ctrls + n)))
        in_idx = ([incoming['ctrl'][i] for i in range(num_ctrls)]
                  + [incoming['q'][i] for i in range(n)])
        out_idx = ([outgoing['ctrl'][i] for i in range(num_ctrls)]
                   + [outgoing['q'][i] for i in range(n)])
        tn.add(qtn.Tensor(data=data, inds=tuple(out_idx + in_idx), tags=[self.pretty_name(), tag]))

    def tensor_contract(self) -> np.ndarray:
        """
        Optimized tensor contraction

        The tensor_contraction_method attribute can force a specific method:
        - None or "auto": use auto-selection logic (default behavior)
        - "incremental": force O(n) incremental contraction
        - "structured": force O(log n) structured contraction (raises error if pattern not detected)
        - "qualtran": contract the decomposition gate by gate with Qualtran

        Under control, the result is block-diagonal: U (or U^dagger) on the block matching
        `ctrl_values` and the identity elsewhere, with the controls as the most significant qubits.

        Returns:
            The full unitary matrix as a numpy array with shape (2^n, 2^n), where n is the number
            of qubits plus the number of controls

        Example:
            >>> trotter = Trotterization.from_method(
            ...     pauli_terms=[("X", 1.0), ("Y", 1.0)],
            ...     method="second order",
            ...     time=1.0,
            ...     num_steps=1000
            ... )
            >>> U = trotter.tensor_contract()  # Uses O(log n) optimization
            >>> U.shape
            (2, 2)
        """
        logger.verbose(f"Computing unitary matrix via tensor contraction...")
        logger.debug(f"-- num_qubits = {self.num_qubits}")
        logger.debug(f"-- num_steps = {self.num_steps}")
        logger.debug(f"-- num_terms = {len(self.pauli_terms)}")
        logger.debug(f"-- expanded_sequence length = {self.num_operations}")

        # The "qualtran" method contracts the decomposition, which already applies the adjoint
        # flag and the controls.
        if self.tensor_contraction_method == "qualtran":
            logger.verbose(f"Using Qualtran's gate-by-gate tensor contraction")
            # Not super().tensor_contract(): that would route back here via add_my_tensors.
            return self.decompose_bloq().tensor_contract()

        # The other methods build the forward matrix, so conjugate-transpose it for the adjoint
        # and then add the controls.
        unitary = self._forward_matrix()
        if self.is_adjoint:
            unitary = unitary.conj().T
        return controls.controlled_matrix(unitary, self.ctrl_values)

    def _forward_matrix(self) -> np.ndarray:
        """The uncontrolled, forward matrix U, by the configured contraction method (not "qualtran")."""
        forced_method = self.tensor_contraction_method
        if forced_method is not None and forced_method != "auto":
            if forced_method == "incremental":
                logger.verbose(f"Using O(n) incremental tensor contraction")
                return self._incremental_contraction()
            elif forced_method == "structured":
                logger.verbose(f"Using O(log n) structured tensor contraction")
                return self._structured_contraction()
            else:
                raise ValueError(
                    f"Invalid tensor_contraction_method '{forced_method}'. "
                    f"Must be None, 'auto', 'incremental', 'structured', or 'qualtran'."
                )

        # Default (non-forced): use structured contraction
        logger.verbose(f"Using default (O(log n) structured) tensor contraction")
        return self._structured_contraction()

    def _structured_contraction(self) -> np.ndarray:
        """O(log n) contraction of the repeating pattern.

        Returns the forward, uncontrolled matrix U: `tensor_contract` applies the adjoint and the
        controls.
        """
        def trm_str(sequence):
            n = len(sequence)
            return f"{n} {'term' if n == 1 else 'terms'}"
        def tms_str(n):
            return f"{n} {'time' if n == 1 else 'times'}"
        logger.debug("Performing tensor contraction with structured method.")
        logger.debug(f"-- prologue contains {trm_str(self._prologue)}")
        logger.debug("-- repeating pattern: core bridge core bridge ... core bridge core")
        logger.debug(f"   -- core contains {trm_str(self._repeat_core)} "
                     f"and repeats {tms_str(self.num_steps)}")
        logger.debug(f"   -- bridge contains {trm_str(self._repeat_bridge)} "
                     f"and repeats {tms_str(self.num_steps-1)}")
        logger.debug(f"-- epilogue contains {trm_str(self._epilogue)}")
        # Repeat section
        result = self._build_matrix_from_terms(self._repeat_core)
        if self.num_steps > 1:
            if self._repeat_bridge:
                temp = self._build_matrix_from_terms(self._repeat_bridge)
                temp = temp @ result
                if self.num_steps > 2:
                    temp = np.linalg.matrix_power(temp, self.num_steps - 1)
                result = result @ temp
            else:
                result = np.linalg.matrix_power(result, self.num_steps)
        # Bookends
        if self._symmetric_bookends:
            if self._prologue:
                temp = self._build_matrix_from_terms(self._prologue)
                result = temp @ result @ temp
        else:
            if self._prologue:
                temp = self._build_matrix_from_terms(self._prologue)
                result = result @ temp
            if self._epilogue:
                temp = self._build_matrix_from_terms(self._epilogue)
                result = temp @ result
        return result

    def _build_matrix_from_terms(self, terms: List[Tuple[int, float]]) -> np.ndarray:
        """
        Build unitary matrix from a sequence of (term_index, coefficient) pairs.

        Args:
            terms: List of (term_idx, coeff) tuples from expanded_sequence

        Returns:
            Unitary matrix for this sequence of terms
        """
        if len(terms) == 0:
            logger.debug(f"    Building identity matrix (empty term list)")
            dim = 2 ** self.num_qubits
            return np.eye(dim, dtype=np.complex128)

        logger.debug(f"Building matrix from {len(terms)} terms")
        dim = 2 ** self.num_qubits
        U = np.eye(dim, dtype=np.complex128)

        dt = self.time / self.num_steps

        for term_idx, coeff in terms:
            pauli_string, h_i = self.pauli_terms[term_idx]

            cpse = CommutingPauliStringEvolution(
                pauli_terms=((pauli_string, h_i * coeff),),
                time=dt,
                hbar=self.hbar
            )

            U_term = cpse.tensor_contract()
            U = U_term @ U

        return U

    def _incremental_contraction(self) -> np.ndarray:
        """
        O(n) incremental contraction

        Iterates through expanded_sequence and multiplies matrices one at a time.
        Used for small num_steps.

        Returns:
            The forward, uncontrolled matrix U: `tensor_contract` applies the adjoint and the
            controls.
        """
        logger.debug(f"Performing incremental contraction over {self.num_operations} operations")
        return self._build_matrix_from_terms(self.expanded_sequence)

    @property
    def num_terms(self) -> int:
        """Number of Pauli string terms in the Hamiltonian."""
        return len(self.pauli_terms)

    @property
    def num_operations(self) -> int:
        """Number of operations in the expanded sequence (after combining)."""
        return len(self._prologue) \
                + len(self._repeat_core) * self.num_steps \
                + len(self._repeat_bridge) * (self.num_steps - 1) \
                + len(self._epilogue)

    def __str__(self) -> str:
        method_name = self._infer_method_name()
        dagger = "†" if self.is_adjoint else ""
        return controls.ctrl_str(self.ctrl_values,
                                 f"Trotterization({method_name}, {self.num_terms} terms, "
                                 f"{self.num_steps} steps, t={self.time}){dagger}")

    def __repr__(self) -> str:
        adjoint = ", is_adjoint=True" if self.is_adjoint else ""
        ctrl = controls.ctrl_repr(self.ctrl_values, self.multi_control_method)
        return (f"Trotterization(pauli_terms={self.pauli_terms}, "
                f"coefficients={self.coefficients}, time={self.time}, "
                f"num_steps={self.num_steps}, hbar={self.hbar}{adjoint}{ctrl})")

    def _infer_method_name(self) -> str:
        """Infer the method name from coefficients for display."""
        if self.coefficients == (1.0,):
            return "first-order"
        elif self.coefficients == (0.5, 0.5):
            return "second-order"
        else:
            return f"{len(self.coefficients)}-coeff"

    def phase_error_from_energy_error(self, energy_error: float) -> float:
        """Eigenphase error, in turns ([0, 1)), corresponding to an energy error.

        U = exp(-i H t / hbar), so an energy E has phase E t / (2 pi hbar) turns. This is a
        magnitude: it uses |time|, so a negative time or an adjoint gives the same (positive) error.
        """
        return energy_error * abs(self.time) / (2 * np.pi * self.hbar)

    def __pow__(self, power: int):
        """Raise the Trotterization to a non-zero integer power.

        For QPE, we need U^k where k is an integer. Since U applies the Trotter
        sequence num_steps times with time step time / num_steps, U^k applies the same
        sequence (num_steps * k) times with the same time step, i.e. for total time (time * k).
        U^(-k) is (U^k).adjoint(). A power of an adjoint keeps the flag: (U^dagger)^k is
        (U^k)^dagger.

        Args:
            power: Non-zero integer power to raise to

        Returns:
            New Trotterization bloq with num_steps and time multiplied by |power|

        Raises:
            ValueError: If power is zero or not an integer
        """
        if isinstance(power, bool) or not isinstance(power, numbers.Integral) or power == 0:
            raise ValueError(f"Power must be a non-zero integer, got {power}")

        k = abs(int(power))
        return attrs.evolve(
            self, num_steps=self.num_steps * k, time=self.time * k,
            is_adjoint=self.is_adjoint != (power < 0))

    def adjoint(self) -> 'Trotterization':
        """Return U^dagger as a Trotterization, so it keeps the optimized tensor contraction.

        Flips `is_adjoint` and leaves `time`, the coefficients and the term order as given, so
        the adjoint of the adjoint equals the original. Negating `time` alone would be wrong for
        methods that are not symmetric; the decomposition and the tensor apply the reversal.
        """
        return attrs.evolve(self, is_adjoint=not self.is_adjoint)

    def get_ctrl_system(
        self, ctrl_spec: Optional[CtrlSpec] = None
    ) -> Tuple['Trotterization', 'AddControlledT']:
        """Controlled version: the same bloq with the new control value prepended to `ctrl_values`.

        Only single-qubit controls on |1> or |0> are supported; others raise ValueError.
        """
        return controls.get_ctrl_system(self, ctrl_spec)

    def my_static_costs(self, cost_key: 'CostKey'):
        """`QubitCount` without decomposing: the qubits, the controls, and the k-1 AND ancillas of
        the "ancilla" method with k >= 2 controls.

        Nothing else in the decomposition allocates qubits, so this is the peak number in use.
        """
        if isinstance(cost_key, QubitCount):
            return controls.peak_qubits(self.num_qubits, self.ctrl_values, self.multi_control_method)
        return NotImplemented


# =================================================================================================
# Compatibility function for old interface
# =================================================================================================

def build_ramped_trotterized_unitary(
    pauli_strings: Union[Sequence[Tuple[str, float]], 'Iterable[Tuple[str, float]]'],
    method,
    timestep: float,
    numsteps: int,
    combine_terms: bool = True,
    tensor_contraction_method: Optional[str] = None,
    multi_control_method: str = controls.DEFAULT_MULTI_CONTROL_METHOD,
):
    """Build a Trotterization using the old interface for backward compatibility.

    This function provides compatibility with the old `common.trotter` module interface
    used by QPE code in the analysis/ directory.

    Args:
        pauli_strings: Iterable of (pauli_string, coefficient) tuples.
            Example: [("XYZ", 0.5), ("ZZI", 0.3)]
        method: Trotterization method name or int (e.g., "second order", 2)
        timestep: Time step for evolution
        numsteps: Number of Trotterization steps
        combine_terms: If True (default), combine adjacent identical terms. If False, keep all terms separate.
        tensor_contraction_method: Optional method to force for tensor contraction.
            None or "auto" = auto-select based on step count
            "incremental" = force O(n) incremental contraction
            "structured" = force O(log n) structured contraction (raises error if pattern not detected)
            "qualtran" = contract the decomposition gate by gate with Qualtran (reference)
        multi_control_method: How 2 or more controls are built, "ancilla" or "direct".

    Returns:
        Trotterization instance

    Example:
        >>> # Old interface (for compatibility)
        >>> from common.trotter_flattened import build_ramped_trotterized_unitary
        >>> bloq = build_ramped_trotterized_unitary(
        ...     [("XY", 0.5), ("YZ", 0.3)],
        ...     "second order",
        ...     1.0,
        ...     10
        ... )
    """
    return Trotterization.from_method(
        pauli_terms=pauli_strings,
        method=method,
        time=timestep,
        num_steps=numsteps,
        hbar=1.0,
        combine_terms=combine_terms,
        tensor_contraction_method=tensor_contraction_method,
        multi_control_method=multi_control_method,
    )
