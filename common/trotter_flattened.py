"""
Trotterization Bloq for time evolution under non-commuting Hamiltonians.

This module implements Trotterization with ramped coefficients, expanding all steps and
ramps into a single flat sequence of Pauli string evolutions. Adjacent terms with the
same Pauli string are combined for efficiency.
"""

import logging
from functools import cached_property
from math import cbrt
from typing import List, Optional, Sequence, Tuple, Union

import attrs
import numpy as np

from qualtran import Bloq, BloqBuilder, Signature, SoquetT
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity

from qhat.common.commuting_pauli_string_evolution import CommutingPauliStringEvolution
from qhat.common.pauli_utils import validate_pauli_string

# Get logger for this module
logger = logging.getLogger(__name__)


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
            s2 = 0.5 / (4.0 - cbrt(4.0))
            k = 0.5 - 4.0 * s2
            return (s2, s2, s2, s2, k, k, s2, s2, s2, s2)
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
            "qualtran" = use Qualtran's inherited Bloq.tensor_contract() method

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

        # Validate each Pauli string
        for pauli_string, coeff in self.pauli_terms:
            validate_pauli_string(pauli_string)

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

        logger.debug("Basic Pattern:")
        logger.debug(f"-- prologue ({len(self._prologue)} terms)")
        for idx, coeff in self._prologue:
            logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
        logger.debug(f"-- repeat core ({len(self._repeat_core)} terms x {self.num_steps} repetitions)")
        for idx, coeff in self._repeat_core:
            logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
        logger.debug(f"-- repeat bridge ({len(self._repeat_bridge)} terms x {self.num_steps-1} repetitions)")
        for idx, coeff in self._repeat_bridge:
            logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
        logger.debug(f"-- epilogue ({len(self._epilogue)} terms)")
        for idx, coeff in self._epilogue:
            logger.debug(f"   -- {coeff:13.6e} {idx:>9}")

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

            logger.debug("Deduplicated Pattern:")
            logger.debug(f"-- prologue ({len(self._prologue)} terms)")
            for idx, coeff in self._prologue:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
            logger.debug(f"-- repeat core ({len(self._repeat_core)} terms x {self.num_steps} repetitions)")
            for idx, coeff in self._repeat_core:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
            logger.debug(f"-- repeat bridge ({len(self._repeat_bridge)} terms x {self.num_steps-1} repetitions)")
            for idx, coeff in self._repeat_bridge:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
            logger.debug(f"-- epilogue ({len(self._epilogue)} terms)")
            for idx, coeff in self._epilogue:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")

            # Combine across steps
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

            logger.debug("Restructured Pattern:")
            logger.debug(f"-- prologue ({len(self._prologue)} terms)")
            for idx, coeff in self._prologue:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
            logger.debug(f"-- repeat core ({len(self._repeat_core)} terms x {self.num_steps} repetitions)")
            for idx, coeff in self._repeat_core:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
            logger.debug(f"-- repeat bridge ({len(self._repeat_bridge)} terms x {self.num_steps-1} repetitions)")
            for idx, coeff in self._repeat_bridge:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")
            logger.debug(f"-- epilogue ({len(self._epilogue)} terms)")
            for idx, coeff in self._epilogue:
                logger.debug(f"   -- {coeff:13.6e} {idx:>9}")


    @classmethod
    def from_method(
        cls,
        pauli_terms: Sequence[Tuple[str, float]],
        method,
        time: float,
        num_steps: int,
        hbar: float = 1.0,
        combine_terms: bool = True,
        tensor_contraction_method: Optional[str] = None
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
            tensor_contraction_method=tensor_contraction_method
        )

    @cached_property
    def num_qubits(self) -> int:
        """Number of qubits the operator acts on."""
        return len(self.pauli_terms[0][0])

    @cached_property
    def signature(self) -> Signature:
        """Return the signature with a register of qubits."""
        from qualtran import Register, QBit, Side
        return Signature([Register('q', QBit(), shape=(self.num_qubits,), side=Side.THRU)])

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
        logger.debug("New Full Sequence ({len(sequence)} terms)")
        for idx, coeff in sequence:
            logger.debug(f"-- {coeff:13.6e} {idx:>9}")
        return sequence

    def build_composite_bloq(self, bb: BloqBuilder, **soqs: SoquetT) -> dict[str, SoquetT]:
        """Decompose into a sequence of CommutingPauliStringEvolution bloqs.

        Uses CommutingPauliStringEvolution to enable future grouping of commuting terms.
        Currently, each instance contains a single term, but the infrastructure is in place
        for future optimization where multiple commuting terms can be grouped together.

        Uses the expanded sequence with combined terms for efficiency.
        """
        # Get the input register
        assert len(soqs) == 1, f"Expected 1 register, got {len(soqs)}"
        reg_name, register = list(soqs.items())[0]

        # Time step for each Trotterization step
        dt = self.time / self.num_steps

        # Apply each term in the expanded sequence
        for term_idx, coeff in self.expanded_sequence:
            pauli_string, h_i = self.pauli_terms[term_idx]

            # Create CommutingPauliStringEvolution for this term
            # The effective evolution is exp(-i * h_i * P_i * coeff * dt / hbar)
            # Currently wraps a single term, but enables future grouping of commuting terms
            cpse = CommutingPauliStringEvolution(
                pauli_terms=((pauli_string, h_i * coeff),),
                time=dt,
                hbar=self.hbar
            )
            register = bb.add(cpse, q=register)

        return {reg_name: register}

    def _t_complexity_(self) -> TComplexity:
        """Return T-complexity by counting term occurrences.

        Since T-complexity is independent of time and coefficients, we count how many
        times each term appears and compute the complexity once per unique term.

        Uses CommutingPauliStringEvolution to match the decomposition used in
        build_composite_bloq.
        """
        # Count occurrences of each term in the expanded sequence
        term_counts = {}
        # -- prologue happens once
        for term_idx, coeff in self._prologue:
            term_counts[term_idx] = term_counts.get(term_idx, 0) + 1
        # -- repeat_core happens num_steps times
        for term_idx, coeff in self._repeat_core:
            term_counts[term_idx] = term_counts.get(term_idx, 0) + self.num_steps
        # -- repeat_bridge happens (num_steps-1) times
        for term_idx, coeff in self._repeat_bridge:
            term_counts[term_idx] = term_counts.get(term_idx, 0) + (self.num_steps - 1)
        # -- epilogue happens once
        for term_idx, coeff in self._epilogue:
            term_counts[term_idx] = term_counts.get(term_idx, 0) + 1

        # Compute complexity once per unique term and multiply by occurrence count
        dt = self.time / self.num_steps
        total_complexity = TComplexity()

        for term_idx, count in term_counts.items():
            pauli_string, h_i = self.pauli_terms[term_idx]
            # Coefficient and time don't affect T-complexity, so use dummy values
            # Use CommutingPauliStringEvolution to match build_composite_bloq
            cpse = CommutingPauliStringEvolution(
                pauli_terms=((pauli_string, 1.0),),  # Dummy coefficient
                time=1.0,         # Dummy value - doesn't affect complexity
                hbar=1.0
            )
            term_complexity = t_complexity(cpse)
            total_complexity += count * term_complexity

        return total_complexity

    def tensor_contract(self) -> np.ndarray:
        """
        Optimized tensor contraction

        The tensor_contraction_method attribute can force a specific method:
        - None or "auto": use auto-selection logic (default behavior)
        - "incremental": force O(n) incremental contraction
        - "structured": force O(log n) structured contraction (raises error if pattern not detected)
        - "qualtran": use Qualtran's inherited Bloq.tensor_contract() method

        Returns:
            The full unitary matrix as a numpy array with shape (2^n_qubits, 2^n_qubits)

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

        # Check if a specific method is forced via configuration
        forced_method = self.tensor_contraction_method
        if forced_method is not None and forced_method != "auto":
            if forced_method == "qualtran":
                logger.verbose(f"Using Qualtran's Bloq.tensor_contract() method")
                return super().tensor_contract()
            elif forced_method == "incremental":
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
        logger.debug("Performing tensor contraction with structured method.")
        logger.debug(f"-- prologue contains {len(self._prologue)} terms")
        logger.debug("-- repeating pattern: core bridge core bridge ... core bridge core")
        logger.debug(f"   -- core contains {len(self._repeat_core)} terms and repeats {self.num_steps} times")
        logger.debug(f"   -- bridge contains {len(self._repeat_bridge)} terms and repeats {self.num_steps-1} times")
        logger.debug(f"-- epilogue contains {len(self._epilogue)} terms")
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

        logger.debug(f"    Building matrix from {len(terms)} terms")
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
            Full unitary matrix
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
        return (f"Trotterization({method_name}, {self.num_terms} terms, "
                f"{self.num_steps} steps, t={self.time})")

    def __repr__(self) -> str:
        return (f"Trotterization(pauli_terms={self.pauli_terms}, "
                f"coefficients={self.coefficients}, time={self.time}, "
                f"num_steps={self.num_steps}, hbar={self.hbar})")

    def _infer_method_name(self) -> str:
        """Infer the method name from coefficients for display."""
        if self.coefficients == (1.0,):
            return "first-order"
        elif self.coefficients == (0.5, 0.5):
            return "second-order"
        else:
            return f"{len(self.coefficients)}-coeff"

    def __pow__(self, power: int):
        """Raise the Trotterization to an integer power.

        For QPE, we need U^k where k is an integer. Since U applies the Trotter
        sequence num_steps times, U^k applies it (num_steps * k) times.

        Args:
            power: Integer power to raise to (must be positive)

        Returns:
            New Trotterization bloq with num_steps multiplied by power

        Raises:
            ValueError: If power is not a positive integer
        """
        if not isinstance(power, int) or power < 1:
            raise ValueError(f"Power must be a positive integer, got {power}")

        return attrs.evolve(self, num_steps=self.num_steps * power)


# =================================================================================================
# Compatibility function for old interface
# =================================================================================================

def build_ramped_trotterized_unitary(
    pauli_strings: Union[Sequence[Tuple[str, float]], 'Iterable[Tuple[str, float]]'],
    method,
    timestep: float,
    numsteps: int,
    combine_terms: bool = True,
    tensor_contraction_method: Optional[str] = None
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
            "qualtran" = use Qualtran's inherited Bloq.tensor_contract() method

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
        tensor_contraction_method=tensor_contraction_method
    )
