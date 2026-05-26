"""
Trotterization Bloq for time evolution under non-commuting Hamiltonians.

This module implements Trotterization with ramped coefficients, expanding all steps and
ramps into a single flat sequence of Pauli string evolutions. Adjacent terms with the
same Pauli string are combined for efficiency.
"""

from functools import cached_property
from math import cbrt
from typing import List, Sequence, Tuple, Union

import attrs
import numpy as np

from qualtran import Bloq, BloqBuilder, Signature, SoquetT
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity

from qhat.common.commuting_pauli_string_evolution import CommutingPauliStringEvolution
from qhat.common.pauli_utils import validate_pauli_string


def get_trotterization_coefficients(method):
    """Get coefficients for a named Trotterization method.

    Args:
        method: Either a string/int naming a method, or a sequence of coefficients.
            Supported names:
            - 1, "first order": Minimal first-order method (1.0,)
            - 2, "second order": Minimal second-order method (0.5, 0.5)
            - 3, "third order", "ruth 1983": Third-order Ruth (1983)
            - "symmetrized ruth 1983": Fourth-order symmetrized Ruth
            - 4, "fourth order", "suzuki 1990": Fourth-order Suzuki (1990)
            - 8, "eighth order", "morales 2022", "morales 2025": Eighth-order Morales et al.

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
            # Minimal first-order method
            return (1.0,)
        elif m in (2, "second order"):
            # Minimal second-order method
            return (0.5, 0.5)
        elif m in (3, "third order", "ruth 1983"):
            # Third-order Ruth (1983)
            return (7./24., 3./8., 3./8., -25./24., 1.0)
        elif m == "symmetrized ruth 1983":
            # Symmetrized version of Ruth (1983) that raises it to fourth order but uses twice as
            # many terms; expected to be less optimal than "suzuki 1990"
            return (7./48., 3./16., 3./16., -25./48., 0.5, 0.5, -25./48., 3./16., 3./16., 7./48.)
        elif m in (4, "fourth order", "suzuki 1990"):
            # Apply the recursion relation from Suzuki (1990) to the standard second-order method
            # to get a fourth-order method; also discussed in Hatano and Suzuki (2005) and
            # Ostmeyer (2023)
            s2 = 0.5 / (4.0 - cbrt(4.0))
            k = 0.5 - 4.0 * s2
            return (s2, s2, s2, s2, k, k, s2, s2, s2, s2)
        elif m in (8, "eighth order", "morales 2022", "morales 2025"):
            # Eighth-order method from Morales et al. (2022), recommended by Ostmeyer (2023);
            # paper modified and updated on arXiv as Morales et al. (2025)
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


def expand_ramped_trotterization(
    num_terms: int,
    coefficients: Sequence[float],
    num_steps: int,
    combine_terms: bool = True
) -> List[Tuple[int, float]]:
    """Expand ramped Trotterization into a flat list of (term_index, coefficient) pairs.

    In ramped Trotterization:
    - Each step consists of multiple ramps with different coefficients
    - Each ramp applies all terms in either ascending (0,1,2,...) or descending (...,2,1,0) order
    - Ramps alternate direction: ascending, descending, ascending, descending, etc.
    - The trailing term of one ramp matches the leading term of the next ramp
    - By default, adjacent occurrences of the same term are combined for efficiency

    Args:
        num_terms: Number of Pauli string terms in the Hamiltonian
        coefficients: Sequence of coefficients for each ramp in a step
        num_steps: Number of Trotterization steps
        combine_terms: If True (default), combine adjacent identical terms. If False, keep all terms separate.

    Returns:
        List of (term_index, coefficient) tuples representing the full sequence

    Example:
        >>> # 3 terms, 2 ramps per step, 1 step, with combining (default)
        >>> expand_ramped_trotterization(3, [0.5, 0.5], 1)
        [(0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)]

        >>> # Same but without combining
        >>> expand_ramped_trotterization(3, [0.5, 0.5], 1, combine_terms=False)
        [(0, 0.5), (1, 0.5), (2, 0.5), (2, 0.5), (1, 0.5), (0, 0.5)]

        Note: term 2 appears once with coefficient 1.0 when combined, twice with 0.5 when not
    """
    if num_terms == 0:
        raise ValueError("Must have at least one term")
    if num_steps == 0:
        raise ValueError("Must have at least one step")
    if len(coefficients) == 0:
        raise ValueError("Must have at least one coefficient")

    result = []
    ascending = True  # Start with ascending direction

    for step in range(num_steps):
        for coeff in coefficients:
            # Generate indices for this ramp
            if ascending:
                indices = range(num_terms)  # 0, 1, 2, ..., n-1
            else:
                indices = range(num_terms - 1, -1, -1)  # n-1, n-2, ..., 1, 0

            # Add each term from this ramp
            for idx in indices:
                result.append((idx, coeff))

            # Alternate direction for next ramp
            ascending = not ascending

    # Optionally combine adjacent terms with the same index
    if not combine_terms:
        return result

    # At this point, result is guaranteed to be non-empty due to validation
    combined = []
    current_idx, current_coeff = result[0]

    for idx, coeff in result[1:]:
        if idx == current_idx:
            # Same term as previous - combine coefficients
            current_coeff += coeff
        else:
            # Different term - save current and start new
            combined.append((current_idx, current_coeff))
            current_idx = idx
            current_coeff = coeff

    # Don't forget the last term
    combined.append((current_idx, current_coeff))

    return combined


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

    @classmethod
    def from_method(
        cls,
        pauli_terms: Sequence[Tuple[str, float]],
        method,
        time: float,
        num_steps: int,
        hbar: float = 1.0,
        combine_terms: bool = True
    ):
        """Create Trotterization using a named method.

        Args:
            pauli_terms: Sequence of (pauli_string, h_i) tuples
            method: Method name (e.g., "second order", "fourth order") or coefficients
            time: Total evolution time
            num_steps: Number of Trotterization steps
            hbar: Reduced Planck constant (default 1.0)
            combine_terms: If True (default), combine adjacent identical terms. If False, keep all terms separate.

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
            combine_terms=combine_terms
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
        return expand_ramped_trotterization(
            num_terms=len(self.pauli_terms),
            coefficients=self.coefficients,
            num_steps=self.num_steps,
            combine_terms=self.combine_terms
        )

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
        for term_idx, coeff in self.expanded_sequence:
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

    @property
    def num_terms(self) -> int:
        """Number of Pauli string terms in the Hamiltonian."""
        return len(self.pauli_terms)

    @property
    def num_operations(self) -> int:
        """Number of operations in the expanded sequence (after combining)."""
        return len(self.expanded_sequence)

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
    combine_terms: bool = True
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
        combine_terms=combine_terms
    )
