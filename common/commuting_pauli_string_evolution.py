"""
Commuting Pauli String Evolution Bloq

This module implements a Qualtran Bloq for exact time evolution under a sum of commuting
Pauli string Hamiltonians:
    U = exp(-i * (c₁P₁ + c₂P₂ + ...) * t / ℏ)
      = exp(-i * c₁P₁ * t / ℏ) * exp(-i * c₂P₂ * t / ℏ) * ...

This exact decomposition is valid when all Pauli strings commute with each other.
"""

import numbers
from functools import cached_property
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Tuple

import attrs
import numpy as np

from qualtran import (Bloq, BloqBuilder, CompositeBloq, CtrlSpec, QBit, Register, Side, Signature,
                      SoquetT)
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity

from qhat.common import controls
from qhat.common.hermiticity import require_real_scalar
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import validate_pauli_string

if TYPE_CHECKING:
    import quimb.tensor as qtn
    from qualtran import AddControlledT


def pauli_strings_commute(pauli_string1: str, pauli_string2: str) -> bool:
    """Check if two Pauli strings commute.

    Two Pauli strings commute if and only if they have an even number of positions
    where both have non-identity Pauli operators that are different.

    Args:
        pauli_string1: First Pauli string (e.g., "XYZI")
        pauli_string2: Second Pauli string (e.g., "IZXY")

    Returns:
        True if the Pauli strings commute, False otherwise.

    Example:
        >>> pauli_strings_commute("XY", "YX")  # 2 differences -> commute
        True
        >>> pauli_strings_commute("XY", "ZI")  # 1 difference -> don't commute
        False
        >>> pauli_strings_commute("XII", "IXI")  # 0 differences -> commute
        True
    """
    if len(pauli_string1) != len(pauli_string2):
        raise ValueError(
            f"Pauli strings must have the same length. "
            f"Got {len(pauli_string1)} and {len(pauli_string2)}"
        )

    # Count positions where both strings have non-identity Paulis that differ
    diff_count = 0
    for p1, p2 in zip(pauli_string1, pauli_string2):
        # Skip if either is identity
        if p1 == 'I' or p2 == 'I':
            continue
        # Count if they're different
        if p1 != p2:
            diff_count += 1

    # Commute if even number of differences
    return diff_count % 2 == 0


def all_commute(pauli_strings: Sequence[str]) -> bool:
    """Check if all Pauli strings in a sequence pairwise commute.

    Args:
        pauli_strings: Sequence of Pauli strings

    Returns:
        True if all pairs commute, False otherwise.
    """
    n = len(pauli_strings)
    for i in range(n):
        for j in range(i + 1, n):
            if not pauli_strings_commute(pauli_strings[i], pauli_strings[j]):
                return False
    return True


@attrs.frozen
class CommutingPauliStringEvolution(Bloq):
    """Exact time evolution for a sum of commuting Pauli string Hamiltonians.

    Implements the unitary operator:
        U = exp(-i * H * time / hbar)

    where H = sum_i (coefficient_i * PauliString_i) and all Pauli strings commute.

    Since all terms commute, we can exactly decompose:
        exp(-i * (A + B + ...) * t / ℏ) = exp(-i * A * t / ℏ) * exp(-i * B * t / ℏ) * ...

    Attributes:
        pauli_terms: Sequence of (pauli_string, coefficient) tuples.
            Example: [("XYZ", 0.5), ("ZZI", 0.3)]
        time: Evolution time parameter.
        hbar: Reduced Planck constant (default 1.0 for natural units).
        is_adjoint: If True, the operator is U^dagger. `adjoint()` flips it; `time` and the term
            order are left as given, and the decomposition applies the terms in reverse order.
        ctrl_values: Controls, as in `PauliStringEvolution`; each term gets the same controls.

    Example:
        >>> terms = [("XZ", 0.5), ("ZX", 0.3)]  # These commute
        >>> cpse = CommutingPauliStringEvolution(pauli_terms=terms, time=1.0)
        >>> U = cpse.tensor_contract()  # Get the unitary matrix
        >>> from qualtran.cirq_interop.t_complexity_protocol import t_complexity
        >>> tc = t_complexity(cpse)  # Get resource estimates
    """

    pauli_terms: Tuple[Tuple[str, float], ...]
    time: float
    hbar: float = 1.0
    is_adjoint: bool = False
    ctrl_values: Tuple[int, ...] = ()

    def __attrs_post_init__(self):
        """Validate that all Pauli strings commute and have the same length."""
        if len(self.pauli_terms) == 0:
            raise ValueError("Must provide at least one Pauli term.")

        # Extract Pauli strings and check they all have the same length
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
        object.__setattr__(self, "time",
                           require_real_scalar("CommutingPauliStringEvolution time", self.time))
        object.__setattr__(self, "hbar",
                           require_real_scalar("CommutingPauliStringEvolution hbar", self.hbar))
        object.__setattr__(self, "ctrl_values", controls.validate_ctrl_values(
            "CommutingPauliStringEvolution", self.ctrl_values))

        # Check that all Pauli strings pairwise commute
        if not all_commute(pauli_strings):
            raise ValueError(
                "Not all Pauli strings commute. CommutingPauliStringEvolution requires "
                "all terms to pairwise commute."
            )

    @classmethod
    def from_pauli_string_iterable(
        cls,
        pauli_strings: Iterable[Tuple[str, float]],
        time: float,
        hbar: float = 1.0
    ):
        """Create from an iterable of (pauli_string, coefficient) tuples.

        This is the format used by hamiltonian.get_all_pauli_strings(return_as='strings').

        Args:
            pauli_strings: Iterable of (string, coefficient) tuples
            time: Evolution time
            hbar: Reduced Planck constant

        Returns:
            CommutingPauliStringEvolution instance
        """
        return cls(pauli_terms=tuple(pauli_strings), time=time, hbar=hbar)

    @cached_property
    def num_qubits(self) -> int:
        """Number of qubits the operator acts on."""
        # Guaranteed to have at least one term due to __attrs_post_init__ validation
        return len(self.pauli_terms[0][0])

    @cached_property
    def signature(self) -> Signature:
        """The `ctrl` register (if controlled), then a register of qubits."""
        return Signature(controls.ctrl_registers(self.ctrl_values)
                         + [Register('q', QBit(), shape=(self.num_qubits,), side=Side.THRU)])

    def build_composite_bloq(self, bb: BloqBuilder, **soqs: SoquetT) -> Dict[str, SoquetT]:
        """Decompose into a sequence of PauliStringEvolution bloqs, each with the same controls.

        Since all terms commute, the exact evolution is the product of individual evolutions.
        """
        # Apply each PauliStringEvolution in sequence (reversed for the adjoint)
        for pse in self._term_evolutions():
            soqs = bb.add_d(pse, **soqs)

        return soqs

    def _term_evolutions(self) -> Tuple[PauliStringEvolution, ...]:
        """One PauliStringEvolution per term, in the order they are applied."""
        terms = reversed(self.pauli_terms) if self.is_adjoint else self.pauli_terms
        return tuple(
            PauliStringEvolution(
                pauli_string=pauli_string,
                coefficient=coefficient,
                time=self.time,
                hbar=self.hbar,
                is_adjoint=self.is_adjoint,
                ctrl_values=self.ctrl_values,
            )
            for pauli_string, coefficient in terms
        )

    def _t_complexity_(self) -> TComplexity:
        """Return the T-complexity by summing individual term complexities."""
        total_complexity = TComplexity()
        for pse in self._term_evolutions():
            total_complexity += t_complexity(pse)
        return total_complexity

    @property
    def num_terms(self) -> int:
        """Number of Pauli string terms in the Hamiltonian."""
        return len(self.pauli_terms)

    def get_ctrl_system(
        self, ctrl_spec: Optional[CtrlSpec] = None
    ) -> Tuple['CommutingPauliStringEvolution', 'AddControlledT']:
        """Controlled version: the same bloq with the new control value prepended to `ctrl_values`.

        Only single-qubit controls on |1> or |0> are supported; others raise ValueError.
        """
        return controls.get_ctrl_system(self, ctrl_spec)

    def adjoint(self) -> 'CommutingPauliStringEvolution':
        """Return U^dagger (flips `is_adjoint`; `time` and the terms are unchanged)."""
        return attrs.evolve(self, is_adjoint=not self.is_adjoint)

    def __pow__(self, power: int) -> 'CommutingPauliStringEvolution':
        """U**k is the same evolution for k times the time; U**(-k) is (U**k).adjoint().

        Raises ValueError for 0 and for non-integer powers.
        """
        if isinstance(power, bool) or not isinstance(power, numbers.Integral) or power == 0:
            raise ValueError(f"Power must be a non-zero integer, got {power}")
        return attrs.evolve(
            self, time=self.time * abs(int(power)), is_adjoint=self.is_adjoint != (power < 0))

    def __str__(self) -> str:
        term_strs = [f"{coeff}*{ps}" for ps, coeff in self.pauli_terms]
        H_str = " + ".join(term_strs)
        text = f"exp(-i * ({H_str}) * {self.time} / {self.hbar})"
        return controls.ctrl_str(self.ctrl_values, f"[{text}]†" if self.is_adjoint else text)

    def __repr__(self) -> str:
        adjoint = ", is_adjoint=True" if self.is_adjoint else ""
        ctrl = f", ctrl_values={self.ctrl_values}" if self.ctrl_values else ""
        return (f"CommutingPauliStringEvolution(pauli_terms={self.pauli_terms}, "
                f"time={self.time}, hbar={self.hbar}{adjoint}{ctrl})")
