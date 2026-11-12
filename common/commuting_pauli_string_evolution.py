"""
Commuting Pauli String Evolution Bloq

This module implements a Qualtran Bloq for exact time evolution under a sum of commuting
Pauli string Hamiltonians:
    U = exp(-i * (c₁P₁ + c₂P₂ + ...) * t / ℏ)
      = exp(-i * c₁P₁ * t / ℏ) * exp(-i * c₂P₂ * t / ℏ) * ...

This exact decomposition is valid when all Pauli strings commute with each other.
"""

import numbers
from collections import Counter
from functools import cached_property
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Set, Tuple

import attrs
import numpy as np

from qualtran import (Bloq, BloqBuilder, CompositeBloq, CtrlSpec, QBit, Register, Side, Signature,
                      SoquetT)
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity
from qualtran.resource_counting import QubitCount

from qhat.common import controls
from qhat.common.hermiticity import require_real_scalar
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import validate_pauli_string

if TYPE_CHECKING:
    import quimb.tensor as qtn
    from qualtran import AddControlledT
    from qualtran.resource_counting import BloqCountT, CostKey, SympySymbolAllocator


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
        multi_control_method: As in `PauliStringEvolution`. With "ancilla" and 2 or more controls,
            the controls are ANDed once for the whole bloq, and each term gets one control.
        controlled_rotation_method: As in `PauliStringEvolution`. With "ancilla", the terms leave
            out their phases on the control, and this bloq applies their sum once.
        omit_control_phase: As in `PauliStringEvolution`: if True, the summed phase is left out
            too, for a parent (`Trotterization`) to apply.

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
    multi_control_method: str = controls.DEFAULT_MULTI_CONTROL_METHOD
    controlled_rotation_method: str = controls.DEFAULT_CONTROLLED_ROTATION_METHOD
    omit_control_phase: bool = False

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
        controls.validate_multi_control_method(
            "CommutingPauliStringEvolution", self.multi_control_method)
        controls.validate_controlled_rotation_method(
            "CommutingPauliStringEvolution", self.controlled_rotation_method, self.ctrl_values,
            self.multi_control_method)
        controls.validate_omit_control_phase(
            "CommutingPauliStringEvolution", self.omit_control_phase, self.ctrl_values,
            self.controlled_rotation_method)

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
        With 2 or more controls and the "ancilla" method, the controls are first ANDed into one
        ancilla, which controls this bloq. With the "ancilla" rotation method, the terms' phases
        on the control are applied last, as one phase (unless `omit_control_phase`).
        """
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.add_with_and_ancilla(bb, self, **soqs)
        # Apply each PauliStringEvolution in sequence (reversed for the adjoint)
        for pse in self._term_evolutions():
            soqs = bb.add_d(pse, **soqs)
        if self._applies_merged_phase:
            ctrl, _ = controls.add_controlled_phase(
                bb, self.control_phase, soqs['ctrl'], self.ctrl_values)
            soqs['ctrl'] = np.array(ctrl)

        return soqs

    @cached_property
    def _merges_control_phases(self) -> bool:
        """True if the terms leave out their phases on the control (the "ancilla" rotation
        method with one control)."""
        return controls.uses_rotation_ancilla(self.ctrl_values, self.controlled_rotation_method)

    @cached_property
    def _applies_merged_phase(self) -> bool:
        """True if this bloq applies the terms' summed phase on the control itself."""
        return self._merges_control_phases and not self.omit_control_phase

    @cached_property
    def control_phase(self) -> float:
        """The sum of the terms' `control_phase`s: the phase exp(-i phi) the "ancilla" rotation
        method leaves on the active control. A parent adds these up and applies them once."""
        return sum(pse.control_phase for pse in self._term_evolutions())

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
                multi_control_method=self.multi_control_method,
                controlled_rotation_method=self.controlled_rotation_method,
                omit_control_phase=self._merges_control_phases,
            )
            for pauli_string, coefficient in terms
        )

    def _merged_phase_call_graph(self) -> Counter:
        """The gates of the summed phase on the control, if this bloq applies it."""
        if not self._applies_merged_phase:
            return Counter()
        return controls.controlled_phase_call_graph(self.control_phase, self.ctrl_values)

    def build_call_graph(self, ssa: 'SympySymbolAllocator') -> Set['BloqCountT']:
        """The decomposition's term evolutions (or, for the "ancilla" method with 2 or more
        controls, the AND, its uncompute and this bloq on one control), without building it."""
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.and_ancilla_call_graph(self)
        return set((Counter(self._term_evolutions()) + self._merged_phase_call_graph()).items())

    def _t_complexity_(self) -> TComplexity:
        """Return the T-complexity by summing individual term complexities (plus the AND ladder
        for the "ancilla" method, and the summed phase for the "ancilla" rotation method)."""
        if controls.uses_and_ancilla(self.ctrl_values, self.multi_control_method):
            return controls.and_ancilla_t_complexity(self)
        total_complexity = TComplexity()
        for pse in self._term_evolutions():
            total_complexity += t_complexity(pse)
        if self._applies_merged_phase:
            total_complexity += controls.controlled_phase_t_complexity(self.ctrl_values)
        return total_complexity

    def my_static_costs(self, cost_key: 'CostKey'):
        """`QubitCount` without decomposing: qubits + controls + AND ancillas (if any), + 1 for the
        "ancilla" rotation method (unless every term is all identities)."""
        if isinstance(cost_key, QubitCount):
            has_rotations = any(set(ps) != {'I'} for ps, _ in self.pauli_terms)
            return controls.peak_qubits(
                self.num_qubits, self.ctrl_values, self.multi_control_method,
                controls.num_rotation_ancillas(self.ctrl_values, self.controlled_rotation_method,
                                               has_rotations))
        return NotImplemented

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
        ctrl = controls.ctrl_repr(self.ctrl_values, self.multi_control_method,
                                  self.controlled_rotation_method)
        omit = ", omit_control_phase=True" if self.omit_control_phase else ""
        return (f"CommutingPauliStringEvolution(pauli_terms={self.pauli_terms}, "
                f"time={self.time}, hbar={self.hbar}{adjoint}{ctrl}{omit})")
