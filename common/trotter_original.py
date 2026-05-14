# TODO: I need to write some tests for this.  The error analysis is suggesting that there may be
#       bugs in this code.  There may well also be bugs in the error analysis code, but I need to
#       verify this works correctly.
#    -- Something like aX + bY + cZ would be very simple, so I could probably compute it
#       analytically for definite comparison.  It has three operators, none of which commute with
#       each other.

# TODO: Check how things are computed for each individual term.  I think the cost of `a^† a` is
#       less than the cost of `a^†` or `a` by themselves.
#       -- This seems too easy and obvious of a gain for Qualtran and/or Cirq to not already have
#          it, but it's worth checking on.

# =================================================================================================

import attrs
from cirq import LineQubit
from common.dense_pauli_exp import DensePauliString
from functools import cached_property, reduce
from math import cbrt
from numpy import exp
from pyLIQTR.utils.pauli_string_manip import convert_to_dense_pauli_string
from qualtran import Bloq, BloqBuilder, Signature, SoquetT
from qualtran.cirq_interop import CirqGateAsBloq
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity
from qualtran.symbolics import SymbolicFloat
from typing import Dict, Sequence

# =================================================================================================

@attrs.frozen
class TrotterRamp(Bloq):
    bloqs: Sequence[Bloq]
    timestep: SymbolicFloat
    direction: bool
    ASCENDING: bool = True
    DESCENDING: bool = False

    # TODO: Possible improvement: try to rearrange/group bloqs together to limit the number of
    #       terms in the Trotterization.  Not sure if this is possible/easy to do here.  There may
    #       also be tools in Cirq and/or Qualtran to do this for us.

    @cached_property
    def signature(self) -> Signature:
        return self.bloqs[0].signature

    def _bloq_iterator(self):
        if self.direction: # True = ascending
            return self.bloqs
        else: # False = descending
            return reversed(self.bloqs)

    def switch(self, direction):
        return not direction

    def build_composite_bloq(self, bb: 'BloqBuilder', **soqs: SoquetT) -> Dict[str, 'SoquetT']:
        # TODO: Check that we got the right coefficients in the right place.  Qualtran and Cirq
        #       sometimes use surprising conventions (e.g., factors of 2, etc), and I also need to
        #       trace my logic to ensure that I didn't drop coefficients somewhere between the
        #       definition of the problem and the final result.
        for bloq in self._bloq_iterator():
            soqs |= bb.add_d(bloq**self.timestep, **soqs)
        return soqs

    def _t_complexity_(self) -> TComplexity:
        return reduce(
                lambda a, b: a + t_complexity(b**self.timestep),
                self.bloqs,
                TComplexity())

# =================================================================================================

@attrs.frozen
class RampedTrotterStep(Bloq):
    bloqs: Sequence[Bloq]
    coefficients: Sequence[SymbolicFloat]
    timestep: SymbolicFloat

    @cached_property
    def signature(self) -> Signature:
        return self.bloqs[0].signature

    def _build_ramp(self, s_dt, direction):
        return TrotterRamp(self.bloqs, s_dt, direction)

    def build_composite_bloq(self, bb: 'BloqBuilder', **soqs: SoquetT) -> Dict[str, 'SoquetT']:
        direction = TrotterRamp.DESCENDING
        for coefficient in self.coefficients:
            direction = not direction # TODO: This is hacky (and may not even work)
            soqs |= bb.add_d(self._build_ramp(coefficient * self.timestep, direction), **soqs)
        return soqs

    def _t_complexity_(self) -> TComplexity:
        # The cost of a ramp is independent of the actual value of the coefficient, so we just
        # simply use coefficient = 1 when estimating resources.
        cost1 = t_complexity(self._build_ramp(self.timestep, TrotterRamp.ASCENDING))
        count = len(self.coefficients)
        return count * cost1
        # TODO: Future improvements could include accounting for the fact that the last term in one
        #       ramp is the same as the first step in the next ramp (because ramps alternate
        #       direction), so they can be combined with a new coefficient.  Since the resource
        #       cost is independent of the coefficient, combining these two terms means you only
        #       have to evaluate the term once instead of twice.

# =================================================================================================

@attrs.frozen
class RampedTrotterizedUnitary(Bloq):
    bloqs: Sequence[Bloq]
    coefficients: Sequence[SymbolicFloat]
    timestep: SymbolicFloat
    numsteps: int

    def __attrs_post_init__(self):
        ref_sig = self.bloqs[0].signature
        for bloq in self.bloqs:
            if bloq.signature != ref_sig:
                raise ValueError(
                    f"Bloqs must have the same signature. Got {ref_sig} and {bloq.signature}"
                )
            if not attrs.has(bloq.__class__):
                raise ValueError("Bloq must be an attrs dataclass.")

    @cached_property
    def signature(self) -> Signature:
        return self.bloqs[0].signature

    def _generate_step_(self):
        return RampedTrotterStep(self.bloqs, self.coefficients, self.timestep / self.numsteps)

    def build_composite_bloq(self, bb: 'BloqBuilder', **soqs: SoquetT) -> Dict[str, 'SoquetT']:
        # quick preliminary implementation based on https://github.com/quantumlib/Qualtran/blob/main/qualtran/bloqs/chemistry/trotter/trotterized_unitary.py#L105
        for n in range(self.numsteps):
            soqs |= bb.add_d(self._generate_step_(), **soqs)
        return soqs

    def __pow__(self, exponent: int):
        return attrs.evolve(self, numsteps = exponent * self.numsteps)

    def _t_complexity_(self) -> TComplexity:
        return self.numsteps * t_complexity(self._generate_step_())
        # TODO: Future improvements could include accounting for the fact that the last term in one
        #       step is the same as the first step in the next step (assuming you have an even
        #       number of coefficients), leading to an optimization like that between ramps (see
        #       RampedTrotterStep).

# =================================================================================================

# TODO: This be a @classmethod of RampedTrotterizedUnitary
def build_bloqs(pauli_strings):
    return tuple(exp(1j * DensePauliString.from_dense_pauli_string(*term))
        for term in pauli_strings)

# -------------------------------------------------------------------------------------------------

def build_coefficients(method):
    if isinstance(method, str) or isinstance(method, int):
        m = method
        if isinstance(method, str):
            m = method.lower()
        if m in (1, "first order"):
            # minimal first-order method
            return (1.0,)
        elif m in (2, "second order"):
            # minimal second-order method
            return (0.5, 0.5)
        elif m in (3, "third order", "ruth 1983"):
            # third-order
            return (7./24., 3./8., 3./8., -25./24., 1.0)
        elif m == "symmetrized ruth 1983":
            # Symmetrized version of Ruth (1983) that raises it to fourth order but uses twice as
            # many terms; expected to be less optimal than "suzuki 2005", hence why that method is
            # the default "fourth order" method
            return (7./48., 3./16., 3./16., -25./48., 0.5, 0.5, -25./48., 3./16., 3./16., 7./48.)
        elif m in (4, "fourth order", "suzuki 1990"):
            # Apply the recursion relation from Suzuki (1990) to the standard second-order method
            # to get a fourth-order method; also discussed in Hatano and Suzuki (2005) and Ostmeyer
            # (2023)
            s2 = 0.5 / (4.0 - cbrt(4.0))
            k = 0.5 - 4.0 * s2
            return (s2, s2, s2, s2, k, k, s2, s2, s2, s2)
        elif m in (8, "eighth order", "morales 2022", "morales 2025"):
            # eighth-order method from Morales et al. (2022), recommended by Ostmeyer (2023); paper
            # modified and updated on arXiv as Morales et al. (2025)
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
        # If it's not a string naming a method, we assume it's a list of coefficients
        return method

# -------------------------------------------------------------------------------------------------

def build_ramped_trotterized_unitary(pauli_strings, method, timestep, numsteps):

    # Each bloq implements exp(i h_n / hbar), where h_n is a single term in the Hamiltonian
    # represented as a "Pauli string" (a set of Pauli operators applied to different qubits)
    bloqs = build_bloqs(pauli_strings)

    coefficients = build_coefficients(method)

    return RampedTrotterizedUnitary(bloqs, coefficients, timestep, numsteps)
