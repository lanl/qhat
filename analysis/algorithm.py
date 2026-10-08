# pyLIQTR provides the following phase estimation tools
# -- PhaseEstimation.pe.PhaseEstimation: Only works for Trotterization, has not been a priority for
#    the pyLIQTR team to develop this to run efficiently
# -- qubitization.phase_estimation.QubitizedPhaseEstimation: Our double-factorization script uses
#    this.  Currently it's unclear to me if this is restricted to "qubitized" methods, and what
#    pyLIQTR consideres a "qubitized" method or not (because I need to understand "qubitization"
#    better).
#    -- QubitizedWalkOperator
#    -- n = 1...precision
#       -- QubitizedReflection
#       -- QubitizedWalkOperator (2^n times)
#       -- QubitizedReflection
# I once traced PhaseEstimation and it does not follow the "usual" phase estimation algorithm
# design.  The QubitizedPhaseEstimation involves quantum walk operators (add that to my list of
# things to learn about), and it's not clear if it follows the usual phase estimation algorithm or
# not.
# 
# This Qualtran issue provides nice links to the "standard" version and the "walk" version
# (although at a quick glance it looks as though both use multiple ancilla qubits, while both
# pyLIQTR versions use a single ancilla qubit): https://github.com/quantumlib/Qualtran/issues/819.
# This brings us to a total of four versions of QPE, and it's not clear to me how much overlap
# there is between the methods.
# 
# I can probably implement the "standard" QPE algorithm in a useful framework, but I will have to
# read more to understand all the variations and figure out how to implement them.
# 
# Addendum: Qualtran provides QPE algorithms, including TextbookQPE.  It may not be the most
# efficient (?) but it may be a reliable starting point.

import logging
import math

from qualtran import CtrlSpec
from qualtran.bloqs.phase_estimation import TextbookQPE, QubitizationQPE
from qualtran.bloqs.qubitization.qubitization_walk_operator import QubitizationWalkOperator

from pyLIQTR.qubitization.phase_estimation import QubitizedPhaseEstimation

from qhat.analysis.config_types import AlgorithmConfiguration
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.qpe_window_state import precision_bits, RectangularWindowState

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

# We define our own version of QubitizationQPE, in order to fix a bug
# in the qualtran 0.4.0 version
class NewQubitizationQPE(QubitizationQPE):
    def build_call_graph(self, ssa: 'SympySymbolAllocator'):
        # Assumes self.unitary is not fast forwardable.
        M = 2**self.m_bits
        return {
            (self.ctrl_state_prep, 1),
            #(self.walk.controlled(control_values=[1]), 1),
            #(self.walk.reflect.controlled(control_values=[0]), 2 * (self.m_bits - 1)),
            # === bugfix, backported from 0.5.0 ===
            (self.walk.controlled(), 1),
            (self.walk.reflect.controlled(ctrl_spec=CtrlSpec(cvs=0)), 2 * (self.m_bits - 1)),
            # === end bugfix ===
            (self.walk, M - 2),
            (self.qft_inv, 1),
        }

# -------------------------------------------------------------------------------------------------

def qpe_requirements(config_algorithm: AlgorithmConfiguration, phase_error):
    """Check and return `(phase_error, probability_of_failure)` for sizing a non-qubitized QPE.

    Used when `algorithm.num_phase_qubits` is not set. `phase_error` is in turns (see
    `qpe_phase_error`).
    """
    if phase_error is None:
        raise ValueError(
            "QPE needs either algorithm.num_phase_qubits or algorithm.energy_error.")
    probability_of_failure = config_algorithm.probability_of_failure
    if probability_of_failure is None:
        raise ValueError(
            "QPE sized from algorithm.energy_error also needs algorithm.probability_of_failure.")
    if not 0 < probability_of_failure < 1:
        raise ValueError(
            f"algorithm.probability_of_failure must be in (0, 1), got {probability_of_failure}.")
    logger.verbose(f"-- target phase error = {phase_error} turns, "
                   f"probability of failure = {probability_of_failure}")
    return phase_error, probability_of_failure

# -------------------------------------------------------------------------------------------------

def build_qpe_qualtran_textbook(
        config_algorithm: AlgorithmConfiguration,
        unitary,
        phase_error):

    logger.verbose("Build a QPE algorithm with Qualtran's \"textbook\" method.")

    # TODO: There is a note in the documentation (see link below) that a fast-forwardable unitary
    #       can lower the cost from (2^m - 1) * cost(C-U) to m * cost(C-U).  If we have a
    #       continuous time-evolution Hamiltonian, we might be able to implement this.  I think it
    #       essentially means, for example, something like if U = e^{i H dt} then instead of using
    #       U^{2^n) you use e^{i H 2^n dt}, and that may be implementable with the same cost (in
    #       terms of gates) as U itself.  If we can demonstrate that a unitary is
    #       "fast-forwardable", then we could get _significant_ improvements in T counts (reduce
    #       from O(2^m) to O(m^2) or even O(m log m) with a faster approximate iQFT).  We'll need
    #       to look into whether or not TextbookQPE already accounts for this, whether or not we're
    #       getting this improvement, and how we could ensure that we do get this improvement.
    # https://qualtran.readthedocs.io/en/latest/bloqs/phase_estimation/text_book_qpe.html#cost-of-textbookqpe
    #       -- U^{2^i} is implemented through cirq.pow:
    #          https://github.com/quantumlib/Qualtran/blob/main/qualtran/bloqs/phase_estimation/text_book_qpe.py#L180C23-L180C25: 
    #       -- cirq.pow will look for U.__pow__ and use that if available; otherwise it will use
    #          some default (I didn't yet read that far):
    #          https://github.com/quantumlib/Cirq/blob/v1.5.0/cirq-core/cirq/protocols/pow_protocol.py#L79
    #       -- So long as we implement `unitary` so that it has a `__pow__` method, then we get the
    #          fast-forwarding improvement.
    #       -- Does pyLIQTR's implementation(s) of QPE get the same enhancement?
    #          -- It looks like pyLIQTR eventually gets down to OpenFermion, which builds on a base
    #             class that does implement __pow__: https://github.com/quantumlib/OpenFermion/blob/master/src/openfermion/ops/operators/symbolic_operator.py#L577
    #          -- pyLIQTR _may_ be getting the fast-forward behavior when using Trotterization and
    #             generating the algorithm:
    #             https://github.com/quantumlib/OpenFermion/blob/master/src/openfermion/ops/operators/symbolic_operator.py#L577
    #          -- When using an arbitrary unitary, pyLIQTR appears to just add the unitary 2^i
    #             times, so it probably is not getting the fast-forward behavior.
    #          -- It doesn't appear that PhaseEstimation is using the _t_complexity_ protocol, so
    #             it's probably just counting the gates manually
    #          -- But QubitizedPhaseEstimation is completely different and does use the
    #             _t_complexity_ protocol, which explicitly includes the (2*i - 1) factor (actually
    #             in an indirect way as a summation, so that it's not as efficient as it could be
    #             even there), so you don't get T-gate estimates with the fast-forward behavior.
    #       -- The best way to confirm the behavior of a given QPE method would be to test it with
    #          the same unitary and all QPE settings the same except for varying the number of
    #          phase qubits, then checking how that scales.
    #       -- Trotterization probably can't get the fast-forward behavior: U^k would basically be
    #          implemented by multiplying the number of steps by k

    # TODO: This uses the default QFT (QFTTextBook(self.m_bits).adjoint()), but we could make it a
    #       user-configurable setting to switch to other iQFT implementations.  See also
    #       https://qualtran.readthedocs.io/en/latest/bloqs/phase_estimation/text_book_qpe.html#cost-of-textbookqpe
    # TODO: This uses the "textbook" state initialization for the phase qubits.  There are other
    #       options, including KaiserWindowState and LPResourceState.  The KaiserWindowState was
    #       added after the version of Qualtran that I'm currently using (0.4.0), but later
    #       versions of Qualtran add a Jupyter notebook to go with the KaiserWindowState that
    #       compares the QPE performance with different window state objects.
    #       -- The RectangularWindowState isn't added until a later version of qualtran than the
    #          one I'm using.  The interface changes in later versions.
    P = config_algorithm.num_phase_qubits
    if P is not None:
        logger.verbose(f"-- using user-specified number of phase qubits ({P})")
        qpe = TextbookQPE(unitary, P)
    else:
        phase_error, probability_of_failure = qpe_requirements(config_algorithm, phase_error)
        qpe = TextbookQPE.from_precision_and_delta(
                unitary, precision_bits(phase_error), probability_of_failure)
    logger.verbose(f"-- number of phase qubits = {qpe.m_bits}")
    return qpe

# -------------------------------------------------------------------------------------------------

# Config names for FlexibleQPE components.  A qft_inv entry maps the phase-register size to a bloq;
# None selects FlexibleQPE's default (the textbook inverse QFT).
ANCILLA_PREPS = {"rectangular": RectangularWindowState}
INVERSE_QFTS = {"textbook": None}

def build_qpe_qhat_flexible(
        config_algorithm: AlgorithmConfiguration,
        unitary,
        phase_error):

    logger.verbose("Build a QPE algorithm with QHAT's FlexibleQPE.")

    ancilla_prep_name = (config_algorithm.ancilla_prep or "rectangular").lower()
    if ancilla_prep_name not in ANCILLA_PREPS:
        raise ValueError(f"Invalid QPE ancilla_prep \"{config_algorithm.ancilla_prep}\".")
    qft_inv_name = (config_algorithm.qft_inv or "textbook").lower()
    if qft_inv_name not in INVERSE_QFTS:
        raise ValueError(f"Invalid QPE qft_inv \"{config_algorithm.qft_inv}\".")

    components = dict(ancilla_prep=ANCILLA_PREPS[ancilla_prep_name],
                      qft_inv=INVERSE_QFTS[qft_inv_name])
    P = config_algorithm.num_phase_qubits
    if P is not None:
        logger.verbose(f"-- using user-specified number of phase qubits ({P})")
        qpe = FlexibleQPE.from_num_phase_qubits(unitary, P, **components)
    else:
        phase_error, probability_of_failure = qpe_requirements(config_algorithm, phase_error)
        qpe = FlexibleQPE.from_requirements(
                unitary, phase_error, probability_of_failure, **components)
    logger.verbose(f"-- number of phase qubits = {qpe.num_ancilla_qubits}")
    return qpe

# -------------------------------------------------------------------------------------------------

def build_qpe_qualtran_qubitized(
        config_algorithm: AlgorithmConfiguration,
        unitary):

    logger.verbose(
            "Build a QPE algorithm with Qualtran's \"QubitizationQPE\" method.")

    P = config_algorithm.num_phase_qubits
    if P is None:
        dE = config_algorithm.energy_error
        alpha = unitary.alpha
        P = int(math.ceil(math.log2(math.pi * alpha / (2 * dE))))

    return NewQubitizationQPE(QubitizationWalkOperator(unitary._select_gate,
                                                       unitary._prepare_gate),
                              P)

# -------------------------------------------------------------------------------------------------

def build_qpe_pyliqtr_qubitized(
        config_algorithm: AlgorithmConfiguration,
        unitary):

    logger.verbose(
            "Build a QPE algorithm with pyLIQTR's \"QubitizedPhaseEstimation\" method.")

    P = config_algorithm.num_phase_qubits
    if P is None:
        dE = config_algorithm.energy_error
        alpha = unitary.alpha
        P = int(math.ceil(math.log2(math.pi * alpha / (2 * dE))))

    # TODO: The name and signature suggest that this may _only_ be valid for block-encoded
    #       unitaries.  Is that true?
    return QubitizedPhaseEstimation(block_encoding=unitary, prec=P)

# -------------------------------------------------------------------------------------------------

def build_time_evolution(
        config_algorithm: AlgorithmConfiguration,
        unitary):

    logger.verbose("Build a time evolution algorithm.")

    return unitary

# -------------------------------------------------------------------------------------------------

def build_controlled_time_evolution(
        config_algorithm: AlgorithmConfiguration,
        unitary):

    logger.verbose("Build a singly-controlled time evolution algorithm.")

    return unitary.controlled()

# -------------------------------------------------------------------------------------------------

def qpe_phase_error(config_algorithm: AlgorithmConfiguration, unitary):
    """Phase error (in turns) corresponding to algorithm.energy_error for this unitary."""
    if config_algorithm.energy_error is None:
        return None
    if not hasattr(unitary, "phase_error_from_energy_error"):
        raise ValueError(
            f"Cannot size QPE from algorithm.energy_error: the unitary ({type(unitary).__name__}) "
            "does not provide phase_error_from_energy_error().  Set algorithm.num_phase_qubits.")
    return unitary.phase_error_from_energy_error(config_algorithm.energy_error)

# -------------------------------------------------------------------------------------------------

def build_algorithm(
        config_algorithm: AlgorithmConfiguration,
        unitary):

    logger.info("Beginning to construct quantum algorithm.")

    if config_algorithm.method.lower() in ("qpe: qualtran textbook",):
        return build_qpe_qualtran_textbook(
                config_algorithm, unitary, qpe_phase_error(config_algorithm, unitary))
    elif config_algorithm.method.lower() in ("qpe: qhat flexible",):
        return build_qpe_qhat_flexible(
                config_algorithm, unitary, qpe_phase_error(config_algorithm, unitary))
    elif config_algorithm.method.lower() in ("qpe: qualtran qubitization",):
        # TODO: This may be more specialized (for LCU only?), but I'm not yet sure of the details.
        return build_qpe_qualtran_qubitized(config_algorithm, unitary)
    elif config_algorithm.method.lower() in ("qpe: pyliqtr qubitized",):
        return build_qpe_pyliqtr_qubitized(config_algorithm, unitary)
    elif config_algorithm.method.lower() in ("time evolution",):
        return build_time_evolution(config_algorithm, unitary)
    elif config_algorithm.method.lower() in ("controlled time evolution",):
        return build_controlled_time_evolution(config_algorithm, unitary)
    else:
        raise ValueError(f"Invalid algorithm method \"{config_algorithm.method}\".")

# -------------------------------------------------------------------------------------------------

def qpe_energy_window(config_algorithm: AlgorithmConfiguration, Elo, Ehi):
    """Energy window [lo, hi) to map onto one full phase turn (t = 2 pi / (hi - lo)).

    With algorithm.energy_error set, the width is rounded up to 2**P0 * energy_error, so each bin
    of a P0-bit phase register is exactly energy_error wide; the extra width is split evenly above
    and below [Elo, Ehi).  Otherwise the window is [Elo, Ehi).  For example, [-1.5, 1.5) with
    energy_error 0.1 gives P0 = 5 and [-1.6, 1.6).
    """
    logger.info("Computing the QPE energy window.")

    dE = config_algorithm.energy_error
    if dE is None:
        return (Elo, Ehi)

    P0 = math.ceil(math.log2((Ehi - Elo) / dE))
    width = 2**P0 * dE
    middle = (Elo + Ehi) / 2
    lo, hi = middle - width / 2, middle + width / 2
    logger.verbose(f"-- P0 = {P0} bits of {dE} each, window = [{lo}, {hi})")
    return (lo, hi)
