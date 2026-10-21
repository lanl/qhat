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
from qhat.common.hermiticity import require_real_scalar
from qhat.common.qpe_window_state import precision_bits, RectangularWindowState
from qhat.common.tapered_qpe_sizing import (
    COST_CEILING,
    fewest_runs,
    multi_run_failure,
    size_tapered_qpe,
)
from qhat.common.tapered_window_state import (
    BandwidthOptimizedWindowState,
    DPSSWindowState,
    TaperedWindowStateBase,
    TruncatedDPSSWindowState,
)

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

def qpe_requirements(config_algorithm: AlgorithmConfiguration, unitary):
    """Check and return `(phase_error, probability_of_failure)` for sizing a non-qubitized QPE.

    Used only when `algorithm.num_phase_qubits` is not set. `phase_error` is in turns (see
    `qpe_phase_error`).
    """
    phase_error = qpe_phase_error(config_algorithm, unitary)
    if phase_error is None:
        raise ValueError(
            "QPE needs either algorithm.num_phase_qubits or algorithm.energy_error.")
    phase_error = require_real_scalar("QPE phase error", phase_error)
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
        unitary):

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
        phase_error, probability_of_failure = qpe_requirements(config_algorithm, unitary)
        qpe = TextbookQPE.from_precision_and_delta(
                unitary, precision_bits(phase_error), probability_of_failure)
    logger.verbose(f"-- number of phase qubits = {qpe.m_bits}")
    return qpe

# -------------------------------------------------------------------------------------------------

# Config names for FlexibleQPE components.  A qft_inv entry maps the phase-register size to a bloq;
# None selects FlexibleQPE's default (the textbook inverse QFT).
ANCILLA_PREPS = {
    "rectangular": RectangularWindowState,
    "dpss": DPSSWindowState,
    "truncated dpss": TruncatedDPSSWindowState,
    "bandwidth-optimized": BandwidthOptimizedWindowState,
}
INVERSE_QFTS = {"textbook": None}

# Multi-run sizing of tapered windows (see qhat.common.tapered_qpe_sizing).
DEFAULT_OVERLAP = 0.5
OVERLAP_TABLE = (0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
MULTI_RUN_OPTIONS = ("overlap", "force_lowest_estimate")
SINGLE_RUN_NOTE = "-- NOTE: the analyses below report on a single run of the algorithm."


def warn_unused_multi_run_options(config_algorithm: AlgorithmConfiguration, where: str):
    for name in MULTI_RUN_OPTIONS:
        if getattr(config_algorithm, name) is not None:
            logger.warning(
                f"algorithm.{name} is ignored with {where}: it applies only to \"QPE: QHAT "
                "flexible\" with a tapered ancilla_prep (\"dpss\", \"truncated dpss\" or "
                "\"bandwidth-optimized\") and algorithm.probability_of_failure.")


def multi_run_options(config_algorithm: AlgorithmConfiguration):
    """Check algorithm.overlap and algorithm.force_lowest_estimate and fill in their defaults.

    Returns `(overlap, overlap_is_default, force_lowest_estimate, force_is_default)`, e.g.
    `(0.5, True, False, True)` when neither is set.
    """
    overlap = config_algorithm.overlap
    overlap_is_default = overlap is None
    if overlap_is_default:
        overlap = DEFAULT_OVERLAP
    else:
        overlap = require_real_scalar("algorithm.overlap", overlap)
        if not 0 < overlap <= 1:
            raise ValueError(f"algorithm.overlap must be in (0, 1], got {overlap}.")
    force = config_algorithm.force_lowest_estimate
    force_is_default = force is None
    if force_is_default:
        force = False
    elif not isinstance(force, bool):
        raise TypeError(
            f"algorithm.force_lowest_estimate must be True or False, got {force!r}.")
    return overlap, overlap_is_default, force, force_is_default


def log_multi_run_options(overlap, overlap_is_default, force, force_is_default):
    def source(is_default):
        return "default" if is_default else "user-specified"
    logger.info(f"-- overlap with target eigenstate = {overlap:g} ({source(overlap_is_default)})")
    if force:
        logger.info("-- report the smallest of N runs (K = 1): forced by "
                    f"algorithm.force_lowest_estimate ({source(force_is_default)})")
    else:
        logger.info("-- report the K-th smallest of N runs: K chosen by sizing "
                    f"({source(force_is_default)})")


def smallest(k: int) -> str:
    """'smallest', '2nd-smallest', '3rd-smallest', '11th-smallest', ..."""
    if k == 1:
        return "smallest"
    suffix = "th" if k % 100 in (11, 12, 13) else {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")
    return f"{k}{suffix}-smallest"


def log_multi_run_plan(window_cls, plan, error_bound: str):
    """Log the failure probability for 1 run and for the plan's N runs, then for N runs at the
    overlaps in OVERLAP_TABLE. `error_bound` is e.g. "energy error <= 0.2"."""
    runs, which = plan.num_runs, smallest(plan.quantile)
    times = "1 time" if runs == 1 else f"{runs} times"
    logger.info("-- if you run the circuit 1 time:")
    logger.info(f"--     probability of failure <= {plan.single_run_failure:.2e}, {error_bound}")
    logger.info(f"-- if you run the circuit {times} and take the {which} estimate:")
    logger.info(f"--     probability of failure <= {plan.multi_run_failure:.2e}, {error_bound}")
    logger.info(f"-- if the true overlap differs from {plan.overlap:g} (same {runs} runs, {which} "
                f"estimate; {error_bound} in every row):")
    logger.info("--     overlap   probability of failure")
    for overlap in OVERLAP_TABLE:
        failure = multi_run_failure(window_cls, plan.num_precision_qubits,
                                    plan.num_confidence_qubits, runs, plan.quantile, overlap)
        logger.info(f"--     {overlap:<9.2f} <= {failure:.2e}")


def counts_text(num_precision, num_confidence):
    return (f"{num_precision} precision + {num_confidence} confidence qubits = "
            f"{num_precision + num_confidence} phase qubits")


def size_tapered_window(config_algorithm, window_cls, phase_error, probability_of_failure):
    """`(num_precision_qubits, num_confidence_qubits)` from the multi-run sizing, logging the
    plan. If no plan exists, logs an error and falls back on `window_cls.from_requirements` (a
    single run, as if the overlap were 1)."""
    overlap, overlap_is_default, force, force_is_default = multi_run_options(config_algorithm)
    log_multi_run_options(overlap, overlap_is_default, force, force_is_default)
    plan = size_tapered_qpe(window_cls, phase_error, probability_of_failure, overlap, force)
    if plan is None:
        fallback = window_cls.from_requirements(phase_error, probability_of_failure)
        logger.error(
            f"Tapered QPE sizing found no plan meeting probability of failure "
            f"{probability_of_failure} at overlap {overlap:g} with cost N x 2^(phase qubits) <= "
            f"{COST_CEILING:g}. FALLBACK: sizing for a single run with overlap 1; the number of "
            "runs and the overlap table are skipped.")
        logger.info("-- FALLBACK tapered QPE sizing: " + counts_text(
            fallback.num_precision_qubits, fallback.num_confidence_qubits))
        logger.info(SINGLE_RUN_NOTE)
        return fallback.num_precision_qubits, fallback.num_confidence_qubits
    logger.info("-- tapered QPE sizing: "
                + counts_text(plan.num_precision_qubits, plan.num_confidence_qubits))
    log_multi_run_plan(window_cls, plan, f"energy error <= {config_algorithm.energy_error:.3g}")
    logger.info(SINGLE_RUN_NOTE)
    return plan.num_precision_qubits, plan.num_confidence_qubits


def report_given_tapered_counts(config_algorithm, unitary, window_cls, window):
    """Log the multi-run plan for user-given tapered qubit counts: the fewest runs meeting
    algorithm.probability_of_failure at the overlap. The error bound is the window's
    2^-(num_precision_qubits + 1) turns, as an energy if the unitary can convert it."""
    num_precision, num_confidence = window.num_precision_qubits, window.num_confidence_qubits
    overlap, overlap_is_default, force, force_is_default = multi_run_options(config_algorithm)
    phase_bound = 2.0**-(num_precision + 1)
    error_bound = f"phase error <= {phase_bound:.3g} turns"
    if hasattr(unitary, "phase_error_from_energy_error"):
        energy_bound = phase_bound / unitary.phase_error_from_energy_error(1.0)
        error_bound = f"energy error <= {energy_bound:.3g}"
        energy_error = config_algorithm.energy_error
        if energy_error is not None and energy_bound > energy_error * (1 + 1e-12):
            logger.warning(
                f"num_precision_qubits = {num_precision} bounds the energy error by "
                f"{energy_bound:.3g}, more than algorithm.energy_error = {energy_error}.")

    logger.info("-- tapered QPE: user-specified " + counts_text(num_precision, num_confidence))
    probability_of_failure = config_algorithm.probability_of_failure
    if probability_of_failure is None:
        warn_unused_multi_run_options(config_algorithm, "no algorithm.probability_of_failure")
        logger.info("-- algorithm.probability_of_failure is not set: skipping the number of runs "
                    "and the overlap table.")
        logger.info(SINGLE_RUN_NOTE)
        return
    if not 0 < probability_of_failure < 1:
        raise ValueError(
            f"algorithm.probability_of_failure must be in (0, 1), got {probability_of_failure}.")
    log_multi_run_options(overlap, overlap_is_default, force, force_is_default)
    plan = fewest_runs(window_cls, num_precision, num_confidence, probability_of_failure, overlap,
                       force)
    if plan is None:
        logger.error(
            f"No number of runs of the tapered QPE with {num_precision} precision and "
            f"{num_confidence} confidence qubits meets probability of failure "
            f"{probability_of_failure} at overlap {overlap:g} with cost N x 2^(phase qubits) <= "
            f"{COST_CEILING:g}; the number of runs and the overlap table are skipped.")
    else:
        log_multi_run_plan(window_cls, plan, error_bound)
    logger.info(SINGLE_RUN_NOTE)


def build_qpe_qhat_flexible(
        config_algorithm: AlgorithmConfiguration,
        unitary):

    logger.verbose("Build a QPE algorithm with QHAT's FlexibleQPE.")

    ancilla_prep_name = (config_algorithm.ancilla_prep or "rectangular").lower()
    if ancilla_prep_name not in ANCILLA_PREPS:
        raise ValueError(f"Invalid QPE ancilla_prep \"{config_algorithm.ancilla_prep}\".")
    qft_inv_name = (config_algorithm.qft_inv or "textbook").lower()
    if qft_inv_name not in INVERSE_QFTS:
        raise ValueError(f"Invalid QPE qft_inv \"{config_algorithm.qft_inv}\".")

    window_cls = ANCILLA_PREPS[ancilla_prep_name]
    tapered = issubclass(window_cls, TaperedWindowStateBase)
    if not tapered:
        warn_unused_multi_run_options(
            config_algorithm, f"algorithm.ancilla_prep = \"{ancilla_prep_name}\"")
    components = dict(ancilla_prep=window_cls, qft_inv=INVERSE_QFTS[qft_inv_name])
    counts = dict(num_phase_qubits=config_algorithm.num_phase_qubits,
                  num_precision_qubits=config_algorithm.num_precision_qubits,
                  num_confidence_qubits=config_algorithm.num_confidence_qubits)
    if any(count is not None for count in counts.values()):
        given = ", ".join(f"{name} = {count}" for name, count in counts.items()
                          if count is not None)
        logger.verbose(f"-- using user-specified qubit counts ({given})")
        try:
            qpe = FlexibleQPE.from_qubit_counts(unitary, **counts, **components)
        except ValueError as error:
            raise ValueError(
                f"algorithm.ancilla_prep = \"{ancilla_prep_name}\" with {given}: {error} Or set "
                "algorithm.energy_error and algorithm.probability_of_failure, and no qubit "
                "counts.") from error
        if tapered:
            report_given_tapered_counts(config_algorithm, unitary, window_cls, qpe.ancilla_prep)
    elif tapered:
        phase_error, probability_of_failure = qpe_requirements(config_algorithm, unitary)
        num_precision, num_confidence = size_tapered_window(
            config_algorithm, window_cls, phase_error, probability_of_failure)
        qpe = FlexibleQPE.from_qubit_counts(
            unitary, num_precision_qubits=num_precision, num_confidence_qubits=num_confidence,
            **components)
    else:
        phase_error, probability_of_failure = qpe_requirements(config_algorithm, unitary)
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

def build_algorithm(
        config_algorithm: AlgorithmConfiguration,
        unitary):

    logger.info("Beginning to construct quantum algorithm.")

    if config_algorithm.method.lower() != "qpe: qhat flexible":
        warn_unused_multi_run_options(
            config_algorithm, f"algorithm.method = \"{config_algorithm.method}\"")
    if config_algorithm.method.lower() in ("qpe: qualtran textbook",):
        return build_qpe_qualtran_textbook(config_algorithm, unitary)
    elif config_algorithm.method.lower() in ("qpe: qhat flexible",):
        return build_qpe_qhat_flexible(config_algorithm, unitary)
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

    Elo = require_real_scalar("lower energy bound", Elo)
    Ehi = require_real_scalar("upper energy bound", Ehi)
    dE = config_algorithm.energy_error
    if dE is None:
        return (Elo, Ehi)
    dE = require_real_scalar("algorithm.energy_error", dE)

    P0 = math.ceil(math.log2((Ehi - Elo) / dE))
    width = 2**P0 * dE
    middle = (Elo + Ehi) / 2
    lo, hi = middle - width / 2, middle + width / 2
    logger.verbose(f"-- P0 = {P0} bits of {dE} each, window = [{lo}, {hi})")
    return (lo, hi)
