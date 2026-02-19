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
# I once traced PhaseEstimation and it does not follow the "usual" phase estimation circuit design.
# The QubitizedPhaseEstimation involves quantum walk operators (add that to my list of things to
# learn about), and it's not clear if it follows the usual phase estimation circuit or not.
# 
# This Qualtran issue provides nice links to the "standard" version and the "walk" version
# (although at a quick glance it looks as though both use multiple ancilla qubits, while both
# pyLIQTR versions use a single ancilla qubit): https://github.com/quantumlib/Qualtran/issues/819.
# This brings us to a total of four versions of QPE, and it's not clear to me how much overlap
# there is between the methods.
# 
# I can probably implement the "standard" QPE circuit in a useful framework, but I will have to
# read more to understand all the variations and figure out how to implement them.
# 
# Addendum: Qualtran provides QPE circuits, including TextbookQPE.  It may not be the most
# efficient (?) but it may be a reliable starting point.

import math

from qre_types import GeneralConfiguration, QPEConfiguration

from pyLIQTR.PhaseEstimation.pe import PhaseEstimation
from pyLIQTR.qubitization.phase_estimation import QubitizedPhaseEstimation
from qualtran.bloqs.phase_estimation import TextbookQPE

# -------------------------------------------------------------------------------------------------

def build_qpe_qualtran_textbook(
        config_general: GeneralConfiguration,
        config_qpe: QPEConfiguration,
        unitary,
        P0):

    config_general.log_verbose("Build a QPE circuit with Qualtran's \"textbook\" method.")

    P = config_qpe.num_phase_qubits
    if P is None:
        assert P0 is not None
        assert config_qpe.probability_of_failure is not None
        Pextra = math.ceil(math.log2(2.0 + 0.5 / config_qpe.probability_of_failure))
        P = P0 + Pextra
        config_general.log_verbose(
                f"-- extending the phase register by {Pextra} qubits (total = {P})")

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
    #             generating the circuit:
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
    return TextbookQPE(unitary, P)

# -------------------------------------------------------------------------------------------------

def build_qpe_pyliqtr_qubitized(
        config_general: GeneralConfiguration,
        config_qpe: QPEConfiguration,
        unitary):

    config_general.log_verbose(
            "Build a QPE circuit with pyLIQTR's \"QubitizedPhaseEstimation\" method.")

    # TODO: The name and signature suggest that this may _only_ be valid for block-encoded
    #       unitaries.  Is that true?
    return QubitizedPhaseEstimation(block_encoding=unitary, prec=config_qpe.num_phase_qubits)

# -------------------------------------------------------------------------------------------------

def build_qpe_circuit(
        config_general: GeneralConfiguration,
        config_qpe: QPEConfiguration,
        unitary,
        P0):

    config_general.log("Beginning to construct quantum phase estimation circuit.")

    if config_qpe.method.lower() in ("qualtran textbook",):
        return build_qpe_qualtran_textbook(config_general, config_qpe, unitary, P0)
    elif config_qpe.method.lower() in ("qualtran qubitization",):
        # TODO: This may be more specialized (for LCU only?), but I'm not yet sure of the details.
        raise NotImplementedError()
    elif config_qpe.method.lower() in ("pyliqtr qubitized",):
        return build_qpe_pyliqtr_qubitized(config_general, config_qpe, unitary)
    else:
        raise ValueError(f"Invalid QPE circuit method \"{config_qpe.method}\".")

# -------------------------------------------------------------------------------------------------

def compute_initial_phase_qubits(
        config_general: GeneralConfiguration,
        config_qpe: QPEConfiguration,
        Elo2, Ehi2):

    config_general.log("Computing initial phase qubits.")

    P0 = math.ceil(math.log2((Ehi2 - Elo2) / config_qpe.energy_error))
    config_general.log_verbose(f"-- initial number of phase qubits = {P0}")
    dE_new = 2**P0 * config_qpe.energy_error
    Elo3 = Elo2
    Ehi3 = Elo3 + dE_new
    config_general.log_verbose(f"-- QPE-optimized bounds = [{Elo3}, {Ehi3})")

    return (P0, Elo3, Ehi3)
