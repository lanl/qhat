import logging

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from qualtran.resource_counting import get_cost_value, QubitCount
from qualtran.cirq_interop.t_complexity_protocol import t_complexity

from qhat.analysis.config_types import AnalysisConfiguration

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def resource_estimation_cirq(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:
    raise NotImplementedError

# -------------------------------------------------------------------------------------------------

def resource_estimation_pyliqtr(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:

    logger.verbose("Estimating resources with pyLIQTR.")

    # TODO: rotation error
    #       -- argument rotation_gate_precision sets the precision for a single rotation gate
    #       -- argument algorithm_precision sets the precision for the whole algorithm (i.e., it
    #          sets rotation_gate_precision to algorithm_precision / number of rotations)
    # TODO: profile?
    #       -- argument profile = True: keep rotations as a separate count
    #       -- argument profile = False: estimate rotations as Clifford+T
    resources = estimate_pyliqtr(algorithm)

    resource_dict = {
        "Clifford_count" : resources["Clifford"],
        "T_count"        : resources["T"],
        }
    if "LogicalQubits" in resources:
        resource_dict["qubit_count"] = resources["LogicalQubits"]
    else:
        get_cost_value(algorithm, QubitCount())
    return resource_dict

# -------------------------------------------------------------------------------------------------

def resource_estimation_qualtran(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:

    logger.verbose("Estimating resources with Qualtran.")

    t_cliff_resources = t_complexity(algorithm)
    qubits = get_cost_value(algorithm, QubitCount())

    # TODO:  add rotation count to this?
    resource_dict = {
        "Clifford_count" : t_cliff_resources.clifford,
        "T_count"        : t_cliff_resources.t,
        "qubit_count"    : qubits
        }
    return resource_dict

# -------------------------------------------------------------------------------------------------

def estimate_resources(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:

    if config_analysis.resource_estimator.lower() == "pyliqtr":
        return resource_estimation_pyliqtr(config_analysis, algorithm)
    elif config_analysis.resource_estimator.lower() == "cirq":
        return resource_estimation_cirq(config_analysis, algorithm)
    elif config_analysis.resource_estimator.lower() == "qualtran":
        return resource_estimation_qualtran(config_analysis, algorithm)
    else:
        raise ValueError(
                f"Invalid resource estimator method \"{config_analysis.resource_estimator}\".")
