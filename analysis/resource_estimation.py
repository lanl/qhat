import logging

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from qualtran.resource_counting import get_cost_value, QubitCount
from qualtran.cirq_interop.t_complexity_protocol import t_complexity

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.utils import normalize_string_or_list_to_list

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def resource_estimation_cirq(algorithm) -> dict:
    raise NotImplementedError

# -------------------------------------------------------------------------------------------------

def resource_estimation_pyliqtr(algorithm) -> dict:

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

def resource_estimation_qualtran(algorithm) -> dict:

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

def estimate_resources(resource_estimator, algorithm) -> dict:
    """
    Estimate quantum algorithm resources using one or more methods.

    Parameters:
        resource_estimator: String or list of strings specifying estimation method(s)
                           Valid values: 'pyliqtr', 'qualtran', 'cirq'
        algorithm: Algorithm bloq to analyze

    Returns:
        Dictionary with resource estimates for each method
    """
    logger.info(f"Computing resource estimates: {resource_estimator}")

    results = {}
    for estimator in normalize_string_or_list_to_list(resource_estimator):
        estimator_normalized = estimator.lower()

        if estimator_normalized == "pyliqtr":
            logger.verbose(f"Estimating resources with pyLIQTR")
            results['pyliqtr'] = resource_estimation_pyliqtr(algorithm)
        elif estimator_normalized == "cirq":
            logger.verbose(f"Estimating resources with Cirq")
            results['cirq'] = resource_estimation_cirq(algorithm)
        elif estimator_normalized == "qualtran":
            logger.verbose(f"Estimating resources with Qualtran")
            results['qualtran'] = resource_estimation_qualtran(algorithm)
        else:
            raise ValueError(
                f"Invalid resource estimator method \"{estimator}\". "
                f"Valid options: 'pyliqtr', 'qualtran', 'cirq'"
            )

    return results
