import logging

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from pyLIQTR.utils.resource_analysis import get_T_counts_from_rotations
from qualtran.resource_counting import get_cost_value, QubitCount
from qualtran.cirq_interop.t_complexity_protocol import t_complexity

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.utils import normalize_string_or_list_to_list

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

DEFAULT_ROTATION_GATE_PRECISION = 1e-10    # pyLIQTR's default

def convert_rotations(t_count, clifford_count, rotations,
                      rotation_gate_precision=None, algorithm_precision=None) -> dict:
    """
    Convert the total rotation count to T and Clifford gates, once, at the end (plan D11).

    Each rotation costs the same number of T gates, from pyLIQTR's fit
    `3.02 log2(1/precision) + 0.77` (with its seeded random offset), plus 2 Cliffords per T and one
    more Clifford overall. The precision is `rotation_gate_precision` for each rotation, or
    `algorithm_precision` divided by the total rotation count. If neither is given, the per-rotation
    default is 1e-10. Giving both raises an error.

    The input counts are returned unconverted, as `T_count_unconverted` and
    `Clifford_count_unconverted`; `T_count` and `Clifford_count` include the rotations.
    E.g. 1000 rotations at 1e-10 each: about 100 T each, so `T_count = t_count + 1000 * ~100`.
    """
    if rotation_gate_precision is not None and algorithm_precision is not None:
        raise ValueError("Set only one of analysis.rotation_gate_precision (per rotation) and "
                         "analysis.algorithm_precision (total, divided by the number of rotations).")
    for name, value in (("rotation_gate_precision", rotation_gate_precision),
                        ("algorithm_precision", algorithm_precision)):
        if value is not None and not 0 < value < 1:
            raise ValueError(f"analysis.{name} must be between 0 and 1, got {value}.")

    result = {
        "rotation_count"           : rotations,
        "T_count_unconverted"      : t_count,
        "Clifford_count_unconverted": clifford_count,
        "T_count"                  : t_count,
        "Clifford_count"           : clifford_count,
        }
    if not rotations:
        return result

    if algorithm_precision is not None:
        precision = algorithm_precision / rotations
        t_total = get_T_counts_from_rotations(rotations, circuit_precision=algorithm_precision)
    else:
        precision = (DEFAULT_ROTATION_GATE_PRECISION if rotation_gate_precision is None
                     else rotation_gate_precision)
        t_total = get_T_counts_from_rotations(rotations, gate_precision=precision)
    logger.info(f"Converting {rotations} rotations at precision {precision:g} each: "
                f"{t_total} T gates ({t_total // rotations} per rotation).")
    result["rotation_precision"] = precision
    result["T_count"] = t_count + t_total
    result["Clifford_count"] = clifford_count + 2 * t_total + 1
    return result

# -------------------------------------------------------------------------------------------------

def resource_estimation_cirq(algorithm) -> dict:
    raise NotImplementedError

# -------------------------------------------------------------------------------------------------

def resource_estimation_pyliqtr(algorithm, rotation_gate_precision=None,
                                algorithm_precision=None) -> dict:

    logger.verbose("Estimating resources with pyLIQTR.")

    # profile=True keeps the rotations as a separate count; `convert_rotations` converts them.
    resources = estimate_pyliqtr(algorithm, profile=True)

    resource_dict = convert_rotations(resources["T"], resources["Clifford"], resources["Rotations"],
                                      rotation_gate_precision, algorithm_precision)
    # The qubit count is the peak in use (plan D3). pyLIQTR's LogicalQubits counts registers only.
    try:
        qubits = get_cost_value(algorithm, QubitCount())
    except Exception as error:    # e.g. known issue K1: "QPE: qualtran textbook"
        if "LogicalQubits" not in resources:
            raise
        qubits = resources["LogicalQubits"]
        logger.warning(f"Qualtran QubitCount failed ({type(error).__name__}); using pyLIQTR's "
                       f"LogicalQubits ({qubits}), which counts registers only.")
    else:
        if "LogicalQubits" in resources and resources["LogicalQubits"] != qubits:
            logger.info(f"pyLIQTR LogicalQubits = {resources['LogicalQubits']} counts registers "
                        f"only; qubit_count uses the peak, {qubits} (Qualtran QubitCount).")
    resource_dict["qubit_count"] = qubits
    return resource_dict

# -------------------------------------------------------------------------------------------------

def resource_estimation_qualtran(algorithm, rotation_gate_precision=None,
                                 algorithm_precision=None) -> dict:

    logger.verbose("Estimating resources with Qualtran.")

    t_cliff_resources = t_complexity(algorithm)
    qubits = get_cost_value(algorithm, QubitCount())

    resource_dict = convert_rotations(t_cliff_resources.t, t_cliff_resources.clifford,
                                      t_cliff_resources.rotations,
                                      rotation_gate_precision, algorithm_precision)
    resource_dict["qubit_count"] = qubits
    return resource_dict

# -------------------------------------------------------------------------------------------------

def estimate_resources(resource_estimator, algorithm, rotation_gate_precision=None,
                       algorithm_precision=None) -> dict:
    """
    Estimate quantum algorithm resources using one or more methods.

    Parameters:
        resource_estimator: String or list of strings specifying estimation method(s)
                           Valid values: 'pyliqtr', 'qualtran', 'cirq'
        algorithm: Algorithm bloq to analyze
        rotation_gate_precision, algorithm_precision: how precisely rotations are converted to
                           T gates (see `convert_rotations`); at most one may be set

    Returns:
        Dictionary with resource estimates for each method
    """
    logger.info(f"Computing resource estimates: {resource_estimator}")

    if rotation_gate_precision is not None and algorithm_precision is not None:
        convert_rotations(0, 0, 0, rotation_gate_precision, algorithm_precision)    # raises

    results = {}
    for estimator in normalize_string_or_list_to_list(resource_estimator):
        estimator_normalized = estimator.lower()

        if estimator_normalized == "pyliqtr":
            logger.verbose(f"Estimating resources with pyLIQTR")
            results['pyliqtr'] = resource_estimation_pyliqtr(
                algorithm, rotation_gate_precision, algorithm_precision)
        elif estimator_normalized == "cirq":
            logger.verbose(f"Estimating resources with Cirq")
            results['cirq'] = resource_estimation_cirq(algorithm)
        elif estimator_normalized == "qualtran":
            logger.verbose(f"Estimating resources with Qualtran")
            results['qualtran'] = resource_estimation_qualtran(
                algorithm, rotation_gate_precision, algorithm_precision)
        else:
            raise ValueError(
                f"Invalid resource estimator method \"{estimator}\". "
                f"Valid options: 'pyliqtr', 'qualtran', 'cirq'"
            )

    return results
