from qre_types import GeneralConfiguration, AnalysisConfiguration

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from qualtran.resource_counting import get_cost_value, QubitCount

# -------------------------------------------------------------------------------------------------

def resource_estimation_cirq(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:
    raise NotImplementedError

# -------------------------------------------------------------------------------------------------

def resource_estimation_pyliqtr(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:

    config_general.log_verbose("Estimating resources with pyLIQTR.")

    # TODO: rotation error
    #       -- argument rotation_gate_precision sets the precision for a single rotation gate
    #       -- argument circuit_precision sets the precision for the whole circuit (i.e., it sets
    #          rotation_gate_precision to circuit_precision / number of rotations)
    # TODO: profile?
    #       -- argument profile = True: keep rotations as a separate count
    #       -- argument profile = False: estimate rotations as Clifford+T
    resources = estimate_pyliqtr(qpe_circuit)

    resource_dict = {
        "Clifford_count" : resources["Clifford"],
        "T_count"        : resources["T"],
        }
    if "LogicalQubits" in resources:
        resource_dict["qubit_count"] = resources["LogicalQubits"]
    else:
        get_cost_value(qpe_circuit, QubitCount())
    return resource_dict

# -------------------------------------------------------------------------------------------------

def estimate_resources(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:

    if config_analysis.resource_estimator.lower() == "pyliqtr":
        return resource_estimation_pyliqtr(config_general, config_analysis, qpe_circuit)
    elif config_analysis.resource_estimator.lower() == "cirq":
        return resource_estimation_cirq(config_general, config_analysis, qpe_circuit)
    else:
        raise ValueError(
                f"Invalid resource estimator method \"{config_analysis.resource_estimator}\".")

# -------------------------------------------------------------------------------------------------

def analyze_circuit(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:

    config_general.log("Beginning circuit analysis.")

    results = {}
    results["resource_estimates"] = estimate_resources(
            config_general, config_analysis, qpe_circuit)
    # TODO: Add error estimation
    # TODO: Add an option for detailed error analysis (explicitly compute the eigenvalues of the
    #       original Hamiltonian and the final unitary, compute ground state energy from both,
    #       compare the results; will only work for small systems)
    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    return results
