import numpy as np
from datetime import datetime
from pathlib import Path

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from qualtran.resource_counting import get_cost_value, QubitCount

from qhat.analysis.config_types import AnalysisConfiguration, GeneralConfiguration
from qhat.analysis.file_io import save_matrix
import logging

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def resource_estimation_cirq(config_analysis: AnalysisConfiguration, algorithm) -> dict:
    raise NotImplementedError

# -------------------------------------------------------------------------------------------------

def resource_estimation_pyliqtr(config_analysis: AnalysisConfiguration, algorithm) -> dict:

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

def estimate_resources(config_analysis: AnalysisConfiguration, algorithm) -> dict:

    if config_analysis.resource_estimator.lower() == "pyliqtr":
        return resource_estimation_pyliqtr(config_analysis, algorithm)
    elif config_analysis.resource_estimator.lower() == "cirq":
        return resource_estimation_cirq(config_analysis, algorithm)
    else:
        raise ValueError(
                f"Invalid resource estimator method \"{config_analysis.resource_estimator}\".")

# -------------------------------------------------------------------------------------------------

def output_unitary_matrix(config_analysis: AnalysisConfiguration, algorithm) -> dict:
    """
    Generate and save the unitary matrix representation of the algorithm.

    Parameters:
        config_analysis: Analysis configuration with matrix_output_format and matrix_output_file
        algorithm: The algorithm bloq to analyze

    Returns:
        Dictionary with matrix metadata: shape, file, format, unitarity_error, norm
    """

    # Check if the algorithm has tensor_contract method
    if not hasattr(algorithm, 'tensor_contract'):
        raise AttributeError(
            f"Algorithm of type {type(algorithm).__name__} does not have a "
            "'tensor_contract()' method. Cannot generate unitary matrix."
        )

    # Extract the unitary matrix
    logger.verbose("Computing unitary matrix via tensor contraction...")
    try:
        unitary_matrix = algorithm.tensor_contract()
    except Exception as e:
        logger.info(
            f"ERROR: Failed to compute unitary matrix: {e}\n"
            "This may indicate a bug in the algorithm's tensor_contract() implementation."
        )
        raise

    # Log basic properties
    logger.verbose(f"Matrix shape: {unitary_matrix.shape}")
    logger.verbose(f"Matrix dtype: {unitary_matrix.dtype}")

    # Compute unitarity check: ||U†U - I||_F
    try:
        matrix_norm = np.linalg.norm(unitary_matrix, ord='fro')
        U_dag_U = np.conj(unitary_matrix.T) @ unitary_matrix
        identity = np.eye(unitary_matrix.shape[0])
        unitarity_error = np.linalg.norm(U_dag_U - identity, ord='fro')
        logger.verbose(f"Matrix Frobenius norm: {matrix_norm:.6e}")
        logger.verbose(f"Unitarity error ||U†U - I||_F: {unitarity_error:.6e}")
    except Exception as e:
        logger.info(f"WARNING: Could not compute unitarity check: {e}")
        matrix_norm = None
        unitarity_error = None

    # Save matrix to file (format auto-detected from extension)
    output_file = config_analysis.matrix_output_file
    save_matrix(
        output_file, unitary_matrix,
        unitarity_error=unitarity_error,
        matrix_norm=matrix_norm
    )

    # Return metadata
    return {
        'matrix_shape': unitary_matrix.shape,
        'matrix_file': str(output_file),
        'matrix_format': Path(output_file).suffix,
        'unitarity_error': float(unitarity_error) if unitarity_error is not None else None,
        'matrix_norm': float(matrix_norm) if matrix_norm is not None else None
    }

# -------------------------------------------------------------------------------------------------

def analyze_algorithm(config_analysis: AnalysisConfiguration, algorithm) -> dict:

    logger.info("Beginning algorithm analysis.")

    # Validate at least one analysis requested
    if (config_analysis.resource_estimator is None and
        config_analysis.matrix_output_file is None):
        raise ValueError(
            "No analyses requested. Set at least one of:\n"
            "  - resource_estimator (e.g., 'pyliqtr', 'cirq')\n"
            "  - matrix_output_file (e.g., 'matrix.npz', 'matrix.h5', 'matrix.txt')"
        )

    results = {}

    # Dispatch to requested analyses
    if config_analysis.resource_estimator is not None:
        logger.info(f"Performing resource estimation using {config_analysis.resource_estimator}.")
        results["resource_estimates"] = estimate_resources(config_analysis, algorithm)

    if config_analysis.matrix_output_file is not None:
        logger.info("Generating unitary matrix output.")
        results["matrix_output"] = output_unitary_matrix(config_analysis, algorithm)

    # TODO: Add error estimation
    # TODO: Add an option for detailed error analysis (explicitly compute the eigenvalues of the
    #       original Hamiltonian and the final unitary, compute ground state energy from both,
    #       compare the results; will only work for small systems)
    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    logger.info("Algorithm analysis complete.")

    return results
