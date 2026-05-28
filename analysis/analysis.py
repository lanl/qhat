import numpy as np
from datetime import datetime
import logging
from pathlib import Path

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from qualtran.resource_counting import get_cost_value, QubitCount

from qhat.analysis.config_types import AnalysisConfiguration, GeneralConfiguration
from qhat.analysis.file_io import save_matrix, load_state, save_state

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

def estimate_resources(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:

    if config_analysis.resource_estimator.lower() == "pyliqtr":
        return resource_estimation_pyliqtr(config_analysis, algorithm)
    elif config_analysis.resource_estimator.lower() == "cirq":
        return resource_estimation_cirq(config_analysis, algorithm)
    else:
        raise ValueError(
                f"Invalid resource estimator method \"{config_analysis.resource_estimator}\".")

# -------------------------------------------------------------------------------------------------

def _compute_unitary_matrix(algorithm):
    """
    Compute the unitary matrix representation of the algorithm.

    Parameters:
        algorithm: The algorithm bloq to analyze

    Returns:
        The unitary matrix as a numpy array

    Raises:
        AttributeError: If algorithm doesn't have tensor_contract() method
        Exception: If tensor_contract() fails
    """
    if not hasattr(algorithm, 'tensor_contract'):
        raise AttributeError(
            f"Algorithm of type {type(algorithm).__name__} does not have a "
            "'tensor_contract()' method. Cannot generate unitary matrix."
        )

    logger.verbose("Computing unitary matrix via tensor contraction...")
    try:
        return algorithm.tensor_contract()
    except Exception as e:
        logger.info(
            f"ERROR: Failed to compute unitary matrix: {e}\n"
            "This may indicate a bug in the algorithm's tensor_contract() implementation."
        )
        raise

# -------------------------------------------------------------------------------------------------

def output_unitary_matrix(
        config_analysis: AnalysisConfiguration,
        algorithm,
        unitary_matrix) -> dict:
    """
    Generate and save the unitary matrix representation of the algorithm.

    Parameters:
        config_analysis: Analysis configuration with matrix_output_format and matrix_output_file
        algorithm: The algorithm bloq to analyze
        unitary_matrix: The unitary matrix to save (pre-computed)

    Returns:
        Dictionary with matrix metadata: shape, file, format, unitarity_error, norm
    """

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

def _compute_exact_matrix(hamiltonian):
    """
    Compute the exact matrix representation of the Hamiltonian.

    This computes the Hamiltonian matrix without any approximations (no Trotter,
    no double-factorization). For small systems (≤15 qubits), returns a dense
    matrix. For larger systems, returns a matrix-free operator.

    Parameters:
        hamiltonian: The Hamiltonian object

    Returns:
        Dense numpy array (small systems) or PauliStringOperator (large systems)

    Raises:
        Exception: If matrix computation fails
    """
    logger.verbose("Computing exact Hamiltonian matrix...")
    try:
        return hamiltonian.to_matrix()
    except Exception as e:
        logger.info(
            f"ERROR: Failed to compute exact Hamiltonian matrix: {e}\n"
            "This may indicate an issue with the Hamiltonian representation."
        )
        raise

# -------------------------------------------------------------------------------------------------

def exact_matrix_output(
        config_analysis: AnalysisConfiguration,
        hamiltonian,
        exact_matrix) -> dict:
    """
    Save the exact Hamiltonian matrix representation.

    Parameters:
        config_analysis: Analysis configuration with exact_matrix_output_file
        hamiltonian: The Hamiltonian object
        exact_matrix: The exact matrix to save (pre-computed)

    Returns:
        Dictionary with matrix metadata: shape, file, format, hermiticity_error, norm

    Note:
        For large systems, exact_matrix may be a matrix-free operator rather than
        a dense array. In that case, certain properties (like saving to file) may
        not be supported or may require special handling.
    """
    from qhat.analysis.matrix_operations import PauliStringOperator

    # Check if this is a matrix-free operator
    is_matrix_free = isinstance(exact_matrix, PauliStringOperator)

    if is_matrix_free:
        logger.verbose(f"Matrix-free operator with shape: {exact_matrix.shape}")
        logger.info(
            "WARNING: Matrix-free operator cannot be directly saved to file. "
            "Skipping matrix output for large system."
        )
        return {
            'matrix_shape': exact_matrix.shape,
            'matrix_file': None,
            'matrix_format': None,
            'hermiticity_error': None,
            'matrix_norm': None,
            'matrix_free': True,
            'note': 'Matrix-free operator not saved (too large)'
        }

    # For dense matrices, proceed with normal output
    logger.verbose(f"Matrix shape: {exact_matrix.shape}")
    logger.verbose(f"Matrix dtype: {exact_matrix.dtype}")

    # Compute Hermiticity check: ||H - H†||_F
    try:
        matrix_norm = np.linalg.norm(exact_matrix, ord='fro')
        H_dag = np.conj(exact_matrix.T)
        hermiticity_error = np.linalg.norm(exact_matrix - H_dag, ord='fro')
        logger.verbose(f"Matrix Frobenius norm: {matrix_norm:.6e}")
        logger.verbose(f"Hermiticity error ||H - H†||_F: {hermiticity_error:.6e}")
    except Exception as e:
        logger.info(f"WARNING: Could not compute Hermiticity check: {e}")
        matrix_norm = None
        hermiticity_error = None

    # Save matrix to file (format auto-detected from extension)
    output_file = config_analysis.exact_matrix_output_file
    save_matrix(
        output_file, exact_matrix,
        hermiticity_error=hermiticity_error,
        matrix_norm=matrix_norm
    )

    # Return metadata
    return {
        'matrix_shape': exact_matrix.shape,
        'matrix_file': str(output_file),
        'matrix_format': Path(output_file).suffix,
        'hermiticity_error': float(hermiticity_error) if hermiticity_error is not None else None,
        'matrix_norm': float(matrix_norm) if matrix_norm is not None else None,
        'matrix_free': False
    }

# -------------------------------------------------------------------------------------------------

def numerical_simulation(
        config_analysis: AnalysisConfiguration,
        algorithm,
        unitary_matrix) -> dict:
    """
    Perform numerical simulation by applying the unitary matrix to input state(s).

    Parameters:
        config_analysis: Analysis configuration with numerical_simulation_inputs
        algorithm: The algorithm bloq to analyze
        unitary_matrix: The unitary matrix to apply (pre-computed)

    Returns:
        Dictionary with simulation metadata: list of {input_file, output_file, output_norm}
    """

    # Log matrix properties
    logger.verbose(f"Matrix shape: {unitary_matrix.shape}")

    # Normalize input to list
    inputs = config_analysis.numerical_simulation_inputs
    if isinstance(inputs, str):
        input_files = [inputs]
    elif isinstance(inputs, list):
        input_files = inputs
    else:
        raise ValueError(
            f"numerical_simulation_inputs must be a string or list of strings, "
            f"got {type(inputs)}"
        )

    logger.info(f"Running numerical simulation on {len(input_files)} input state(s)")

    results = []

    for input_file in input_files:
        logger.verbose(f"Processing {input_file}")

        # Load initial state
        try:
            initial_state = load_state(input_file)
        except Exception as e:
            logger.info(f"ERROR: Failed to load state from {input_file}: {e}")
            raise

        # Validate dimensions
        if initial_state.shape[0] != unitary_matrix.shape[1]:
            raise ValueError(
                f"Dimension mismatch: state vector has dimension {initial_state.shape[0]} "
                f"but matrix expects {unitary_matrix.shape[1]}"
            )

        # Perform matrix-vector multiplication
        logger.verbose("Performing matrix-vector multiplication")
        final_state = unitary_matrix @ initial_state

        # Compute norm
        final_norm = np.linalg.norm(final_state)
        logger.verbose(f"Final state norm: {final_norm:.6e}")

        # Generate output filename: input.npy -> input_final.npy
        input_path = Path(input_file)
        output_file = str(input_path.parent / f"{input_path.stem}_final{input_path.suffix}")

        # Save final state
        try:
            save_state(output_file, final_state)
        except Exception as e:
            logger.info(f"ERROR: Failed to save state to {output_file}: {e}")
            raise

        logger.info(f"Simulation complete: {input_file} -> {output_file}")

        results.append({
            'input_file': input_file,
            'output_file': output_file,
            'output_norm': float(final_norm)
        })

    return {'simulations': results}

# -------------------------------------------------------------------------------------------------

def analyze_algorithm(
        config_analysis: AnalysisConfiguration,
        algorithm,
        hamiltonian=None) -> dict:

    logger.info("Beginning algorithm analysis.")

    # Validate at least one analysis requested
    if (config_analysis.resource_estimator is None and
        config_analysis.matrix_output_file is None and
        config_analysis.numerical_simulation_inputs is None and
        config_analysis.exact_matrix_output_file is None):
        raise ValueError(
            "No analyses requested. Set at least one of:\n"
            "  - resource_estimator (e.g., 'pyliqtr', 'cirq')\n"
            "  - matrix_output_file (e.g., 'matrix.npz', 'matrix.h5', 'matrix.txt')\n"
            "  - numerical_simulation_inputs (e.g., 'state.npy' or ['state1.npy', 'state2.npy'])\n"
            "  - exact_matrix_output_file (e.g., 'exact_hamiltonian.npz')"
        )

    results = {}

    # Check if any analysis needs the unitary matrix
    needs_matrix = (config_analysis.matrix_output_file is not None or
                    config_analysis.numerical_simulation_inputs is not None)

    # Check if any analysis needs the exact Hamiltonian matrix
    needs_exact_matrix = (config_analysis.exact_matrix_output_file is not None)

    # Compute matrices once if needed
    unitary_matrix = None
    if needs_matrix:
        unitary_matrix = _compute_unitary_matrix(algorithm)

    exact_matrix = None
    if needs_exact_matrix:
        if hamiltonian is None:
            raise ValueError(
                "Exact matrix computation requires hamiltonian parameter. "
                "Pass hamiltonian to analyze_algorithm()."
            )
        exact_matrix = _compute_exact_matrix(hamiltonian)

    # Dispatch to requested analyses
    if config_analysis.resource_estimator is not None:
        logger.info(f"Performing resource estimation using {config_analysis.resource_estimator}.")
        results["resource_estimates"] = estimate_resources(config_analysis, algorithm)

    if config_analysis.matrix_output_file is not None:
        logger.info("Generating unitary matrix output.")
        results["matrix_output"] = output_unitary_matrix(config_analysis, algorithm, unitary_matrix)

    if config_analysis.exact_matrix_output_file is not None:
        logger.info("Generating exact Hamiltonian matrix output.")
        results["exact_matrix_output"] = exact_matrix_output(config_analysis, hamiltonian, exact_matrix)

    if config_analysis.numerical_simulation_inputs is not None:
        logger.info("Performing numerical simulation.")
        results["numerical_simulation"] = numerical_simulation(config_analysis, algorithm, unitary_matrix)

    # TODO: Add error estimation
    # TODO: Add an option for detailed error analysis (explicitly compute the eigenvalues of the
    #       original Hamiltonian and the final unitary, compute ground state energy from both,
    #       compare the results; will only work for small systems)
    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    logger.info("Algorithm analysis complete.")

    return results
