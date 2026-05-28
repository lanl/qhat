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

def _compute_exact_matrix(hamiltonian, config_analysis):
    """
    Compute the exact matrix representation of the Hamiltonian.

    This computes the Hamiltonian matrix without any approximations (no Trotter,
    no double-factorization). The choice between dense and sparse/matrix-free
    representation is based on the memory threshold in config_analysis.

    Parameters:
        hamiltonian: The Hamiltonian object
        config_analysis: Analysis configuration with matrix_memory_threshold_gb

    Returns:
        Dense numpy array (small systems) or PauliStringOperator (large systems)

    Raises:
        Exception: If matrix computation fails
    """
    logger.verbose("Computing exact Hamiltonian matrix...")
    try:
        return hamiltonian.to_matrix(
            memory_threshold_gb=config_analysis.matrix_memory_threshold_gb
        )
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

def eigendecomposition_analysis(
        config_analysis: AnalysisConfiguration,
        hamiltonian,
        algorithm,
        exact_matrix=None,
        unitary_matrix=None) -> dict:
    """
    Compute eigendecompositions of exact and/or approximate matrices.

    Parameters:
        config_analysis: Analysis configuration with eigendecomposition settings
        hamiltonian: Hamiltonian object (needed if exact_matrix not provided and required)
        algorithm: Algorithm bloq (needed if unitary_matrix not provided and required)
        exact_matrix: Pre-computed exact matrix (optional)
        unitary_matrix: Pre-computed unitary matrix (optional)

    Returns:
        Dictionary with eigendecomposition results and file paths
    """
    from qhat.analysis.file_io import save_eigendecomposition
    from qhat.analysis.matrix_operations import PauliStringOperator
    import scipy.linalg
    import scipy.sparse.linalg

    num_eigenvalues = config_analysis.num_eigenvalues
    which_matrices = config_analysis.eigendecomposition_matrices
    which_eigs = config_analysis.which_eigenvalues

    logger.info(f"Starting eigendecomposition analysis")
    logger.verbose(f"  num_eigenvalues: {num_eigenvalues}")
    logger.verbose(f"  eigendecomposition_matrices: {which_matrices}")
    logger.verbose(f"  which_eigenvalues: {which_eigs}")

    # Determine if full or partial eigendecomposition
    is_full = isinstance(num_eigenvalues, str) and num_eigenvalues.lower() == "all"

    # Determine which matrices we need
    need_exact = which_matrices in ['exact', 'both']
    need_approx = which_matrices in ['approximate', 'both']

    results = {}

    # Helper function to perform eigendecomposition
    def _eigendecompose(matrix, matrix_type, num_qubits):
        dimension = 2 ** num_qubits
        is_matrix_free = isinstance(matrix, PauliStringOperator) or hasattr(matrix, 'matvec')

        if is_full:
            # Full eigendecomposition
            if is_matrix_free or dimension > 32768:
                raise ValueError(
                    f"Full eigendecomposition not supported for systems with >15 qubits "
                    f"(dimension={dimension}). Use num_eigenvalues=k for partial decomposition."
                )
            logger.verbose(f"Computing full eigendecomposition for {matrix_type} matrix")
            eigenvalues, eigenvectors = scipy.linalg.eigh(matrix)
            num_computed = len(eigenvalues)

        else:
            # Partial eigendecomposition using sparse methods
            k = int(num_eigenvalues)
            if k <= 0:
                raise ValueError(f"num_eigenvalues must be positive, got {k}")
            if k >= dimension:
                raise ValueError(
                    f"num_eigenvalues ({k}) must be less than dimension ({dimension}). "
                    f"Use num_eigenvalues='all' for full decomposition."
                )

            # Map user-friendly values to scipy's 'which' parameter
            which_map = {
                'smallest': 'SA',  # Smallest Algebraic (most negative)
                'largest': 'LA',   # Largest Algebraic (most positive)
            }

            if which_eigs == 'both':
                # Compute both smallest and largest
                logger.verbose(f"Computing {k} smallest and {k} largest eigenvalues for {matrix_type} matrix")
                eigs_small, vecs_small = scipy.sparse.linalg.eigsh(
                    matrix, k=k, which='SA', return_eigenvectors=True
                )
                eigs_large, vecs_large = scipy.sparse.linalg.eigsh(
                    matrix, k=k, which='LA', return_eigenvectors=True
                )
                # Concatenate results
                eigenvalues = np.concatenate([eigs_small, eigs_large])
                eigenvectors = np.concatenate([vecs_small, vecs_large], axis=1)
                num_computed = len(eigenvalues)
            else:
                # Compute only one set
                which_scipy = which_map.get(which_eigs)
                if which_scipy is None:
                    raise ValueError(
                        f"which_eigenvalues must be 'smallest', 'largest', or 'both', "
                        f"got '{which_eigs}'"
                    )
                logger.verbose(f"Computing {k} {which_eigs} eigenvalues for {matrix_type} matrix")
                eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(
                    matrix, k=k, which=which_scipy, return_eigenvectors=True
                )
                num_computed = len(eigenvalues)

        logger.info(f"Computed {num_computed} eigenvalues for {matrix_type} matrix")
        logger.verbose(f"  Eigenvalue range: [{eigenvalues.min():.6e}, {eigenvalues.max():.6e}]")

        return eigenvalues, eigenvectors, num_computed

    # Get git hash for provenance
    from qhat.analysis.config_types import _get_git_hash
    git_hash = _get_git_hash()

    # Compute exact eigendecomposition if requested
    if need_exact:
        logger.info("Computing exact matrix eigendecomposition")

        # Get exact matrix if not provided
        if exact_matrix is None:
            if hamiltonian is None:
                raise ValueError("hamiltonian required for exact eigendecomposition")
            exact_matrix = _compute_exact_matrix(hamiltonian)

        num_qubits = hamiltonian.num_qubits()
        eigs, vecs, num_computed = _eigendecompose(exact_matrix, 'exact', num_qubits)

        # Save to file
        output_file = 'exact_eigendecomposition.npz'
        save_eigendecomposition(
            output_file, eigs, vecs, 'exact',
            num_eigenvalues, which_eigs, git_hash=git_hash
        )

        results['exact_eigendecomposition'] = {
            'file': output_file,
            'num_eigenvalues': num_computed,
            'eigenvalue_range': [float(eigs.min()), float(eigs.max())],
            'which': which_eigs
        }

    # Compute approximate eigendecomposition if requested
    if need_approx:
        logger.info("Computing approximate matrix eigendecomposition")

        # Get unitary matrix if not provided
        if unitary_matrix is None:
            if algorithm is None:
                raise ValueError("algorithm required for approximate eigendecomposition")
            unitary_matrix = _compute_unitary_matrix(algorithm)

        # Infer num_qubits from matrix dimension
        dimension = unitary_matrix.shape[0]
        num_qubits = int(np.log2(dimension))

        eigs, vecs, num_computed = _eigendecompose(unitary_matrix, 'approximate', num_qubits)

        # Save to file
        output_file = 'approximate_eigendecomposition.npz'
        save_eigendecomposition(
            output_file, eigs, vecs, 'approximate',
            num_eigenvalues, which_eigs, git_hash=git_hash
        )

        results['approximate_eigendecomposition'] = {
            'file': output_file,
            'num_eigenvalues': num_computed,
            'eigenvalue_range': [float(eigs.min()), float(eigs.max())],
            'which': which_eigs
        }

    return results

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
    num_eigenvalues = config_analysis.num_eigenvalues
    eigendecomposition_requested = (
        isinstance(num_eigenvalues, int) and num_eigenvalues > 0
    ) or (
        isinstance(num_eigenvalues, str) and num_eigenvalues.lower() == "all"
    )

    if (config_analysis.resource_estimator is None and
        config_analysis.matrix_output_file is None and
        config_analysis.numerical_simulation_inputs is None and
        config_analysis.exact_matrix_output_file is None and
        not eigendecomposition_requested):
        raise ValueError(
            "No analyses requested. Set at least one of:\n"
            "  - resource_estimator (e.g., 'pyliqtr', 'cirq')\n"
            "  - matrix_output_file (e.g., 'matrix.npz', 'matrix.h5', 'matrix.txt')\n"
            "  - numerical_simulation_inputs (e.g., 'state.npy' or ['state1.npy', 'state2.npy'])\n"
            "  - exact_matrix_output_file (e.g., 'exact_hamiltonian.npz')\n"
            "  - num_eigenvalues (e.g., 5 or 'all')"
        )

    results = {}

    # Check if any analysis needs the unitary matrix
    needs_matrix = (
        config_analysis.matrix_output_file is not None or
        config_analysis.numerical_simulation_inputs is not None or
        (eigendecomposition_requested and
         config_analysis.eigendecomposition_matrices in ['approximate', 'both'])
    )

    # Check if any analysis needs the exact Hamiltonian matrix
    needs_exact_matrix = (
        config_analysis.exact_matrix_output_file is not None or
        (eigendecomposition_requested and
         config_analysis.eigendecomposition_matrices in ['exact', 'both'])
    )

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
        exact_matrix = _compute_exact_matrix(hamiltonian, config_analysis)

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

    if eigendecomposition_requested:
        logger.info("Performing eigendecomposition analysis.")
        results["eigendecomposition"] = eigendecomposition_analysis(
            config_analysis, hamiltonian, algorithm,
            exact_matrix=exact_matrix,
            unitary_matrix=unitary_matrix
        )

    # TODO: Add error estimation
    # TODO: Add an option for detailed error analysis (explicitly compute the eigenvalues of the
    #       original Hamiltonian and the final unitary, compute ground state energy from both,
    #       compare the results; will only work for small systems)
    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    logger.info("Algorithm analysis complete.")

    return results
